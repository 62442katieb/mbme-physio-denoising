import bids
import json
import json
import gzip
import shutil
#import errno
import os
import argparse

import pandas as pd
#import numpy as np

#import heartpy as hp
from systole.plots import plot_raw
from systole.detection import interpolate_clipping, ecg_peaks
from systole.utils import heart_rate

import scipy.signal as signal
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from filtering import comb_band_stop, fourier_freq, plot_signal_fourier

pd.options.mode.chained_assignment = None

parser = argparse.ArgumentParser(description='Accept BIDS directory, specify # slices if slice timing isn\' specified in BOLD sidecar.')
parser.add_argument('dset', type=str,
                    help='Valid BIDS dataset containing physiological data.')
parser.add_argument('--slices', 
                    type=int, 
                    help='''The number of slices acquired in BOLD sequences concurrent with physiological data acquisition.
                            If multiple BOLD sequences were used with different numbers of slices, this option is invalid, please
                            make sure slice timing is specified in the json sidecar for each BOLD sequence.
                    ''')
parser.add_argument('--mb', type=int, default=1,
                    help='Multiband factor of fMRI scan sequence (if single band, --mb=1).')
parser.add_argument('--biopac', action='store_true',
                    help='Run only BIOPAC-recommended filtering at single-band slice collection frequency.')
parser.add_argument('--verbose', action='store_true', 
                    help='Print information as the cleaning script runs.')

args = parser.parse_args()

bids_dir = args.dset
deriv_dir = os.path.join(bids_dir, 'derivatives', 'combed_physio')
if not os.path.exists(deriv_dir):
    os.makedirs(deriv_dir)
dset = bids.BIDSLayout(bids_dir)



dataset_description = {
    "Name": "Combed Physio Data",
    "BIDSVersion": "1.4.0",
    "DatasetType": "derivative",
    "GeneratedBy": [
        {
            "Name": "physio_comb",
            "Version": "0.2.0"
        }
    ],
}

if len(dset.description["DatasetDOI"]) > 0:
    dataset_description["SourceDatasets"] = [
            {
                "DOI": dset.description['DatasetDOI'],
            }
        ]
else:
    pass

physio_jsons = dset.get(suffix='physio', extension='json')
slices = args.slices

if args.mb:
    mb = args.mb
else:
    mb = 1

# check json for cardiac column
# extract sampling frequency if cardiac is present
# and then unzip the accompanying physio tsv.gz 
# and then load the physio.tsv (maybe as a separate step)
physio_files = []
for file in physio_jsons[:5]:
    path = file.path
    with open(path) as json_file:
        data = json.load(json_file)
    if 'cardiac' in data['Columns']:
        data_file = file.get_associations()
        assert len(data_file) == 1, f"Found {len(data_file)} physio files instead of the expected 1."
        data_file = data_file[0].path
        if 'gz' in data_file[-2:]:
            data_tsv = data_file[:-3]
            if os.path.exists(data_tsv):
                pass
            else:
                with gzip.open(data_file, 'rb') as f_in:
                    with open(data_tsv, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
        elif 'tsv' in data_file[:-3]:
            data_tsv = data_file
        else:
            raise FileNotFoundError

        # save entities as variables to pull associated BOLD meta-data
        task = file.entities['task']
        subj = file.entities['subject']
        dtype = file.entities['datatype']

        if 'session' in file.entities.keys():
            sesh = file.entities['session']
        else:
            sesh = None
        if 'run' in file.entities.keys():
            run = file.entities['run']
        else:
            run = None
        
        bold_json = dset.get(task=task, 
                                subject=subj, 
                                session=sesh,
                                datatype=dtype, 
                                run=run, 
                                suffix='bold',
                                extension='json')
        assert len(bold_json) == 1, f"Looking for one associated BOLD file, found {len(bold_json)}."
        bold_json = bold_json[0]
                               
        physio_dict = file.get_dict()
        
        bold_dict = bold_json.get_dict()
        

        out_path = os.path.join(deriv_dir, 
                                f'sub-{subj}')
        if sesh:
            out_path = os.path.join(out_path, 
                                    f'ses-{sesh}',
                                    'func')
            if not os.path.exists(out_path):
                os.makedirs(out_path)
            if run:
                out_path = os.path.join(out_path, 
                                        f'sub-{subj}_ses-{sesh}_task-{task}_run-{run}_desc-filtered_physio')
        else:
            pass
        if run:
            out_path = os.path.join(deriv_dir, 
                                f'sub-{subj}', 
                                f'ses-{sesh}',
                                'func')
            if not os.path.exists(out_path):    
                os.makedirs(out_path)
            out_path = os.path.join(out_path, 
                                    f'sub-{subj}_task-{task}_run-{run}_desc-filtered_physio')
        else:
            out_path = os.path.join(deriv_dir, 
                                    f'sub-{subj}', 
                                    'func')
            if not os.path.exists(out_path):
                os.makedirs(out_path)
            out_path = os.path.join(out_path, 
                                    f'sub-{subj}_task-{task}_desc-filtered_physio')
        
        notches = {}
        # params needed for physio denoising
        tr = bold_dict['RepetitionTime']
        if len(bold_dict['SliceTiming']) < 1:
            pass
        else:
            slices = len(bold_dict['SliceTiming'])
        if 'MultibandAccelerationFactor' in bold_dict.keys():
            mb = bold_dict['MultibandAccelerationFactor']
        else:
            mb = 1
        if 'echo' in bold_json.entities.keys():
            notches['tr'] = 1/tr
            tr_filter = True
        else:
            tr_filter = False
        cutoff = 120

        fs = physio_dict['SamplingFrequency']

        nyquist = fs/2
        Q = 100 
        print(f'tr: {tr}\nmb: {mb}\nslices: {slices}\nfs: {fs}')
        
        # load the data
        dat = pd.read_csv(data_tsv, sep='\t', header=0)
        dat.columns = physio_dict['Columns']
        dat['seconds'] = dat.index / fs
        # and then run the denoising filters 
        notches['slices'] = slices / mb / tr

        if args.biopac:
            notches = {'slices': slices/tr}
        
        # compute power spectrum of raw signal
        fft_ecg, _, freq, flimit = fourier_freq(dat['cardiac'], 1/fs, 60)
        # first, plot the raw signal and its power spectrum
        # how many samples in six seconds?
        lim = 6 * fs
        downsample = 10
        fig = plot_signal_fourier(time=dat['seconds'], 
                    data=dat['cardiac'], 
                    downsample=downsample, 
                    low_lim=0,
                    high_lim=int(lim * downsample), 
                    fft=fft_ecg, 
                    freq=freq, 
                    lim_fmax=flimit, 
                    annotate=False,
                    peaks=None,
                    slice_peaks=None,
                    title='Raw cardiac', 
                    save=True)
        fig.savefig(f'{out_path}.png', dpi=400, bbox_inches='tight')

        # let the filtering begin
        for notch in notches.keys():
            filtered = comb_band_stop(notch, dat['cardiac'], Q, fs)

    else:
        pass