import bids
import json
#import gzip
#import shutil
#import errno
import os
import gc
import argparse

import pandas as pd
#import numpy as np

#import heartpy as hp
#from systole.plots import plot_raw
#from systole.detection import interpolate_clipping, ecg_peaks
#from systole.utils import heart_rate

#import scipy.signal as signal
#import numpy as np
import pandas as pd
#import seaborn as sns
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
parser.add_argument('--tr', type=float, default=None,
                    help='Repetition time of fMRI scan sequence.')
parser.add_argument('--multiecho', type=bool, default=None,
                    help='Does the fMRI scan sequence collect multiple echoes.')
parser.add_argument('--biopac', action='store_true',
                    help='Run only BIOPAC-recommended filtering at single-band slice collection frequency.')
parser.add_argument('--progress', action='store_true',
                    help='Display progress bar (note: requires enligten).')
parser.add_argument('--multicomb', type=bool, default=None,
                    help='Run multiple comb notch filtering strategies and save all outputs.')
parser.add_argument('--verbose', action='store_true', 
                    help='Print information as the cleaning script runs.')

args = parser.parse_args()

bids_dir = args.dset
deriv_dir = os.path.join(bids_dir, 'derivatives', 'PhysioComb')
if not os.path.exists(deriv_dir):
    os.makedirs(deriv_dir)
dset = bids.BIDSLayout(bids_dir)

dataset_description = {
    "Name": "Combed Physio Data",
    "BIDSVersion": "1.4.0",
    "DatasetType": "derivative",
    "GeneratedBy": [
        {
            "Name": "PhysioComb",
            "Version": "0.2.0"
        }
    ],
}

if args.multiecho:
    me = args.multiecho
else:
    me = False
if args.tr:
    tr = args.tr
else:
    tr = None
if args.slices:
    slices = args.slices
else:
    slices = None
if args.mb:
    mb = args.mb
else:
    mb = 1

CANDIDATES = ['cardiac', 'ecg', 'ekg', # are there any ecg signals? note: doesn't distinguish between ppg and ecg
              'electrodermal', 'electrodermal activity', 'eda', 'scr']

LIMITS = {'cardiac': 60, 
          'ecg': 60, 
          'ekg': 60,
          'electrodermal': 1, 
          'electrodermal activity': 1, 
          'eda': 1, 
          'scr': 1}

if len(dset.description["DatasetDOI"]) > 0:
    dataset_description["SourceDatasets"] = [
            {
                "DOI": dset.description['DatasetDOI'],
            }
        ]
else:
    pass

with open(os.path.join(deriv_dir, 'dataset_description.json'), 'w') as fp:
    json.dump(dataset_description, fp)

physio_jsons = dset.get(suffix='physio', extension='json')
slices = args.slices

# check json for cardiac column
# extract sampling frequency if cardiac is present
# and then load the physio.tsv (maybe as a separate step)
physio_files = []
for file in physio_jsons:
    path = file.path
    with open(path) as json_file:
        data = json.load(json_file)
    # need to replace this if statement with something that covers all ecg and eda options
    check_cols = any(item in CANDIDATES for item in data['Columns'])
    if check_cols:
        intersection = [value for value in CANDIDATES if value in data['Columns']]
        physio_file = file.get_associations()
        assert len(physio_file) == 1, f"Found {len(physio_file)} physio files instead of the expected 1."
        data_file = physio_file[0].path
        fs = physio_file[0].get_metadata()['SamplingFrequency']

        

        print(data_file)
        # save entities as variables to pull associated BOLD meta-data
        task = file.entities['task']
        subj = file.entities['subject']
        #dtype = file.entities['datatype']

        if 'session' in file.entities.keys():
            sesh = file.entities['session']
        else:
            sesh = None
        if 'run' in file.entities.keys():
            run = file.entities['run']
        else:
            run = None
                               
        physio_dict = file.get_dict()

        out_path = os.path.join(deriv_dir, 
                                f'sub-{subj}')
        if sesh:
            out_path = os.path.join(out_path, 
                                    f'ses-{sesh}',
                                    'func')
            if not os.path.exists(out_path):
                os.makedirs(out_path)
            if run:
                basename = f'sub-{subj}_ses-{sesh}_task-{task}_run-{run}_'
            else:
                basename = f'sub-{subj}_ses-{sesh}_task-{task}_'
        else:
            out_path = os.path.join(deriv_dir, 
                                    f'sub-{subj}', 
                                    'func')
            if not os.path.exists(out_path):    
                os.makedirs(out_path)
            if run:
                basename =  f'sub-{subj}_task-{task}_run-{run}_'
            else:
                basename =  f'sub-{subj}_task-{task}_'
        
        sources = [path]
        
        notches = {}
        # params needed for physio denoising
        if tr is not None and slices is not None:
            pass
        else:
            bold_json = dset.get(task=task, 
                                    subject=subj, 
                                    session=sesh,
                                    #datatype=dtype, 
                                    run=run, 
                                    suffix='bold',
                                    extension='json')
            assert len(bold_json) == 1, f"Looking for one associated BOLD file, found {len(bold_json)}.\nEither specify --slices, --tr, --mb, --multiecho or verify that associated BOLD json sidecar exists."
            bold_json = bold_json[0]
            bold_dict = bold_json.get_dict()
            sources.append(bold_json.path)

            tr = bold_dict['RepetitionTime']
            if len(bold_dict['SliceTiming']) < 1:
                if args.slices is None:
                    raise ValueError("SliceTiming not found in BOLD sidecar, please provide using --slices")
            else:
                slices = len(bold_dict['SliceTiming'])
            if args.biopac == True: # ignore data's mb factor and filter like biopac
                mb = 1
            elif 'MultibandAccelerationFactor' in bold_dict.keys():
                mb = bold_dict['MultibandAccelerationFactor']
            elif args.mb:
                mb = args.mb
            else: # mb factor not specified, assume single band
                mb = 1

        fs = physio_dict['SamplingFrequency']

        # this is going to have to approzimate the comb 
        # so that fs is divisible by f0 ( = slices / mb / tr)
        #print(f'tr: {tr}\nmb: {mb}\nslices: {slices}\nfs: {fs}')
        
        # load the data
        dat = pd.read_table(data_file,)
        dat.columns = physio_dict['Columns']
        dat['seconds'] = dat.index / fs
        # and then run the denoising filters 
        notches['slices'] = slices / mb / tr
        
        for column in intersection:
            # compute power spectrum of raw signal
            fft_ecg, _, freq, flimit = fourier_freq(dat[column], 1/fs, LIMITS[column])
            # first, plot the raw signal and its power spectrum
            # how many samples in six seconds?
            samples = int(6 * fs)
            downsample = 10
            fig = plot_signal_fourier(time=dat['seconds'], 
                        data=dat[column], 
                        downsample=downsample, 
                        samples=samples, 
                        fft=fft_ecg, 
                        freq=freq, 
                        lim_fmax=flimit, 
                        annotate=False,
                        peaks=None,
                        slice_peaks=None,
                        title=f'Raw {column}', 
                        save=True)
            fig.savefig(os.path.join(out_path, f'{basename}desc-raw{column.capitalize()}_physio.png'), dpi=400, bbox_inches='tight')
            plt.clf()
            # let the filtering begin
            #print(notches)
            if args.multicomb:
                combs = {'biopac': {'slices': slices / tr},
                         'bottenhorn': {'slices': slices / mb / tr, 
                                        'tr': 1 / tr}}
                if column in ['eda', 'scr', 'electrodermal']:
                    Qs = {'slices': 10, 'tr': 10}
                else:
                    Qs = {'slices': 10, 'tr': 100}
                for comb in combs.keys():
                    notches = combs[comb]
                    print(column, notches)
                    filtered = dat[column]
                    for notch in notches:
                        print(notch, notches[notch],'\n',len(filtered))
                        filtered = comb_band_stop(notches[notch], filtered, Qs[notch], fs)
                        print(len(filtered))
                    dat[f'{column}_{comb}-filtered'] = filtered
                
                    # fourier transform filtered data
                    fft_ecg, _, freq, flimit = fourier_freq(filtered, 1/fs, 60)
                    fig = plot_signal_fourier(time=dat['seconds'],
                                data=filtered, 
                                downsample=downsample, 
                                samples=samples, 
                                fft=fft_ecg, 
                                freq=freq, 
                                lim_fmax=flimit, 
                                annotate=False,
                                peaks=None,
                                slice_peaks=None,
                                title=f'Filtered {column}', 
                                save=True)
                    fig.savefig(os.path.join(out_path, f'{basename}desc-{comb}Filtered{column.capitalize()}_physio.png'), 
                                dpi=400, 
                                bbox_inches='tight')
                    plt.clf()
            else:
                Q = 10
                filtered = dat[column]
                for notch in notches:
                        filtered = comb_band_stop(notches[notch], filtered, Q, fs)
                dat[f'{column}_filtered'] = filtered
                
                # fourier transform filtered data
                fft_ecg, _, freq, flimit = fourier_freq(filtered, 1/fs, 60)
                fig = plot_signal_fourier(time=dat['seconds'],
                            data=filtered, 
                            downsample=downsample, 
                            samples=samples, 
                            fft=fft_ecg, 
                            freq=freq, 
                            lim_fmax=flimit, 
                            annotate=False,
                            peaks=None,
                            slice_peaks=None,
                            title=f'Filtered {column}', 
                            save=True)
                fig.savefig(os.path.join(out_path, f'{basename}desc-filtered{column.capitalize()}_physio.png'), dpi=400, bbox_inches='tight')
                plt.clf()
        
        out_json = {
            'Description': 'raw and filtered electrophysiological measures (i.e., cardiac and skin conductance).',
            'Sources': sources,
            'RawSources': data_file,
            'SamplingFrequency': fs,
            'Columns': list(dat.columns),
        }

        if args.multicomb:
            out_json['Note'] = f'''Infinite Impulse Response comb notch filters were applied to raw data to generate `_filtered` data, 
                        with the following notch frequencies: {combs}. Slices {slices} / MB factor {mb} / TR {tr}; 1 / TR {tr}'''
        else:
            out_json['Note'] = f'''Infinite Impulse Response comb notch filters were applied to raw data to generate `_filtered` data, 
                        with the following notch frequencies: {notches}. Slices {slices} / MB factor {mb} / TR {tr}'''

        out_json_path = os.path.join(out_path, f'{basename}desc-filtered_physio.json')
        with open(out_json_path, 'w') as fp:
            json.dump(out_json, fp)
        dat.to_csv(os.path.join(out_path, f'{basename}desc-filtered_physio.tsv'), sep='\t')
        del dat
        gc.collect()
    else:
        pass