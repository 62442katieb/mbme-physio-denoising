import bids

import bioread as br
import pandas as pd
import seaborn as sns

from os.path import join, exists
from filtering import plot_signal_fourier, fourier_freq

diva = '/home/data/nbc/Laird_DIVA/dset'
raw = '/home/data/nbc/Laird_DIVA/sourcedata/biopac'

layout = bids.BIDSLayout(diva, derivatives=True)

subjects = layout.get_subjects()
sessions = layout.get_sessions()

for subject in subjects:
    print(subject)
    if not 'PILOT' in subject:
        for session in sessions:
            print(session)
            biopac = join(raw, f'sub-{subject}', f'sub-{subject}_ses-{session}.acq')
            out_path = join(diva, 
                            'derivatives', 
                            'PhysioComb', 
                            f'sub-{subject}', 
                            f'ses-{session}', 
                            f'sub-{subject}_ses-{session}_desc-noMR_physio.tsv')
            if exists(biopac):
                dat = br.read_file(biopac)

                timeseries = pd.DataFrame(columns=['cardiac', 'eda', 'trigger', 'seconds'])

                timeseries['trigger'] = dat.named_channels['Trigger'].data
                timeseries['cardiac'] = dat.named_channels['ECG100C - ECG100C'].data
                timeseries['eda'] = dat.named_channels['EDA100C - EDA100C-MRI'].data
                timeseries['seconds'] = dat.time_index

                no_epi = timeseries[timeseries['trigger'] != 5]

                fft_ecg_bt, _, freq, limit = fourier_freq(no_epi['cardiac'], 0.0005, 60)
                fig = plot_signal_fourier(time=no_epi['seconds'][100000:],
                                    data=no_epi['cardiac'][100000:],
                                    downsample=10,
                                    samples=12000,
                                    fft=fft_ecg_bt,
                                    freq=freq,
                                    lim_fmax=60,
                                    annotate=False,
                                    peaks=None,
                                    slice_peaks=None,
                                    title=f'No-MR ECG - {subject}, session {session}',
                                    save=True)
                fig.savefig(f'{out_path[:-4]}_ecg.png', dpi=400, bbox_inches='tight')

                fft_eda_bt, _, freq, limit = fourier_freq(no_epi['eda'], 0.0005, 60)
                fig = plot_signal_fourier(time=no_epi['seconds'][100000:],
                                    data=no_epi['eda'][100000:],
                                    downsample=10,
                                    samples=12000,
                                    fft=fft_eda_bt,
                                    freq=freq,
                                    lim_fmax=60,
                                    annotate=False,
                                    peaks=None,
                                    slice_peaks=None,
                                    title=f'No-MR EDA - {subject}, session {session}',
                                    save=True)
                fig.savefig(f'{out_path[:-4]}_eda.png', dpi=400, bbox_inches='tight')

                
                no_epi.to_csv(out_path, sep='\t') 
            else:
                print('couldn\'t find', biopac)
            