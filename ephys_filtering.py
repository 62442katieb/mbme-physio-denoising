import argparse
import os
import bioread as br
import scipy.signal as signal
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

pd.options.mode.chained_assignment = None

def consecutive(data, stepsize=0.000501):
    return np.split(data, np.where(np.diff(data) != stepsize)[0]+1)

def comb_band_stop(notches, filtered, Q, fs):
    for notch in notches:
        print('Cleaning', notch, '@', notches[notch])
        for i in np.arange(1, (nyquist / notches[notch])):
            #print(notches[notch] * i)
            f0 = notches[notch] * i
            w0 = f0/nyquist
            b,a = signal.iirnotch(w0, Q)
            filtered = signal.filtfilt(b, a, filtered)
    return filtered

parser = argparse.ArgumentParser(description='Accept physio + sidecar input for cleaning EDA/ECG collected simultaneously with fMRI.')
parser.add_argument('in_file', type=str, 
                    help='AcqKnowledge file containing physio measurements from a single scan session.')
parser.add_argument('out_dir', type=str, help='Absolute or relative path to directory where output (figures and cleaned data) will be saved.')
parser.add_argument('tr', type=float, 
                    help='The TR (repetition time) in seconds of the fMRI scan sequence used during colleciton of these data.')
parser.add_argument('--mb', type=int, 
                    help='Multiband factor of fMRI scan sequence (if single band, --mb=1).')
parser.add_argument('--verbose', action='store_true', 
                    help='Print information as the cleaning script runs.')
parser.add_argument('slices', type=int, 
                    help='Number of slices collected by fMRI scan sequence.')
args = parser.parse_args()


data_fname = args.in_file
ext = data_fname.split('.')[-1]
if 'acq' in ext:
    data = br.read_file(data_fname)
if 'csv' in ext:
    data = pd.read_csv(data_fname)
if 'tsv' in ext:
    data = pd.read_csv(date_fname, sep='\t')
basename = data_fname.split('/')[-1][:-4]
out_dir = args.out_dir

if not os.path.exists('{0}/data'.format(out_dir)):
    os.mkdir('{0}/data'.format(out_dir))
    os.mkdir('{0}/data/clean'.format(out_dir))
    os.mkdir('{0}/data/raw'.format(out_dir))
    os.mkdir('{0}/figures'.format(out_dir))

if args.verbose:
    print(basename)

#arranging noisy frequencies to be filtered out
slices = args.slices
tr = args.tr

if args.mb:
    mb = args.mb
else:
    mb = 1

#print the scan parameters:
if args.verbose:
    print('======== fMRI sequence parameters ========\nTR = {0} s\n# Slices = {1}\nMB Factor = {2}'.format(tr, slices, mb))
    print('==========================================\n')

notches = {'slices': slices / mb / tr, 
           'tr': 1 / tr}      #turns out those evenly-spaced interfering frequencies are 0.66Hz apart

fs = data.samples_per_second


#I don't know if this is right, been playing around with the volume of Q
Q = 100

nyquist = fs/2

for channel in data.named_channels:
    #print(channel)
    if 'ECG' in channel:
        ecg_channel = channel
        if args.verbose:
            print('ECG channel:', channel)
    elif 'Digital' in channel:
        trigger = channel
        print('Trigger channel:', channel)
    elif 'EDA' in channel:
        eda_channel = channel
        if args.verbose:
            print('EDA channel:', channel)
    elif 'Respiration' in channel:
        resp_channel = channel
        if args.verbose:
            print('Resp channel:', channel)

timeseries = pd.DataFrame(columns=['ECG', 'EDA', 'Trigger', 'Resp', 'seconds'])

timeseries['Trigger'] = data.named_channels[trigger].data
timeseries['ECG'] = data.named_channels[ecg_channel].data
timeseries['EDA'] = data.named_channels[eda_channel].data
timeseries['Resp'] = data.named_channels[resp_channel].data
timeseries['seconds'] = data.time_index
timeseries.to_csv('{0}/data/raw/{1}-raw.csv'.format(out_dir, basename))

data = None

fives = timeseries[timeseries['Trigger'] == 5].index.values
scan_idx = consecutive(fives, stepsize=1)

#separating timeseries collected during BOLD scan
#(where trigger channel = 5V)
fft_ecg = np.fft.fft(timeseries['ECG'][:1000000])
freq = np.fft.fftfreq(timeseries['ECG'][:1000000].shape[-1], d=0.0005)
fft_ecg_db = 10 * np.log10(abs(fft_ecg))

if args.verbose:
    print('Plotting ECG example...')
gridkw = dict(width_ratios=[2,1])
fig, (ax1, ax2) = plt.subplots(1, 2, gridspec_kw=gridkw, figsize=(20,4))
plt.tight_layout(pad=2)

sns.lineplot(signal.decimate(timeseries['seconds'][50000:80000],10)[:8000], 
             signal.decimate(timeseries['ECG'][50000:80000], 10)[:8000], 
             linewidth=1, ax=ax1) #array, top subplot
sns.lineplot(freq[:30000], 
             fft_ecg.real[:30000], 
             ax=ax2, linewidth=1)#bottom subplot
ax1.set_xlabel('seconds')
ax1.set_ylabel('mV')
#ax3.set_xlabel('Hz')
ax2.set_xlabel('Hz')
ax2.set_yticklabels([0,-10, -5, 0, 5])
ax2.set_ylabel('Power (x10,000)')
fig.savefig('{0}/figures/{1}-ecg_raw-noepi.png'.format(out_dir, basename), dpi=300)


fft_eda = np.fft.fft(timeseries['EDA'][:1000000])
freq = np.fft.fftfreq(timeseries['EDA'][:1000000].shape[-1], d=0.0005)
fft_eda_db = 10 * np.log10(abs(fft_ecg))

if args.verbose:
    print('Plotting EDA example...')
gridkw = dict(width_ratios=[2,1])
fig, (ax1, ax2) = plt.subplots(1, 2, gridspec_kw=gridkw, figsize=(20,4))
plt.tight_layout(pad=2)

sns.lineplot(signal.decimate(timeseries['seconds'][50000:80000],10)[:80000], 
             signal.decimate(timeseries['EDA'][50000:80000], 10)[:80000], 
             linewidth=1, ax=ax1) #array, top subplot
sns.lineplot(freq[:80000], 
             fft_eda_db[:80000], 
             ax=ax2, linewidth=1)#bottom subplot
ax1.set_xlabel('seconds')
ax1.set_ylabel('mV')
#ax3.set_xlabel('Hz')
ax2.set_xlabel('Hz')
ax2.set_yticklabels([0,-10, -5, 0, 5])
ax2.set_ylabel('Power (x10,000)')
fig.savefig('{0}/figures/{1}-eda_raw-noepi.png'.format(out_dir, basename), dpi=300)


i = 1
for scan in scan_idx:
    if args.verbose:
        print('Now cleaning scan #', i)
    scan1 = timeseries[timeseries.index.isin(scan)]

    fft_ecg = np.fft.fft(scan1['ECG'].values)
    freq = np.fft.fftfreq(scan1['ECG'].values.shape[-1], d=0.0005)
    fft_ecg_db = 10 * np.log10(abs(fft_ecg))
    gridkw = dict(width_ratios=[2,1])
    fig, (ax1, ax2) = plt.subplots(1, 2, gridspec_kw=gridkw, figsize=(20,4))
    plt.tight_layout(pad=2)

    sns.lineplot(signal.decimate(scan1['seconds'][:954000], 10)[:8000], 
                 signal.decimate(scan1['ECG'][:954000], 10)[:8000], 
                 linewidth=1, ax=ax1) #array, top subplot
    #sns.lineplot(freq[:100000], fft_ecg_db[:100000], linewidth=1, ax=ax2)
    sns.lineplot(freq[:30000], 
                 fft_ecg.real[:30000], 
                 ax=ax2, linewidth=1)#bottom subplot
    ax1.set_xlabel('seconds')
    ax1.set_ylabel('mV')
    ax2.set_xlabel('Hz')
    ax2.set_yticklabels([0,-10, -5, 0, 5])
    ax2.set_ylabel('Power (x10,000)')
    
    fig.savefig('{0}/figures/{1}_scan-{2}_ecg-raw.png'.format(out_dir, basename, i), dpi=300)
    plt.close()

    filtered = comb_band_stop(notches, scan1['ECG'], Q, fs)

    fft_filt = np.fft.fft(filtered)
    freq = np.fft.fftfreq(filtered.shape[-1], d=0.0005)
    fft_ecg_db = 10 * np.log10(abs(filtered))

    gridkw = dict(width_ratios=[2,1])
    fig, (ax1, ax2) = plt.subplots(1, 2, gridspec_kw=gridkw, figsize=(20,4))

    plt.tight_layout(pad=2)

    sns.lineplot(signal.decimate(scan1['seconds'][:954000], 10)[:8000], 
                 signal.decimate(filtered, 10)[:8000], 
                 linewidth=1, ax=ax1) #array, top subplot
    #sns.lineplot(freq[:100000], fft_ecg_db[:100000], linewidth=1, ax=ax2)
    sns.lineplot(freq[:30000], 
                 fft_filt.real[:30000], 
                 ax=ax2, linewidth=1)#bottom subplot
    ax1.set_xlabel('seconds')
    ax1.set_ylabel('mV')
    #ax3.set_xlabel('Hz')
    ax2.set_xlabel('Hz')
    ax2.set_yticklabels([0,-2.5, -5, 0, 2.5, 5])
    ax2.set_ylabel('Power (x10,000)')
    fig.savefig('{0}/figures/{1}-scan{2}-ecg_clean.png'.format(out_dir, basename, i), dpi=300)
    plt.close()

    scan1.at[scan, 'Clean ECG'] = filtered

    fft_eda = np.fft.fft(scan1['EDA'].values)
    freq = np.fft.fftfreq(scan1['EDA'].values.shape[-1], d=0.0005)
    fft_eda_db = 10 * np.log10(abs(fft_eda))
    
    #plotting raw data and frequency spectrum 
    gridkw = dict(width_ratios=[2,1])
    fig, (ax1, ax2) = plt.subplots(1, 2, gridspec_kw=gridkw, figsize=(20,4))
    plt.tight_layout(pad=2)

    sns.lineplot(signal.decimate(scan1['seconds'], 10), 
                 signal.decimate(scan1['EDA'], 10), 
                 linewidth=1, ax=ax1) #array, top subplot
    #sns.lineplot(freq[:100000], fft_ecg_db[:100000], linewidth=1, ax=ax2)
    sns.lineplot(freq[:80000], 
                 fft_eda_db[:80000], 
                 ax=ax2, linewidth=1)#bottom subplot
    ax1.set_xlabel('Seconds')
    ax1.set_ylabel('microsiemens')
    ax2.set_xlabel('Hz')
    ax2.set_ylabel('dB')

    fig.savefig('{0}/figures/{1}-scan{2}-eda_raw.png'.format(out_dir, basename, i), dpi=300)
    plt.close()
    
    filtered = comb_band_stop(notches, scan1['EDA'], 100, fs)

    fft_filt = np.fft.fft(filtered)
    freq = np.fft.fftfreq(filtered.shape[-1], d=0.0005)
    fft_eda_db = 10 * np.log10(abs(filtered))

    gridkw = dict(width_ratios=[2,1])
    fig, (ax1, ax2) = plt.subplots(1, 2, gridspec_kw=gridkw, figsize=(20,4))

    plt.tight_layout(pad=2)

    sns.lineplot(signal.decimate(scan1['seconds'], 10), 
                 signal.decimate(filtered, 10), 
                 linewidth=1, ax=ax1) #array, top subplot
    #sns.lineplot(freq[:100000], fft_ecg_db[:100000], linewidth=1, ax=ax2)
    sns.lineplot(freq[:80000], 
                 fft_eda_db[:80000], 
                 ax=ax2, linewidth=1)#bottom subplot
    ax1.set_xlabel('Seconds')
    ax1.set_ylabel('microsiemens')
    ax2.set_xlabel('Hz')
    ax2.set_ylabel('dB')
    #ax2.set_yticklabels([0,-2.5, -5, 0, 2.5, 5])
    fig.savefig('{0}/figures/{1}-scan{2}-eda_clean.png'.format(out_dir, basename, i), dpi=300)
    plt.close()
    scan1.at[scan, 'Clean EDA'] = filtered
    scan1.drop(['Trigger'], axis=1, inplace=True)
    scan1.to_csv('{0}/data/clean/{1}-scan{2}-clean.csv'.format(out_dir, basename, i))
    i += 1
