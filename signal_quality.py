import bids
import argparse

import numpy as np
import pandas as pd
import neurokit2 as nk
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.signal import welch
from os.path import join, exists
#from scipy.stats import kurtosis, wilcoxon, ttest_rel
sns.set_context('talk')

# these are the command line imports

#bids_dir = '/Users/katherine.b/Dropbox/Data/ds001242'
parser = argparse.ArgumentParser(description='Accept BIDS directory, specify # slices if slice timing isn\' specified in BOLD sidecar.')
parser.add_argument('dset', type=str,
                    help='Valid BIDS dataset containing physiological data with PhysioComb derivatives (i.e., run physioComb.py first).')
parser.add_argument('--verbose', action='store_true', 
                    help='Print filename as the script runs.')
args = parser.parse_args()
bids_dir = args.dset

if not exists(join(bids_dir, 'derivatives', 'PhysioComb')):
    raise FileNotFoundError("Cannot find PhysioComb derivaties, please run physioComb.py.")

layout = bids.BIDSLayout(bids_dir, derivatives=True)
files = layout.derivatives['PhysioComb'].get(#return_type='filename', 
                                     extension='tsv', 
                                     desc='filtered',
                                     invalid_filters='allow')

# read exp design from bids layout
runs = layout.get_runs()
tasks = layout.get_tasks()
subjects = layout.get_subjects()
sessions = layout.get_session()

# set up multiindex for comparing measures across filtered and unfiltered data
measures = ['kurtosis', 'zhao_quality', 'bpm_mean', 'bpm_sdev', 'snr', 'noise']
cols = ['']
if len(sessions) > 1:
    index = pd.MultiIndex.from_product([subjects, sessions, tasks, runs])
else:
    sessions = [1]
    index = pd.MultiIndex.from_product([subjects, sessions, tasks, runs])
columns = pd.MultiIndex.from_product([measures, cols])
ktdf = pd.DataFrame(index=index, columns=columns)

# loop through filtered physio data and compute relevant outcomes
for file in files:
    filename = file.filename
    if args.verbose:
        print(filename)
    subject = file.entities['subject']
    task = file.entities['task']
    fs = file.get_metadata()['SamplingFrequency']
    try:
        run = file.entities['run']
    except:
        run = 1
    try: 
        session = file.entities['session']
    except:
        session = 1
    df = pd.read_table(file.path)
    temp = df.filter(regex='cardiac.*', axis=1).kurtosis()
    for col in temp.keys():
        ktdf.loc[(subject, session, task, run), ('kurtosis', col)] = temp[col]
        ktdf.loc[(subject, session, task, run), ('snr', col)] = np.var(np.abs(df[col])) / np.var(df[col])
        try:
            rpeaks, info = nk.ecg_peaks(df[col], sampling_rate=fs)
            # Compute rate
            ecg_rate = nk.ecg_rate(rpeaks, sampling_rate=fs, desired_length=len(df[col]))
            q = nk.ecg_quality(df[col],
                                sampling_rate=fs,
                                method="zhao2018",
                                approach="fuzzy")
            
            ktdf.at[(subject, session, task, run), ('zhao_quality', col)] = q
            ktdf.at[(subject, session, task, run), ('bpm_mean', col)] = np.mean(ecg_rate)
            ktdf.at[(subject, session, task, run), ('bpm_sdev', col)] = np.std(ecg_rate)
        except Exception as e:
            #print(e)
            ktdf.at[(subject, session, task, run), ('zhao_quality', col)] = np.nan
            ktdf.at[(subject, session, task, run), ('bpm_mean', col)] = np.nan
            ktdf.at[(subject, session, task, run), ('bpm_sdev', col)] = np.nan
    # and now do EDA
    temp = df.filter(regex='scr.*', axis=1).kurtosis()
    for col in temp.keys():
        # compute noise power
        f, Pxx = welch(df[col], fs, nperseg=8192)
        noise_pct = np.trapz(Pxx[f > 0.4], f[f > 0.4]) / (np.trapz(Pxx, f))
        ktdf.loc[(subject, session, task, run), ('noise', col)] = noise_pct

# not all tasks will have all runs, sessions, etc.
# remove blank col created for purpose of multiindex
ktdf.dropna(how='all', inplace=True, axis=1)
ktdf.dropna(how='all', inplace=True, axis=0)

# save outcomes
ktdf.to_csv(join(bids_dir, 'derivatives', 'PhysioComb', 'cardiac_outcomes.tsv'), sep='\t')

filenames = layout.derivatives['PhysioComb'].get(return_type='filename', 
                                     extension='tsv', 
                                     desc='filtered',
                                     invalid_filters='allow')
filenames = layout.get(return_type='filename', 
                                     extension='tsv', 
                                     suffix='physio',
                                     invalid_filters='allow')
# write a corresponding sidecar
ktjson = {
    'Description': '''Measures computed from cardiac data, both filtered and unfiltered, 
                      for use in future analysis and to diagnose filtration efficacy.''',
    'Sources': filenames,
    }

cardiac_cols = list(ktdf.columns.get_level_values(1).unique())

# plots and comparisons
good_peak_long = ktdf['zhao_quality'].melt(value_vars=ktdf['zhao_quality'].columns,
                           value_name='quality', 
                           var_name='data')
bpm_long = ktdf['bpm_mean'].melt(value_vars=ktdf['bpm_mean'].columns, 
                     value_name='bpm', 
                     var_name='data')
kurt_long = ktdf['kurtosis'].melt(value_vars=ktdf['kurtosis'].columns, 
                     value_name='kurtosis', 
                     var_name='data')
snr_long = ktdf['snr'].melt(value_vars=ktdf['snr'].columns, 
                     value_name='snr', 
                     var_name='data')

noise_long = ktdf['noise'].melt(value_vars=ktdf['noise'].columns, 
                     value_name='noise', 
                     var_name='data')


# Use this to position the legend for this plot
loc = {'Unnacceptable': 'upper left',
       'Barely acceptable': 'upper center',
       'Excellent': 'upper right'}
# Find quality level with least # of runs
x = {}
for qual in ['Unnacceptable', 'Barely acceptable', 'Excellent']:
    x[qual] = len(good_peak_long[good_peak_long['quality'] == qual].index)
least = list(dict(sorted(x.items(), key=lambda item: item[1])).keys())[0]

fig,ax = plt.subplots(figsize=(10,7))
sns.countplot(x='quality', 
            data=good_peak_long, 
            hue='data', 
            order=['Unnacceptable', 'Barely acceptable', 'Excellent'],
            palette='cubehelix')
ax.legend(loc=loc[least])
fig.savefig(join(bids_dir, 'derivatives', 'PhysioComb', 'zhao_quality.png'), dpi=400, bbox_inches='tight')

fig,ax = plt.subplots(figsize=(10,7))
ax.set_xlim(left=bpm_long['bpm'].min() * 0.9, right=bpm_long['bpm'].max() * 1.1)
g = sns.histplot(x='bpm', data=bpm_long, hue='data', fill=True, alpha=0.5, palette='cubehelix')
fig.savefig(join(bids_dir, 'derivatives', 'PhysioComb', 'bpm.png'), dpi=400, bbox_inches='tight')

fig,ax = plt.subplots(figsize=(10,7))
ax.set_xlim(left=kurt_long['kurtosis'].min() * 0.9, right=kurt_long['kurtosis'].max() * 1.1)
g = sns.kdeplot(x='kurtosis', 
                data=kurt_long, 
                hue='data', 
                fill=True, 
                alpha=0.5,
                bw_adjust=0.5,
                palette='cubehelix')
fig.savefig(join(bids_dir, 'derivatives', 'PhysioComb', 'kurtosis.png'), dpi=400, bbox_inches='tight')

fig,ax = plt.subplots(figsize=(10,7))
ax.set_xlim(left=snr_long['snr'].min() * 0.9, right=snr_long['snr'].max() * 1.1)
g = sns.kdeplot(x='snr', 
                data=snr_long, 
                hue='data', 
                fill=True, 
                alpha=0.5,
                bw_adjust=0.5,
                palette='cubehelix')
fig.savefig(join(bids_dir, 'derivatives', 'PhysioComb', 'snr.png'), dpi=400, bbox_inches='tight')

fig,ax = plt.subplots(figsize=(10,7))
ax.set_xlim(0,1)
g = sns.kdeplot(x='noise', 
                data=noise_long, 
                hue='data', 
                fill=True, 
                alpha=0.5,
                bw_adjust=0.5,
                palette='cubehelix')
fig.savefig(join(bids_dir, 'derivatives', 'PhysioComb', 'noise.png'), dpi=400, bbox_inches='tight')