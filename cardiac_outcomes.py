import bids
import argparse

import numpy as np
import pandas as pd
import heartpy as hp
import seaborn as sns
import matplotlib.pyplot as plt

from os.path import join, exists
#from scipy.stats import kurtosis, wilcoxon, ttest_rel
sns.set_context('talk')

# these are the command line imports

#bids_dir = '/Users/katherine.b/Dropbox/Data/ds001242'
parser = argparse.ArgumentParser(description='Accept BIDS directory, specify # slices if slice timing isn\' specified in BOLD sidecar.')
parser.add_argument('dset', type=str,
                    help='Valid BIDS dataset containing physiological data with PhysioComb derivatives (i.e., run physioComb.py first).')
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
measures = ['kurtosis', 'good_peaks', 'bpm_mean', 'bpm_sdev', 'snr']
cols = ['']
index = pd.MultiIndex.from_product([subjects, tasks, runs])
columns = pd.MultiIndex.from_product([measures, cols])
ktdf = pd.DataFrame(index=index, columns=columns)

# loop through filtered physio data and compute relevant outcomes
for file in files:
    filename = file.filename
    subject = file.entities['subject']
    task = file.entities['task']
    fs = file.get_metadata()['SamplingFrequency']
    try:
        run = file.entities['run']
    except:
        run = 1
    df = pd.read_table(file, index_col=0)
    temp = df.filter(regex='cardiac.*', axis=1).kurtosis()
    for col in temp.keys():
        ktdf.loc[(subject, task, run), ('kurtosis', col)] = temp[col]
        ktdf.loc[(subject, task, run), ('snr', col)] = np.var(np.abs(df[col])) / np.var(df[col])
        try:
            working_data, measures = hp.process(df[col], fs)
            ktdf.at[(subject, task, run), ('good_peaks', col)] = np.mean(working_data['binary_peaklist'])
            ktdf.at[(subject, task, run), ('bpm_mean', col)] = measures['bpm']
            ktdf.at[(subject, task, run), ('bpm_sdev', col)] = working_data['hr'].std()
        except Exception as e:
            #print(e)
            ktdf.at[(subject, task, run), ('good_peaks', col)] = np.nan
            ktdf.at[(subject, task, run), ('bpm_mean', col)] = np.nan
            ktdf.at[(subject, task, run), ('bpm_sdev', col)] = np.nan

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
good_peak_long = ktdf['good_peaks'].melt(value_vars=cardiac_cols,
                           value_name='good peaks', 
                           var_name='data')
bpm_long = ktdf['bpm_mean'].melt(value_vars=cardiac_cols, 
                     value_name='bpm', 
                     var_name='data')
kurt_long = ktdf['kurtosis'].melt(value_vars=cardiac_cols, 
                     value_name='kurtosis', 
                     var_name='data')
snr_long = ktdf['snr'].melt(value_vars=cardiac_cols, 
                     value_name='snr', 
                     var_name='data')

fig,ax = plt.subplots(figsize=(10,7))
ax.set_xlim(good_peak_long['good peaks'].min(), good_peak_long['good peaks'].max() * 1.1)
sns.kdeplot(x='good peaks', 
            data=good_peak_long, 
            hue='data', 
            fill=True, 
            alpha=0.5, 
            bw_adjust=0.5,
            palette='cubehelix')
fig.savefig(join(bids_dir, 'derivatives', 'PhysioComb', 'good_peaks.png'), dpi=400, bbox_inches='tight')

fig,ax = plt.subplots(figsize=(10,7))
ax.set_xlim(left=0, right=600)
g = sns.histplot(x='bpm', data=bpm_long, hue='data', fill=True, alpha=0.5, palette='cubehelix')
fig.savefig(join(bids_dir, 'derivatives', 'PhysioComb', 'bpm.png'), dpi=400, bbox_inches='tight')

fig,ax = plt.subplots(figsize=(10,7))
ax.set_xlim(left=kurt_long['kurtosis'].min(), right=100)
g = sns.kdeplot(x='kurtosis', 
                data=kurt_long, 
                hue='data', 
                fill=True, 
                alpha=0.5,
                bw_adjust=0.5,
                palette='cubehelix')
fig.savefig(join(bids_dir, 'derivatives', 'PhysioComb', 'kurtosis.png'), dpi=400, bbox_inches='tight')

fig,ax = plt.subplots(figsize=(10,7))
ax.set_xlim(left=snr_long['snr'].min(), right=100)
g = sns.kdeplot(x='snr', 
                data=snr_long, 
                hue='data', 
                fill=True, 
                alpha=0.5,
                bw_adjust=0.5,
                palette='cubehelix')
fig.savefig(join(bids_dir, 'derivatives', 'PhysioComb', 'snr.png'), dpi=400, bbox_inches='tight')