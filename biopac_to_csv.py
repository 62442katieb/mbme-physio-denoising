import os
import numpy as np
import pandas as pd
import bioread as br
import argparse

def consecutive(data, stepsize=0.000501):
    return np.split(data, np.where(np.diff(data) != stepsize)[0]+1)

parser = argparse.ArgumentParser(description='Accept BIOPAC .acq files as input, save them as .csv files.')
parser.add_argument('in_files', nargs='+', type=str, 
                    help='AcqKnowledge file(s) containing physio measurements from a single scan session.')
parser.add_argument('out_dir', type=str, help='Absolute or relative path to directory where output (figures and cleaned data) will be saved.')
args = parser.parse_args()

if args.out_dir:
    out_dir = args.out_dir
else:
    out_dir = '.'

files = args.in_files
for file in files:
    print(file)
    data_fname = file
    data = br.read_file(data_fname)
    basename = data_fname.split('/')[-1][:-4]
    

    if not os.path.exists('{0}/data'.format(out_dir)):
        os.mkdir('{0}/data'.format(out_dir))
        os.mkdir('{0}/data/clean'.format(out_dir))
        os.mkdir('{0}/data/raw'.format(out_dir))
        os.mkdir('{0}/figures'.format(out_dir))

    for channel in data.named_channels:
        #print(channel)
        if 'ECG' in channel:
            print('ECG channel:', channel)
            ecg_channel = channel
        elif 'Trigger' in channel:
            print('Trigger channel:', channel)
            trigger = channel
        elif 'Digital' in channel:
            print('Trigger channel:', channel)
            trigger = channel
        elif 'EDA' in channel:
            print('EDA channel:', channel)
            eda_channel = channel
        elif 'Respiration' in channel:
            print('Resp channel:', channel)
            resp_channel = channel

    timeseries = pd.DataFrame(columns=['ECG', 'EDA', 'Trigger', 'Resp', 'seconds'])

    timeseries['Trigger'] = data.named_channels[trigger].data
    timeseries['ECG'] = data.named_channels[ecg_channel].data
    timeseries['EDA'] = data.named_channels[eda_channel].data
    timeseries['Resp'] = data.named_channels[resp_channel].data
    timeseries['seconds'] = data.time_index


    #separating timeseries collected during BOLD scan
    #(where trigger channel = 5V)
    fives = timeseries[timeseries['Trigger'] == 5].index.values
    scan_idx = consecutive(fives, stepsize=1)

    scans = {}
    for i in range(len(scan_idx)):
        duration = len(scan_idx[i])/2000./60
        #print(i, np.round(duration, 1), 'minutes')
        if duration > 2.:
            scans[i] = duration
    print('Scan index and duration (in minutes):\n',scans)

    for i in range(0, len(scan_idx)):
        for j in scan_idx[i]:
            timeseries.at[j,'scan'] = i

    timeseries.to_csv('{0}/data/raw/{1}-raw.csv'.format(out_dir,basename))
    print('File saved at:{0}/data/raw/{1}-raw.csv'.format(out_dir,basename))