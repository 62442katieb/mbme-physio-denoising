#import os
#import bioread as br
import scipy.signal as signal
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

pd.options.mode.chained_assignment = None

sns.set_context('paper', font_scale=3)

pal = sns.crayon_palette(['Gray', 'Razzmatazz', 'Tropical Rain Forest'])
sns.set_palette(pal)

def consecutive(data, stepsize=0.000501):
    return np.split(data, np.where(np.diff(data) != stepsize)[0]+1)

#def comb_band_stop(notches, filtered, Q, fs):
#    for notch in notches:
#        print('Cleaning', notch, '@', notches[notch])
#        for i in np.arange(1, (nyquist / notches[notch])):
#            #print(notches[notch] * i)
#            f0 = notches[notch] * i
#            w0 = f0/nyquist
#            b,a = signal.iirnotch(w0, Q)
#            filtered = signal.filtfilt(b, a, filtered)
#    return filtered

def butter_highpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = signal.butter(order, normal_cutoff, btype='high', analog=False)
    return b, a

def butter_highpass_filter(data, cutoff, fs, order=5):
    b, a = butter_highpass(cutoff, fs, order=order)
    y = signal.filtfilt(b, a, data)
    return y

def comb_band_stop(notches, data, Q, fs):
        nyquist = fs / 2
        for notch in notches:
            f0 = notches[notch]
            w0 = f0 / nyquist
            if fs % w0 > 0:
                # do it the old way
                filtered = data.copy()
                max_harmonic = int(nyquist / f0)
                for i in np.arange(1, max_harmonic):
                    
                    f0 = notches[notch] * i
                    w0 = f0 / nyquist
                    b,a = signal.iirnotch(w0, Q, fs=fs)
                    filtered = signal.filtfilt(b, a, filtered)
            else:
                b,a = signal.iircomb(w0, Q, ftype='notch', fs=fs)
                filtered = signal.filtfilt(b, a, data)
        return filtered

def fourier_freq(timeseries, d, fmax):
    fft = np.fft.fft(timeseries)
    freq = np.fft.fftfreq(timeseries.shape[-1], d=d)
    fft_db = 10 * np.log10(abs(fft))
    limit = np.where(freq >= fmax)[0][0]
    return fft, fft_db, freq, limit

def plot_signal_fourier(time, data, downsample, samples, fft, freq, lim_fmax, 
                        annotate=False, peaks=None, slice_peaks=None, title=None, save=True):
    gridkw = dict(width_ratios=[2,1])
    fig, (ax1, ax2) = plt.subplots(1, 2, gridspec_kw=gridkw, figsize=(20,5))
    plt.tight_layout(h_pad=4)

    time_dec = signal.decimate(time, downsample)
    data_dec = signal.decimate(data, downsample)

    limit = int(samples / downsample)

    time_domain = pd.DataFrame()
    time_domain['x'] = time_dec[0: limit]
    time_domain['y'] = data_dec[0: limit]
    freq_domain = pd.DataFrame()
    freq_domain['x'] = freq[:lim_fmax]
    freq_domain['y'] = fft.real[:lim_fmax]
    sns.lineplot(x='x', y='y', data=time_domain, 
                 linewidth=1, ax=ax1) #array, left subplot
    sns.lineplot(x='x', y='y', data=freq_domain, 
                 ax=ax2, linewidth=1)#right subplot
    if annotate:
        ax2.plot(freq[peaks][:50], fft.real[peaks][:50], "^", ms=5)
        ax2.plot(freq[slice_peaks][:4], fft.real[slice_peaks][:4], "o", ms=7)
    else:
        pass
    ax1.set_xlabel('seconds')
    ax1.set_ylabel('mV')
    #ax1.set_yticks([-1,0,1,2])
    ax2.set_xlabel('Hz')
    ax2.ticklabel_format(axis='y', style='sci', scilimits=(-1,1))
    ax2.set_ylabel('Power')
    ax1.set_title(title, pad=40)
    ax2.set_title('{0} frequencies'.format(title),pad=40)
    plt.close()
    #if save:
    #    
    #    fig.savefig('../figures/{title}.svg'.format(title=title))
    #else:
    #    pass
    return fig