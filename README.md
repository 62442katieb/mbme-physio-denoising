# Removing MR-related noise from concurrently collected electrophysiological recordings

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/62442katieb/mbme-physio-denoising/binder-live?filepath=notebooks%2Fdenoising_ecg.ipynb)
## Table of Contents
1. [Directory](https://github.com/62442katieb/mbme-physio-denoising#Directory)
- [Workflow](https://github.com/62442katieb/mbme-physio-denoising#Workflow)
- [Outdated](https://github.com/62442katieb/mbme-physio-denoising#Outdated)
- [Notebooks](https://github.com/62442katieb/mbme-physio-denoising#Notebooks)
2. [Code for Manuscript](https://github.com/62442katieb/mbme-physio-denoising/edit/main/README.md#code-for-manuscript)

Collecting physiological data during fMRI experiments can improve fMRI data cleaning and 
contribute to our understanding of psychophysiological processes; however, these recordings are 
frequently fraught with artifacts from the MRI pulse sequence. Here, we look at manufacturer 
recommendations for filtering such artifacts from physiological data collected during 
single-band, single-echo fMRI sequences and extend these recommendations to address artifacts 
associated with multiband, multi-echo fMRI sequences. While the magnitude and frequencies of 
artifacts differ with these aspects of pulse sequences, their effects can be mitigated via 
application of digital filters focused on the slice collection and repetition time.

## Directory
### Workflow
The main program in this package is `physioComb.py `, which is a BIDS-friendly command-line Python script that performs digital filtering of MR-related 
artifacts from eletrophysiological data (i.e., electrocardiogram and electrodermal activity), using IIR notch filters centered at the slice collection frequency and (optionally) TR frequency. NOTE: current BIDS specifications recommend using `cardiac` to name columns in physio containing cardiac data, regardless of recording source (i.e., ECG or PPG) and `physioComb.py` just looks for a `cardiac` column. Proceed with caution and/or discard `cardiac_filtered*` columns if your cardiac data is from a photoplethysmogram.
```
usage: physioComb.py [-h] [--slices SLICES] [--mb MB] [--tr TR] [--multiecho MULTIECHO] [--biopac] [--progress] [--multicomb MULTICOMB]
                     [--verbose]
                     dset

Accept BIDS directory, specify # slices if slice timing isn't specified in BOLD sidecar.

positional arguments:
  dset                  Valid BIDS dataset containing physiological data.

optional arguments:
  -h, --help            show this help message and exit
  --slices SLICES       The number of slices acquired in BOLD sequences concurrent with physiological data acquisition. If multiple BOLD
                        sequences were used with different numbers of slices, this option is invalid, please make sure slice timing is
                        specified in the json sidecar for each BOLD sequence.
  --mb MB               Multiband factor of fMRI scan sequence (if single band, --mb=1).
  --tr TR               Repetition time of fMRI scan sequence.
  --multiecho MULTIECHO
                        Does the fMRI scan sequence collect multiple echoes.
  --biopac              Run only BIOPAC-recommended filtering at single-band slice collection frequency.
  --progress            Display progress bar (note: requires enligten).
  --multicomb MULTICOMB
                        Run multiple comb notch filtering strategies and save all outputs.
  --verbose             Print information as the cleaning script runs.
```
To compare signal quality indices across signals (i.e., filtered and raw), `signalQuality.py` reads in a BIDS dataset (with PhysioComb derivatives) and calculates kurtosis, signal-to-noise ratio, and an integrated signal quality heuristic based on work from [Zhao & Zhang (2018)](https://www.frontiersin.org/articles/10.3389/fphys.2018.00727).
```
usage: signal_quality.py [-h] [--verbose] dset

Accept BIDS directory, specify # slices if slice timing isn't specified in BOLD sidecar.

positional arguments:
  dset        Valid BIDS dataset containing physiological data with PhysioComb derivatives (i.e., run physioComb.py first).

optional arguments:
  -h, --help  show this help message and exit
  --verbose   Print filename as the script runs.
```

### Outdated
These scripts were written for a previous version of a [manuscript-in-progress](https://www.biorxiv.org/content/10.1101/2021.04.01.437293) comparing comb filters for MR-related noise removal.

`ephys_filtering.py` is a command-line Python script that performs digital filtering of MR-related artifacts from eletrophysiological data, 
using IIR notch filters centered at slice collection and TR frequencies.
```
usage: ephys_filtering.py [-h] [--mb MB] [--verbose] in_file out_dir tr slices

Accept physio + sidecar input for cleaning EDA/ECG collected simultaneously
with fMRI.

positional arguments:
  in_file     AcqKnowledge file containing physio measurements from a single
              scan session.
  out_dir     Absolute or relative path to directory where output (figures and
              cleaned data) will be saved.
  tr          The TR (repetition time) in seconds of the fMRI scan sequence
              used during colleciton of these data.
  slices      Number of slices collected by fMRI scan sequence.

optional arguments:
  -h, --help  show this help message and exit
  --mb MB     Multiband factor of fMRI scan sequence (if single band, --mb=1).
  --verbose   Print information as the cleaning script runs.
```

`biopac_to_csv.py` is a command-line Python script that converts BIOPAC AcqKnowledge files (extension `.acq`) to comma-separated value files (`.csv`) 
with a common time index and standardized headers (currently only supports ECG, EDA, and Respiration recordings).
```
usage: biopac_to_csv.py [-h] in_files [in_files ...] out_dir

Accept BIOPAC .acq files as input, save them as .csv files.

positional arguments:
  in_files    AcqKnowledge file(s) containing physio measurements from a
              single scan session.
  out_dir     Absolute or relative path to directory where output (figures and
              cleaned data) will be saved.

optional arguments:
  -h, --help  show this help message and exit
```

`requirements.txt` specifies Python packages (and versions thereof) used in the creation and application of these filters.

### Notebooks
Contains interactive Jupyter notebooks (`denoising_ecg.ipynb` and `denoising_eda.ipynb`) that apply the filters to electrocardiogram and electrodermal
activity data, respectively, collected during multiband, single-echo and multi-echo EPI sequences, plots raw and cleaned data, 
and compares ferquency spectra across the filtering strategies using magnitude squared coherence.

## Code for Manuscript
To perform the analyses detailed in [Bottenhorn et al., (preprint)](https://doi.org/10.1101/2021.04.01.437293), the following workflow was run:
1. `physioComb.py --multicomb=True` to filter ECG and EDA signals according to BIOPAC recommendations and updated filters
2. `msc_and_plot.py` to calculates magnitude squared coherence, pairwise, between signals per supported modality (i.e., ECG, EDA)
3. `signal_quality.py` to calculate signal quality indices across filtering strategies.
