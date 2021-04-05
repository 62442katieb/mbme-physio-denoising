# Removing MR-related noise from concurrently collected electrophysiological recordings

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/62442katieb/mbme-physio-denoising/binder-live?filepath=notebooks%2Fdenoising_ecg.ipynb)


Collecting physiological data during fMRI experiments can improve fMRI data cleaning and 
contribute to our understanding of psychophysiological processes; however, these recordings are 
frequently fraught with artifacts from the MRI pulse sequence. Here, we look at manufacturer 
recommendations for filtering such artifacts from physiological data collected during 
single-band, single-echo fMRI sequences and extend these recommendations to address artifacts 
associated with multiband, multi-echo fMRI sequences. While the magnitude and frequencies of 
artifacts differ with these aspects of pulse sequences, their effects can be mitigated via 
application of digital filters focused on the slice collection and repetition time.

## Directory
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

### Figures
Contains output of Jupyter notebooks: png files of signals, their frequency spectra, and magnitude squared coherence plots.

### Example Data
Contains electrocardiogram and electrodermal activity recordings collected during one run of each a single-echo and multi-echo EPI sequence on a Siemens PRISMA 3T MRI scanner.
