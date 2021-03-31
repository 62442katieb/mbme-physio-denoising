# Removing MR-related noise from concurrently collected electrophysiological recordings
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

`requirements.txt` specifies Python packages (and versions thereof) used in the creation and application of these filters.

### Examples
Contains interactive Jupyter notebooks (`denoising_ecg.ipynb` and `denoising_eda.ipynb`) that apply the filters to electrocardiogram and electrodermal
activity data, respectively, collected during multiband, single-echo and multi-echo EPI sequences, plots raw and cleaned data, 
and compares ferquency spectra across the filtering strategies using magnitude squared coherence.

### Figures
Contains output of Jupyter notebooks: png files of signals, their frequency spectra, and magnitude squared coherence plots.

### Example Data
Contains electrocardiogram and electrodermal activity recordings collected during one run of each a single-echo and multi-echo EPI sequence on a Siemens PRISMA 3T MRI scanner.
