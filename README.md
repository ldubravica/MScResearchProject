### 1. Preprocessing

**Code:** Depression_Study\phase_one_raw.ipynb

Preprocceses the raw EEG data in Python. Produces clean and filtered EEG per subject in *Depression_Study\export_mat\XXX_Depression_REST_processed.mat*.

### 2. Source Reconstruction & CSER Calculation

**Code:** EEGtoCSER.m

Performs source reconstruction and CSER calculations in MATLAB on the electrode-based preprocessed EEG. Produces CSER values per subject / source / band in *cser_values.mat*.

### 3. Complexity & PSD Analysis 

**Code:** analyse_cser.ipynb, analyse_lzc.ipynb, analyse_psd.ipynb

Analyses CSER results and performs LZC and PSD (power spectral density) calculations in Python. Produces significant observations and graphs for CSER, LZC and PSD in depressed subjects vs healthy control.

### 4. Transfer Entropy Analysis

**Code:** EEGtoMVGC.m (and related files)

Perform transfer entropy analysis using four different methods and various correction methods.

### 5. Predictive Processing Simulations

**Code:** analyse_pp.ipynb

Model and simulate predictive processing models of a depressed brain based on healthy data.

### TO-DO

The code is aimed to be cleaned up, further separated and grouped into units, and turned into a library for easier reusability.
