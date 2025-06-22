### Phase One: Depression_Study\phase_one_raw.ipynb

**Input:** depression EEG study files - Depression_Study\depression_data\matlab_files\XXX_Depression_REST.mat

**Function:** preprocceses the EEG in Python

**Output:** clean and filtered EEG per subject - Depression_Study\export_mat\XXX_Depression_REST_processed.mat

### Phase Two: EEGtoCSER.m

**Input:** preprocessed EEG files - Depression_Study\export_mat\XXX_Depression_REST_processed.mat

**Function:** performs source reconstruction and CSER calculation in MatLab

**Output:** CSER values per subject / source / band - cser_values.mat (as well as source reconstructed EEG)

### Phase Three: analyse_cser.ipynb, analyse_lzc.ipynb, analyse_psd.ipynb

**Input:** CSER values (cser_values.mat) or reconstructed EEG source (eeg_source.mat)

**Function:** analyses CSER results and perform LZC and PSD calculations in Python

**Output:** significant observations regarding CSER, LZC and PSD in depression vs control
