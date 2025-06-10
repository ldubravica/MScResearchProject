### Phase One: Depression_Study\phase_one_raw.ipynb

**Input:** depression EEG study files - Depression_Study\depression_data\matlab_files\###_Depression_REST.mat

**Function:** preprocceses the EEG in Python

**Output:** clean and filtered EEG per subject - Depression_Study\export_mat\###_Depression_REST_processed.mat

### Phase Two: EEGtoCSER.m (replacing PhaseTwo.m and CalcAvgER.m combination)

**Input:** preprocessed EEG files - Depression_Study\export_mat\###_Depression_REST_processed.mat

**Function:** performs source reconstruction and CSER analysis in MatLab

**Output:** CSER values per subject / source / band - cser_values.mat

### Phase Three: cser_analysis.ipynb

**Input:** CSER values - cser_values.mat

**Function:** analyses CSER results in Python

**Output:** significant observations regarding CSER in depression vs control
