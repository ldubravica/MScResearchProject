function [cser_values, eeg_source, aals] = EEGtoCSER()

    %% Load .mat and .m files

    dataDir = fullfile('Depression_Study', 'export_mat'); % TEST directory
    files = dir(fullfile(dataDir, '*.mat')); % Get all .mat files in the directory

    addpath('Source_Reconstruction'); % Add path to Source_Reconstruction folder
    addpath('EntRate'); % Add path to EntRate folder

    %% Prepare structs

    cser_open = struct();
    cser_closed = struct();
    cser_source_open = struct();
    cser_source_closed = struct();
    cser_band_open = struct();
    cser_band_closed = struct();
    cser_source_band_open = struct();
    cser_source_band_closed = struct();

    eeg_source_open = struct();
    eeg_source_closed = struct();

    source_names = [];

    %% Perform source reconstruction, rotation and CSER calculation for each subject

    % Loop through each file
    for i = 1:length(files)
        % Load the .mat file
        fprintf('\n*** PROCESSING FILE: %s *** \n', files(i).name);
        filePath = fullfile(dataDir, files(i).name);

        % Reconstruct 60 sources from the EEG data
        fprintf('\n*** RECONSTRUCTING SOURCE *** \n\n');
        [source_ts_open, source_ts_closed, aals] = SourceReconMatlab(filePath);

        % Rotate [trials x sources x time] to [sources x time x trials]
        fprintf('\n*** ROTATING DATA *** \n');
        rotated_open = permute(source_ts_open, [2, 3, 1]);
        rotated_closed = permute(source_ts_closed, [2, 3, 1]);

        % Calculate entropy rate per source & per source per band
        fprintf('\n*** CALCULATING ENTROPY RATE %s - (%i/%i) *** \n', files(i).name(1:3), i, length(files));
        Fs = 500; % Sampling frequency in Hz
        basic_bands = [1, 4; 4, 8; 8, 12; 12, 30; 30, 100];
        custom_bands = [1, 12];
        bands = [basic_bands; custom_bands];
        [H_src_open, bH_src_open] = StateSpaceEntropyRatePerSource(rotated_open, Fs, 'yes', bands);
        [H_src_closed, bH_src_closed] = StateSpaceEntropyRatePerSource(rotated_closed, Fs, 'yes', bands);

        % Calculate entropy rate per subject & per subject per band
        H_open = mean(H_src_open, 'omitnan');
        H_closed = mean(H_src_closed, 'omitnan');
        bH_open = mean(bH_src_open, 2, 'omitnan');
        bH_closed = mean(bH_src_closed, 2, 'omitnan');

        % Store results
        key = sprintf('x%s', files(i).name(1:3));

        cser_open.(key) = H_open;  % CSER / subject (eyes open)
        cser_closed.(key) = H_closed;  % CSER / subject (eyes closed)
        cser_source_open.(key) = H_src_open;  % CSER / subject / source (eyes open)
        cser_source_closed.(key) = H_src_closed;  % CSER / subject / source (eyes closed)
        cser_band_open.(key) = bH_open';  % CSER / subject / band (eyes open)
        cser_band_closed.(key) = bH_closed';  % CSER / subject / band (eyes closed)
        cser_source_band_open.(key) = bH_src_open;  % CSER / subject / source / band (eyes open)
        cser_source_band_closed.(key) = bH_src_closed;  % CSER / subject / source / band (eyes closed)

        % EEG reconstructed from channels to sources
        eeg_source_open.(key) = source_ts_open;
        eeg_source_closed.(key) = source_ts_closed;

        if isempty(source_names)
            source_names = aals;
        end

        fprintf('\n*** RESULTS SAVED: %s *** \n', files(i).name);
    end

    %% Save results into .mat files

    save(fullfile('cser_values.mat'), ...
        'cser_open', 'cser_closed', 'cser_source_open', 'cser_source_closed', ...
        'cser_band_open', 'cser_band_closed', 'cser_source_band_open', 'cser_source_band_closed', ...
        'source_names');

    % save(fullfile('eeg_source.mat'), ...
    %     'eeg_source_open', 'eeg_source_closed', 'source_names', '-v7.3');

    %% Save results in a structured format for return
    
    cser_values = struct(...
        'cser_open', cser_open, ...
        'cser_closed', cser_closed, ...
        'cser_source_open', cser_source_open, ...
        'cser_source_closed', cser_source_closed, ...
        'cser_band_open', cser_band_open, ...
        'cser_band_closed', cser_band_closed, ...
        'cser_source_band_open', cser_source_band_open, ...
        'cser_source_band_closed', cser_source_band_closed);
    
    % eeg_source = struct(...
    %     'eeg_source_open', eeg_source_open, ...
    %     'eeg_source_closed', eeg_source_closed);

    fprintf('\n*** FILES SAVED *** \n');

end