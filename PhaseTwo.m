function [cser_source_open, cser_source_closed, cser_source_band_open, cser_source_band_closed] = PhaseTwo()
    % Define the directory containing the .mat files
    dataDir = fullfile('Depression_Study', 'export_mat'); % TEST directory
    files = dir(fullfile(dataDir, '*.mat')); % Get all .mat files in the directory

    cser_source_open = struct();
    cser_source_closed = struct();
    cser_source_band_open = struct();
    cser_source_band_closed = struct();

    % Loop through each file
    for i = 1:length(files)
        % Load the .mat file
        fprintf('\n*** PROCESSING FILE: %s *** \n', files(i).name);
        filePath = fullfile(dataDir, files(i).name);

        % Reconstruct 60 sources from the EEG data
        fprintf('\n*** RECONSTRUCTING SOURCE *** \n\n');
        addpath('Source_Reconstruction'); % Add path to Source_Reconstruction folder
        [source_ts_open, source_ts_closed, aals] = SourceRecon_matlab(filePath);

        % Rotating [trials x sources x time] to [sources x time x trials]
        fprintf('\n*** ROTATING DATA *** \n');
        rotated_open = permute(source_ts_open, [2, 3, 1]);
        rotated_closed = permute(source_ts_closed, [2, 3, 1]);

        % Calculate entropy rate per source for open and closed epochs
        fprintf('\n*** CALCULATING ENTROPY RATE *** \n');
        addpath('EntRate'); % Add path to EntRate folder
        Fs = 500; % Sampling frequency in Hz
        bands = [1, 4; 4, 8; 8, 12; 12, 30; 30, 100];
        % H_open = StateSpaceEntropyRatePerSource(rotated_open, Fs);
        % H_closed = StateSpaceEntropyRatePerSource(rotated_closed, Fs);
        [H_open, bH_open] = StateSpaceEntropyRatePerSource(rotated_open, Fs, 'yes', bands);
        [H_closed, bH_closed] = StateSpaceEntropyRatePerSource(rotated_closed, Fs, 'yes', bands);

        % Store results
        key = sprintf('x%s', files(i).name(1:3));
        cser_source_open.(key) = H_open;
        cser_source_closed.(key) = H_closed;
        cser_source_band_open.(key) = bH_open;
        cser_source_band_closed.(key) = bH_closed;
        fprintf('\n*** RESULTS SAVED: %s *** \n', files(i).name);
    end

    % Save results to .mat files
    % save(fullfile('entropy_source_open.mat'), 'entropy_source_open');
    % save(fullfile('entropy_source_closed.mat'), 'entropy_source_closed');
    % save(fullfile('entropy_source_band_open.mat'), 'entropy_source_band_open');
    % save(fullfile('entropy_source_band_closed.mat'), 'entropy_source_band_closed');

    % Save all results in a single file
    save(fullfile('cser_values.mat'), ...
        'cser_source_open', 'cser_source_closed', 'cser_source_band_open', 'cser_source_band_closed');

    fprintf('\n*** FILES SAVED *** \n');
end