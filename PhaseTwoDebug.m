function [entropy_open, entropy_closed, entropy_bands_open, entropy_bands_closed] = PhaseTwoDebug()
    % Define the directory containing the .mat files
    dataDir = fullfile('Depression_Study', 'export_mat'); % TEST directory
    files = dir(fullfile(dataDir, '*.mat')); % Get all .mat files in the directory

    entropy_open = struct();
    entropy_closed = struct();
    entropy_bands_open = struct();
    entropy_bands_closed = struct();

    %%% REGULAR %%%

    % Loop through each file
    % for i = 1:length(files)
    %     [key, H_open, H_closed] = processfile(dataDir, files(i).name);
    %     entropy_open.(key) = H_open;
    %     entropy_closed.(key) = H_closed;
    %     fprintf('\n*** RESULTS SAVED: %s *** \n', files(i).name);
    % end

    %%% DEBUG %%%

    filename = '527_Depression_REST_processed.mat';
    if isfile(fullfile(dataDir, filename))
        [key, H_open, H_closed, bH_open, bH_closed] = processfile(dataDir, filename);
        entropy_open.(key) = H_open;
        entropy_closed.(key) = H_closed;
        entropy_bands_open.(key) = bH_open;
        entropy_bands_closed.(key) = bH_closed;
        fprintf('\n*** RESULTS SAVED: %s *** \n', filename);
    else
        fprintf('\n*** SPECIFIC FILE NOT FOUND: %s *** \n', filename);
    end

    %%% EXPORT FILES %%%

    % Save results to .mat files
    % save(fullfile('entropy_open.mat'), 'entropy_open');
    % save(fullfile('entropy_closed.mat'), 'entropy_closed');
    % fprintf('\n*** FILES SAVED *** \n');
end

% PhaseTwo Pipeline function
function [key, H_open, H_closed, bH_open, bH_closed] = processfile(dataDir, filename)
    % Load the .mat file
    fprintf('\n*** PROCESSING FILE: %s *** \n', filename);
    filePath = fullfile(dataDir, filename);

    % Reconstruct 60 sources from the EEG data
    fprintf('\n*** RECONSTRUCTING SOURCE *** \n\n');
    addpath('Source_Reconstruction'); % Add path to Source_Reconstruction folder
    [source_ts_open, source_ts_closed, aals] = SourceReconMatlab(filePath);

    % Rotating [trials x channels x time] to [channels x time x trials]
    fprintf('\n*** ROTATING DATA *** \n');
    rotated_open = permute(source_ts_open, [2, 3, 1]);
    rotated_closed = permute(source_ts_closed, [2, 3, 1]);

    % Calculate entropy rate per channel for open and closed epochs
    fprintf('\n*** CALCULATING ENTROPY RATE *** \n');
    addpath('EntRate'); % Add path to EntRate folder
    Fs = 500; % Sampling frequency in Hz
    bands = [1, 4; 4, 8; 8, 12; 12, 30; 30, 100];
    % H_open = StateSpaceEntropyRatePerChannel(rotated_open, Fs);
    % H_closed = StateSpaceEntropyRatePerChannel(rotated_closed, Fs);
    [H_open, bH_open] = StateSpaceEntropyRatePerChannel(rotated_open, Fs, 'yes', bands);
    [H_closed, bH_closed] = StateSpaceEntropyRatePerChannel(rotated_closed, Fs, 'yes', bands);

    key = sprintf('x%s', filename(1:3));
end
