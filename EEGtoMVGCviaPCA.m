function [mvgc_open, mvgc_closed, ss_info_open, ss_info_closed] = EEGtoMVGCviaPCA()

    %% Load .mat and .m files

    dataDir = fullfile('Depression_Study', 'export_mat'); % TEST directory
    files = dir(fullfile(dataDir, '*.mat')); % Get all .mat files in the directory

    addpath('Source_Reconstruction'); % Add path to Source_Reconstruction folder
    addpath('EntRate'); % Add path to EntRate folder

    % Load MVGC toolbox - silently adding MVGC to current path
    % p = mfilename('fullpath');
    % addpath(strrep(p, 'EEGtoMVGCviaPCA', 'EntRate/private/mvgc_v2.0'));
    % evalc('mvgc_startup;');

    %% Initialize variables

    mvgc_open = struct();
    mvgc_closed = struct();

    ss_info_open = struct();
    ss_info_closed = struct();

    frontal = [3:16, 19:24];
    occipital = [25:36];
    parietal = [39:50];
    sensorimotor = [1,2,17,18,37,38];
    temporal = [51:60];

    region_names = {'Frontal', 'Occipital', 'Parietal', 'Sensorimotor', 'Temporal'};
    region_indices = {frontal, occipital, parietal, sensorimotor, temporal};
    n_regions = 5;

    LOAD_MVGC_FILE = false;

    %% Perform source reconstruction, rotation and CSER calculation for each subject

    if ~LOAD_MVGC_FILE

        % Loop through each file
        % for i = 1:length(files)
        % for i = 55:length(files)
        for i = 1:1
            % Load the .mat file
            fprintf('\n*** PROCESSING FILE: %s *** \n', files(i).name);
            filePath = fullfile(dataDir, files(i).name);
            key = sprintf('x%s', files(i).name(1:3));

            % Reconstruct 60 sources from the EEG data
            fprintf('\n*** RECONSTRUCTING SOURCE *** \n\n');
            [source_ts_open, source_ts_closed, aals] = SourceReconMatlab(filePath);
            
            % Rotate [trials x sources x time] to [sources x time x trials]
            fprintf('\n*** ROTATING DATA *** \n');
            subj_open = permute(source_ts_open, [2, 3, 1]);
            subj_closed = permute(source_ts_closed, [2, 3, 1]);

            % Calculate and store MVGC for each subject's data
            fprintf('\n*** CALCULATING MVGC %s - (%i/%i) *** \n', files(i).name(1:3), i, length(files));
            [pwcgc_values_open, info_open] = calculateMVGC(subj_open);
            [pwcgc_values_closed, info_closed] = calculateMVGC(subj_closed);

            mvgc_open.(key) = pwcgc_values_open;
            mvgc_closed.(key) = pwcgc_values_closed;

            ss_info_open.(key) = info_open;
            ss_info_closed.(key) = info_closed;

            fprintf('\n*** RESULTS SAVED: %s *** \n', files(i).name);
        end
    
    end

    %% Function to calculate MVGC for a given subject's data

    function [pwcgc_values, info] = calculateMVGC(subj_data)

        % Downsampling the data
        fprintf('\nDownsampling data to 250 Hz...\n');
        disp(size(subj_data));
        subj_data = subj_data(:, 1:2:end, :);  % Downsample to 250 Hz (from 500 Hz)
        disp(size(subj_data));

        % Z-score each region across time for each trial
        fprintf('\nZ-scoring data...\n');
        [n_sources, n_samples, n_trials] = size(subj_data);
        for trial = 1:n_trials
            for source = 1:n_sources
                try
                    subj_data(source, :, trial) = zscore(subj_data(source, :, trial));
                catch
                    fprintf('Error in z-scoring source %d, trial %d. Adding small noise.\n', source, trial);
                    subj_data(source, :, trial) = zscore(subj_data(source, :, trial)) + randn(1, n_samples) * 0.001;  % Add small noise
                end 
            end
        end

        % Create empty 3D matrix of 10 PCAs x samples x trials
        pca_data = zeros(10, n_samples, n_trials);

        % Acquire top 2 PCA components for each region
        fprintf('\nPerforming PCA on each region...\n');
        n_components = 2; % Number of components to retain
        for region_idx = 1:n_regions
            region_set = region_indices{region_idx};
            fprintf('Region %s: %s\n', region_names{region_idx}, mat2str(region_set));

            % Extract data for the current region
            region_data = subj_data(region_set, :, :);
            [n_sources_region, n_samples, n_trials] = size(region_data);

            % Reshape to [sources x (samples * trials)]
            reshaped_data = reshape(region_data, n_sources_region, n_samples * n_trials);
            % Perform PCA
            [coeff, score, latent] = pca(reshaped_data', 'NumComponents', n_components);
            % Store the PCA scores in the pca_data matrix
            pca_data(2*region_idx-1:2*region_idx, :, :) = reshape(score', n_components, n_samples, n_trials);
        end

        %% Calculate MVGC between sources

        % Initialize MVGC matrix for region pairs
        pwcgc_values = zeros(n_regions, n_regions);

        % VAR model order estimation
        varmomax = 10;
        varmosel = 'AIC';
        fprintf('\nEstimating VAR model order (max order = %d)...\n', varmomax);
        [varmoaic, varmobic, varmohqc, varmolrt] = tsdata_to_varmo(pca_data, varmomax, 'LWR', [], []);
        varmo = moselect(sprintf('VAR model order selection (max = %d)', varmomax), ... 
                varmosel, 'AIC', varmoaic, 'BIC', varmobic, 'HQC', varmohqc, 'LRT', varmolrt);
        fprintf('Selected VAR model order (%s): %d\n', varmosel, varmo);

        % State-space model order estimation using SVC (faster)
        fprintf('\nEstimating state-space model order using SVC...\n');
        [ssmo_svc, ssmo_max] = tsdata_to_sssvc(pca_data, 2*varmo, [], []);
        ssmo = ssmo_svc;
        fprintf('Selected state-space model order (SVC): %d (max = %d)\n', ssmo, ssmo_max);
        
        % Estimate SS model
        fprintf('Estimating state-space model parameters...\n');
        [A, C, K, V] = tsdata_to_ss(pca_data, 2*varmo, ssmo);
        info = ss_info(A, C, K, V, 1);
        fprintf('  Sigma SPD: %d\n', info.sigspd);
        if info.error
            fprintf('State-space model estimation encountered errors.\n');
        else
            fprintf('State-space model estimation successful (no errors).\n\n');
        end

        % Loop through all pairs of regions
        for from_idx = 1:n_regions
            for to_idx = 1:n_regions

                if from_idx == to_idx
                    continue;  % Skip self-connections
                end

                from_indices = 2*from_idx-1:2*from_idx;
                to_indices = 2*to_idx-1:2*to_idx;

                % Compute MVGC from from_set to to_set
                try
                    mvgc_val = ss_to_mvgc(A, C, K, V, to_indices, from_indices);
                    fprintf('\nMVGC value (%d-%s -> %d-%s):\n', ...
                        from_idx, region_names{from_idx}, to_idx, region_names{to_idx});
                    disp(mvgc_val);
                    % Average MVGC value from all from->to connections
                    pwcgc_values(from_idx, to_idx) = mean(mvgc_val(:));
                catch
                    pwcgc_values(from_idx, to_idx) = NaN;
                    fprintf('\nMVGC (%d-%s -> %d-%s) computation failed\n', ...
                        from_idx, region_names{from_idx}, to_idx, region_names{to_idx});
                end
            end
        end

        fprintf('\nMVGC computation complete.\n');

    end

    %% Load MVGC results from mvgc_values.mat

    if LOAD_MVGC_FILE

        fprintf('\n*** LOADING MVGC RESULTS *** \n');
        if exist('mvgc_values.mat', 'file')
            load('mvgc_values.mat', 'mvgc_open', 'mvgc_closed');
            fprintf('MVGC results loaded successfully.\n');
        else
            fprintf('No MVGC results found. Please run EEGtoMVGCviaPCA first.\n');
            return;
        end

    end

    %% Split mvgc_open and mvgc_closed based on healthy/depressed subjects

    fprintf('\n*** SPLITTING MVGC RESULTS *** \n');

    % Initialize structures for healthy and depressed subjects
    mvgc_open_healthy = struct();
    mvgc_open_depressed = struct();
    mvgc_closed_healthy = struct();
    mvgc_closed_depressed = struct();

    % Load healthy and depressed subject lists
    healthy_depressed_samples = load('healthy_depressed_samples.mat');
    healthy_sample = healthy_depressed_samples.healthy_sample;
    depressed_sample = healthy_depressed_samples.depressed_sample;

    % Split MVGC results based on whether subject name appears in healthy or depressed lists
    for i = 1:length(files)
        key = sprintf('x%s', files(i).name(1:3));
        % Define variable key_int as files(i).name(1:3), but converted to an integer
        key_int = str2double(files(i).name(1:3));

        if isfield(mvgc_open, key)
            if ismember(key_int, healthy_sample)
                mvgc_open_healthy.(key) = mvgc_open.(key);
            elseif ismember(key_int, depressed_sample)
                mvgc_open_depressed.(key) = mvgc_open.(key);
            end
        end
        if isfield(mvgc_closed, key)
            if ismember(key_int, healthy_sample)
                mvgc_closed_healthy.(key) = mvgc_closed.(key);
            elseif ismember(key_int, depressed_sample)
                mvgc_closed_depressed.(key) = mvgc_closed.(key);
            end
        end
    end

    %% Average MVGC across open/closed healthy/depressed subjects

    fprintf('\n*** AVERAGING MVGC ACROSS SUBJECTS *** \n');

    % Helper function to average MVGC matrices in a struct
    function avg_mat = average_mvgc_struct(mvgc_struct)
        avg_mat = zeros(length(region_names), length(region_names));
        count = 0;
        for i = 1:length(files)
            key = sprintf('x%s', files(i).name(1:3));
            if isfield(mvgc_struct, key)
                avg_mat = avg_mat + mvgc_struct.(key);
                count = count + 1;
            end
        end
        if count > 0
            avg_mat = avg_mat / count;
        end
    end

    mvgc_open_healthy_avg = average_mvgc_struct(mvgc_open_healthy);
    mvgc_open_depressed_avg = average_mvgc_struct(mvgc_open_depressed);
    mvgc_closed_healthy_avg = average_mvgc_struct(mvgc_closed_healthy);
    mvgc_closed_depressed_avg = average_mvgc_struct(mvgc_closed_depressed);

    mvgc_open_avg = (mvgc_open_healthy_avg + mvgc_open_depressed_avg) / 2;
    mvgc_closed_avg = (mvgc_closed_healthy_avg + mvgc_closed_depressed_avg) / 2;
    mvgc_healthy_avg = (mvgc_open_healthy_avg + mvgc_closed_healthy_avg) / 2;
    mvgc_depressed_avg = (mvgc_open_depressed_avg + mvgc_closed_depressed_avg) / 2;

    fprintf('\n*** AVERAGE MVGC OPEN HEALTHY *** \n');
    disp(mvgc_open_healthy_avg);
    fprintf('\n*** AVERAGE MVGC OPEN DEPRESSED *** \n');
    disp(mvgc_open_depressed_avg);
    fprintf('\n*** AVERAGE MVGC CLOSED HEALTHY *** \n');
    disp(mvgc_closed_healthy_avg);
    fprintf('\n*** AVERAGE MVGC CLOSED DEPRESSED *** \n');
    disp(mvgc_closed_depressed_avg);
    fprintf('\n*** AVERAGE MVGC OPEN *** \n');
    disp(mvgc_open_avg);
    fprintf('\n*** AVERAGE MVGC CLOSED *** \n');
    disp(mvgc_closed_avg);
    fprintf('\n*** AVERAGE MVGC HEALTHY *** \n');
    disp(mvgc_healthy_avg);
    fprintf('\n*** AVERAGE MVGC DEPRESSED *** \n');
    disp(mvgc_depressed_avg);

    %% Plot average MVGC matrices across open/closed healthy/depressed subjects
    % plot_mvgc_matrix_grid(...
    %     {mvgc_open_healthy_avg, mvgc_open_depressed_avg, mvgc_closed_healthy_avg, mvgc_closed_depressed_avg}, ...
    %     region_names, ...
    %     {'Open Healthy', 'Open Depressed', 'Closed Healthy', 'Closed Depressed'}, ...
    %     'Average MVGC Matrices');

    % plot_mvgc_matrix_grid(...
    %     {mvgc_open_avg, mvgc_closed_avg, mvgc_healthy_avg, mvgc_depressed_avg}, ...
    %     region_names, ...
    %     {'Open', 'Closed', 'Healthy', 'Depressed'}, ...
    %     'Average MVGC Matrices');

    function plot_mvgc_matrix_grid(mvgc_matrices, region_names, titles, super_title)
        % mvgc_matrices: cell array of 4 matrices
        % region_names: cell array of region names
        % titles: cell array of 4 subplot titles
        % super_title: string for sgtitle

        % Find global color limits
        all_vals = [];
        for i = 1:numel(mvgc_matrices)
            all_vals = [all_vals; mvgc_matrices{i}(:)];
        end
        clim = [min(all_vals), max(all_vals)];

        figure;
        for i = 1:4
            subplot(2,2,i);
            imagesc(mvgc_matrices{i}, clim);
            colorbar;
            title(sprintf('Average MVGC - %s', titles{i}));
            set(gca, 'XTick', 1:length(region_names), 'XTickLabel', region_names, ...
                     'YTick', 1:length(region_names), 'YTickLabel', region_names);
            xlabel('To Region');
            ylabel('From Region');
            axis square;
        end
        sgtitle(super_title);
    end

    %% Save results into .mat files

    % save(fullfile('mvgc_values.mat'), ...
    %      'mvgc_open', 'mvgc_closed', ...
    %      'mvgc_open_healthy', 'mvgc_open_depressed', ...
    %      'mvgc_closed_healthy', 'mvgc_closed_depressed', ...
    %      'mvgc_open_healthy_avg', 'mvgc_open_depressed_avg', ...
    %      'mvgc_closed_healthy_avg', 'mvgc_closed_depressed_avg');

    fprintf('\n*** FILES SAVED *** \n');

end