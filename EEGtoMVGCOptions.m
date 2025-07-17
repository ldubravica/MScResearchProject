function [mvgc_open, mvgc_closed, ss_info_open, ss_info_closed] = EEGtoMVGCOptions(mvgc_format, load_mvgc_file, save_mvgc_file, plot_data, save_plots)
    
    % EEGtoMVGCOptions - Main function for EEG to MVGC analysis
    % Usage:
    %   EEGtoMVGCOptions('region')

    if nargin < 1 || isempty(mvgc_format)
        mvgc_format = 'source';  % Options: 'source', 'region', 'pca', 'brain_source'
    end
    if nargin < 2 || isempty(load_mvgc_file)
        load_mvgc_file = false;
    end
    if nargin < 3 || isempty(save_mvgc_file)
        save_mvgc_file = false;  % Save results to .mat file
    end
    if nargin < 4 || isempty(plot_data)
        plot_data = false;  % Plot results
    end
    if nargin < 5 || isempty(save_plots)
        save_plots = false;  % Save plots to files
    end

    % filename = 'mvgc_values.mat';
    filename = sprintf('mvgc_results_%s.mat', mvgc_format);

    %% Load .mat and .m files

    dataDir = fullfile('Depression_Study', 'export_mat'); % TEST directory
    files = dir(fullfile(dataDir, '*.mat')); % Get all .mat files in the directory

    addpath('./Source_Reconstruction'); % Add path to Source_Reconstruction folder
    addpath('./EntRate'); % Add path to EntRate folder
    addpath('./fdr_bh');

    % Load MVGC toolbox - silently adding MVGC to current path
    % p = mfilename('fullpath');
    % addpath(strrep(p, 'EEGtoMVGC', 'EntRate/private/mvgc_v2.0'));
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

    n_subjects = length(files);

    %% Perform source reconstruction, rotation, prep and MVGC calculation for each subject

    if ~load_mvgc_file

        % Loop through each file
        for file_idx = 1:length(files)
        % for file_idx = 1:1
            % Load the .mat file
            fprintf('\n*** PROCESSING FILE: %s *** \n', files(file_idx).name);
            filePath = fullfile(dataDir, files(file_idx).name);
            key = sprintf('x%s', files(file_idx).name(1:3));

            % Reconstruct 60 sources from the EEG data
            fprintf('\n*** RECONSTRUCTING SOURCE *** \n\n');
            [source_ts_open, source_ts_closed, ~] = SourceReconMatlab(filePath);
            
            % Rotate [trials x sources x time] to [sources x time x trials]
            fprintf('\n*** ROTATING DATA *** \n');
            subj_open = permute(source_ts_open, [2, 3, 1]);
            subj_closed = permute(source_ts_closed, [2, 3, 1]);

            % Downsample and standardize the data
            fprintf('\n*** DOWNSAMPLING AND STANDARDIZING DATA *** \n');
            subj_open = downsampleAndStandardize(subj_open);
            subj_closed = downsampleAndStandardize(subj_closed);

            % Calculate and store MVGC for each subject's data
            fprintf('\n*** CALCULATING MVGC %s - (%i/%i) *** \n', files(file_idx).name(1:3), file_idx, length(files));

            if strcmp(mvgc_format, 'source')
                [pwcgc_values_open, info_open] = calculateMVGCSource(subj_open);
                [pwcgc_values_closed, info_closed] = calculateMVGCSource(subj_closed);
            elseif strcmp(mvgc_format, 'region')
                [pwcgc_values_open, info_open] = calculateMVGCRegion(subj_open);
                [pwcgc_values_closed, info_closed] = calculateMVGCRegion(subj_closed);
            elseif strcmp(mvgc_format, 'pca')
                [pwcgc_values_open, info_open] = calculateMVGCviaPCA(subj_open);
                [pwcgc_values_closed, info_closed] = calculateMVGCviaPCA(subj_closed);
            elseif strcmp(mvgc_format, 'brain_source')
                [pwcgc_values_open, info_open] = calculateMVGCBrainSource(subj_open);
                [pwcgc_values_closed, info_closed] = calculateMVGCBrainSource(subj_closed);
            else
                error('Invalid mvgc_format specified. Choose from: source, region, pca, brain_source');
            end

            % Display results
            fprintf('\nMVGC computation for %s completed.\n', files(file_idx).name);
            fprintf('\nOpen condition MVGC values:\n');
            disp(pwcgc_values_open);
            fprintf('\nClosed condition MVGC values:\n');
            disp(pwcgc_values_closed);

            % fprintf('\nMVGC computation complete.\n');

            mvgc_open.(key) = pwcgc_values_open;
            mvgc_closed.(key) = pwcgc_values_closed;

            ss_info_open.(key) = info_open;
            ss_info_closed.(key) = info_closed;

            fprintf('\n*** RESULTS SAVED: %s *** \n', files(file_idx).name);
        end
    
    end

    function [subj_data] = downsampleAndStandardize(subj_data)

        % Downsampling the data
        subj_data = subj_data(:, 1:2:end, :);  % Downsample to 250 Hz (from 500 Hz)

        % Z-score each region across time for each trial
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

        % disp(eig(cov(subj_data(:, :, 1)')));  % Display eigenvalues of covariance matrix for open condition

    end

    %% Function to calculate MVGC for a given subject's data

    function [A, C, K, V, info] = estimateVARSS(X)
        varmomax = 10;
        varmosel = 'AIC';

        fprintf('\nEstimating VAR order (max=%d)...\n', varmomax);
        [aic, bic, hqc, lrt] = tsdata_to_varmo(X, varmomax, 'LWR', [], []);
        varmo = moselect('VAR Order Selection', varmosel, 'AIC', aic, 'BIC', bic, 'HQC', hqc, 'LRT', lrt);
        % fprintf('Selected VAR order: %d (%s)\n', varmo, varmosel);

        fprintf('Estimating SS model order (SVC) and parameters...\n');
        [ssmo, ~] = tsdata_to_sssvc(X, 2*varmo, [], []);
        [A, C, K, V] = tsdata_to_ss(X, 2*varmo, ssmo);
        info = ss_info(A, C, K, V, 1);
        % fprintf('  Sigma SPD: %d\n', info.sigspd);

        if info.error
            % fprintf('⚠️ SS model estimation encountered errors.\n');
            fprintf('⚠️ SS model estimation encountered errors: %s\n', dec2bin(info.error));
        else
            fprintf('SS model estimation successful.\n');
        end
    end

    % Fits SS to the region pair subset of the data - runs ss_to_mvgc with indeces of the data
    function [pwcgc_values, info] = calculateMVGCSource(subj_data)

        % Initialize MVGC matrix for region pairs
        pwcgc_values = zeros(n_regions, n_regions);

        % Loop through all pairs of regions
        for from_idx = 1:n_regions
            for to_idx = 1:n_regions

                if from_idx == to_idx
                    continue; % Skip self-connections
                end

                fprintf('\n** MVGC (%d-%s -> %d-%s) **\n', ...
                        from_idx, region_names{from_idx}, to_idx, region_names{to_idx});

                from_set = region_indices{from_idx};
                to_set = region_indices{to_idx};
                combined_set = [from_set, to_set];
                X = subj_data(combined_set, :, :);

                % Model order estimation
                [A, C, K, V, info] = estimateVARSS(X);

                % Calculate from_indices and to_indices
                from_set_X = 1:length(from_set);
                to_set_X = (length(from_set) + 1):(length(from_set) + length(to_set));
                
                % Compute MVGC from from_set to to_set
                try
                    mvgc_val = ss_to_mvgc(A, C, K, V, to_set_X, from_set_X);
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

    end

    % Fits SS to the entire data of 60 channels - runs ss_to_mvgc with indeces of the data 
    function [pwcgc_values, info] = calculateMVGCBrainSource(subj_data)

        % Initialize MVGC matrix for region pairs
        pwcgc_values = zeros(n_regions, n_regions);

        % Model order estimation
        [A, C, K, V, info] = estimateVARSS(subj_data);

        % Loop through all pairs of regions
        for from_idx = 1:n_regions
            for to_idx = 1:n_regions

                if from_idx == to_idx
                    continue; % Skip self-connections
                end

                fprintf('\n* MVGC (%d-%s -> %d-%s) *\n', ...
                        from_idx, region_names{from_idx}, to_idx, region_names{to_idx});

                from_set = region_indices{from_idx};
                to_set = region_indices{to_idx};
                combined_set = [from_set, to_set];

                % Find the indices of from_set and to_set within combined_set
                % combined_set = sort(combined_set);
                [~, from_indices_sub] = ismember(from_set, combined_set);
                [~, to_indices_sub] = ismember(to_set, combined_set);
                
                % Compute MVGC from from_set to to_set
                try
                    mvgc_val = ss_to_mvgc(A, C, K, V, to_indices_sub, from_indices_sub);
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

    end

    % Fits SS to the 5 means of regional EEGs - runs ss_to_pwgc with pairwise of 5 regions
    function [pwcgc_values, info] = calculateMVGCRegion(subj_data)

        [~, n_samples, n_trials] = size(subj_data);
        regional_eeg = zeros(n_regions, n_samples, n_trials);  % [regions x samples x trials]

        for r = 1:n_regions
            regional_data = subj_data(region_indices{r}, :, :);
            regional_eeg(r, :, :) = squeeze(mean(regional_data, 1));  % Average across sources
        end

        % Z-score each region across time for each trial
        fprintf('\nZ-scoring each region across time for each trial...\n');
        for trial = 1:n_trials
            for region = 1:n_regions
                regional_eeg(region, :, trial) = zscore(regional_eeg(region, :, trial));
            end
        end

        % Model order estimation
        [A, C, K, V, info] = estimateVARSS(regional_eeg);

        % Compute time-domain pairwise conditional GC
        pwcgc_values = ss_to_pwcgc(A, C, K, V);

        % Validate result
        assert(~isbad(pwcgc_values, false), 'pwcgc_values estimation failed');

    end

    % Fits SS to the 10 PCA components (top 2 per 5 regions) - runs ss_to_mvgc with 2 pairs of components
    function [pwcgc_values, info] = calculateMVGCviaPCA(subj_data)

        n_components = 2;  % Number of components to retain
        [~, n_samples, n_trials] = size(subj_data);
        pca_data = zeros(n_components * n_regions, n_samples, n_trials);

        % Acquire top 2 PCA components for each region
        fprintf('\nPerforming PCA on each region...\n');
        for r = 1:n_regions
            % Extract data for the current region
            region_data = subj_data(region_indices{r}, :, :);
            [n_sources_region, n_samples, n_trials] = size(region_data);

            % Reshape to [sources x (samples * trials)]
            reshaped = reshape(region_data, n_sources_region, n_samples * n_trials);
            % Perform and store PCA
            [coeff, score, latent] = pca(reshaped', 'NumComponents', n_components);
            pca_data((2*r-1):(2*r), :, :) = reshape(score', n_components, n_samples, n_trials);
        end

        % Initialize MVGC matrix for region pairs
        pwcgc_values = zeros(n_regions, n_regions);

        % Model order estimation
        [A, C, K, V, info] = estimateVARSS(pca_data);

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

    end

    %% Load MVGC results from mvgc_values.mat

    if load_mvgc_file

        fprintf('\n*** LOADING MVGC RESULTS *** \n');
        if exist(filename, 'file')
            load(filename, 'mvgc_open', 'mvgc_closed');
            fprintf('MVGC results loaded successfully.\n');
        else
            fprintf('%s not found. Please run such MVGC format first.\n', filename);
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
    for file_idx = 1:length(files)
        key = sprintf('x%s', files(file_idx).name(1:3));
        key_int = str2double(files(file_idx).name(1:3));

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

    % Create a binary list of whether each subject is healthy or depressed
    % disp(depressed_sample)
    depressed_binary_labels = zeros(1, length(files));
    for file_idx = 1:length(files)
        key = sprintf('x%s', files(file_idx).name(1:3));
        key_int = str2double(files(file_idx).name(1:3));

        if isfield(mvgc_open, key)
            if ismember(key_int, depressed_sample)
                depressed_binary_labels(file_idx) = 1;
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

    %% Load Depression_Study\depression_data\Data_4_Import_REST.xlsx

    excelFile = fullfile('Depression_Study', 'depression_data', 'Data_4_Import_REST.xlsx');
    T = readtable(excelFile);
    % T.id, T.MDD, T.sex, T.age, T.BDI, T.BDI_Anh, T.BDI_Mel, T.TAI

    %% Calculate p-values and t-values for each MVGC matrix

    fprintf('\n*** CALCULATING P-VALUES AND T-VALUES *** \n\n');

    % Convert MVGC struct to 3D matrix
    mvgc_open_3d   = struct_to_3d(mvgc_open);
    mvgc_closed_3d = struct_to_3d(mvgc_closed);

    function mvgc_3d = struct_to_3d(mvgc_struct)
        fields = fieldnames(mvgc_struct);
        mvgc_3d = zeros(n_regions, n_regions, n_subjects);
        for i = 1:length(fields)
            mvgc_3d(:, :, i) = mvgc_struct.(fields{i});
        end
    end

    % Initialize results
    mvgc_open_pvals   = zeros(n_regions, n_regions);
    mvgc_open_tvals   = zeros(n_regions, n_regions);
    mvgc_closed_pvals = zeros(n_regions, n_regions);
    mvgc_closed_tvals = zeros(n_regions, n_regions);

    n_tests = n_regions^2 - n_regions;
    open_p_list   = zeros(n_tests, 1);
    closed_p_list = zeros(n_tests, 1);
    index_map     = zeros(n_tests, 2);

    map_idx = 1;

    for r1 = 1:n_regions
        for r2 = 1:n_regions
            if r1 == r2
                continue;  % Skip self-connections
            end

            index_map(map_idx, :) = [r1, r2];

            mvgc_r1_r2_open = squeeze(mvgc_open_3d(r1, r2, :))';
            [p_open, t_open] = get_p_t_value(mvgc_r1_r2_open);
            open_p_list(map_idx) = p_open;
            mvgc_open_pvals(r1, r2) = p_open;
            mvgc_open_tvals(r1, r2) = t_open;

            mvgc_r1_r2_closed = squeeze(mvgc_closed_3d(r1, r2, :))';
            [p_closed, t_closed] = get_p_t_value(mvgc_r1_r2_closed);
            closed_p_list(map_idx) = p_closed;
            mvgc_closed_pvals(r1, r2) = p_closed;
            mvgc_closed_tvals(r1, r2) = t_closed;

            map_idx = map_idx + 1;

            fprintf('Region %d-%s to %d-%s: p_open = %.3f, t_open = %.3f, p_closed = %.3f, t_closed = %.3f\n', ...
                    r1, region_names{r1}, r2, region_names{r2}, p_open, t_open, p_closed, t_closed);
        end
    end

    function [p, t] = get_p_t_value(x)
        healthy   = x(depressed_binary_labels == 0);
        depressed = x(depressed_binary_labels == 1);
        [~, p, ~, stats] = ttest2(healthy, depressed);
        t = stats.tstat;
    end

    %% Calculate p-values and t-values using Robust Linear Model

    fprintf('\n*** CALCULATING ROBUST LINEAR MODEL P-VALUES AND T-VALUES *** \n\n');

    % Prepare subject IDs from file names
    subject_ids = zeros(n_subjects, 1);
    for i = 1:n_subjects
        subject_ids(i) = str2double(files(i).name(1:3));
    end
    % Find indices in T corresponding to subject_ids
    [~, T_idx] = ismember(subject_ids, T.id);

    % Extract covariates
    sex = T.sex(T_idx) - 1; % shift to 0-1 binary
    age = T.age(T_idx);
    depressed = depressed_binary_labels';

    % Initialize matrices
    mvgc_open_rlm_pvals_included = zeros(n_regions, n_regions);
    mvgc_open_rlm_tvals_included = zeros(n_regions, n_regions);
    mvgc_closed_rlm_pvals_included = zeros(n_regions, n_regions);
    mvgc_closed_rlm_tvals_included = zeros(n_regions, n_regions);

    mvgc_open_rlm_pvals_residuals = zeros(n_regions, n_regions);
    mvgc_open_rlm_tvals_residuals = zeros(n_regions, n_regions);
    mvgc_closed_rlm_pvals_residuals = zeros(n_regions, n_regions);
    mvgc_closed_rlm_tvals_residuals = zeros(n_regions, n_regions);

    for r1 = 1:n_regions
        for r2 = 1:n_regions
            if r1 == r2
                continue;
            end

            fprintf('\n** RLM (%d-%s -> %d-%s) **\n', r1, region_names{r1}, r2, region_names{r2});

            % Open condition
            y_open = squeeze(mvgc_open_3d(r1, r2, :)); % 119 x 1

            % Open condition - residuals
            tbl_open = table(y_open, sex, age, 'VariableNames', {'MVGC','Sex','Age'});
            mdl_open = fitlm(tbl_open, 'MVGC ~ Sex + Age', 'RobustOpts', 'on');
            residuals = mdl_open.Residuals.Raw;  % 119 x 1

            healthy = residuals(depressed_binary_labels == 0);
            depressed = residuals(depressed_binary_labels == 1);
            [~, p, ~, stats] = ttest2(healthy, depressed);
            t = stats.tstat;
            mvgc_open_rlm_pvals_residuals(r1, r2) = p;
            mvgc_open_rlm_tvals_residuals(r1, r2) = t;

            % Open condition - with Depressed covariate
            tbl_open = table(y_open, sex, age, depressed, 'VariableNames', {'MVGC','Sex','Age','Depressed'});
            mdl_open = fitlm(tbl_open, 'MVGC ~ Sex + Age + Depressed', 'RobustOpts', 'on');
            coef_open = mdl_open.Coefficients;
            mvgc_open_rlm_pvals_included(r1, r2) = coef_open.pValue(4); % Depressed p-value
            mvgc_open_rlm_tvals_included(r1, r2) = coef_open.tStat(4);  % Depressed t-value

            % Closed condition
            y_closed = squeeze(mvgc_closed_3d(r1, r2, :));  % 119 x 1

            % Closed condition - residuals
            tbl_closed = table(y_closed, sex, age, 'VariableNames', {'MVGC','Sex','Age'});
            mdl_closed = fitlm(tbl_closed, 'MVGC ~ Sex + Age', 'RobustOpts', 'on');
            residuals = mdl_closed.Residuals.Raw; % 119 x 1

            healthy = residuals(depressed_binary_labels == 0);
            depressed = residuals(depressed_binary_labels == 1);
            [~, p, ~, stats] = ttest2(healthy, depressed);
            t = stats.tstat;
            mvgc_closed_rlm_pvals_residuals(r1, r2) = p;
            mvgc_closed_rlm_tvals_residuals(r1, r2) = t;

            % Closed condition - with Depressed covariate
            tbl_closed = table(y_closed, sex, age, depressed, 'VariableNames', {'MVGC','Sex','Age','Depressed'});
            mdl_closed = fitlm(tbl_closed, 'MVGC ~ Sex + Age + Depressed', 'RobustOpts', 'on');
            coef_closed = mdl_closed.Coefficients;
            mvgc_closed_rlm_pvals_included(r1, r2) = coef_closed.pValue(4); % Depressed p-value
            mvgc_closed_rlm_tvals_included(r1, r2) = coef_closed.tStat(4);  % Depressed t-value

        end
    end

    %% Correct p-values and mask t-values based on significance

    fprintf('\n*** CORRECTING P-VALUES AND MASKING T-VALUES *** \n\n');

    mvgc_open_pvals_corrected = NaN(n_regions, n_regions);
    mvgc_closed_pvals_corrected = NaN(n_regions, n_regions);

    % OPTION 1 - Apply FDR correction using fdr_bh
    % [~, ~, ~, open_fdr_pvals] = fdr_bh(open_p_list, 0.05, 'pdep', 'yes');
    % [~, ~, ~, closed_fdr_pvals] = fdr_bh(closed_p_list, 0.05, 'pdep', 'yes');

    % OPTION 2 - Apply FDR correction using mafdr (same results as fdr_bh(open_p_list, 0.05, 'pdep', 'yes'))
    open_fdr_pvals = mafdr(open_p_list, 'BHFDR', true);
    closed_fdr_pvals = mafdr(closed_p_list, 'BHFDR', true);

    for i = 1:(n_regions * (n_regions - 1))
        r1 = index_map(i,1);
        r2 = index_map(i,2);
        mvgc_open_pvals_corrected(r1, r2) = open_fdr_pvals(i);
        mvgc_closed_pvals_corrected(r1, r2) = closed_fdr_pvals(i);
    end

    % Mask t-values based on corrected significance (e.g., α = 0.05)

    alpha = 0.05;
    mvgc_open_tvals_masked = mvgc_open_tvals;
    mvgc_open_tvals_masked(mvgc_open_pvals_corrected > alpha) = 0;

    mvgc_closed_tvals_masked = mvgc_closed_tvals;
    mvgc_closed_tvals_masked(mvgc_closed_pvals_corrected > alpha) = 0;

    %% Display p-value & t-value results

    disp('mvgc_open_pvals:');
    disp(mvgc_open_pvals);
    disp('mvgc_open_pvals_corrected:');
    disp(mvgc_open_pvals_corrected);
    disp('mvgc_open_rlm_pvals_included:');
    disp(mvgc_open_rlm_pvals_included);
    disp('mvgc_open_rlm_pvals_residuals:');
    disp(mvgc_open_rlm_pvals_residuals);

    disp('mvgc_closed_pvals:');
    disp(mvgc_closed_pvals);
    disp('mvgc_closed_pvals_corrected:');
    disp(mvgc_closed_pvals_corrected);
    disp('mvgc_closed_rlm_pvals_included:');
    disp(mvgc_closed_rlm_pvals_included);
    disp('mvgc_closed_rlm_pvals_residuals:');
    disp(mvgc_closed_rlm_pvals_residuals);

    disp('mvgc_open_tvals:');
    disp(mvgc_open_tvals);
    disp('mvgc_open_tvals_masked:');
    disp(mvgc_open_tvals_masked);
    disp('mvgc_open_rlm_tvals_included:');
    disp(mvgc_open_rlm_tvals_included);
    disp('mvgc_open_rlm_tvals_residuals:');
    disp(mvgc_open_rlm_tvals_residuals);

    disp('mvgc_closed_tvals:');
    disp(mvgc_closed_tvals);
    disp('mvgc_closed_tvals_masked:');
    disp(mvgc_closed_tvals_masked);
    disp('mvgc_closed_rlm_tvals_included:');
    disp(mvgc_closed_rlm_tvals_included);
    disp('mvgc_closed_rlm_tvals_residuals:');
    disp(mvgc_closed_rlm_tvals_residuals);

    %% Plotting matrices
    
    % Plot t-test p-values and t-values

    % plot_mvgc(...
    %     {mvgc_open_pvals_corrected, mvgc_closed_pvals_corrected, mvgc_open_tvals_masked, mvgc_closed_tvals_masked}, ...
    %     {'Open Corrected P-Values', 'Closed Corrected P-Values', 'Open Masked T-Values', 'Closed Masked T-Values'}, ...
    %     sprintf('MVGC P-Values and T-Values (%s)', strrep(mvgc_format, '_', ' ')));

    % Plot RLM p-values and t-values

    plot_mvgc(...
        {mvgc_open_pvals, mvgc_open_rlm_pvals_included, mvgc_open_rlm_pvals_residuals, ...
         mvgc_closed_pvals, mvgc_closed_rlm_pvals_included, mvgc_closed_rlm_pvals_residuals}, ...
        {'Open (ttest2)', 'Open RLM (Sex+Age+Depression)', 'Open RLM (Residuals)', ...
         'Closed (ttest2)', 'Closed RLM (Sex+Age+Depression)', 'Closed RLM (Residuals)'}, ...
        sprintf('MVGC P-Values (%s)', strrep(mvgc_format, '_', ' ')));

    % Plot differences between MVGC matrices

    mvgc_open_avg_diff = mvgc_open_healthy_avg - mvgc_open_depressed_avg;
    mvgc_closed_avg_diff = mvgc_closed_healthy_avg - mvgc_closed_depressed_avg;
    mvgc_eyes_diff = mvgc_open_avg - mvgc_closed_avg;
    mvgc_health_diff = mvgc_healthy_avg - mvgc_depressed_avg;

    % plot_mvgc(...
    %     {mvgc_open_avg_diff, mvgc_closed_avg_diff, mvgc_eyes_diff, mvgc_health_diff}, ...
    %     {'Open Healthy vs Open Depressed', 'Closed Healthy vs Closed Depressed', 'Open vs Closed', 'Healthy vs Depressed'}, ...
    %     sprintf('MVGC Differences (%s)', strrep(mvgc_format, '_', ' ')));

    % Plot Function

    function plot_mvgc(mvgc_matrices, titles, super_title, signif_threshs)
        if ~plot_data
            return;
        end

        if nargin < 4
            signif_threshs = [];
        end

        % Find global color limits
        num_matrices = numel(mvgc_matrices);
        mat_size = numel(mvgc_matrices{1});
        all_vals = zeros(num_matrices * mat_size, 1);
        for i = 1:num_matrices
            idx_start = (i-1)*mat_size + 1;
            idx_end = i*mat_size;
            all_vals(idx_start:idx_end) = mvgc_matrices{i}(:);
        end
        clim = [min(all_vals), max(all_vals)];

        figure;
        for i = 1:numel(mvgc_matrices)
            if numel(mvgc_matrices) == 2
                subplot(1, 2, i)
            end
            if numel(mvgc_matrices) == 4
                subplot(2, 2, i)
            end
            if numel(mvgc_matrices) == 6
                subplot(2, 3, i)
            end
            imagesc(mvgc_matrices{i}, clim);
            colorbar;
            title(titles{i});
            set(gca, 'XTick', 1:length(region_names), 'XTickLabel', region_names, ...
                     'YTick', 1:length(region_names), 'YTickLabel', region_names);
            xlabel('To Region');
            ylabel('From Region');
            axis square;

            % Write the values in the center of each cell
            for r = 1:length(region_names)
                for c = 1:length(region_names)
                    text(c, r, sprintf('%.3f', mvgc_matrices{i}(r, c)), ...
                         'Color', 'k', 'FontSize', 10, 'HorizontalAlignment', 'center');
                end
            end

            % Mark significant differences with stars
            if ~isempty(signif_threshs)
                hold on;
                for r = 1:length(region_names)
                    for c = 1:length(region_names)
                        if abs(mvgc_matrices{i}(r, c)) >= signif_threshs(i)  % Example threshold
                            text(c, r, '*', 'Color', 'k', 'FontSize', 14, 'HorizontalAlignment', 'center');
                        end
                    end
                end
            end
        end
        sgtitle(super_title);

        if (save_plots)
            if ~exist(fullfile('output', 'images'), 'dir')
                mkdir(fullfile('output', 'images'));
            end
            filename = sprintf('%s - %s.png', super_title, datetime('now', 'Format', 'yyyy_MM_dd'));
            saveas(gcf, fullfile('output', 'images', filename));
            % close(gcf);  % Close the figure after saving
        end
        
    end

    %% Save results into .mat files

    if save_mvgc_file
        fprintf('\n*** SAVING MVGC RESULTS TO %s *** \n', filename);

        save(fullfile(filename), ...
         'mvgc_open', 'mvgc_closed', ...
         'mvgc_open_healthy', 'mvgc_closed_healthy', ...
         'mvgc_open_depressed', 'mvgc_closed_depressed', ...
         'mvgc_open_healthy_avg', 'mvgc_closed_healthy_avg', ...
         'mvgc_open_depressed_avg', 'mvgc_closed_depressed_avg', ...
         'mvgc_open_avg', 'mvgc_closed_avg', ...
         'mvgc_healthy_avg', 'mvgc_depressed_avg', ...
         'mvgc_open_pvals', 'mvgc_closed_pvals', ...
         'mvgc_open_pvals_corrected', 'mvgc_closed_pvals_corrected', ...
         'mvgc_open_tvals', 'mvgc_closed_tvals', ...
         'mvgc_open_tvals_masked', 'mvgc_closed_tvals_masked', ...
         'ss_info_open', 'ss_info_closed');

        fprintf('\n*** FILES SAVED *** \n');
    else
        fprintf('\n*** NOT SAVING MVGC RESULTS *** \n');
    end

end