function [mvgc_open, mvgc_closed, ss_info_open, ss_info_closed] = EEGtoMVGC2()

    %% Load .mat and .m files

    dataDir = fullfile('Depression_Study', 'export_mat'); % TEST directory
    files = dir(fullfile(dataDir, '*.mat')); % Get all .mat files in the directory

    addpath('Source_Reconstruction'); % Add path to Source_Reconstruction folder
    addpath('EntRate'); % Add path to EntRate folder

    % Load MVGC toolbox - silently adding MVGC to current path
    % p = mfilename('fullpath');
    % addpath(strrep(p, 'EEGtoMVGC', 'EntRate/private/mvgc_v2.0'));
    % evalc('mvgc_startup;');

    %% Initialize variables

    mvgc_open = struct();
    mvgc_closed = struct();

    pwcgc_open = struct();
    pwcgc_closed = struct();

    spwcgc_open = struct();
    spwcgc_closed = struct();

    bands_open = struct();
    bands_closed = struct();

    ss_info_open = struct();
    ss_info_closed = struct();

    region_names = {'Frontal', 'Occipital', 'Parietal', 'Sensorimotor', 'Temporal'};

    frontal = [3:16, 19:24];
    occipital = [25:36];
    parietal = [39:50];
    sensorimotor = [1,2,17,18,37,38];
    temporal = [51:60];

    LOAD_MVGC_FILE = false;

    %% Perform source reconstruction, rotation and CSER calculation for each subject

    if ~LOAD_MVGC_FILE

        % Loop through each file
        for i = 55:length(files)
        % for i = 1:5
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
            [pwcgc_values_open, spwcgc_values_open, bands_values_open, info_open] = calculateMVGC(subj_open);
            [pwcgc_values_closed, spwcgc_values_closed, bands_values_closed, info_closed] = calculateMVGC(subj_closed);

            pwcgc_open.(key) = pwcgc_values_open;
            pwcgc_closed.(key) = pwcgc_values_closed;

            spwcgc_open.(key) = spwcgc_values_open;
            spwcgc_closed.(key) = spwcgc_values_closed;

            bands_open.(key) = bands_values_open;
            bands_closed.(key) = bands_values_closed;

            ss_info_open.(key) = info_open;
            ss_info_closed.(key) = info_closed;

            fprintf('\n*** RESULTS SAVED: %s *** \n', files(i).name);
        end
    
    end

    mvgc_open = pwcgc_open;  % Store pairwise conditional GC
    mvgc_closed = pwcgc_closed;  % Store pairwise conditional GC

    %% Function to calculate MVGC for a given subject's data

    function [pwcgc_values, spwcgc_values, bands_values, info] = calculateMVGC(subj_data)

        % Reconstruct 5 regions out of 60 sources
        frontal_data = subj_data(frontal, :, :);  % 20 srcs * 2000 samples * trials
        occipital_data = subj_data(occipital, :, :);  % 12 srcs * 2000 samples * trials
        parietal_data = subj_data(parietal, :, :);  % 12 srcs * 2000 samples * trials
        sensorimotor_data = subj_data(sensorimotor, :, :);  % 6 srcs * 2000 samples * trials
        temporal_data = subj_data(temporal, :, :);  % 10 srcs * 2000 samples * trials

        num_regions = 5;
        [~, num_samples, num_trials] = size(frontal_data);
        regional_eeg = zeros(num_regions, num_samples, num_trials);  % [regions x samples x trials]

        regional_eeg(1, :, :) = squeeze(mean(frontal_data, 1));  % [samples x trials]
        regional_eeg(2, :, :) = squeeze(mean(occipital_data, 1));
        regional_eeg(3, :, :) = squeeze(mean(parietal_data, 1));
        regional_eeg(4, :, :) = squeeze(mean(sensorimotor_data, 1));
        regional_eeg(5, :, :) = squeeze(mean(temporal_data, 1));

        %% Calculate MVGC between regions

        % Downsampling the data
        fprintf('\nDownsampling data to 250 Hz...\n');
        disp(size(subj_data));
        subj_data = subj_data(:, 1:2:end, :);  % Downsample to 250 Hz (from 500 Hz)
        disp(size(subj_data));

        % X = regional_eeg;  % [regions x samples x trials]
        X = subj_data;
        X = subj_data(:, :, :);  % Use only frontal region for now [frontal, occipital]
        disp(size(X));

        % Z-score each region across time for each trial
        for trial = 1:size(X, 3)
            for region = 1:size(X, 1)
                X(region, :, trial) = zscore(X(region, :, trial));
                % TODO: catch if error, and add noise
                % X(region, :, trial) = zscore(X(region, :, trial)) + randn(1, size(X, 2)) * 0.001;  % Add small noise to avoid numerical issues  
            end
        end

        disp(eig(cov(X(:, :, 1)')));  % Display eigenvalues of covariance matrix for open condition


        % VAR model order estimation
        varmomax = 10;  % has to be less than nobs
        varmosel = 'AIC';
        fprintf('\nEstimating VAR model order (max order = %d)...\n', varmomax);
        [varmoaic, varmobic, varmohqc, varmolrt] = tsdata_to_varmo(X, varmomax, 'LWR', [], []);
        varmo = moselect(sprintf('VAR model order selection (max = %d)', varmomax), varmosel, 'AIC', varmoaic, 'BIC', varmobic, 'HQC', varmohqc, 'LRT', varmolrt);
        fprintf('Selected VAR model order (AIC): %d\n', varmo);

        % State-space model order estimation using SVC (faster)
        fprintf('Estimating state-space model order using SVC...\n');
        [ssmo_svc, ssmo_max] = tsdata_to_sssvc(X, 2*varmo, [], []);
        % [ssmo_svc, ssmo_max] = tsdata_to_ssmo(X, 2*varmo);  % Incredibly slow
        ssmo = ssmo_svc;
        fprintf('Selected state-space model order (SVC): %d (max = %d)\n', ssmo, ssmo_max);

        % Estimate SS model
        fprintf('Estimating state-space model parameters...\n');
        [A, C, K, V] = tsdata_to_ss(X, 2*varmo, ssmo);
        info = ss_info(A, C, K, V, 1);
        % fprintf('\nState-space model information:\n');
        % fprintf('  Observables: %d\n', info.observ);
        % fprintf('  Model order: %d\n', info.morder);
        % fprintf('  Rho A: %.4f\n', info.rhoA);
        % fprintf('  Rho B: %.4f\n', info.rhoB);
        % fprintf('  Autocovariance decay: %d\n', info.acdec);
        fprintf('  Sigma SPD: %d\n', info.sigspd);
        % fprintf('  Multi-information: %.4f\n', info.mii);
        % fprintf('  Multi-information (uniform): %.4f\n', info.mmii);
        if info.error
            fprintf('State-space model estimation encountered errors.\n');
            % fprintf('Error code: %d\n', info.error);
            % return;
        else
            fprintf('State-space model estimation successful (no errors).\n\n');
        end

        % Compute time-domain pairwise conditional GC
        fprintf('Computing time-domain pairwise conditional MVGC... (1/2)\n');
        % pwcgc_values = ss_to_pwcgc(A, C, K, V);
        disp(frontal);
        disp(occipital);
        from_indices = 1:length(frontal);
        to_indices = length(frontal)+1:length(frontal)+length(occipital);
        disp(from_indices);
        disp(to_indices);
        pwcgc_values = ss_to_mvgc(A, C, K, V, to_indices, from_indices);

        % Compute spectral pairwise conditional GC
        % fprintf('Computing spectral pairwise conditional GC... (2/2)\n');
        % fres = 200;  % frequency resolution
        % spwcgc_values = ss_to_spwcgc(A, C, K, V, fres);  % [nvars x nvars x fres]
        % bands_values = compute_band_gc(spwcgc_values, 200, fres);
        spwcgc_values = struct();
        bands_values = struct();
        
        fprintf('\nMVGC computation complete.\n');

        % Validate result
        if isbad(pwcgc_values, false)
            error('pwcgc_values estimation failed');
        end

        % if isbad(spwcgc_values, false)
        %     error('spwcgc_values estimation failed');
        % end

        % Print mvgc contents
        fprintf('\npwcgc_values contents:\n\n');
        disp(pwcgc_values);

        % fprintf('\nbands_values.delta contents:\n\n');
        % disp(bands_values.delta);

        % fprintf('\nbands_values.theta contents:\n\n');
        % disp(bands_values.theta);

        % fprintf('\nbands_values.alpha contents:\n\n');
        % disp(bands_values.alpha);

        % fprintf('\nbands_values.beta contents:\n\n');
        % disp(bands_values.beta);

        % fprintf('\nbands_values.gamma contents:\n\n');
        % disp(bands_values.gamma);

    end

    function band_gc = compute_band_gc(gc_spectral, fs, fres)

        % Frequency vector (0 to Nyquist)
        nyq = fs / 2;
        freqs = linspace(0, nyq, fres + 1);

        % Define frequency bands (in Hz)
        bands = struct( ...
            'delta', [0.5 4], ...
            'theta', [4 8], ...
            'alpha', [8 13], ...
            'beta',  [13 30], ...
            'gamma', [30 80] ...
        );

        band_names = fieldnames(bands);

        % Loop through each band and average across frequency indices
        for b = 1:numel(band_names)
            band = band_names{b};
            range = bands.(band);
            idx = freqs >= range(1) & freqs <= range(2);
            band_gc.(band) = mean(gc_spectral(:, :, idx), 3);
        end
    end

    %% Load MVGC results from mvgc_values.mat

    if LOAD_MVGC_FILE

        fprintf('\n*** LOADING MVGC RESULTS *** \n');
        if exist('mvgc_values.mat', 'file')
            load('mvgc_values.mat', 'mvgc_open', 'mvgc_closed');
            fprintf('MVGC results loaded successfully.\n');
        else
            fprintf('No MVGC results found. Please run EEGtoMVGC first.\n');
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
    plot_mvgc_matrix_grid(...
        {mvgc_open_healthy_avg, mvgc_open_depressed_avg, mvgc_closed_healthy_avg, mvgc_closed_depressed_avg}, ...
        region_names, ...
        {'Open Healthy', 'Open Depressed', 'Closed Healthy', 'Closed Depressed'}, ...
        'Average MVGC Matrices');

    plot_mvgc_matrix_grid(...
        {mvgc_open_avg, mvgc_closed_avg, mvgc_healthy_avg, mvgc_depressed_avg}, ...
        region_names, ...
        {'Open', 'Closed', 'Healthy', 'Depressed'}, ...
        'Average MVGC Matrices');

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

    save(fullfile('mvgc_values.mat'), ...
         'mvgc_open', 'mvgc_closed', ...
         'mvgc_open_healthy', 'mvgc_open_depressed', ...
         'mvgc_closed_healthy', 'mvgc_closed_depressed', ...
         'mvgc_open_healthy_avg', 'mvgc_open_depressed_avg', ...
         'mvgc_closed_healthy_avg', 'mvgc_closed_depressed_avg');
    
    save(fullfile('spwcgc_values.mat'), ...
         'spwcgc_open', 'spwcgc_closed', ...
         'bands_open', 'bands_closed');

    fprintf('\n*** FILES SAVED *** \n');

end