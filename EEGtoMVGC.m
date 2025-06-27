function [mvgc_open, mvgc_closed, mvgc_open_avg, mvgc_closed_avg] = EEGtoMVGC()

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

    region_names = {'Frontal', 'Occipital', 'Parietal', 'Sensorimotor', 'Temporal'};

    frontal = [3:16, 19:24];
    occipital = [25:36];
    parietal = [39:50];
    sensorimotor = [1,2,17,18,37,38];
    temporal = [51:60];

    LOAD_MVGC_FILE = true;

    %% Perform source reconstruction, rotation and CSER calculation for each subject

    if ~LOAD_MVGC_FILE

        % Loop through each file
        for i = 1:length(files)
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
            mvgc_open.(key) = calculateMVGC(subj_open);
            mvgc_closed.(key) = calculateMVGC(subj_closed);

            fprintf('\n*** RESULTS SAVED: %s *** \n', files(i).name);
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

    % Initialize matrices to accumulate MVGC values for healthy and depressed groups
    mvgc_open_healthy_avg = zeros(length(region_names), length(region_names));
    mvgc_open_depressed_avg = zeros(length(region_names), length(region_names));
    mvgc_closed_healthy_avg = zeros(length(region_names), length(region_names));
    mvgc_closed_depressed_avg = zeros(length(region_names), length(region_names));

    % Count subjects in each group
    open_healthy_count = 0;
    open_depressed_count = 0;
    closed_healthy_count = 0;
    closed_depressed_count = 0;

    % Loop through each subject's MVGC results for open
    for i = 1:length(files)
        key = sprintf('x%s', files(i).name(1:3));
        if isfield(mvgc_open_healthy, key)
            mvgc_open_healthy_avg = mvgc_open_healthy_avg + mvgc_open_healthy.(key);
            open_healthy_count = open_healthy_count + 1;
        end
        if isfield(mvgc_open_depressed, key)
            mvgc_open_depressed_avg = mvgc_open_depressed_avg + mvgc_open_depressed.(key);
            open_depressed_count = open_depressed_count + 1;
        end
        if isfield(mvgc_closed_healthy, key)
            mvgc_closed_healthy_avg = mvgc_closed_healthy_avg + mvgc_closed_healthy.(key);
            closed_healthy_count = closed_healthy_count + 1;
        end
        if isfield(mvgc_closed_depressed, key)
            mvgc_closed_depressed_avg = mvgc_closed_depressed_avg + mvgc_closed_depressed.(key);
            closed_depressed_count = closed_depressed_count + 1;
        end
    end

    % Average the MVGC values for each group (avoid division by zero)
    if open_healthy_count > 0
        mvgc_open_healthy_avg = mvgc_open_healthy_avg / open_healthy_count;
    end
    if open_depressed_count > 0
        mvgc_open_depressed_avg = mvgc_open_depressed_avg / open_depressed_count;
    end
    if closed_healthy_count > 0
        mvgc_closed_healthy_avg = mvgc_closed_healthy_avg / closed_healthy_count;
    end
    if closed_depressed_count > 0
        mvgc_closed_depressed_avg = mvgc_closed_depressed_avg / closed_depressed_count;
    end

    % Display average MVGC matrices
    fprintf('\n*** AVERAGE MVGC OPEN HEALTHY *** \n');
    disp(mvgc_open_healthy_avg);
    fprintf('\n*** AVERAGE MVGC OPEN DEPRESSED *** \n');
    disp(mvgc_open_depressed_avg);
    fprintf('\n*** AVERAGE MVGC CLOSED HEALTHY *** \n');
    disp(mvgc_closed_healthy_avg);
    fprintf('\n*** AVERAGE MVGC CLOSED DEPRESSED *** \n');
    disp(mvgc_closed_depressed_avg);

    %% Average MVGC across open/closed subjects and healthy/depressed subjects separately

    fprintf('\n*** AVERAGING MVGC ACROSS OPEN/CLOSED AND HEALTHY/DEPRESSED SUBJECTS SEPARATELY *** \n');
    mvgc_open_avg = (mvgc_open_healthy_avg + mvgc_open_depressed_avg) / 2;
    mvgc_closed_avg = (mvgc_closed_healthy_avg + mvgc_closed_depressed_avg) / 2;
    mvgc_healthy_avg = (mvgc_open_healthy_avg + mvgc_closed_healthy_avg) / 2;
    mvgc_depressed_avg = (mvgc_open_depressed_avg + mvgc_closed_depressed_avg) / 2;

    % Display average MVGC matrices
    fprintf('\n*** AVERAGE MVGC OPEN *** \n');
    disp(mvgc_open_avg);
    fprintf('\n*** AVERAGE MVGC CLOSED *** \n');
    disp(mvgc_closed_avg);
    fprintf('\n*** AVERAGE MVGC HEALTHY *** \n');
    disp(mvgc_healthy_avg);
    fprintf('\n*** AVERAGE MVGC DEPRESSED *** \n');
    disp(mvgc_depressed_avg);
    

    %% Plot average MVGC matrices across open/closed healthy/depressed subjects

    % Find global min and max for color scaling
    all_mvgc = [mvgc_open_healthy_avg(:); mvgc_open_depressed_avg(:); ...
                mvgc_closed_healthy_avg(:); mvgc_closed_depressed_avg(:)];
    clim = [min(all_mvgc), max(all_mvgc)];

    figure;

    subplot(2, 2, 1);
    imagesc(mvgc_open_healthy_avg, clim);
    colorbar;
    title('Average MVGC - Open Healthy');
    set(gca, 'XTick', 1:length(region_names), 'XTickLabel', region_names, ...
             'YTick', 1:length(region_names), 'YTickLabel', region_names);
    xlabel('To Region');
    ylabel('From Region');
    axis square;

    subplot(2, 2, 2);
    imagesc(mvgc_open_depressed_avg, clim);
    colorbar;
    title('Average MVGC - Open Depressed');
    set(gca, 'XTick', 1:length(region_names), 'XTickLabel', region_names, ...
             'YTick', 1:length(region_names), 'YTickLabel', region_names);
    xlabel('To Region');
    ylabel('From Region');
    axis square;

    subplot(2, 2, 3);
    imagesc(mvgc_closed_healthy_avg, clim);
    colorbar;
    title('Average MVGC - Closed Healthy');
    set(gca, 'XTick', 1:length(region_names), 'XTickLabel', region_names, ...
             'YTick', 1:length(region_names), 'YTickLabel', region_names);
    xlabel('To Region');
    ylabel('From Region');
    axis square;

    subplot(2, 2, 4);
    imagesc(mvgc_closed_depressed_avg, clim);
    colorbar;
    title('Average MVGC - Closed Depressed');
    set(gca, 'XTick', 1:length(region_names), 'XTickLabel', region_names, ...
             'YTick', 1:length(region_names), 'YTickLabel', region_names);
    xlabel('To Region');
    ylabel('From Region');
    axis square;

    sgtitle('Average MVGC Matrices');

    %% Plot average MVGC matrices across open/closed and healthy/depressed subjects separately

    % Find global min and max for color scaling
    all_mvgc = [mvgc_open_avg(:); mvgc_closed_avg(:); ...
                mvgc_healthy_avg(:); mvgc_depressed_avg(:)];
    clim = [min(all_mvgc), max(all_mvgc)];

    figure;

    subplot(2, 2, 1);
    imagesc(mvgc_open_avg, clim);
    colorbar;
    title('Average MVGC - Open');
    set(gca, 'XTick', 1:length(region_names), 'XTickLabel', region_names, ...
             'YTick', 1:length(region_names), 'YTickLabel', region_names);
    xlabel('To Region');
    ylabel('From Region');
    axis square;

    subplot(2, 2, 2);
    imagesc(mvgc_closed_avg, clim);
    colorbar;
    title('Average MVGC - Closed');
    set(gca, 'XTick', 1:length(region_names), 'XTickLabel', region_names, ...
             'YTick', 1:length(region_names), 'YTickLabel', region_names);
    xlabel('To Region');
    ylabel('From Region');
    axis square;

    subplot(2, 2, 3);
    imagesc(mvgc_healthy_avg, clim);
    colorbar;
    title('Average MVGC - Healthy');
    set(gca, 'XTick', 1:length(region_names), 'XTickLabel', region_names, ...
             'YTick', 1:length(region_names), 'YTickLabel', region_names);
    xlabel('To Region');
    ylabel('From Region');
    axis square;

    subplot(2, 2, 4);
    imagesc(mvgc_depressed_avg, clim);
    colorbar;
    title('Average MVGC - Depressed');
    set(gca, 'XTick', 1:length(region_names), 'XTickLabel', region_names, ...
             'YTick', 1:length(region_names), 'YTickLabel', region_names);
    xlabel('To Region');
    ylabel('From Region');
    axis square;

    sgtitle('Average MVGC Matrices');

    %% Save results into .mat files

    save(fullfile('mvgc_values.mat'), ...
         'mvgc_open', 'mvgc_closed', ...
         'mvgc_open_healthy', 'mvgc_open_depressed', ...
         'mvgc_closed_healthy', 'mvgc_closed_depressed', ...
         'mvgc_open_healthy_avg', 'mvgc_open_depressed_avg', ...
         'mvgc_closed_healthy_avg', 'mvgc_closed_depressed_avg');

    fprintf('\n*** FILES SAVED *** \n');

    %% Function to calculate MVGC for a given subject's data

    function [mvgc_values] = calculateMVGC(subj_data)

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

        % %% Calculate MVGC between regions (needs a further breakdown)

        % % Select VAR model order
        % [varmoaic, varmobic, varmohqc, varmolrt] = tsdata_to_varmo(y, varmomax, 'LWR', [], false);

        % % Select SS model order
        % if varmohqc < 1
        %     % Force the algorithm to fit a SS model even if varmo gives up
        %     varmohqc = 1;
        %     ssmo = 2;
        % else
        %     [ssmo, ~] = tsdata_to_sssvc(y, 2*varmohqc, [], []);
        % end

        % % Fit SS model
        % [A, C, K, V] = tsdata_to_ss(y, 2*varmohqc, ssmo);
        % info = ss_info(A, C, K, V, 0);
        % F = ss_to_mvgc(A,C,K,V,x,y);

        X = regional_eeg;  % [regions x samples x trials]

        % Set parameters
        nvars = size(X, 1);       % should be 5
        nobs = size(X, 2);        % should be 2000
        ntrials = size(X, 3);
        fs = 200;                 % Sampling rate

        % VAR model order estimation
        max_var_order = 20;       % reasonable maximum (has to be less than nobs)
        fprintf('\nEstimating VAR model order (max order = %d)...\n', max_var_order);
        [varmo_aic, ~, ~, ~] = tsdata_to_varmo(X, max_var_order, 'LWR');
        varmo = find(varmo_aic == min(varmo_aic), 1); % Select best AIC order
        fprintf('Selected VAR model order (AIC): %d\n', varmo);

        % State-space model order estimation using SVC (faster)
        fprintf('Estimating state-space model order using SVC...\n');
        [ssmo_svc, ~] = tsdata_to_sssvc(X, 2*varmo, [], []);
        ssmo = ssmo_svc;
        fprintf('Selected state-space model order (SVC): %d\n', ssmo);

        % Estimate SS model
        fprintf('Estimating state-space model parameters...\n');
        [A, C, K, V] = tsdata_to_ss(X, 2*varmo, ssmo);

        % Compute time-domain pairwise conditional GC
        fprintf('Computing time-domain pairwise conditional MVGC...\n');
        mvgc_values = ss_to_pwcgc(A, C, K, V);
        fprintf('\nMVGC computation complete.\n');

        % Validate result
        if isbad(mvgc_values, false)
            error('GC estimation failed');
        end

        % Print mvgc_values's size and contents
        fprintf('\nmvgc_values contents:\n\n');
        disp(mvgc_values);

    end

end