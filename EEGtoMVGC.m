function [mvgc_open, mvgc_closed] = EEGtoMVGC()

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
        subj_open = permute(source_ts_open, [2, 3, 1]);
        subj_closed = permute(source_ts_closed, [2, 3, 1]);

        % Calculate MVGC for each subject's data
        fprintf('\n*** CALCULATING MVGC %s - (%i/%i) *** \n', files(i).name(1:3), i, length(files));
        subj_open_mvgc = calculateMVGC(subj_open);
        subj_closed_mvgc = calculateMVGC(subj_closed);

        % Store MVGC per subject per region combination
        key = sprintf('x%s', files(i).name(1:3));

        mvgc_open.(key) = subj_open_mvgc;
        mvgc_closed.(key) = subj_closed_mvgc;

        fprintf('\n*** RESULTS SAVED: %s *** \n', files(i).name);
    end

    %% Save results into .mat files

    % save(fullfile('mvgc_values.mat'), 'mvgc_open', 'mvgc_closed');

    fprintf('\n*** FILES SAVED *** \n');

    %% Function to calculate MVGC for a given subject's data

    function [mvgc_values] = calculateMVGC(subj_data)

        mvgc_values = struct();

        % Reconstruct 5 regions out of 60 sources
        frontal_data = subj_data(frontal, :, :);  % 20 srcs * 2000 samples * trials
        occipital_data = subj_data(occipital, :, :);  % 12 srcs * 2000 samples * trials
        parietal_data = subj_data(parietal, :, :);  % 12 srcs * 2000 samples * trials
        sensorimotor_data = subj_data(sensorimotor, :, :);  % 6 srcs * 2000 samples * trials
        temporal_data = subj_data(temporal, :, :);  % 10 srcs * 2000 samples * trials

        % Calculate MVGC between regions (needs a further breakdown)

        % Select VAR model order
        [varmoaic, varmobic, varmohqc, varmolrt] = tsdata_to_varmo(y, varmomax, 'LWR', [], false);

        % Select SS model order
        if varmohqc < 1
            % Force the algorithm to fit a SS model even if varmo gives up
            varmohqc = 1;
            ssmo = 2;
        else
            [ssmo, ~] = tsdata_to_sssvc(y, 2*varmohqc, [], []);
        end

        % Fit SS model
        [A, C, K, V] = tsdata_to_ss(y, 2*varmohqc, ssmo);
        info = ss_info(A, C, K, V, 0);
        F = ss_to_mvgc(A,C,K,V,x,y);

    end

end