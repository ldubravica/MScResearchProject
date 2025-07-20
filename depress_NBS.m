%% Depression tests using MVGC data

% Set up and load MVGC data from EEGtoMVGCOptions
addpath('./NBS_test/'); % Add parent directory to access EEGtoMVGCOptions

% Load MVGC data 5x5 connectivity matrices
fprintf('Loading MVGC data...\n');
[mvgc_open, mvgc_closed, ~, ~] = EEGtoMVGCOptions('brain_source', true); % Load from saved file

% Load demographic data
excelFile = fullfile('Depression_Study', 'depression_data', 'Data_4_Import_REST.xlsx');
T = readtable(excelFile);

% Drop rows with no associated file or bad data (already done in EEGtoMVGCOptions)
T(T.id == 544 | T.id == 571 | T.id == 572 | T.id == 527 | T.id == 535, :) = [];

% Get subject IDs from MVGC data
subject_names = fieldnames(mvgc_open);
subject_ids = cellfun(@(x) str2double(x(2:end)), subject_names); % Remove 'x' prefix and convert to number

% Find matching rows in demographic table
[~, T_idx] = ismember(subject_ids, T.id);
T_matched = T(T_idx, :);

nb_subjects = length(subject_names);
assert(height(T_matched) == nb_subjects, 'Mismatch between MVGC data and demographic data');

% Create binary depression labels (MDD <= 2 indicates depression)
is_depressed = double(T_matched.MDD <= 2);

% Convert MVGC structures to cell arrays for NBS
mvgc_open_nets = cell([nb_subjects, 1]);
mvgc_closed_nets = cell([nb_subjects, 1]);
for i = 1:nb_subjects
    mvgc_open_nets{i} = mvgc_open.(subject_names{i});
    mvgc_closed_nets{i} = mvgc_closed.(subject_names{i});
end

fprintf('*** RUNNING DEPRESSION NBS TESTS ***\n');

% Initialize UI structure for NBS analysis
% Default UI structure based on NBS requirements
UI = struct();
UI.method.ui = 'Connectivity';
UI.test.ui = 't-test';
UI.perms.ui = '5000';
UI.thresh.ui = '3.1';
UI.alpha.ui = '0.05';
UI.size.ui = 'Extent';
UI.exchange.ui = '';

% Set up design matrix: [intercept, depressed, age, age^2, sex]
depress_UI_open = UI;
depress_UI_open.design.ui = [ones(nb_subjects, 1), is_depressed, T_matched.age, T_matched.age.^2, T_matched.sex];
depress_UI_open.matrices.ui = mvgc_open_nets;

depress_UI_closed = UI;
depress_UI_closed.design.ui = [ones(nb_subjects, 1), is_depressed, T_matched.age, T_matched.age.^2, T_matched.sex];
depress_UI_closed.matrices.ui = mvgc_closed_nets;

% Define ROI names for region connectivity matrices  
ROI = {'Frontal', 'Occipital', 'Parietal', 'Sensorimotor', 'Temporal'};

% Test for depression-related changes in open eyes condition
fprintf('\nTesting for depression-related increases in open eyes condition...\n');
depress_UI_open.contrast.ui = [0 +1 0 0 0]; % depressed > healthy
depress_NBS_open_pos = CompositeNBStest({depress_UI_open});

fprintf('\nTesting for depression-related decreases in open eyes condition...\n');
depress_UI_open.contrast.ui = [0 -1 0 0 0]; % healthy > depressed
depress_NBS_open_neg = CompositeNBStest({depress_UI_open});

% Test for depression-related changes in closed eyes condition
fprintf('\nTesting for depression-related increases in closed eyes condition...\n');
depress_UI_closed.contrast.ui = [0 +1 0 0 0]; % depressed > healthy
depress_NBS_closed_pos = CompositeNBStest({depress_UI_closed});

fprintf('\nTesting for depression-related decreases in closed eyes condition...\n');
depress_UI_closed.contrast.ui = [0 -1 0 0 0]; % healthy > depressed
depress_NBS_closed_neg = CompositeNBStest({depress_UI_closed});

% Test for eye condition differences (closed > open)
fprintf('\nTesting for closed > open eyes differences...\n');
depress_UI_diff = UI;
depress_UI_diff.design.ui = [ones(nb_subjects, 1), is_depressed, T_matched.age, T_matched.age.^2, T_matched.sex];
% Create difference matrices (closed - open)
mvgc_diff_nets = cell([nb_subjects, 1]);
for i = 1:nb_subjects
    mvgc_diff_nets{i} = mvgc_closed_nets{i} - mvgc_open_nets{i};
end
depress_UI_diff.matrices.ui = mvgc_diff_nets;
depress_UI_diff.contrast.ui = [0 +1 0 0 0]; % depressed difference > healthy difference
depress_NBS_diff = CompositeNBStest({depress_UI_diff});

% Save significant positive clusters for open eyes condition
fprintf('\nSaving open eyes condition results...\n');
T_open = array2table(zeros([0,4]), 'VariableNames', {'From', 'To', 'Tvalue', 'Significant'});
for i=1:5
  for j=1:5
    if i == j, continue; end
    T_open = [T_open; {ROI{i}, ROI{j}, depress_NBS_open_pos.test_stat{1}(i,j), depress_NBS_open_pos.con_mat{1}(i,j)}];
  end
end
writetable(T_open, 'depress_network_results_open.csv');

% Save significant positive clusters for closed eyes condition
fprintf('\nSaving closed eyes condition results...\n');
T_closed = array2table(zeros([0,4]), 'VariableNames', {'From', 'To', 'Tvalue', 'Significant'});
for i=1:5
  for j=1:5
    if i == j, continue; end
    T_closed = [T_closed; {ROI{i}, ROI{j}, depress_NBS_closed_pos.test_stat{1}(i,j), depress_NBS_closed_pos.con_mat{1}(i,j)}];
  end
end
writetable(T_closed, 'depress_network_results_closed.csv');

% Save condition difference results
fprintf('\nSaving condition difference results...\n');
T_diff = array2table(zeros([0,4]), 'VariableNames', {'From', 'To', 'Tvalue', 'Significant'});
for i=1:5
  for j=1:5
    if i == j, continue; end
    T_diff = [T_diff; {ROI{i}, ROI{j}, depress_NBS_diff.test_stat{1}(i,j), depress_NBS_diff.con_mat{1}(i,j)}];
  end
end
writetable(T_diff, 'depress_network_results_diff.csv');

% Display summary results
fprintf('\n========== RESULTS SUMMARY ==========\n');
fprintf('Open eyes condition (depressed > healthy):\n');
fprintf('  Number of significant clusters: %d\n', depress_NBS_open_pos.n);
if depress_NBS_open_pos.n > 0
    fprintf('  P-values: ');
    fprintf('%.4f ', depress_NBS_open_pos.pval);
    fprintf('\n');
end

fprintf('\nClosed eyes condition (depressed > healthy):\n');
fprintf('  Number of significant clusters: %d\n', depress_NBS_closed_pos.n);
if depress_NBS_closed_pos.n > 0
    fprintf('  P-values: ');
    fprintf('%.4f ', depress_NBS_closed_pos.pval);
    fprintf('\n');
end

fprintf('\nCondition difference (closed - open, depressed > healthy):\n');
fprintf('  Number of significant clusters: %d\n', depress_NBS_diff.n);
if depress_NBS_diff.n > 0
    fprintf('  P-values: ');
    fprintf('%.4f ', depress_NBS_diff.pval);
    fprintf('\n');
end

fprintf('\nNumber of negative clusters (should be 0):\n');
fprintf('  Open eyes: %d\n', depress_NBS_open_neg.n);
fprintf('  Closed eyes: %d\n', depress_NBS_closed_neg.n);
fprintf('=====================================\n');

