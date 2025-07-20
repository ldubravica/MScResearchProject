%% Load subject excel metadata

excelFile = fullfile('Depression_Study', 'depression_data', 'Data_4_Import_REST.xlsx');
T = readtable(excelFile);
T(T.id == 544 | T.id == 571 | T.id == 572 | T.id == 527 | T.id == 535, :) = [];
% Remove 570 due to no age data
T(T.id == 570, :) = [];

%% Load MVGC matrices

addpath('./NBS_test/')
addpath('./NBS_test/NBS1.2/');
[mvgc_open, mvgc_closed, ~, ~] = EEGtoMVGCOptions('brain_source', true);

% Remove 570 due to no age data
mvgc_open = rmfield(mvgc_open, 'x570');
mvgc_closed = rmfield(mvgc_closed, 'x570');

% Get subject IDs from MVGC data (and remove 'x' prefix and convert to number)
subject_names = fieldnames(mvgc_open);
subject_ids = cellfun(@(x) str2double(x(2:end)), subject_names);

% Find matching rows in demographic table
[~, T_idx] = ismember(subject_ids, T.id);
T_matched = T(T_idx, :);

nb_subjects = length(subject_names);
assert(height(T_matched) == nb_subjects, 'Mismatch between MVGC data and demographic data');

% Convert MVGC structures to cell arrays for NBS
mvgc_open_nets = cell([nb_subjects, 1]);
mvgc_closed_nets = cell([nb_subjects, 1]);
for i = 1:nb_subjects
    mvgc_open_nets{i} = mvgc_open.(subject_names{i});
    mvgc_closed_nets{i} = mvgc_closed.(subject_names{i});
end

%% TEMP DIAGNOSIS

is_depressed = double(T_matched.MDD <= 2);

% Add this diagnostic code after creating your design matrix
design_matrix = [is_depressed, T_matched.age, T_matched.sex];

% disp('Design matrix:');
% disp(design_matrix);

% Check for problems
fprintf('Design matrix diagnostics:\n');
fprintf('Matrix size: %d x %d\n', size(design_matrix));

% Check correlations between predictors
corr_matrix = corrcoef(design_matrix);
[row, col] = find(abs(corr_matrix) > 0.9 & corr_matrix ~= 1);
for i = 1:length(row)
    fprintf('High correlation (>0.9): Column %d vs %d: %.3f\n', row(i), col(i), corr_matrix(row(i), col(i)));
end

% Check for constant columns
for i = 1:size(design_matrix, 2)
    if std(design_matrix(:,i)) < 1e-10
        fprintf('Column %d is nearly constant\n', i);
    end
end

fprintf('\n');

%% Run depression NBS tests

is_depressed = double(T_matched.MDD <= 2);
age = T_matched.age;
sex = T_matched.sex - 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
% depress_UI_open.design.ui = [ones(nb_subjects, 1), is_depressed, T_matched.age, T_matched.age.^2, T_matched.sex];
depress_UI_open.design.ui = [is_depressed, T_matched.age, T_matched.sex];
depress_UI_open.matrices.ui = mvgc_open_nets;

depress_UI_closed = UI;
% depress_UI_closed.design.ui = [ones(nb_subjects, 1), is_depressed, T_matched.age, T_matched.age.^2, T_matched.sex];
depress_UI_closed.design.ui = [is_depressed, T_matched.age, T_matched.sex];
depress_UI_closed.matrices.ui = mvgc_closed_nets;

% Test for depression-related changes in open eyes condition
fprintf('\nTesting for depression-related increases in open eyes condition...\n');
% depress_UI_open.contrast.ui = [0 +1 0 0 0]; % depressed > healthy
depress_UI_open.contrast.ui = [+1 0 0]; % depressed > healthy
depress_NBS_open_pos = CompositeNBStest({depress_UI_open});

fprintf('\nTesting for depression-related decreases in open eyes condition...\n');
% depress_UI_open.contrast.ui = [0 -1 0 0 0]; % healthy > depressed
depress_UI_open.contrast.ui = [-1 0 0]; % healthy > depressed
depress_NBS_open_neg = CompositeNBStest({depress_UI_open});

% Test for depression-related changes in closed eyes condition
fprintf('\nTesting for depression-related increases in closed eyes condition...\n');
% depress_UI_closed.contrast.ui = [0 +1 0 0 0]; % depressed > healthy
depress_UI_closed.contrast.ui = [+1 0 0]; % depressed > healthy
depress_NBS_closed_pos = CompositeNBStest({depress_UI_closed});

fprintf('\nTesting for depression-related decreases in closed eyes condition...\n');
% depress_UI_closed.contrast.ui = [0 -1 0 0 0]; % healthy > depressed
depress_UI_closed.contrast.ui = [-1 0 0]; % healthy > depressed
depress_NBS_closed_neg = CompositeNBStest({depress_UI_closed});

%%%%%%%%%%%%%%%%%%%%%%%%%%%

% depress_UI = UI;
% depress_UI.design.ui = 1*[is_depressed, ~is_depressed, age, age.^2, sex];
% depress_UI.exchange.ui = '';

% % Test for OPEN EYES
% depress_UI.matrices.ui = mvgc_open_nets;

% % Test for increases under Depression
% depress_UI.contrast.ui = [+1 -1 0 0 0 0];
% depress_NBS_open_pos = CompositeNBStest({depress_UI})

% % Test for decreases under Depression
% depress_UI.contrast.ui = [-1 +1 0 0 0 0];
% depress_NBS_open_neg = CompositeNBStest({depress_UI})

% % Test for CLOSED EYES
% depress_UI.matrices.ui = mvgc_closed_nets;

% % Test for increases under Depression
% depress_UI.contrast.ui = [+1 -1 0 0 0 0];
% depress_NBS_closed_pos = CompositeNBStest({depress_UI})

% % Test for decreases under Depression
% depress_UI.contrast.ui = [-1 +1 0 0 0 0];
% depress_NBS_closed_neg = CompositeNBStest({depress_UI})

%% Display results

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

fprintf('\nNumber of negative clusters (should be 0):\n');
fprintf('  Open eyes: %d\n', depress_NBS_open_neg.n);
fprintf('  Closed eyes: %d\n', depress_NBS_closed_neg.n);
fprintf('=====================================\n');

%% Export results

% Assert there are no negative clusters
assert(depress_NBS_neg.n == 0);

%%%%%%%%%%%%%%%%%%%%

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