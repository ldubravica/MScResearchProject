%% Load subject excel metadata

excelFile = fullfile('Depression_Study', 'depression_data', 'Data_4_Import_REST.xlsx');
T = readtable(excelFile);
T(T.id == 544 | T.id == 571 | T.id == 572 | T.id == 527 | T.id == 535, :) = [];
% Remove 570 due to no age data
T(T.id == 570, :) = [];

%% Load MVGC matrices

addpath('./NBS_test/')
addpath('./NBS_test/NBS1.2/');
% addpath('./NBS_test/NBS_v1.2/');
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

%% Diagnostics - Check MVGC data quality

% fprintf('\n==== MVGC Data Diagnostics ====\n\n');

% fprintf('Number of subjects: %d\n', nb_subjects);
% fprintf('Depressed subjects: %d\n', sum(is_depressed));
% fprintf('Healthy subjects: %d\n', sum(~is_depressed));

% fprintf('MVGC data diagnostics:\n');
% for i = 1:min(3, nb_subjects)  % Check first 3 subjects
%     fprintf('Subject %d - Open eyes matrix:\n', i);
%     fprintf('  Min: %.4f, Max: %.4f, Mean: %.4f\n', ...
%         min(mvgc_open_nets{i}(:)), max(mvgc_open_nets{i}(:)), mean(mvgc_open_nets{i}(:)));
%     fprintf('  Has NaN: %d, Has Inf: %d\n', ...
%         any(isnan(mvgc_open_nets{i}(:))), any(isinf(mvgc_open_nets{i}(:))));
% end

% fprintf('\n===============================\n');

%% Run depression NBS tests

is_depressed = double(T_matched.MDD <= 2);
age = T_matched.age;
sex = T_matched.sex - 1;

mvgc_diff_nets = cell([nb_subjects, 1]);
for i = 1:nb_subjects
	mvgc_diff_nets{i} = mvgc_closed_nets{i} - mvgc_open_nets{i};
end

% Initialize UI structure for NBS analysis
% Default UI structure based on NBS requirements
UI = struct();
UI.method.ui = 'Run NBS'; % 'Run NBS' | 'Run FDR'
UI.test.ui = 't-test'; % statistical test - 'One Sample' | 't-test' | 'F-test'
UI.perms.ui = '5000';
UI.thresh.ui = '1.7'; % example was set to '3.1'
UI.alpha.ui = '0.05';
UI.size.ui = 'Extent'; % measure size of a network component - 'Extent' | 'Intensity'
UI.exchange.ui = '';

% pedro
% UI.design.ui = 1*[is_depressed, ~is_depressed, age, age.^2, sex];
% contrast = [+1 -1 0 0 0];

% claude
% UI.design.ui = [ones(nb_subjects, 1), is_depressed, age, age.^2, sex];
% contrast = [0 +1 0 0 0];

% short
UI.design.ui = [is_depressed, age, sex];
contrast = [+1 0 0];

% == DIAGNOSTICS ==

% Check for problems
fprintf('\nDesign matrix diagnostics (%d x %d):\n', size(UI.design.ui));

% Check correlations between predictors
corr_matrix = corrcoef(UI.design.ui);
[row, col] = find(abs(corr_matrix) > 0.9 & corr_matrix ~= 1);
for i = 1:length(row)
	fprintf('High correlation (>0.9): Column %d vs %d: %.3f\n', row(i), col(i), corr_matrix(row(i), col(i)));
end

% Check for constant columns
for i = 1:size(UI.design.ui, 2)
	if std(UI.design.ui(:,i)) < 1e-10
		fprintf('Column %d is nearly constant\n', i);
	end
end

fprintf('\n');

% == Composite NBS tests ==

fprintf('Testing open eyes condition...\n');
UI.matrices.ui = mvgc_open_nets;
UI.contrast.ui = contrast;
NBS_open_pos = CompositeNBStest({UI})
UI.contrast.ui = contrast * -1;
NBS_open_neg = CompositeNBStest({UI})

fprintf('Testing closed eyes condition...\n');
UI.matrices.ui = mvgc_closed_nets;
UI.contrast.ui = contrast;
NBS_closed_pos = CompositeNBStest({UI})
UI.contrast.ui = contrast * -1;
NBS_closed_neg = CompositeNBStest({UI})

% fprintf('Testing difference (closed - open) condition...\n');
% UI.matrices.ui = mvgc_diff_nets;
% UI.contrast.ui = contrast;
% NBS_diff_pos = CompositeNBStest({UI})
% UI.contrast.ui = contrast * -1;
% NBS_diff_neg = CompositeNBStest({UI})

%% Display results

fprintf('\n========== RESULTS SUMMARY ==========\n\n');
fprintf('Open eyes condition (depressed > healthy):\n');
fprintf('  Number of significant clusters: %d\n', NBS_open_pos.n);
if NBS_open_pos.n > 0
    fprintf('  P-values: ');
    fprintf('%.4f ', NBS_open_pos.pval);
    fprintf('\n');
end

fprintf('\nClosed eyes condition (depressed > healthy):\n');
fprintf('  Number of significant clusters: %d\n', NBS_closed_pos.n);
if NBS_closed_pos.n > 0
    fprintf('  P-values: ');
    fprintf('%.4f ', NBS_closed_pos.pval);
    fprintf('\n');
end

% fprintf('\nCondition difference (closed - open, depressed > healthy):\n');
% fprintf('  Number of significant clusters: %d\n', NBS_diff_pos.n);
% if NBS_diff_pos.n > 0
%     fprintf('  P-values: ');
%     fprintf('%.4f ', NBS_diff_pos.pval);
%     fprintf('\n');
% end

fprintf('\nNumber of negative clusters (should be 0):\n');
fprintf('  Open eyes: %d\n', NBS_open_neg.n);
fprintf('  Closed eyes: %d\n', NBS_closed_neg.n);
fprintf('\n=====================================\n');

%%%%%%%%%%%%%%%%%

fprintf('\n');
fprintf('NBS_open_pos.test_stat:\n');
disp(NBS_open_pos.test_stat{1});
fprintf('NBS_closed_pos.test_stat:\n');
disp(NBS_closed_pos.test_stat{1});

significant_connections_open = abs(NBS_open_pos.test_stat{1}) > 1.7;
fprintf('\nSignificant connections open (|t-statistic| > 1.7):\n');
disp(significant_connections_open);

significant_connections_closed = abs(NBS_closed_pos.test_stat{1}) > 1.7;
fprintf('\nSignificant connections closed (|t-statistic| > 1.7):\n');
disp(significant_connections_closed);

% Calculate Cohen's d for each connection
mean_depressed_open = mean(cat(3, mvgc_open_nets{is_depressed==1}), 3);
mean_healthy_open = mean(cat(3, mvgc_open_nets{is_depressed==0}), 3);
pooled_std_open = sqrt(var(cat(3, mvgc_open_nets{:}), [], 3));
cohens_d_open = (mean_depressed_open - mean_healthy_open) ./ pooled_std_open;

mean_depressed_closed = mean(cat(3, mvgc_closed_nets{is_depressed==1}), 3);
mean_healthy_closed = mean(cat(3, mvgc_closed_nets{is_depressed==0}), 3);
pooled_std_closed = sqrt(var(cat(3, mvgc_closed_nets{:}), [], 3));
cohens_d_closed = (mean_depressed_closed - mean_healthy_closed) ./ pooled_std_closed;

fprintf('\nCohen''s d for each open connection:\n');
disp(cohens_d_open);
fprintf('\nCohen''s d for each closed connection:\n');
disp(cohens_d_closed);
fprintf('\n');

%% Export results

ROI = {'Frontal', 'Occipital', 'Parietal', 'Sensorimotor', 'Temporal'};

% Save significant positive clusters for open eyes condition
fprintf('\nSaving open eyes condition results...\n');
T_open = array2table(zeros([0,4]), 'VariableNames', {'From', 'To', 'Tvalue', 'Significant'});
for i=1:5
  for j=1:5
    if i == j, continue; end
    T_open = [T_open; {ROI{i}, ROI{j}, NBS_open_pos.test_stat{1}(i,j), NBS_open_pos.con_mat(i,j)}];
  end
end
writetable(T_open, 'depress_network_results_open.csv');

% Save significant positive clusters for closed eyes condition
fprintf('\nSaving closed eyes condition results...\n');
T_closed = array2table(zeros([0,4]), 'VariableNames', {'From', 'To', 'Tvalue', 'Significant'});
for i=1:5
  for j=1:5
    if i == j, continue; end
    T_closed = [T_closed; {ROI{i}, ROI{j}, NBS_closed_pos.test_stat{1}(i,j), NBS_closed_pos.con_mat(i,j)}];
  end
end
writetable(T_closed, 'depress_network_results_closed.csv');