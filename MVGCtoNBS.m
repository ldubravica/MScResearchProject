%% Load subject excel metadata

excelFile = fullfile('Depression_Study', 'depression_data', 'Data_4_Import_REST.xlsx');
T = readtable(excelFile);
T(T.id == 544 | T.id == 571 | T.id == 572 | T.id == 527 | T.id == 535, :) = [];
% Remove 570 due to no age data
T(T.id == 570, :) = [];

%% Load MVGC matrices

mvgc_format = 'brain_source';  % 'source' | 'brain source' | 'pca' | 'region'
save_plots = false;

addpath('./NBS_test/')
addpath('./NBS_test/NBS_v1.2/');
[mvgc_open, mvgc_closed, ~, ~] = EEGtoMVGCOptions(mvgc_format, true, false, true, false);

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

%% Run depression NBS tests

is_depressed = double(T_matched.MDD <= 2);
age = T_matched.age;
sex = T_matched.sex - 1;

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
% UI.design.ui = [is_depressed, ~is_depressed, age, age.^2, sex];
% contrast = [+1 -1 0 0 0];

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

fprintf('\nNumber of negative clusters (should be 0):\n');
fprintf('  Open eyes: %d\n', NBS_open_neg.n);
fprintf('  Closed eyes: %d\n', NBS_closed_neg.n);
fprintf('\n=====================================\n');

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

%% Plot

disp(NBS_open_pos);

plot_mvgc(...
    {NBS_open_pos.test_stat{1}, NBS_closed_pos.test_stat{1}, NBS_open_pos.test_stat{1}, NBS_closed_pos.test_stat{1}}, ...
    {'Open NBS T-Values', 'Closed NBS T-Values', 'Open NBS T-Values', 'Closed NBS T-Values'}, ...
    sprintf('MVGC NBS-Corrected T-Values (method: %s)', strrep(mvgc_format, '_', '-')), ...
    [1.7, 1.7, 1.7, 1.7]);  % Significance threshold for both matrices

function plot_mvgc(mvgc_matrices, titles, super_title, signif_threshs)
    plot_data = true;
    save_plots = true;
    region_names = {'Frontal', 'Occipital', 'Parietal', 'Sensorimotor', 'Temporal'};

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
                    if abs(mvgc_matrices{i}(r, c)) >= signif_threshs(i)
                        % Place the star slightly above center (y-0.25)
                        text(c, r-0.1, '*', 'Color', 'k', 'FontSize', 14, 'HorizontalAlignment', 'center');
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