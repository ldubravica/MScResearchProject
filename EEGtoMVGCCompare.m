function [] = EEGtoMVGCCompare(output_csv, output_figures)

    if nargin < 1 || isempty(output_csv)
        output_csv = false;
    end
    if nargin < 2 || isempty(output_figures)
        output_figures = false;
    end

    method_names = {'brain_source', 'source', 'pca', 'region'};

    values_to_compare = { ...
        'mvgc_open_healthy_avg', 'mvgc_closed_healthy_avg', ...
        'mvgc_open_depressed_avg', 'mvgc_closed_depressed_avg', ...
        'mvgc_open_avg', 'mvgc_closed_avg', ...
        'mvgc_healthy_avg', 'mvgc_depressed_avg', ...
        'mvgc_open_pvals', 'mvgc_closed_pvals', ...
        'mvgc_open_pvals_corrected', 'mvgc_closed_pvals_corrected', ...
        'mvgc_open_tvals', 'mvgc_closed_tvals', ...
        'mvgc_open_tvals_masked', 'mvgc_closed_tvals_masked' ...
    };

    % values_to_compare = {'mvgc_open_healthy_avg'};

    % add struct fields to values_to_compare - TODO
    % values_to_compare = [values_to_compare, ...
    %     'mvgc_open_healthy', 'mvgc_closed_healthy', ...
    %     'mvgc_open_depressed', 'mvgc_closed_depressed', ...
    %     'mvgc_open', 'mvgc_closed'];

    % Compare MVGC matrices and derived metrics
    compare_mvgc_methods(method_names, values_to_compare, output_csv, output_figures);

    % Compare ss_info parameters across methods
    info_to_compare = {'morder', 'rhoA', 'rhoB', 'acdec', 'sigspd', 'mii', 'mmii'};
    info_descriptions = {'Model Order', 'Rho A', 'Rho B', 'Autocovariance Decay', ...
                         'Sigma SPD (the covariance matrix is symmetric positive definite)', ...
                         'Mutual Information Index', 'Multivariate Mutual Information Index'};
    compare_ss_info(method_names, info_to_compare, info_descriptions, output_csv);

end

function compare_mvgc_methods(method_names, fields_to_compare, output_csv, output_figures)

    n_methods = length(method_names);
    similarity_results = struct();
    all_results = table('Size', [length(fields_to_compare), 6], ...
        'VariableTypes', {'string', 'string', 'string', 'double', 'double', 'double'}, ...
        'VariableNames', {'Field', 'Method1', 'Method2', 'Correlation', 'Cosine', 'MSE'});

    row_idx = 1;
    for f = 1:length(fields_to_compare)
        field = fields_to_compare{f};
        fprintf('\n=== Comparing field: %s ===\n\n', field);

        matrix_vectors = cell(n_methods, 1);
        for i = 1:n_methods
            data = load(sprintf('mvgc_results_%s.mat', method_names{i}));
            M = data.(field);

            fprintf('%s:\n', method_names{i});
            disp(M);

            M(logical(eye(size(M)))) = NaN;  % Exclude self-connections
            matrix_vectors{i} = M(:);        % Flatten
        end

        % Pairwise similarity matrices
        corr_mat = NaN(n_methods);
        cosine_mat = NaN(n_methods);
        mse_mat = NaN(n_methods);

        for i = 1:n_methods
            for j = i+1:n_methods
                m1 = matrix_vectors{i};
                m2 = matrix_vectors{j};
                mask = ~isnan(m1) & ~isnan(m2);  % Mask to ignore NaNs

                if nnz(mask) == 0
                    continue;
                end

                % Correlation
                corr_mat(i,j) = corr(m1(mask), m2(mask), 'Rows', 'pairwise');
                corr_mat(j,i) = corr_mat(i,j);  % Symmetric

                % Cosine similarity
                cosine_mat(i,j) = dot(m1(mask), m2(mask)) / (norm(m1(mask)) * norm(m2(mask)));
                cosine_mat(j,i) = cosine_mat(i,j);  % Symmetric

                % Mean Squared Error
                mse_mat(i,j) = mean((m1(mask) - m2(mask)).^2, 'omitnan');
                mse_mat(j,i) = mse_mat(i,j);  % Symmetric
            end
        end

        similarity_results.(field).correlation = corr_mat;
        similarity_results.(field).cosine = cosine_mat;
        similarity_results.(field).mse = mse_mat;

        % --- Display Results ---

        fprintf('\nResults for field "%s":\n\n', field);
        fprintf('Correlation Matrix:\n');
        disp(corr_mat);
        fprintf('Cosine Similarity Matrix:\n');
        disp(cosine_mat);
        fprintf('MSE Matrix:\n');
        disp(mse_mat);

        % --- Visualization ---

        figure('Name', field, 'NumberTitle', 'off');
        sgtitle(sprintf('MVGC Similarity - %s', strrep(field, '_', ' ')), 'FontSize', 16);

        % Define colors
        n = 128;
        neg_colors = [linspace(0,1,n)', linspace(0,1,n)', ones(n,1)]; % blue to white
        pos_colors = [ones(n,1), linspace(1,0,n)', linspace(1,0,n)']; % white to red
        diverging_cmap = [neg_colors; pos_colors];

        subplot_similarity_matrix(1, corr_mat, 'Correlation Matrix', [-1 1], diverging_cmap);
        subplot_similarity_matrix(2, cosine_mat, 'Cosine Similarity Matrix', [0 1], pos_colors);
        subplot_similarity_matrix(3, mse_mat, 'MSE Matrix', [], neg_colors);

        if (output_figures)
            if ~exist(fullfile('output', 'images'), 'dir')
                mkdir(fullfile('output', 'images'));
            end
            filename = sprintf('mvgc_similarity - %s - %s.png', field, datetime('now', 'Format', 'yyyy_MM_dd'));
            saveas(gcf, fullfile('output', 'images', filename));
            % close(gcf);  % Close the figure after saving
        end

        % --- Store results ---

        for i = 1:n_methods
            for j = i+1:n_methods
                all_results(row_idx, :) = { ...
                    string(field), ...
                    string(method_names{i}), ...
                    string(method_names{j}), ...
                    corr_mat(i,j), ...
                    cosine_mat(i,j), ...
                    mse_mat(i,j)};
                row_idx = row_idx + 1;
            end
        end
    end

    % --- Export CSV ---

    if output_csv
        writetable(all_results, 'mvgc_similarity.csv');
    end

    function [] = subplot_similarity_matrix(idx, similarity_mat, title_str, range, colormap_name)
        subplot(1,3,idx);
        if range % Check if range is provided
            imagesc(similarity_mat, range);
        else
            imagesc(similarity_mat);
        end
        colorbar;
        title(title_str);
        xticks(1:n_methods); yticks(1:n_methods);
        xticklabels(strrep(method_names, '_', ' ')); 
        yticklabels(strrep(method_names, '_', ' '));
        xtickangle(45);
        axis square;
        colormap(subplot(1,3,idx), colormap_name);
        set(gca, 'Color', [0 0 0]);  % Set background to black
        set(get(gca, 'Children'), 'AlphaData', ~isnan(similarity_mat));  % Set all NaNs to transparent

        % Write the values in the center of each cell
        for m1 = 1:length(method_names)
            for m2 = 1:length(method_names)
                text(m2, m1, sprintf('%.2f', similarity_mat(m1, m2)), ...
                        'Color', 'k', 'FontSize', 11, 'HorizontalAlignment', 'center');
            end
        end
    end

end

function compare_ss_info(method_names, ss_fields, ss_descriptions, output_csv)
    n_methods = length(method_names);
    ss_data_open = struct();
    ss_data_closed = struct();

    fprintf('\n');

    % Load all data for open/closed
    for i = 1:n_methods
        filename = ['mvgc_results_' method_names{i} '.mat'];
        data = load(filename, 'ss_info_open', 'ss_info_closed');
        ss_data_open.(method_names{i}) = data.ss_info_open;
        ss_data_closed.(method_names{i}) = data.ss_info_closed;
    end

    values_open = struct();
    values_closed = struct();

    calcs_open = struct();
    calcs_closed = struct();

    for method_idx = 1:n_methods
        method = method_names{method_idx};
        calcs_open.(method) = struct();
        calcs_closed.(method) = struct();

        for field_idx = 1:length(ss_fields)
            field = ss_fields{field_idx};

            subj_names = fieldnames(ss_data_open.(method));
            n_subj = numel(subj_names);

            vals_open = zeros(1, n_subj);
            vals_closed = zeros(1, n_subj);

            for s = 1:n_subj
                subj = subj_names{s};
                vals_open(s) = ss_data_open.(method).(subj).(field);
                vals_closed(s) = ss_data_closed.(method).(subj).(field);
            end

            values_open.(method).(field) = vals_open;
            values_closed.(method).(field) = vals_closed;

            if ~isempty(vals_open)
                calcs_open.(method).(field) = mean(vals_open);
                calcs_open.(method).([field '_std']) = std(vals_open);
            else
                calcs_open.(method).(field) = NaN;
                calcs_open.(method).([field '_std']) = NaN;
            end

            if ~isempty(vals_closed)
                calcs_closed.(method).(field) = mean(vals_closed);
                calcs_closed.(method).([field '_std']) = std(vals_closed);
            else
                calcs_closed.(method).(field) = NaN;
                calcs_closed.(method).([field '_std']) = NaN;
            end
        end
    end

    % Display results
    for field_idx = 1:length(ss_fields)
        field = ss_fields{field_idx};
        % description = ss_descriptions{field_idx};

        fprintf('\n\033[1mSS Info Field: %s\033[0m\n\n', field); % Bold in most terminals
        fprintf('%-12s | %-30s | %-30s\n', 'Method', 'Open (mean ± std)', 'Closed (mean ± std)');
        fprintf('%s\n', repmat('-', 1, 80));
        for method_idx = 1:n_methods
            method = method_names{method_idx};
            fprintf('%-12s | %-30s | %-30s\n', method, ...
                sprintf('%.4f ± %.4f', calcs_open.(method).(field), calcs_open.(method).([field '_std'])), ...
                sprintf('%.4f ± %.4f', calcs_closed.(method).(field), calcs_closed.(method).([field '_std'])));
        end
        fprintf('\n');
    end

    % Export CSV
    if output_csv
        all_results = table('Size', [length(ss_fields), 6], ...
            'VariableTypes', {'string', 'string', 'string', 'double', 'double', 'double'}, ...
            'VariableNames', {'Field', 'Method', 'Open Mean', 'Open Std', 'Closed Mean', 'Closed Std'});

        row_idx = 1;
        for field_idx = 1:length(ss_fields)
            field = ss_fields{field_idx};
            for method_idx = 1:n_methods
                method = method_names{method_idx};
                all_results(row_idx, :) = { ...
                    string(field), ...
                    string(method), ...
                    calcs_open.(method).(field), ...
                    calcs_open.(method).([field '_std']), ...
                    calcs_closed.(method).(field), ...
                    calcs_closed.(method).([field '_std'])};
                row_idx = row_idx + 1;
            end
        end

        writetable(all_results, 'ss_info_comparison.csv');
    end

end
