function [] = EEGtoMVGCCompare()

    method_names = {'brain_source', 'source', 'pca', 'region'};

    fields_to_compare = { ...
        'mvgc_open_healthy_avg', 'mvgc_closed_healthy_avg', ...
        'mvgc_open_depressed_avg', 'mvgc_closed_depressed_avg', ...
        'mvgc_open_avg', 'mvgc_closed_avg', ...
        'mvgc_healthy_avg', 'mvgc_depressed_avg', ...
        'mvgc_open_pvals', 'mvgc_closed_pvals', ...
        'mvgc_open_pvals_corrected', 'mvgc_closed_pvals_corrected', ...
        'mvgc_open_tvals', 'mvgc_closed_tvals', ...
        'mvgc_open_tvals_masked', 'mvgc_closed_tvals_masked' ...
    };

    % fields_to_compare = { ...
    %     'mvgc_open_healthy_avg', 'mvgc_closed_healthy_avg', ...
    %     'mvgc_open_depressed_avg', ...
    %     'mvgc_open_avg', 'mvgc_closed_avg', ...
    %     'mvgc_open_pvals', 'mvgc_closed_pvals', ...
    %     'mvgc_open_tvals_masked' ...
    % };

    % fields_to_compare = {'mvgc_open_healthy_avg'};

    % add struct fields to fields_to_compare - TODO
    % fields_to_compare = [fields_to_compare, ...
    %     'mvgc_open_healthy', 'mvgc_closed_healthy', ...
    %     'mvgc_open_depressed', 'mvgc_closed_depressed', ...
    %     'mvgc_open', 'mvgc_closed'];

    compare_mvgc_methods(method_names, fields_to_compare, false);

end

function compare_mvgc_methods(method_names, fields_to_compare, output_csv)

    n_methods = length(method_names);
    similarity_results = struct();
    all_results = table(length(fields_to_compare));

    % Define table() of size fields_to_compare x (6) to store results
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
                v1 = matrix_vectors{i};
                v2 = matrix_vectors{j};
                mask = ~(isnan(v1) | isnan(v2));
                if nnz(mask) > 0
                    r = corr(v1(mask), v2(mask), 'type', 'Pearson');
                    cosine_sim = dot(v1(mask), v2(mask)) / (norm(v1(mask)) * norm(v2(mask)));
                    mse = mean((v1(mask) - v2(mask)).^2);
                else
                    r = NaN; cosine_sim = NaN; mse = NaN;
                end

                corr_mat(i, j) = r; corr_mat(j, i) = r;
                cosine_mat(i, j) = cosine_sim; cosine_mat(j, i) = cosine_sim;
                mse_mat(i, j) = mse; mse_mat(j, i) = mse;
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

        % Define colors
        n = 128;
        neg_colors = [linspace(0,1,n)', linspace(0,1,n)', ones(n,1)]; % blue to white
        pos_colors = [ones(n,1), linspace(1,0,n)', linspace(1,0,n)']; % white to red
        diverging_cmap = [neg_colors; pos_colors];

        subplot_similarity_matrix(1, corr_mat, 'Correlation Matrix', [-1 1], diverging_cmap);
        subplot_similarity_matrix(2, cosine_mat, 'Cosine Similarity Matrix', [0 1], pos_colors);
        subplot_similarity_matrix(3, mse_mat, 'MSE Matrix', [], neg_colors);

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
