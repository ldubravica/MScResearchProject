function [n_cnt, con_mat, pval] = NBSmultistats(STATS)
%NBSMULTISTATS Network-Based Statistic for multiple test statistics
%
%   [N_CNT, CON_MAT, PVAL] = NBSMULTISTATS(STATS) computes network 
%   components where ALL test statistics in the cell array STATS exceed 
%   their respective thresholds. This implements a composite hypothesis 
%   test where clusters must be significant across multiple conditions.
%
%   Input:
%       STATS:      Cell array of STATS structures, each containing:
%                   - test_stat: Test statistics matrix
%                   - thresh: Primary threshold
%                   - alpha: Significance level
%                   - N: Number of nodes
%                   - size: 'Extent' or 'Intensity'
%
%   Output:
%       N_CNT:      Number of significant components
%       CON_MAT:    Cell array of adjacency matrices for significant components
%       PVAL:       P-values for each significant component
%
%   This function extends NBSstats to handle composite hypotheses where
%   edges must survive thresholding in ALL test statistics.

%% Initialize
nb_tests = length(STATS);
if nb_tests == 0
    error('STATS cell array cannot be empty');
end

% Get dimensions from first test
N = STATS{1}.N;
J = N * (N - 1) / 2;
ind_upper = find(triu(ones(N, N), 1));

% Check all tests have same dimensions and parameters
for i = 2:nb_tests
    if STATS{i}.N ~= N
        error('All tests must have the same number of nodes');
    end
    if ~strcmp(STATS{1}.size, STATS{i}.size)
        error('All tests must use the same size measure');
    end
end

% Get number of permutations (should be same for all tests)
K = size(STATS{1}.test_stat, 1) - 1;

%% Find edges that survive ALL thresholds (composite hypothesis)
% For observed data (first row)
surviving_edges = true(1, J);
for i = 1:nb_tests
    test_stat_row = STATS{i}.test_stat(1, :);
    thresh_val = STATS{i}.thresh;
    % Ensure thresh is scalar or same size as test_stat_row
    if isscalar(thresh_val)
        surviving_edges = surviving_edges & (test_stat_row > thresh_val);
    else
        % Ensure thresh_val is same size as test_stat_row
        if length(thresh_val) == length(test_stat_row)
            surviving_edges = surviving_edges & (test_stat_row > thresh_val(:)');
        else
            error('Threshold dimensions must match test statistic dimensions');
        end
    end
end
ind_observed = ind_upper(surviving_edges);

%% Set up component size measurement
Intensity = strcmp(STATS{1}.size, 'Intensity');

if Intensity
    % Create composite test statistic matrix for intensity measurement
    test_stat_mat_obs = zeros(N, N);
    for i = 1:nb_tests
        temp_mat = zeros(N, N);
        temp_mat(ind_upper) = STATS{i}.test_stat(1, :) - STATS{i}.thresh;
        test_stat_mat_obs = test_stat_mat_obs + temp_mat;
    end
    test_stat_mat_obs = (test_stat_mat_obs + test_stat_mat_obs') / nb_tests;
end

%% Find components in observed data
adj_obs = sparse(N, N);
adj_obs(ind_observed) = 1;
adj_obs = adj_obs + adj_obs';

% Find connected components
if exist('components', 'file') == 2
    [a, sz] = components(adj_obs);
else
    [a, sz] = get_components(adj_obs);
end

% Only consider components with more than one node (at least one edge)
ind_sz = find(sz > 1);
sz_links_obs = zeros(1, length(ind_sz));
max_sz_obs = 0;

for i = 1:length(ind_sz)
    nodes = find(ind_sz(i) == a);
    if Intensity
        sz_links_obs(i) = sum(sum(adj_obs(nodes, nodes) .* test_stat_mat_obs(nodes, nodes))) / 2;
    else
        sz_links_obs(i) = sum(sum(adj_obs(nodes, nodes))) / 2;
    end
    adj_obs(nodes, nodes) = adj_obs(nodes, nodes) * (i + 1);
    if max_sz_obs < sz_links_obs(i)
        max_sz_obs = sz_links_obs(i);
    end
end

% Adjust adjacency matrix indexing
    % Find edges surviving ALL thresholds in this permutation
    surviving_edges_perm = true(1, J);
    for i = 1:nb_tests
        test_stat_row = STATS{i}.test_stat(perm, :);
        thresh_val = STATS{i}.thresh;
        % Ensure thresh is scalar or same size as test_stat_row
        if isscalar(thresh_val)
            surviving_edges_perm = surviving_edges_perm & (test_stat_row > thresh_val);
        else
            surviving_edges_perm = surviving_edges_perm & (test_stat_row > thresh_val(:)');
        end
    end
for perm = 2:K+1
    % Find edges surviving ALL thresholds in this permutation
    surviving_edges_perm = true(1, J);
    for i = 1:nb_tests
        surviving_edges_perm = surviving_edges_perm & ...
            (STATS{i}.test_stat(perm, :) > STATS{i}.thresh);
    end
    ind_perm = ind_upper(surviving_edges_perm);
    
    if Intensity
        % Create composite test statistic matrix for this permutation
        test_stat_mat_perm = zeros(N, N);
        for i = 1:nb_tests
            temp_mat = zeros(N, N);
            temp_mat(ind_upper) = STATS{i}.test_stat(perm, :) - STATS{i}.thresh;
            test_stat_mat_perm = test_stat_mat_perm + temp_mat;
        end
        test_stat_mat_perm = (test_stat_mat_perm + test_stat_mat_perm') / nb_tests;
    end
    
    % Create adjacency matrix for this permutation
    adj_perm = sparse(N, N);
    adj_perm(ind_perm) = 1;
    adj_perm = adj_perm + adj_perm';
    
    % Find components
    if exist('components', 'file') == 2
        [a_perm, sz_perm] = components(adj_perm);
    else
        [a_perm, sz_perm] = get_components(adj_perm);
    end
    
    % Find maximum component size in this permutation
    ind_sz_perm = find(sz_perm > 1);
    max_sz_perm = 0;
    
    for j = 1:length(ind_sz_perm)
        nodes = find(ind_sz_perm(j) == a_perm);
        if Intensity
            tmp = sum(sum(adj_perm(nodes, nodes) .* test_stat_mat_perm(nodes, nodes))) / 2;
        else
            tmp = sum(sum(adj_perm(nodes, nodes))) / 2;
        end
        if tmp > max_sz_perm
            max_sz_perm = full(tmp);
        end
    end
    
    null_dist(perm - 1) = max_sz_perm;
    
    if mod(perm - 1, 500) == 0
        fprintf('Completed %d/%d permutations\n', perm - 1, K);
    end
end

%% Determine significant components
n_cnt = 0;
con_mat = {};
pval = [];

alpha = STATS{1}.alpha; % Use alpha from first test (should be same for all)

for i = 1:length(sz_links_obs)
    p_val = sum(null_dist >= sz_links_obs(i)) / K;
    if p_val <= alpha
        n_cnt = n_cnt + 1;
        
        % Create adjacency matrix for this significant component
        ind_component = find(adj_obs == i);
        con_mat{n_cnt} = sparse(N, N);
        con_mat{n_cnt}(ind_component) = 1;
        con_mat{n_cnt} = triu(con_mat{n_cnt}, 1);
        
        pval(n_cnt) = p_val;
    end
end

if n_cnt == 0
    con_mat = {};
    pval = [];
end

fprintf('Found %d significant components\n', n_cnt);

end


%% Helper function for connected components (if MatlabBGL not available)
function [a, sz] = get_components(adj)
    N = size(adj, 1);
    a = zeros(N, 1);
    sz = [];
    component_id = 0;
    
    for i = 1:N
        if a(i) == 0  % Unvisited node
            component_id = component_id + 1;
            component_nodes = dfs(adj, i, N);
            a(component_nodes) = component_id;
            sz(component_id) = length(component_nodes);
        end
    end
end

function component_nodes = dfs(adj, start_node, N)
    visited = false(N, 1);
    stack = start_node;
    component_nodes = [];
    
    while ~isempty(stack)
        node = stack(end);
        stack(end) = [];
        
        if ~visited(node)
            visited(node) = true;
            component_nodes = [component_nodes, node];
            
            % Add unvisited neighbors to stack
            neighbors = find(adj(node, :));
            for neighbor = neighbors
                if ~visited(neighbor)
                    stack = [stack, neighbor];
                end
            end
        end
    end
end