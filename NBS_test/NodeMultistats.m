function [n_cnt, con_mat, pval] = NodeMultistats(STATS, adj, GLM)
%%NODEMULTISTATS Cluster statistics for node-level composite cluster tests
%
%   [N_CNT, CON_MAT, PVAL] = NBSMULTISTATS(STATS, A) computes cluster statistics
%   and null distribution of composite test given by 1D cell array STATS. Each
%   test represents a node in a network, and binary adjacency matrix A is used
%   to cluster significant nodes together.
%
%   [...] = NBSMULTISTATS(STATS, A, GLM) uses GLM to compute test statistics on
%   the fly, as in NBSSTATS. Slower, but saves memory.
%
%   This code is based on NBSstats.m from the NBS library (tested with v1.2).
%   Reference:
%
%     https://www.nitrc.org/projects/nbs/
%
%   WARNING: this function assumes that data has already been checked by
%   CompositeNodeClusterTest() in terms of network size, clustering criterion, etc.
%
% Pedro Mediano, Jul 2021


%% Parameter checks and initialisation
%Number of nodes
N = STATS{1}.N; 

% Number of tests
nb_tests = length(STATS);


%% Step 1: find nodes where all test statistics are above threshold

% Determine whether test statistics have been precomputed and determine
% index of nodes exceeding the primary threshold for each dataset
ind = true(1,N);
for k=1:nb_tests
  if ~isempty(STATS{k}.test_stat)
      %Precomputed test statistics
      ind = ind & (STATS{k}.test_stat(1,:) > STATS{k}.thresh); 
      %Number of permutations
      K = size(STATS{k}.test_stat,1)-1; 
  else
      %Compute test statistics on the fly
      %Get desired number of permutations
      K = GLM{k}.perms;
      %Set to 1, since NBSglm will be called separately for each permutation
      GLM{k}.perms = 1; 
      test_stat{k} = NBSglm(GLM{k});  
      ind = ind & (test_stat{k}(1,:) > STATS{k}.thresh);
  end
end

%Size of a component measured using extent or intensity? 
Intensity=0;
if strcmp(STATS{1}.size,'Intensity')    
    %If size measure using intensity, create a 1 x N matrix cotaining the 
    %test statistic for each node minus the test statistic threshold
    %(primary threshold)
    Intensity = 1; 

    % Test statistic is the minimum across all individual tests
    test_stat_mat = inf*ones(1,N);
    for k=1:nb_tests
      if ~isempty(STATS{k}.test_stat)
          %Precomputed
          test_stat_mat = min(test_stat_mat, STATS{k}.test_stat(1,:) - STATS{k}.thresh);
      else
          %Not precomputed
          test_stat_mat = min(test_stat_mat, test_stat{k}(1,:) - STATS{k}.thresh);
      end
    end
end


%% Step 2: compute cluster stats for real (unshuffled) data
cluster_labels = 1:N;
cluster_labels(~ind) = 0;
clusters = combineClusters(uint32(cluster_labels'), logical(adj | adj'), uint32(N));
ind_sz = setdiff(unique(clusters), [0]);

sz_links=zeros(1,length(ind_sz));
max_sz=0; 
for i=1:length(ind_sz)
    nodes=find(ind_sz(i)==clusters);
    if Intensity
        %Measure size as intensity
        sz_links(i)=sum(test_stat_mat(nodes));
    else
        %Measure size as extent
        sz_links(i)=length(nodes);
    end
    if max_sz<sz_links(i)
        max_sz=sz_links(i);
    end
end



%% Step 3: repeat for shuffled data

%Repeat above for each permutation
%Empirical null distribution of maximum component size
null_dist=zeros(K,1); 
p_approx=0;
%First row of test_stat is the observed test statistics, so start at the
%second row
for i=2:K+1


    % Randomly pick one dataset to shuffle, keep the others intact. The point
    % behind shuffling exactly one dataset is that it is the strongest null
    % hypothesis in which not all tests are significant.
    surr_idx = ones([1, nb_tests]);
    shuf_idx = randi(nb_tests);
    if ~isempty(STATS{shuf_idx}.test_stat)
        surr_idx(shuf_idx) = i;
    else
        surr_idx(shuf_idx) = 2;
    end


    % Compute test statistic matrix for composite test
    ind = true(1,N);
    % test_stat_mat = zeros(N,N); 
    test_stat = {};
    for k=1:nb_tests
        if ~isempty(STATS{k}.test_stat)
            %Precomputed test statistics 
            ind = ind & (STATS{k}.test_stat(surr_idx(k),:) > STATS{k}.thresh); 
        else
            %Compute on the fly
            test_stat{k} = NBSglm(GLM{k});
            ind = ind & (test_stat{k}(surr_idx(k),:) > STATS{k}.thresh);
        end
    end


    % Test statistic is the minimum across all individual tests
    % test_stat_mat(ind_upper) = inf;
    test_stat_mat = inf*ones(1,N);
    for k=1:nb_tests
        %Compute a test statistic matrix
        if Intensity 
            if ~isempty(STATS{k}.test_stat)
                test_stat_mat = min(test_stat_mat, STATS{k}.test_stat(surr_idx(k),:) - STATS{k}.thresh);
            else
                test_stat_mat = min(test_stat_mat, test_stat{k}(surr_idx(k),:) - STATS{k}.thresh);
            end
        end
    end


    % Cluster, compute cluster statistics and save to null distribution
    cluster_labels = 1:N;
    cluster_labels(~ind) = 0;
    clusters_perm = combineClusters(uint32(cluster_labels'), logical(adj | adj'), uint32(N));
    ind_sz_perm = setdiff(unique(clusters_perm), [0]);

    max_sz_perm=0; 
    for j=1:length(ind_sz_perm)
        nodes=find(ind_sz_perm(j)==clusters_perm);
        if Intensity
            tmp=sum(test_stat_mat(nodes));
        else
            tmp=length(nodes);
        end
        if tmp>max_sz_perm
            max_sz_perm=full(tmp);
        end   
    end
    null_dist(i-1)=max_sz_perm; 
    if max_sz_perm>=max_sz
        p_approx=p_approx+1;
    end
end


%% Step 4: compare observed cluster against null distribution

%Determine components satisfying alpha significance threshold
n_cnt=0; 
for i=1:length(sz_links)
    tmp=sum(null_dist>=sz_links(i))/K;
    if tmp<=STATS{1}.alpha
        n_cnt=n_cnt+1;
        con_mat{n_cnt} = zeros(1,N);
        con_mat{n_cnt}(clusters == ind_sz(i)) = 1;
        pval(n_cnt) = tmp;
    end
end
if n_cnt==0
    pval=[]; con_mat=[]; 
end

