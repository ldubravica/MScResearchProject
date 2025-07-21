function [n_cnt, con_mat, pval] = NBSmultistats(STATS, GLM)
%%NBSMULTISTATS Cluster statistics for composite NBS tests
%
%   [N_CNT, CON_MAT, PVAL] = NBSMULTISTATS(STATS) computes cluster statistics
%   and null distribution of composite test given by 1D cell array STATS.
%
%   [...] = NBSMULTISTATS(STATS, GLM) uses GLM to compute test statistics on
%   the fly, as in NBSSTATS. Slower, but saves memory.
%
%   This code is based on NBSstats.m from the NBS library (tested with v1.2).
%   Reference:
%
%     https://www.nitrc.org/projects/nbs/
%
%   WARNING: this function assumes that data has already been checked by
%   CompositeNBStest() in terms of network size, clustering criterion, etc.
%
% Pedro Mediano, Oct 2020


%% Parameter checks and initialisation

%Is BGL available?
bgl=0;
if exist('components','file')==2
    %Use components.m provided by MatlabBGL, otherwise use get_components.m
    bgl=1;
end

%Number of nodes
N = STATS{1}.N; 

%Number of edges
J = N*(N-1)/2;

% Number of tests
nb_tests = length(STATS);


% Check if network is undirected by comparing number of dependent variables and
% number of nodes.
is_net_undirected = all(cellfun(@(s) size(s.test_stat, 2) == J, STATS));


% Index of elements in the matrix to test: upper triangular if undirected, all
% except diagonal if directed
if is_net_undirected
  ind_fun = @(n) triu(ones(n), 1);
else
  ind_fun = @(n) ones(n) - eye(n);
end

ind_upper = find(ind_fun(N))';


%% Step 1: find edges where all test statistics are above threshold

% Determine whether test statistics have been precomputed and determine
% index of edges exceeding the primary threshold for each dataset
ind = 1:(N*N);
for k=1:nb_tests
  if ~isempty(STATS{k}.test_stat)
      %Precomputed test statistics
      ind = intersect(ind, ind_upper(STATS{k}.test_stat(1,:) > STATS{k}.thresh)); 
      %Number of permutations
      K = size(STATS{k}.test_stat,1)-1; 
  else
      %Compute test statistics on the fly
      %Get desired number of permutations
      K = GLM{k}.perms;
      %Set to 1, since NBSglm will be called separately for each permutation
      GLM{k}.perms = 1; 
      test_stat{k} = NBSglm(GLM{k});  
      ind = intersect(ind, ind_upper(test_stat{k}(1,:) > STATS{k}.thresh));
  end
end

%Size of a component measured using extent or intensity? 
Intensity=0;
if strcmp(STATS{1}.size,'Intensity')    
    %If size measure using intensity, create an N x N matrix cotaining the 
    %test statistic for each edge minus the test statistic threshold
    %(primary threshold)
    Intensity = 1; 
    %Compute a test statistic matrix
    test_stat_mat = zeros(N,N); 

    % Test statistic is the minimum across all individual tests
    test_stat_mat(ind_upper) = inf;
    for k=1:nb_tests
      if ~isempty(STATS{k}.test_stat)
          %Precomputed
          test_stat_mat(ind_upper) = min(test_stat_mat(ind_upper), STATS{k}.test_stat(1,:) - STATS{k}.thresh);
      else
          %Not precomputed
          test_stat_mat(ind_upper) = min(test_stat_mat(ind_upper), test_stat{k}(1,:) - STATS{k}.thresh);
      end
    end

    if is_net_undirected
      test_stat_mat = (test_stat_mat+test_stat_mat');
    end
end


%% Step 2: compute cluster stats for real (unshuffled) data

adj=spalloc(N,N,length(ind)*2);
adj(ind) = 1;
if is_net_undirected
  adj = adj + adj';
end
%Only consider components comprising more than one node, equivalent to at
%least one edge
if bgl==1
    [a,sz]=components(adj | adj');
else
    [a,sz]=get_components(adj | adj');
end
ind_sz=find(sz>1);
sz_links=zeros(1,length(ind_sz));
max_sz=0; 
for i=1:length(ind_sz)
    nodes=find(ind_sz(i)==a);
    if Intensity
        %Measure size as intensity
        sz_links(i)=sum(sum(adj(nodes,nodes).*test_stat_mat(nodes,nodes)))/2;
    else
        %Measure size as extent
        sz_links(i)=sum(sum(adj(nodes,nodes)))/2;
    end
    adj(nodes,nodes)=adj(nodes,nodes)*(i+1);
    if max_sz<sz_links(i)
        max_sz=sz_links(i);
    end
end

%Subtract one to remove edges not part of a component
%Although one is also subtracted from edges comprising a component, this is 
%compensated by the (i+1) above
adj(~~adj)=adj(~~adj)-1;


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
    ind = 1:(N*N);
    test_stat_mat = zeros(N,N); 
    test_stat = {};
    for k=1:nb_tests
        if ~isempty(STATS{k}.test_stat)
            %Precomputed test statistics 
            ind = intersect(ind, ind_upper(STATS{k}.test_stat(surr_idx(k),:) > STATS{k}.thresh)); 
        else
            %Compute on the fly
            test_stat{k} = NBSglm(GLM{k});
            ind = intersect(ind, ind_upper(test_stat{k}(surr_idx(k),:) > STATS{k}.thresh));
        end
    end


    % Test statistic is the minimum across all individual tests
    test_stat_mat(ind_upper) = inf;
    for k=1:nb_tests
        %Compute a test statistic matrix
        if Intensity 
            if ~isempty(STATS{k}.test_stat)
                test_stat_mat(ind_upper) = min(test_stat_mat(ind_upper), STATS{k}.test_stat(surr_idx(k),:) - STATS{k}.thresh);
            else
                test_stat_mat(ind_upper) = min(test_stat_mat(ind_upper), test_stat{k}(surr_idx(k),:) - STATS{k}.thresh);
            end
        end
    end

    if is_net_undirected
      test_stat_mat = (test_stat_mat+test_stat_mat');
    end

    % Cluster, compute cluster statistics and save to null distribution
    adj_perm=spalloc(N,N,length(ind)*2);
    adj_perm(ind) = 1;
    if bgl==1
        [a,sz]=components(adj_perm | adj_perm');
    else
        [a,sz]=get_components(adj_perm | adj_perm');
    end
    ind_sz=find(sz>1);
    max_sz_perm=0; 
    for j=1:length(ind_sz)
        nodes=find(ind_sz(j)==a);
        if Intensity
            tmp=sum(sum(adj_perm(nodes,nodes).*test_stat_mat(nodes,nodes)))/2;
        else
            tmp=sum(sum(adj_perm(nodes,nodes)))/2;
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
        ind=find(adj==i);
        con_mat{n_cnt}=spalloc(N,N,length(ind)*2);
        con_mat{n_cnt}(ind)=1; 
        if is_net_undirected
            con_mat{n_cnt}=triu(con_mat{n_cnt},1);
        else
            con_mat{n_cnt} = con_mat{n_cnt} - diag(diag(con_mat{n_cnt}));
        end
        pval(n_cnt)=tmp;
    end
end
if n_cnt==0
    pval=[]; con_mat=[]; 
end
