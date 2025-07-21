function [ NBS ] = CompositeNBStest(UI_structs)
%%COMPOSITENBSTEST NBS-corrected hypothesis test for composite hypotheses
%
%   NBS = COMPOSITENBSTEST(UI_STRUCTS), where UI_STRUCTS is a 1D cell array of
%   NBS-ready UI structures, detects clusters where the hypotheses encoded by
%   all UI structs are all true.  For a single test, it reduces to the usual
%   NBS test -- i.e. COMPOSITENBSTEST({UI}) yields (statistically) the same
%   result as NBSRUN(UI).
%
%   This code depends on the NBS library (tested with v1.2). Reference:
%
%     https://www.nitrc.org/projects/nbs/
%
% Pedro Mediano, Oct 2020

%% Parameter checks and initialisation
if isempty(UI_structs) || ~iscell(UI_structs) || ~isvector(UI_structs)
  error('Input must be a 1D cell array of NBS-ready UI structures.');
end
nb_tests = length(UI_structs);

% Silently add NBS to current path
% p = mfilename('fullpath');
% addpath(strrep(p, 'CompositeNBStest', 'private/NBS_v1.2'));

% Check that stats configuration and input dimension match across datasets
for j=2:nb_tests
  if ~all([strcmp(UI_structs{1}.perms.ui,  UI_structs{j}.perms.ui),
           strcmp(UI_structs{1}.method.ui, UI_structs{j}.method.ui),
           strcmp(UI_structs{1}.alpha.ui,  UI_structs{j}.alpha.ui),
           strcmp(UI_structs{1}.size.ui,   UI_structs{j}.size.ui)])
    error('Parameters perms, method, alpha, and size must match for all tests.');
  end
end

% Check if any of the networks are undirected
is_net_undirected = true;
for j=1:nb_tests
  data = readUInum(UI_structs{j}.matrices.ui);
  if ~iscell(data)
    data = squeeze(num2cell(data, [1,2]));
  end

  is_net_undirected = all(cellfun(@(m) all(all(isnan(m - m') | (abs(m - m') < 1e-8))), data));

  if ~is_net_undirected
    break;
  end

end


%% Get and reshape input
GLM   = cell([nb_tests, 1]);
STATS = cell([nb_tests, 1]);
DIMS  = cell([nb_tests, 1]);
for j=1:nb_tests
  [GLM{j}, STATS{j}, DIMS{j}] = PreprocessUIStruct(UI_structs{j}, is_net_undirected);

  % Normalise dependent variables by their second moment to ensure variances
  % of all null distributions are comparable.
  GLM{j}.y = GLM{j}.y./repmat(sqrt(mean(GLM{j}.y.^2, 1)), [DIMS{j}.observations, 1]);

  if DIMS{1}.nodes ~= DIMS{j}.nodes
    error('All inputs must have the same number of nodes.');
  end
end


%% Run GLMs for each dataset and save test stats
for j=1:nb_tests
  STATS{j}.test_stat = NBSglm(GLM{j});
end


%% Run modified version of main NBS algorithm and find significant components
[NBS.n, NBS.con_mat, NBS.pval] = NBSmultistats(STATS);


%% Reshape results back and return
% Copy test statistics to NBS strucutre so that they can be displayed with
% each link
if is_net_undirected
  ind = find(triu(ones(DIMS{1}.nodes, DIMS{1}.nodes),1));
else
  ind = find(ones(DIMS{1}.nodes) - eye(DIMS{1}.nodes));
end
NBS.test_stat = cell([nb_tests, 1]);
for k=1:nb_tests
  NBS.test_stat{k} = zeros(STATS{k}.N, STATS{k}.N);
  NBS.test_stat{k}(ind) = STATS{k}.test_stat(1,:);
  if is_net_undirected
    NBS.test_stat{k} = NBS.test_stat{k} + NBS.test_stat{k}';
  end
end

end


%%%%%%%%%%%%%%%%%%%%%  AUXILIARY FUNCTIONS  %%%%%%%%%%%%%%%%%%%%%%%%%

function [GLM,STATS,DIMS] = PreprocessUIStruct(UI, is_undirected)

if nargin < 2 || isempty(is_undirected)
  is_undirected = true;
end

%Assume UI is valid to begin with
%Can be set to zero after reading UI or performing error checking
UI.method.ok=1;
UI.design.ok=1;
UI.contrast.ok=1;
UI.thresh.ok=1;
UI.test.ok=1;
UI.matrices.ok=1;
UI.perms.ok=1;
UI.alpha.ok=1;
UI.size.ok=1;
UI.exchange.ok=1;
UI.node_coor.ok=1;
UI.node_label.ok=1;
    
%Read UI and assign to appropriate structure
%Connectivity matrices

if isstr(UI.matrices.ui) && exist(fileparts(UI.matrices.ui),'dir')
    [GLM.y,UI.matrices.ok,DIMS]=read_matrices([fileparts(UI.matrices.ui),filesep], is_undirected);
    if ~UI.matrices.ok
        [GLM.y,UI.matrices.ok,DIMS]=read_matrices(UI.matrices.ui, is_undirected);
    end
else
    [GLM.y,UI.matrices.ok,DIMS]=read_matrices(UI.matrices.ui, is_undirected);
end


%Design matrix
[GLM.X,UI.design.ok,DIMS]=read_design(UI.design.ui,DIMS); 
%Contrast
[GLM.contrast,UI.contrast.ok]=read_contrast(UI.contrast.ui,DIMS, UI.test.ui);
%Exchange blocks for permutation [optional]
[tmp,UI.exchange.ok]=read_exchange(UI.exchange.ui,DIMS);
if UI.exchange.ok
    GLM.exchange=tmp; 
elseif isfield(GLM,'exchange')
    GLM=rmfield(GLM,'exchange');
end
%Test statistic
try GLM.test=UI.test.ui; 
    if strcmp(GLM.test,'One Sample')
        GLM.test='onesample';
    elseif strcmp(GLM.test,'t-test')
        GLM.test='ttest';
    elseif strcmp(GLM.test,'F-test')
        GLM.test='ftest';
    end
catch; UI.test.ok=0; end
%Number of permutations
try GLM.perms=str2num(UI.perms.ui); catch; UI.perms.ok=0; end 
try if ~isnumeric(GLM.perms) || ~(GLM.perms>0)
    UI.perms.ok=0; end
catch; UI.perms.ok=0; end
%Test statistic threshold
try STATS.thresh=str2num(UI.thresh.ui); catch; UI.thresh.ok=0; end
try if ~isnumeric(STATS.thresh) || ~(STATS.thresh>0)
    UI.thresh.ok=0; end
catch; UI.thresh.ok=0; end
%Corrected p-value threshold
try STATS.alpha=str2num(UI.alpha.ui); catch; UI.alpha.ok=0; end 
try if ~isnumeric(STATS.alpha) || ~(STATS.alpha>0)
    UI.alpha.ok=0; end
catch; UI.alpha.ok=0; end
%Component size 
try STATS.size=UI.size.ui; catch; UI.size.ok=0; end 
%Number of nodes
STATS.N=DIMS.nodes;

%Do error checking on user inputs
[msg, stop] = errorcheck(UI, DIMS, []);

if stop
  error(msg{1});
end

end


%% Small wrapper around NBS's readUI function which can simply pass numeric
% data (to avoid having to write it to a file and then read it again)
function [ data ] = readUInum(UI)
  if isnumeric(UI) || iscell(UI)
    data = UI;
  else
    data = readUI(UI);
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Code below has been copied with edits from NBSrun.m in NBS v1.2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Read connectivity matrices and vectorize the upper triangle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y,ok,DIMS]=read_matrices(Name, is_undirected)
    if nargin < 2 || isempty(is_undirected)
      is_undirected = true;
    end
    if is_undirected
      ind_fun = @(n) triu(ones(n), 1);
    else
      ind_fun = @(n) ones(n) - eye(n);
    end
    ok=1;
    data=readUInum(Name);
    if ~isempty(data)
        [nr,nc,ns]=size(data);
        if nr==nc && ns>0 && ~iscell(data) && isnumeric(data)
            ind_upper=find(ind_fun(nr));
            y=zeros(ns,length(ind_upper));
            %Collapse matrices
            for i=1:ns
                tmp=data(:,:,i);
                y(i,:)=tmp(ind_upper);
            end
        elseif iscell(data)
            [nr,nc]=size(data{1});
            ns=length(data);
            if nr==nc && ns>0
                ind_upper=find(ind_fun(nr));
                y=zeros(ns,length(ind_upper));
                %Collapse matrices
                for i=1:ns
                    [nrr,ncc]=size(data{i});
                    if nrr==nr && ncc==nc && isnumeric(data{i})
                        y(i,:)=data{i}(ind_upper);
                    else
                        ok=0; y=[]; 
                        break
                    end
                end
            else
                ok=0; y=[];
            end
        end
    else
        ok=0; y=[];
    end
    if ok==1
        %Number of nodes
        DIMS.nodes=nr;
        %Number of matrices
        DIMS.observations=ns;
    else
        %Number of nodes
        DIMS.nodes=0;
        %Number of matrices
        DIMS.observations=0;
    end
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Read design matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,ok,DIMS]=read_design(Name,DIMS)
ok=1;
data=readUInum(Name);
if ~isempty(data)
    [nr,nc,ns]=size(data);
    if nr==DIMS.observations && nc>0 && ns==1 && isnumeric(data) 
        X=data; 
    else
        ok=0; X=[];
    end
else
    ok=0; X=[];
end
clear data
if ok==1
    %Number of predictors
    DIMS.predictors=nc;
else
    DIMS.predictors=0;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Read contrast
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [contrast,ok]=read_contrast(Name,DIMS,testname)
ok=1; 
data=readUInum(Name);
if ~isempty(data)
    [nr,nc,ns]=size(data); 
    if strcmp(testname, 'Composite F-test') && nc==DIMS.predictors && ns==1 && isnumeric(data) 
        contrast=data;
    elseif nr==1 && nc==DIMS.predictors && ns==1 && isnumeric(data) 
        contrast=data; 
    else
        ok=0; contrast=[];
    end
else
    ok=0; contrast=[];
end     
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Read node coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [node_coor,ok]=read_node_coor(Name,DIMS)
ok=1;
data=readUInum(Name);
if ~isempty(data)
    [nr,nc,ns]=size(data);
    if nr==DIMS.nodes && nc==3 && ns==1 && isnumeric(data)
        node_coor=data; 
    else
        ok=0; node_coor=[];
    end
else
    ok=0; node_coor=[];
end        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Read node labels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [node_label,ok]=read_node_label(Name,DIMS)
ok=1;
data=readUInum(Name);
if ~isempty(data)
    [nr,nc,ns]=size(data);
    if nr==DIMS.nodes && nc==1 && ns==1
        node_label=data; 
    else
        ok=0; node_label=[];
    end
else
    ok=0; node_label=[]; 
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Read permutation exchange blocks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [exchange,ok]=read_exchange(Name,DIMS)
ok=1;
data=readUInum(Name);
if ~isempty(data)
    [nr,nc,ns]=size(data);
    if nr==DIMS.observations && nc==1 && ns==1
        exchange=data'; 
    else
        if nc==DIMS.observations && nr==1 && ns==1
            exchange=data;
        else
        ok=0; exchange=[];
        end
    end
else
    ok=0; exchange=[];
end
    %Set up exchange blocks
    blks=unique(exchange); 
    %Number of blocks
    n_blks=length(blks);
    %Number of observations per block
    sz_blk=length(exchange)/n_blks;
    if rem(sz_blk,1)>0
        ok=0;
        exchange=[];
    end
end
    
function [msg,stop]=errorcheck(UI,DIMS,S)
    stop=1;
    %Mandatroy UI
    %UI.method.ok %no need to check
    if ~UI.matrices.ok
        msg={'Stop: Connectivity Matrices not found or inconsistent'};
        try set(S.DATA.matrices.text,'ForegroundColor','red');
        catch; end
        return;
    end
    if ~UI.design.ok
        msg={'Stop: Design Matrix not found or inconsistent'};
        try set(S.STATS.design.text,'ForegroundColor','red');
        catch; end
        return;
    end
    if ~UI.contrast.ok 
        msg={'Stop: Contrast not found or inconsistent'};
        try set(S.STATS.contrast.text,'ForegroundColor','red');
        catch; end
        return;
    end
    if ~UI.thresh.ok
        msg={'Stop: Threshold not found or inconsistent'};
        try set(S.STATS.thresh.text,'ForegroundColor','red');
        catch; end
        return;
    end
    if ~UI.test.ok 
        msg={'Stop: Statistical Test not found or inconsistent'};
        try set(S.STATS.test.text,'ForegroundColor','red');
        catch; end
        return;
    end
    if ~UI.perms.ok
        msg={'Stop: Permutations not found or inconsistent'};
        try set(S.ADV.perms.text,'ForegroundColor','red');
        catch; end
        return;
    end        
    if ~UI.alpha.ok
        msg={'Stop: Significance not found or inconsistent'};
        try set(S.ADV.alpha.text,'ForegroundColor','red');
        catch; end
        return;
    end
    if ~UI.size.ok
        msg={'Stop: Component Size not found or inconsistent'};
        try set(S.ADV.size.text,'ForegroundColor','red');
        catch; end
        return;
    end
    stop=0;
    
    msg=[{sprintf('Nodes: %d',DIMS.nodes)};...
             {sprintf('Observations: %d',DIMS.observations)};...
             {sprintf('Predictors: %d',DIMS.predictors)}]; 
    
    %Optional, but mandatory for NBSview
    if ~UI.node_coor.ok
        msg=[msg;{'Node Coordinates: No'}];
        try set(S.DATA.node_coor.text,'ForegroundColor','red');
        catch; end
    end
    
    %Optional
    if ~UI.exchange.ok
        msg=[msg;{'Exchange Blocks: No'}];
        try set(S.ADV.exchange.text,'ForegroundColor','red');
        catch; end
    else
        msg=[msg;{'Exchange Blocks: Yes'}];
    end
    if ~UI.node_label.ok
        msg=[msg;{'Node Labels: No'}];
        try set(S.DATA.node_label.text,'ForegroundColor','red');
        catch; end
    end
end
