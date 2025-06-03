function [source_ts,aals] = SourceRecon(filename,filepath)
%% SourceTimeSeries
%
% Get the time series (i.e. 'virtual sensors') of sources at the centroid of
% outer AAL regions for a .set EEG dataset located at `filename`.
%
% Inputs:
%   filename -- string
%
% Outputs:
%   source_ts -- N-by-T array with beamformed source time series
%
% Pedro Mediano, Hardik Rajpal 2019
%% Defining the ROIs to select
Frontal = [3:16, 19:20, 23:26];
Occipital = [43:54];
Parietal = [59:70];
Sensorimotor = [1,2,17,18,57,58];
Temporal = [81:90];
select_aal_idx = sort(cat(2,Frontal, Occipital, Parietal, Sensorimotor, Temporal));

%% Add paths and load relevant Fieldtrip files
fieldtrip_folder = 'E:\PhD\fieldtrip-20190724\fieldtrip-20190724'; %Change to the relavant FieldTrip Folder
addpath(fieldtrip_folder);
ft_defaults;
eeglab_path = 'E:\PhD\eeglab_current\eeglab2019_1\'; % Not needed if not loading data in EEGLab format
addpath(genpath(eeglab_path));

addpath([fieldtrip_folder, '/external/eeglab']); % Not needed if not loading data in EEGLab format

if ~exist('AAL.mat', 'file')
  error('AAL file missing.');
end
load('AAL.mat');

source_pos = AAL_centroids(select_aal_idx,:);
aals = AAL_Labels(select_aal_idx);

%% Load data
data = pop_loadset('filename',filename,'filepath',filepath); % Loading data in EEGLab Format
data = eeg_regepochs(data,'recurrence',2,'limits',[0,2],'rmbase',NaN); % Epoching data into 2 second epochs in EEGLab
data = eeglab2fieldtrip(data,'raw','none') % Converting EEGLab data to fieldtrip format
% If loading data in a different format, skip above and convert the epoched
% data directly into fieldtrip format.

% cfg = struct('dataset', filename, 'channel', {'eeg'});
%data = ft_preprocessing(cfg);

% Calculate means and covariances for source analysis
cfg                   = [];
cfg.covariance        = 'yes';
cfg.covariancewindow  = 'all';
cfg.channel           = 'EEG';
cfg.keeptrials        = 'yes';
tlock                 = ft_timelockanalysis(cfg, data);

%% Load Fieldtrip's template headmodel
vol = load([fieldtrip_folder, '/template/headmodel/standard_bem.mat']);
template_vol = vol.vol;

% Run source analysis on the provided positions
cfg                    = [];
cfg.method             = 'lcmv';
cfg.grid               = [];
cfg.grid.pos           = source_pos;
cfg.vol                = template_vol;
cfg.lcmv.keepfilter    = 'yes';
cfg.lcmv.lambda        = '5%';
cfg.lcmv.fixedori      = 'yes';
source  = ft_sourceanalysis(cfg, tlock);

% Use filters to reconstruct the source time series from EEG time series
source_filter = cell2mat(source.avg.filter);
[nsource,nelec] = size(source_filter);
[ntrials,nelec,T] = size(tlock.trial);
source_ts = zeros(ntrials,nsource,T);
for i=1:ntrials
    source_ts(i,:,:) = source_filter*squeeze(tlock.trial(i,:,:));
end

