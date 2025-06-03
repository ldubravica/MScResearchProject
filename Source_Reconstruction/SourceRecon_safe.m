function [source_ts,aals] = SourceRecon_safe(filename,filepath)
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
fieldtrip_folder = 'C:\code\MScResearchProject\fieldtrip-20250423'; % Change to the relavant FieldTrip Folder
addpath(fieldtrip_folder); % Adds FieldTrip to the MATLAB path
ft_defaults; % Initializes FieldTrip with ft_defaults

if ~exist('AAL.mat', 'file')
  error('AAL file missing.');
end
load('AAL.mat');

source_pos = AAL_centroids(select_aal_idx,:); % 3D coordinates (x,y,z) for each region
aals = AAL_Labels(select_aal_idx); % text labels (e.g., "Frontal_Sup_L")

%% Load data - Epoched MATLAB Format

% Load MATLAB file containing the pre-epoched data - assuming filename = 'your_file.mat'
load(fullfile(filepath, filename), 'epochs_open', 'epochs_closed', 'channel_names');
epochs = epochs_open; % or epochs_closed

% Create FieldTrip-compatible data structure manually
[n_trials, n_channels, n_times] = size(epochs);
fsample = 500
label = channel_names

% Build the FieldTrip 'raw' data structure manually
data = [];
data.label     = label;               % Cell array of channel labels
data.fsample   = fsample;
data.trial     = cell(1, n_trials);    % FieldTrip needs trial to be a cell array
data.time      = cell(1, n_trials);    % Time vector for each trial
data.sampleinfo = zeros(n_trials,2);   % (optional) sample start/end per trial

for i = 1:n_trials
    data.trial{i} = squeeze(epochs(i,:,:)); % [channels x times]
    data.time{i} = (0:n_times-1)/data.fsample;
    data.sampleinfo(i,:) = [1 n_times]; % could also adjust if needed
end

%% Prepare Beamforming
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

%% Source Localization
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

%% Reconstruction
% Use filters to reconstruct the source time series from EEG time series
source_filter = cell2mat(source.avg.filter);
[nsource,nelec] = size(source_filter);
[ntrials,nelec,T] = size(tlock.trial);
source_ts = zeros(ntrials,nsource,T);
for i=1:ntrials
    source_ts(i,:,:) = source_filter*squeeze(tlock.trial(i,:,:));
end

