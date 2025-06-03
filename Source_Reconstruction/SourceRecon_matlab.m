function [source_ts_open, source_ts_closed, aals] = SourceRecon_matlab(filename)

%% SourceRecon_matlab
% 
% Perform source reconstruction (beamforming) from epoched data stored in a MATLAB file.
%
% Inputs:
%   filename -- string, path to the MATLAB file containing:
%               - epochs_open [n_epochs × n_channels × n_times]
%               - epochs_closed [n_epochs × n_channels × n_times]
%               - channel_names {n_channels×1 cell}
%
% Outputs:
%   source_ts_open   -- [n_epochs_open × n_sources × n_times]
%   source_ts_closed -- [n_epochs_closed × n_sources × n_times]
%   aals             -- {n_sources×1 cell}, names of selected AAL regions
%
% Pedro Mediano, Hardik Rajpal, Modified 2025

%% Defining the ROIs to select
Frontal = [3:16, 19:20, 23:26];
Occipital = [43:54];
Parietal = [59:70];
Sensorimotor = [1,2,17,18,57,58];
Temporal = [81:90];
select_aal_idx = sort(cat(2,Frontal, Occipital, Parietal, Sensorimotor, Temporal));

%% Load AAL info
if ~exist('AAL.mat', 'file')
    error('AAL.mat file not found.');
end
load('AAL.mat'); % Should contain AAL_centroids and AAL_Labels
source_pos = AAL_centroids(select_aal_idx,:); % 3D coordinates (x,y,z) for each region
aals = AAL_Labels(select_aal_idx); % text labels (e.g., "Frontal_Sup_L")

%% Prepare Fieldtrip
fieldtrip_folder = 'C:\code\MScResearchProject\fieldtrip-20250423'; % Adjust path!
addpath(fieldtrip_folder); % Adds FieldTrip to the MATLAB path
ft_defaults; % Initializes FieldTrip with ft_defaults

%% Setup timelock and source analysis configuration
cfg_timelock = [];
cfg_timelock.covariance = 'yes';
cfg_timelock.covariancewindow = 'all';
cfg_timelock.channel = 'EEG';
cfg_timelock.keeptrials = 'yes';

% Load standard headmodel
vol = load([fieldtrip_folder, '/template/headmodel/standard_bem.mat']);
template_vol = vol.vol;

cfg_source = [];
cfg_source.method = 'lcmv';
cfg_source.sourcemodel.pos = source_pos;
cfg_source.headmodel = template_vol;
cfg_source.lcmv.keepfilter = 'yes';
cfg_source.lcmv.lambda = '5%';
cfg_source.lcmv.fixedori = 'yes';
cfg_source.elec = ft_read_sens('/template/electrode/standard_1020.elc'); 

%% Load data
data_struct = load(filename);
epochs_open = data_struct.epochs_open;
epochs_closed = data_struct.epochs_closed;
channel_names = data_struct.channel_names;
fs = 500; % Sampling frequency

%% Process Data
if ~isempty(epochs_open)
    % Convert epochs to FieldTrip format
    data_open = array_to_ft(epochs_open, channel_names, fs);
    % Compute covariance (timelock analysis)
    tlock_open = ft_timelockanalysis(cfg_timelock, data_open);
    % Source Analysis
    source_open = ft_sourceanalysis(cfg_source, tlock_open);
    % Reconstruct source time series
    source_ts_open = apply_spatial_filter(source_open, tlock_open);
else
    warning('epochs_open is empty. Skipping.');
    source_ts_open = [];
end

if ~isempty(epochs_closed)
    % Convert epochs to FieldTrip format
    data_closed = array_to_ft(epochs_closed, channel_names, fs);
    % Compute covariance (timelock analysis)
    tlock_closed = ft_timelockanalysis(cfg_timelock, data_closed);
    % Source Analysis
    source_closed = ft_sourceanalysis(cfg_source, tlock_closed);
    % Reconstruct source time series
    source_ts_closed = apply_spatial_filter(source_closed, tlock_closed);
else
    warning('epochs_closed is empty. Skipping.');
    source_ts_closed = [];
end

end

%% Helper function: Convert 3D array into FieldTrip data
function ft_data = array_to_ft(data, channel_names, fs)
    [n_epochs, n_channels, n_times] = size(data);
    ft_data = [];
    ft_data.label = channel_names;
    ft_data.fsample = fs;
    ft_data.trial = cell(1,n_epochs);
    ft_data.time = cell(1,n_epochs);
    ft_data.sampleinfo = zeros(n_epochs, 2);

    for i = 1:n_epochs
        ft_data.trial{i} = squeeze(data(i,:,:));
        ft_data.time{i} = (0:(n_times-1))/fs;
        ft_data.sampleinfo(i,:) = [(i-1)*n_times + 1, i*n_times]; % fabricate sample info
    end
end

%% Helper function: Apply spatial filter
function source_ts = apply_spatial_filter(source, tlock)
    source_filter = cell2mat(source.avg.filter);
    [n_sources, n_channels] = size(source_filter);
    [n_trials, ~, n_times] = size(tlock.trial);
    source_ts = zeros(n_trials, n_sources, n_times);
    for i = 1:n_trials
        source_ts(i,:,:) = source_filter * squeeze(tlock.trial(i,:,:));
    end
end
