clear; close all; clc;

% Setting up paths for data
curr_dir = pwd;
path_parts = strsplit(curr_dir, filesep);
fieldtrip_path = fullfile(filesep, path_parts{1:end-1}, 'fieldtrip');
dataPath = fullfile(filesep, 'scratch', 'work', 'courses', ...
                'PSYCH-GA-3405-2024fa');

meg_path = fullfile(dataPath, 'MEG');
eeg_path = fullfile(dataPath, 'EEG');

% Add fieldtrip to path
addpath(fullfile(fieldtrip_path));
ft_defaults;

%% TODO: this needs to be BIDS
groupName = 'GroupD';
recRoot = fullfile(eeg_path, groupName, 'Recording');

lstFiles = dir(recRoot);
lstFiles = {lstFiles.name};

% Filter the files based on their extensions
bdfFiles = lstFiles(endsWith(lstFiles, '.bdf'));

% Load data
cfg = [];
cfg.dataset = fullfile(recRoot, bdfFiles{1});
data = ft_preprocessing(cfg);

% Visualize the data
% cfg = [];
% cfg.viewmode = 'vertical';
% cfg.preproc.hpfilter = 'yes';
% cfg.preproc.hpfreq = 0.5;
% ft_databrowser(cfg, data)

%% Epoch the data
cfg = [];
cfg.trialfun = 'ft_trialfun_general';
cfg.dataset = fullfile(recRoot, bdfFiles{1});
cfg_trl = ft_definetrial(cfg);
