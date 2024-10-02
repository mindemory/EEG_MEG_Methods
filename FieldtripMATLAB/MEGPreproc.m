clear; close all; clc;

% Initialize Fieldtrip
addpath('/Users/mrugankdake/Documents/MATLAB/fieldtrip-20220104');
ft_defaults;

%% Load data
% Define the paths
rootPath = '/Users/mrugankdake/Library/CloudStorage/GoogleDrive-mdd9787@nyu.edu/My Drive/Coursework/EEG MEG methods/ClassData/MEG';
groupName = 'GroupB';
recRoot = fullfile(rootPath, groupName, 'Recording');

lstFiles = dir(recRoot);
lstFiles = {lstFiles.name};

% Filter the files based on their extensions
sqdFiles = lstFiles(endsWith(lstFiles, '_NR.fif'));
% mrkFiles = lstFiles(endsWith(lstFiles, '.mrk'));
% elpFiles = lstFiles(endsWith(lstFiles, 'Points.txt'));
% hsFiles = lstFiles(endsWith(lstFiles, 'HS.txt'));

% In this case, second file is the one we want to load
ff = sqdFiles{1};
% mrk = cellfun(@(x) fullfile(recRoot, x), mrkFiles, 'UniformOutput', false);
% elp = fullfile(recRoot, elpFiles{1});
% hsp = fullfile(recRoot, hsFiles{1});
kitpath = fullfile(recRoot, ff);

%%
% Load the raw data using FieldTrip
cfg = [];
cfg.datafile = kitpath;
rawData = ft_preprocessing(cfg);

data = rawData.trial{1};
chLabel = rawData.label;

%% Visualize data and screen for artifacts
megChans = find(contains(chLabel, 'MEG'));
cfg = [];
cfg.viewmode = 'vertical';
cfg.channel = megChans;
cfg.ylim = [-1e-13, 1e-13];
ft_databrowser(cfg, rawData)