clear; close all; clc;

[ret, hostname] = system('hostname');
if ret ~= 0
    hostnme = getenv('HOSTNAME');
end
hostname = strtrim(hostname);

if contains(hostname, 'hpc')
    % You are running your code on HPC
    curr_dir = pwd;
    path_parts = strsplit(curr_dir, filesep);
    fieldtrip_path = fullfile(filesep, path_parts{1:end-1}, 'fieldtrip');
    dataPath = fullfile(filesep, 'scratch', 'work', 'courses', ...
                    'PSYCH-GA-3405-2024fa');
elseif strcmp(hostname, 'sebastian_mac')
    % You are Sebastian
    % Setting up paths for data
    mydir = pwd;
    path_parts = strsplit(mydir, filesep);
    fieldtrip_path = fullfile(filesep, path_parts{1:end-1}, 'fieldtrip');
    idcs   = strfind(mydir, filesep);
    assert(strcmp(mydir(end-14:end), 'FieldtripMATLAB'), ...
        'Please cd to folder');
    dataPath = mydir(1:idcs(end-2)-1);
else
    % You are running on your device with Google drive path
    % Define the paths and initialize Fieldtrip
    my_user_id = 'mdd9787'; % change this to your netID
    curr_dir = pwd;
    path_parts = strsplit(curr_dir, filesep);
    base_dir = fullfile(filesep, path_parts{1:3});
    fieldtrip_path = fullfile(base_dir, 'Documents', 'MATLAB', 'fieldtrip-20220104');
    dataPath = fullfile(base_dir, 'Library', 'CloudStorage', ...
        strcat('GoogleDrive-', my_user_id, '@nyu.edu'), ...
        'My Drive', 'Coursework', 'EEG MEG methods', 'ClassData');
end

megRoot = fullfile(dataPath, 'MEGBids');

% Add fieldtrip to path
addpath(fullfile(fieldtrip_path));
ft_defaults;

%% Load data
groupName = 'GroupB'; % Change group name to analyze different dataset (valid IDs: GroupA, GroupC, GroupD)
userName = 'mrugank'; % Make sure to ensure that you are writing to your derivatives 
taskName = 'oddball';

% Get subject code for your group
groupMap = struct('GroupA', '001', 'GroupB', '002', 'GroupC', '003', 'GroupD', '004');
subjCode = groupMap.(groupName);

dataPath = fullfile(megRoot, ['sub-' subjCode], 'meg');
derivPath = fullfile(megRoot, 'derivatives', userName, 'preprocessing', ...
                    ['sub-' subjCode]);

if ~exist(derivPath,'dir'); mkdir(derivPath); end
saveRoot = ['sub-' subjCode '_task-' taskName '_'];

%%
cfg = [];
cfg.dataset = fullfile(dataPath, ...
                       ['sub-' subjCode '_task-' taskName '_eventmodel.sqd']);
cfg.continuous = 'yes';
data_raw = ft_preprocessing(cfg);

% save the data_raw (because it takes forever to load!)
% save(fullfile(derivPath, [saveRoot 'raw.mat']), 'data_raw', '-v7.3');

%% Visualize data
cfg = [];
cfg.channel = data_raw.label(1:30); % First 157 channels are data
cfg.blocksize = 40;
cfg.preproc.hpfilter = 'yes';
cfg.preproc.hpfreq = 1;
cfg.ylim = [ -10-13  10-13 ];
cfg_art = ft_databrowser(cfg, data_raw);

% save the cfg_art
% save(fullfile(derivPath, [saveRoot '_art.mat']), 'cfg_art', '-v7.3');

%% Bad channels
if strcmp(subjCode, '002')
    bad_channels = {'AG015', 'AG079'};
end

%% prepare a neighborhood structure for interpolation
%(we need to know what channels to use)
cfg = [];
cfg.method = 'triangulation';
cfg.feedback = 'yes';
cfg.grad = data_raw.grad;
neighbours = ft_prepare_neighbours(cfg);
%% Interpolate bad channels
cfg = [];
cfg.method = 'weighted';
cfg.neighbours = neighbours;
cfg.badchannel = bad_channels; % List your bad channels
data_interp = ft_channelrepair(cfg, data_raw);


%% compute one big ICA ...
n_sample = size(data_raw.trial{1},2);
artifact_def = cfg_art.artfctdef.visual.artifact;
art_vec = zeros(1, n_sample);
for aa = 1 : size(artifact_def,1)
    art_vec(artifact_def(aa,1):...
        artifact_def(aa,2)) = 1;
end
% get the rank of the data (is 240 here because we 
% excluded -15 channels and did avg reference -1)


%% Select only MEG data
cfg = [];
cfg.channel = data_interp.label(1:157);
data_interp = ft_selectdata(cfg, data_interp);
%% ... or with eye channels
rnk = rank(data_interp.trial{1});

% highpass filter for ica
cfg = [];
cfg.hpfilter = 'yes';
cfg.hpfreq = 1;
data_nan = ft_preprocessing(cfg, data_interp);

% mask artifacts for ICA
data_nan.trial{1}(:,art_vec==1)  = nan;

% compute the ICA on masked+filtered data
cfg = [];
cfg.method = 'runica';
cfg.numcomponent= rnk;
% cfg.unmixing = unmixing_eye; % if you have your unmixing matrix handy.
% cfg.topolabel = data_nan.label;
data_comp = ft_componentanalysis(cfg, data_nan);

%% add eye channels to layout

% add eye positions to layout
layeye = lay;
layeye.label = [lay.label; {'EOG-1'; 'EOG-2'; 'EOG-3'; 'EOG-4'}];
% layeye.pos = [lay.pos; [-0.4 0.28]; [-0.2  0.4]; [0.2 0.4]; [0.4 0.28]]; % Example positions
layeye.pos = [lay.pos; [-0.49 0.49]; [-0.48 0.48];  [0.48 0.48]; [0.49 0.49]]; % Example positions

layeye.width = [lay.width;ones(4, 1).*lay.width(1)];
layeye.height = [lay.width;ones(4, 1).*lay.width(1)];
cfg = []; cfg.layout = layeye;
ft_layoutplot(cfg)
%%
cfg = [];
cfg.layout = layeye;
% cfg.preproc.hpfilter = 'yes';
% cfg.preproc.hpfreq = 0.5;

cfg.viewmode = 'component';
cfg.ylim = [-100, 100];
ft_databrowser(cfg, data_comp);

%% reject components:
cfg           = [];
cfg.component = [1 3 10, 13];
data_clean    = ft_rejectcomponent(cfg,data_comp, ...
    data_combined);

% save(fullfile(derivPath, [saveRoot '_cleaned.mat']), 'data_clean', '-v7.3');
%% 
cfg = [];
cfg.layout = layeye;
cfg.viewmode = 'vertical';
cfg.viewmode = 'vertical';
cfg.preproc.hpfilter = 'yes';
cfg.preproc.hpfreq = 0.5;
cfg.ylim = [-10 10];
ft_databrowser(cfg, data_clean)

%% read in trigger info
event_info = ...
    readtable(...
    fullfile(dataPath,...
    ['sub-' subjCode '_task-oddball_events.tsv']), 'FileType', 'text');
n_events = size(event_info,1);

sum(event_info.trialNumber ~=0)
n_trials = sum(event_info.trialNumber ~=0 & ...
    event_info.value_fixed ~=16 ) ;
samp = nan(n_trials, 1);
info = cell(n_trials,1);
inf.group = '';
inf.type = '';
trial_idx = 0;

for ee = 1: n_events
    switch event_info.value_fixed(ee)
        case 1
            inf.group = 'tone';
        case 2
            inf.group = 'vis';
        case 4
            inf.group = 'aud';
        case 8
            inf.group = 'story_start';
        case 16
            inf.group = 'story_seg';
        case 64
            trial_idx = trial_idx+1;
            samp(trial_idx) =  ...
                event_info.sampleOnset(ee);

            info{trial_idx} = inf;
            info{trial_idx}.blockType = ...
                event_info.blockTypeActual(ee);
            info{trial_idx}.trialnumber = ...
                event_info.trialNumber(ee);
            info{trial_idx}.sampleinfo = ...
                event_info.sampleOnset(ee);
            info{trial_idx}.type  = 'even';
            info{trial_idx}.blocknumber = ...
                event_info.blockNumber(ee);
            info{trial_idx}.trigger =  ...
               event_info.value_fixed(ee);
        case 128
            trial_idx = trial_idx+1;
            samp(trial_idx) =  ...
                event_info.sampleOnset(ee);
            info{trial_idx} = inf;
            info{trial_idx}.blockType = ...
                event_info.blockTypeActual(ee);
            info{trial_idx}.trialnumber = ...
                event_info.trialNumber(ee);
            info{trial_idx}.sampleinfo = ...
                event_info.sampleOnset(ee);
            info{trial_idx}.type  = 'odd';
            info{trial_idx}.blocknumber = ...
                event_info.blockNumber(ee);
            info{trial_idx}.trigger =  ...
                event_info.value_fixed(ee);
            
    end

end
 
%% Epoch storySegments
cfg = [];
cfg.trials = find(...
    cellfun(@(x) ...
    strcmp('story_start', string(x.group)), ...
    data_clean.trialinfo));

tl_vis = ft_timelockanalysis(cfg, data_bl);

