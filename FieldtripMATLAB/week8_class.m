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

eegRoot = fullfile(dataPath, 'EEGBids');

% Add fieldtrip to path
addpath(fullfile(fieldtrip_path));
ft_defaults;

%% Load data
groupName = 'GroupA'; % Change group name to analyze different dataset (valid IDs: GroupA, GroupC, GroupD)
userName = 'mrugank'; % Make sure to ensure that you are writing to your derivatives 
taskName = 'oddball';

% Get subject code for your group
groupMap = struct('GroupA', '001', 'GroupC', '003', 'GroupD', '004');
subjCode = groupMap.(groupName);

dataPath = fullfile(eegRoot, ['sub-' subjCode], 'eeg');
derivPath = fullfile(eegRoot, 'derivatives', userName, 'preprocessing', ...
                    ['sub-' subjCode]);

if ~exist(derivPath,'dir'); mkdir(derivPath); end
saveRoot = ['sub-' subjCode '_task-' taskName '_'];

%% Load data
cfg = [];
cfg.dataset = fullfile(dataPath, ...
                       ['sub-' subjCode '_task-' taskName '_eeg.bdf']);
cfg.continuous = 'yes';
data_raw = ft_preprocessing(cfg);


%% read in electrodes with known coordinates
elecPath = fullfile(dataPath, ['sub-' subjCode '_space-CapTrak_electrodes.tsv']);
elecTbl = readtable(elecPath, 'FileType', 'text');
elecTbl = rmmissing(elecTbl); % Remove NaN (for EXG channels)


%% Keep the electrode positions
elec = [];
elec.elecpos = [elecTbl.x, elecTbl.y, elecTbl.z];
elec.chanpos = elec.elecpos;
elec.unit = 'mm';
n_elec = size(elec.elecpos,1);
n_sample = size(data_raw.trial{1},2);
SR = data_raw.fsample;
elec.label = data_raw.label(1:size(elec.elecpos,1));
data_raw.elec = elec;
%% prepare the standard layout (we can do this from elec too but need to adjust scale and rotation)
cfg = [];
cfg.layout = 'BioSemi256';
cfg.feedback = 'yes';
lay = ft_prepare_layout(cfg);
%%
% Visualize the data to check for rogue channels
cfg = [];
cfg.layout = lay;
cfg.viewmode = 'vertical';
cfg.preproc.hpfilter = 'yes';
cfg.preproc.hpfreq = 0.5;
cfg.ylim = [-100 100];
cfg.blocksize = 10;
ft_databrowser(cfg, data_raw);
%% First step is we re-reference the data
cfg = [];
cfg.channel = data_raw.label(1:n_elec); % this is the default
cfg.reref = 'yes';
cfg.refmethod = 'avg';
cfg.refchannel = data_raw.label(1:n_elec);
data_avg = ft_preprocessing(cfg, data_raw);

%%
cfg = [];
cfg.layout = lay;
cfg.viewmode = 'vertical';
cfg.preproc.hpfilter = 'yes';
cfg.preproc.hpfreq = 0.5;
cfg.ylim = [-10 10];
cfg.blocksize = 10;
ft_databrowser(cfg, data_avg);
%% now identify bad channels
if strcmp(subjCode, '004')
    %'B7, 9, 25?, C18?, 'G30'? !F2 E5
    bad_channels = {'A26', 'B2', 'B3', 'B4', 'B7' 'B8',...
        'B19', 'C9', 'C30','D15','E5', 'F2' ,'F22',...
        'G25', 'G29'};
elseif strcmp(subjCode, '001')
    % Added by Mrugank
    % D22-32 & E7-16 & E23-32 & F6-15 & F25-32 (blinks)
    bad_channels = {'B1', 'B2', 'B8', 'B29', 'B30', ...
        'C6', 'C7', 'C11', 'C18', 'C19', 'C32', 'D7', 'D8', 'D9', 'D10', ...
        'D13', 'D14', 'D21', 'D27', 'G30', 'H1', 'H5', 'H11', 'H12', 'H14', 'H15', 'H16', ...
        };
end
% go back one step and interpolate
% then re-reference again!
% then epoch the data

%% prepare a neighborhood structure for interpolation
%(we need to know what channels to use)
cfg = [];
cfg.method = 'triangulation';
cfg.feedback = 'yes';
cfg.elec = elec;
neighbours = ft_prepare_neighbours(cfg);
%% Interpolate bad channels
cfg = [];
cfg.method = 'weighted';
cfg.neighbours = neighbours;
cfg.badchannel = bad_channels; % List your bad channels
data_interp = ft_channelrepair(cfg, data_raw);

% Save data interp
save(fullfile(derivPath, [saveRoot 'interp.mat']), ...
    'data_interp', '-v7.3');

%% redo the average reference
cfg = [];
cfg.channel = data_raw.label(1:n_elec); 
cfg.reref = 'yes';
cfg.refmethod = 'avg';
cfg.refchannel = data_raw.label(1:n_elec);
data_avg = ft_preprocessing(cfg, data_interp);

% compute the average here
avg_vec = mean(data_interp.trial{1}(1:n_elec,:));
% clean up
clear data_interp;
%% check the data again to see if we should exclude 
% more channels
% if so.. go back to previous step
cfg = [];
cfg.layout = lay;
cfg.viewmode = 'vertical';
cfg.preproc.hpfilter = 'yes';
cfg.preproc.hpfreq = 0.5;
cfg.ylim = [-10 10];
cfg.blocksize = 10;
ft_databrowser(cfg, data_avg);
%% concatenate with the EOG;
cfg = [];
cfg.channel = data_raw.label(n_elec+1:n_elec+4);

data_eye = ft_selectdata(cfg, data_raw);

% because eog channels are recorded with amplifier refrence, we need to
% re-reference them as well
data_eye.trial{1} = data_eye.trial{1} - avg_vec;
% EOG channels are recognized by databrowser via ft_channelselection
data_eye.label = {'EOG-1'; 'EOG-2'; 'EOG-3';'EOG-4'};
data_combined = ft_appenddata([],data_avg, data_eye);


%% now manually mark artifacts
cfg = [];
cfg.layout = lay;
cfg.viewmode = 'vertical';
cfg.preproc.hpfilter = 'yes';
cfg.preproc.hpfreq = 0.5;
cfg.ylim = [-10 10];
cfg.blocksize = 10;
cfg.eogscale = 0.8; % Four channels are eye channels
cfg_art = ft_databrowser(cfg, data_combined);

save(fullfile(derivPath, ['sub-' subjCode '_task-' taskName '_cfg_art.mat']), 'cfg_art', '-v7.3')
%% compute one big ICA ...
artifact_def = cfg_art.artfctdef.visual.artifact;
art_vec = zeros(1, n_sample);
for aa = 1 : size(artifact_def,1)
    art_vec(artifact_def(aa,1):...
        artifact_def(aa,2)) = 1;
end
% get the rank of the data (is 240 here because we 
% excluded -15 channels and did avg reference -1)

%% ... without eye channels

rnk = rank(data_avg.trial{1});

% highpass filter for ica
cfg = [];
cfg.hpfilter = 'yes';
cfg.hpfreq = 1;
data_nan = ft_preprocessing(cfg, data_avg);

% mask artifacts for ICA
data_nan.trial{1}(:,art_vec==1)  = nan;

% compute the ICA on masked+filtered data
cfg = [];
cfg.method = 'runica';
cfg.numcomponent= rnk;
data_comp = ft_componentanalysis(cfg, data_nan);

cfg = [];
cfg.layout = lay;
cfg.viewmode = 'component';
ft_databrowser(cfg, data_comp);
% remove 1 5? 7

%% ... or with eye channels
rnk = rank(data_combined.trial{1});

% highpass filter for ica
cfg = [];
cfg.hpfilter = 'yes';
cfg.hpfreq = 1;
data_nan = ft_preprocessing(cfg, data_combined);

% mask artifacts for ICA
data_nan.trial{1}(:,art_vec==1)  = nan;

% compute the ICA on masked+filtered data
cfg = [];
% cfg.method = 'runica';
cfg.method = 'fastica';
cfg.numcomponent= rnk;
data_comp = ft_componentanalysis(cfg, data_nan);

unmixing_eye = data_comp.unmixing;  
% 
% % Save the unmixing matrix
save(fullfile(derivPath, ['sub-' subjCode '_task-' taskName '_unmixing_eye.mat']), 'unmixing_eye', '-v7.3')
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

% reject 1, 3(pulse),  10?hz-eye and 13
%% use unmixing matrix on unfiltered data
cfg = [];
cfg.unmixing = unmixing_eye;
cfg.topolabel = data_combined.label;
data_comp = ft_componentanalysis(cfg, data_combined);
%% reject components:
cfg           = [];
cfg.component = [1 3 10, 13];
data_clean    = ft_rejectcomponent(cfg,data_comp, ...
    data_combined);
%% 
cfg = [];
cfg.layout = layeye;
cfg.viewmode = 'vertical';
cfg.viewmode = 'vertical';
cfg.preproc.hpfilter = 'yes';
cfg.preproc.hpfreq = 0.5;
cfg.ylim = [-10 10];
ft_databrowser(cfg, data_clean)

save(fullfile(derivPath, [saveRoot 'data_clean.mat']), 'data_clean', '-v7.3');
clear data_comp data_combined data_raw data_avg
%%
%% read in trigger info
event_info = ...
    getEvents(...
    fullfile(dataPath,...
    ['sub-' subjCode '_task-oddball_events.tsv']));
n_events = size(event_info,1);

sum(event_info.trialNumber ~=0)
n_trials = sum(event_info.trialNumber ~=0 & ...
    event_info.trigger_codes ~=16 ) ;
samp = nan(n_trials, 1);
info = cell(n_trials,1);
inf.group = '';
inf.type = '';
trial_idx = 0;

for ee = 1: n_events
    switch event_info.trigger_codes_fixed(ee)
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
                event_info.blockType(ee);
            info{trial_idx}.trialnumber = ...
                event_info.trialNumber(ee);
            info{trial_idx}.sampleinfo = ...
                event_info.sampleOnset(ee);
            info{trial_idx}.type  = 'even';
            info{trial_idx}.blocknumber = ...
                event_info.blockNumber(ee);
            info{trial_idx}.trigger =  ...
               event_info.trigger_codes_fixed(ee);
        case 128
            trial_idx = trial_idx+1;
            samp(trial_idx) =  ...
                event_info.sampleOnset(ee);
            info{trial_idx} = inf;
            info{trial_idx}.blockType = ...
                event_info.blockType(ee);
            info{trial_idx}.trialnumber = ...
                event_info.trialNumber(ee);
            info{trial_idx}.sampleinfo = ...
                event_info.sampleOnset(ee);
            info{trial_idx}.type  = 'odd';
            info{trial_idx}.blocknumber = ...
                event_info.blockNumber(ee);
            info{trial_idx}.trigger =  ...
                event_info.trigger_codes_fixed(ee);
            
    end

end
% NOTE: blocktype is correct
% legth of story is 742s 

%% overwrite trialinfo:

%%

% ?? maybe filter?
% cfg = [];
% cfg.hpfilter = 'yes';
% cfg.hpfreq = 0.5;
% data_clean = ft_preprocessing(cfg, data_clean);
% 
% cut into 3 second trials (lots of ovelap)
cfg = [];
cfg.trl = [samp-1*SR, samp+2*SR, ...
    repmat(-1*SR, [n_trials,1])];
data_cut = ft_redefinetrial(cfg, data_clean);
clear data_clean;
%%
% also remove eye channels
cfg = [];
cfg.channel = data_cut.label(1:n_elec);
data_cut = ft_selectdata(cfg, data_cut);
%% inspect one more time
cfg = [];
cfg.preproc.detrend = 'yes'; % now we detrend instead of highpass filtering! 
cfg.viewmode = 'vertical';
cfg.ylim = [-10 10];
ft_databrowser(cfg, data_cut);
%% overwrite trialingo
data_cut.trialinfo = info;

%% now we do "rough" preprocessing

% instead, we should define artifacts manually and check for overlap with
% the sampleinfo!

cfg = [];
cfg.preproc.detrend = 'yes'; % now we detrend instead of highpass filtering! 
cfg.method = 'summary';
data_pruned = ft_rejectvisual(cfg, data_cut);


% automatic muscle artifact detection
cfg = [];
cfg.metric = 'zvalue'
cfg.preproc.bpfilter    = 'yes';
cfg.preproc.bpfreq      = [70 90];% [110 140] this is EU setting
cfg.preproc.bpfiltord   =  8;
cfg.preproc.bpfilttype  = 'but';
cfg.preproc.rectify     = 'yes';
cfg.preproc.boxcar      = 0.2;
data_pruned = ft_rejectvisual(cfg, data_pruned);

%% processing for ERP

% TODO select only 
cfg = [];
cfg.detrend = 'no';
cfg.lpfilter = 'yes';
cfg.lpfreq = 15;
cfg.detrend = 'yes';
cfg.demean              = 'yes';    % we demean (baseline correct) ...
cfg.baselinewindow      = [-0.5 0]; % using the mean activity in this window
%      = ft_preprocessing(cfg, data_pruned);
%%
cfg = [];
cfg.trials = find(...
    cellfun(@(x) ...
    strcmp('semanticVis', string(x.blockType)), ...
    data_bl.trialinfo));

tl_vis = ft_timelockanalysis(cfg, data_bl);
cfg.trials = find(...
    cellfun(@(x) ...
    strcmp('semanticAud', string(x.blockType)), ...
    data_bl.trialinfo));
tl_aud = ft_timelockanalysis(cfg, data_bl);


% 
% cfg.keeptrials = 'yes';
% tl_kt = ft_timelockanalysis(cfg, data_bl);;

avg_new = trimmean(tl_kt.trial, 20);

tl.avg      = squeeze(avg_new);
% keep as template

% do manual trimmean and overwrite

% compute lda between conditions


figure;
cfg =[];
cfg.layout = lay;
ft_multiplotER(cfg, tl_vis, tl_aud);


cfg =[];
cfg.layout = lay;
ft_multiplotER(cfg, tl_odd, tl_even);