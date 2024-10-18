clear; close all; clc;

% Setting up paths for data
mydir = pwd;
path_parts = strsplit(mydir, filesep);
fieldtrip_path = fullfile(filesep, path_parts{1:end-1}, 'fieldtrip');

if ~ismac && ~ispc
    dataPath = fullfile(filesep, 'scratch', 'work', 'courses', ...
        'PSYCH-GA-3405-2024fa');
    % TODO
    % reorganize such that: dataPath = mydir(1:idcs(end-1)-1);

else

    idcs   = strfind(mydir, filesep);
    assert(strcmp(mydir(end-14:end), 'FieldtripMATLAB'), ...
        'Please cd to folder');
    dataPath = mydir(1:idcs(end-2)-1);
end


meg_path = fullfile(dataPath, 'MEG');
eeg_path = fullfile(dataPath, 'EEG');

% Add fieldtrip to path
addpath(fullfile(fieldtrip_path));
ft_defaults;
%% list of subjects (couls be cell array)

sjs = ...
    ['001']; % todo change to make BIDS compliant
%%

acq = '';
code_ = sjs(1,:); % pick one subject
subj_pth_ = ...
    [dataPath filesep 'sub-' code_ filesep 'eeg'];
deriv_pth = [dataPath filesep ...
    'derivatives' filesep 'preprocessing', filesep ...
    'sub-' code_ ];
if ~exist(deriv_pth,'dir'); mkdir(deriv_pth); end

%% 

cfg = [];
cfg.dataset = fullfile(subj_pth_, ...
    ['sub-' code_ '.bdf']);
cfg.continuous = 'yes';
data_raw = ft_preprocessing(cfg);


%% read in electrodes with known coordinates
opts = delimitedTextImportOptions("NumVariables", 7);
% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = "\t";
% Specify column names and types
opts.VariableNames = ["name1", "x1", "y1", "z1", "Var5", "Var6", "Var7"];
opts.SelectedVariableNames = ["name1", "x1", "y1", "z1"];
opts.VariableTypes = ["string", "double", "double", "double", "string", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["name1", "Var5", "Var6", "Var7"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["name1", "Var5", "Var6", "Var7"], "EmptyFieldRule", "auto");
tbl = readtable(fullfile(subj_pth_, ['sub-' code_ '_electrodes.tsv']), opts);

%% Keep the electrode positions
elec = [];
elec.elecpos = [ tbl.x1,tbl.y1,tbl.z1];
n_elec = size(elec.elecpos,1);
elec.label = data_raw.label(1:size(elec.elecpos,1));
data_raw.elec = elec;
%% prepare the standard layout (we can do this from elec too but need to adjust scale and rotation)
cfg = [];
cfg.layout = 'BioSemi256'
cfg.feedback = 'yes';
lay = ft_prepare_layout(cfg);
%%
% Visualize the data to check for rogue channels
cfg = [];
cfg.layout = lay;
cfg.viewmode = 'vertical';
cfg.preproc.hpfilter = 'yes';
cfg.preproc.hpfreq = 0.5;
ft_databrowser(cfg, data_raw);
%% First step is we re-reference the data
cfg = [];
cfg.channel = data_raw.label(1:n_elec); % this is the default
cfg.reref = 'yes';
cfg.refmethod = 'avg';
cfg.refchannel = data_raw.label(1:n_elec);
data_avg = ft_preprocessing(cfg, data_raw);

%% now identify bad channels
% go back one step and interpolate
% then re-reference again!
% then epoch the data
%%
cfg = [];
cfg.layout = lay;
cfg.viewmode = 'vertical';
cfg.preproc.hpfilter = 'yes';
cfg.preproc.hpfreq = 0.5;
ft_databrowser(cfg, data_avg);
%% Epoch the data
cfg = [];
cfg.trialfun = 'ft_trialfun_general';
cfg.dataset = fullfile(recRoot, bdfFiles{1});
cfg_trl = ft_definetrial(cfg);
