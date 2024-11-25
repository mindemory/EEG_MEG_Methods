clear; close all; clc;

% Setting up paths for data
mydir = pwd;
path_parts = strsplit(mydir, filesep);
fieldtrip_path = fullfile(filesep, path_parts{1:end-1}, 'fieldtrip');

if ~ismac & ~ispc
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


% Add fieldtrip to path
addpath(fullfile(fieldtrip_path));
ft_defaults;
%% list of subjects (could be cell array)

meg_sjs = ...
    ['002']; %
eeg_sjs = ...
    ['004']; % todo change to make BIDS compliant
%% set all paths

acq = '';
eeg_code_ = eeg_sjs(1,:); % pick one subject
meg_code_ = meg_sjs(1,:); % pick one subject

anat_pth_eeg = ...
    [dataPath filesep 'sub-' eeg_code_ filesep 'anat'];

anat_pth_meg = ...
    [dataPath filesep 'sub-' meg_code_ filesep 'anat'];

meg_pth_ = ...
    [dataPath filesep 'sub-' meg_code_ filesep 'meg'];
eeg_pth_ = ...
    [dataPath filesep 'sub-' eeg_code_ filesep 'eeg'];

deriv_pth_eeg = [dataPath filesep ...
    'derivatives' filesep 'eeg_preprocessing', filesep ...
    'sub-' eeg_code_ ];

deriv_pth_meg = [dataPath filesep ...
    'derivatives' filesep 'meg_preprocessing', filesep ...
    'sub-' meg_code_ ];

if ~exist(deriv_pth_meg,'dir'); mkdir(deriv_pth_meg); end
if ~exist(deriv_pth_eeg,'dir'); mkdir(deriv_pth_eeg); end

hdm_pth_eeg = [dataPath filesep ...
    'derivatives' filesep 'eeg_headmodeling', filesep ...
    'sub-' eeg_code_ ];

hdm_pth_meg = [dataPath filesep ...
    'derivatives' filesep 'meg_headmodeling', filesep ...
    'sub-' meg_code_ ];

%% prepare for plots..MEG

% layout with YOKOGAWA is not defined in FT
% we can hack one together:

hdr = ft_read_header(fullfile(meg_pth_, ...
    ['sub-' meg_code_ '_task-oddball_eventmodel.sqd']));
% read the position of the sensors from the data
grad                        = ft_read_sens(fullfile(meg_pth_, ...
    ['sub-' meg_code_ '_task-oddball_eventmodel.sqd'])); % this can be inspected with ft_plot_sens(grad)

grad_mm = ft_convert_units(grad, 'mm');
% prepare the custom channel layout
cfg                         = [];
cfg.grad                    = grad;
lay                      = ft_prepare_layout(cfg);
sel                         = 1:(length(lay.label)-2); % the last two are COMNT and SCALE

% scale & stretch the position of the sensors
lay.pos(sel,:)           = lay.pos(sel,:) * 1.05;
lay.pos(sel,2)           = lay.pos(sel,2) * 1.08 + 0.02;

% load the CTF151 helmet and mask
cfg                         = [];
cfg.layout                  = 'CTF151_helmet';
ctf151                      = ft_prepare_layout(cfg);

% use the CTF151 outlint and mask instead of the circle
lay.outline              = ctf151.outline;
lay.mask                 = ctf151.mask;


% %% load all information
% load([deriv_pth_eeg filesep 'sub-'...
%     eeg_code_ '_task-oddball_epochedStory.mat']);
% data_eeg = epochStory;
% 
% % data_eeg.time{1} = data_eeg.time{1} - 8;
load([deriv_pth_meg filesep 'sub-'...
    meg_code_ '_task-oddball_epochedStory.mat']);
data_meg = epochStory;
% % data_meg.time{1} = data_meg.time{1} - 8;
% %% load artifact information
% load([deriv_pth_eeg filesep 'sub-'...
%     eeg_code_ '_task-oddball__art.mat']);
% art_eeg = cfg_art.artfctdef.visual.artifact;
% art_eeg(art_eeg(:,2) < data_eeg.sampleinfo(1) |...
%     art_eeg(:,1) > data_eeg.sampleinfo(2),:) = [];
load([deriv_pth_meg filesep 'sub-'...
    meg_code_ '_task-oddball__art.mat']);
art_meg = cfg_art.artfctdef.visual.artifact;
art_meg(art_meg(:,2) < data_meg.sampleinfo(1) |...
    art_meg(:,1) > data_meg.sampleinfo(2),:) = [];
% 
% 
% %% data browse this to check
cfg = [];
cfg.blocksize = 10;
cfg.ylim = [-10 10];
cfg.preproc.hpfilter ='yes';
cfg.preproc.hpfreq = 1;
cfg.viewmode = 'vertical';
% ft_databrowser(cfg, data_eeg);
% 
cfg.ylim = [ -9.6e-13     9.6e-13 ];
ft_databrowser(cfg, data_meg);

%% prepare the data (add a lowpass filter)
cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = [1 55];
cfg.detrend = 'yes';
data_meg = ft_preprocessing(cfg, data_meg);  
% %%

%% now compare to the MEG story reconstruction

% load([hdm_pth_meg filesep 'grid_indiv.mat']);
load([hdm_pth_meg filesep 'vol.mat']);
% vol = vol_improved;
load([hdm_pth_meg filesep 'mri_reslice.mat']);
%% make a new very coarse grid:

% cfg = [];
% cfg.grid.resolution = 35; % in mm
% cfg.grid.unit = 'mm';
% cfg.tight = 'yes'; % ensures grid points are inside the brain
% cfg.inwardshift = 5; % shift grid points slightly inward to ensure they are inside the brain
% cfg.headmodel = vol;
% grid = ft_prepare_sourcemodel(cfg);

%% actually make a coarse standard grid
ftpath = fileparts(which('ft_defaults'));
%prepare a template grid with inward shift and warp it to this persons space
template = load(fullfile(ftpath, 'template/headmodel/standard_singleshell'));
template_vol = ft_convert_units(template.vol, 'mm');
% Prepare the source model grid
cfg = [];
cfg.grid.resolution = 30; % in mm
cfg.grid.unit = 'mm';
cfg.tight = 'yes'; % ensures grid points are inside the brain
cfg.inwardshift = 5; % shift grid points 5mm inward
cfg.headmodel = template.vol;
template_grid = ft_prepare_sourcemodel(cfg);
%% plot template grid
ft_plot_mesh(template_vol.bnd, 'facecolor', 'k',...
    'facealpha', 0.2); % brain
hold on;
ft_plot_mesh(template_grid.pos(template_grid.inside,:), ...
    'vertexcolor', 'r');
pos_mat = template_grid.pos(template_grid.inside,:);
%% now warp the new template grid to individual head space:
% use the resliced MRI and warp the grid to it. 
cfg = [];
cfg.warpmni             = 'yes';
cfg.template            = template_grid;
cfg.nonlinear           = 'yes'; % use non-linear normalization
cfg.mri                 = mri_reslice;
cfg.sourcemodel.unit    = 'mm';
grid = ft_prepare_sourcemodel(cfg);   % ras + 
%% alignment plot
figure,ft_plot_ortho(mri_reslice.anatomy, 'style','intersect','transform',mri_reslice.transform)
ft_plot_mesh(grid.pos(grid.inside,:),'vertexcolor',[1,0,0])
view([80,-15,0])
hold on;
ft_plot_mesh(vol.bnd, 'facecolor', ...
    'y', 'facealpha', 0.2); % brain

ft_plot_sens(grad_mm,'facecolor', 'g');
%% compute sources but mask artifacts with nan

data_nan = data_meg;
samplevec1 = data_nan.sampleinfo(1,1):...
    data_nan.sampleinfo(1,2);
samplevec2 = data_nan.sampleinfo(2,1):...
    data_nan.sampleinfo(2,2);
sample_vec_all = [samplevec1, samplevec2];


artifact_all = zeros(size(sample_vec_all));
for art = 1: size(art_meg,1)
    artifact_all(find(sample_vec_all==art_meg(art,1)):...
        find(sample_vec_all==art_meg(art,2))) = 1;
end


data_nan.trial{1}(:, find(...
    artifact_all(1:length(samplevec1)))) = nan;

data_nan.trial{2}(:, find(...
    artifact_all(length(samplevec2)+1:end))) = nan;

cfg = [];
cfg.covariance = 'yes';
cfg.covariancewindow = 'all'; % Use the entire time window
timelock = ft_timelockanalysis(cfg, data_nan);
clear data_nan
%% LCMV computes leadfield on the fly (use one common filter)
cfg                 = [];
cfg.method          = 'lcmv';
cfg.sourcemodel     = grid;
cfg.headmodel       = vol;
cfg.grad            = grad_mm;
cfg.keepleadfield   = 'yes';
cfg.lcmv.keepfilter = 'yes';
cfg.lcmv.fixedori   = 'yes';
cfg.lcmv.lambda     = '1%';  % must be '1%' not 1
cfg.normalize  = 'yes'; % normalize LF to attenuate depth bias
cfg.normalizeparam  = 0.5; % default = 0.5
source = ft_sourceanalysis(cfg, timelock);
inside_pos = find(source.inside);

% Extract all filters for inside voxels
W_meg = cell2mat(source.avg.filter(source.inside));
%% create sourcedata
cfg =[];
cfg.resamplefs = 100;
data_meg_rs = ft_resampledata(cfg, data_meg);

source_trial1 = W_meg*data_meg_rs.trial{1};
source_trial2 = W_meg*data_meg_rs.trial{2};
%% make a source data structure 
sourcedata.label = cell(numel(inside_pos), 1); % Initialize with empty cells

% Assign labels only to inside sources
for i = 1:numel(inside_pos)

    sourcedata.label{i} = sprintf('S_%d',...
        inside_pos(i));
end

sourcedata.trial{1} = zscore(source_trial1,0,2);
sourcedata.trial{2} =  zscore(source_trial2,0,2);
sourcedata.time = data_meg_rs.time;
cfg = [];
cfg.latency = [-1 800];

sourcedata = ft_selectdata(cfg, sourcedata);

%%
% rmpath(genpath('MVGC1'));

regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'OLS';  % information criteria regression mode ('OLS', 'LWR' or empty for default)
momax     = 20;     % maximum model order for model order is 20, which corresponds to 2s


% concatenate the trials
X = cat(3,sourcedata.trial{1},sourcedata.trial{2});


% maximally 1 second
[AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X,momax,icregmode);
% 20 and 8
morder = moBIC;

% finally estimate Granger causality in the traditional way i.e in
% the time domain
[F1,A1,SIG1, E1] = ...
    GCCA_tsdata_to_pwcgc(sourcedata.trial{1},morder,regmode); % use same model order for reduced as for full regressions

[F2,A2,SIG2, E2] = ...
    GCCA_tsdata_to_pwcgc(sourcedata.trial{2},morder,regmode); % use same model order for reduced as for full regressions

% plot the matrices
subplot(131);
imagesc(F1);

subplot(132);
imagesc(F2)

F3 = F2-F1;

subplot(133);
imagesc(F3);

%%

%%  make the colormap with the right number of entries
% Define the original blue to red colormap
originalColormap = [0 0 1; 1 0 0]; % Blue to Red

% Define the number of entries you want in the new colormap
numEntries = sum(grid.inside);

% Generate the new colormap
x = linspace(1, size(originalColormap, 1), numEntries);
newColormap = interp1(1:size(originalColormap, 1), originalColormap, x);

% Apply the new colormap
colormap(newColormap);

% Display the colormap
colorbar;
%% save
save F1 F1
save F2 F2
save F3 F3
save A1 A1
save A2 A2
save SIG1 SIG1
save SIG2 SIG2
pos_mat = template_grid.pos(template_grid.inside,:);
save pos_mat pos_mat
save newColormap newColormap

