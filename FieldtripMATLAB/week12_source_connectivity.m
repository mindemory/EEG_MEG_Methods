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
%% prepare for plots..EEG
cfg = [];
cfg.layout = 'BioSemi256';
lay_eeg = ft_prepare_layout(cfg);

cfg_topo= {'mask' lay_eeg.mask, 'datmask' [], 'interplim' 'mask'};
cfg_lay = { 'point' 'no' 'box' 'no' 'label' 'no' };
n_eeg_channels = 256;
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

% 
% %% correlate the first and second story in source space
% source_trial1 = W_eeg*data_eeg.trial{1};
% source_trial2 = W_eeg*data_eeg.trial{2};
% 
% % get the trial data from pre and post
% datpre = [source_trial1(:, data_eeg.time{1}<0),source_trial2(:, data_eeg.time{2}<0)];
% % get the trial data from pre and post
% datpost = [source_trial1(:, data_eeg.time{1}>0 & ...
%     data_eeg.time{1}<700), source_trial2(:, data_eeg.time{2}>0 & ...
%     data_eeg.time{2}<700)];
% 
% 
% post_var = var(datpost, 0,2);
% pre_var = var(datpre, 0,2);
% % 
% 
% effect  = (post_var-pre_var)./(post_var+pre_var);
% source.story = (nan(size(source.inside)));
% source.story(source.inside) = effect;
% 
% %% plot
% cfg = [];
% cfg.parameter = {'story'};
% [interp] = ft_sourceinterpolate(cfg, source, mri_reslice);
% % plot the correlation
% 
% cfg = [];
% cfg.method        = 'ortho';
% cfg.crosshair = 'no';
% cfg.funparameter  = 'story';
% cfg.funcolormap = 'jet';
% cfg.funcolorlim   = [-0.4 0.4];
% % cfg.location = [0 -59 15];
% ft_sourceplot(cfg, interp);
%% now compare to the MEG story reconstruction

load([hdm_pth_meg filesep 'grid_indiv.mat']);
load([hdm_pth_meg filesep 'vol.mat']);
% vol = vol_improved;
load([hdm_pth_meg filesep 'mri_reslice.mat']);

%% alignment plot
figure,ft_plot_ortho(mri_reslice.anatomy, 'style','intersect','transform',mri_reslice.transform)
ft_plot_mesh(grid.pos(grid.inside,:),'vertexcolor',[1,0,0])
view([80,-15,0])
hold on;
ft_plot_mesh(vol.bnd, 'facecolor', ...
    'r', 'facealpha', 0.2); % brain

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
%% localize the story variance increase
cfg =[];
cfg.resamplefs = 250;
data_meg_rs = ft_resampledata(cfg, data_meg);
% source_trial1 = sourcemodel.trial{1};
% source_trial2 = sourcemodel.trial{2};
seg_length = 10*60*250 + 4*250; % 10 minutes of data?
source_trial1 = W_meg*data_meg_rs.trial{1};
source_trial2 = W_meg*data_meg_rs.trial{2};

% get the trial data from pre and post
datpre =...
    [source_trial1(...
    :, data_meg_rs.time{1}<0),source_trial2(...
    :, data_meg_rs.time{2}<0)];
% get the trial data from pre and post
% datpost = [source_trial1(:, data_meg_rs.time{1}>0 & ...
%     data_meg_rs.time{1}<700),...
%     source_trial2(:, data_meg_rs.time{2}>0 & ...
%     data_meg_rs.time{2}<700)];
sel = data_meg_rs.time{1}>0;
datpost = [source_trial1(:, sel(1:seg_length)),...
    source_trial2(:, sel(1:seg_length))];

post_var = var(datpost, 0,2);
pre_var = var(datpre, 0,2);
% 

effect  = (post_var-pre_var)./(post_var+pre_var);
source.story = (nan(size(source.inside)));
source.story(source.inside) = effect;

%% plot
cfg = [];
cfg.parameter = {'story'};
[interp] = ft_sourceinterpolate(cfg, source, mri_reslice);
% plot the correlation

cfg = [];
cfg.method        = 'ortho';
cfg.crosshair = 'yes';
cfg.funparameter  = 'story';
cfg.funcolormap = 'jet';
% cfg.funcolorlim   = [0 0.06];
cfg.location = [2 38 48];
ft_sourceplot(cfg, interp);
%% make a source data structure and filter in the delta band
sourcedata.label = cell(numel(inside_pos), 1); % Initialize with empty cells

% Assign labels only to inside sources
for i = 1:numel(inside_pos)

    sourcedata.label{i} = sprintf('S_%d',...
        inside_pos(i));
end

sourcedata.trial{1} = source_trial1;
sourcedata.trial{2} = source_trial2;
sourcedata.time = data_meg_rs.time;

cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = [2 5];
sourcedata = ft_preprocessing(cfg, sourcedata);

for tt = 1: numel(sourcedata.label)
    sourcedata.trial{1}(tt,:) = (hilbert(...
        sourcedata.trial{1}(tt,:)));
    sourcedata.trial{2}(tt,:) = (hilbert(...
        sourcedata.trial{2}(tt,:)));
end
%% try the variance increase again
sourceD_trial1 = abs(sourcedata.trial{1});
sourceD_trial2 = abs(sourcedata.trial{2});

% get the trial data from pre and post
datpre =...
    [sourceD_trial1(...
    :, data_meg_rs.time{1}<0),sourceD_trial2(...
    :, data_meg_rs.time{2}<0)];

datpost = [sourceD_trial1(:, sel(1:seg_length)),...
    sourceD_trial2(:, sel(1:seg_length))];

post_var = var(datpost, 0,2);
pre_var = var(datpre, 0,2);
% 

effect  = (post_var-pre_var)./(post_var+pre_var);
source.story = (nan(size(source.inside)));
source.story(source.inside) = effect;

%% plot
cfg = [];
cfg.parameter = {'story'};
[interp] = ft_sourceinterpolate(cfg, source, mri_reslice);
% plot the correlation

cfg = [];
cfg.method        = 'ortho';
cfg.crosshair = 'yes';
cfg.funparameter  = 'story';
cfg.funcolormap = 'jet';
% cfg.funcolorlim   = [0 0.06];
cfg.location = [2 38 48];
ft_sourceplot(cfg, interp);
%% try again with delta power phase coherence over time
% in first 10 minutes

source_phase1 = (sourcedata.trial{1}(:,1:seg_length))./...
    abs(sourcedata.trial{1}(:,1:seg_length));
source_phase2 = (sourcedata.trial{2}(:,1:seg_length))./...
    abs(sourcedata.trial{2}(:,1:seg_length));


delta_pc = abs(mean(source_phase1.*conj(source_phase2),2));


source.retest = (nan(size(source.inside)));
source.retest(source.inside) = delta_pc;

%% plot
close all;
cfg = [];
cfg.parameter = {'retest'};
[interp] = ft_sourceinterpolate(cfg, source, mri_reslice);
% plot the correlation

cfg = [];
cfg.method        = 'ortho';
cfg.crosshair = 'yes';
cfg.funparameter  = 'retest';
cfg.funcolormap = 'jet';
% cfg.funcolorlim   = [0 0.06];
% cfg.location = [2 38 48];
ft_sourceplot(cfg, interp);

%% pick peak source position in auditory cortex and compute simple correlation

[~, aud_index] =...
    max(source.retest(source.inside) );
% 
% [~, aud_index] =...
%     min(sum(abs(source.pos(inside_pos,:) - [9 59 29]),2));

source_tc1 = source_trial1(aud_index, :);
source_tc2 = source_trial2(aud_index, :);

%% compute correlation of auditory source with all others
source_corrs = zeros(numel(inside_pos),1);

for ii = 1: size(source_corrs, 1)
    ii
     source_corrs(ii) = (corr(source_trial1(ii,:)', ...
        source_tc1','Type', 'Spearman')+...
        corr(source_trial2(ii,:)', ...
        source_tc2','Type', 'Spearman'))/2;
end


source.corrs = nan(size(source.inside));
source.corrs(inside_pos) = source_corrs;
%% plot
cfg = [];
cfg.parameter = {'corrs'};
[interp] = ft_sourceinterpolate(cfg, source, mri_reslice);
% plot the correlation

cfg = [];
cfg.method        = 'ortho';
cfg.crosshair = 'yes';
cfg.funparameter  = 'corrs';
cfg.funcolormap = 'jet';
cfg.funcolorlim   = [0 0.8];
ft_sourceplot(cfg, interp);

%% coherence in delta band
% pick 10 minutes of story
source_trial1c = sourcedata.trial{1}(:,1:seg_length);
source_trial2c = sourcedata.trial{2}(:,1:seg_length);

% phase differences with the selected source
coh1 = source_trial1c * ...
    (source_trial1c(aud_index, :)');% complex conjugate transpose!!!
coh2 = source_trial2c *...
    (source_trial2c(aud_index, :)');% complex conjugate transpose!!!

% for correct scaling we should divide by n
coh1 = coh1./length(source_trial1c);
coh2 = coh2./length(source_trial2c);

source.coh = nan(size(source.inside));
source.coh(inside_pos) = (abs(coh1) + abs(coh2))./2;
% imaginary part:
source.icoh = nan(size(source.inside));
source.icoh(inside_pos) = (imag(coh1) + imag(coh2))./2;

%% compute phase coherence as well:
% normalize to unit amplitude
clear source_phase*
source_trial1c = source_trial1c./abs(source_trial1c);
source_trial2c = source_trial2c./abs(source_trial2c);

% coherence is now automatically phase coherence
pcoh1 = (source_trial1c * ...
    (source_trial1c(aud_index, :)'));% complex conjugate transpose!!!
pcoh2 = (source_trial2c *...
    (source_trial2c(aud_index, :)'));% complex conjugate transpose!!!


% for correct scaling we should divide by n
pcoh1 = pcoh1./length(source_trial1c);
pcoh2 = pcoh2./length(source_trial2c);

source.pcoh = nan(size(source.inside));
source.pcoh(inside_pos) = (abs(pcoh1) + abs(pcoh2))./2;
% imaginary part:
source.ipcoh = nan(size(source.inside));
source.ipcoh(inside_pos) = (imag(pcoh1) + imag(pcoh2))./2;

%% plot
cfg = [];
cfg.parameter = {'coh','pcoh', 'icoh', 'ipcoh'};
[interp] = ft_sourceinterpolate(cfg, source, mri_reslice);
% plot the coherence

cfg = [];
cfg.method        = 'ortho';
cfg.crosshair = 'yes';
cfg.funparameter  = 'coh';
cfg.figurename = 'coh';
cfg.funcolormap = 'jet';
ft_sourceplot(cfg, interp);

cfg = [];
cfg.method        = 'ortho';
cfg.crosshair = 'yes';
cfg.funparameter  = 'pcoh';
cfg.figurename = 'pcoh';
cfg.funcolormap = 'jet';

ft_sourceplot(cfg, interp);


cfg = [];
cfg.method        = 'ortho';
cfg.funparameter  = 'icoh';
cfg.funcolormap = 'jet';
cfg.location =  [9 59 29];
cfg.figurename = 'icoh';
ft_sourceplot(cfg, interp);


cfg = [];
cfg.method        = 'ortho';
cfg.funparameter  = 'ipcoh';
cfg.location =  [9 59 29];

cfg.figurename  = 'ipcoh';

cfg.funcolormap = 'jet';
ft_sourceplot(cfg, interp);
%% REPEAT WITH EEG
clear data* source*
load([deriv_pth_eeg filesep 'sub-'...
    eeg_code_ '_task-oddball_epochedStory.mat']);
data_eeg = epochStory;

load([deriv_pth_eeg filesep 'sub-'...
    eeg_code_ '_task-oddball__art.mat']);
art_eeg = cfg_art.artfctdef.visual.artifact;
art_eeg(art_eeg(:,2) < data_eeg.sampleinfo(1) |...
    art_eeg(:,1) > data_eeg.sampleinfo(2),:) = [];
%%
cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = [1 55];
cfg.channel = data_eeg.label(1:end-4);
cfg.detrend = 'yes';
data_eeg = ft_preprocessing(cfg, data_eeg);  

%% repeat the EEG source reconstruction with story data
load([hdm_pth_eeg filesep 'elec_fin.mat']);
load([hdm_pth_eeg filesep 'grid_indiv.mat']);
load([hdm_pth_eeg filesep 'vol_dipoli.mat']);
% vol = vol_improved;
load([hdm_pth_eeg filesep 'mri_reslice.mat']);

%% scale outwards a bit:
elec_fin.elecpos = elec_fin.elecpos.*1.05;
elec_fin.elecpos(:,2) = elec_fin.elecpos(:,2) + 8;
elec_fin.chanpos = elec_fin.elecpos;

%% alignment plot
figure,ft_plot_ortho(mri_reslice.anatomy, 'style','intersect','transform',mri_reslice.transform)
ft_plot_mesh(grid.pos(grid.inside,:),'vertexcolor',[1,0,0])
view([80,-15,0])
hold on;
ft_plot_mesh(vol.bnd(1), 'facecolor', ...
    'k', 'facealpha', 0.2); % brain
view([0 -1 0]); % from the right side
ft_plot_mesh(vol.bnd(2), 'facecolor', ...
    'g', 'facealpha', 0.2); % skull
view([0 -1 0]); % from the right side

ft_plot_mesh(vol.bnd(3), 'facecolor', ...
    'b', 'facealpha', 0.2); 
ft_plot_sens(elec_fin,'facecolor', 'g');

%%

data_nan = data_eeg;
samplevec1 = data_nan.sampleinfo(1,1):...
    data_nan.sampleinfo(1,2);
samplevec2 = data_nan.sampleinfo(2,1):...
    data_nan.sampleinfo(2,2);
sample_vec_all = [samplevec1, samplevec2];


artifact_all = zeros(size(sample_vec_all));
for art = 1: size(art_eeg,1)
    artifact_all(find(sample_vec_all==art_eeg(art,1)):...
        find(sample_vec_all==art_eeg(art,2))) = 1;
end


data_nan.trial{1}(:, find(...
    artifact_all(1:length(samplevec1)))) = nan;

data_nan.trial{2}(:, find(...
    artifact_all(length(samplevec1)+1:end))) = nan;

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
cfg.elec            = elec_fin;
cfg.keepleadfield   = 'yes';
cfg.normalize  = 'yes'; % normalize LF to attenuate depth bias
cfg.normalizeparam  = 0.5; % default = 0.5
cfg.lcmv.keepfilter = 'yes';
cfg.lcmv.fixedori   = 'yes';
cfg.lcmv.lambda     = '1%';  % must be '1%' not 1
source = ft_sourceanalysis(cfg, timelock);
inside_pos = find(source.inside);

% Extract all filters for inside voxels
W_eeg = cell2mat(source.avg.filter(source.inside));

%% localize the story variance increase (in delta directly)
cfg =[];
cfg.resamplefs = 250;
data_eeg_rs = ft_resampledata(cfg, data_eeg);
% source_trial1 = sourcemodel.trial{1};
% source_trial2 = sourcemodel.trial{2};
seg_length = 10*60*250 + 4*250; % 10 minutes of data?
source_trial1 = W_eeg*data_eeg_rs.trial{1};
source_trial2 = W_eeg*data_eeg_rs.trial{2};

%% make a source data structure and filter in the delta band
sourcedata.label = cell(numel(inside_pos), 1); % Initialize with empty cells

% Assign labels only to inside sources
for i = 1:numel(inside_pos)

    sourcedata.label{i} = sprintf('S_%d',...
        inside_pos(i));
end

sourcedata.trial{1} = source_trial1;
sourcedata.trial{2} = source_trial2;
sourcedata.time = data_eeg_rs.time;

cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = [2 5];
sourcedata = ft_preprocessing(cfg, sourcedata);

for tt = 1: numel(sourcedata.label)
    sourcedata.trial{1}(tt,:) = (hilbert(...
        sourcedata.trial{1}(tt,:)));
    sourcedata.trial{2}(tt,:) = (hilbert(...
        sourcedata.trial{2}(tt,:)));
end
%% try the variance increase again
sourceD_trial1 = abs(sourcedata.trial{1});
sourceD_trial2 = abs(sourcedata.trial{2});
sel = data_eeg_rs.time{1}>0;
% get the trial data from pre and post
datpre =...
    [sourceD_trial1(...
    :, data_eeg_rs.time{1}<0),sourceD_trial2(...
    :, data_eeg_rs.time{2}<0)];

datpost = [sourceD_trial1(:, sel(1:seg_length)),...
    sourceD_trial2(:, sel(1:seg_length))];

post_var = var(datpost, 0,2);
pre_var = var(datpre, 0,2);
% 

effect  = (post_var-pre_var)./(post_var+pre_var);
source.story = (nan(size(source.inside)));
source.story(source.inside) = effect;

%% plot
cfg = [];
cfg.parameter = {'story'};
[interp] = ft_sourceinterpolate(cfg, source, mri_reslice);
% plot the correlation

cfg = [];
cfg.method        = 'ortho';
cfg.crosshair = 'yes';
cfg.funparameter  = 'story';
cfg.funcolormap = 'jet';
% cfg.funcolorlim   = [0 0.06];
ft_sourceplot(cfg, interp);
%% try again with delta power phase coherence over time
% in first 10 minutes

source_phase1 = (sourcedata.trial{1}(:,1:seg_length))./...
    abs(sourcedata.trial{1}(:,1:seg_length));
source_phase2 = (sourcedata.trial{2}(:,1:seg_length))./...
    abs(sourcedata.trial{2}(:,1:seg_length));


delta_pc = abs(mean(source_phase1.*conj(source_phase2),2));


source.retest = (nan(size(source.inside)));
source.retest(source.inside) = delta_pc;

%% plot
close all;
cfg = [];
cfg.parameter = {'retest'};
[interp] = ft_sourceinterpolate(cfg, source, mri_reslice);
% plot the correlation

cfg = [];
cfg.method        = 'ortho';
cfg.crosshair = 'yes';
cfg.funparameter  = 'retest';
cfg.funcolormap = 'jet';
% cfg.funcolorlim   = [0 0.06];
% cfg.location = [2 38 48];
ft_sourceplot(cfg, interp);

%% pick peak source position in auditory cortex and compute simple correlation

[~, aud_index] =...
    max(source.retest(source.inside) );
% 
% [~, aud_index] =...
%     min(sum(abs(source.pos(inside_pos,:) - [9 59 29]),2));

source_tc1 = source_trial1(aud_index, :);
source_tc2 = source_trial2(aud_index, :);

%% compute correlation of auditory source with all others
source_corrs = zeros(numel(inside_pos),1);

for ii = 1: size(source_corrs, 1)
    ii
     source_corrs(ii) = (corr(source_trial1(ii,:)', ...
        source_tc1','Type', 'Spearman')+...
        corr(source_trial2(ii,:)', ...
        source_tc2','Type', 'Spearman'))/2;
end


source.corrs = nan(size(source.inside));
source.corrs(inside_pos) = source_corrs;
%% plot
cfg = [];
cfg.parameter = {'corrs'};
[interp] = ft_sourceinterpolate(cfg, source, mri_reslice);
% plot the correlation

cfg = [];
cfg.method        = 'ortho';
cfg.crosshair = 'yes';
cfg.funparameter  = 'corrs';
cfg.funcolormap = 'jet';
cfg.funcolorlim   = [0 0.8];
ft_sourceplot(cfg, interp);

%% coherence in delta band
% pick 10 minutes of story
source_trial1c = sourcedata.trial{1}(:,1:seg_length);
source_trial2c = sourcedata.trial{2}(:,1:seg_length);

% phase differences with the selected source
coh1 = source_trial1c * ...
    (source_trial1c(aud_index, :)'); % complex conjugate transpose!!!
coh2 = source_trial2c *...
    (source_trial2c(aud_index, :)');% complex conjugate transpose!!!

% for correct scaling we should divide by n
coh1 = coh1./length(source_trial1c);
coh2 = coh2./length(source_trial2c);

source.coh = nan(size(source.inside));
source.coh(inside_pos) = (abs(coh1) + abs(coh2))./2;
% imaginary part:
source.icoh = nan(size(source.inside));
source.icoh(inside_pos) = (imag(coh1) + imag(coh2))./2;

%% compute phase coherence as well:
% normalize to unit amplitude
clear source_phase*
source_trial1c = source_trial1c./abs(source_trial1c);
source_trial2c = source_trial2c./abs(source_trial2c);

% coherence is now automatically phase coherence
pcoh1 = (source_trial1c * ...
    (source_trial1c(aud_index, :)'));% complex conjugate transpose!!!
pcoh2 = (source_trial2c *...
    (source_trial2c(aud_index, :)'));% complex conjugate transpose!!!


% for correct scaling we should divide by n
pcoh1 = pcoh1./length(source_trial1c);
pcoh2 = pcoh2./length(source_trial2c);

source.pcoh = nan(size(source.inside));
source.pcoh(inside_pos) = (abs(pcoh1) + abs(pcoh2))./2;
% imaginary part:
source.ipcoh = nan(size(source.inside));
source.ipcoh(inside_pos) = (imag(pcoh1) + imag(pcoh2))./2;

%% plot
cfg = [];
cfg.parameter = {'coh','pcoh', 'icoh', 'ipcoh'};
[interp] = ft_sourceinterpolate(cfg, source, mri_reslice);
% plot the coherence

cfg = [];
cfg.method        = 'ortho';
cfg.crosshair = 'yes';
cfg.funparameter  = 'coh';
cfg.title = 'coh';
cfg.funcolormap = 'jet';
ft_sourceplot(cfg, interp);

cfg = [];
cfg.method        = 'ortho';
cfg.crosshair = 'yes';
cfg.funparameter  = 'pcoh';
cfg.title = 'pcoh';
cfg.funcolormap = 'jet';
cfg.title = 'pcoh';
ft_sourceplot(cfg, interp);

cfg = [];
cfg.method        = 'ortho';
cfg.crosshair = 'yes';
cfg.funparameter  = 'icoh';
cfg.title = 'icoh';
cfg.funcolormap = 'jet';
ft_sourceplot(cfg, interp);

cfg = [];
cfg.method        = 'ortho';
cfg.crosshair = 'yes';
cfg.funparameter  = 'ipcoh';
cfg.title = 'ipcoh';
cfg.funcolormap = 'jet';
cfg.title = 'pcoh';
ft_sourceplot(cfg, interp);