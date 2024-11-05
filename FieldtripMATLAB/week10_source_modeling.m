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

sjs = ...
    ['004']; % todo change to make BIDS compliant
%% set all paths

acq = '';
code_ = sjs(1,:); % pick one subject

anat_pth_ = ...
    [dataPath filesep 'sub-' code_ filesep 'anat'];
eeg_pth_ = ...
    [dataPath filesep 'sub-' code_ filesep 'eeg'];
deriv_pth = [dataPath filesep ...
    'derivatives' filesep 'preprocessing', filesep ...
    'sub-' code_ ];
if ~exist(deriv_pth,'dir'); mkdir(deriv_pth); end

hdm_pth = [dataPath filesep ...
    'derivatives' filesep 'headmodeling', filesep ...
    'sub-' code_ ];
%% prepare for plots..
cfg = [];
cfg.layout = 'BioSemi256';
lay = ft_prepare_layout(cfg);

cfg_topo= {'mask' lay.mask, 'datmask' [], 'interplim' 'mask'};
cfg_lay = { 'point' 'no' 'box' 'no' 'label' 'no' };
n_channels = 256;
%% load all information
load([deriv_pth filesep 'data_pruned.mat']);

load([hdm_pth filesep 'elec_fin_scaled.mat']);

load([hdm_pth filesep 'grid_indiv.mat']);
load([hdm_pth filesep 'vol_dipoli.mat']);
% vol = vol_improved;
load([hdm_pth filesep 'mri_reslice.mat']);
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
%% prepare only the leadfied first
cfg                  = [];
cfg.elec             = elec_fin;  % electrode positions
cfg.headmodel        = vol;   % volume conduction headmodel
cfg.sourcemodel      = grid;
lf                   = ft_prepare_leadfield(cfg);
%% example position
inside_pos = find(grid.inside);

figure,ft_plot_ortho(mri_reslice.anatomy, ...
    'style','intersect','transform',mri_reslice.transform)
hold on;
ft_plot_mesh(grid.pos(inside_pos(910),:), ...
    'vertexcolor', 'r');
%% plot the leadfield at position 910

lf_example = lf.leadfield(inside_pos(910));
figure
subplot(1,3,1), 
    ft_plot_topo(lay.pos(1:n_channels,1), lay.pos(1:n_channels,2), lf_example{1}(:,1), cfg_topo{:});
    ft_plot_layout(lay, cfg_lay{:})
    title('LF X direction')
subplot(1,3,2), 
    ft_plot_topo(lay.pos(1:n_channels,1), lay.pos(1:n_channels,2), lf_example{1}(:,2), cfg_topo{:});
    ft_plot_layout(lay, cfg_lay{:})
    title('LF Y direction')
subplot(1,3,3), 
    ft_plot_topo(lay.pos(1:n_channels,1), lay.pos(1:n_channels,2), lf_example{1}(:,3), cfg_topo{:});
    ft_plot_layout(lay, cfg_lay{:})
    title('LF Z direction');
%% now pull out the whole leadfield
lf_cell = lf.leadfield(inside_pos);
 
% count the number of channels and leadfield components
Nchan   = size(lf_cell{1},1);
Nsource = size(lf_cell{1},2)* numel(lf_cell);

% concatenate the leadfield components of all sources into one large matrix
H = zeros(Nchan, Nsource);
n = 1;
for i=1:size(inside_pos,1)
    cbeg = n;
    cend = n + size(lf_cell{i}, 2) - 1;
    H(:,cbeg:cend) = lf_cell{i};
    n = n + size(lf_cell{i}, 2);
end
%% plot the whole LF matrix
figure;imagesc(H); axis xy;
xlabel('source XYZ');
ylabel('channel');
title('leadfield matrix')

%% prepare the data (add a lowpass filter)
cfg = [];
cfg.lpfilter = 'yes';
cfg.lpfreq = 30;
cfg.detrend = 'yes';
data_fil = ft_preprocessing(cfg, data_pruned);  
%% compute the visual ERP
cfg = [];
cfg.trials = find(...
    cellfun(@(x) ...
    strcmp('semanticVis', string(x.blockType)), ...
    data_pruned.trialinfo));
cfg.covariance = 'yes';
cfg.covariancewindow = [-0.5 0]; %it will calculate the covariance matrix
% cfg.demean              = 'yes';    % we demean (baseline correct) ...
% cfg.baselinewindow      = [-0.5 0]; % using the mean activity in this window

tl_vis = ft_timelockanalysis(cfg, data_fil);  
dat = tl_vis.avg;
%% pull out covariance matrix for later
C = tl_vis.cov;
% regularize the covariance matrix
nu = trace(C)/n_channels;
I  = nu*eye(n_channels);
gamma = 0.0001;
% The regularized covariance matrix is a convex combination of the empirical
% covariance matrix and a diagonal covariance matrix with equal total variance.
C = C*(1-gamma) + I*gamma;
%% MNE manual
t_axis = tl_vis.time;

Wmne = pinv(H);

%% MNE manual (Lin paper)

H2 = inv(C)*H;
dat2 = inv(C)*dat;
% estimate lambda? 
lambda = 0.05;

Wmne2 = H2'*inv(H2*H2'+ (lambda^2)*eye(size(H2,1)));

% Wmne = Wmne2;
%%
%% plot leadfield and filter
indx = 3*910;
figure
subplot(2,3,1), 
    ft_plot_topo(lay.pos(1:n_channels,1), lay.pos(1:n_channels,2), H(:,indx), cfg_topo{:});
    ft_plot_layout(lay, cfg_lay{:})
    title('LF X direction')
subplot(2,3,2), 
    ft_plot_topo(lay.pos(1:n_channels,1), lay.pos(1:n_channels,2), H(:,indx+1), cfg_topo{:});
    ft_plot_layout(lay, cfg_lay{:})
    title('LF Y direction')
subplot(2,3,3), 
    ft_plot_topo(lay.pos(1:n_channels,1), lay.pos(1:n_channels,2), H(:,indx+2), cfg_topo{:});
    ft_plot_layout(lay, cfg_lay{:})
    title('LF Z direction');
subplot(2,3,4), 
    ft_plot_topo(lay.pos(1:n_channels,1), lay.pos(1:n_channels,2), Wmne(indx,:), cfg_topo{:});
    ft_plot_layout(lay, cfg_lay{:})
    title('LF X direction')
subplot(2,3,5), 
    ft_plot_topo(lay.pos(1:n_channels,1), lay.pos(1:n_channels,2), Wmne(indx+1,:), cfg_topo{:});
    ft_plot_layout(lay, cfg_lay{:})
    title('LF Y direction')
subplot(2,3,6), 
    ft_plot_topo(lay.pos(1:n_channels,1), lay.pos(1:n_channels,2), Wmne(indx+2,:), cfg_topo{:});
    ft_plot_layout(lay, cfg_lay{:})
    title('LF Z direction');

%% create a template data structure for plotting

source_template = ft_datatype_source(grid);
source_template.erp = ...
    nan(size(source_template.inside));

%% project the ERP into source space
source_dat = (Wmne * dat);

source1D = nan(size(source_dat, 1)/3, ...
    size(source_dat, 2));
%estimate the source strength across directions
n = 1;
for i=1:size(inside_pos,1)
    cbeg = n;
    cend = n +  2;
    source1D(i,:) = sum(source_dat(cbeg:cend,:).^2,1);
    n = n + 3;
end

%% baseline correct the ERP
source_template.erp(source_template.inside) = ...
    (mean(source1D(:,nearest(t_axis, 0.12): ...
    nearest(t_axis, 0.15)),2) - ...
    mean(source1D(:,nearest(t_axis, -0.5): ...
    nearest(t_axis, 0)),2))./...
    mean(source1D(:,nearest(t_axis, -0.5): ...
    nearest(t_axis, 0)),2);


cfg = [];
cfg.parameter = {'erp'};
[interp] = ft_sourceinterpolate(cfg, ...
    source_template, mri_reslice);

cfg = [];
cfg.method        = 'ortho';
cfg.crosshair = 'no';
cfg.funparameter  = 'erp';
cfg.funcolormap = 'jet';
% cfg.funcolorlim   = [-12 12];
cfg.location = [0 -59 15];
ft_sourceplot(cfg, interp);

%% MNE fieldtrip
cfg = [];
cfg.trials = find(...
    cellfun(@(x) ...
    strcmp('semanticVis', string(x.blockType)), ...
    data_pruned.trialinfo));
cfg.covariance = 'yes';
cfg.covariancewindow = [-0.7 0]; %it will calculate the covariance matrix
% on the timepoints that are before the zero-time point in the trials
tl_vis = ft_timelockanalysis(cfg, data_fil);
%%
cfg = [];
cfg.trials = find(...
    cellfun(@(x) ...
    strcmp('semanticVis', string(x.blockType)), ...
    data_pruned.trialinfo));
cfg.covariance = 'yes';
cfg.covariancewindow = [-0.5 0]; %it will calculate the covariance matrix
% cfg.demean              = 'yes';    % we demean (baseline correct) ...
% cfg.baselinewindow      = [-0.5 0]; % using the mean activity in this window

tl_vis = ft_timelockanalysis(cfg, data_fil);  

%% mne source reconstruction
cfg                 = [];
cfg.trials = find(...
    cellfun(@(x) ...
    strcmp('semanticVis', string(x.blockType)), ...
    data_pruned.trialinfo));
cfg.method          = 'mne';
cfg.sourcemodel     = grid;
cfg.headmodel       = vol;
cfg.elec            = elec_fin;
cfg.mne.prewhiten = 'yes';
cfg.mne.lambda    = 0.05;
cfg.mne.scalesourcecov = 'yes';
sourceMNEvis = ft_sourceanalysis(cfg, tl_vis);
%% plot
t_axis = sourceMNEvis.time;
sourceMNEvis.diff = ...
    (mean(sourceMNEvis.avg.pow(:,nearest(t_axis, 0.12): ...
    nearest(t_axis, 0.15)),2) - ...
    mean(sourceMNEvis.avg.pow(:,nearest(t_axis, -0.5):...
    nearest(t_axis, 0)),2))./...
     mean(sourceMNEvis.avg.pow(:,nearest(t_axis, -0.5):...
    nearest(t_axis, 0)),2);


cfg = [];
cfg.parameter = {'diff'};
[interp] = ft_sourceinterpolate(cfg, ...
    sourceMNEvis, mri_reslice);


cfg = [];
cfg.method        = 'ortho';
cfg.crosshair = 'no';
cfg.funparameter  = 'diff';
cfg.funcolormap = 'jet';
% cfg.funcolorlim   = [-7.5 7.5];
cfg.location = [0 -59 15];
ft_sourceplot(cfg, interp);

%% LCMV manual
%% compute data covariance
cfg = [];
% cfg.trials = find(...
%     cellfun(@(x) ...
%     strcmp('semanticVis', string(x.blockType)), ...
%     data_pruned.trialinfo));
cfg.covariance = 'yes';
cfg.covariancewindow = 'all'; %it will calculate the covariance matrix
% on the timepoints that are before the zero-time point in the trials
tl_all = ft_timelockanalysis(cfg, data_fil);
%% regularize !
C = tl_all.cov;

nu = trace(C)/n_channels;
I  = nu*eye(n_channels);
gamma = 0.0001;
% The regularized covariance matrix is a convex combination of the empirical
% covariance matrix and a diagonal covariance matrix with equal total variance.
C = C*(1-gamma) + I*gamma;
%%
% Van Veen formula 23 (use pinv for badly scaled matrix
Wlcmv = pinv(H' * inv(C) * H) * H' * inv(C);
%% project on the main axis
Wlcmv2 = nan(size(Wlcmv,1)/3, size(Wlcmv,2));
n = 1;
for i=1:size(inside_pos,1)
    
    cfilt = ...
        Wlcmv(cbeg:cend,:) * C * Wlcmv(cbeg:cend,:)';

    cbeg = n;
    cend = n +  2;

    [u, s, v] = svd(cfilt);
    Wlcmv2(i,:)  = u(:,1)' *  Wlcmv(cbeg:cend,:);
    n = n + 3;
end

%% compute the increase in variance and the contrast:
cfg = []; 
cfg.trials = find(...
    cellfun(@(x) ...
    strcmp('semanticVis', string(x.blockType)), ...
    data_pruned.trialinfo));
cfg.latency            = [-0.5 0]; 

data_vis_pre = ft_selectdata(cfg, data_fil);
cfg.latency            = [0 1]; 
data_vis_post = ft_selectdata(cfg, data_fil);

% get the trial data from pre and post
datpre = cell2mat(data_vis_pre.trial);
datpost = cell2mat(data_vis_post.trial);

post_var = var(Wlcmv2*datpost, 0,2);
pre_var = var(Wlcmv2*datpre, 0,2);
% 

effect  = (post_var-pre_var)./(post_var+pre_var);
source_template.vis = (nan(size(source_template.inside)));
source_template.vis(source_template.inside) = effect;

cfg = [];
cfg.parameter = {'vis'};
[interp] = ft_sourceinterpolate(cfg,...
    source_template, mri_reslice);


cfg = [];
cfg.method        = 'ortho';
cfg.crosshair = 'no';
cfg.funparameter  = 'vis';
cfg.funcolormap = 'jet';
cfg.funcolorlim   = [-0 0.35];
cfg.location = [0 -59 15];
ft_sourceplot(cfg, interp);

%% LCMV computes leadfield on the fly (use one common filter)
cfg                 = [];
cfg.method          = 'lcmv';
cfg.sourcemodel     = grid;
cfg.headmodel       = vol;
cfg.elec            = elec_fin;
cfg.keepleadfield   = 'yes';
cfg.lcmv.keepfilter = 'yes';
cfg.lcmv.fixedori   = 'yes';
cfg.lcmv.lambda     = '1%';  % must be '1%' not 1
source = ft_sourceanalysis(cfg, data_fil);

% Extract all filters for inside voxels
W = cell2mat(source.avg.filter(source.inside));

%% make nice plots of leadfield for an example position

figure,ft_plot_ortho(mri_reslice.anatomy, ...
    'style','intersect','transform',mri_reslice.transform)
hold on;
ft_plot_mesh(grid.pos(inside_pos(910),:), ...
    'vertexcolor', 'r');

inside_pos = find(grid.inside);
grid.pos(inside_pos(910),:)

lf_example = source.leadfield(inside_pos(910));



%% plot the leadfield and filter for position 910
n_channels = 256;
figure
subplot(1,4,1), 
    ft_plot_topo(lay.pos(1:n_channels,1), lay.pos(1:n_channels,2), lf_example{1}(:,1), cfg_topo{:});
    ft_plot_layout(lay, cfg_lay{:})
    title('LF X direction')
subplot(1,4,2), 
    ft_plot_topo(lay.pos(1:n_channels,1), lay.pos(1:n_channels,2), lf_example{1}(:,2), cfg_topo{:});
    ft_plot_layout(lay, cfg_lay{:})
    title('LF Y direction')
subplot(1,4,3), 
    ft_plot_topo(lay.pos(1:n_channels,1), lay.pos(1:n_channels,2), lf_example{1}(:,3), cfg_topo{:});
    ft_plot_layout(lay, cfg_lay{:})
    title('LF Z direction')

subplot(1,4,4), 
    ft_plot_topo(lay.pos(1:n_channels,1), lay.pos(1:n_channels,2), W(910,:), cfg_topo{:});
    ft_plot_layout(lay, cfg_lay{:})
    title('Filter')

%% compute the increase in variance and the contrast:
cfg = []; 
cfg.trials = find(...
    cellfun(@(x) ...
    strcmp('semanticVis', string(x.blockType)), ...
    data_pruned.trialinfo));
cfg.latency            = [-1 0]; 

data_vis_pre = ft_selectdata(cfg, data_fil);
cfg.latency            = [0 2]; 
data_vis_post = ft_selectdata(cfg, data_fil);

% get the trial data from pre and post
datpre = cell2mat(data_vis_pre.trial);
datpost = cell2mat(data_vis_post.trial);

post_var = var(W*datpost, 0,2);
pre_var = var(W*datpre, 0,2);
% 

effect  = (post_var-pre_var)./(post_var+pre_var);
source.vis = (nan(size(source.inside)));
source.vis(source.inside) = effect;


%%
cfg.trials = find(...
    cellfun(@(x) ...
    strcmp('semanticAud', string(x.blockType)), ...
    data_fil.trialinfo)|...
    ...
    cellfun(@(x) ...
    strcmp('classicalAud', string(x.blockType)), ...
    data_fil.trialinfo));
cfg.latency            = [-1 0]; 

data_aud_pre = ft_selectdata(cfg, data_fil);
cfg.latency            = [0 2]; 
data_aud_post = ft_selectdata(cfg, data_fil);

% get the trial data from pre and post
datpre = cell2mat(data_aud_pre.trial);
datpost = cell2mat(data_aud_post.trial);

post_var = var(W*datpost, 0,2);
pre_var = var(W*datpre, 0,2);
% 

effect  = (post_var-pre_var)./(post_var+pre_var);
source.aud = (nan(size(source.inside)));
source.aud(source.inside) = effect;
%%
source.vis_min_aud = source.vis - source.aud;

% %% with shared grid we can overwrite the grid interpolate onto mni mri
% ftpath = fileparts(which('ft_defaults'));
% ft_read_mri([ftpath filesep, 'template', filesep ...
%     'anatomy',filesep,  'single_subj_T1_1mm.nii' ])
%%
cfg = [];
cfg.parameter = {'aud', 'vis', 'vis_min_aud'};
[interp] = ft_sourceinterpolate(cfg, source, mri_reslice);


cfg = [];
cfg.method        = 'ortho';
cfg.crosshair = 'no';
cfg.funparameter  = 'vis';
cfg.funcolormap = 'jet';
cfg.funcolorlim   = [0 0.06];
cfg.location = [0 -59 15];
ft_sourceplot(cfg, interp);

%% plot the audio

cfg = [];
cfg.method        = 'ortho';
cfg.crosshair = 'no';
% cfg.funcolorlim   = [-0.2 0.2];
cfg.funcolormap = 'jet';
cfg.location = [0 0 0];
cfg.funparameter  = 'aud';
ft_sourceplot(cfg, interp);



%% split in l/r for better separation of correlated sources
right_selection = find(elec_fin.elecpos(:,1)>0);
left_selection = find(elec_fin.elecpos(:,1)<=0);

elec_right = elec_fin;
elec_right.elecpos = ...
    elec_fin.elecpos(right_selection,:);
elec_right.chanpos = elec_right.elecpos;
elec_right.chantype = elec_fin.chantype(right_selection);
elec_right.chanunit = elec_fin.chanunit(right_selection);
elec_right.label = elec_fin.label(right_selection);

elec_left = elec_fin;
elec_left.elecpos = ...
    elec_fin.elecpos(left_selection,:);
elec_left.chanpos = elec_left.elecpos;
elec_left.chantype = elec_fin.chantype(left_selection);
elec_left.chanunit = elec_fin.chanunit(left_selection);
elec_left.label = elec_fin.label(left_selection);

cfg = [];
cfg.channel = data_fil.label(right_selection);
data_right = ft_selectdata(cfg, data_fil);
cfg.channel = data_fil.label(left_selection);
data_left = ft_selectdata(cfg, data_fil);

cfg                 = [];
cfg.method          = 'lcmv';
cfg.sourcemodel     = grid;
cfg.headmodel       = vol;
cfg.elec            = elec_right;
cfg.keepleadfield   = 'yes';
cfg.lcmv.keepfilter = 'yes';
cfg.lcmv.fixedori   = 'yes';
cfg.lcmv.lambda     = '1%';  % must be '1%' not 1
source_right = ft_sourceanalysis(cfg, data_right);
cfg.elec            = elec_left;
source_left = ft_sourceanalysis(cfg, data_left);

all_inside_pos = grid.pos(grid.inside,:);
right_inside = find(all_inside_pos(:,1)>0);
left_inside = find(all_inside_pos(:,1)<=0);

Wr = cell2mat(source_right.avg.filter(source_right.inside));
Wl = cell2mat(source_left.avg.filter(source_left.inside));

Wr = Wr(right_inside,:);
Wl = Wl(left_inside,:);

%%
cfg = [];
cfg.trials = find(...
    cellfun(@(x) ...
    strcmp('semanticAud', string(x.blockType)), ...
    data_fil.trialinfo)|...
    ...
    cellfun(@(x) ...
    strcmp('classicalAud', string(x.blockType)), ...
    data_fil.trialinfo));
cfg.latency            = [-1 0]; 

% right side
data_aud_pre = ft_selectdata(cfg, data_right);
cfg.latency            = [0 2]; 
data_aud_post = ft_selectdata(cfg, data_right);

% get the trial data from pre and post
datpre = cell2mat(data_aud_pre.trial);
datpost = cell2mat(data_aud_post.trial);

post_var_r = var(Wr*datpost, 0,2);
pre_var_r = var(Wr*datpre, 0,2);

% left side
cfg.latency            = [-1 0]; 
data_aud_pre = ft_selectdata(cfg, data_left);
cfg.latency            = [0 2]; 
data_aud_post = ft_selectdata(cfg, data_left);

% get the trial data from pre and post
datpre = cell2mat(data_aud_pre.trial);
datpost = cell2mat(data_aud_post.trial);

post_var_l = var(Wl*datpost, 0,2);
pre_var_l = var(Wl*datpre, 0,2);

% 
% post_var = (post_var_r + post_var_l) ./2;
% pre_var = (pre_var_r + pre_var_l) ./2;
% 
% effect  = (post_var-pre_var)./pre_var;

source.aud = (nan(size(grid.inside)));
inside_idx = find(grid.inside);
source.aud(inside_idx(left_inside)) ...
    = (post_var_l-pre_var_l)./(post_var_l+pre_var_l);
source.aud(inside_idx(right_inside))...
    = (post_var_r-pre_var_r)./(post_var_r+pre_var_r);

cfg = [];
cfg.parameter = 'aud';
[interp] = ft_sourceinterpolate(cfg, source, mri_reslice);


cfg = [];
cfg.method        = 'ortho';
cfg.funparameter  = 'aud';
cfg.funcolormap = 'jet';
cfg.location = [40 8 -6];
% cfg.funcolorlim   = [-0.1 0.1];
ft_sourceplot(cfg, interp);