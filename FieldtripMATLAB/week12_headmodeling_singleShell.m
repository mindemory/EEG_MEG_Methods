
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
    ['002']; % todo change to make BIDS compliant
%% set all paths

acq = '';
code_ = sjs(1,:); % pick one subject

anat_pth_ = ...
    [dataPath filesep 'sub-' code_ filesep 'anat'];
meg_pth_ = ...
    [dataPath filesep 'sub-' code_ filesep 'meg'];
deriv_pth = [dataPath filesep ...
    'derivatives' filesep 'meg_preprocessing', filesep ...
    'sub-' code_ ];
if ~exist(deriv_pth,'dir'); mkdir(deriv_pth); end

hdm_pth = [dataPath filesep ...
    'derivatives' filesep 'meg_headmodeling', filesep ...
    'sub-' code_ ];
%% layout with YOKOGAWA is not defined in FT
% we can hack one together:

hdr = ft_read_header(fullfile(meg_pth_, ...
    ['sub-' code_ '_task-oddball_eventmodel.sqd']));
% read the position of the sensors from the data
grad                        = ft_read_sens(fullfile(meg_pth_, ...
    ['sub-' code_ '_task-oddball_eventmodel.sqd'])); % this can be inspected with ft_plot_sens(grad)

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

% plot the custom layout
figure;
ft_plot_layout(lay, 'box', 1);
%% read in the headhape too (it's useless for this subject because fiducials are not aligned)
% Specify variable properties
tbl = readtable(fullfile(anat_pth_, ['sub-' code_ '_headshape.tsv']), 'FileType', 'text');
tbl.Properties.VariableNames = {'Label', 'X', 'Y', 'Z'};
%  headshape information + fiducials
% Step 2: Extract the relevant data
tmp_pos = [tbl.X, tbl.Y, tbl.Z] ;
%% something is off here... ?? fixed now..
ft_plot_mesh(tmp_pos(1:8,:), 'vertexcolor', 'r', 'vertexsize', 20);
hold on;
ft_plot_mesh(tmp_pos(9:end,:), 'vertexcolor', 'k',  'vertexsize', 5)

% pos_new = tmp_pos(4:end,:);
% pos_new(1:3, :) = (tmp_pos(1:3,:) + tmp_pos(4:6,:))./2;
headshape.pos = tmp_pos; % 
% Step 3: Create the FieldTrip headshape structure
headshape.unit = 'mm'; % or 'cm', depending on your data
headshape.label = tbl.Label(1:end);

%% first align the headshape to the MEG fiducials!
elec_dummy = headshape;
elec_dummy.elecpos = headshape.pos;
elec_dummy.chanpos = elec_dummy.elecpos;
%% get the HPI coil positions:

hpi_shape1 = ft_read_headshape(fullfile(meg_pth_, ...
    '240926-1.mrk'));
hpi_shape2 = ft_read_headshape(fullfile(meg_pth_, ...
    '240926-2.mrk'));
% average the 2 recordis
fids = 0.5.*(hpi_shape1.fid.pos(:,:) + ...
    hpi_shape2.fid.pos(:,:)).*10; % convert to mm

%% make an unaligned plot:
ft_plot_mesh(tmp_pos(1:8,:), 'vertexcolor', 'r', 'vertexsize', 20);
hold on;
ft_plot_mesh(tmp_pos(9:end,:), 'vertexcolor', 'k',  'vertexsize', 5);
ft_plot_mesh(fids, 'vertexcolor', 'g',  'vertexsize', 20);
ft_plot_mesh(tmp_pos(8,:), 'vertexcolor', 'y',  'vertexsize', 5)

% meg-fids go: middle (nas), lpa, rpa, l_extra, r_extra)
% head-fids go: 1_nas (bottom), 2_lpa_tragus, 3_rpa_tragus, 
% 4_lpa_coil, 5_rpa_coil, 6_nas (top/coil), 7_l_extra, 8_r_extra

ft_plot_sens(grad_mm);
set(gcf, 'color','w');
%% now we can align based on the fiducial positions

elec_coil.elecpos = fids; % convert to mm
elec_coil.label = {'nas', 'lpa', 'rpa', 'l_extra', 'r_extra'};
elec_coil.unit  = 'mm';
elec_dummy.label{4} = 'lpa';
elec_dummy.label{5} = 'rpa';
elec_dummy.label{6} = 'nas';
elec_dummy.label{7} = 'l_extra';
elec_dummy.label{8} = 'r_extra';

% coregister the electrodes to the MRI using fiducials
cfg = [];
cfg.method   = 'fiducial';
cfg.template = elec_coil;
cfg.elec     = elec_dummy;
cfg.feedback = 'yes';
cfg.fiducial = {'nas', 'lpa', 'rpa'};%, 'l_extra', 'r_extra'};
%red before /magenta after /blue mri reference
elec_aligned = ft_electroderealign(cfg);
%% overwrite headshape.pos and plot aligned grad + hs
headshape.pos = elec_aligned.chanpos;
%% make an aligned plot:
ft_plot_mesh(headshape.pos(1:8,:), 'vertexcolor', 'r', 'vertexsize', 5);
hold on;
ft_plot_headshape(headshape,'vertexcolor', 'k');
% ft_plot_sens(grad_mm);
hold on;

ft_plot_mesh(fids, 'vertexcolor', 'g',  'vertexsize', 30);
hold on;
ft_plot_sens(grad_mm, 'style', '*b', ...
    'facecolor' , 'y', 'facealpha' , 0.5);
%% get the MRI
mriFile = ...
    [anat_pth_ filesep 'sub-' code_ '_T1w.nii'];


% subject specific folder
save_pth = fullfile(deriv_pth);
if ~ exist(save_pth, 'dir'); mkdir(save_pth); end

% read mri 
mri = ft_read_mri(mriFile);
fprintf('...reading mri \n');

coordsys = ft_determine_coordsys(mri);
mri.coordsys = coordsys.coordsys;
mri_orig = mri;
%% now align the MRI to the headshape!
%% start by putting in this transform
% 0 15 -90
% 1 1 1
% 12 0 30
cfg = [];
cfg.method = 'headshape';  % Use the 'fiducial' method
cfg.headshape.headshape = headshape;
cfg.headshape.icp = 'yes'; 
mri_aligned = ft_volumerealign(cfg, mri);
%% double check the alignment after icp:
ft_volumerealign(cfg, mri_aligned);
%% now save the aligned mri

ft_plot_sens(grad_mm, 'style', '*b', ...
    'facecolor' , 'y', 'facealpha' , 0.5);
hold on;
ft_plot_mesh(elec_coil.elecpos, 'vertexcolor', 'r');
ft_plot_ortho(mri_aligned.anatomy,'transform',mri_aligned.transform,'style','intersect')
%% needs to be resliced so VOL is aligned to grad
mri_reslice = ft_volumereslice([],mri_aligned);

%% now segment the mri (we only need the brain)

cfg = [];
cfg.output = {'brain'}; % but singleshell only uses brain anyways
segmentedmri = ft_volumesegment(cfg, mri_reslice);
% OR map is brain
% segmentedmri.brain = segmentedmri.gray | segmentedmri.white;
%% plot again in 2 figures to check

figure;
ft_plot_ortho(segmentedmri.brain, 'style', 'intersect','transform',segmentedmri.transform);
% %% first manually rotate and translate to help the icp algorithm. After that ICP is always run

ft_plot_sens(grad_mm, 'style', '*b', ...
    'facecolor' , 'y', 'facealpha' , 0.5);


%% now overwrite the transform of the segmentation and prepare the vol
% segmentedmri.transform = mri.transform;
cfg = [];
cfg.method='singleshell';
vol = ft_prepare_headmodel(cfg, segmentedmri);
close all
 %note: the vol is now in headspace and doesn't have a transform anymore.
ft_plot_mesh(vol.bnd)
save vol vol;
save segmentedmri segmentedmri;

%%
load('mri_reslice.mat');
load('vol.mat');
figure;
ft_plot_mesh(vol.bnd(1), 'facecolor', 'r', 'facealpha', 0.5);
hold on;
ft_plot_ortho(mri_reslice.anatomy,'transform',mri_reslice.transform,'style','intersect')

hold on;
ft_plot_sens(grad_mm, 'style', '*b', 'facecolor' , 'y', 'facealpha' , 0.5);
view(25, 10)
set(gcf, 'color', 'w')
saveas(gcf, 'position_in_dewar.png');
%% this is for later:
%% plot everything one more time
figure

ft_plot_ortho(mri_reslice.anatomy, 'style', 'intersect', ...
    'transform', mri_reslice.transform)
hold on;
ft_plot_mesh(vol.bnd, 'facecolor', 'r','facealpha', 0.2); % scalp
view([0 -1 0]); % from the right side

hold on;
ft_plot_sens(grad_mm, 'style', '*b', 'facecolor' , 'y', 'facealpha' , 0.5);
%% create boundary element method forward model

% Create a volume conduction model
cfg        = [];
cfg.method = 'dipoli'; % You can also specify 'openmeeg', 'bemcp', or another method
vol  = ft_prepare_headmodel(cfg, vol);

save vol_dipoli vol -v7.3

%% also prepare the grid
% %make a subject specific grid from the inverse warped template
% %% use the warped template grid
% % we use a grid defined in MNI space. It is inverse-warped onto the subject
% % MRI. We can use a predefined standard sourcemodel: http://www.fieldtriptoolbox.org/template/sourcemodel/
% ftpath = fileparts(which('ft_defaults'));
% load(fullfile(ftpath, 'template/sourcemodel/standard_sourcemodel3d10mm'))
% template_grid = ft_convert_units(sourcemodel,'mm');
% clear sourcemodel;
% 
% % use the resliced MRI and warp the grid to it. 
% cfg = [];
% cfg.warpmni             = 'yes';
% cfg.template            = template_grid;
% cfg.nonlinear           = 'yes'; % use non-linear normalization
% cfg.mri                 = mri_reslice;
% cfg.sourcemodel.unit    = 'mm';
% grid = ft_prepare_sourcemodel(cfg);   % ras + n


%%
% save grid grid


%% or prepare an individual grid for this person only
% Prepare the source model grid
cfg = [];
cfg.grid.resolution = 10; % in mm
cfg.grid.unit = 'mm';
cfg.tight = 'yes'; % ensures grid points are inside the brain
cfg.inwardshift = 2; % shift grid points slightly inward to ensure they are inside the brain
cfg.headmodel = vol;
grid = ft_prepare_sourcemodel(cfg);
%% plot individual grid
ft_plot_mesh(vol.bnd, 'facecolor', 'k',...
    'facealpha', 0.2); % scalp
hold on;
ft_plot_mesh(grid.pos(grid.inside,:), ...
    'vertexcolor', 'r');
%% or prepare a template grid with inward shift and warp it to this persons space
ftpath = fileparts(which('ft_defaults'));
template = load(fullfile(ftpath, 'template/headmodel/standard_singleshell'));
template_vol = ft_convert_units(template.vol, 'mm');
% Prepare the source model grid
cfg = [];
cfg.grid.resolution = 10; % in mm
cfg.grid.unit = 'mm';
cfg.tight = 'yes'; % ensures grid points are inside the brain
cfg.inwardshift = 5; % shift grid points 5mm inward
cfg.headmodel = template_vol;
template_grid = ft_prepare_sourcemodel(cfg);
%% plot template grid
ft_plot_mesh(template_vol.bnd, 'facecolor', 'k',...
    'facealpha', 0.2); % brain
hold on;
ft_plot_mesh(template_grid.pos(template_grid.inside,:), ...
    'vertexcolor', 'r');

%% now warp the new template grid to individual head space:
% use the resliced MRI and warp the grid to it. 
cfg = [];
cfg.warpmni             = 'yes';
cfg.template            = template_grid;
cfg.nonlinear           = 'yes'; % use non-linear normalization
cfg.mri                 = mri_reslice;
cfg.sourcemodel.unit    = 'mm';
grid = ft_prepare_sourcemodel(cfg);   % ras + 
%% plot the shared individual grid
ft_plot_mesh(vol.bnd, 'facecolor', 'g',...
    'facealpha', 0.2); % brain
hold on;
ft_plot_mesh(grid.pos(grid.inside,:), ...
    'vertexcolor', 'r');
%%

ft_plot_ortho(mri_reslice.anatomy, 'style', 'intersect', ...
    'transform', mri_reslice.transform)
hold on;
ft_plot_mesh(vol.bnd, 'facecolor', 'k', 'facealpha', 0.2); % brain
view([0 -1 0]); % from the right side


ft_plot_mesh(grid.pos(grid.inside,:), ...
    'vertexcolor', [1 0 0]);

%% SLIDES

%% %END END

