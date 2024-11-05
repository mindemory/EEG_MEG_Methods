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
%%

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


cfg = [];
cfg.dataset = fullfile(eeg_pth_, ...
    ['sub-' code_ '_task-oddball_eeg.bdf']);
cfg.continuous = 'yes';
% data_raw = ft_preprocessing(cfg);
hdr = ft_read_header(cfg.dataset);

%% LOAD IN ELEC AND HEADSHAPE
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
tbl = readtable(fullfile(eeg_pth_, ['sub-' code_ '_electrodes.tsv']), opts);

%% Keep the electrode positions
n_elec = 256;
elec = [];
elec.label = hdr.label(1:n_elec);
elec.elecpos = ...
    [tbl.x1(1:n_elec), ...
    tbl.y1(1:n_elec), ...
    tbl.z1(1:n_elec)];
elec.chanpos = elec.elecpos;
elec.unit = 'mm';
n_elec = size(elec.elecpos,1);
%% read in the headhape too
opts = delimitedTextImportOptions("NumVariables", 6);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["X", "y", "z", "label", "VarName5", "VarName6"];
opts.VariableTypes = ["double", "double", "double", "categorical", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["VarName5", "VarName6"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["label", "VarName5", "VarName6"], "EmptyFieldRule", "auto");

% Specify variable properties
tbl = readtable(fullfile(anat_pth_, ['sub-' code_ '_headshape.tsv']), opts);
%%  headshape information + fiducials
% Step 2: Extract the relevant data
headshape.pos = [tbl.X, tbl.y, tbl.z];

% Step 3: Create the FieldTrip headshape structure
headshape.unit = 'mm'; % or 'cm', depending on your data
headshape.label = cell(size(headshape.pos,1),1); % add labels if available
n_hs = numel(headshape.label);
for ll = 1: numel(headshape.label)
    if ll < 4; headshape.label{ll}; % skip... = string(tbl.label(ll));
    else
        headshape.label{ll} = 'extra';
    end
end
headshape.label{1} = 'na';
headshape.label{2} = 'lpa';
headshape.label{3} = 'rpa';
%% Visualize the headshape and electrodes
ft_plot_headshape(headshape);
hold on;
ft_plot_sens(elec);
%% read in the anatomical information
% individual MRI scan
mri_file = ...
    [anat_pth_ filesep 'sub-' code_ '_T1w.nii.gz'];
% read MRI image
mri_orig = ft_read_mri(mri_file);

%% determine the coordinate system
coordsys = ft_determine_coordsys(mri_orig);
% and save
fprintf(['coordsys defined as: ' ...
    coordsys.coordsys '\n']);
mri_orig.coordsys = coordsys.coordsys;
%% do a first segmentation keeping only the scalp for 
% alignment and plotting
cfg =[];
cfg.output = 'scalp';
segmentedmri = ft_volumesegment(cfg, mri_orig);
%% initial plot of segmentation
% add anatomical information to the segmentation
segmentedmri.transform = mri_orig.transform;
segmentedmri.anatomy   = mri_orig.anatomy;

% typically axes are A L S and N in this dataset
% plot the segmented mri 
cfg              = [];
cfg.funparameter = 'scalp';
ft_sourceplot(cfg, segmentedmri);
%%  make a 3D model of the scalp for plotting 
% and alignment
cfg        = [];
cfg.tissue = 'scalp';
cfg.method = 'singleshell';
hdm        = ft_prepare_headmodel(cfg, segmentedmri);

%% combine the elec-structure with fiducials and 
% headshape positions for optimal alignment
elec_unaligned = elec;
elec_unaligned.elecpos = cat(1, headshape.pos(1:3,:),...
    elec.elecpos(1:n_elec,:), headshape.pos(4:end,:));
elec_unaligned.label = cat(1, headshape.label(1:3), ...
    elec.label(1:n_elec), headshape.label(4:end));
elec_unaligned.chanpos = elec_unaligned.elecpos;
elec_unaligned

%% make an unaligned plot:
figure;
ft_plot_mesh(hdm.bnd,'facecolor','none'); 

hold on;
ft_plot_sens(elec_unaligned, 'facecolor', ...
    [0 1 0], 'edgecolor', [1 0 0]);

%% first align the electrodes to the MRI

% best to read these out from MRICroGL
nas_vox  =  [96, 219, 135];%
lpa_vox =  [12 122 93];% 
rpa_vox =  [178 120 96];%
transm = mri_orig.transform;
% bring the fiducials from voxel-space to headspace
nas = ft_warp_apply(transm, nas_vox, 'homogenous');
lpa = ft_warp_apply(transm, lpa_vox, 'homogenous');
rpa = ft_warp_apply(transm, rpa_vox, 'homogenous');

%-- alternatively read out with fieldtrip function
% cfg = [];
% cfg.method = 'interactive';
% %cfg.method = 'headshape';
% %cfg.headshape = headshape;
% % NOTE: this is just to interactively pick the positions in voxel space
% mri_trash = ft_volumerealign(cfg, mri_orig); % don't use the output
% % save at this point?
% nas_vox = mri_trash.cfg.fiducial.nas;
% lpa_vox = mri_trash.cfg.fiducial.lpa;
% rpa_vox = mri_trash.cfg.fiducial.rpa;
% 
% % bring the fiducials from voxel-space to headspace
% nas = ft_warp_apply(transm, nas_vox, 'homogenous');
% lpa = ft_warp_apply(transm, lpa_vox, 'homogenous');
% rpa = ft_warp_apply(transm, rpa_vox, 'homogenous');
%% 
%% check if the mri has been scaled
% Compute the scaling factors
scaling_factors = sqrt(sum(transm(1:3, 1:3).^2, 1));

% Display the scaling factors
disp('Scaling factors:');
disp(scaling_factors);
%% scale the elec to MRI?
elec_unaligned.elecpos = ...
    elec_unaligned.elecpos./scaling_factors;
elec_unaligned.chanpos = elec_unaligned.elecpos;
%% now we can align based on the fiducial positions

elec_mri.elecpos = [
  nas
  lpa
  rpa
  ];
elec_mri.label = {'na', 'lpa', 'rpa'};
elec_mri.unit  = 'mm';

% coregister the electrodes to the MRI using fiducials
cfg = [];
cfg.method   = 'fiducial';
cfg.template = elec_mri;
cfg.elec     = elec_unaligned;
cfg.feedback = 'yes';
cfg.fiducial = {'na', 'lpa', 'rpa'};
%red before /magenta after /blue mri reference
elec_aligned = ft_electroderealign(cfg);


%% plot together

figure;
ft_plot_mesh(hdm.bnd,'facecolor','none'); 

hold on;
ft_plot_sens(elec_aligned, 'facecolor', [0 1 0], 'edgecolor', [1 0 0]);

elec_fin = elec_aligned;
elec_fin.label = elec_aligned.label(...
    (ismember(elec_aligned.label, hdr.label)));
elec_fin.elecpos     = elec_aligned.elecpos(...
    (ismember(elec_aligned.label, hdr.label)),:);
elec_fin.chanpos = elec_fin.elecpos;
elec_fin.chantype     = elec_aligned.chantype(...
    (ismember(elec_aligned.label, hdr.label)));
elec_fin.chanunit     = elec_aligned.chanunit(...
    (ismember(elec_aligned.label, hdr.label)));
%%
% alternative
% % interactively coregister the electrodes to the BEM head model
% % this is a visual check and refinement step
% cfg = [];
% cfg.method    = 'interactive';
% cfg.elec      = elec_aligned;
% cfg.headshape = hdm.bnd(1);
% elec_new = ft_electroderealign(cfg);

% alternative: align mri to electrodes
% elec_unaligned.pos = elec_unaligned.elecpos;
% cfg = [];
% cfg.coordsys = 'ras'
% cfg.method = 'headshape';
% cfg.headshape.icp            = 'yes';
% cfg.headshape.headshape = elec_new;
% mri_aligned = ft_volumerealign(cfg, mri_orig,elec_new); % don't use the output

%% NOTE the reslicing may be necessary if the voxels  are not aligned with head axes
% 
cfg = [];
cfg.dim = mri_orig.dim;
mri_reslice = ft_volumereslice(cfg, mri_orig);

%% segmentation of MRI into scalp, skull, brain surfaces
cfg           = [];
% FEM ONLY
% cfg.output    = {'gray','white','csf','skull','scalp'};
cfg.output = {'scalp', 'skull', 'brain'};
segmented_mri  = ft_volumesegment(cfg, mri_reslice);

%% % plot segmented tissues
% create indexed representations
seg_i = ft_datatype_segmentation(segmented_mri,'segmentationstyle','indexed');
seg_i.seg = seg_i.tissue;

% FEM  only
%seg_i.seglabel  = {'gray','white','csf','skull','scalp'};

seg_i.seglabel = {'scalp', 'skull', 'brain'};
cfg              = [];
cfg.funparameter = 'seg';
cfg.funcolormap  = lines(numel(seg_i.seglabel)+1); % distinct color per tissue
% cfg.location     = 'center';
cfg.atlas        = seg_i;    % the segmentation can also be used as atlas
cfg.location = 'center';
ft_sourceplot(cfg, seg_i);
% saveas(gcf, [fig_path filesep sj_code '_segmentation.png'])
%% create a BEM 

cfg = [];
cfg.tissue      = {'scalp', 'skull', 'brain'}; %{'brain', 'skull', 'scalp'};
cfg.numvertices = [3000 6000 9000];
mesh = ft_prepare_mesh(cfg, segmented_mri);

% save mesh_BEM mesh
disp(mesh(1))

%% plot everything one more time
figure

ft_plot_ortho(mri_reslice.anatomy, 'style', 'intersect', ...
    'transform', mri_reslice.transform)
hold on;
ft_plot_mesh(mesh(1), 'facecolor', 'k','facealpha', 0.2); % scalp
view([0 -1 0]); % from the right side
ft_plot_mesh(mesh(2), 'facecolor', 'g', 'facealpha', 0.5); % skull
view([0 -1 0]); % from the right side

ft_plot_mesh(mesh(3), 'facecolor', 'b', 'facealpha', 0.8); %brain
ft_plot_sens(elec_fin, 'facecolor', ...
    [1 0 0], 'edgecolor', [1 0 0])
%% create boundary element method forward model

% Create a volume conduction model
cfg        = [];
cfg.method = 'dipoli'; % You can also specify 'openmeeg', 'bemcp', or another method
vol  = ft_prepare_headmodel(cfg, mesh);

save vol_dipoli vol -v7.3

%% fine tune the conductivity based on Meta-Analysis? 
% The weighted average mean and standard deviation (in S/m) 
% for each tissue type were: 
% --> scalp = 0.41 ± 0.18, 
% --> whole skull = 0.02 ± 0.02, 
% --> spongiform skull layer = 0.048 ± 0.07, 
% ---> whole compact skull layer = 0.005 ± 0.002, 
% ---> outer compact = 0.005 ± 0.003, 
% ----> inner compact = 0.007 ± 0.004, 
% ---> CSF = 1.71 ± 0.3, 
% ---> GM = 0.47 ± 0.24, WM = 0.22 ± 0.17, 
% ---> WM perpendicular = 0.12 ± 0.05,
% ---> WM parallel = 0.12 ± 0.09, 
% ---> Blood = 0.57 ± 0.11 and BSCR = 50.4 ± 39.

% dipoli says skin skull brain [0.3300 0.0042 0.3300]
cfg        = [];
cfg.conductivity = [0.41 ,0.02, 0.33];
cfg.method = 'dipoli'; % You can also specify 'openmeeg', 'bemcp', or another method
vol_improved  = ft_prepare_headmodel(cfg, mesh);
save vol_dipoli_improved vol_improved -v7.3
%% %% (FEM) mesh (3d model of each tissue layer in the brain)
% cfg             = [];
% cfg.shift       = 0.3;
% cfg.method      = 'hexahedral';
% mesh = ft_prepare_mesh(cfg, segmented_mri);
% save the mesh for the server

%% final check on the elec alignment in interactive mode
% 
cfg          = [];
cfg.method   = 'interactive';
cfg.elec     = elec_aligned;
cfg.headshape = mesh;
elec_trash = ft_electroderealign(cfg);

%% doesn't work at the moment
% %% create the forward model (model of the volume via Finite Element )
% %load the big mesh
% 
% cfg                 = [];
% cfg.method          = 'simbio';
% cfg.conductivity    = [0.33 0.14 1.79 0.01 0.43];   % order follows mesh.tissuelabel
% vol = ft_prepare_headmodel(cfg, mesh);   
%
%% also prepare the grid
%make a subject specific grid from the inverse warped template
%% use the warped template grid
% we use a grid defined in MNI space. It is inverse-warped onto the subject
% MRI. We can use a predefined standard sourcemodel: http://www.fieldtriptoolbox.org/template/sourcemodel/
ftpath = fileparts(which('ft_defaults'));
load(fullfile(ftpath, 'template/sourcemodel/standard_sourcemodel3d10mm'))
template_grid = ft_convert_units(sourcemodel,'mm');
clear sourcemodel;

% use the resliced MRI and warp the grid to it. 
cfg = [];
cfg.warpmni             = 'yes';
cfg.template            = template_grid;
cfg.nonlinear           = 'yes'; % use non-linear normalization
cfg.mri                 = mri_reslice;
cfg.sourcemodel.unit    = 'mm';
grid = ft_prepare_sourcemodel(cfg);   % ras + n

%%
ft_plot_mesh(mesh(3), 'facecolor', 'k',...
    'facealpha', 0.2); % scalp
hold on;
ft_plot_mesh(grid.pos(grid.inside,:), ...
    'vertexcolor', 'r');
%%
% save grid grid


%% or prepare an individual grid for this person only
% Prepare the source model grid
cfg = [];
cfg.grid.resolution = 10; % in mm
cfg.grid.unit = 'mm';
cfg.tight = 'yes'; % ensures grid points are inside the brain
cfg.inwardshift = 2; % shift grid points slightly inward to ensure they are inside the brain
cfg.headmodel = vol_improved;
grid = ft_prepare_sourcemodel(cfg);
%% plot individual grid
ft_plot_mesh(mesh(3), 'facecolor', 'k',...
    'facealpha', 0.2); % scalp
hold on;
ft_plot_mesh(grid.pos(grid.inside,:), ...
    'vertexcolor', 'r');
%% or prepare a template grid with inward shift and warp it to this persons space
template = load(fullfile(ftpath, 'template/headmodel/standard_bem'));

% Prepare the source model grid
cfg = [];
cfg.grid.resolution = 10; % in mm
cfg.grid.unit = 'mm';
cfg.tight = 'yes'; % ensures grid points are inside the brain
cfg.inwardshift = 5; % shift grid points 5mm inward
cfg.headmodel = template.vol;
template_grid = ft_prepare_sourcemodel(cfg);
%% plot template grid
ft_plot_mesh(template.vol.bnd(3), 'facecolor', 'k',...
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
ft_plot_mesh(vol_improved.bnd(3), 'facecolor', 'k',...
    'facealpha', 0.2); % brain
hold on;
ft_plot_mesh(grid.pos(grid.inside,:), ...
    'vertexcolor', 'r');
%%

ft_plot_ortho(mri_reslice.anatomy, 'style', 'intersect', ...
    'transform', mri_reslice.transform)
hold on;
ft_plot_mesh(mesh(1), 'facecolor', 'k', 'facealpha', 0.2); % brain
view([0 -1 0]); % from the right side
ft_plot_mesh(mesh(2), 'facecolor', 'g', 'facealpha', 0.2); % skull
view([0 -1 0]); % from the right side

ft_plot_mesh(mesh(3), 'facecolor', 'b', 'facealpha', 0.2); 
ft_plot_sens(elec_fin,'facecolor', 'g');

ft_plot_mesh(grid.pos(grid.inside,:), ...
    'vertexcolor', [1 0 0]);

%% SLIDES