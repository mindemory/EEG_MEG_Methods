clear; close all; clc;

% Initialize Fieldtrip
addpath('/Users/mrugankdake/Documents/MATLAB/fieldtrip-20220104');
ft_defaults;
addpath('/Users/mrugankdake/Documents/MATLAB/fieldtrip-20220104/external/freesurfer') 
addpath('/Users/mrugankdake/Documents/MATLAB/fieldtrip-20220104/external/spm12') 

% addpath([toolpth_ filesep 'fieldtrip' filesep 'external' ...
%     filesep 'spm12']);

%% Load sample dataset
% For EEG dataset, we will use CHB-MIT Scalp EEG Dataset from here: 
% https://physionet.org/content/chbmit/1.0.0/chb01/#files-panel

% The original dataset contains 22 subjects. However, here in this 
% tutorial we will explore data from Subject 01

% You can download the dataset from Brightspace.

%% Load EEG data

% Path to EEG file
eegPath = '../../Datasets/EEG/sub-01/eeg/sub-01_task-daf_eeg_filtered.vhdr';
% 
% elecPosFile = '../../Datasets/iEEG/sub-01/eeg/sub-01_electrodes.tsv';
% elecPos = readtable(elecPosFile, 'FileType', 'text', 'Delimiter', '\t');
% % Add fiducials
% fiducials.name      = {'Nz', 'LPA', 'RPA'}';
% fiducials.x         = [-4.129838157917329e-18, -0.0729282673627754, 0.08278152042487033]';
% fiducials.y         = [0.10011015398430487, 3.008505424862354e-18, -3.414981080487009e-18]';
% fiducials.z         = [-5.7777898331617076e-33, 3.851859888774472e-34, 3.4666738998970245e-33]';
% fiducials.impedance = {'n/a', 'n/a', 'n/a'}';
% fiducialsTable      = struct2table(fiducials);
% elecPos             = [elecPos; struct2table(fiducials)];
% 
% elec                = [];
% elec.label          = elecPos.name;
% elec.elecpos        = [elecPos.x, elecPos.y, elecPos.z];
% elec.chanpos        = elec.elecpos;
% elec.unit           = 'm';
% Load the EEG file using Fieldtrip
% Fieldtrip requires 'cfg' structure to load the data
cfg                 = [];
cfg.dataset         = eegPath;
% Add electrode positions
% cfg.elec            = elec;
data                = ft_preprocessing(cfg);

%% Plot EEG data (single channel plotting)
% Let us try plotting an EEG channel using default MATLAB plot function
% Let us plot the first channel
chToPlot = 10; % Channel to plot
timeToPlot = 5; % Time to plot in seconds
samplesToPlot = timeToPlot * data.fsample;
figure();
plot(data.time{1}(1:samplesToPlot), ...
    data.trial{1}(chToPlot, 1:samplesToPlot))
xlabel('Time (s)');
ylabel('Amplitude (uV)');
title('EEG Channel 1');

%% Plot EEG Data (databrowser)
% Now let's try using the interactive plot function from Fieldtrip
% This is done using ft_databrowser
% This will open a new window with interactive plot
% This is useful for exploring the data, checking for any artifacts, finding bad channels, or sections of recordings that need to be removed
cfg                 = [];
cfg.viewmode        = 'vertical';
ft_databrowser(cfg, data);

%% Visualizing Montage of electrodes
% The EEG data is recorded using multiple electrodes placed on the scalp.
% The arrangement of electrodes is called a montage.
% The montage can be visualized using ft_plot_sens
% But let us try visualizing the montage as a 3D plot usint MATLAB's plot3 function
figure();
plot3(data.elec.elecpos(:, 1), data.elec.elecpos(:, 2), data.elec.elecpos(:, 3), 'bo', 'MarkerFaceColor','b');
hold on;
text(data.elec.elecpos(:, 1), data.elec.elecpos(:, 2), data.elec.elecpos(:, 3), data.elec.label, 'FontSize', 8);

plot3(fiducialsTable.x, fiducialsTable.y, fiducialsTable.z, 'ro', 'MarkerFaceColor', 'r');
xlabel('X');
ylabel('Y');
zlabel('Z');
title('EEG Electrode Montage');
view(180, 0)

%% Rereferencing data
% EEG data is usually recorded with respect to a reference electrode. The reference electrode can be different for different EEG systems.
% The data can be re-referenced to a common reference electrode, such as average reference, linked mastoids, etc.
% Here we are going to re-reference the data using the average reference and the mastoids as reference electrodes
% The average refernce is calculated by taking the average of all the electrodes
dataRerefAvg           = data;
dataRerefAvg.trial{1}  = data.trial{1} - mean(data.trial{1}, 1);


% We can also reference data using a specific channel as reference
% Here we are going to use the mastoids as reference
mastoidIdx            = find(ismember(data.label, 'TP9') | ismember(data.label, 'TP10'));
dataRerefMastoids     = data;
dataRerefMastoids.trial{1} = dataRerefMastoids.trial{1} - (dataRerefMastoids.trial{1}(mastoidIdx(1), :) + ...
                                                           dataRerefMastoids.trial{1}(mastoidIdx(2), :))/2;

% Plot the data before and after re-referencing
figure();
chToPlot = 10; % Channel to plot
timeToPlot = 5; % Time to plot in seconds
samplesToPlot = timeToPlot * data.fsample;

subplot(3,1,1);
plot(data.time{1}(1:samplesToPlot), ...
    data.trial{1}(chToPlot,1:samplesToPlot));
xlabel('Time (s)');
ylabel('Amplitude (uV)');
title(['EEG Channel ' num2str(chToPlot) ' - Original']);

subplot(3,1,2);
plot(dataRerefAvg.time{1}(1:samplesToPlot), ...
    dataRerefAvg.trial{1}(chToPlot, 1:samplesToPlot));
xlabel('Time (s)');
ylabel('Amplitude (uV)');
title(['EEG Channel ' num2str(chToPlot) ' - Average Reference']);

subplot(3,1,3);
plot(dataRerefMastoids.time{1}(1:samplesToPlot), ...
    dataRerefMastoids.trial{1}(chToPlot, 1:samplesToPlot));
xlabel('Time (s)');
ylabel('Amplitude (uV)');
title(['EEG Channel ' num2str(chToPlot) ' - Mastoids Reference']);


%% Load MEG data
% Path to MEG file
megPath              = '../../Datasets/MEG/Subject01/Subject01.ds';

% Load the MEG file using Fieldtrip
cfg                 = [];
cfg.dataset         = megPath;
megData             = ft_preprocessing(cfg);

% Extract the channel positions and the coil positions
% The channel positions refer to the location of the physical sensors in
% MEG
% And the coil positions are the locations of the gradiometers in the MEG
% setup
chanpos             = megData.grad.chanpos; % Channel positions
coilpos             = megData.grad.coilpos; % Coil positions

figure;
hold on;
plot3(chanpos(:,1), chanpos(:,2), chanpos(:,3), 'bo', 'MarkerSize', 8, 'DisplayName', 'Channels');
plot3(coilpos(:,1), coilpos(:,2), coilpos(:,3), 'ro', 'MarkerSize', 4, 'DisplayName', 'Coils');
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Visualization of Channel and Coil Positions');
legend;
view(3);
axis equal;
grid on;


% We can also visualize these channels with respect to the subject head if
% we have the anatomicals of the subject.
% Load and segment MRI (anatomical)
mriPath              = '../../Datasets/MEG/Subject01/Subject01.mri';
mri                  = ft_read_mri(mriPath);

output_nii_path = '../../Datasets/MEG/Subject01/Subject01.nii';  % Define the output NIfTI file name
ft_write_mri(output_nii_path, mri, 'dataformat', 'nifti');

% Visualize MRI 
cfg                  = [];
cfg.method           = 'ortho';
cfg.interactive      = 'yes';
ft_sourceplot(cfg, mri);

% First we segment the anatomical to identify different layers
segmentedmri         = ft_volumesegment([], mri);
% add anatomical information to the segmentation
segmentedmri.transform = mri.transform;
segmentedmri.anatomy   = mri.anatomy;

% Next we construct a headmodel using the segmented MRI
cfg        = [];
cfg.method = 'singleshell';
hdm        = ft_prepare_headmodel(cfg, segmentedmri);
% Convert the units to the units present in megData.grad
hdm        = ft_convert_units(hdm, 'cm');

% Now we can visualize the mesh of subject head along with grad
figure;
ft_plot_mesh(hdm.bnd,'facecolor','none'); 
hold on;
ft_plot_sens(megData.grad, 'facecolor', [0 1 0], 'edgecolor', [1 0 0]);
view([0 0]);