clear; close all; clc;

% Initialize Fieldtrip
addpath('/Users/mrugankdake/Documents/MATLAB/fieldtrip-20220104');
ft_defaults;

%% Load sample dataset
% For EEG dataset, we will use EEG Motor Movement/Imagery Dataset from
% here: https://physionet.org/content/eegmmidb/1.0.0/#files-panel

% The original dataset contains 109 volunteers. However, here in this
% tutorial we will explore data from subject 42

% You can download the dataset from Brightspace.

%% Load EEG data

% Path to EEG file
eegPath = '../../Datasets/EEG/S042/S042R01.edf';

% Load the EEG file using Fieldtrip
% Fieldtrip requires 'cfg' structure to load the data
cfg                 = [];
cfg.dataset         = eegPath;
cfg.continuous      = 'yes';
cfg.channel         = 'all';
data                = ft_preprocessing(cfg);

%% Plot EEG data
% Let us try plotting an EEG channel using default MATLAB plot function
% Let us plot the first channel
figure();
plot(data.time{1}, data.trial{1}(1,:));
xlabel('Time (s)');
ylabel('Amplitude (uV)');
title('EEG Channel 1');

% Now let's try using the interactive plot function from Fieldtrip
% This is done using ft_databrowser
% This will open a new window with interactive plot
% This is useful for exploring the data, checking for any artifacts, finding bad channels, or sections of recordings that need to be removed
cfg                 = [];
cfg.viewmode        = 'vertical';
ft_databrowser(cfg, data);

%% Rereferencing data
% EEG data is usually recorded with respect to a reference electrode. The reference electrode can be different for different EEG systems.
% The data can be re-referenced to a common reference electrode, such as average reference, linked mastoids, etc.
% Here we are going to re-reference the data using the average reference and the mastoids as reference electrodes
% The average refernce is calculated by taking the average of all the electrodes
cfg                 = [];
cfg.reref           = 'yes';
cfg.refchannel      = 'all'; % 'all' means all channels will be re-referenced
cfg.refmethod       = 'avg'; % 'avg' means average reference
dataRerefAvg        = ft_preprocessing(cfg, data);

% We can also reference data using a specific channel as reference
% Here we are going to use the mastoids as reference
cfg                 = [];
cfg.reref           = 'yes';
cfg.refchannel      = {'T7..', 'T8..'}; % Mastoids
cfg.refmethod       = 'avg';
dataRerefMastoids   = ft_preprocessing(cfg, data);

% Plot the data before and after re-referencing
figure();
subplot(3,1,1);
plot(data.time{1}, data.trial{1}(1,:));
xlabel('Time (s)');
ylabel('Amplitude (uV)');
title('EEG Channel 1 - Original');

subplot(3,1,2);
plot(dataRerefAvg.time{1}, dataRerefAvg.trial{1}(1,:));
xlabel('Time (s)');
ylabel('Amplitude (uV)');
title('EEG Channel 1 - Average Reference');

subplot(3,1,3);
plot(dataRerefMastoids.time{1}, dataRerefMastoids.trial{1}(1,:));
xlabel('Time (s)');
ylabel('Amplitude (uV)');
title('EEG Channel 1 - Mastoids Reference');