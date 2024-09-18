clear; close all; clc;

% Initialize Fieldtrip
addpath('/Users/mrugankdake/Documents/MATLAB/fieldtrip-20220104');
ft_defaults;
addpath('/Users/mrugankdake/Documents/MATLAB/fieldtrip-20220104/external/freesurfer') 
addpath('/Users/mrugankdake/Documents/MATLAB/fieldtrip-20220104/external/spm12') 

%% Load sample dataset
% For EEG dataset, we will use CHB-MIT Scalp EEG Dataset from here: 
% https://physionet.org/content/chbmit/1.0.0/chb01/#files-panel

% The original dataset contains 22 subjects. However, here in this 
% tutorial we will explore data from Subject 01

% You can download the dataset from Brightspace.

%% Filters as matrix multiplications
%% Load EEG data

% Path to EEG file
eegPath = '../../Datasets/EEG/sub-01/eeg/sub-01_task-daf_eeg_filtered.vhdr';

elecPosFile = '../../Datasets/EEG/sub-01/eeg/sub-01_electrodes.tsv';
elecPos = readtable(elecPosFile, 'FileType', 'text', 'Delimiter', '\t');
% Add fiducials
fiducials.name      = {'Nz', 'LPA', 'RPA'}';
fiducials.x         = [-4.129838157917329e-18, -0.0729282673627754, 0.08278152042487033]';
fiducials.y         = [0.10011015398430487, 3.008505424862354e-18, -3.414981080487009e-18]';
fiducials.z         = [-5.7777898331617076e-33, 3.851859888774472e-34, 3.4666738998970245e-33]';
fiducials.impedance = {'n/a', 'n/a', 'n/a'}';
fiducialsTable      = struct2table(fiducials);
elecPos             = [elecPos; struct2table(fiducials)];

elec                = [];
elec.label          = elecPos.name;
elec.elecpos        = [elecPos.x, elecPos.y, elecPos.z];
elec.chanpos        = elec.elecpos;
elec.unit           = 'm';
% Load the EEG file using Fieldtrip
% Fieldtrip requires 'cfg' structure to load the data
cfg                 = [];
cfg.dataset         = eegPath;
% Add electrode positions
cfg.elec            = elec;
rawData             = ft_preprocessing(cfg);

data                = rawData.trial{1};
n_channels          = size(data, 1);
n_samples           = size(data, 2);

%% Select an electrode via vectro multiplication
% Select electrode index
selectElecIdx = 2;

% Create a vector with a 1 at the selected index
selectElecVec = zeros(1, n_channels);
selectElecVec(selectElecIdx) = 1;

% Select the data for the chosen electrode
selectElecData = selectElecVec * data;

disp(selectElecVec);

%% Select a subset of electrodes via matrix multiplication
% Let's select channels 1, 3, and 5
selected_channels = [1, 3, 5]; 

% Create a matrix with selected channels
ChannelMultiplier = zeros(n_channels, n_channels);
for i = selected_channels
    ChannelMultiplier(i, i) = 1;
end

disp(ChannelMultiplier);

% Apply the ChannelMultiplier to the data
data_selected = ChannelMultiplier * data;

figure();
imagesc(ChannelMultiplier)


%% Average an electrode-ROI via vector multiplication
% Let's average across occipital electrodes
elecList = {'O1', 'Oz', 'O2', 'PO3', 'POz', 'PO4'};

% Assuming elecPos is a table with electrode names and their indices
% Get indices for the specified electrodes
elecIdx = zeros(1, numel(elecList));
for i = 1:numel(elecList)
    elecIdx(i) = find(strcmp(elecPos.name, elecList{i}));
end

% Create a vector with 1s at the selected indices
elecVec = zeros(1, n_channels);
elecVec(elecIdx) = 1;

% Average the data across the specified electrodes
elecData = elecVec * data;

disp(elecVec);
disp(size(elecData));

%% Channel interpolation via matrix multiplication
% Let's interpolate channel 2 using channels 1 and 3
interp_channel = 2; 

InterpMultiplier = eye(n_channels);
InterpMultiplier(interp_channel, interp_channel - 1) = 0.5;
InterpMultiplier(interp_channel, interp_channel + 1) = 0.5;
InterpMultiplier(interp_channel, interp_channel) = 0;

disp(InterpMultiplier);

data_interp = InterpMultiplier * data;

figure();
imagesc(InterpMultiplier)

% Plot the interpolated EEG data
% n_rows = ceil(n_channels / 4); % Adjust number of rows and columns as needed
% n_cols = min(4, n_channels); % Adjust number of columns as needed

% figure;
% for i = 1:n_channels
%     subplot(n_rows, n_cols, i);
%     plot(data_interp(i, :));
%     title(['Channel ', num2str(i)]);
%     xlabel('Time');
%     ylabel('Amplitude');
% end

%% Fixing swapped electrodes
% Channel 1 and 3 have been plugged into wrong position in the cap
% Swap the channels using matrix multiplication
SwapMultiplier = eye(n_channels);
SwapMultiplier(1, 3) = 1;
SwapMultiplier(3, 1) = 1;
SwapMultiplier(1, 1) = 0;
SwapMultiplier(3, 3) = 0;

disp(SwapMultiplier);

data_swapped = SwapMultiplier * data;

figure;
imagesc(SwapMultiplier);
colorbar;


%% Re-referencing data to a single channel using matrix multiplication
% Let's use channel 10 as the reference channel
ref_channel = 10;

ReferenceMultiplier = eye(n_channels);
ReferenceMultiplier(:, ref_channel) = -1;
ReferenceMultiplier(ref_channel, :) = 0;
ReferenceMultiplier(ref_channel, ref_channel) = 1;

disp(ReferenceMultiplier);

% Re-reference the data
data_single_ref = ReferenceMultiplier * data;

figure();
imagesc(ReferenceMultiplier)

% n_rows = ceil(n_channels / 4); % Adjust number of rows and columns as needed
% n_cols = min(4, n_channels); % Adjust number of columns as needed
% 
% figure;
% for i = 1:n_channels
%     subplot(n_rows, n_cols, i);
%     plot(data_single_ref(i, :));
%     title(['Channel ', num2str(i)]);
%     xlabel('Time');
%     ylabel('Amplitude');
% end


%% Average re-referencing via matrix multiplication
MeanMultiplier = eye(n_channels) - ones(n_channels, n_channels) / n_channels;

disp(MeanMultiplier);

data_re_ref = MeanMultiplier * data;

% % Plot the re-referenced EEG data
% n_rows = ceil(n_channels / 4); % Adjust number of rows and columns as needed
% n_cols = min(4, n_channels); % Adjust number of columns as needed
% 
% figure;
% for i = 1:n_channels
%     subplot(n_rows, n_cols, i);
%     plot(data_re_ref(i, :));
%     title(['Channel ', num2str(i)]);
%     xlabel('Time');
%     ylabel('Amplitude');
% end

% Plot the MeanMultiplier matrix
figure;
imagesc(MeanMultiplier);

