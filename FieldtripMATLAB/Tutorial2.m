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

% Try extracting data from rawData
data                = ...
n_channels          = size(data, 1);
n_samples           = size(data, 2);

%% Select an electrode via vectro multiplication
selectElecIdx = 2;

selectElecVec = zeros(1, n_channels);
% Try creating vector for the selected electrode

selectElecData = selectElecVec * data;

disp(selectElecVec);

%% Select a subset of electrodes via matrix multiplication
% Let's select channels 1, 3, and 5
selected_channels = [1, 3, 5]; 

ChannelMultiplier = zeros(n_channels, n_channels);
% Make the ChannelMultiplier select only the selected channels
...

disp(ChannelMultiplier);

data_selected = ChannelMultiplier * data;

figure();
imagesc(ChannelMultiplier)


%% Average an electrode-ROI via vector multiplication
% Let's average across occipital electrodes
elecList = {'O1', 'Oz', 'O2', 'PO3', 'POz', 'PO4'};

% Get indices for the specified electrodes
elecIdx = zeros(1, numel(elecList));
for i = 1:numel(elecList)
    elecIdx(i) = find(strcmp(elecPos.name, elecList{i}));
end

elecVec = zeros(1, n_channels);
% Try creating vector for the selected electrodes to average
...

elecData = elecVec * data;

disp(elecVec);
disp(size(elecData));

%% Channel interpolation via matrix multiplication
% Let's interpolate channel 2 using channels 1 and 3
interp_channel = 2; 

InterpMultiplier = eye(n_channels);
% Make the InterpMultiplier interpolate the selected channel
...

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
% Make the SwapMultiplier swap the selected channels
...

disp(SwapMultiplier);

data_swapped = SwapMultiplier * data;

figure;
imagesc(SwapMultiplier);
colorbar;

%% Re-referencing data to a single channel using matrix multiplication
% Let's use channel 10 as the reference channel
ref_channel = 10;

ReferenceMultiplier = eye(n_channels);
% Make the ReferenceMultiplier re-reference the data to the selected channel
...

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
MeanMultiplier = ...

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

