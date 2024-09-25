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

%% Select an electrode via vector multiplication
selectElecIdx = 2;

% Create a vector with a 1 at the selected index
selectElecVec = zeros(1, n_channels);
selectElecVec(selectElecIdx) = 1;
disp(size(selectElecVec))
disp(size(data))

% Select the data for the chosen electrode
selectElecData = selectElecVec * data;

disp(size(selectElecData))
disp(selectElecData(1:10))

%% Select a subset of electrodes via matrix multiplication
% Let's select channels 1, 3, and 5
selected_channels = [1, 3, 5]; 

% Create a matrix with selected channels
ChannelMultiplier = zeros(length(selected_channels), n_channels);
for i = 1:length(selected_channels)
    chIdx = selected_channels(i);
    ChannelMultiplier(i, chIdx) = 1;
end
disp(size(ChannelMultiplier));
disp(size(data));

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
elecVec = zeros(length(elecIdx), n_channels);
for i = 1:length(elecIdx)
    elec = elecIdx(i);  
    elecVec(i, elec) = 1 / length(elecIdx);  
end

% Average the data across the specified electrodes
elecData = elecVec * data;

disp(size(elecVec));
disp(size(data));

figure();
imagesc(elecVec)

%% Channel interpolation via matrix multiplication
% Let's interpolate channel 31 using channels 30 and 32
interp_channel = 31;  
n_channels = size(data, 1);  % Number of channels

% Initialize InterpMultiplier as an identity matrix
InterpMultiplier = eye(n_channels);

% Set the weights for the interpolation
InterpMultiplier(interp_channel, interp_channel-1) = 1/2;  % Channel 30
InterpMultiplier(interp_channel, interp_channel+1) = 1/2;  % Channel 32
InterpMultiplier(interp_channel, interp_channel) = 0;      % Zero out the original channel 31

disp(size(InterpMultiplier));
disp(size(data));

% Apply the interpolation via matrix multiplication
data_interp = InterpMultiplier * data;

% Plot the InterpMultiplier matrix
figure;
imagesc(InterpMultiplier);
axis image;
title('Interpolation Multiplier Matrix');
colorbar;

% Plot the original EEG data for channels 30, 31, and 32
figure;
for i = 1:3
    chan = interp_channel - 2 + i;  % Channels 30, 31, 32 in MATLAB (1-based indexing)
    subplot(3, 1, i);
    plot(data(chan, 1:2000));  % Plot first 2000 time points
    title(['Channel ', num2str(chan)]);
    xlabel('Time');
    ylabel('Amplitude');
end

% Plot the interpolated EEG data for channels 30, 31, and 32
figure;
for i = 1:3
    chan = interp_channel - 2 + i;  % Channels 30, 31, 32
    subplot(3, 1, i);
    plot(data_interp(chan, 1:2000));  % Plot first 2000 time points
    title(['Channel ', num2str(chan)]);
    xlabel('Time');
    ylabel('Amplitude');
end

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

%% Average Re-referencing via matrix multiplication
n_channels = size(data, 1);  % Number of channels

% Create averages_mat which is a row vector with each element as 1/n_channels
averages_mat = ones(1, n_channels) / n_channels;

% Perform the average re-referencing operation
data_re_ref = (eye(n_channels) * data) - (averages_mat * data);  
% Matrix multiplication is distributive
% and hence we can do referencing as: 
% A * X - B * X
% OR
% (A - B) * X
data_re_ref = (eye(n_channels) - averages_mat) * data;  

% Create the MeanMultiplier matrix for later use
MeanMultiplier = eye(n_channels) - ones(n_channels, n_channels) / n_channels;
disp(MeanMultiplier);

% Plot the re-referenced EEG data
% figure;
% nrows = ceil(sqrt(n_channels));  
% ncols = ceil(n_channels / nrows); 
% 
% for i = 1:n_channels
%     subplot(nrows, ncols, i);  
%     plot(data_re_ref(i, :));  % Plot the re-referenced data for channel i
%     title(['Channel ', num2str(i)]);
%     xlabel('Time');
%     ylabel('Amplitude');
% end

% Plot the MeanMultiplier matrix as an image
figure;
imagesc(MeanMultiplier);
axis image;
title('Mean Multiplier Matrix');
colorbar;

%% Creating a fake channel with a polynomial trend
fake_data = data(1, 1:10000);  % Assuming 'raw' has similar data structure
t = linspace(0, 1, length(fake_data));
rng(42);  
trend = randn(1, 3);  % Random coefficients for polynomial trend
fake_data = fake_data + 4e+3 * polyval(trend, t);

% Remove trend by fitting polynomials of different orders

% 0th order polynomial (mean-centering)
trend0 = polyfit(t, fake_data, 0);
fake_data_centered = fake_data - polyval(trend0, t);

% 1st order polynomial (linear detrending)
trend1 = polyfit(t, fake_data, 1);
fake_data_detrended = fake_data - polyval(trend1, t);

% 2nd order polynomial (quadratic detrending)
trend2 = polyfit(t, fake_data, 2);
fake_data_detrended2 = fake_data - polyval(trend2, t);

% Plotting the original and detrended data
figure;
subplot(2, 2, 1);
plot(t, fake_data);
title('Original Data');
xlabel('Time');
ylabel('Amplitude');

subplot(2, 2, 2);
plot(t, fake_data_centered);
title('Mean-Centered Data');
xlabel('Time');
ylabel('Amplitude');

subplot(2, 2, 3);
plot(t, fake_data_detrended);
title('Linear Detrended Data');
xlabel('Time');
ylabel('Amplitude');

subplot(2, 2, 4);
plot(t, fake_data_detrended2);
title('Polynomial Detrended Data');
xlabel('Time');
ylabel('Amplitude');

%% Compute power at 10 Hz for occipital electrode O2
timeToPlot = 50; % Time in seconds to plot
tArr = rawData.time{1};
t = tArr(tArr < timeToPlot);
N = length(t);
chName = 'O2';

% Assuming elecPos is a table or structure where the electrode positions are stored
elecIdx = find(strcmp(elecPos.name, chName));

x = (0:N-1) / N;
X = data(elecIdx, tArr < timeToPlot);

SR = rawData.fsample;  % Sampling rate
freqs = (0:(N-1)) * SR / N;
freqs = freqs(1:N/2);
k = 10;

% Create a complex sine wave at frequency k
W = exp(-1i * 2 * pi * k * x);

% Compute the inner product (DFT coefficient for frequency k)
X_k = W * X.';

% Compute power at frequency k
power = real(X_k)^2 + imag(X_k)^2;
fprintf('Power at %d Hz: %f uV^2\n', k, power);


%% Let us compute power at all frequencies for occipital electrode O2
timeToPlot = 50; % Time in seconds to plot
tArr = rawData.time{1};
t = tArr(tArr < timeToPlot);
N = length(t);
chName = 'O2';

% Assuming elecPos is a table or structure where the electrode positions are stored
elecIdx = find(strcmp(elecPos.name, chName));

x = (0:N-1) / N;
X = data(elecIdx, tArr < timeToPlot);

SR = rawData.fsample;  % Sampling rate
freqs = (0:(N-1)) * SR / N;
freqs = freqs(1:N/2);
powerVec = zeros(length(freqs));
for i = 1:length(freqs)
    k = freqs(i);
    % Create a complex sine wave at frequency k
    W = exp(-1i * 2 * pi * k * x);

    % Compute the inner product (DFT coefficient for frequency k)
    X_k = W * X.';

    powerVec(i) =  real(X_k).^2 + imag(X_k).^2;
end

figure();
subplot(2, 1, 1)
plot(t, X)
xlabel('Time (s)')
ylabel('Amplitude (uV)')
title(['Channel: ' chName])

subplot(2, 1, 2)
plot(freqs, powerVec)
xlabel('Frequency (Hz)')
ylabel('Power')

%% Compute the 8-12Hz power at all channels using matrix multiplication
timeToPlot = 500; % Time in seconds to plot
tArr = rawData.time{1};

t = tArr(tArr < timeToPlot);
N = length(t);
x = (0:N-1) / N;
X = data(:, tArr < timeToPlot);

SR = rawData.fsample;
freq_res = SR / length(X);
freqs = (0:N-1)' * (1/SR);
freqs = freqs(1:N/2);
k_range = 8:12;

powerMat = zeros(n_channels, length(k_range));
for i = 1:length(k_range)
    k = k_range(i);
    % Create a complex sine wave at frequency k
    W = exp(-1i * 2 * pi * k * x);

    % Compute the inner product (DFT coefficient for frequency k)
    % We use matrix multiplication!
    X_k = (W * X')'; 

    powerMat(:, i) = real(X_k).^2 + imag(X_k).^2;
end

powerVec = mean(powerMat, 2);

% Plotting
figure(); 
plot(powerVec(1:end-4), 'o--')
xlabel('Channel');
ylabel('Power');
xticklabels(rawData.label(1:end-4))
title(['Power at ', num2str(k_range), ' Hz']);
xtickangle(90); 
grid on; 

%% We can visualize this data as a topolplot in Fieldtrip
% Exclude the non-EEG channels
nonEEG_channels = {'leog', 'reog', 'egg', 'audio'};
picks = find(~ismember(rawData.label, nonEEG_channels)); % Get indices of EEG channels

% Select alpha power for only EEG channels
alpha_power = powerVec; 
alpha_power_eeg = alpha_power(picks); 

% Prepare data for topoplot
cfg = [];
cfg.layout = rawData.elec; % Use the electrode layout from rawData.elecs
cfg.zlim = 'maxabs'; % Color limits
cfg.comment = 'no'; 
cfg.colormap = 'RdBu_r';

% Create a structure for the topoplot
topoData = [];
topoData.label = rawData.label(picks); 
topoData.avg = alpha_power_eeg; 
topoData.dimord = 'chan'; 
topoData.time = 0; 

% Plot the topomap for alpha power (only using the EEG channels)
figure;
ft_topoplotER(cfg, topoData);
colorbar; 
title('Alpha Power Topomap (8-12 Hz)');



