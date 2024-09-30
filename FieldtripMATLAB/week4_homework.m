clear; close all; clc;

% Initialize Fieldtrip
addpath('/Users/mrugankdake/Documents/MATLAB/fieldtrip-20220104');
ft_defaults;
addpath('/Users/mrugankdake/Documents/MATLAB/fieldtrip-20220104/external/freesurfer') 
addpath('/Users/mrugankdake/Documents/MATLAB/fieldtrip-20220104/external/spm12') 

%%        MORLET WAVELETS    
% Define the parameters for the Morlet wavelet
f0_10 = 10; % Central frequency of the wavelet
f0_15 = 15; % Central frequency of the wavelet

SR = 100;
sigma = 1; % Standard deviation of the Gaussian window
t = -5:1/SR:5; % Time vector
N_ml = length(t);

% Create the Morlet wavelet (%frequency based on time)
A10 = 1 / (sqrt(sigma * sqrt(pi)) * sqrt(f0_10)); % amplitude scaling for the frequency-specific Morlet wavelet
morlet_wavelet10 = A10*exp(2 * 1i * pi * f0_10 * t) .* exp(-t.^2 / (2 * sigma^2));

% Create the Morlet wavelet (%frequency based on time)
A15 = 1 / (sqrt(sigma * sqrt(pi)) * sqrt(f0_15)); % amplitude scaling for the frequency-specific Morlet wavelet
morlet_wavelet15 = A15*exp(2 * 1i * pi * f0_15 * t) .* exp(-t.^2 / (2 * sigma^2));

figure;
subplot(1, 2, 1)
title('Morlet Wavelet (10Hz)')
xlabel('Sample Index')
ylabel('Amplitude')
plot(t, real(morlet_wavelet10))
xlim([t(1) t(end)])

subplot(1, 2, 2)
title('Morlet Wavelet (15Hz)')
xlabel('Sample Index')
ylabel('Amplitude')
plot(t, real(morlet_wavelet15))
xlim([t(1) t(end)])

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

%% convolve your EEG data with a Morlet wavelet
signal = data(3, 1:10000);  

conv10 = conv(signal, morlet_wavelet10, 'same');
conv15 = conv(signal, morlet_wavelet15, 'same');

figure;
subplot(3, 1, 1);
plot(real(conv10)); 
title('Convolution with Morlet Wavelet (10Hz)');
xlabel('Sample Index');
ylabel('Amplitude');
xlim([0 length(signal)]);

subplot(3, 1, 2);
plot(real(conv15)); 
title('Convolution with Morlet Wavelet (15Hz)');
xlabel('Sample Index');
ylabel('Amplitude');
xlim([0 length(signal)]);

subplot(3, 1, 3);
plot(signal); 
title('Original Signal');
xlabel('Sample Index');
ylabel('Amplitude');
xlim([0 length(signal)]);

%% add a spike to your EEG data by setting a random value to the value 10000
% what happens after you convolve your EEG with the morlet wavelet?
signal(5000) = 10000;

conv10_spike = conv(signal, morlet_wavelet10, 'same');
conv15_spike = conv(signal, morlet_wavelet15, 'same');

figure;
subplot(3, 1, 1);
plot(real(conv10_spike)); 
title('Convolution with Morlet Wavelet (10Hz)');
xlabel('Sample Index');
ylabel('Amplitude');
xlim([0 length(signal)]);

subplot(3, 1, 2);
plot(real(conv15_spike)); 
title('Convolution with Morlet Wavelet (15Hz)');
xlabel('Sample Index');
ylabel('Amplitude');
xlim([0 length(signal)]);

subplot(3, 1, 3);
plot(signal); 
title('Original Signal with Spike');
xlabel('Sample Index');
ylabel('Amplitude');
xlim([0 length(signal)]);

