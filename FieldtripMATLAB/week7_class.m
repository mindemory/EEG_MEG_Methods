clear; close all; clc;

% You are running your code on HPC
curr_dir = pwd;
path_parts = strsplit(curr_dir, filesep);
fieldtrip_path = fullfile(filesep, path_parts{1:end-1}, 'fieldtrip');
dataPath = fullfile(filesep, 'scratch', 'work', 'courses', ...
    'PSYCH-GA-3405-2024fa');

meg_path = fullfile(dataPath, 'MEG');
eeg_path = fullfile(dataPath, 'EEG');

% Add fieldtrip to path
addpath(fullfile(fieldtrip_path));
ft_defaults;
%%
groupName = 'GroupB';
recRoot = fullfile(meg_path, groupName, 'Recording');

lstFiles = dir(recRoot);
lstFiles = {lstFiles.name};

% Filter the files based on their extensions
sqdFiles = lstFiles(endsWith(lstFiles, '_NR.sqd'));
mrkFiles = lstFiles(endsWith(lstFiles, '.mrk'));
elpFiles = lstFiles(endsWith(lstFiles, 'Points.txt'));
hsFiles = lstFiles(endsWith(lstFiles, 'HS.txt'));
matFiles = lstFiles(endsWith(lstFiles, '.mat'));

% In this case, second file is the one we want to load
% this will vary depending on which group's data you are analyzing and how
% many files you have.
ff = sqdFiles{2};
rec_filepath = fullfile(recRoot, ff);

% % Place holder to load marker, fiducials and headshape files
% mrk = cellfun(@(x) fullfile(recRoot, x), mrkFiles, 'UniformOutput', false);
% elp = fullfile(recRoot, elpFiles{1});
% hsp = fullfile(recRoot, hsFiles{1});

cfg = [];
cfg.dataset = fullfile(rec_filepath);
all_data = ft_preprocessing(cfg);
%%  identify bad channels here
%@mrugank / todo

cfg = [];
cfg.channel = all_data.label(1:157);
all_data = ft_selectdata(cfg, all_data);

cfg = [];
cfg.resamplefs = 256;
all_data = ft_resampledata(cfg, all_data);

% load in artifactual sampls
dummy = [1:10; 158]
%% inspect the data
cfg = [];
cfg.blocksize = 10;
cfg.preproc.hpfilter = 'yes';
cfg.preproc.hpfreq = 0.5;
ft_databrowser(cfg, all_data);
%% do a pca manually
data = all_data.trial{1}(1:157,:);

N = size(data,2);
% center the data
data = bsxfun(@minus, data, nanmean(data, 2));
% Compute the covariance matrix
covMatrix = (data * data')./(N-1);
covMatrix2 = cov(data');

%% plot the covariance matrix:
figure;
imagesc(covMatrix);
caxis([-max(abs(covMatrix(:))) max(abs(covMatrix(:)))]);colorbar;
colormap(jet(256));
%% PCA
% Perform eigenvalue decomposition
[eigenVectors, eigenValues] = eig(covMatrix);

% Sort the eigenvalues and corresponding eigenvectors in descending order
[eigenValuesSorted, sortOrder] = sort(diag(eigenValues), 'descend');
eigenVectorsSorted = eigenVectors(:, sortOrder);

unmixing = eigenVectorsSorted';
mixing = inv(unmixing); % orthogonal matrices: the transpose is the inverse

% Project the data onto the principal components
principalComponents =  unmixing'* data;

% Display the results
% disp('Eigenvalues:');
% disp(eigenValuesSorted);
% disp('Principal Components:');
% disp(principalComponents);


%% how much would each channel get from component 1 if we were to backproject?

cfg = [];
cfg.grad = all_data.grad;
lay = ft_prepare_layout(cfg);
% Prepare data for topoplot
cfg = [];
cfg.layout = lay; % Use the electrode layout from rawData.elecs
cfg.zlim = 'maxabs'; % Color limits
cfg.comment = 'no'; 

% Create a structure for the topoplot
topoData = [];
topoData.label = all_data.grad.label; 
topoData.avg = mixing(:,1); 
topoData.dimord = 'chan'; 
topoData.time = 0; 

% Plot the topomap for alpha power (only using the EEG channels)
figure;
ft_topoplotER(cfg, topoData);
colorbar; 
title('Component Topomap');
%%
% Extract the first principal component
component1 = principalComponents(1, :);

% Perform FFT
N2 = length(component1);
fftResult = fft(component1);

% Compute the frequency axis
Fs = all_data.fsample; % Sampling frequency in Hz (adjust as needed)
f = (0:N2-1)*(Fs/N2);

% Plot the magnitude of the FFT result
figure;
plot(f, abs(fftResult));
title('FFT of the First Principal Component');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([0 Fs/2]); % Plot up to the Nyquist frequency


figure;
plot(all_data.time{1}(1:100*Fs), component1(1:100*Fs));

cfg = [];
cfg.channel = all_data.label(1:157);
data = ft_selectdata(cfg, all_data);
cfg = [];
cfg.method = 'pca';
data_comp = ft_componentanalysis(cfg, data);

cfg = [];
cfg.layout = lay;
cfg.viewmode = 'component';
ft_databrowser(cfg, data_comp); 

%%
% high pass filter the data
cfg = [];
cfg.hpfilter = 'yes';
cfg.hpfreq = 1;
data = ft_preprocessing(cfg, data);


%%  run an ica (EEG lab algo runica)
cfg = [];
cfg.method = 'runica';
%cfg.unmixing = unmixing; % use precomputed unmixing matrix
data_comp = ft_componentanalysis(cfg, data);