clear; close all; clc;

[ret, hostname] = system('hostname');
if ret ~= 0
    hostnme = getenv('HOSTNAME');
end
hostname = strtrim(hostname);

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
%% TODO: this needs to be BIDS


data_path = fullfile(meg_path, 'GroupB', 'Recording');

cfg = [];
cfg.dataset =fullfile(data_path, 'R2470_EventModel_9.26.24_NR.sqd');
all_data = ft_preprocessing(cfg);

%% inspect the data
cfg = [];
cfg.channel = all_data.label(1:157);
cfg.blocksize = 10;
cfg.preproc.hpfilter = 'yes';
cfg.preproc.hpfreq = 1;
ft_databrowser(cfg, all_data);
%%  reject bad channels here
cfg = [];
cfg.channel = all_data.label([1:14,16:157]);
all_data = ft_selectdata(cfg, all_data);



% load in artifactual samples here and use them

% sampling info in raw file corresponding to data
% si_vector = [eeg.sampleinfo(1,1):eeg.sampleinfo(1,2), ...
%     eeg.sampleinfo(2,1):eeg.sampleinfo(2,2)];
% 
% 
% dummy = [1:10, 158:190 ];

%% inspect the data
cfg = [];
cfg.blocksize = 10;
cfg.preproc.hpfilter = 'yes';
cfg.preproc.hpfreq = 0.5;
ft_databrowser(cfg, all_data);
% cfg = [];
% cfg.channel = all_data.label([1:14, 16:157]);
% all_data = ft_selectdata(cfg, all_data);
%% do a pca manually


% resample only for speed (usually use 512 Hz)
% note that artifact annotation will be relative to sampleinfo. This
% changes with re-sampling
% cfg = [];
% cfg.resamplefs = 256;
% all_data = ft_resampledata(cfg, all_data);

% let's not use artifacts in this stepe
%data(:,artifact_vector) = [];
data = all_data.trial{1};
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
topoData.label = all_data.grad.label([1:14,16:end]); 
topoData.avg = mixing(:,1); 
topoData.dimord = 'chan'; 
topoData.time = 0; 

% Plot the topomap for alpha power (only using the EEG channels)
figure;
ft_topoplotER(cfg, topoData);
colorbar; 
title('Component Topomap');
%% plot the component 
% Extract the first principal component
component1 = principalComponents(1, :);

% plot the component
figure;
plot(all_data.time{1}(1:100*Fs), component1(1:100*Fs));

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


%% Use fieldtrip to compute the pca

cfg = [];
cfg.method = 'pca';
data_comp = ft_componentanalysis(cfg, all_data);

cfg = [];
cfg.layout = lay;
cfg.viewmode = 'component';
ft_databrowser(cfg, data_comp); 

% high pass filter the data
data = all_data;
cfg = [];
cfg.hpfilter = 'yes';
cfg.hpfreq = 1;
data = ft_preprocessing(cfg, data);

%% here we could overwrite data with nan for ICA computation
% Load in the cfg from artifact annotation 
% TODO...
artifact_vector = zeros(1, size(data.trial{1},2));
for aa = 1: size(cfg_artifacts.artfctdef.visual.artifact,1)
    
    artifact_vector(cfg_artifacts.artfctdef.visual.artifact(aa,1)...
        : cfg_artifacts.artfctdef.visual.artifact(aa,2)) = 1;
end
data.trial{1}(:, artifact_vec) = nan;


%%  run an ica (EEG lab algo runica)
cfg = [];
cfg.method = 'runica';
data_comp = ft_componentanalysis(cfg, data);

% save your unmixing matrix here

% decompose the original data as it was prior to downsampling, filtering, 
% and artifact masking

cfg = [];
cfg.method = 'runica';
% load pre-computed solution
% TODO
%cfg.unmixing = unmixing; % use precomputed unmixing matrix
%cfg.topolabel = all_data.label
data_comp = ft_componentanalysis(cfg, all_data);

%% now plot FFT of the component, and use databrowser as componentbrowser
% this only works on local machine
% mainfig = figure(); 
% tabgroup = uitabgroup(mainfig, 'Position', [0 .05 1 .95]);
% cnt = 0;
% ax_ = nan(size(data_comp.label)); %init handles to the fft axes
% minsize_ = min(cellfun(@(x) size(x,2), data_comp.trial));
% for step_ =1 : 10: 10% size(data_comp.label,1); %loop through components in steps of 10...
%     cnt = cnt+1; vec = [step_:step_+9];   
%     tab(cnt)=uitab(tabgroup,'Title', ['Comps_', num2str(vec(1)) '-' num2str(vec(end))]);
%     axes('parent',tab(cnt))    
%     ind_ = 0;
%     for comp = vec; % loop through the subcomponents in order to fft and draw
%         if comp> size(data_comp.trial{1,1},1); continue; end %do not plot if theres no comp
%         ind_ = ind_+1;
%         % ----- here we calculate the fft to get the powerspectrum
%         pspctrm = 0; cnt = 0;
%         for tr = 1 : numel(data_comp.trial)
% %             if strcmp(data_comp.trialinfo{tr}.type, 'retrieval'); continue; end;
%             signal  = data_comp.trial{1,tr}(comp,1:minsize_);
%             N       = length(signal);
%             nyquist = data_comp.fsample/2;
%             fourierCoefs = zeros(size(signal));
%             frequencies = linspace(0,nyquist,floor(N/2)+1);
%             fourierCoefsF = fft(signal) / N;
%             pspctrm = pspctrm + abs(fourierCoefsF(1:length(frequencies)))*2;
%             cnt = cnt +1;
%         end
%         pspctrm = pspctrm/cnt;
%         % ----------------------------------------------------------------
% 
%         subplot(5,4,2*(ind_-1)+1),hold on %open a new subplot for the topo
% 
%         cfg           = [];
%         cfg.component = [comp:comp];       % specify the component(s) that should be plotted
%         cfg.layout = lay; % specify the layout file that should be used for plotting
%         cfg.comment   = 'no';
%         ft_topoplotIC(cfg, data_comp)
%         ax_(ind_+step_-1)=subplot(5,4,2*ind_);
%         set(ax_(ind_),'xlim', [0 60]);
%         hold on;
%         plot(frequencies,pspctrm)
% 
%     end
%     ax_ = ax_(~isnan(ax_));
%     hax = @(src, ax_) set(ax_, 'xlim', [0 get( src, 'Value' )] );
%     uicontrol('Style', 'slider', 'Units','normalized', 'Min',20,'Max',300,'Value',60, 'Position', [0 0 1 0.05], 'Callback', @(src,evt) hax( src, ax_ ) );
% end

cfg = [];
cfg.layout =  lay;
cfg.channel = [1:10]; % components to be plotted
cfg.viewmode = 'component';
ft_databrowser(cfg, data_comp)

%% reject the components
cfg           = [];
cfg.component = [1 2 15 20];
data_clean    = ft_rejectcomponent(cfg,data_comp, data);
% the second input is important to adjust data.grad for future source
% reconstruction

% also store which components were rejected for this subject
% component indices + unmixing matrix can reproduce the preprocessing
% pipeline quickly
