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


meg_path = fullfile(dataPath, 'MEG');
eeg_path = fullfile(dataPath, 'EEG');

% Add fieldtrip to path
addpath(fullfile(fieldtrip_path));
ft_defaults;
%% list of subjects (could be cell array)

sjs = ...
    {'004'}; % todo change to make BIDS compliant

derivPath = fullfile(dataPath, 'derivatives', ...
    'preprocessing', ['sub-' sjs{1}]);
load(fullfile(derivPath, "data_pruned.mat"))

%% TF decomposition via wavelet


cfg            = [];
cfg.trials = find(...
    cellfun(@(x) ...
    strcmp('semanticVis', string(x.blockType)), ...
    data_pruned.trialinfo));
cfg.output     = 'pow';
cfg.channel    = 'all';
cfg.method     = 'wavelet';
cfg.toi        = -0.5 : 0.01 : 1;
cfg.foi        = 2:30;

% wavelet parameter (width in cycles)
cfg.width      = 5; % compare with 3 and 7??

tfr       = ft_freqanalysis(cfg, data_pruned);
%% inspect the power spectrum for a single channel:
figure;
imagesc(tfr.time, tfr.freq, squeeze(...
    tfr.powspctrm(20, :,:)));
axis xy
%% baseline correction
cfg = [];
cfg.baselinetype = 'relchange';
cfg.baseline     = [-0.5 0];
tfr_bl = ft_freqbaseline(cfg, tfr);

% save this as a dummy for later!
tfr_dummy = tfr_bl;

%% from the subfunction:
% 
% if (strcmp(baselinetype, 'absolute'))
%   data = data - meanVals;
% elseif (strcmp(baselinetype, 'relative'))
%   data = data ./ meanVals;
% elseif (strcmp(baselinetype, 'relchange'))
%   data = (data - meanVals) ./ meanVals;
% elseif (strcmp(baselinetype, 'normchange')) || (strcmp(baselinetype, 'vssum'))
%   data = (data - meanVals) ./ (data + meanVals);
% elseif (strcmp(baselinetype, 'db'))
%   data = 10*log10(data ./ meanVals);
% elseif (strcmp(baselinetype,'zscore'))
%     stdVals = repmat(nanstd(data(:,:,baselineTimes),1, 3), [1 1 size(data, 3)]);
%     data=(data-meanVals)./stdVals;
% else
%   ft_error('unsupported method for baseline normalization: %s', baselinetype);
% end

%%
%% plot the power spectrum
imagesc(tfr_bl.time, tfr_bl.freq, squeeze(...
    tfr_bl.powspctrm(20, :,:)));
axis xy

%% prepare the standard layout (we can do this from elec too but need to adjust scale and rotation)
cfg = [];
cfg.layout = 'BioSemi256';
cfg.feedback = 'yes';
lay = ft_prepare_layout(cfg);
%% plot the distrbution of power

figure;
cfg =[];
cfg.layout = lay;
ft_multiplotTFR(cfg, tfr_bl);

%% let's use hanning window

cfg            = [];
cfg.trials = find(...
    cellfun(@(x) ...
    strcmp('semanticVis', string(x.blockType)), ...
    data_pruned.trialinfo));
cfg.output     = 'pow';
cfg.toi        = -0.5 : 0.01 : 1;
cfg.foi        = 2:30;

cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
% Hanning parameters 
cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec

% alternative 
%cfg.t_ftimwin    = 7./cfg.foi;  % 7 cycles per time window

tfr       = ft_freqanalysis(cfg, data_pruned);

%% baseline correction
cfg = [];
cfg.baselinetype = 'relchange';
cfg.baseline     = [-0.5 0];
tfr_bl = ft_freqbaseline(cfg, tfr);
%% inspect the power spectrum for a single channel:
figure;
imagesc(tfr_bl.time, tfr_bl.freq, squeeze(...
    tfr_bl.powspctrm(20, :,:)));
axis xy

%% plot the distrbution of power

figure;
cfg =[];
cfg.layout = lay;
ft_multiplotTFR(cfg, tfr_bl);

%%%%%%%%%
% also show the different setting
%%%%%%%%%
% test different wavelet n-cycles
% mention frequency vs. time smoothing

%% Multitapers
% Multitapers are typically used in order to achieve 
% better control over the frequency smoothing. 
% More tapers for a given time window will result in 
% stronger smoothing. For frequencies above 30 Hz, 
% smoothing has been shown to be advantageous, 
% increasing sensitivity thanks to reduced variance in 
% the estimates despite reduced effective spectral 
% resolution. 
%%
cfg = [];
cfg.trials = find(...
    cellfun(@(x) ...
    strcmp('semanticVis', string(x.blockType)), ...
    data_pruned.trialinfo));
cfg.output     = 'pow';
cfg.toi        = -0.5 : 0.01 : 1;
cfg.foi        = 2:2:30;
cfg.method     = 'mtmconvol';
cfg.taper      =  'dpss'; % discrete prolate spheroidal sequences
cfg.t_ftimwin  = 5./cfg.foi;
cfg.tapsmofrq  = 0.4 *cfg.foi;
% number, the amount of spectral smoothing through
% multi-tapering. Note that 4 Hz smoothing means
% plus-minus 4 Hz, i.e. a 8 Hz smoothing box.
tfrmult = ft_freqanalysis(cfg, data_pruned);
%% baseline correction and plot
cfg = [];
cfg.baselinetype = 'relchange';
cfg.baseline     = [-0.5 0];
tfr_bl = ft_freqbaseline(cfg, tfrmult);
figure;
cfg =[];
cfg.layout = lay;
ft_multiplotTFR(cfg, tfr_bl);
%% SLIDES


%% now let's compute inter trial phase coherence:
cfg            = [];
cfg.trials = find(...
    cellfun(@(x) ...
    strcmp('semanticVis', string(x.blockType)), ...
    data_pruned.trialinfo));
cfg.output     = 'fourier'; % return fourier coefficients
cfg.channel    = 'all';
cfg.method     = 'wavelet';
cfg.width      = 5;
cfg.toi        = -0.5 : 0.01 : 1;
cfg.foi        = 2:30;
tfr       = ft_freqanalysis(cfg, data_pruned);

size(tfr.fourierspctrm);
fspcrtm_norm = tfr.fourierspctrm./abs(tfr.fourierspctrm);

tfr_dummy.powspctrm = squeeze(abs(mean(fspcrtm_norm)));
%% plot the inter-trial-phase-coherence

figure;
cfg =[];
cfg.layout = lay;
ft_multiplotTFR(cfg, tfr_dummy);


%% (homework?) demonstrate trial bias


%%(homework?) compute pairwise phase consistency on a subset of trials


