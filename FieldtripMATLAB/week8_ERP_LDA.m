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

%% processing for ERP
cfg = [];
cfg.detrend = 'no';
cfg.lpfilter = 'yes';
cfg.lpfreq = 15;
cfg.detrend = 'yes';
cfg.demean              = 'yes';    % we demean (baseline correct) ...
cfg.baselinewindow      = [-0.5 0]; % using the mean activity in this window
data_bl     = ft_preprocessing(cfg, data_pruned);
clear data_pruned
%% prepare the standard layout (we can do this from elec too but need to adjust scale and rotation)
cfg = [];
cfg.layout = 'BioSemi256';
cfg.feedback = 'yes';
lay = ft_prepare_layout(cfg);
%% compute ERPs via fieldtrip function

cfg = [];
cfg.latency            = [-0.3 1]; 
cfg.trials = find(...
    cellfun(@(x) ...
    strcmp('semanticVis', string(x.blockType)), ...
    data_bl.trialinfo));

tl_vis = ft_timelockanalysis(cfg, data_bl);
cfg.trials = find(...
    cellfun(@(x) ...
    strcmp('semanticAud', string(x.blockType)), ...
    data_bl.trialinfo));
tl_aud = ft_timelockanalysis(cfg, data_bl);

%% plot the ERPs separately

figure;
cfg =[];
cfg.layout = lay;
ft_multiplotER(cfg, tl_vis);
%%
cfg =[];
cfg.layout = lay;
ft_multiplotER(cfg, tl_aud);

%% trimmean is a robust measure that discards a percentage of extreme observations
% we have to compute this manually
tl = tl_vis;
cfg =[];
cfg.latency            = [-0.3 1]; 
cfg.keeptrials = 'yes';
cfg.trials = find(...
    cellfun(@(x) ...
    strcmp('semanticVis', string(x.blockType)), ...
    data_bl.trialinfo));
tl_kt = ft_timelockanalysis(cfg, data_bl);;
avg_new = trimmean(tl_kt.trial, 10);

%% overwrite and plot the new ERP
tl.avg      = squeeze(avg_new);
figure;
cfg =[];
cfg.layout = lay;
ft_multiplotER(cfg, tl);

%% compare the two

ft_multiplotER(cfg, tl, tl_vis);

%% prepare for LDA
cfg =[];
cfg.latency            = [-0.3 1]; 
cfg.keeptrials = 'yes';
cfg.trials = find(...
    cellfun(@(x) ...
    strcmp('semanticVis', string(x.blockType)), ...
    data_bl.trialinfo));
tl_vis = ft_timelockanalysis(cfg, data_bl);


cfg.trials = find(...
    cellfun(@(x) ...
    strcmp('semanticAud', string(x.blockType)), ...
    data_bl.trialinfo));
tl_aud = ft_timelockanalysis(cfg, data_bl);

% Combine data for LDA
data_combined = cat(1, tl_vis.trial, tl_aud.trial);
labels = [ones(size(tl_vis.trial, 1), 1); 2 * ones(size(tl_aud.trial, 1), 1)];


% Reshape data for LDA
[n_trials, n_channels, n_timepoints] = size(data_combined);

% Initialize accuracy array
accuracy_timecourse = zeros(n_timepoints, 1);


%% 10-Fold Cross-Validation
cv = cvpartition(labels, 'KFold', 10);

tic;
% Loop through each time point
for t = 1:n_timepoints
    t/n_timepoints
    data_t = data_combined(:, :, t);
    
    fold_accuracies = zeros(cv.NumTestSets, 1);
    % Loop over folds
    for i = 1:cv.NumTestSets
        
        train_idx = training(cv, i);
        test_idx = test(cv, i);
        
        train_data = data_t(train_idx, :);
        train_labels = labels(train_idx);
        test_data = data_t(test_idx, :);
        test_labels = labels(test_idx);
        
        % Perform LDA
        mean_vis = mean(train_data(train_labels == 1, :), 1);
        mean_aud = mean(train_data(train_labels == 2, :), 1);
        % within class covariance
        Sw = cov(train_data(train_labels == 1, :)) + cov(train_data(train_labels == 2, :));
        % between class covariance
        Sb = (mean_vis - mean_aud)' * (mean_vis - mean_aud);

        % get the eigenvectors
        [V, D] = eig(Sb, Sw);
        % and sort them
        [~, idx] = sort(diag(D), 'descend');
        V = V(:, idx);
        
        % Project test data onto first LDA component
        test_projection = test_data * V(:, 1);
        % ... explain multi-class LDA here

        % compute the mean of projected training data
        mean_vis_proj = mean(train_data(train_labels == 1, :) * V(:, 1));
        mean_aud_proj = mean(train_data(train_labels == 2, :) * V(:, 1));
        
        % test if projected data is closer to mean of class 1 or 2
        predictions = zeros(size(test_labels));
        for j = 1:length(test_labels)
            if abs(test_projection(j) - mean_vis_proj) < abs(test_projection(j) - mean_aud_proj)
                predictions(j) = 1;
            else
                predictions(j) = 2;
            end
        end
        
        % Calculate accuracy for this fold
        fold_accuracies(i) = mean(predictions == test_labels);
    end
    
    % Average accuracy across folds for this time point
    accuracy_timecourse(t) = mean(fold_accuracies);
end

%%
% Plot the time-course of classification accuracy
time = linspace(cfg.latency(1), cfg.latency(2), n_timepoints);
figure;
plot(time, accuracy_timecourse, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Classification Accuracy');
title('Time-Course of Classification Accuracy');
grid on;
%% LDA beamforming -----
% to extract sources of ERPs (Treder et al.)


%% get the reference data
cfg = [];
cfg.latency = [-0.5 1];
cfg.trials = find(...
    cellfun(@(x) ...
    strcmp('semanticVis', string(x.blockType)), ...
    data_bl.trialinfo));
data_reference = ft_selectdata(cfg, data_bl);

clear data_bl;
%% compute the ERP
cfg =[];
cfg.latency            = [-0.5 1]; % with a bit more padding 
cfg.keeptrials = 'no'; % AVERAGING
cfg.trials = find(...
    cellfun(@(x) ...
    strcmp('semanticVis', string(x.blockType)), ...
    data_reference.trialinfo));
tl_vis = ft_timelockanalysis(cfg, data_reference);
%% the interval of interest
ival = [0.12 0.15];

% To extract the spatial pattern and store it in a vector, we first find
% the corresponding time samples in the 0.45-0.75s window, and then 
% calculate the average spatial pattern across these time samples
time_idx = (tl_vis.time >= ival(1) & tl_vis.time<= ival(2));
P = mean( tl_vis.avg(:, time_idx),2);

%% Use the low-level plotting function for plotting arbitrary vectors
cfg_topo= {'mask' lay.mask, 'datmask' [], 'interplim' 'mask'};
cfg_lay = { 'point' 'no' 'box' 'no' 'label' 'no' };

% plot the pattern
ft_plot_topo(lay.pos(1:n_channels,1),...
    lay.pos(1:n_channels,2), P, cfg_topo{:});
 ft_plot_layout(lay, cfg_lay{:})
    title('Spatial pattern P')


%% compare with manual computation

% compute the Covariance matrix
C2 = zeros(n_channels);
n_trials = numel(data_reference.trial);
for ii=1:n_trials
    % center
    X = data_reference.trial{ii}...
        - mean(data_reference.trial{ii}, 2);
    C2 = C2 + cov(squeeze(X'));
end
%% regularize the covariance matrix
% Define shrinkage target as diagonal matrix with with entries equal to
% trace/(number sensors) (=mean variance) of the covariance matrix. The
% scaling factor nu puts it on equal footing (same total variance)
% with the empirical covariance matrix
nu = trace(C2)/n_channels;
I  = nu*eye(n_channels);
gamma = 0.0001;
% The regularized covariance matrix is a convex combination of the empirical
% covariance matrix and a diagonal covariance matrix with equal total variance.
C2 = C2*(1-gamma) + I*gamma;

C2 = C2 / n_trials;

%% Calculate spatial filter w

w2 = inv(C2)*P*inv(P'*inv(C2)*P);

%% apply spatial filter to the data:
data_lda = data_reference;
for tt = 1: n_trials
    data_lda.trial{tt} = w2'*data_reference.trial{tt};
end
data_lda.label = {'ERP_component'};
data_lda = rmfield(data_lda, 'elec');
data_lda = rmfield(data_lda, 'hdr');
data_lda = rmfield(data_lda, 'cfg');
data_lda = rmfield(data_lda, 'sampleinfo');

lda_tl = ft_timelockanalysis([], data_lda);
figure; 
plot(lda_tl.time, lda_tl.avg)
%% compare this to Treder et al. toolbox
[w,lda,C] = LDAbeamformer(P,data_reference);
%%
figure
subplot(1,3,1), 
    imagesc(C),title('Covariance matrix');
    clim([-10 10]);
    colorbar;
subplot(1,3,2), 
    ft_plot_topo(lay.pos(1:n_channels,1), lay.pos(1:n_channels,2), P, cfg_topo{:});
    ft_plot_layout(lay, cfg_lay{:})
    title('Spatial pattern P')
subplot(1,3,3), 
    ft_plot_topo(lay.pos(1:n_channels,1), lay.pos(1:n_channels,2), w, cfg_topo{:});
    ft_plot_layout(lay, cfg_lay{:})
    title('Spatial filter w');

 %% Source-level ERP
figure
avg_LDA = ft_timelockanalysis([], lda);

cfg= [];
cfg.xlim        = [-0.3 1.0];
cfg.channel     = 'B17';
ft_singleplotER(cfg,tl_vis);
hold on;
% Scale LDA source output to fit the channel plot
yl= ylim;
hold all
plot(avg_LDA.time, avg_LDA.avg./max(avg_LDA.avg(:)) .* yl(2), 'r')
ylim(1.1.*ylim);
set(gca,'ylim',[yl(1) max(yl(2), max(avg_LDA.avg * yl(2)) )])

legend({cfg.channel, avg_LDA.label{1}})
title('Event-related potential (best electrode vs. spatial filter)');

%% we can use this to analyze single trial responses
plot(data_lda.time{1}, cell2mat(data_best.trial(1:50)')', 'k')

%% select best electrode data
cfg = [];
cfg.channel = 'B17';
data_best = ft_selectdata(cfg,data_reference);
%% plot them together
subplot(211);
plot(data_best.time{1}, cell2mat(data_best.trial(1:end)')', 'k')
title('best electrode single trials')
subplot(212);
plot(data_lda.time{1}, cell2mat(data_lda.trial(1:end)')', 'k')
title('spatially filtered single trials (noise attenuated via LDA)')



