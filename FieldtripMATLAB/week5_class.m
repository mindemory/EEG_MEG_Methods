%% Butterworth filter

% Sample time series data
fs = 1024; % Sampling frequency (Hz)

% Design a low-pass Butterworth IIR filter
order = 4; % Filter order
cutoff = 100; % Cutoff frequency (Hz)
[b, a] = butter(order, cutoff/(fs/2)); % Butterworth filter design


% Generate an impulse function

impulseX = [1, zeros(1, 199)];


% Filter the impulse function
impulse_response = filter(b, a, impulseX);


% Plot the impulse response
figure;
plot(impulse_response);
title('Impulse Response');
xlabel('Sample');
ylabel('Amplitude');

%% SHOW THAT IT DOESN't go to zero!

%% Compute the frequency response using DFT
N = 1024;
H = fft(impulse_response, N);

% how we know it: 
figure;
subplot(211);
plot((0:N-1), real(H).^2 + imag(H).^2);
xlim([0 512]);  % Nyquist
title('Frequency Response');
xlabel('frequency (hz)')
% Compute and plot the phase response
subplot(2,1,2);
plot((0:N-1),  angle(H));
xlim([0 512]);  % Nyquist
title('Phase Response');
xlabel('frequency (hz)')

% note: mention the discontinuity!

%% normalizations...

w = 2*pi*(0:N-1)./N; % Frequency vector w = 2*pi*f
decibel_power = 20*log10(abs(H)); % compute power in decibel
% Plot the magnitude response
figure;
subplot(2,1,1);
plot(w/pi, decibel_power);
title('Frequency Response');
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)');
xlim([0 1])

subplot(2,1,2);
plot(w/pi,  unwrap(angle(H)));
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Phase (radians)');
title('Phase Response');
xlim([0 1])

%% using freqz and phasez

[H2, w] = freqz(b, a, fs, 'whole');
figure;
subplot(211);
plot(w/pi,20*log10(abs(H2)));
title('Frequency Response');
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)');
xlim([0 1])

[PHI, w] = phasez(b, a, fs, 'whole');

subplot(212);
plot(w/pi,PHI);
title('Phase Response');
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Angle (radians)');
xlim([0 1])

%% filtfilt function

% Generate an impulse function
impulseX = [1, zeros(1, 199)];

% Filter the impulse function using filter
impulse_response = filter(b, a, impulseX);

% Filter the impulse function using filtfilt
impulse_response_filtfilt = filtfilt(b, a, impulseX);

% Plot the impulse responses
figure;
subplot(2,1,1);
plot(impulse_response);
title('Impulse Response using filter');
xlabel('Sample');
ylabel('Amplitude');

subplot(2,1,2);
plot(impulse_response_filtfilt);
title('Impulse Response using filtfilt');
xlabel('Sample');
ylabel('Amplitude');

% Show that the filtfilt response does not shift in time
figure;
plot(impulse_response, 'b', 'DisplayName', 'filter');
hold on;
plot(impulse_response_filtfilt, 'r', 'DisplayName', 'filtfilt');
title('Comparison of Impulse Responses');
xlabel('Sample');
ylabel('Amplitude');
legend;
hold off;

%% filter a signal

N = 1024; % Number of points for DFT
x = 0:1/fs:(N-1)/fs; % Time vector
signal = sin(2*pi*50*x) + sin(2*pi*120*x) + 2*rand(1,length(x)); % Sample signal with two frequencies (50 Hz and 120 Hz)

% Apply the IIR filter to the signal
y = filter(b, a, signal);

% Apply the filtfilt function to the signal
y_filtfilt = filtfilt(b, a, signal);

% Plot the original and filtered signals
figure;
subplot(3,1,1);
plot(x, signal);
title('Original Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(3,1,2);
plot(x, y);
title('Filtered Signal using filter');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(3,1,3);
plot(x, y_filtfilt);
title('Filtered Signal using filtfilt');
xlabel('Time (s)');
ylabel('Amplitude');
%%
% Plot the original and filtered signals on top of each other
figure;
plot(x, signal, 'k', 'DisplayName', 'Original Signal'); % Original signal in black
hold on;
plot(x, y, 'b', 'DisplayName', 'Filtered Signal (filter)'); % Filtered signal using filter in blue
plot(x, y_filtfilt, 'r', 'DisplayName', 'Filtered Signal (filtfilt)'); % Filtered signal using filtfilt in red
title('Comparison of Original and Filtered Signals');
xlabel('Time (s)');
ylabel('Amplitude');
legend;
hold off;
xlim([0 0.2]);

%% homework (start in class)
% create a signal and filter it with a  6th order low-pass butterworth filter

% characterize the impulse response of the filter

% filter EEG with a 6th order low-pass butterworth filter

% repeat this using fieldtrip's ft_preprocessing function/mne-python
% mne.filter.filter_data

%%
% read in meg 
% compute the spectrogram per channel with FFT
% screen continuous data for first quality check
% cut generously into trials
% mark moments of artifacts in trials and remaining data

% artifact definition should be in sample points of the recording

% have one dataset with trials as words
% keep track of the trialinfo (what word/odd or even/ modality/ category)

