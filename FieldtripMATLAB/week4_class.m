clear;
clc;

%%  note: mention t (from 0 to T) vs. x (from 0 to 1)
%% -- let's continue with ALL frequencies and use Matrix multiplication
% Create the Signal again:
% Define the parameters
SR = 512;             % Sampling frequency
dT = 1/SR;            % Sampling period
N = 2048;             % Length of signal 
n = 0:N-1;            % Sample indices
t = n*dT;             % Time vector
x = (0:N-1) ./ N;      % Sample vector runs from 0 to 1
% Create a simple signal (sum of sine waves)
X = 2*sin(2 * pi * 2 * n / SR) + sin(2 * pi * 10 * n / SR) + ...
    3*cos(2 * pi * 6 * n / SR + pi) ... 
    + cos(2 * pi * 10 * n / SR + pi) + rand(1, N); % cosine has a phase offset of pi

% Plot the signal
figure;
plot(t, X);
title('Original Signal');
xlabel('t (seconds)');
ylabel('X(t)');

%% Here we compute the DFT Manually:
% Compute the DFT Matrix
% initialize a zero-matrix 
W_dft = zeros(N,N);
% fill the matrix by looping over frequencies
for k = 0:N-1
    % note that using the sample vector x instead of t defines frequencies
    % as repetitions in the whole segment
    W_dft(k+1,:) =  exp(-1i * 2 * pi * k * x); 
end
% here we compute all the inner products in one go (each row is a frequency)
X_dft = W_dft * X';


% now we compute the power of our complex Fourier coefficients
power = real(X_dft).^2 + imag(X_dft).^2; % NOTE: this is elementwise .^

% for the correct frequency axis, we need to rescale to 1/second (Hertz)
frequencies = [0:N-1] .* (SR / N);

figure
% Plot the barplot of power versus frequency (we only plot up to Nyquist
% frequency. Everything above is mirrored.)
bar(frequencies(1:N/2), power(1:N/2))
xlabel('Frequency (Hz)')
ylabel('Power')
title('Power vs Frequency');
xlim([0 20])
% ===> Note: Explain the zero-frequency
%% Next, we compute the inverse DFT Manually:

% initialize matrix
W_idft = zeros(N,N);
% Fill (now we will sum across frequencies to get time). Note that this is
% almost the same formula, only i becomes positive now (i.e., this is the
% conjugate)
for n = 0:N-1
    W_idft(n+1,:) =  exp(1i * 2 * pi * (0:N-1) *  n/ N);
end

X_idft = W_idft*X_dft;

% we note that W_idft is simply the conjugate of W_dft so instead, we 
% multiply with the conjugate from the left (summing over frequencies)

X_idft = (conj(W_dft)*X_dft)/N;
%%
%Compare the Original and Reconstructed Signals:
% Plot both signals for comparison. We can now shift between representing
% our signal in the time domain and in the frequency domain. 
figure;
subplot(3,1,1);
plot(t, X, 'r');
title('Original Signal');
xlabel('t (seconds)');
ylabel('X(t)');

subplot(3,1,2);
plot(t, X_idft, 'b');
title('Reconstructed Signal');
xlabel('t (seconds)');
ylabel('X\_reconstructed(t)');

subplot(3,1,3);
plot(t, X, 'r');
hold on;
plot(t, X_idft, 'b--');
title('Both Signals');
xlabel('t (seconds)');
xlabel('t (seconds)');

%% HOMEWORK FROM LAST WEEK
%% compute the Inverse DFT, but set the low frequencies to 0:
W_idft = conj(W_dft);
X_dft_filt = X_dft;
X_dft_filt(1:36) = 0;
X_dft_filt(N-36:N) = 0; % set the negative frequencies to 0 as well
% 
% alternatively, we could set our W_inv to zero and not consider those 
% frequencies in the reconstruction 
%-------
% W_idft(:, 1:36) = 0 + 1i*0;
% % the high frequencies are mirrored
% W_idft(:, N-36:N) = 0 + 1i*0;
%-------

X_idft_filtered = (W_idft*X_dft_filt)/N;

%% Compare the Original and Reconstructed Signals:
% Plot both signals for comparison
figure;
subplot(3,1,1);
plot(t, X, 'r');
title('Original Signal');
xlabel('t (seconds)');
ylabel('X(t)');

subplot(3,1,2);
plot(t, X_idft_filtered, 'b');
title('Filtered Signal');
xlabel('t (seconds)');
ylabel('X\_reconstructed(t)');


subplot(3,1,3);
plot(t, X, 'r');
hold on;
plot(t, X_idft_filtered, 'b--');
title('Both Signals');
xlabel('t (seconds)');

%%  ==> we have now high-pass filtered our data.
% We will soon see that this is not the ideal way to filter our data

%% how does non-stationarity affect the DFT?
SR = 512;             % Sampling frequency
dT = 1/SR;            % Sampling period
N = 2048;             % Length of signal 
n = 0:N-1;            % Sample indices
t = n*dT;             % Time vector
x = (0:N-1) ./ N;      % Sample vector runs from 0 to 1



% make an oscillation with changing amplitude
f_ns = 6; % Frequency of the non-stationary signal
decay = exp(-n/(N/4)); % Exponential decay function
X_6 = 4 * 3*cos(2 * pi * f_ns * n / SR);
ns_X = (decay .* X_6);
% Plot the signal
figure;
plot(t, ns_X);
title('Non-stationary Signal');
xlabel('t (seconds)');
ylabel('X(t)');

% For comparison, we create the same signal as before
X = 2*sin(2 * pi * 2 * n / SR) + sin(2 * pi * 10 * n / SR) + ...
    3*cos(2 * pi * 6 * n / SR + pi) ... 
    + cos(2 * pi * 10 * n / SR + pi) + rand(1, N); % cosine has a phase offset of pi


% here we change our 6Hz cosine for the non-stationary data
X_ns = 2*sin(2 * pi * 2 * n / SR) + sin(2 * pi * 10 * n / SR) + ...
    ns_X ... 
    + cos(2 * pi * 10 * n / SR + pi) + rand(1, N); % cosine has a phase offset of pi

%%
% Plot the signals
figure;
subplot(211)
plot(t, X);
title('Previous Signal');
xlabel('t (seconds)');
ylabel('X(t)');

subplot(212)
plot(t, X_ns);
title('Non-stationary Signal');
xlabel('t (seconds)');
ylabel('X(t)');
%% Here we compute the DFT Manually:
% Compute the DFT Matrix
% initialize a zero-matrix 
W_dft = zeros(N,N);
% fill the matrix by looping over frequencies
for k = 0:N-1
    % note that using the sample vector x instead of t defines frequencies
    % as repetitions in the whole segment
    W_dft(k+1,:) =  exp(-1i * 2 * pi * k * x); 
end
% here we compute all the inner products in one go (each row is a frequency)
X_dft = W_dft * X';

% here we compute the power of the non-stationary signal
X_dft_ns = W_dft * X_ns';

% now we compute the power of our complex Fourier coefficients
power = real(X_dft).^2 + imag(X_dft).^2; % NOTE: this is elementwise .^

power_ns = real(X_dft_ns).^2 + imag(X_dft_ns).^2; % NOTE: this is elementwise .^

% for the correct frequency axis, we need to rescale to 1/second (Hertz)
frequencies = [0:N-1] .* (SR / N);

figure;
subplot(211);
% Plot the barplot of power versus frequency (we only plot up to Nyquist
% frequency. Everything above is mirrored.)
bar(frequencies(1:N/2), power(1:N/2))
xlabel('Frequency (Hz)')
ylabel('Power')
title('Power vs Frequency');
xlim([0 20])

subplot(212);
% Plot the barplot of power versus frequency (we only plot up to Nyquist
% frequency. Everything above is mirrored.)
bar(frequencies(1:N/2), power_ns(1:N/2))
xlabel('Frequency (Hz)')
ylabel('Power')
title('Power vs Frequency');
xlim([0 20])

%% SLIDES

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%       CONVOLUTION              %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% EDGE DETECTION (important for finding triggers)
% Define the signal and the kernel
signal = [zeros(1,512), ones(1,512), zeros(1,512)];
kernel = [1, 0, -1]; % Example of a simple edge-detection kernel

% Perform convolution
convolved_signal = conv(signal, kernel, 'same');

% Plot the original signal and the convolved signal
figure;
subplot(2,1,1);
plot(signal, 'r');
title('Original Signal');
xlabel('Sample Index');
ylabel('Amplitude');

subplot(2,1,2);
plot(convolved_signal, 'b');
title('Convolved Signal');
xlabel('Sample Index');
ylabel('Amplitude');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%       FILTER KERNELS           %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Play with different filter kernels
low_pass_kernel = ones(1, 5) / 5; % Simple moving average filter
high_pass_kernel = [-1, -1, 4, -1, -1]; % Example of a simple high-pass filter

% Create a Gaussian filter manually
sigma = 1; % Standard deviation
kernel_size = 51; % Size of the kernel
t = linspace(-2, 2, kernel_size);
gaussian_kernel = exp(-t.^2 / (2 * sigma^2));
gaussian_kernel = gaussian_kernel / sum(gaussian_kernel); % Normalize

% Apply the low-pass filter
filtered_signal_low = conv(signal, low_pass_kernel, 'same');

% Apply the high-pass filter
filtered_signal_high = conv(signal, high_pass_kernel, 'same');

% Apply the Gaussian filter
filtered_signal_gaussian = conv(signal, gaussian_kernel, 'same');

% Plot the original and filtered signals
figure;
subplot(5,1,1);
plot(signal, 'r');
title('Original Signal');
xlabel('Sample Index');
ylabel('Amplitude');
ylim([-1 2]);
xlim([1 length(signal)]);

subplot(5,1,2);
plot(filtered_signal_low, 'b');
title('Low-Pass Filtered Signal');
xlabel('Sample Index');
ylabel('Amplitude');
ylim([-1 2])
xlim([1 length(signal)])
subplot(5,1,3);
plot(filtered_signal_high, 'g');
title('High-Pass Filtered Signal');
xlabel('Sample Index');
ylabel('Amplitude');
ylim([-1 2])
xlim([1 length(signal)])

subplot(5,1,4);

plot(gaussian_kernel, 'm');
title('Gaussian Kernel');
xlabel('Sample Index');
ylabel('Amplitude');
xlim([-100 length(gaussian_kernel)+100])

subplot(5,1,5);
plot(filtered_signal_gaussian, 'm');
title('Gaussian Filtered Signal');
xlabel('Sample Index');
ylabel('Amplitude');
ylim([-1 2])
xlim([1 length(signal)])
%% convolution in time is multiplication in frequency:

conv_length = length(signal) + length(gaussian_kernel) -1;
x=fft(signal,conv_length);
y=fft(gaussian_kernel,conv_length);
y_ifft = ifft(x.*y);
y_conv = conv(signal,gaussian_kernel, 'same');


subplot(5,1,1);
plot(signal, 'r');
title('Original Signal');
xlabel('Sample Index');
ylabel('Amplitude');
ylim([-1 2]);
xlim([1 length(signal)]);

subplot(5,1,2);
plot(gaussian_kernel, 'm');
title('Gaussian Kernel');
xlabel('Sample Index');
ylabel('Amplitude');
xlim([-100 length(gaussian_kernel)+100])

subplot(5,1,3);
plot(y_conv, 'r');
title('Gaussian Filtered Signal (convolution)');
xlabel('Sample Index');
ylabel('Amplitude');
ylim([-1 2])
xlim([1 length(signal)]);

subplot(5,1,4);
plot(y_ifft(...
    ceil(length(...
    gaussian_kernel)/2)+1 ... 
    :length(y_conv)+floor(length(gaussian_kernel)/2)), 'b');
title('Gaussian Filtered Signal (multiplication in frequency)');
xlabel('Sample Index');
ylabel('Amplitude');
ylim([-1 2])
xlim([1 length(signal)]);

subplot(5,1,5);
plot(y_conv, 'r'); hold on;
title('Gaussian Filtered Signal (convolution)');
plot(y_ifft(...
    ceil(length(...
    gaussian_kernel)/2)+1 ... 
    :length(y_conv)+floor(length(gaussian_kernel)/2)), 'b--');
title('Comparison');
xlabel('Sample Index');
ylabel('Amplitude');
ylim([-1 2])
xlim([1 length(signal)])


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%       MORLET WAVELETS          %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the parameters for the Morlet wavelet
f0 = 6; % Central frequency of the wavelet
SR = 100;
sigma = 1; % Standard deviation of the Gaussian window
t = -5:1/SR:5; % Time vector
N_ml = length(t);

% Create the Morlet wavelet (%frequency based on time)
A = 1 / (sqrt(sigma * sqrt(pi)) * sqrt(f0)); % amplitude scaling for the frequency-specific Morlet wavelet
morlet_wavelet = exp(2 * 1i * pi * f0 * t) .* exp(-t.^2 / (2 * sigma^2));

figure;
subplot(3,1,1);
title('Cosine (real part of complex exponential)');
xlabel('Sample Index');
ylabel('Amplitude');
plot( t,real(exp(2 * 1i * pi * f0 * t)))
xlim([t(1) t(end)]);
subplot(3,1,2);
title('Gaussian');
xlabel('Sample Index');
ylabel('Amplitude');
plot( t,exp(-t.^2 / (2 * sigma^2)));
xlim([t(1) t(end)]);
subplot(3,1,3);
title('Wavelet (real part)');
xlabel('Sample Index');
ylabel('Amplitude');
plot( t,real(morlet_wavelet));
xlim([t(1) t(end)]);

%%
% Plot the Morlet wavelet and the transformed signal
figure;
subplot(3,1,1);

plot(signal, 'r');
title('Original Signal');
xlabel('Sample Index');
ylabel('Amplitude');
ylim([-1 2]);
xlim([1 length(signal)]);


% Apply the Morlet wavelet to the signal (convolution)
wavelet_transformed_signal = conv(signal, morlet_wavelet, 'same');

subplot(3,1,2);
plot3(t, real(morlet_wavelet), imag(morlet_wavelet), 'r');
title('Morlet Wavelet');
xlabel('Time (seconds)');
ylabel('Amplitude');
legend('Real Part', 'Imaginary Part');

subplot(3,1,3);
plot3([1:length(signal)], real(wavelet_transformed_signal), ...
    imag(wavelet_transformed_signal), 'b');

title('Wavelet Transformed Signal');
xlabel('Sample Index');
ylabel('Amplitude');
%% NOW WE WANT TO CONVOLVE an impulse function with OUR WAVELET 

impulseX = zeros([1 5000]);

impulseX(2501) = 1;

y_conv = conv(impulseX, morlet_wavelet, 'same');

figure;
subplot(211);
plot(real(morlet_wavelet));
subplot(212);
plot(real(y_conv));

%% create X again
N = 5000;
n = 0:N-1;            % Sample indices

f_ns = 6; % Frequency of the non-stationary signal 

% narrow gaussian window shifted by 10 seconds 
gaussian_window = exp(-0.5 * (([-(N-1)/2-10*SR:(N-1)/2-10*SR] ./ (N/12)).^2));
%normalized to Area of 1
gaussian_window_norm = gaussian_window./sum(gaussian_window);
% modulated 6 Hz oscillaiton
X_6 = 800 * cos(2 * pi * f_ns * n / SR);
ns_X = (gaussian_window_norm .* X_6);
% x = (0:N-1) ./ N;     
X = 1*sin(2 * pi * 2 * n / SR) + sin(2 * pi * 10 * n / SR) + ...
    ns_X ... 
    + cos(2 * pi * 10 * n / SR + pi) + rand(1, N); % cosine has a phase offset of pi
figure;
subplot(211);
title('non-stationary part')
plot(ns_X);
subplot(212);
plot(X);
title('combined signal')

xlim([0 N]);
%% what is the frequency spectrum of our wavelet?

% get the correct frequency axis
frequencies = SR*(0:(length(morlet_wavelet)-1)/2)/length(morlet_wavelet); %[0:N-1] .* (SR / N);
morlet_fft = fft(morlet_wavelet);

% Compute the power spectrum
power = real(morlet_fft).^2 + imag(morlet_fft).^2; % 


figure;
% Plot the barplot of power versus frequency (we only plot up to Nyquist
% frequency. Everything above is mirrored.)
bar(frequencies, power(1:(length(morlet_wavelet)-1)/2+1))
xlabel('Frequency (Hz)')
ylabel('Power')
title('Power spectrum of Morlet Wavelet');
xlim([0 20]);
%% what is the frequency spectrum of our data?

% get the correct frequency axis
frequencies = SR*(0:(N-1)/2)/N; %[0:N-1] .* (SR / N);
X_fft = fft(X);

% Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1.
powerX = real(X_fft).^2 + imag(X_fft).^2; % 


figure;
% Plot the barplot of power versus frequency (we only plot up to Nyquist
% frequency. Everything above is mirrored.)
bar(frequencies, powerX(1:(N/2)))
xlabel('Frequency (Hz)')
ylabel('Power')
title('Power spectrum of original data');
xlim([0 15]);
%% now filter

conv_length = length(X) + length(morlet_wavelet) -1;

X_morlet_filtered = ifft(fft(morlet_wavelet, conv_length).*fft(X,conv_length)); % each fft is padded with zeros
X_morlet_convolved = conv(X, morlet_wavelet, 'same'); % left fft is padded with zeros
%%
figure;
subplot(411);
plot(n / SR,X,'r');
xlabel('time');
title('combined signal')
xlim([0 N/SR])
subplot(412);
plot(n / SR,real(X_morlet_filtered(501:end-500)),'r');
xlim([0 N/SR])
xlabel('time');
title('Morlet filtered (FFT)')
subplot(413);
plot(n / SR,real(X_morlet_convolved),'b');
xlim([0 N/SR])
xlabel('time');
title('Morlet filtered (Conv)')
subplot(414);
plot(n / SR,real(X_morlet_filtered(501:end-500)),'r'); hold on;
plot(n / SR,real(X_morlet_convolved),'b--');
xlim([0 N/SR])
xlabel('time');
title('Both')
%% what is the frequency spectrum of our filtered data?

% get the correct frequency axis
frequencies = SR*(0:(N-1)/2)/N; %[0:N-1] .* (SR / N);
X_fft = fft(X_morlet_convolved);

% Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1.
powerX = real(X_fft).^2 + imag(X_fft).^2; % 


figure;
% Plot the barplot of power versus frequency (we only plot up to Nyquist
% frequency. Everything above is mirrored.)
bar(frequencies, powerX(1:(N/2)))
xlabel('Frequency (Hz)')
ylabel('Power')
title('Power Spectrum of Wavelet-filtered data');
xlim([0 15]);

%% Practice (HOMEWORK) 
% create a Morlet wavelet at 10Hz and at 15Hz.
% convolve your EEG data with a Morlet wavelet
% add a spike to your EEG data by setting a random value to the value
% 10,000
% what happens after you convolve your EEG with the morlet wavelet?