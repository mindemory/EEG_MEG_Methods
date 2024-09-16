clear;
clc;
%% SLIDES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%             THE DOT PRODUCT                     %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define two vectors
A = [1; 2; 3];
B = [4; 5; 6];
%% Plot the vectors
figure;
quiver3(0, 0, 0, A(1), A(2), A(3), 'r', 'LineWidth', 2);
hold on;
quiver3(0, 0, 0, B(1), B(2), B(3), 'b', 'LineWidth', 2);
grid on;
xlabel('X');
ylabel('Y');
zlabel('Z');
legend('Vector A', 'Vector B');
axis equal;
hold off;
%% Calculate with a loop 
dot_product = 0;
for idx = 1 : size(A,1)
    dot_product = dot_product + A(idx)*B(idx);
end

fprintf('Dot product: %f\n', dot_product);
%% Calculate as sum of elementwise produt
dot_product = sum(A .* B);
fprintf('Dot product: %f\n', dot_product);

%% Calculate via vector multiplication
dot_product = A' * B;
fprintf('Dot product: %f\n', dot_product);

%% Calculate the dot product using the built-in function
dot_product = dot(A, B);
fprintf('Dot product: %f\n', dot_product);

%% SLIDES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%             MATRIX MULTIPLICATION               %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define two matrices
A = [1, 2; 2, 3; 3, 4];
B = [4, 5; 5, 6; 6, 7];

%% plot the matrices
figure;
subplot(121);
imagesc(A);  clim([0 8]);
title('3X2 Matrix')
axis off
subplot(122);
imagesc(B);  clim([0 8]);
title('3X2 Matrix');
axis off

%% compute a matrix multiplication
C = A' *B;
fprintf('Pairwise dot products: %f, %f, %f, %f\n', C(1, 1), C(1, 2), ...
    C(2, 1), C(2, 2));

%% plot the matrices
figure;
subplot(131);
imagesc(A');  clim([0 8]);
title('2X3 Matrix A^T');colorbar;
axis equal;
axis off
pbaspect([3 2 1]); % Set aspect ratio to 3:2
subplot(132);
imagesc(B);  clim([0 8]);
title('3X2 Matrix B');colorbar;
axis equal;
axis off
pbaspect([2 3 1]); % Set aspect ratio to 3:2
subplot(133);
imagesc(C);  clim([0 60]);
title('2X2 Matrix C'); colorbar;
axis equal;
axis off
pbaspect([1 1 1]); 
%% Exercises
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%   Exercises with EEG (Mrugank) %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - select an electrode via vector multiplication
% - elect a subset of eletrodes via matrix multiplication
% - average an electrode-ROI via vector multiplication
% - do a channel interpolation via matrix multiplication
% - channel 1 and 3 have been plugged into the wrong position in the cap. 
%   Reorder them in the recording via matrix multiplication.
% - re-reference your data to a different electrode via matrix multiplication
% - compute an average reference via matrix multiplication
%% SLIDES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%   From SPATIAL to TEMPORAL FILTERS    %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sample time series data
t = linspace(0, 10, 100); % Time vector
y = 5 + 2*t + 0.5*t.^2 + randn(size(t)); % Sample data with a polynomial trend

% Plot original data
figure;
subplot(6,1,1);
plot(t, y);
title('Original Data');
xlabel('Time');
ylabel('Value');

% Remove zero-order trend (mean)
p0 = polyfit(t, y, 0); % Fit a zero-order polynomial (mean)
y_zero_order_trend = polyval(p0, t); % Evaluate the polynomial
y_detrended_zero_order = y - y_zero_order_trend; % Subtract the zero-order trend
subplot(6,1,2);
plot(t, y_detrended_zero_order);
title('Zero-Order Detrended Data (Mean Removed)');
xlabel('Time');
ylabel('Value');

% Remove first-order trend (linear)
p1 = polyfit(t, y, 1); % Fit a first-order polynomial (linear)
y_first_order_trend = polyval(p1, t); % Evaluate the polynomial
y_detrended_first_order = y - y_first_order_trend; % Subtract the first-order trend

subplot(6,1,3);
plot(t, y_first_order_trend);
title('Linear Trend');
xlabel('Time');
ylabel('Value');

subplot(6,1,4);
plot(t, y_detrended_first_order);
title('First-Order Detrended Data (Linear Trend Removed)');
xlabel('Time');
ylabel('Value');

% Remove polynomial trend (e.g., second-order)
p2 = polyfit(t, y, 2); % Fit a second-order polynomial
y_poly_trend = polyval(p2, t); % Evaluate the polynomial


subplot(6,1,5);
plot(t, y_poly_trend);
title('Second order Trend');
xlabel('Time');
ylabel('Value');

y_detrended_poly = y - y_poly_trend; % Subtract the polynomial trend
subplot(6,1,6);
plot(t, y_detrended_poly);
title('Polynomial Detrended Data (Second-Order Trend Removed)');
xlabel('Time');
ylabel('Value');

%% EXERCISE (Mrugank) 
% detrend EEG channels and/or artificial channels



%% SLIDES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%   PLAYING WITH SINE WAVES      %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_cycles = 5;
SR = 512;
% Define the range of x (cycle 5 times)
t = linspace(0, n_cycles*2*pi,n_cycles*SR);

% Generate the sine and cosine waves
y_cos = cos(t);
y_sin = sin(t);

%% plot them separately and together
figure; 
subplot(311);
plot(t, y_cos, 'r'); xlim([t(1) t(end)])
title('cosine')
subplot(312);
plot(t, y_sin, 'b'); xlim([t(1) t(end)])
title('sine')
subplot(313);
plot(t, y_cos, 'r');
hold on;
plot(t, y_sin, 'b');xlim([t(1) t(end)])
title('cosine and sine')

%% We want to visualize COSINE and SINE together at the same time.
% To this end, we give sine its own AXIS

% make a polar plot of sine and cosine on the UNIT CIRCLE
% Define the phase
theta = pi/8;

%Calculate cosine and sine values
cos_val = cos(theta);
sin_val = sin(theta);

% Create the polar plot
figure;
polarplot([0 theta], [0 1], 'k'); % Unit circle
hold on;

% Highlight the cosine and sine on the unit circle
polarplot([0 0], [0 cos_val], 'r--'); % Cosine vector
polarplot([0 sin_val], [cos_val 1], 'b--'); % Sine vector

% Add labels and title
text(theta/2, cos_val/2, 'cos(\pi/8)', 'Color', 'r', 'FontSize', 12);
text(pi/16, 1, 'sin(\pi/8)', 'Color', 'b', 'FontSize', 12);
title('Polar Plot of Cosine and Sine at \pi/8');
hold off;
grid off;
pax = gca; % Get current axes
pax.ThetaAxisUnits = 'radians'; % Set theta axis to radians
%% Let's visualize COSINE and SINE together over time.
% Sine still gets its own axix; over time this becomes a 3D plot
close all;
% Create the figure
figure;
% Subplot for the 3D plot
subplot(2, 1, 1);
h3d = plot3(t, y_cos, y_sin, 'k', 'LineWidth', 2);
hold on;
dot1 = plot3(t(1), y_cos(1), y_sin(1), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
dot2 = plot3(t(1), y_cos(1), y_sin(1), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
xlabel('x');
ylabel('Real Part (cos(x))');
zlabel('Imaginary Part (sin(x))');
title('Complex Function: cos(x) + i*sin(x)');
grid on;

view([30 10]);

% Subplot for the polar plot
subplot(2, 1, 2);
hpolar = compass(y_cos(1), y_sin(1));
hold on;
hcos = quiver(0, 0, y_cos(1), 0, 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);
hsin = quiver(0, 0, 0, y_sin(1), 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);
cos_label = text(1.1, 0, 'cos', 'Color', 'r', 'FontSize', 12);
sin_label = text(0, 1.1, 'sin', 'Color', 'b', 'FontSize', 12);
title('Polar Coordinates');
hold off;

% Animation loop
for k = 1:length(t)
    % Update the position of the dots
    set(dot1, 'XData', t(k), 'YData', y_cos(k), 'ZData', y_sin(k));
    set(dot2, 'XData', t(k), 'YData', y_cos(k), 'ZData', y_sin(k));
    
    % Update the compass plot
    set(hpolar(1), 'XData', [0 y_cos(k)], 'YData', [0 y_sin(k)]);    
    
    % Update the quivers
    set(hcos, 'UData', y_cos(k), 'VData', 0);
    set(hsin, 'UData', 0, 'VData', y_sin(k));
    
    % Update the labels
    set(cos_label, 'Position', [y_cos(k) + 0.1, 0, 0]);
    set(sin_label, 'Position', [0, y_sin(k) + 0.1, 0]);
    
    % Pause to create animation effect
    pause(0.01);
end
%% WE CALL ONE AXIS THE REAL AND THE OTHER THE IMAGINARY AXIS
%% SLIDES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%  COMPLEX NUMBERS AND EULER'S FORMULA     %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% We convince ourselves that exp(1i * x) = cos(x) + i*sin(x)

% Define the range of t
n_cycles = 5;
t = linspace(0, n_cycles*2*pi, n_cycles*100);

% Generate the trigonometric form
y_cos = cos(t);
y_sin = sin(t);

% Generate the exponential form
y_exp = exp(1i * t);

% Create the figure
figure;
% Plot the trigonometric form
subplot(3, 1, 1);
plot3(t, cos(t), sin(t), 'r', 'LineWidth', 2);

legend('cos(x) + i*sin(x)', 'Location', 'northeast');
xlabel('x');
ylabel('Real Part');
zlabel('Imaginary Part');

title('Trigonometric Form');
grid on; view([30 10]);

subplot(3, 1, 2);
plot3(t, real(y_exp), imag(y_exp), 'b', 'LineWidth', 2);
title('Exponential Form');
xlabel('x');
ylabel('Real Part');
zlabel('Imaginary Part');
legend('e^{ix}', 'Location', 'northeast');
grid on; view([30 10]);

subplot(3, 1, 3);
plot3(t, cos(t), sin(t), 'r', 'LineWidth', 2);
hold on
plot3(t, real(y_exp), imag(y_exp), 'b--', 'LineWidth', 2);
title('Trigonometric Form vs Exponential Form');
xlabel('x');
ylabel('Real Part');
zlabel('Imaginary Part');
legend('cos(x) + i*sin(x)', 'e^{ix}', 'Location', 'northeast');
grid on; view([30 10]);

%% We create waves with different frequencies 
% Define the range of x
n_cycles = 1;
SR = 512;
t = linspace(0, n_cycles*2*pi, n_cycles*SR);

% Generate the exponential form
freq_a = 1;
freq_b = 5;

y_exp = exp(1i * t * freq_a);

plot3(t, real(y_exp), imag(y_exp), 'r', 'LineWidth', 2);
hold on;
y_exp = exp(1i * t * freq_b);

plot3(t, real(y_exp), imag(y_exp), 'b', 'LineWidth', 2);

title('Exponential Form - different frequencies');
xlabel('x');
ylabel('Real Part');
zlabel('Imaginary Part');
legend({'Frequency = 1','Frequency = 5'}, 'Location', 'northeast');
grid on;
%% We create waves with different amplitudes
% Define the range of x
figure;
n_cycles = 1;
SR = 512;
t = linspace(0, n_cycles*2*pi, n_cycles*SR);

% Generate the exponential form
freq_a = 1;
freq_b = 5;
amp_a = 1;
amp_b = 1.5;
y_exp = amp_a*exp(1i * t * freq_a);

% plot in 3D
plot3(t, real(y_exp), imag(y_exp), 'r', 'LineWidth', 2);
hold on;
y_exp = amp_b*exp(1i * t * freq_b);

plot3(t, real(y_exp), imag(y_exp), 'b', 'LineWidth', 2);

title('Exponential Form - different frequencies');
xlabel('x');
ylabel('Real Part');
zlabel('Imaginary Part');
legend({'Frequency = 1  (amp = 1)','Frequency = 5 (amp = 1.5)'}, ...
    'Location', 'northeast');
grid on;


%% We create waves with different phase offsets

% Define the range of x
n_cycles = 1;
SR = 512;
t = linspace(0, n_cycles*2*pi, n_cycles*SR);

% Generate the exponential form
freq = 1;
phase_offset = pi;
y_exp = exp(1i * (t * freq));

plot3(t, real(y_exp), imag(y_exp), 'r', 'LineWidth', 2);
hold on;
y_exp = exp(1i * (t * freq + phase_offset));

plot3(t, real(y_exp), imag(y_exp), 'b', 'LineWidth', 2);

title('Exponential Form - different frequencies');
xlabel('x');
ylabel('Real Part');
zlabel('Imaginary Part');
e_legend = legend({'Frequency = 1','Frequency = 1 (phase offset = pi)'}, 'Location', 'northeast');
grid on;

%% Exercises
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%   Exercises (Mrugank) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -Play with sine waves: Create waves with different frequencies, phase and
% amplitude

%% SLIDES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Discrete Fourier Transform    %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% let's make a signal that contains different frequencies and noise
% Define the parameters
N = 512; % Number of samples
n = 0:N-1; % Sample indices
t = n./N; % instead of using linspace we define t explicitly
% Generate an example signal (FREQUENCIES are 2 and 10) with some noise
X = sin(2 * pi * 2 * t) + sin(2 * pi * 10 * t) + ...
    cos(2 * pi * 10 * t + pi) + rand(1, N); % cosine has phase offset pi
%% plot the signal
figure; 
plot(X); xlim([0 N-1])
%% We want to get a fourier series -- let's start with a single frequency
k = 2; % Frequency index

% Generate the complex exponential (basis function) for frequency k
W = exp(-1i * 2 * pi * k * n / N);

% Compute the inner product (DFT coefficient for frequency k)
X_k = sum(X .* W);

% Compute the inner product (DFT coefficient for frequency k) via vector
% multiplication!
X_k_vm = W * X';

power = real(X_k)^2 + imag(X_k)^2;
% Display the results
fprintf('Power at frequency 2: %f', power);


%% -- let's continue with more than one frequency
% We also: 
% (2) use a longer signal,
% (3) define everything in Hz (1/s)
% (4) add a high amplitude signal at 6 hz

% Define the parameters
N = 2048; % Number of samples (extended signal length)
frequencies = 1:12; % Frequencies from 1 to 12 Hz
n = 0:N-1; % Sample indices
SR = 512;  % Sample Rate
dT = 1/SR; % Sample interval (size of a time-step)
t = n.*dT; % time vector
% Generate a sample signal (FREQUENCIES sin2, sin10, 3Xsin6, and cos10 hz
X = 2*sin(2 * pi * 2 * t) + sin(2 * pi * 10 * t) + ...
    3*cos(2 * pi * 6 * t + pi) ... 
    + cos(2 * pi * 10 * t + pi)+ rand(1, N); % cosine has a phase offset of pi
%% plot the signal
figure; 
plot(X); xlim([0 N-1])
%% Compute the power in the same way as before, however, 
% we loop through more frequencies

% Initialize an array to store the power for each frequency
power = zeros(size(frequencies));

% Compute the power for each frequency
for k = frequencies
    % Generate the complex exponential (basis function) for frequency k
    W = exp(-1j * 2 * pi * k * t);
    
    % Compute the inner product (DFT coefficient for frequency k)
    % we use vector multiplication! 
    X_k =  W * X';
    
    % Compute the power
    power(k) = real(X_k)^2 + imag(X_k)^2;
end

% Display the power for each frequency
for k = frequencies
    fprintf('Power at frequency %d Hz: %f\n', k, power(k));
end

%% Plot the barplot of power versus frequency
bar(frequencies, power)
xlabel('Frequency (Hz)')
ylabel('Power')
title('Power vs Frequency')
%% EXERCISE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%   Exercises with EEG (Mrugank) %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - select an electrode via vector multiplication
% - compute the power of the 60Hz frequency at that channel by taking the
% - inner product with a complex sine wave
% - compute the 60Hz power at all channels at once via matrix multiplication 

% - compare the 60Hz power to 50 Hz and 70 Hz power
% - do this for the alpha frequency
% - make a topoplot of alpha power
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% -- let's continue with ALL frequencies and use Matrix multiplication
% Create the Signal again:
% Define the parameters
SR = 512;             % Sampling frequency
dT = 1/SR;            % Sampling period
N = 2048;             % Length of signal 
n = 0:N-1;            % Sample indices
t = n*dT;             % Time vector
x = (0:N-1) ./ N;      % Sample vector
% Create a simple signal (sum of sine waves)
X = 2*sin(2 * pi * 2 * n / SR) + sin(2 * pi * 10 * n / SR) + ...
    3*cos(2 * pi * 6 * n / SR + pi) ... 
    + cos(2 * pi * 10 * n / SR + pi) + rand(1, N); % cosine has a phase offset of pi

% ----- this is could be interesting to show later. 
% % to add non-stationary noise:
% % 5 Hz  noise
% f_noise = 5; % Frequency of the noise
% decay = exp(-n/(N/4)); % Exponential decay function
% ns_noise = 10*(decay .* sin(2 * pi * f_noise * n / SR));
% 
% % Add the line noise to the original signal
% X_noisy = X + ns_noise;
% X = X_noisy;

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

% --- alternatively we could compute the DFT with a loop over frequencies
% X_dft = zeros(1, N);
% for k = 0:N-1
%     for n = 0:N-1
%         X_dft(k+1) = X_dft(k+1) + X(n+1) * exp(-1i * 2 * pi * k * n / N);
%     end
% end

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

% if we wanted to do this with a loop...
% X_idft = zeros(1, N);
% 
% % Compute the inverse DFT
% for n = 0:N-1
%     for k = 0:N-1
%         X_idft(n+1) = X_idft(n+1) + X_dft(k+1) * exp(1i * 2 * pi * k * n / N);
%     end
%     X_idft(n+1) = X_idft(n+1) / N;
% end

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

% ==> we have now high-pass filtered our data.
%% SLIDES

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%       CONVOLUTION              %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% Define different filter kernels
low_pass_kernel = ones(1, 5) / 5; % Simple moving average filter
high_pass_kernel = [-1, -1, 4, -1, -1]; % Example of a high-pass filter

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
%%%%%%%%%   Redo our highpass filter via convolution %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% recreate the signal from before
SR = 512;             % Sampling frequency
dT = 1/SR;            % Sampling period
N = 2048;             % Length of signal 
n = 0:N-1;            % Sample indices
t = n*dT;             % Time vector
% Create a simple signal (sum of sine waves)
X = 2*sin(2 * pi * 2 * n / SR) + sin(2 * pi * 10 * n / SR) + ...
    3*cos(2 * pi * 6 * n / SR + pi) ... 
    + cos(2 * pi * 10 * n / SR + pi) + rand(1, N); % cosine has a phase offset of pi

figure;
plot(t, X);
title('Original Signal');
xlabel('t (seconds)');
ylabel('X(t)');
%% Compute an FFT this time
% here we compute all the dot products in one go (eeach row is a frequency)
X_fft = fft(X);

% compute the power
power = real(X_fft).^2 + imag(X_fft).^2; % NOTE: this is elementwise .^

% for the correct axis, we need to rescale to 1/second
frequencies = [0:N-1] .* (SR / N);

figure
% Plot the barplot of power versus frequency (we only plot up to Nyquist
% frequency. Everything above is mirrored.
bar(frequencies(1:N/2), power(1:N/2))
xlabel('Frequency (Hz)')
ylabel('Power')
title('Power vs Frequency');
xlim([0 20])

%% Recompute our LP filtered data:

X_fft_filt = X_fft;
X_fft_filt(:, 1:36) = 0 ;
% the high frequencies are mirrored
X_fft_filt(:, N-36:N) = 0 ;
X_ifft_filtered = ifft(X_fft_filt);

%% Let's compute the filter kernel from our previous high-pass filter

% also ifft the kernel
filter_spec  = ones(1, N);
filter_spec(:, 1:36) = 0 ;
% the high frequencies are mirrored
filter_spec(:, N-36:N) = 0 ;
filter_kernel = fftshift(ifft(filter_spec));
X_convfiltered = conv(X, filter_kernel, 'same');
%% plot the filter kernel
figure; plot(real(filter_kernel));

%% Compare the Filtered Signals:
% Plot both signals for comparison
figure;
subplot(3,1,1);
plot(real(X_ifft_filtered), 'r');
title('ifft filtered Signal');
xlabel('samples');
ylabel('X\_filtered(t)');

subplot(3,1,2);
plot(real(X_convfiltered), 'b');
title('Convolution filtered Signal');
xlabel('samples');
ylabel('X\_filtered(t)');
subplot(3,1,3);
plot(real(X_ifft_filtered), 'r');
hold on;
plot(real(X_convfiltered), 'b--');

title('Both Signals');
xlabel('samples');
ylabel('X\_filtered(t)');

