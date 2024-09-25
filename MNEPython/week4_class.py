#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 10:45:42 2024

@author: sebastianmichelmann
"""

import numpy as np
import matplotlib.pyplot as plt

# Define the parameters
SR = 512             # Sampling frequency
dT = 1/SR            # Sampling period
N = 2048             # Length of signal 
n = np.arange(N)     # Sample indices
t = n * dT           # Time vector
x = np.arange(N) / N # Sample vector runs from 0 to 1

# Create a simple signal (sum of sine waves)
X = 2 * np.sin(2 * np.pi * 2 * n / SR) + np.sin(2 * np.pi * 10 * n / SR) + \
    3 * np.cos(2 * np.pi * 6 * n / SR + np.pi) + \
    np.cos(2 * np.pi * 10 * n / SR + np.pi) + np.random.rand(N)

# Plot the signal
plt.figure()
plt.plot(t, X)
plt.title('Original Signal')
plt.xlabel('t (seconds)')
plt.ylabel('X(t)')
plt.show()
#%%
# Compute the DFT Matrix
W_dft = np.zeros((N, N), dtype=complex)
for k in range(N):
    W_dft[k, :] = np.exp(-1j * 2 * np.pi * k * x)

# Compute the DFT
X_dft = np.dot(W_dft, X)

# Compute the power of the complex Fourier coefficients
power = np.real(X_dft)**2 + np.imag(X_dft)**2

# Rescale to 1/second (Hertz)
frequencies = np.arange(N) * (SR / N)

# Plot the power vs frequency
plt.figure()
plt.bar(frequencies[:N//2], power[:N//2])
plt.xlabel('Frequency (Hz)')
plt.ylabel('Power')
plt.title('Power vs Frequency')
plt.xlim([0, 20])
plt.show()

# Compute the inverse DFT Matrix
W_idft = np.zeros((N, N), dtype=complex)
for n in range(N):
    W_idft[n, :] = np.exp(1j * 2 * np.pi * np.arange(N) * n / N)

# Compute the inverse DFT
X_idft = np.dot(W_idft, X_dft) / N

# Compare the Original and Reconstructed Signals
plt.figure()
plt.subplot(3, 1, 1)
plt.plot(t, X, 'r')
plt.title('Original Signal')
plt.xlabel('t (seconds)')
plt.ylabel('X(t)')

plt.subplot(3, 1, 2)
plt.plot(t, X_idft, 'b')
plt.title('Reconstructed Signal')
plt.xlabel('t (seconds)')
plt.ylabel('X_reconstructed(t)')

plt.subplot(3, 1, 3)
plt.plot(t, X, 'r')
plt.plot(t, X_idft, 'b--')
plt.title('Both Signals')
plt.xlabel('t (seconds)')
plt.show()
#%%
# Compute the Inverse DFT, but set the low frequencies to 0
W_idft = np.conj(W_dft)
X_dft_filt = X_dft.copy()
X_dft_filt[:36] = 0
X_dft_filt[-36:] = 0

X_idft_filtered = np.dot(W_idft, X_dft_filt) / N

# Compare the Original and Filtered Signals
plt.figure()
plt.subplot(3, 1, 1)
plt.plot(t, X, 'r')
plt.title('Original Signal')
plt.xlabel('t (seconds)')
plt.ylabel('X(t)')

plt.subplot(3, 1, 2)
plt.plot(t, X_idft_filtered, 'b')
plt.title('Filtered Signal')
plt.xlabel('t (seconds)')
plt.ylabel('X_reconstructed(t)')

plt.subplot(3, 1, 3)
plt.plot(t, X, 'r')
plt.plot(t, X_idft_filtered, 'b--')
plt.title('Both Signals')
plt.xlabel('t (seconds)')
plt.show()
#%%
# Non-stationary signal
n = np.arange(N)     # Sample indices
f_ns = 6  # Frequency of the non-stationary signal
decay = np.exp(-n / (N / 4))  # Exponential decay function
X_6 = 4 * 3 * np.cos(2 * np.pi * f_ns * n / SR)
ns_X = decay * X_6

# Plot the non-stationary signal
plt.figure()
plt.plot(t, ns_X)
plt.title('Non-stationary Signal')
plt.xlabel('t (seconds)')
plt.ylabel('X(t)')
plt.show()
#%%
# Create the same signal as before
X = 2 * np.sin(2 * np.pi * 2 * n / SR) + np.sin(2 * np.pi * 10 * n / SR) + \
    3 * np.cos(2 * np.pi * 6 * n / SR + np.pi) + \
    np.cos(2 * np.pi * 10 * n / SR + np.pi) + np.random.rand(N)

# Change the 6Hz cosine for the non-stationary data
X_ns = 2 * np.sin(2 * np.pi * 2 * n / SR) + np.sin(2 * np.pi * 10 * n / SR) + \
    ns_X + np.cos(2 * np.pi * 10 * n / SR + np.pi) + np.random.rand(N)

# Plot the signals
plt.figure()
plt.subplot(2, 1, 1)
plt.plot(t, X)
plt.title('Previous Signal')
plt.xlabel('t (seconds)')
plt.ylabel('X(t)')

plt.subplot(2, 1, 2)
plt.plot(t, X_ns)
plt.title('Non-stationary Signal')
plt.xlabel('t (seconds)')
plt.ylabel('X(t)')
plt.show()

# Compute the DFT of the non-stationary signal
X_dft_ns = np.dot(W_dft, X_ns)

# Compute the power of the non-stationary signal
power_ns = np.real(X_dft_ns)**2 + np.imag(X_dft_ns)**2

# Plot the power vs frequency for the non-stationary signal
plt.figure()
plt.subplot(2, 1, 1)
plt.bar(frequencies[:N//2], power[:N//2])
plt.xlabel('Frequency (Hz)')
plt.ylabel('Power')
plt.title('Power vs Frequency')
plt.xlim([0, 20])

plt.subplot(2, 1, 2)
plt.bar(frequencies[:N//2], power_ns[:N//2])
plt.xlabel('Frequency (Hz)')
plt.ylabel('Power')
plt.title('Power vs Frequency')
plt.xlim([0, 20])
plt.show()

# Convolution for edge detection
signal = np.concatenate([np.zeros(512), np.ones(512), np.zeros(512)])
kernel = [1, 0, -1]  # Example of a simple edge-detection kernel

# Perform convolution
convolved_signal = np.convolve(signal, kernel, 'same')

# Plot the original signal and the convolved signal
plt.figure()
plt.subplot(2, 1, 1)
plt.plot(signal, 'r')
plt.title('Original Signal')
plt.xlabel('Sample Index')
plt.ylabel('Amplitude')

plt.subplot(2, 1, 2)
plt.plot(convolved_signal, 'b')
plt.title('Convolved Signal')
plt.xlabel('Sample Index')
plt.ylabel('Amplitude')
plt.show()
#%% DIFFERENT FILTER KERNELS

# Play with different filter kernels
low_pass_kernel = np.ones(5) / 5  # Simple moving average filter
high_pass_kernel = np.array([-1, -1, 4, -1, -1])  # Example of a simple high-pass filter

# Create a Gaussian filter manually
sigma = 1  # Standard deviation
kernel_size = 51  # Size of the kernel
t = np.linspace(-2, 2, kernel_size)
gaussian_kernel = np.exp(-t**2 / (2 * sigma**2))
gaussian_kernel = gaussian_kernel / np.sum(gaussian_kernel)  # Normalize

# Apply the low-pass filter
filtered_signal_low = np.convolve(signal, low_pass_kernel, 'same')

# Apply the high-pass filter
filtered_signal_high = np.convolve(signal, high_pass_kernel, 'same')

# Apply the Gaussian filter
filtered_signal_gaussian = np.convolve(signal, gaussian_kernel, 'same')

# Plot the original and filtered signals
plt.figure()
plt.subplot(5, 1, 1)
plt.plot(signal, 'r')
plt.title('Original Signal')
plt.xlabel('Sample Index')
plt.ylabel('Amplitude')
plt.ylim([-1, 2])
plt.xlim([1, len(signal)])

plt.subplot(5, 1, 2)
plt.plot(filtered_signal_low, 'b')
plt.title('Low-Pass Filtered Signal')
plt.xlabel('Sample Index')
plt.ylabel('Amplitude')
plt.ylim([-1, 2])
plt.xlim([1, len(signal)])

plt.subplot(5, 1, 3)
plt.plot(filtered_signal_high, 'g')
plt.title('High-Pass Filtered Signal')
plt.xlabel('Sample Index')
plt.ylabel('Amplitude')
plt.ylim([-1, 2])
plt.xlim([1, len(signal)])

plt.subplot(5, 1, 4)
plt.plot(gaussian_kernel, 'm')
plt.title('Gaussian Kernel')
plt.xlabel('Sample Index')
plt.ylabel('Amplitude')
plt.xlim([-100, len(gaussian_kernel) + 100])

plt.subplot(5, 1, 5)
plt.plot(filtered_signal_gaussian, 'm')
plt.title('Gaussian Filtered Signal')
plt.xlabel('Sample Index')
plt.ylabel('Amplitude')
plt.ylim([-1, 2])
plt.xlim([1, len(signal)])
plt.show()

# Convolution in time is multiplication in frequency
conv_length = len(signal) + len(gaussian_kernel) - 1
x = np.fft.fft(signal, conv_length)
y = np.fft.fft(gaussian_kernel, conv_length)
y_ifft = np.fft.ifft(x * y)
y_conv = np.convolve(signal, gaussian_kernel, 'same')

plt.figure()
plt.subplot(5, 1, 1)
plt.plot(signal, 'r')
plt.title('Original Signal')
plt.xlabel('Sample Index')
plt.ylabel('Amplitude')
plt.ylim([-1, 2])
plt.xlim([1, len(signal)])

plt.subplot(5, 1, 2)
plt.plot(gaussian_kernel, 'm')
plt.title('Gaussian Kernel')
plt.xlabel('Sample Index')
plt.ylabel('Amplitude')
plt.xlim([-100, len(gaussian_kernel) + 100])

plt.subplot(5, 1, 3)
plt.plot(y_conv, 'r')
plt.title('Gaussian Filtered Signal (convolution)')
plt.xlabel('Sample Index')
plt.ylabel('Amplitude')
plt.ylim([-1, 2])
plt.xlim([1, len(signal)])

plt.subplot(5, 1, 4)
plt.plot(y_ifft[int(np.ceil(len(gaussian_kernel) / 2)) + 1:len(y_conv) + int(np.floor(len(gaussian_kernel) / 2))], 'b')
plt.title('Gaussian Filtered Signal (multiplication in frequency)')
plt.xlabel('Sample Index')
plt.ylabel('Amplitude')
plt.ylim([-1, 2])
plt.xlim([1, len(signal)])

plt.subplot(5, 1, 5)
plt.plot(y_conv, 'r')
plt.plot(y_ifft[int(np.ceil(len(gaussian_kernel) / 2)) + 1:len(y_conv) + int(np.floor(len(gaussian_kernel) / 2))], 'b--')
plt.title('Comparison')
plt.xlabel('Sample Index')
plt.ylabel('Amplitude')
plt.ylim([-1, 2])
plt.xlim([1, len(signal)])
plt.show()
#%% MORLET WAVELETS
# Define the parameters for the Morlet wavelet
f0 = 6  # Central frequency of the wavelet
SR = 100
sigma = 1  # Standard deviation of the Gaussian window
t = np.arange(-5, 5 + 1/SR, 1/SR)  # Time vector
N_ml = len(t)

# Create the Morlet wavelet (frequency based on time)
A = 1 / (np.sqrt(sigma * np.sqrt(np.pi)) * np.sqrt(f0))  # amplitude scaling for the frequency-specific Morlet wavelet
morlet_wavelet = A * np.exp(2 * 1j * np.pi * f0 * t) * np.exp(-t**2 / (2 * sigma**2))

plt.figure()
plt.subplot(3, 1, 1)
plt.title('Cosine (real part of complex exponential)')
plt.xlabel('Sample Index')
plt.ylabel('Amplitude')
plt.plot(t, np.real(np.exp(2 * 1j * np.pi * f0 * t)))
plt.xlim([t[0], t[-1]])

plt.subplot(3, 1, 2)
plt.title('Gaussian')
plt.xlabel('Sample Index')
plt.ylabel('Amplitude')
plt.plot(t, np.exp(-t**2 / (2 * sigma**2)))
plt.xlim([t[0], t[-1]])

plt.subplot(3, 1, 3)
plt.title('Wavelet (real part)')
plt.xlabel('Sample Index')
plt.ylabel('Amplitude')
plt.plot(t, np.real(morlet_wavelet))
plt.xlim([t[0], t[-1]])
plt.show()

# Plot the Morlet wavelet and the transformed signal
plt.figure()
plt.subplot(3, 1, 1)
plt.plot(signal, 'r')
plt.title('Original Signal')
plt.xlabel('Sample Index')
plt.ylabel('Amplitude')
plt.ylim([-1, 2])
plt.xlim([1, len(signal)])

# Apply the Morlet wavelet to the signal (convolution)
wavelet_transformed_signal = np.convolve(signal, morlet_wavelet, 'same')

plt.subplot(3, 1, 2)
plt.plot(t, np.real(morlet_wavelet), 'r', label='Real Part')
plt.plot(t, np.imag(morlet_wavelet), 'b', label='Imaginary Part')
plt.title('Morlet Wavelet')
plt.xlabel('Time (seconds)')
plt.ylabel('Amplitude')
plt.legend()

plt.subplot(3, 1, 3)
plt.plot(np.arange(1, len(signal) + 1), np.real(wavelet_transformed_signal), 'r', label='Real Part')
plt.plot(np.arange(1, len(signal) + 1), np.imag(wavelet_transformed_signal), 'b', label='Imaginary Part')
plt.title('Wavelet Transformed Signal')
plt.xlabel('Sample Index')
plt.ylabel('Amplitude')
plt.legend()
plt.show()

# Convolve an impulse function with the Morlet wavelet
impulseX = np.zeros(5000)
impulseX[2501] = 1

y_conv = np.convolve(impulseX, morlet_wavelet, 'same')

plt.figure()
plt.subplot(2, 1, 1)
plt.plot(np.real(morlet_wavelet))
plt.subplot(2, 1, 2)
plt.plot(np.real(y_conv))
plt.show()

# Create X again
N = 5000
n = np.arange(N)  # Sample indices

f_ns = 6  # Frequency of the non-stationary signal

# Narrow Gaussian window shifted by 10 seconds
gaussian_window = np.exp(-0.5 * ((np.arange(-(N-1)/2-10*SR, (N-1)/2-10*SR + 1)) / (N/12))**2)
# Normalize to area of 1
gaussian_window_norm = gaussian_window / np.sum(gaussian_window)
# Modulated 6 Hz oscillation
X_6 = 800 * np.cos(2 * np.pi * f_ns * n / SR)
ns_X = gaussian_window_norm * X_6

X = np.sin(2 * np.pi * 2 * n / SR) + np.sin(2 * np.pi * 10 * n / SR) + ns_X + \
    np.cos(2 * np.pi * 10 * n / SR + np.pi) + np.random.rand(N)

plt.figure()
plt.subplot(2, 1, 1)
plt.title('Non-stationary part')
plt.plot(ns_X)
plt.subplot(2, 1, 2)
plt.plot(X)
plt.title('Combined Signal')
plt.xlim([0, N])
plt.show()

# Frequency spectrum of the wavelet
frequencies = SR * np.arange((len(morlet_wavelet) - 1) / 2 + 1) / len(morlet_wavelet)
morlet_fft = np.fft.fft(morlet_wavelet)

# Compute the power spectrum
power = np.real(morlet_fft)**2 + np.imag(morlet_fft)**2

plt.figure()
plt.bar(frequencies, power[:int((len(morlet_wavelet) - 1) / 2 + 1)])
plt.xlabel('Frequency (Hz)')
plt.ylabel('Power')
plt.title('Power Spectrum of Morlet Wavelet')
plt.xlim([0, 20])
plt.show()

# Frequency spectrum of the data
frequencies = SR * np.arange(N / 2) / N
X_fft = np.fft.fft(X)

# Compute the power spectrum
powerX = np.real(X_fft)**2 + np.imag(X_fft)**2

plt.figure()
plt.bar(frequencies, powerX[:N//2])
plt.xlabel('Frequency (Hz)')
plt.ylabel('Power')
plt.title('Power Spectrum of Original Data')
plt.xlim([0, 15])
plt.show()

# Filter the data
conv_length = len(X) + len(morlet_wavelet) - 1

X_morlet_filtered = np.fft.ifft(np.fft.fft(morlet_wavelet, conv_length) * np.fft.fft(X, conv_length))
X_morlet_convolved = np.convolve(X, morlet_wavelet, 'same')

plt.figure()
plt.subplot(4, 1, 1)
plt.plot(n / SR, X, 'r')
plt.xlabel('Time')
plt.title('Combined Signal')
plt.xlim([0, N / SR])

plt.subplot(4, 1, 2)
plt.plot(n / SR, np.real(X_morlet_filtered[500:len(X_morlet_filtered)-500]), 'r')
plt.xlabel('Time')
plt.title('Morlet Filtered (FFT)')
plt.xlim([0, N / SR])

plt.subplot(4, 1, 3)
plt.plot(n / SR, np.real(X_morlet_convolved), 'b')
plt.xlabel('Time')
plt.title('Morlet Filtered (Conv)')
plt.xlim([0, N / SR])

plt.subplot(4, 1, 4)
plt.plot(n / SR, np.real(X_morlet_filtered[500:len(X_morlet_filtered)-500]), 'r', label='FFT')
plt.plot(n / SR, np.real(X_morlet_convolved), 'b--', label='Conv')
plt.xlabel('Time')
plt.title('Both')
plt.legend()
plt.xlim([0, N / SR])
plt.show()

# Frequency spectrum of the filtered data
frequencies = SR * np.arange(N / 2) / N
X_fft = np.fft.fft(X_morlet_convolved)

# Compute the power spectrum
powerX = np.real(X_fft)**2 + np.imag(X_fft)**2

plt.figure()
plt.bar(frequencies, powerX[:N//2])
plt.xlabel('Frequency (Hz)')
plt.ylabel('Power')
plt.title('Power Spectrum of Wavelet-filtered Data')
plt.xlim([0, 15])
plt.show()
