clc; 
clear all; 
close all;

% Reading the original signal
fs = 44100; % Sampling frequency
[y, fs] = audioread('nokia.mp3'); % Load the audio signal
y = y(:, 1); % Use only one channel if stereo
t = (0:length(y)-1) / fs; % Time vector

% Add noise to the signal
noise_level_dB = 10; % Noise level in dB
signal_power = mean(y.^2); % Power of the original signal
noise_power = signal_power / (10^(noise_level_dB / 10)); % Noise power calculation
noise = sqrt(noise_power) * randn(size(y)); % Generate Gaussian noise
y_noisy = y + noise; % Add noise to the signal



% Analyze frequency spectrum before filtering
N = length(y); % Number of samples
f = (0:N-1) * (fs / N); % Frequency vector
fft_y = fft(y); % Compute FFT of the original signal
fft_y_noisy = fft(y_noisy); % Compute FFT of the noisy signal

% Plot the frequency spectrum of original and noisy signals
figure;

% Original signal spectrum
subplot(2, 1, 1);
plot(f(1:N/2), abs(fft_y(1:N/2)));
title('Frequency Spectrum of Original Signal');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
grid on;

% Noisy signal spectrum
subplot(2, 1, 2);
plot(f(1:N/2), abs(fft_y_noisy(1:N/2)));
title('Frequency Spectrum of Noisy Signal');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
grid on;


% Remove noise using a low-pass filter
fc = 20000; % Cutoff frequency (Hz)
%[b, a] = butter(6, fc / (fs / 2), 'low'); % 6th order Butterworth low-pass filter
%y_filtered = filter(b, a, y_noisy); % Apply the filter
[b, a] = butter(10, fc / (fs / 2), 'low'); % 10th-order Butterworth filter
y_filtered = filter(b, a, y_noisy);


% Calculate SNR
snr_noisy = 10 * log10(sum(y.^2) / sum((y_noisy - y).^2)); % SNR of noisy signal
snr_filtered = 10 * log10(sum(y.^2) / sum((y_filtered - y).^2)); % SNR of filtered signal

% Display SNR values
fprintf('SNR of before filter signal: %.2f dB\n', snr_noisy);
fprintf('SNR of after filter signal: %.2f dB\n', snr_filtered);

% Plot the waveforms
figure;

% Original signal
subplot(3, 1, 1);
plot(t, y);
title('Waveform of Nokia MP3 Sound Original Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Noisy signal
subplot(3, 1, 2);
plot(t, y_noisy);
title('Waveform of Nokia MP3 Sound with adding noise so Noisy Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Filtered signal
subplot(3, 1, 3);
plot(t, y_filtered);
title('Waveform of Nokia MP3 Sound with Filtered Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;
