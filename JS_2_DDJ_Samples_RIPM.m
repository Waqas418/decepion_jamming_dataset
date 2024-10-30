clc;
clear all;

% Define the path to save the spectrogram images
folderPath = 'C:\Users\ASD LAB 21\Desktop\New folder\STFT';

% Create the folder if it doesn't exist
if ~exist(folderPath, 'dir')
    mkdir(folderPath);
end

% Parameters for the original chirp waveform
start_frequency = -6e6;          % Start frequency in Hz
end_frequency = 6e6;             % End frequency in Hz
pulse_duration = 20e-6;          % Pulse duration in seconds
total_duration = 100e-6;         % Total duration to visualize
start_time = 52e-6;              % Start time for the chirp waveform

% Sampling frequency
fs = 75e6;  % Adjust as needed

% Number of samples to generate
num_samples = 2;

% Initialize an array to store the delay_additional values
delay_additional_values = zeros(1, num_samples);

% Generate delay_additional values within the range [1e-6, 10e-6]
for i = 1:num_samples
    delay_additional_values(i) = (10 - 1) * rand() * 1e-6 + 1e-6; % Randomly sample within the range [1e-6, 10e-6]
end

% Loop to generate samples with different delay_additional values
for i = 1:num_samples
    % Time vector
    t = linspace(0, total_duration, total_duration * fs);
    
    % Generate original chirp waveform with pulse duration and additional delay
    chirp_waveform_original = chirp(t - start_time - delay_additional_values(i), start_frequency, pulse_duration, end_frequency, 'linear','complex');
    
    % Ensure original chirp waveform is zero before 0 microseconds and after 20 microseconds
    chirp_waveform_original(t < (start_time + delay_additional_values(i))) = 0;
    chirp_waveform_original(t > (start_time + delay_additional_values(i) + pulse_duration)) = 0;
    
    % Calculate the modulus (magnitude) and phase of the complex chirp waveform
    modulus_chirp = abs(chirp_waveform_original);
    phase_chirp = angle(chirp_waveform_original);

    % Plot the spectrogram of the real part of the original chirp waveform with the additional delay
    figure;
    spectrogram(chirp_waveform_original), hann(256), 250, 1024, fs, 'centered', 'yaxis');
    title(sprintf('Distance Deception Jamming (Spectrogram - Real Part) - Sample %d', i));
    set(gca, 'YDir', 'reverse');
    % % Save the spectrogram image
    % filename = fullfile(folderPath, sprintf('Spectrogram_Real_Part_Sample_%d.png', i));
    % saveas(gcf, filename);
    % 
    % % Plot the spectrogram of the imaginary part of the original chirp waveform with the additional delay
    % figure;
    % spectrogram(imag(chirp_waveform_original), hann(256), 250, 1024, fs, 'centered', 'yaxis');
    % title(sprintf('Distance Deception Jamming (Spectrogram - Imaginary Part) - Sample %d', i));
    % set(gca, 'YDir', 'reverse');
    % % Save the spectrogram image
    % filename = fullfile(folderPath, sprintf('Spectrogram_Imaginary_Part_Sample_%d.png', i));
    % saveas(gcf, filename);
    % 
    % % Plot the spectrogram of the modulus (magnitude) of the original chirp waveform
    % figure;
    % spectrogram(modulus_chirp, hann(256), 250, 1024, fs, 'centered', 'yaxis');
    % title(sprintf('Distance Deception Jamming (Spectrogram - Modulus) - Sample %d', i));
    % set(gca, 'YDir', 'reverse');
    % % Save the spectrogram image
    % filename = fullfile(folderPath, sprintf('Spectrogram_Modulus_Sample_%d.png', i));
    % saveas(gcf, filename);
    % 
    % % Plot the spectrogram of the phase of the original chirp waveform
    % figure;
    % spectrogram(phase_chirp, hann(256), 250, 1024, fs, 'centered', 'yaxis');
    % title(sprintf('Distance Deception Jamming (Spectrogram - Phase) - Sample %d', i));
    % set(gca, 'YDir', 'reverse');
    % % Save the spectrogram image
    % filename = fullfile(folderPath, sprintf('Spectrogram_Phase_Sample_%d.png', i));
    % saveas(gcf, filename);
end
