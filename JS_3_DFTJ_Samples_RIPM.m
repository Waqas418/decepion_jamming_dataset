clc;
clear all;

% Define the path to save the spectrogram images
folderPath = 'C:\Users\waq\Desktop\Saved Images\New folder';

% Create the folder if it doesn't exist
if ~exist(folderPath, 'dir')
    mkdir(folderPath);
end

% Parameters for the original chirp waveform
start_frequency = -6e6;          % Start frequency in Hz
end_frequency = 6e6;             % End frequency in Hz
pulse_duration = 20e-6;          % Pulse duration in seconds
total_duration = 100e-6;         % Total duration to visualize

% Sampling frequency
fs = 70e6;  % Adjust as needed

% Number of samples
num_samples = 2;

for sample = 1:num_samples
    % Generate random delay for the first pulse between 1 and 10 microseconds
    delay_first_pulse = randi([1, 10]) * 1e-6;

    % Time vector
    t = linspace(0, total_duration, total_duration * fs);

    % Generate original chirp waveform with pulse duration and random delay for the first pulse
    start_time_original = 25e-6 + delay_first_pulse; % Start time for the original chirp waveform with delay
    chirp_waveform_original = chirp(t - start_time_original, start_frequency, pulse_duration, end_frequency, 'linear', 'complex');

    % Ensure original chirp waveform is zero before 0 microseconds and after 20 microseconds
    chirp_waveform_original(t < start_time_original) = 0;
    chirp_waveform_original(t > (start_time_original + pulse_duration)) = 0;

    % Generate additional chirp waveforms at specified intervals with the same delay
    start_times_additional = [25, 30, 35, 40, 45] * 1e-6 + delay_first_pulse;  % Start times for additional chirps with delay
    chirp_waveforms_additional = zeros(length(start_times_additional), length(t));
    for i = 1:length(start_times_additional)
        chirp_waveforms_additional(i, :) = chirp(t - start_times_additional(i), start_frequency, pulse_duration, end_frequency, 'linear', 'complex');
        % Ensure additional chirp waveform is zero before start time and after end time
        chirp_waveforms_additional(i, t < start_times_additional(i)) = 0;
        chirp_waveforms_additional(i, t > (start_times_additional(i) + pulse_duration)) = 0;
    end

    % Combine all chirp waveforms
    composite_signal = chirp_waveform_original + sum(chirp_waveforms_additional, 1);

    % Calculate the modulus (magnitude) and phase of the complex chirp waveform
    modulus_chirp = abs(composite_signal);
    phase_chirp = angle(composite_signal);

    % Plot the spectrogram of the real part of the composite signal
    figure;
    spectrogram(real(composite_signal), hann(256), 250, 1024, fs, 'centered', 'yaxis');
    title(sprintf('Dense False Target Jamming (Spectrogram - Real Part) - Sample %d', sample));
    set(gca, 'YDir', 'reverse');
    % Save the spectrogram image
    filename = fullfile(folderPath, sprintf('Spectrogram_Real_Part_Sample_%d.png', sample));
    saveas(gcf, filename);
    % 
    % % Plot the spectrogram of the imaginary part of the composite signal
    % figure;
    % spectrogram(imag(composite_signal), hann(256), 250, 1024, fs, 'centered', 'yaxis');
    % title(sprintf('Dense False Target Jamming (Spectrogram - Imaginary Part) - Sample %d', sample));
    % set(gca, 'YDir', 'reverse');
    % % Save the spectrogram image
    % filename = fullfile(folderPath, sprintf('Spectrogram_Imaginary_Part_Sample_%d.png', sample));
    % saveas(gcf, filename);
    % 
    % % Plot the spectrogram of the modulus (magnitude) of the composite signal
    % figure;
    % spectrogram(modulus_chirp, hann(256), 250, 1024, fs, 'centered', 'yaxis');
    % title(sprintf('Dense False Target Jamming (Spectrogram - Modulus) - Sample %d', sample));
    % set(gca, 'YDir', 'reverse');
    % % Save the spectrogram image
    % filename = fullfile(folderPath, sprintf('Spectrogram_Modulus_Sample_%d.png', sample));
    % saveas(gcf, filename);
    % 
    % % Plot the spectrogram of the phase of the composite signal
    % figure;
    % spectrogram(phase_chirp, hann(256), 250, 1024, fs, 'centered', 'yaxis');
    % title(sprintf('Dense False Target Jamming (Spectrogram - Phase) - Sample %d', sample));
    % set(gca, 'YDir', 'reverse');
    % % Save the spectrogram image
    % filename = fullfile(folderPath, sprintf('Spectrogram_Phase_Sample_%d.png', sample));
    % saveas(gcf, filename);
end
