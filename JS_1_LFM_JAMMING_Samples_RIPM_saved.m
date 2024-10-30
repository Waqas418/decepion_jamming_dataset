clc;
clear all;

% % Define the path to save the spectrogram images
% folderPath = 'L:\JS_1_LFM_JAMMING_Samples_RIPM\Real';
% 
% % Create the folder if it doesn't exist
% % if ~exist(folderPath, 'dir')
% %     mkdir(folderPath);
% end

% Parameters for the original chirp waveform
start_frequency = -6e6;          % Start frequency in Hz
end_frequency = 6e6;             % End frequency in Hz
pulse_duration = 20e-6;          % Pulse duration in seconds
total_duration = 100e-6;         % Total duration to visualize
start_time = 55e-6;              % Start time for the chirp waveform

% Sampling frequency
fs = 30e6;  % Adjust as needed

% Time vector
t = linspace(0, total_duration, total_duration * fs);

% Generate original chirp waveform with pulse duration
chirp_waveform_original = chirp(t - start_time, start_frequency, pulse_duration, end_frequency, 'linear','complex');

% Ensure original chirp waveform is zero before 0 microseconds and after 20 microseconds
chirp_waveform_original(t < start_time) = 0;
chirp_waveform_original(t > (start_time + pulse_duration)) = 0;

% Generate 5 samples with random noise levels between 0.25 and 1
num_samples = 1500;
for i = 1:num_samples
    % Generate pure noise with random amplitude between 0.25 and 1
    noise_amplitude = 0.25 + (1 - 0.25) * rand();
    noise = (rand(size(t)) - 0.5) * 2 * noise_amplitude;
    
    % Add noise to the chirp waveform
    chirp_with_noise = chirp_waveform_original + noise;

    % Calculate the modulus (magnitude) and phase of the complex chirp waveform
    modulus_chirp = abs(chirp_with_noise);
    phase_chirp = angle(chirp_with_noise);

    % Plot the spectrogram of the real part of the original chirp waveform with the additional delay
    figure;
    spectrogram(real(chirp_with_noise), hann(256), 250, 1024, fs, 'centered', 'yaxis');
    title(sprintf('LFM Jamming Signal (Spectrogram - Real Part) - Sample %d', i));
    set(gca, 'YDir', 'reverse');
    % % Save the spectrogram image
    % filename = fullfile(folderPath, sprintf('Spectrogram_Real_Part_Sample_%d.png', i));
    % saveas(gcf, filename);

    % Plot the spectrogram of the imaginary part of the original chirp waveform with the additional delay
    figure;
    spectrogram(imag(chirp_with_noise), hann(256), 250, 1024, fs, 'centered', 'yaxis');
    title(sprintf('LFM Jamming Signal (Spectrogram - Imaginary Part) - Sample %d', i));
    set(gca, 'YDir', 'reverse');
    % % Save the spectrogram image
    % filename = fullfile(folderPath, sprintf('Spectrogram_Imaginary_Part_Sample_%d.png', i));
    % saveas(gcf, filename);
    % 
    % Plot the spectrogram of the modulus (magnitude) of the original chirp waveform
    figure;
    spectrogram(modulus_chirp, hann(256), 250, 1024, fs, 'centered', 'yaxis');
    title(sprintf('LFM Jamming Signal (Spectrogram - Modulus) - Sample %d', i));
    set(gca, 'YDir', 'reverse');
    % % Save the spectrogram image
    % filename = fullfile(folderPath, sprintf('Spectrogram_Modulus_Sample_%d.png', i));
    % saveas(gcf, filename);
    % 
    % Plot the spectrogram of the phase of the original chirp waveform
    figure;
    spectrogram(phase_chirp, hann(256), 250, 1024, fs, 'centered', 'yaxis');
    title(sprintf('LFM Jamming Signal (Spectrogram - Phase) - Sample %d', i));
    set(gca, 'YDir', 'reverse');
    % % Save the spectrogram image
    % filename = fullfile(folderPath, sprintf('Spectrogram_Phase_Sample_%d.png', i));
    % saveas(gcf, filename);
end
