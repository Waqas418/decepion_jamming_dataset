clc;
clear all;

% Parameters for the original chirp waveform
start_frequency = -6e6;          % Start frequency in Hz
end_frequency = 6e6;             % End frequency in Hz
pulse_duration = 20e-6;          % Pulse duration in seconds
total_duration = 100e-6;         % Total duration to visualize
start_time = 52e-6;              % Start time for the chirp waveform

% Sampling frequency
fs = 75e6;  % Adjust as needed

% Number of samples to generate
num_samples =1;

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
    
    % % Plot the original chirp waveform with the additional delay
    % figure;
    % plot(t * 1e6, real(chirp_waveform_original));
    % title(sprintf('Distance Deception Jamming (Time Domain Waveform) - Sample %d', i));
    % xlabel('Time (\mu s)');
    % ylabel('Normalized Amplitude');

    % Plot the spectrogram
    figure;
    spectrogram(chirp_waveform_original, hann(256), 250, 1024, fs,'centered', 'yaxis');
    title(sprintf('Distance Deception Jamming (Time Frequency Spectrogram) - Sample %d', i));
    set(gca, 'YDir','reverse');

%     % Save the image with the new destination path
% saveas(gcf, ['L:\Simulation Framework\Jamming Signals\JS_2_DDJ\Dataset_Samples\New folder', sprintf('Spectrogram_Sample_%d.png', i)]);

end
