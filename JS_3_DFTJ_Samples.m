clc;
clear all;

% Parameters for the original chirp waveform
start_frequency = -6e6;          % Start frequency in Hz
end_frequency = 6e6;             % End frequency in Hz
pulse_duration = 20e-6;          % Pulse duration in seconds
total_duration = 100e-6;         % Total duration to visualize

% Sampling frequency
fs = 70e6;  % Adjust as needed

% Number of samples
num_samples = 1;

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

    

    % % Plot composite signal waveform
    % figure;
    % plot(t * 1e6, real(composite_signal), 'r', 'LineWidth', 1.5);
    % title(sprintf('Sample %d: Dense False Target Jamming (Time Domain Waveform)', sample));
    % xlabel('Time (\mu s)');
    % ylabel('Normalized Amplitude');

    % Plot spectrum of the composite signal
    figure;
    spectrogram(composite_signal, hann(256), 250, 1024, fs, 'centered', 'yaxis');
    title(sprintf('Sample %d: Dense False Target Jamming (Time Frequency Waveform)', sample));
    set(gca, 'YDir','reverse');
end
