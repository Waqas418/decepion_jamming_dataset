clc;
clear all;

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
num_samples = 5;
for i = 1:num_samples
    % Generate pure noise with random amplitude between 0.25 and 1
    noise_amplitude = 0.25 + (1 - 0.25) * rand();
    noise = (rand(size(t)) - 0.5) * 2 * noise_amplitude;
    
    % Add noise to the chirp waveform
    chirp_with_noise = chirp_waveform_original + noise;

    % Plot the time waveform
    figure;
    plot(t * 1e6, real(chirp_with_noise) / max(abs(chirp_waveform_original)) * 0.5);
    title(sprintf('LFM Jamming Signal (Time Domain Waveform) - Random Noise Level %.2f', noise_amplitude));
    xlabel('Time (\mu s)');
    ylabel('Normalized Amplitude');
    ylim([-0.5, 0.5]);

%     % Plot the time-frequency spectrogram
%    figure;
%     spectrogram(chirp_with_noise, hann(256), 250, 1024, fs, 'centered', 'yaxis');
%     title(sprintf('LFM Jamming Signal (Time Frequency Spectrogram) - Random Noise Level %.2f', noise_amplitude));
%     set(gca, 'YDir','reverse');
end
