clc;
clear;

% Parameters for TS(t) pulse train
Ts = 15e-6; % Pulse width in seconds
PRF_TS = 15e-6; % Pulse repetition frequency in seconds
A = 1; % Amplitude of TS(t) pulse train

% Parameters for TW(t) pulse train
Tw = 1.25e-6; % Pulse width in seconds
PRF_TW = 1.25e-6; % Pulse repetition frequency in seconds
B = 1; % Amplitude of TW(t) and TW1(t) pulse trains

% Parameters for the original chirp waveform
start_frequency = -6e6;          % Start frequency in Hz
end_frequency = 6e6;           % End frequency in Hz
pulse_duration = 20e-6;         % Pulse duration in seconds
total_duration = 100e-6;        % Total duration to visualize
start_time = 5e-6;             % Start time for the chirp waveform

% Sampling frequency
fs = 75e6;  % Adjust as needed

% Time vector
t = linspace(0, total_duration, total_duration * fs);

% Generate original chirp waveform with pulse duration
chirp_waveform_original = chirp(t-start_time, start_frequency, pulse_duration, end_frequency, 'linear', 'complex');

% Ensure original chirp waveform is zero before 0 microseconds and after 20 microseconds
chirp_waveform_original(t < start_time) = 0;
chirp_waveform_original(t > (start_time+pulse_duration)) = 0;

% Generate the TS(t) pulse train
TS_pulse_train = zeros(size(t));
TS_pulse_train(mod(t, PRF_TS) <= Ts) = A; % Generate rectangular pulses

% Generate the TW(t) pulse train
TW_pulse_train = zeros(size(t));
TW_pulse_train(mod(t, PRF_TW) <= Tw) = B; % Generate rectangular pulses

% Generate the delayed TW1(t) pulse train
delayed_TW_pulse_train = zeros(size(t));
for i = 1:length(t)
    % Check if the current time is within the TW(t) pulse width
    if mod(t(i), PRF_TW) <= Tw
        % Check if it's within the first Tw microseconds
        if mod(t(i), PRF_TS) <= Tw
            delayed_TW_pulse_train(i) = 0; % Set to zero
        else
            delayed_TW_pulse_train(i) = B; % Set to amplitude B
        end
    end
end

% Multiply TS(t) with TW(t) and TW1(t)
Pulse_train_multiplied = TS_pulse_train .* TW_pulse_train;
Pulse_train_multiplied_delayed = TS_pulse_train .* delayed_TW_pulse_train;

% Multiply P(t) with X(t) to get J(t)
J_signal = Pulse_train_multiplied_delayed .* chirp_waveform_original;

% Generate white Gaussian noise
sigma_squared = 0.1; % Variance of the white Gaussian noise
white_gaussian_noise = sqrt(sigma_squared) * randn(size(t)); % Scale by square root of variance

% Multiply J(t) with white Gaussian noise
J_with_noise = 0.0005 * (J_signal .* white_gaussian_noise);

% Generate additional chirp waveforms at randomly selected intervals
additional_start_times = [20, 25, 30] * 1e-6;  % Start times for additional chirps
num_samples = 2;
combined_signals = cell(num_samples, 1);
for sample_idx = 1:num_samples
    % Randomly select start times from additional_start_times
    random_indices = randi(length(additional_start_times), [1, randi(3)]);
    selected_start_times = additional_start_times(random_indices);
    
    % Generate chirp waveforms at selected start times
    additional_chirp_samples = zeros(length(selected_start_times), length(t));
    for i = 1:length(selected_start_times)
        % Generate chirp waveform with the same length as the time vector
        additional_chirp_samples(i, :) = chirp(t - selected_start_times(i), start_frequency, pulse_duration, end_frequency, 'linear', 'complex');
        % Ensure additional chirp waveform is zero before start time and after end time
        additional_chirp_samples(i, t < selected_start_times(i)) = 0;
        additional_chirp_samples(i, t > (selected_start_times(i) + pulse_duration)) = 0;
    end
    
    % Combine all chirp waveforms for each sample
    combined_signals{sample_idx} = (J_with_noise + 0.00005 * (chirp_waveform_original + sum(additional_chirp_samples, 1)));
end

% Define the folder path
folderPath = 'L:\JS_1_LFM_JAMMING_Samples_RIPM';

% Create the folder if it doesn't exist
if ~isfolder(folderPath)
    mkdir(folderPath);
end

% Plot J(t) and its spectrogram for each combined signal
for sample_idx = 1:num_samples
    % Calculate the modulus (magnitude) and phase of the complex chirp waveform
    modulus_chirp = abs(combined_signals{sample_idx});
    phase_chirp = angle(combined_signals{sample_idx});

    % Plot the spectrogram of the real part of the original chirp waveform with the additional delay
    figure;
    spectrogram(real(combined_signals{sample_idx}), hann(256), 250, 1024, fs, 'centered', 'yaxis');
    title(sprintf('(DFTJ+SNJ) (Spectrogram - Real Part) - Sample %d', sample_idx));
    set(gca, 'YDir', 'reverse');
    saveas(gcf, fullfile(folderPath, sprintf('Real_Part_Sample_%d.png', sample_idx)));

    % % Plot the spectrogram of the imaginary part of the original chirp waveform with the additional delay
    % figure;
    % spectrogram(imag(combined_signals{sample_idx}), hann(256), 250, 1024, fs, 'centered', 'yaxis');
    % title(sprintf('(DFTJ+SNJ) (Spectrogram - Imaginary Part) - Sample %d', sample_idx));
    % set(gca, 'YDir', 'reverse');
    % saveas(gcf, fullfile(folderPath, sprintf('Imaginary_Part_Sample_%d.png', sample_idx)));
    % 
    % % Plot the spectrogram of the modulus (magnitude) of the original chirp waveform
    % figure;
    % spectrogram(modulus_chirp, hann(256), 250, 1024, fs, 'centered', 'yaxis');
    % title(sprintf('(DFTJ+SNJ) (Spectrogram - Modulus) - Sample %d', sample_idx));
    % set(gca, 'YDir', 'reverse');
    % saveas(gcf, fullfile(folderPath, sprintf('Modulus_Sample_%d.png', sample_idx)));
    % 
    % % Plot the spectrogram of the phase of the original chirp waveform
    % figure;
    % spectrogram(phase_chirp, hann(256), 250, 1024, fs, 'centered', 'yaxis');
    % title(sprintf('(DFTJ+SNJ) (Spectrogram - Phase) - Sample %d', sample_idx));
    % set(gca, 'YDir', 'reverse');
    % saveas(gcf, fullfile(folderPath, sprintf('Phase_Sample_%d.png', sample_idx)));
end
