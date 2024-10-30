clc;
clear all;

% Define the folder path
folderPath = 'L:\JS_1_LFM_JAMMING_Samples_RIPM';

% Create the folder if it doesn't exist
if ~isfolder(folderPath)
    mkdir(folderPath);
end

% Parameters for the original chirp waveform
start_frequency = -6e6;          % Start frequency in Hz
end_frequency = 6e6;             % End frequency in Hz
pulse_duration = 25e-6;          % Pulse duration in seconds
total_duration = 100e-6;         % Total duration to visualize
start_time = 10e-6;              % Start time for the chirp waveform

% Sampling frequency
fs = 100e6;  % Adjust as needed

% Number of samples to generate
num_samples = 2;

% Loop to generate samples
for sample = 1:num_samples
    % Random delay in the range of 1 to 10 microseconds
    delay_additional = (10 - 1) * rand() * 1e-6 + 1e-6;
    
    % Random slice period between 5 to 10 microseconds
    slice_period = (10e-6 - 5e-6) * rand() + 5e-6;
    
    % Time vector
    t = linspace(0, total_duration, total_duration * fs);

    % Generate original chirp waveform with pulse duration and additional delay
    chirp_waveform_original = chirp(t - start_time - delay_additional, start_frequency, pulse_duration, end_frequency, 'linear','complex');
    
    % Ensure original chirp waveform is zero before 0 microseconds and after 20 microseconds
    chirp_waveform_original(t < (start_time + delay_additional)) = 0;
    chirp_waveform_original(t > (start_time + delay_additional + pulse_duration)) = 0;

    % Parameters for TS(t) pulse train
    Ts = 5e-6; % Pulse width in seconds
    PRF_TS = 5e-6; % Pulse repetition frequency in seconds
    A = 1; % Amplitude of TS(t) pulse train

    % Parameters for TW(t) pulse train
    Tw = 1.25e-6; % Pulse width in seconds
    PRF_TW = 1.25e-6; % Pulse repetition frequency in seconds
    B = 1; % Amplitude of TW(t) and TW1(t) pulse trains

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
    J_signal = (Pulse_train_multiplied_delayed .* chirp_waveform_original);

    % Generate white Gaussian noise
    sigma_squared = 0.1; % Variance of the white Gaussian noise
    white_gaussian_noise = 1 * (sqrt(sigma_squared) * randn(size(t))); % Scale by square root of variance

    % Multiply J(t) with white Gaussian noise
    J_with_noise = 0.0005*(J_signal .* white_gaussian_noise);

    % Combine the DDJ and SNJ signals
    J_combined = (0.00005*chirp_waveform_original) + J_with_noise;

    % Random number of forwards between 1 to 4
    num_forwards = randi([1, 4]);

    % Random jamming duty cycle between 20 to 50
    jamming_duty_cycle = (50 - 20) * rand + 20;

    % Generate the jamming pulse
    jamming_pulse_duration = slice_period * jamming_duty_cycle / 100;

    % Initialize J(t) to accumulate ISRJ signals
    J_signal_accumulated = zeros(size(t));

    % Main loop for iterating over forwards
    for forward = 1:num_forwards
        % Generate the TW(t) pulse train
        TW_pulse_train = zeros(size(t));
        TW_pulse_train(mod(t, PRF_TW) <= Tw) = B; % Generate rectangular pulses

        % Generate the delayed TW1(t) pulse train
        delayed_TW_pulse_train = zeros(size(t));
        for i = 1:length(t)
            % Check if the current time is within the TW(t) pulse width
            if mod(t(i), PRF_TW) <= Tw
                % Check if it's within the jamming pulse duration
                if mod(t(i), slice_period) <= jamming_pulse_duration
                    delayed_TW_pulse_train(i) = B; % Set to amplitude B (jamming)
                else
                    delayed_TW_pulse_train(i) = 0; % Set to zero (idle)
                end
            end
        end

        % Multiply TS(t) with TW(t) and TW1(t)
        Pulse_train_multiplied = TS_pulse_train .* TW_pulse_train;
        Pulse_train_multiplied_delayed = TS_pulse_train .* delayed_TW_pulse_train;

        % Multiply P(t) with X(t) to get J(t)
        J_signal = Pulse_train_multiplied_delayed .* chirp_waveform_original;

        % Accumulate ISRJ signals
        J_signal_accumulated = J_signal_accumulated + J_signal;
    end

    % Multiply J(t) with white Gaussian noise
    J_with_noise = 0.0005*(J_signal_accumulated .* white_gaussian_noise);

    % Combine the DDJ and SNJ signals
    J_combined = J_combined + J_with_noise;

    % % Plot J(t) and its spectrogram
    % figure;
    % plot(t * 1e6, real(J_combined));
    % xlabel('Time (\mu s)');
    % ylabel('Normalized Amplitude');
    % title(['(DDJ+SNJ) Time Domain Waveform - Sample ', num2str(sample)]);

    % figure;
    % spectrogram(J_combined, hann(256), 250, 1024, fs, 'centered', 'yaxis');
    % set(gca, 'YDir', 'reverse');
    % title(['(DDJ+SNJ) Time Frequency Spectrogram - Sample ', num2str(sample)]);


    % Calculate the modulus (magnitude) and phase of the complex chirp waveform
    modulus_chirp = abs(J_combined);
    phase_chirp = angle(J_combined);

    % Plot the spectrogram of the real part of the original chirp waveform with the additional delay
    figure;
    spectrogram(real(J_combined), hann(256), 250, 1024, fs, 'centered', 'yaxis');
    title(sprintf('(DDJ+SNJ) (Spectrogram - Real Part) - Sample %d', sample));
    set(gca, 'YDir', 'reverse');
    
    % Save the image
    saveas(gcf, fullfile(folderPath, sprintf('Sample_%d_Real_Part.png', sample)));

    % % Plot the spectrogram of the imaginary part of the original chirp waveform with the additional delay
    % figure;
    % spectrogram(imag(J_combined), hann(256), 250, 1024, fs, 'centered', 'yaxis');
    % title(sprintf('(DDJ+SNJ) (Spectrogram - Imaginary Part) - Sample %d', sample));
    % set(gca, 'YDir', 'reverse');
    % 
    % % Save the image
    % saveas(gcf, fullfile(folderPath, sprintf('Sample_%d_Imaginary_Part.png', sample)));
    % 
    % % Plot the spectrogram of the modulus (magnitude) of the original chirp waveform
    % figure;
    % spectrogram(modulus_chirp, hann(256), 250, 1024, fs, 'centered', 'yaxis');
    % title(sprintf('(DDJ+SNJ) (Spectrogram - Modulus) - Sample %d', sample));
    % set(gca, 'YDir', 'reverse');
    % 
    % % Save the image
    % saveas(gcf, fullfile(folderPath, sprintf('Sample_%d_Modulus.png', sample)));
    % 
    % % Plot the spectrogram of the phase of the original chirp waveform
    % figure;
    % spectrogram(phase_chirp, hann(256), 250, 1024, fs, 'centered', 'yaxis');
    % title(sprintf('(DDJ+SNJ) (Spectrogram - Phase) - Sample %d', sample));
    % set(gca, 'YDir', 'reverse');
    % 
    % % Save the image
    % saveas(gcf, fullfile(folderPath, sprintf('Sample_%d_Phase.png', sample)));
end
