clc;
clear all;

% Define the folder path
folderPath = 'D:\SubDataSet _Saved_Samples\Extraction\JS_12_DDJ_ISRJ_Samples_RIPM\Phase';



% Sampling frequency
fs = 100e6;  % Adjust as needed

% Number of samples to generate
num_samples = 500;

% Parameters for the original chirp waveform
start_frequency = -6e6;          % Start frequency in Hz
end_frequency = 6e6;             % End frequency in Hz
pulse_duration = 20e-6;          % Pulse duration in seconds
total_duration = 100e-6;         % Total duration to visualize
start_time = 5e-6;               % Start time for the chirp waveform
delay_additional = 10e-6;        % Additional delay in seconds

% Time vector for the original chirp waveform
t = linspace(0, total_duration, total_duration * fs);

% Generate original chirp waveform with pulse duration and additional delay
chirp_waveform_original1 = chirp(t - start_time - delay_additional, start_frequency, pulse_duration, end_frequency, 'linear','complex');

% Ensure original chirp waveform is zero before 0 microseconds and after 20 microseconds
chirp_waveform_original1(t < (start_time + delay_additional)) = 0;
chirp_waveform_original1(t > (start_time + delay_additional + pulse_duration)) = 0;

% Loop to generate samples
for sample = 1:num_samples
    % Parameters for TS(t) pulse train
    Ts = 5e-6; % Pulse width in seconds
    PRF_TS = 5e-6; % Pulse repetition frequency in seconds
    A = 1; % Amplitude of TS(t) pulse train

    % Parameters for TW(t) pulse train
    Tw = 2e-6; % Pulse width in seconds
    PRF_TW = 2e-6; % Pulse repetition frequency in seconds
    B = 1; % Amplitude of TW(t) and TW1(t) pulse trains

    % Generate the TS(t) pulse train
    TS_pulse_train = zeros(size(t));
    TS_pulse_train(mod(t, PRF_TS) <= Ts) = A; % Generate rectangular pulses

    % Random slice period between 5 to 10 microseconds
    slice_period = (10e-6 - 5e-6) * rand() + 5e-6;

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
        J_signal = Pulse_train_multiplied_delayed .* chirp_waveform_original1;

        % Accumulate ISRJ signals
        J_signal_accumulated = J_signal_accumulated + J_signal;
    end

    % Combine the DDJ and ISRJ signals
    J_combined = J_signal_accumulated + chirp_waveform_original1;

    % Extract real, imaginary, phase, and modulus parts
    real_part = real(J_combined);
    imag_part = imag(J_combined);
    phase = angle(J_combined);
    modulus = abs(J_combined);

    % Plot the spectrogram of the real part of the original chirp waveform with the additional delay
    % figure;
    % spectrogram(real_part,hann(256), 250, 1024, fs, 'centered', 'yaxis');
    % title(sprintf('DDJ+ISRJ(Spectrogram - Real Part) - Sample %d', sample));
    % set(gca, 'YDir', 'reverse');
    % % Save the figure
    % saveas(gcf, fullfile(folderPath, sprintf('DDJ+ISRJ_Real_Part_Sample_%d.png', sample)));

 % % Plot the spectrogram of the imaginary part of the original chirp waveform with the additional delay
    % figure;
    % spectrogram(imag_part, hann(256), 250, 1024, fs, 'centered', 'yaxis');
    % title(sprintf('DDJ+ISRJ (Spectrogram - imaginary Part) - Sample %d', sample));
    % set(gca, 'YDir', 'reverse');
    % % Save the figure
    % saveas(gcf, fullfile(folderPath, sprintf('DDJ+ISRJ_Imaginary_Part_Sample_%d.png', sample)));

    % 
    % % Plot the spectrogram of the modulus (magnitude) of the original chirp waveform
    % figure;
    % spectrogram(modulus, hann(256), 250, 1024, fs, 'centered', 'yaxis');
    % title(sprintf('DDJ+ISRJ (Spectrogram - Modulus) - Sample %d', sample));
    % set(gca, 'YDir', 'reverse');
    % % Save the figure
    % saveas(gcf, fullfile(folderPath, sprintf('DDJ+ISRJ_Modulus_Sample_%d.png', sample)));

    % 
    % % Plot the spectrogram of the phase of the original chirp waveform
    figure;
    spectrogram(phase, hann(256), 250, 1024, fs, 'centered', 'yaxis');
    title(sprintf('DDJ+ISRJ (Spectrogram - Phase) - Sample %d', sample));
    set(gca, 'YDir', 'reverse');
    % Save the figure
    saveas(gcf, fullfile(folderPath, sprintf('DDJ+ISRJ_Phase_Sample_%d.png', sample)));
end
