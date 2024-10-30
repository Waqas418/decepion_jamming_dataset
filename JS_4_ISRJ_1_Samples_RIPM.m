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
start_time = 25e-6;              % Start time for the chirp waveform

% Sampling frequency
fs = 55e6;  % Adjust as needed

% Number of samples
num_samples = 2;

for sample = 1:num_samples
    disp(['Sample ', num2str(sample), ':']);
    
    % Parameters for TS(t) pulse train
    Ts = 5e-6; % Pulse width in seconds
    PRF_TS = 5e-6; % Pulse repetition frequency in seconds
    A = 1; % Amplitude of TS(t) pulse train

    % Parameters for TW(t) pulse train
    Tw = 2e-6; % Pulse width in seconds
    PRF_TW = 2e-6; % Pulse repetition frequency in seconds
    B = 1; % Amplitude of TW(t) and TW1(t) pulse trains

    % Number of forwards
    num_forwards = 1;

    % Slice period
    slice_period = (10e-6 - 1e-6) * rand + 1e-6; % in microseconds

    %  select the jamming duty cycle between 50% 
    jamming_duty_cycle = 50; %

    % Calculate the duration of the jamming pulse
    jamming_pulse_duration = slice_period * jamming_duty_cycle / 100;

    % Generate original chirp waveform with pulse duration
    t = linspace(0, total_duration, total_duration * fs); % Time vector
    chirp_waveform_original = chirp(t-start_time, start_frequency, pulse_duration, end_frequency, 'linear', 'complex');

    % Ensure original chirp waveform is zero before 0 microseconds and after 20 microseconds
    chirp_waveform_original(t < start_time) = 0;
    chirp_waveform_original(t > (start_time+pulse_duration)) = 0;

    % Initialize J(t) to accumulate ISRJ signals
    J_signal_accumulated = zeros(size(t));

    % Main loop for iterating over forwards
    for forward = 1:num_forwards
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

    % % Plot J(t) and its spectrogram
    % figure;
    % plot(t * 1e6, real(J_signal_accumulated));
    % xlabel('Time (\mu s)');
    % ylabel('Normalized Amplitude');
    % title('Accumulated ISRJ Time Domain Waveform');

    % figure;
    % spectrogram(J_signal_accumulated, hann(256), 250, 1024, fs, 'centered', 'yaxis');
    % set(gca, 'YDir', 'reverse');
    % title('Accumulated ISRJ Time Frequency Spectrogram');

    % Calculate the modulus (magnitude) and phase of the complex chirp waveform
    modulus_chirp = abs(J_signal_accumulated);
    phase_chirp = angle(J_signal_accumulated);

    % Plot the spectrogram of the real part of the composite signal
    figure;
    spectrogram(real(J_signal_accumulated), hann(256), 250, 1024, fs, 'centered', 'yaxis');
    title(sprintf('ISRJ1 (Spectrogram - Real Part) - Sample %d', sample));
    set(gca, 'YDir', 'reverse');
    % Save the spectrogram image
    filename = fullfile(folderPath, sprintf('ISRJ1_Spectrogram_Real_Part_Sample_%d.png', sample));
    saveas(gcf, filename);

    % % Plot the spectrogram of the imaginary part of the composite signal
    % figure;
    % spectrogram(imag(J_signal_accumulated), hann(256), 250, 1024, fs, 'centered', 'yaxis');
    % title(sprintf('ISRJ1 (Spectrogram - Imaginary Part) - Sample %d', sample));
    % set(gca, 'YDir', 'reverse');
    % % Save the spectrogram image
    % filename = fullfile(folderPath, sprintf('ISRJ1_Spectrogram_Imaginary_Part_Sample_%d.png', sample));
    % saveas(gcf, filename);
    % 
    % % Plot the spectrogram of the modulus (magnitude) of the composite signal
    % figure;
    % spectrogram(modulus_chirp, hann(256), 250, 1024, fs, 'centered', 'yaxis');
    % title(sprintf('ISRJ1 (Spectrogram - Modulus) - Sample %d', sample));
    % set(gca, 'YDir', 'reverse');
    % % Save the spectrogram image
    % filename = fullfile(folderPath, sprintf('ISRJ1_Spectrogram_Modulus_Sample_%d.png', sample));
    % saveas(gcf, filename);
    % 
    % % Plot the spectrogram of the phase of the composite signal
    % figure;
    % spectrogram(phase_chirp, hann(256), 250, 1024, fs, 'centered', 'yaxis');
    % title(sprintf('ISRJ1 (Spectrogram - Phase) - Sample %d', sample));
    % set(gca, 'YDir', 'reverse');
    % % Save the spectrogram image
    % filename = fullfile(folderPath, sprintf('ISRJ1_Spectrogram_Phase_Sample_%d.png', sample));
    % saveas(gcf, filename);
end
