clc;
clear all;

% Parameters for the original waveform
start_frequency = -6e6;        % Start frequency in Hz
end_frequency = 6e6;           % End frequency in Hz
pulse_duration = 20e-6;        % Pulse duration in seconds
total_duration = 100e-6;       % Total duration to visualize

% Sampling frequency
fs = 55e6;                    

% Number of samples for each type of signal
num_samples = 1500;

% LFM
for sample = 1:num_samples
    % Time vector
    t = linspace(0, total_duration, total_duration * fs);
    
    % Generate lfm waveform with pulse duration of 20 microseconds
    lfm = chirp(t - 55e-6, start_frequency, pulse_duration, end_frequency, 'linear', 'complex');
    
    % Ensure original chirp waveform is zero before 0 microseconds and after 20 microseconds
    lfm(t < 55e-6) = 0;
    lfm(t > (55e-6 + pulse_duration)) = 0;
    
    % Generate pure noise with random amplitude between 0.25 and 1 
    noise_amplitude = 0.25 + (1 - 0.25) * rand();
    noise = (rand(size(t)) - 0.5) * 2 * noise_amplitude;
    
    % Add noise to the origional waveform
    rad_wavform = lfm + noise;

    %Plot the time-frequency spectrogram for LFM
    figure;
    spectrogram(rad_wavform, hann(256), 250, 1024, fs, 'centered', 'yaxis');
    title(sprintf('LFM Signal - Sample %d - Random Noise Level %.2f', sample, noise_amplitude));
    set(gca, 'YDir','reverse');

    % Save the image enter the path to save samples  for example(C:\Users\waq\Desktop\Saved Images\New folder\)
    saveas(gcf, ['your path to save the samples ', sprintf('LFM_%d.png', sample)]);
end


% DDJ
for sample = 1:num_samples
    % Time vector
    t = linspace(0, total_duration, total_duration * fs);

    % Generate distance deception jamming  with  additional delay 
    ddj = chirp(t - 52e-6 - rand() * 9e-6, start_frequency, pulse_duration, end_frequency, 'linear','complex');

    % Ensure waveform is zero before 0 microseconds and after 20 microseconds for DDJ
    ddj(t < 52e-6) = 0;
    ddj(t > (52e-6 + pulse_duration)) = 0;

    % Plot the range-frequency spectrogram for DDJ
    figure;
    spectrogram(ddj, hann(256), 250, 1024, fs,'centered', 'yaxis');
    title(sprintf('Distance Deception Jamming (Time Frequency Spectrogram) - Sample %d', sample));
    set(gca, 'YDir','reverse');

    % Save the image
    saveas(gcf, ['C:\Users\waq\Desktop\Saved Images\New folder\', sprintf('DDJ_%d.png', sample)]);
end


% DFTJ
for sample = 1:num_samples
    % Generate DFTJ Jamming waveform with pulse duration and random delay for the first pulse
    delay_first_pulse = randi([1, 10]) * 1e-6;
    start_time_original = 25e-6 + delay_first_pulse;
    dftj = chirp(t - start_time_original, start_frequency, pulse_duration, end_frequency, 'linear', 'complex');

    % Ensure  waveform is zero before 0 microseconds and after 20 microseconds
    dftj(t < start_time_original) = 0;
    dftj(t > (start_time_original + pulse_duration)) = 0;

    % Generate additional  waveforms at specified intervals with the same delay
    start_times_additional = [25, 30, 35, 40, 45] * 1e-6 + delay_first_pulse;
    dftj_additional = zeros(length(start_times_additional), length(t));
    for i = 1:length(start_times_additional)
        dftj_additional(i, :) = chirp(t - start_times_additional(i), start_frequency, pulse_duration, end_frequency, 'linear', 'complex');
        dftj_additional(i, t < start_times_additional(i)) = 0;
        dftj_additional(i, t > (start_times_additional(i) + pulse_duration)) = 0;
    end

    % Combine all dftj waveforms
    composite_signal = dftj + sum(dftj_additional, 1);

    % Plot spectrum of DFTJ
    figure;
    spectrogram(composite_signal, hann(256), 250, 1024, fs, 'centered', 'yaxis');
    title(sprintf('Sample %d: Dense False Target Jamming (Time Frequency Waveform)', sample));
    set(gca, 'YDir','reverse');

    % Save the image
    saveas(gcf, ['C:\Users\waq\Desktop\Saved Images\New folder\', sprintf('DFTJ_%d.png', sample)]);
end

% ISRJ1 
for sample = 1:num_samples
    % Parameters for ISRJ1
    Ts = 5e-6;                    % Pulse width in seconds
    PRF_TS = 5e-6;                % Pulse repetition frequency in seconds
    Tw = 2e-6;                    % Pulse width in seconds
    PRF_TW = 2e-6;                % Pulse repetition frequency in seconds
    slice_period_1 = (10e-6 - 1e-6) * rand + 1e-6; % in microseconds
    jamming_duty_cycle_3 = 50;      % Jamming duty cycle in percent

    % Calculate the duration of the jamming pulse
    jamming_pulse_duration_2 = slice_period_1 * jamming_duty_cycle_3 / 100;

    % Time vector
    t = linspace(0, total_duration, total_duration * fs);

    % Generate ISRJ1  waveform with pulse duration
    ISRJ  = chirp(t - 25e-6, start_frequency, pulse_duration, end_frequency, 'linear', 'complex');

    % Ensure  waveform is zero before 0 microseconds and after 20 microseconds
    ISRJ(t < 25e-6) = 0;
    ISRJ(t > (25e-6 + pulse_duration)) = 0;

    % Initialize J(t) to accumulate ISRJ signals
    J_signal_accumulated_1 = zeros(size(t));

    % Main loop for iterating over forwards
    for forward = 1:num_samples
        % Generate the TS(t) pulse train
        TS_pulse_train = zeros(size(t));
        TS_pulse_train(mod(t, PRF_TS) <= Ts) = 1;

        % Generate the TW(t) pulse train
        TW_pulse_train = zeros(size(t));
        TW_pulse_train(mod(t, PRF_TW) <= Tw) = 1;

        % Generate the delayed TW1(t) pulse train
        delayed_TW_pulse_train_9 = zeros(size(t));
        for i = 1:length(t)
            % Check if the current time is within the TW(t) pulse width
            if mod(t(i), PRF_TW) <= Tw
                % Check if it's within the jamming pulse duration
                if mod(t(i), slice_period_1) <= jamming_pulse_duration_2
                    delayed_TW_pulse_train_9(i) = 1;
                else
                    delayed_TW_pulse_train_9(i) = 0;
                end
            end
        end

        % Multiply TS(t) with TW(t) and TW1(t)
        Pulse_train_multiplied_delayed = TS_pulse_train .* delayed_TW_pulse_train_9;

        % Generate J(t)
        J_signal = Pulse_train_multiplied_delayed .* ISRJ;

        % Accumulate ISRJ1 signals
        J_signal_accumulated_1 = J_signal_accumulated_1 + J_signal;
    end

    % Plot spectrogram ISRJ1
    figure;
    spectrogram(J_signal_accumulated_1, hann(256), 250, 1024, fs, 'centered', 'yaxis');
    set(gca, 'YDir', 'reverse');
    title(sprintf('Sample %d: ISRJ1 (Time Frequency Waveform)', sample));

    % Save the samples
    saveas(gcf, ['C:\Users\waq\Desktop\Saved Images\New folder\', sprintf('ISRJ1_%d.png', sample)]);
end

% ISRJ2 
for sample = 1:num_samples
    % Parameters for ISRJ2
    num_forwards_3 = 3;
    slice_period_2= (20e-6 - 1e-6) * rand + 1e-6; % in microseconds
    jamming_duty_cycle_3 = 75;      % Jamming duty cycle in percent

    % Calculate the duration of the jamming pulse
    jamming_pulse_duration_2 = slice_period_2 * jamming_duty_cycle_3 / 100;

    % Time vector
    t = linspace(0, total_duration, total_duration * fs);

    % % Generate  isrj2 waveform with pulse duration
    % chirp_waveform_original = chirp(t - 25e-6, start_frequency, pulse_duration, end_frequency, 'linear', 'complex');
    % 
    % % Ensure original chirp waveform is zero before 0 microseconds and after 20 microseconds
    % chirp_waveform_original(t < 25e-6) = 0;
    % chirp_waveform_original(t > (25e-6 + pulse_duration)) = 0;

    % Initialize J(t) to accumulate ISRJ signals
    J_signal_accumulated_1 = zeros(size(t));

    % Main loop for iterating over forwards
    for forward = 1:num_forwards_3
        % Generate the TS(t) pulse train
        TS_pulse_train = zeros(size(t));
        TS_pulse_train(mod(t, 5e-6) <= 5e-6) = 1;

        % Generate the TW(t) pulse train
        TW_pulse_train = zeros(size(t));
        TW_pulse_train(mod(t, 2e-6) <= 2e-6) = 1;

        % Generate the delayed TW1(t) pulse train
        delayed_TW_pulse_train_9 = zeros(size(t));
        for i = 1:length(t)
            % Check if the current time is within the TW(t) pulse width
            if mod(t(i), 2e-6) <= 2e-6
                % Check if it's within the jamming pulse duration
                if mod(t(i), slice_period_2) <= jamming_pulse_duration_2
                    delayed_TW_pulse_train_9(i) = 1;
                else
                    delayed_TW_pulse_train_9(i) = 0;
                end
            end
        end

        % Multiply TS(t) with TW(t) and TW1(t)
        Pulse_train_multiplied_delayed = TS_pulse_train .* delayed_TW_pulse_train_9;

        % Generate J(t)
        J_signal = Pulse_train_multiplied_delayed .*  ISRJ;

        % Accumulate ISRJ signals
        J_signal_accumulated_1 = J_signal_accumulated_1 + J_signal;
    end

    % Plot spectrogram ISRJ2
    figure;
    spectrogram(J_signal_accumulated_1, hann(256), 250, 1024, fs, 'centered', 'yaxis');
    set(gca, 'YDir', 'reverse');
    title(sprintf('Sample %d: ISRJ2 (Time Frequency Waveform)', sample));

    % Save the samples
    saveas(gcf, ['C:\Users\waq\Desktop\Saved Images\New folder\', sprintf('ISRJ2_%d.png', sample)]);
end


% ISRJ3
for sample = 1:num_samples
    % Parameters for ISRJ3 simulation
    Ts = 5e-6;
    PRF_TS = 5e-6;
    A = 1;
    Tw = 2e-6;
    PRF_TW = 2e-6;
    B = 1;
    start_time_original = 25e-6;
    num_forwards_3 = 6;
    slice_period_3 = (14e-6 - 1e-6) * rand + 1e-6;
    jamming_duty_cycle_3 = 86;
    jamming_pulse_duration_3 = slice_period_3 * jamming_duty_cycle_3 / 100;

    t = linspace(0, total_duration, total_duration * fs);
    ISRJ3 = chirp(t-start_time_original, start_frequency, pulse_duration, end_frequency, 'linear', 'complex');
    ISRJ3(t < start_time_original) = 0;
    ISRJ3(t > (start_time_original+pulse_duration)) = 0;

    J_signal_accumulated = zeros(size(t));

    for forward = 1:num_forwards_3
        TS_pulse_train = zeros(size(t));
        TS_pulse_train(mod(t, PRF_TS) <= Ts) = A;
        TW_pulse_train = zeros(size(t));
        TW_pulse_train(mod(t, PRF_TW) <= Tw) = B;
        delayed_TW_pulse_train_9 = zeros(size(t));
        for i = 1:length(t)
            if mod(t(i), PRF_TW) <= Tw
                if mod(t(i), slice_period_3) <= jamming_pulse_duration_3
                    delayed_TW_pulse_train_9(i) = B;
                else
                    delayed_TW_pulse_train_9(i) = 0;
                end
            end
        end
        Pulse_train_multiplied = TS_pulse_train .* TW_pulse_train;
        Pulse_train_multiplied_delayed = TS_pulse_train .* delayed_TW_pulse_train_9;
        J_signal = Pulse_train_multiplied_delayed .* ISRJ3;
        J_signal_accumulated = J_signal_accumulated + J_signal;
    end

    figure;
    spectrogram(J_signal_accumulated, hann(256), 250, 1024, fs, 'centered', 'yaxis');
    set(gca, 'YDir', 'reverse');
    title(sprintf('Sample %d: ISRJ3 (Time Frequency Waveform)', sample));
    saveas(gcf, ['C:\Users\waq\Desktop\Saved Images\New folder\', sprintf('ISRJ3_%d.png', sample)]);
end


% ISRJ4


for sample = 1:num_samples
    % Parameters for ISRJ4 simulation
    Ts = 5e-6;
    PRF_TS = 5e-6;
    A = 1;
    Tw = 2e-6;
    PRF_TW = 2e-6;
    B = 1;
    start_time_original = 25e-6;
    num_forwards_4 = 10;
    slice_period_4 = (22e-6 - 1e-6) * rand + 1e-6;
    jamming_duty_cycle_4 = 91;
    jamming_pulse_duration_4= slice_period_4 * jamming_duty_cycle_4 / 100;

    t = linspace(0, total_duration, total_duration * fs);
    ISRJ4 = chirp(t-start_time_original, start_frequency, pulse_duration, end_frequency, 'linear', 'complex');
    ISRJ4(t < start_time_original) = 0;
    ISRJ4(t > (start_time_original+pulse_duration)) = 0;

    J_signal_accumulated = zeros(size(t));

    for forward = 1:num_forwards_4
        TS_pulse_train = zeros(size(t));
        TS_pulse_train(mod(t, PRF_TS) <= Ts) = A;
        TW_pulse_train = zeros(size(t));
        TW_pulse_train(mod(t, PRF_TW) <= Tw) = B;
        delayed_TW_pulse_train_9 = zeros(size(t));
        for i = 1:length(t)
            if mod(t(i), PRF_TW) <= Tw
                if mod(t(i), slice_period_4) <= jamming_pulse_duration_4
                    delayed_TW_pulse_train_9(i) = B;
                else
                    delayed_TW_pulse_train_9(i) = 0;
                end
            end
        end
        Pulse_train_multiplied = TS_pulse_train .* TW_pulse_train;
        Pulse_train_multiplied_delayed = TS_pulse_train .* delayed_TW_pulse_train_9;
        J_signal = Pulse_train_multiplied_delayed .* ISRJ4;
        J_signal_accumulated = J_signal_accumulated + J_signal;
    end

    figure;
    spectrogram(J_signal_accumulated, hann(256), 250, 1024, fs, 'centered', 'yaxis');
    set(gca, 'YDir', 'reverse');
    title(sprintf('Sample %d: ISRJ4 (Time Frequency Waveform)', sample));
    saveas(gcf, ['C:\Users\waq\Desktop\Saved Images\New folder\', sprintf('ISRJ4_%d.png', sample)]);
end

% ISRJ5 

for sample = 1:num_samples
    % Parameters for ISRJ5 simulation
    Ts = 5e-6;
    PRF_TS = 5e-6;
    A = 1;
    Tw = 2e-6;
    PRF_TW = 2e-6;
    B = 1;
    start_time_original = 25e-6;
    num_forwards_5 = 20;
    slice_period_5 = (40e-6 - 1e-6) * rand + 1e-6;
    jamming_duty_cycle_5 = 95;
    jamming_pulse_duration_5 = slice_period_5 * jamming_duty_cycle_5 / 100;

    t = linspace(0, total_duration, total_duration * fs);
    ISRJ5 = chirp(t-start_time_original, start_frequency, pulse_duration, end_frequency, 'linear', 'complex');
    ISRJ5(t < start_time_original) = 0;
    ISRJ5(t > (start_time_original+pulse_duration)) = 0;

    J_signal_accumulated = zeros(size(t));

    for forward = 1:num_forwards_5
        TS_pulse_train = zeros(size(t));
        TS_pulse_train(mod(t, PRF_TS) <= Ts) = A;
        TW_pulse_train = zeros(size(t));
        TW_pulse_train(mod(t, PRF_TW) <= Tw) = B;
        delayed_TW_pulse_train_9 = zeros(size(t));
        for i = 1:length(t)
            if mod(t(i), PRF_TW) <= Tw
                if mod(t(i), slice_period_5) <= jamming_pulse_duration_5
                    delayed_TW_pulse_train_9(i) = B;
                else
                    delayed_TW_pulse_train_9(i) = 0;
                end
            end
        end
        Pulse_train_multiplied = TS_pulse_train .* TW_pulse_train;
        Pulse_train_multiplied_delayed = TS_pulse_train .* delayed_TW_pulse_train_9;
        J_signal = Pulse_train_multiplied_delayed .* ISRJ5;
        J_signal_accumulated = J_signal_accumulated + J_signal;
    end

    figure;
    spectrogram(J_signal_accumulated, hann(256), 250, 1024, fs, 'centered', 'yaxis');
    set(gca, 'YDir', 'reverse');
    title(sprintf('Sample %d: ISRJ5 (Time Frequency Waveform)', sample));
    saveas(gcf, ['C:\Users\waq\Desktop\Saved Images\New folder\', sprintf('ISRJ5_%d.png', sample)]);
end


% SNJ
for sample = 1:num_samples
    % Parameters for SNJ
    Ts = 5e-6;                    % Pulse width in seconds
    PRF_TS = 5e-6;                % Pulse repetition frequency in seconds
    Tw = 2e-6;                    % Pulse width in seconds
    PRF_TW = 2e-6;                % Pulse repetition frequency in seconds
    start_time = 60e-6;           % Start time for the chirp waveform
    num_forwards_6 = randi([1, 4]); % Random number of forwards
    slice_period_SNJ = (10e-6 - 5e-6) * rand + 5e-6; % Random slice period between 5 and 10 microseconds
    jamming_duty_cycle_SNJ = (50 - 20) * rand + 20; % Random duty cycle between 20% and 50%
    jamming_pulse_duration_SNJ = slice_period_SNJ * jamming_duty_cycle_SNJ / 100;

    % Time vector
    t = linspace(0, total_duration, total_duration * fs);

    % Generate SNJ waveform with pulse duration
    SNJ = chirp(t - start_time, start_frequency, pulse_duration, end_frequency, 'linear', 'complex');

    % Ensure original chirp waveform is zero before 0 microseconds and after 20 microseconds
    SNJ(t < start_time) = 0;
    SNJ(t > (start_time + pulse_duration)) = 0;

    % Initialize J(t) to accumulate SNJ signals
    J_signal_accumulated = zeros(size(t));

    % Main loop for iterating over forwards
    for forward = 1:num_forwards_6
        % Generate the TS(t) pulse train
        TS_pulse_train = zeros(size(t));
        TS_pulse_train(mod(t, PRF_TS) <= Ts) = 1;

        % Generate the TW(t) pulse train
        TW_pulse_train = zeros(size(t));
        TW_pulse_train(mod(t, PRF_TW) <= Tw) = 1;

        % Generate the delayed TW1(t) pulse train
        delayed_TW_pulse_train_9 = zeros(size(t));
        for i = 1:length(t)
            % Check if the current time is within the TW(t) pulse width
            if mod(t(i), PRF_TW) <= Tw
                % Check if it's within the jamming pulse duration
                if mod(t(i), slice_period_SNJ) <= jamming_pulse_duration_SNJ
                    delayed_TW_pulse_train_9(i) = 1;
                else
                    delayed_TW_pulse_train_9(i) = 0;
                end
            end
        end

        % Multiply TS(t) with TW(t) and TW1(t)
        Pulse_train_multiplied_delayed = TS_pulse_train .* delayed_TW_pulse_train_9;

        % Generate J(t)
        J_signal = Pulse_train_multiplied_delayed .* SNJ;

        % Accumulate SNJ signals
        J_signal_accumulated = J_signal_accumulated + J_signal;
    end

    % Generate white Gaussian noise
    sigma_squared = 0.1; % Variance of the white Gaussian noise
    white_gaussian_noise = sqrt(sigma_squared) * randn(size(t)); % Scale by square root of variance

    % Addition of Gaussian noise
    J_with_noise = 0.0005 * (J_signal_accumulated .* white_gaussian_noise);

    % Plot spectrogram SNJ
    figure;
    spectrogram(J_with_noise, hann(256), 250, 1024, fs, 'centered', 'yaxis');
    set(gca, 'YDir', 'reverse');
    title(sprintf('Sample %d: SNJ (Time Frequency Waveform)', sample));
    saveas(gcf, ['C:\Users\waq\Desktop\Saved Images\New folder\', sprintf('SNJ_%d.png', sample)]);
end
%DDJ_SNJ
% Number of samples to generate
num_samples = 3;

% Loop to generate samples
for sample = 1:num_samples
    % Random delay in the range of 1 to 10 microseconds
    delay_additional = (10 - 1) * rand() * 1e-6 + 1e-6;
    
    % Random slice period between 5 to 10 microseconds
    slice_period = (10e-6 - 5e-6) * rand() + 5e-6;
    
    % Time vector
    t = linspace(0, total_duration, total_duration * fs);

    % Generate original chirp waveform with pulse duration and additional delay
    dftj_snj = chirp(t - start_time - delay_additional, start_frequency, pulse_duration, end_frequency, 'linear','complex');
    
    % Ensure original chirp waveform is zero before 0 microseconds and after 20 microseconds
    dftj_snj(t < (start_time + delay_additional)) = 0;
    dftj_snj(t > (start_time + delay_additional + pulse_duration)) = 0;

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
    delayed_TW_pulse_train(mod(t, PRF_TW) <= Tw & mod(t, PRF_TS) > Tw) = B; % Jamming condition

    % Multiply TS(t) with TW(t) and TW1(t)
    Pulse_train_multiplied_delayed = TS_pulse_train .* delayed_TW_pulse_train;

    % Generate J(t)
    J_signal = Pulse_train_multiplied_delayed .* dftj_snj;

    % Generate white Gaussian noise
    sigma_squared = 0.1; % Variance of the white Gaussian noise
    white_gaussian_noise = (sqrt(sigma_squared) * randn(size(t))); % Scale by square root of variance

    % Multiply J(t) with white Gaussian noise
    J_with_noise = 0.0005*(J_signal .* white_gaussian_noise);

    % Combine the DDJ and SNJ signals
    J_combined = (0.00005*dftj_snj) + J_with_noise;

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
        TW_pulse_train(mod(t, PRF_TW) <= Tw & mod(t, PRF_TS) > Tw) = B; % Jamming condition

        % Generate the delayed TW1(t) pulse train
        delayed_TW_pulse_train = zeros(size(t));
        delayed_TW_pulse_train(mod(t, PRF_TW) <= Tw & mod(t, slice_period) <= jamming_pulse_duration) = B; % Jamming condition

        % Multiply TS(t) with TW(t) and TW1(t)
        Pulse_train_multiplied_delayed = TS_pulse_train .* delayed_TW_pulse_train;

        % Generate J(t)
        J_signal = Pulse_train_multiplied_delayed .* dftj_snj;

        % Accumulate ISRJ signals
        J_signal_accumulated = J_signal_accumulated + J_signal;
    end

    % Multiply J(t) with white Gaussian noise
    J_with_noise = 0.0005*(J_signal_accumulated .* white_gaussian_noise);

    % Combine the DDJ and SNJ signals
    J_combined = J_combined + J_with_noise;

    % Plot spectrogram
    figure;
    spectrogram(J_combined, hann(256), 250, 1024, fs, 'centered', 'yaxis');
    set(gca, 'YDir', 'reverse');
    title(sprintf('Sample %d: (DDJ+SNJ) (Time Frequency Waveform)', sample));
    saveas(gcf, ['C:\Users\waq\Desktop\Saved Images\New folder\', sprintf('DDJ + SNJ_%d.png',  sample)]);
end

%DFTJ_SNJ

% Parameters for TS(t) pulse train
Ts = 15e-6; % Pulse width in seconds
PRF_TS = 15e-6; % Pulse repetition frequency in seconds
A = 1; % Amplitude of TS(t) pulse train

% Parameters for TW(t) pulse train
Tw = 1.25e-6; % Pulse width in seconds
PRF_TW = 1.25e-6; % Pulse repetition frequency in seconds
B = 1; % Amplitude of TW(t) and TW1(t) pulse trains

% Parameters for the original chirp waveform

start_time = 5e-6;          

% Time vector
t = linspace(0, total_duration, total_duration * fs);

% Generate original chirp waveform with pulse duration
dftj_snj = chirp(t-start_time, start_frequency, pulse_duration, end_frequency, 'linear', 'complex');

% Ensure original chirp waveform is zero before 0 microseconds and after 20 microseconds
dftj_snj(t < start_time) = 0;
dftj_snj(t > (start_time+pulse_duration)) = 0;

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
J_signal = Pulse_train_multiplied_delayed .* dftj_snj;

% Generate white Gaussian noise
sigma_squared = 0.1; % Variance of the white Gaussian noise
white_gaussian_noise = sqrt(sigma_squared) * randn(size(t)); % Scale by square root of variance

% Multiply J(t) with white Gaussian noise
J_with_noise = 0.0005 * (J_signal .* white_gaussian_noise);

% Generate additional chirp waveforms at randomly selected intervals
additional_start_times = [20, 25, 30] * 1e-6;  % Start times for additional chirps
num_samples = 3;
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
    combined_signals{sample_idx} = (J_with_noise + 0.00005 * (dftj_snj + sum(additional_chirp_samples, 1)));
end

% Plot J(t) and its spectrogram for each combined signal
for sample_idx = 1:num_samples
    % figure;
    % 
    % plot(t * 1e6, real(combined_signals{sample_idx}));
    % xlabel('Time (\mu s)');
    % ylabel('Normalized Amplitude');
    % title(['(DFTJ+SNJ) Time Domain Waveform - Sample ', num2str(sample_idx)]);

    figure;
    spectrogram(combined_signals{sample_idx}, hann(256), 250, 1024, fs, 'centered', 'yaxis');
    set(gca, 'YDir', 'reverse');

     title(sprintf('Sample %d: (DFTJ+SNJ) (Time Frequency Waveform)', sample_idx));
        saveas(gcf, ['C:\Users\waq\Desktop\Saved Images\New folder\', sprintf('DFTJ + SNJ_%d.png',  sample_idx)]);
end
%DDJ_ISRJ

% Parameters for  waveform
start_frequency = -6e6;          % Start frequency in Hz
end_frequency = 6e6;             % End frequency in Hz
pulse_duration = 20e-6;          % Pulse duration in seconds
total_duration = 100e-6;         % Total duration to visualize
start_time = 5e-6;               % Start time for the chirp waveform
delay_additional = 10e-6;        % Additional delay in seconds

% Time vector for (DDJ_ISRJ) waveform
t = linspace(0, total_duration, total_duration * fs);

% Generate  waveform with pulse duration and additional delay
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

    % Plot and save the spectrogram
    figure;
    spectrogram(J_combined, hann(256), 250, 1024, fs, 'centered', 'yaxis');
    set(gca, 'YDir', 'reverse');
     
    title(sprintf('Sample %d: (DDJ+ISRJ) (Time Frequency Waveform)', sample));
    saveas(gcf, ['C:\Users\waq\Desktop\Saved Images\New folder\', sprintf('DDJ + ISRJ_%d.png',  sample)]);
end

      
   
