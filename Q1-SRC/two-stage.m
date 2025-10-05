% Two-stage SRC example (96 kHz -> 44.1 kHz) — cleaned and commented
% - Produces a two-tone test signal
% - Performs first-stage upsample->filter->downsample (L/M)
% - Performs second-stage upsample->filter->downsample (L1/M1)


%% 0) Basic parameters and test signal
Fs = 96e3;                       % original sampling rate (Hz)
t  = 0 : 1/Fs : 0.002;           % short time vector for plotting (0..2 ms)

% Two test tones (exercise passband/aliasing behaviour)
f1 = 20e3;  % tone 1 (Hz)
f2 = 30e3;  % tone 2 (Hz)
x  = sin(2*pi*f1.*t) + sin(2*pi*f2.*t);  % composite test signal

%% -------------------- Stage 1 (coarse-ratio SRC) ------------------------
% Integer factors for stage 1
L  = 3;      % interpolation factor (upsample)
M  = 4;      % decimation factor (downsample)

% 1. Upsample by zero-stuffing (no filtering yet)
x_up = upsample(x, L);

% 2. Design anti-imaging / anti-alias lowpass filter for Stage 1
Fs_up_stage1 = Fs * L;           % sampling rate after upsampling (Hz)
Fpass = 20e3;                    % passband edge (Hz)
Fstop = 22e3;                    % stopband start (Hz)
Nyq_stage1 = Fs_up_stage1 / 2;   % Nyquist after upsampling

% Design normalized cutoff frequencies relative to Nyquist for designfilt
lpFilt = designfilt('lowpassfir', ...
    'PassbandFrequency', Fpass / Nyq_stage1, ...
    'StopbandFrequency', Fstop / Nyq_stage1, ...
    'PassbandRipple', 0.01, ...
    'StopbandAttenuation', 100, ...
    'DesignMethod', 'kaiserwin');

% Extract coefficient vector and print filter length for diagnostics
b = lpFilt.Coefficients;
fprintf('Stage 1 filter length: %d taps\n', length(b));

% 3. Filter the upsampled signal (removes imaging)
x_filt = filter(lpFilt, x_up);

% 4. Downsample by M to produce the stage-1 output
y = downsample(x_filt, M);

% 5. Effective sampling rate after stage 1
Fs_stage1 = Fs * (L / M);
fprintf('Sampling rate after Stage 1: %.2f Hz (%.3f kHz)\n', Fs_stage1, Fs_stage1/1e3);

%% -------------------- Stage 2 (fine-ratio SRC) --------------------------
% Integer factors for stage 2
L1 = 49;    % interpolation factor
M1 = 80;    % decimation factor

% 1. Upsample stage-1 output (zero-stuff)
x_up_1 = upsample(y, L1);

% 2. Design anti-imaging / anti-alias filter for Stage 2
Fs_up_stage2 = Fs_stage1 * L1;    % sampling rate after stage-2 upsampling
Nyq_stage2 = Fs_up_stage2 / 2;

% Use same pass/stop edges as the original design (relative to stage2 Nyquist)
lpFilt1 = designfilt('lowpassfir', ...
    'PassbandFrequency', Fpass / Nyq_stage2, ...
    'StopbandFrequency', Fstop / Nyq_stage2, ...
    'PassbandRipple', 0.01, ...
    'StopbandAttenuation', 100, ...
    'DesignMethod', 'kaiserwin');

b1 = lpFilt1.Coefficients;
fprintf('Stage 2 filter length: %d taps\n', length(b1));

% 3. Filter the upsampled signal for stage 2
x_filt1 = filter(lpFilt1, x_up_1);

% 4. Downsample by M1 to produce final output
y1 = downsample(x_filt1, M1);

% 5. Final sampling rate after the two-stage SRC
Fs_final = Fs_stage1 * (L1 / M1);
fprintf('Final sampling rate after Stage 2: %.2f Hz (%.3f kHz)\n', Fs_final, Fs_final/1e3);

%% ----------------------------- Plotting --------------------------------
% Short-time stem plots to visualise intermediate signals.
% Titles use the actual factor variables to avoid confusion.

subplot(5,1,1);
stem(0:100, x(1:101), 'filled', 'MarkerSize', 3);
title('Original signal (96 kHz)');
xlabel('Sample index'); ylabel('Amplitude');

subplot(5,1,2);
stem(0:443, x_up(1:444), 'filled', 'MarkerSize', 3);
title(sprintf('Stage 1 — After Upsampling (L = %d)', L));
xlabel('Sample index'); ylabel('Amplitude');

subplot(5,1,3);
stem(0:80, y(1:81), 'filled', 'MarkerSize', 3);
title(sprintf('Stage 1 — After Downsampling (M = %d), Fs = %.1f kHz', M, Fs_stage1/1e3));
xlabel('Sample index'); ylabel('Amplitude');

subplot(5,1,4);
stem(0:100, x_up_1(1:101), 'filled', 'MarkerSize', 3);
title(sprintf('Stage 2 — After Upsampling (L1 = %d)', L1));
xlabel('Sample index'); ylabel('Amplitude');

subplot(5,1,5);
stem(0:80, y1(1:81), 'filled', 'MarkerSize', 3);
title(sprintf('Stage 2 — After Downsampling (M1 = %d), Fs = %.3f kHz', M1, Fs_final/1e3));
xlabel('Sample index'); ylabel('Amplitude');