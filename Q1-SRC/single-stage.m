#% Single-stage polyphase-style SRC example (96 kHz -> 44.1 kHz)
% - Generates a two-tone test signal
% - Performs zero-stuff upsampling by L
% - Designs an anti-image / anti-alias FIR via 'designfilt'
% - Filters the upsampled signal, then downsamples by M


%% 1) Basic parameters and test signal
Fs = 96e3;                    % original sampling rate (Hz)
t  = 0 : 1/Fs : 0.002;        % short time vector for plotting (0..2 ms)

% Two test tones (used to exercise the passband and aliasing behaviour)
f1 = 20e3;  % tone 1 (Hz)
f2 = 30e3;  % tone 2 (Hz)
x  = sin(2*pi*f1.*t) + sin(2*pi*f2.*t);  % composite test signal

%% 2) SRC integer factors (interpolation / decimation)
L = 147;    % interpolation (upsampling) factor
M = 320;    % decimation (downsampling) factor

%% 3) Upsample (zero stuffing)
% Note: upsample() inserts L-1 zeros between input samples (no filtering)
x_up = upsample(x, L);

%% 4) Anti-imaging / anti-alias filter design (FIR)
% After upsampling the new sample rate is Fs_up. Design a lowpass FIR
Fs_up = Fs * L;                 % sampling rate after zero-stuff upsampling
Fpass = 20e3;                   % desired passband edge (Hz)
Fstop = 22e3;                   % stopband start (Hz)
F_nyq = Fs_up / 2;              % Nyquist frequency after upsampling

% Design normalized frequencies for designfilt (PassbandFrequency and StopbandFrequency
% expect normalized frequencies in the range (0,1) referring to Nyquist).
lpFilt = designfilt('lowpassfir', ...
    'PassbandFrequency', Fpass / F_nyq, ...
    'StopbandFrequency', Fstop / F_nyq, ...
    'PassbandRipple', 0.01, ...
    'StopbandAttenuation', 100, ...
    'DesignMethod','kaiserwin');

% Visualise filter (optional). Comment out if running batch/non-interactive.
fvtool(lpFilt)

% Extract FIR coefficients (numerical) and report filter length
b = lpFilt.Coefficients;
fprintf('Filter length: %d taps\n', length(b));

%% 5) Filter the upsampled signal (remove images introduced by upsampling)
% Use the same system object returned by designfilt; calling filter applies
% the FIR to the (zero-stuffed) upsampled signal.
x_filt = filter(lpFilt, x_up);

%% 6) Downsample to the new rate
% Keep every M-th sample (decimation). This produces the final output y.
y = downsample(x_filt, M);

%% 7) Report resulting sampling rate
Fs_new = Fs * (L / M);    % resulting sampling rate after SRC
fprintf('New sampling rate after conversion: %.2f Hz (%.3f kHz)\n', Fs_new, Fs_new/1e3);

%% 8) Simple plotting for quick visual checks
% The plotting indices and stem styles are preserved from the original script
subplot(3,1,1);
stem(0:100, x(1:101), 'filled', 'MarkerSize', 3);
title('Original signal (96 kHz)');
xlabel('Sample index'); ylabel('Amplitude');

subplot(3,1,2);
stem(0:443, x_up(1:444), 'filled', 'MarkerSize', 3);
title(sprintf('Upsampled (L = %d) — zero-stuffed', L));
xlabel('Sample index'); ylabel('Amplitude');

subplot(3,1,3);
stem(0:80, y(1:81), 'filled', 'MarkerSize', 3);
title(sprintf('After downsampling (M = %d) — Fs = %.1f kHz', M, Fs_new/1e3));
xlabel('Sample index'); ylabel('Amplitude');