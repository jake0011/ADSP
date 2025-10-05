% MATLAB_OFDM_64QAM_Simulation_CP.m
% Comprehensive OFDM chain simulation (64-QAM, CP handling, channel + sync + equalization)
% - Implements Schmidl & Cox timing sync, fractional+integer CFO estimation,
%   pilot-assisted LS vs MMSE channel estimation, per-tone equalizers (ZF/MMSE),
%   BER vs Eb/N0, and BER vs Doppler (ICI).
% - Uses an 802.11a-style preamble and pilot layout.
% - Authors: Jonathan Edem Akuaku, Felicia Archeresuah, Grace Wiredu - 2025-09-22
%
% Notes:
% - Sampling units are seconds; frequency units are Hz; subcarrier units are in
%   subcarrier spacing increments unless otherwise noted.
% - The script is organized into: parameters, preamble/pilots, experiments (timing,
%   CFO, channel-est MSE, BER vs Eb/N0, BER vs Doppler) and helper functions.
%
clear; clc; close all; rng(0);

%% ------------------------ System parameters ------------------------------
Fs   = 20e6;      % sampling (Hz)
N    = 64;        % FFT size (number of subcarriers)
Ncp  = 16;        % cyclic prefix length (samples)
df   = Fs / N;    % subcarrier spacing (Hz)
Tfft = N / Fs;    % FFT duration (s)
Tcp  = Ncp / Fs;  % CP duration (s)
Tsym = Tfft + Tcp; % OFDM symbol duration (s)
k    = 6;         % bits per 64-QAM symbol
M    = 64;        % modulation order (64-QAM)

% Schmidl & Cox / preamble parameters
L           = 32;   % half-length for long training (two identical halves -> 64)
short_len   = 16;   % short symbol length (samples)
num_short   = 10;   % repeat count for short preamble (10 * 16 = 160 samples)

% Pilot & data allocation consistent with 802.11a
active_pos = [-26:-1, 1:26];        % indices of active data carriers (DC excluded)
idx_active  = 33 + active_pos;      % MATLAB indexing (1-based)
pilot_rel   = [-21 -7 7 21];        % pilot tone relative positions
pilot_idx   = 33 + pilot_rel;       % MATLAB indexing for pilots
pilot_pattern = [1 1 1 -1].';       % BPSK pilot pattern (col vector)
pilot_symbol_freq = zeros(N,1);
pilot_symbol_freq(pilot_idx) = pilot_pattern;
data_idx = setdiff(idx_active, pilot_idx).';   % data carrier indices (column)

% Simulation control: SNR ranges and number of frames
EbN0_timing = 5:1:15;         % Eb/N0 points for timing variance experiment
EbN0_ber    = 0:2:20;         % Eb/N0 points for BER experiment
nFrames_timing = 2000;        % frames for timing experiment (statistical)
nFrames_ber    = 2000;        % frames for BER measurement (statistical)

% Channel PDP target (RMS delay and tap delays)
tau_rms_target = 200e-9;                      % target RMS delay (seconds)
tap_delays_s   = [0 1 2 4 8] * (1/Fs);        % tap delays (s): 0, 50ns, 100ns, 200ns, 400ns

% Compute exponential PDP parameters to match desired RMS delay
alpha = find_alpha_for_rms(tap_delays_s, tau_rms_target);
p = exp(-alpha*tap_delays_s); p = p / sum(p); % normalized tap power profile

fprintf('Assumed carrier frequency (for Doppler calcs) = 2.4 GHz (GUESS).\n');
fprintf('Tap delays (s): %s\n', mat2str(tap_delays_s));
fprintf('Tap powers (normalized): %s\n', mat2str(p.',6));

%% --------------------------- Helpers / short aliases ---------------------
ifft64 = @(X) ifft(X, N);
fft64  = @(x) fft(x, N);

% Note: 'firwin' style is not used in MATLAB; mapping/demapping use built-ins
map64qam   = @(bits) qammod(bits, M, 'InputType','bit', 'UnitAveragePower',true);
demap64qam = @(sym)  qamdemod(sym, M, 'OutputType','bit','UnitAveragePower',true);

add_cp    = @(x) [x(end-Ncp+1:end); x];   % append CP to time-domain OFDM symbol
remove_cp = @(x) x(Ncp+1 : Ncp+N);       % remove CP from received symbol

%% --------------------------- Preamble & Pilot construction --------------
% Short preamble: synthesize a 16-sample QPSK tile and repeat it 10 times
short_sym_time      = exp(1j*pi/2 * (randi([0 3], short_len, 1)));
short_preamble_time = repmat(short_sym_time, num_short, 1); % 160 samples total

% Long preamble: two identical halves for S&C timing/CFO fractional estimate
long_half        = exp(1j * 2*pi * rand(L,1));  % random complex half
long_symbol_time = [long_half; long_half];      % 64 samples total

% Complete preamble for a frame: short + long
preamble_time = [short_preamble_time; long_symbol_time];

% Block pilot OFDM symbol (frequency domain)
block_pilot_freq = zeros(N,1);
block_pilot_freq(data_idx)  = 1;               % +1 on all data tones
block_pilot_freq(pilot_idx) = pilot_pattern;   % pilots in known pattern
pilot_time_no_cp = ifft64(block_pilot_freq);
pilot_time       = add_cp(pilot_time_no_cp);   % add CP for TX pilot symbol

%% ----------------------- 1) Timing variance (Schmidl & Cox) ------------
fprintf('Running timing variance experiment (Schmidl & Cox with 802.11a-style preamble)...\n');
timing_var = zeros(size(EbN0_timing));

for idx = 1:length(EbN0_timing)
    EbN0dB = EbN0_timing(idx);
    sigma2 = compute_noise_variance(EbN0dB, k, N, Ncp);
    errors = zeros(nFrames_timing,1);

    for f = 1:nFrames_timing
        % Build a realistic transmit frame (preamble + pilot + one data symbol)
        Xd = zeros(N,1);
        bits_tmp = randi([0 1], numel(data_idx) * k, 1);
        Xd(data_idx) = map64qam(bits_tmp);
        data_time_no_cp = ifft64(Xd);
        data_time = add_cp(data_time_no_cp);

        tx_frame = [preamble_time; pilot_time; data_time];

        % AWGN channel
        rx = tx_frame + sqrt(sigma2/2) * (randn(size(tx_frame)) + 1j * randn(size(tx_frame)));

        % Evaluate S&C metric in a local search window around where the long symbol should be
        search_start = max(1, length(short_preamble_time) - 20);
        search_end   = search_start + 200;
        Mmet = zeros(search_end,1);
        for d = search_start:search_end - 2*L
            r1 = rx(d : d + L - 1);
            r2 = rx(d + L : d + 2*L - 1);
            P  = sum(conj(r1) .* r2);
            R  = sum(abs(r2).^2);
            Mmet(d) = abs(P)^2 / (R^2 + eps);
        end
        [~, d_est] = max(Mmet);
        true_long_start = length(short_preamble_time) + 1;
        errors(f) = d_est - true_long_start;
    end

    timing_var(idx) = var(errors);
    fprintf('Eb/N0 = %2.1f dB: timing variance = %.4f samples^2\n', EbN0dB, timing_var(idx));
end

figure; plot(EbN0_timing, timing_var, '-o'); grid on;
xlabel('Eb/N0 (dB)'); ylabel('Timing estimation variance (samples^2)');
title('Schmidl & Cox timing detection variance vs Eb/N0 (802.11a-style preamble)');

%% ------------- 2) CFO estimator: fractional (long) + integer (pilot) ----------
fprintf('Testing CFO estimator (fractional + integer) with 802.11a-style preamble...\n');
EbN0_test = [5 10 15];
cfo_true_list = [-0.4 -0.1 0.05 0.3]; % in subcarrier units
residuals = zeros(length(EbN0_test), length(cfo_true_list));

for ie = 1:length(EbN0_test)
    EbN0dB = EbN0_test(ie);
    sigma2 = compute_noise_variance(EbN0dB, k, N, Ncp);

    for ic = 1:length(cfo_true_list)
        eps_true_sc = cfo_true_list(ic);     % CFO in subcarrier units
        eps_true_cs = eps_true_sc / N;       % cycles per sample
        ntrials = 500; 
        resid = zeros(ntrials,1);

        for t = 1:ntrials
            tx_frame = [preamble_time; pilot_time];
            n = (0:length(tx_frame)-1).';
            rx = tx_frame .* exp(1j*2*pi*eps_true_cs*n);
            rx = rx + sqrt(sigma2/2) * (randn(size(rx)) + 1j*randn(size(rx)));

            % Assume timing is known: estimate fractional CFO using long training halves
            d = length(short_preamble_time) + 1;
            r1 = rx(d : d + L - 1);
            r2 = rx(d + L : d + 2*L - 1);
            P  = sum(conj(r1) .* r2);
            eps_frac_hat_cs = angle(P) / (2*pi*L); % cycles/sample

            % Correct fractional CFO on received samples
            rx_corr = rx .* exp(-1j * 2*pi * eps_frac_hat_cs * (0:length(rx)-1).');

            % Integer CFO: detect integer subcarrier shift using pilot correlation
            pilot_pos = d + 2*L; % pilot starts immediately after long training
            pilot_rx_cp = rx_corr(pilot_pos : pilot_pos + Ncp + N - 1);
            pilot_rx_td = remove_cp(pilot_rx_cp);
            Y = fft64(pilot_rx_td);
            ref = block_pilot_freq;

            corr_vals = zeros(N,1);
            for shift = 0:N-1
                ref_shift = circshift(ref, shift);
                corr_vals(shift+1) = abs(ref_shift' * conj(Y));
            end
            [~, shift_hat] = max(corr_vals);
            shift_hat = shift_hat - 1; % MATLAB 1-based -> 0-based shift

            eps_hat_total_cs = eps_frac_hat_cs + shift_hat/N;   % cycles/sample
            % Compute residual error in subcarrier units
            resid(t) = abs(eps_true_sc - eps_hat_total_cs * N);
        end

        residuals(ie, ic) = mean(resid);
        fprintf('Eb/N0=%d dB, true eps=%.3f subcarrier: residual (avg)=%.4f subcarrier units\n', ...
                 EbN0dB, eps_true_sc, residuals(ie,ic));
    end
end

fprintf('Fractional+Integer residuals (avg):\n');
disp(residuals);

%% ---- 3) Channel estimation MSE comparison: LS vs MMSE (block pilot) ----
fprintf('Running channel estimation MSE comparison (LS vs MMSE) with 802.11a pilot layout...\n');
EbN0_list = 0:2:20;
mse_ls   = zeros(size(EbN0_list));
mse_mmse = zeros(size(EbN0_list));

% Frequency-domain channel covariance computed from PDP
Rhh = channel_freq_covariance_from_pdp(p, tap_delays_s, N, Fs);

for ii = 1:length(EbN0_list)
    EbN0dB = EbN0_list(ii);
    sigma2 = compute_noise_variance(EbN0dB, k, N, Ncp);
    ntrials = 1000; se_ls = 0; se_mmse = 0;

    for t = 1:ntrials
        % Draw a random channel impulse response (time domain) from PDP
        h_time = (randn(length(p),1) + 1j*randn(length(p),1)) / sqrt(2) .* sqrt(p(:));
        H_freq = fft([h_time; zeros(N-length(h_time),1)], N);

        % Transmit pilot OFDM symbol (with CP), apply channel (linear conv), add AWGN
        tx = pilot_time;
        rx = conv(tx, [h_time; zeros(N-length(h_time),1)]);
        rx = rx(1:length(tx));
        rx = rx + sqrt(sigma2/2) * (randn(size(rx)) + 1j*randn(size(rx)));

        % Receiver: remove CP and FFT pilot symbol
        y_td = remove_cp(rx(1:Ncp+N));
        Y = fft64(y_td);

        Xp = block_pilot_freq;   % known frequency-domain pilot symbol
        H_true = H_freq;

        % LS estimate on pilot tones and linear interpolation to full band
        Hls_pilots = Y(pilot_idx) ./ Xp(pilot_idx);
        Hls_full = interp1(pilot_idx, Hls_pilots, 1:N, 'linear', 'extrap').';

        % MMSE (Wiener) interpolation leveraging channel covariance Rhh
        Rp = Rhh(pilot_idx, pilot_idx);
        Rhp = Rhh(:, pilot_idx);
        P = numel(pilot_idx);
        Xp_mat = diag(Xp(pilot_idx));
        A = Xp_mat * Rp * Xp_mat' + sigma2 * eye(P);
        W = Rhp / A * Xp_mat;       % N x P Wiener interpolation matrix
        Hmmse = (W * Hls_pilots).'; % mmse estimate across all tones

        se_ls   = se_ls   + mean(abs(H_true - Hls_full).^2);
        se_mmse = se_mmse + mean(abs(H_true - Hmmse.').^2);
    end

    mse_ls(ii)   = se_ls / ntrials;
    mse_mmse(ii) = se_mmse / ntrials;
    fprintf('EbN0=%d dB: MSE_LS=%.4e, MSE_MMSE=%.4e\n', EbN0dB, mse_ls(ii), mse_mmse(ii));
end

figure; semilogy(EbN0_list, mse_ls, '-o', EbN0_list, mse_mmse, '-x'); grid on;
legend('LS','MMSE'); xlabel('Eb/N0 (dB)'); ylabel('MSE');
title('Channel estimation MSE: LS vs MMSE (802.11a pilots, CP handled)');

%% ----- 4) BER simulations (MMSE vs ZF) — with CP handled at RX ------
fprintf('Running BER simulations (MMSE vs ZF). This section may be time-consuming...\n');
BER_mmse = zeros(size(EbN0_ber));
BER_zf   = zeros(size(EbN0_ber));

for ii = 1:length(EbN0_ber)
    EbN0dB = EbN0_ber(ii);
    sigma2 = compute_noise_variance(EbN0dB, k, N, Ncp);
    nErr_mmse = 0; nErr_zf = 0; nBitsTotal = 0;

    for f = 1:nFrames_ber
        % Data symbol generation on data subcarriers
        nData = numel(data_idx);
        bits = randi([0 1], nData * k, 1);
        symbols = map64qam(bits);

        X = zeros(N,1);
        X(data_idx)  = symbols;
        X(pilot_idx) = pilot_pattern;

        data_time_no_cp = ifft64(X);
        data_time = add_cp(data_time_no_cp);
        tx_frame = [preamble_time; pilot_time; data_time];

        % Rayleigh channel draw per frame from PDP and linear convolution
        h_time = (randn(length(p),1) + 1j*randn(length(p),1)) / sqrt(2) .* sqrt(p(:));
        rx = conv(tx_frame, [h_time; zeros(N-length(h_time),1)]);
        rx = rx(1:length(tx_frame));

        % Add small random CFO (subcarrier units) and AWGN
        eps_sc = (randn * 0.02);  % random small CFO (subcarrier units)
        n = (0:length(rx)-1).';
        rx = rx .* exp(1j * 2*pi * (eps_sc / N) * n);
        rx = rx + sqrt(sigma2/2) * (randn(size(rx)) + 1j * randn(size(rx)));

        % Assume frame sync found: locate pilot symbol and process
        long_start = length(short_preamble_time) + 1;
        pilot_pos  = long_start + 2*L;   % pilot symbol position in frame

        % Pilot: remove CP, FFT, LS estimate on pilots and interpolate
        pilot_rx_cp = rx(pilot_pos : pilot_pos + Ncp + N - 1);
        pilot_rx_td = remove_cp(pilot_rx_cp);
        Y_pilot = fft64(pilot_rx_td);

        Hls_pilots = Y_pilot(pilot_idx) ./ block_pilot_freq(pilot_idx);
        Hls_full = interp1(pilot_idx, Hls_pilots, 1:N, 'linear', 'extrap').';

        % Per-tone equalizer designs (ZF and simple MMSE)
        H_est = Hls_full;
        G_zf = 1 ./ (H_est + 1e-12);
        G_mmse = conj(H_est) ./ (abs(H_est).^2 + sigma2);

        % Data: remove CP, FFT, equalize and demap
        data_pos = pilot_pos + (Ncp + N);
        data_rx_cp = rx(data_pos : data_pos + Ncp + N - 1);
        data_rx_td = remove_cp(data_rx_cp);
        Y_data = fft64(data_rx_td);

        eq_zf = G_zf .* Y_data;
        eq_mmse = G_mmse .* Y_data;

        rx_symbols_zf = eq_zf(data_idx);
        rx_symbols_mmse = eq_mmse(data_idx);

        bits_hat_zf = demap64qam(rx_symbols_zf);
        bits_hat_mmse = demap64qam(rx_symbols_mmse);

        nErr_zf = nErr_zf + sum(bits_hat_zf ~= bits);
        nErr_mmse = nErr_mmse + sum(bits_hat_mmse ~= bits);
        nBitsTotal = nBitsTotal + length(bits);
    end

    BER_zf(ii) = nErr_zf / nBitsTotal;
    BER_mmse(ii) = nErr_mmse / nBitsTotal;
    fprintf('Eb/N0=%d dB: BER_ZF=%.3e, BER_MMSE=%.3e\n', EbN0dB, BER_zf(ii), BER_mmse(ii));
end

% AWGN reference (approximate uncoded 64-QAM bound)
EbN0_lin = 10.^(EbN0_ber / 10);
EsN0_lin = EbN0_lin * k;
ser_awgn = 4 * (1 - 1/sqrt(M)) .* qfunclocal(sqrt(3/(M-1) * EsN0_lin));
ber_awgn_approx = ser_awgn / k;

figure; semilogy(EbN0_ber, BER_mmse, '-o', EbN0_ber, BER_zf, '-x', EbN0_ber, ber_awgn_approx, '--'); grid on;
legend('MMSE (sim)','ZF (sim)','AWGN (approx)');
xlabel('Eb/N0 (dB)'); ylabel('BER');
title('BER: uncoded 64-QAM OFDM (802.11a pilots/preamble, CP handled)');

%% -------- 5) BER vs Doppler (ICI impact) at Eb/N0 = 15 dB -------------
fprintf('Testing ICI impact across Doppler values...\n');
fc = 2.4e9;  % assumed carrier (Hz)
fD_list = [0 100 500 1000 2000]; % Doppler values to test (Hz)
BER_dop  = zeros(size(fD_list));
EbN0dB = 15; sigma2 = compute_noise_variance(EbN0dB, k, N, Ncp);

for ii = 1:length(fD_list)
    fD = fD_list(ii);
    nErr = 0; nBits = 0; nRuns = 500;
    rho = besselj(0, 2*pi*fD*Tsym); % symbol-to-symbol correlation (Clarke model approx)

    for r = 1:nRuns
        nData = numel(data_idx);
        bits = randi([0 1], nData * k, 1);
        symbols = map64qam(bits);

        X = zeros(N,1); X(data_idx) = symbols; X(pilot_idx) = pilot_pattern;
        data_time_no_cp = ifft64(X);
        data_time = add_cp(data_time_no_cp);
        tx_frame = [preamble_time; pilot_time; data_time];

        % Time-varying channel modeled as first-order Gauss-Markov (rho)
        h_prev = (randn(length(p),1) + 1j*randn(length(p),1)) / sqrt(2) .* sqrt(p(:));
        h_time = rho * h_prev + sqrt(max(0, 1 - rho^2)) * (randn(length(p),1) + 1j*randn(length(p),1)) / sqrt(2) .* sqrt(p(:));

        rx = conv(tx_frame, [h_time; zeros(N-length(h_time),1)]);
        rx = rx(1:length(tx_frame));
        rx = rx + sqrt(sigma2/2) * (randn(size(rx)) + 1j * randn(size(rx)));

        long_start = length(short_preamble_time) + 1;
        pilot_pos = long_start + 2*L;

        % Pilot-based LS channel estimation (interpolate to all tones)
        pilot_rx_cp = rx(pilot_pos : pilot_pos + Ncp + N - 1);
        pilot_rx_td = remove_cp(pilot_rx_cp);
        Y_pilot = fft64(pilot_rx_td);
        Hls_pilots = Y_pilot(pilot_idx) ./ block_pilot_freq(pilot_idx);
        H_est = interp1(pilot_idx, Hls_pilots, 1:N, 'linear', 'extrap').';

        % Data reception and MMSE equalization
        data_pos = pilot_pos + (Ncp + N);
        data_rx_cp = rx(data_pos : data_pos + Ncp + N - 1);
        data_rx_td = remove_cp(data_rx_cp);
        Y_data = fft64(data_rx_td);

        G_mmse = conj(H_est) ./ (abs(H_est).^2 + sigma2);
        eq_mmse = G_mmse .* Y_data;
        rx_symbols_mmse = eq_mmse(data_idx);
        bits_hat = demap64qam(rx_symbols_mmse);

        nErr = nErr + sum(bits_hat ~= bits);
        nBits = nBits + length(bits);
    end

    BER_dop(ii) = nErr / nBits;
    fprintf('Doppler=%d Hz: BER=%.3e\n', fD, BER_dop(ii));
end

figure; semilogy(fD_list, BER_dop, '-o'); grid on;
xlabel('Doppler (Hz)'); ylabel('BER at Eb/N0=15 dB');
title('BER vs Doppler (ICI floor behavior)');

%% -------------------- Local helper functions (bottom) -------------------
function alpha = find_alpha_for_rms(tau_vec, tau_rms_target)
    % Find exponential decay alpha such that exponential PDP has desired RMS delay
    obj = @(a) (compute_rms_for_alpha(a, tau_vec) - tau_rms_target).^2;
    a_low = 1e2; a_high = 1e10;
    try
        alpha = fminbnd(obj, a_low, a_high, optimset('TolX',1e-6,'Display','off'));
        if ~isfinite(alpha) || alpha <= 0
            error('bad alpha');
        end
    catch
        % Fallback grid search if optimizer fails
        grid = logspace(log10(a_low), log10(a_high), 200);
        vals = arrayfun(@(a) obj(a), grid);
        [~, idxmin] = min(vals);
        alpha = grid(idxmin);
    end
end

function rms = compute_rms_for_alpha(alpha, tau_vec)
    p = exp(-alpha * tau_vec);
    p = p / sum(p);
    tau_mean = sum(p .* tau_vec);
    rms = sqrt(sum(p .* (tau_vec - tau_mean).^2));
end

function sigma2 = compute_noise_variance(EbN0dB, k, N, Ncp)
    % Convert Eb/N0 (dB) to complex-sample noise variance assuming unit-power symbols
    EbN0 = 10^(EbN0dB / 10);
    EsN0 = EbN0 * k;
    useful_frac = N / (N + Ncp);      % account for CP overhead
    EsN0_eff = EsN0 * useful_frac;
    sigma2 = 1 / max(EsN0_eff, 1e-12); % noise variance per complex sample
end

function Rhh = channel_freq_covariance_from_pdp(p, delays, N, Fs)
    % Build frequency-domain covariance matrix R_HH from PDP p and tap delays
    Rhh = zeros(N, N);
    for i = 1:N
        for j = 1:N
            df_ij = (i - j) * Fs / N; % frequency difference (Hz)
            s = 0;
            for t = 1:length(p)
                s = s + p(t) * exp(-1j * 2 * pi * df_ij * delays(t));
            end
            Rhh(i, j) = s;
        end
    end
end

function y = qfunclocal(x)
    % Q-function (tail probability of standard normal)
    y = 0.5 * erfc(x ./ sqrt(2));
end