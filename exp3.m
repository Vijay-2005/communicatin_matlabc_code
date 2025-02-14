%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clean, Error-Free Complete Code for PAM System Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;

%-----------------------------
% 1) System Parameters
%-----------------------------
N        = 1e5;       % Number of symbols to transmit
M        = 4;         % Modulation order (4-PAM)
L        = 4;         % Oversampling factor (samples per symbol)
beta     = 0.3;       % Roll-off factor for SRRC filter
Nsym     = 8;         % SRRC filter span (in symbol durations)
EbN0dB   = 1000;      % Eb/N0 in dB (very large => nearly no noise)
snr      = 10*log10(log2(M)) + EbN0dB;  % Convert Eb/N0 to SNR in dB

%-----------------------------
% 2) Generate Random Symbols
%-----------------------------
d = ceil(M .* rand(1, N));    % Random symbols in [1, M]
u = pammod(d - 1, M);         % 4-PAM modulation (shift to [0..M-1])

figure(1);
stem(real(u));
title('1) PAM Modulated Symbols u(k)');
xlim([0 20]); ylim([-5 5]);

%-----------------------------
% 3) Oversample (Zero-Stuffing)
%-----------------------------
v = [u; zeros(L-1, length(u))];  % Insert L-1 zeros after each symbol
v = v(:).';                      % Convert to single stream

figure(2);
stem(real(v));
title('2) Oversampled Symbols v(n)');
xlim([0 150]); ylim([-5 5]);

%-----------------------------
% 4) Transmit Filter (SRRC)
%-----------------------------
[p, t, filtDelay] = srrcFunction(beta, L, Nsym);  % Design SRRC filter
s = conv(v, p, 'full');                           % TX pulse shaping

figure(3);
plot(real(s), 'r');
title('3) Pulse-Shaped Symbols s(n)');
xlim([0 150]); ylim([-5 5]);

%-----------------------------
% 5) AWGN Channel
%-----------------------------
r = add_awgn_noise(s, snr, L);   % Add noise at specified SNR

figure(4);
plot(real(r), 'r');
title('4) Received Signal r(n)');
xlim([0 150]); ylim([-5 5]);

%-----------------------------
% 6) Matched Filtering (SRRC)
%-----------------------------
vCap = conv(r, p, 'full');   % RX matched filter

figure(5);
plot(real(vCap), 'r');
title('5) After Matched Filtering \hat{v}(n)');
xlim([0 150]); ylim([-20 20]);

%-----------------------------
% 7) Symbol-Rate Sampling
%-----------------------------
%   - We skip the first 2*filtDelay samples,
%   - then pick every L-th sample,
%   - stop 2*filtDelay from the end,
%   - and divide by L to compensate the filter gain.
uCap = vCap(2*filtDelay + 1 : L : end - (2*filtDelay)) / L;

figure(6);
stem(real(uCap));
title('6) After Symbol-Rate Sampler \hat{u}(n)');
xlim([0 20]); ylim([-5 5]);

% Demodulate
dCap = pamdemod(uCap, M);

%-----------------------------
% 8) Eye Diagram
%-----------------------------
%   - We will overlay 100 traces,
%   - Each covering 3 symbol durations,
%   - Starting at sample index 2*filtDelay.
figure(7); 
hold on;
plotEyeDiagram(vCap, L, 100, 2*filtDelay, 3);  % (signal, L, nTraces, offset, nSymbolsPerTrace)
xlim([0 3]); ylim([-15 15]);
title('7) Eye Diagram (Matched Filter Output)');
xlabel('Time (in symbol durations)');
ylabel('Amplitude');
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Supporting Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-----------------------------------------------------------
% Square-Root Raised Cosine (SRRC) Filter
%-----------------------------------------------------------
function [p, t, filtDelay] = srrcFunction(beta, L, Nsym)
    % beta  : roll-off factor
    % L     : Oversampling factor
    % Nsym  : Filter span in symbol durations
    % p     : SRRC filter coefficients
    % t     : Time vector
    % filtDelay : Filter delay in samples

    T = 1;                 % Symbol duration (normalized)
    Ts = T / L;            % Sampling period
    t = -Nsym/2 : Ts : Nsym/2;  % Time vector

    % SRRC formula
    p = (sinc(t/T) .* cos(pi*beta*t/T)) ./ (1 - (2*beta*t/T).^2);

    % Fix potential divide-by-zero/singularity at t = ±T/(2β)
    singularity = find(abs(2*beta*t/T - 1) < 1e-9);
    if ~isempty(singularity)
        p(singularity) = pi/4 * sinc(1/(2*beta));
    end

    % Normalize filter energy
    p = p / sqrt(sum(p.^2));

    % Filter delay in samples
    filtDelay = (Nsym * L) / 2;
end

%-----------------------------------------------------------
% Add AWGN Noise
%-----------------------------------------------------------
function r = add_awgn_noise(s, snr_dB, L)
    % s       : Transmitted signal
    % snr_dB  : SNR in dB
    % L       : Oversampling factor (not strictly needed for basic AWGN)
    
    snr_linear   = 10^(snr_dB / 10);
    signal_power = mean(abs(s).^2);
    noise_power  = signal_power / snr_linear;

    % Complex white Gaussian noise
    noise = sqrt(noise_power/2) .* (randn(size(s)) + 1i*randn(size(s)));

    % Noisy signal
    r = s + noise;
end

%-----------------------------------------------------------
% Eye Diagram Plotter
%-----------------------------------------------------------
function plotEyeDiagram(signal, L, nTraces, offset, nSymbolsPerTrace)
    % Plots an eye diagram by overlaying multiple segments
    % of 'signal' after matched filtering.
    %
    %   signal           : The oversampled signal (vector).
    %   L                : Oversampling factor (samples per symbol).
    %   nTraces          : Number of traces (symbol intervals) to overlay.
    %   offset           : Starting sample index (e.g., 2*filtDelay).
    %   nSymbolsPerTrace : How many symbol durations each trace covers.
    
    samplesPerTrace = nSymbolsPerTrace * L;
    cmap = lines(nTraces);   % Distinct colors for each trace

    for i = 1 : nTraces
        startIdx = offset + (i - 1)*L;
        stopIdx  = startIdx + samplesPerTrace - 1;
        
        % Only plot if we have enough samples
        if stopIdx <= length(signal)
            traceData = signal(startIdx : stopIdx);
            plot((0:samplesPerTrace-1)/L, real(traceData), ...
                 'Color', cmap(mod(i-1, size(cmap,1)) + 1, :));
        else
            break;  % No more data to plot
        end
    end
end
