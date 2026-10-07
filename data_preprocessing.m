function [Ymat,A_true_clip,S,U,cfg]=data_preprocessing(data,s,cfg)

impresp= data.impulseRespRec;
delaycomp=data.delayCompRec;

y=data.outputSignal;
s = s(:);
y = y(:);
if ~isreal(s) || ~isreal(y) || ~isreal(impresp)
    error('Real-passband processing requires real waveform, outputSignal, and impulseRespRec data.');
end
cfg.signal_domain = 'real-passband';

Zpad_sim   = 0.1*cfg.Timp_min;   



[~, Nwm] = size(impresp);

dCompRef = delaycomp(1);
A_true = zeros(size(impresp), 'like', impresp);
for m = 1:Nwm
    dRel = delaycomp(m) - dCompRef;
    A_true(:,m) = fractional_delay_real(impresp(:,m), dRel * cfg.Fs);
end

%% ---------------------- LFM pulse & linearization ---------------------------
tp=(0:cfg.Nlfm-1)'/cfg.Fs;
s=s(1:length(tp));
t0 = (cfg.Nlfm-1)/(2*cfg.Fs);
ds = gradient(s,tp);
u  = (tp-t0) .* ds;

S = convmtx(s(:), cfg.L);
U = convmtx(u(:), cfg.L);
[cfg.N,cfg.L]=size(S);




%% Reshaping the full measurment vector into N x K (samples per ping x number of pings)
cfg.tau_excess = cfg.tap_delays - cfg.tap_delays(1);
Zpad_extra = 0.005;    % EXTRA zero-pad (e.g., 5 ms)
% 3) Preallocate
Ymat = zeros(cfg.N, cfg.Np, 'like', y); % preserve type 
% 4) Starting index for first ping (account for propagation delay)
startIdx = floor((cfg.initial_delay-Zpad_extra)*cfg.Fs)+1;   % MATLAB 


k_filled = 0;
for k = 1:cfg.Np
    idxStart = startIdx + round(((k-1)*cfg.PRI)*cfg.Fs);
    idxEnd   = idxStart + cfg.N - 1;
    if idxStart < 1 || idxEnd > length(y)
        break;
    end

    Ymat(:,k) = y(idxStart:idxEnd);
    k_filled = k_filled + 1;
end
cfg.Np = k_filled;
Ymat = Ymat(:,1:cfg.Np);    % trim



if cfg.add_noise
    if ~isfield(cfg,'SNRdB') || isempty(cfg.SNRdB)
        error('cfg.SNRdB must be set before adding noise.');
    end

    Y_clean = Ymat;

    % per-ping average power
    Ps_col = mean(abs(Y_clean).^2, 1);                 % 1 x Np

    % noise std per ping so that Ps/sigma^2 = 10^(SNR/10)
    sigma_e_col = sqrt(Ps_col ./ (10^(cfg.SNRdB/10))); % 1 x Np

    Z = randn(size(Ymat));
    E = Z .* repmat(sigma_e_col, size(Ymat,1), 1);

    Ymat = Ymat + E;

    cfg.sigma_e = sigma_e_col;

    measSNRdB = 10*log10(mean(abs(Y_clean(:)).^2) / mean(abs(E(:)).^2));
    fprintf('Measured SNR ≈ %.2f dB (requested %.2f dB)\n', measSNRdB, cfg.SNRdB);
end




%% True CIR preprocesing


cfg.t_wm  = (0:Nwm-1) * cfg.waymarkPeriod;            % Waymark times (s)
cfg.t_ping = cfg.initial_delay + (0:cfg.Np-1) * cfg.PRI;         % Your ping times (s)

A_true_interp = interp1(cfg.t_wm, A_true.', cfg.t_ping, 'linear', 'extrap').';

row0_true = floor(Zpad_sim*cfg.Fs) + 1;           % account for simulator's built-in 0.1% zero-pad

pad_samp  = max(0, floor(Zpad_extra * cfg.Fs));
A_true_pad = [zeros(pad_samp, cfg.Np, 'like', A_true_interp); A_true_interp];

r0 = row0_true;
r1 = r0 + cfg.L - 1;
if r1 > size(A_true_pad,1)
    A_true_pad = [A_true_pad; zeros(r1 - size(A_true_pad,1), cfg.Np, 'like', A_true_pad)];
end
A_true_clip = A_true_pad(r0:r1, :);           % L x cfg.Np, aligned with Ymat but with extra leading zeros
cfg.tau_axis = (0:cfg.L-1).' / cfg.Fs + (cfg.initial_delay - Zpad_extra);

end

function y = fractional_delay_real(x, delay_samples)

N = numel(x);
query = (0:N-1).' - delay_samples;
idx0 = floor(query);
offsets = -7:8;
idx = idx0 + offsets;
d = query - idx;

window_arg = d / 8.5;
window = 0.5 + 0.5*cos(pi*window_arg);
window(abs(window_arg) >= 1) = 0;

z = pi*d;
sinc_weights = ones(size(z));
nonzero = z ~= 0;
sinc_weights(nonzero) = sin(z(nonzero)) ./ z(nonzero);
weights = window .* sinc_weights;

valid = idx >= 0 & idx < N;
idx = min(max(idx, 0), N-1) + 1;
samples = x(idx);
samples(~valid) = 0;
y = sum(weights .* samples, 2);

end

