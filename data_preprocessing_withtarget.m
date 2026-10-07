function [Ymat,A_true_clip,S,U, S_tgt_mat, cfg]=data_preprocessing_withtarget(data,s,cfg)

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




%%
cfg.tau_excess = cfg.tap_delays - cfg.tap_delays(1);
Zpad_extra = 0.005;
Ymat     = zeros(cfg.N, cfg.Np, 'like', y);
S_tgt_mat = zeros(cfg.N, cfg.Np, 'like', y);



kchirp = (cfg.f1 - cfg.f0) / cfg.Tp;     % Hz/s
t_fast = (0:cfg.N-1).' / cfg.Fs;         % [N x 1] fast-time for a ping window



startIdx = floor((cfg.initial_delay - Zpad_extra) * cfg.Fs) + 1;
k_filled = 0;
for k = 1:cfg.Np
    idxStart = startIdx + round(((k-1) * cfg.PRI) * cfg.Fs);
    idxEnd   = idxStart + cfg.N - 1;

    if idxStart < 1 || idxEnd > length(y)
        break;
    end

    k_filled = k_filled + 1;
    Ymat(:,k_filled) = y(idxStart:idxEnd);

    if cfg.add_target && isfinite(cfg.beta_target(k)) && isfinite(cfg.tau_target(k))

        beta_k = cfg.beta_target(k);
        tau_k  = cfg.tau_target(k);

        tau_rel = tau_k - (cfg.initial_delay - Zpad_extra);

        tprime = beta_k * (t_fast - tau_rel);

        gate = (tprime >= 0) & (tprime <= cfg.Tp);

        x = zeros(cfg.N,1,'like',y);
        if any(gate)
            tp = tprime(gate);                 

            w = hann(numel(tp),'periodic');

            xseg = w .* cos(2*pi*(cfg.f0*tp + 0.5*kchirp*tp.^2));
            E = sum(xseg.^2);
            xseg = xseg ./ sqrt(E);

            x(gate) = xseg;
        end

        S_tgt_mat(:,k_filled) = x;
    end
end

cfg.Np = k_filled;
Ymat = Ymat(:,1:cfg.Np);
if cfg.add_target
    S_tgt_mat = S_tgt_mat(:,1:cfg.Np);
end

%% Adding Noise and target
Y_clean = Ymat;
Pclutter = mean(abs(Y_clean).^2, 1);
sigmae2  = Pclutter ./ 10^(cfg.CNRdB/10); % 1 x Np
cfg.sigma_e = sqrt(sigmae2);     % 1 x Np
if cfg.add_target
    Ptgt_des = sigmae2 .* 10^(cfg.SNRdB/10);        % 1xNp
    Ptgt0 = mean(abs(S_tgt_mat).^2, 1);                    %1 x Np
    alpha = sqrt(Ptgt_des ./Ptgt0);                   % 1 x Np
    alpha(~isfinite(alpha)) = 0;

    S_tgt_mat = S_tgt_mat .* reshape(alpha, 1, []);

    % Add target to the clean measurement windows
    Ymat = Ymat + S_tgt_mat;
end
% ---- Add noise (after target) to achieve CNR wrt original clutter power ----
if cfg.add_noise

    Z = randn(size(Ymat));
    sigma_e_col = sqrt(sigmae2);        % 1 x Np
    E = Z .* reshape(sigma_e_col, 1, []);  % 1 x Np
    Ymat = Ymat + E;
    measCNRdB = 10*log10( mean(abs(Y_clean(:)).^2) / mean(abs(E(:)).^2) );
    fprintf('Measured CNR ≈ %.2f dB (supplied %.2f dB)\n', measCNRdB, cfg.CNRdB);


end
if cfg.add_target
    Ptgt_meas = mean(abs(S_tgt_mat(:,cfg.Ntrain+1:end)).^2, 1);
    SNR_ping_dB = 10*log10(Ptgt_meas/ mean(sigmae2(cfg.Ntrain+1:end)));
fprintf('Target SNR (power-based): mean %.2f dB \n', mean(SNR_ping_dB));
end


%% True CIR preprocesing


cfg.t_wm  = (0:Nwm-1) * cfg.waymarkPeriod;            % Waymark times (s)
% cfg.t_ping = cfg.initial_delay + (0:cfg.Np-1) * cfg.PRI;         % Your ping times (s)

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
