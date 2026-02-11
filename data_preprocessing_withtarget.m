function [Ymat,A_true_clip,S,U, S_tgt_mat, cfg]=data_preprocessing_withtarget(data,s,cfg)

impresp= data.impulseRespRec;
delaycomp=data.delayCompRec;

y=data.outputSignal;


s = hilbert(s(:));
y = hilbert(y(:));

if isfield(cfg,'baseband') && cfg.baseband
    if ~isfield(cfg,'fc')
        error('cfg.fc (Hz) must be provided when cfg.baseband == true.');
    end
    ns = (0:numel(s)-1).';         % sample indices
    ny = (0:numel(y)-1).';
    s  = s .* exp(-1j*2*pi*cfg.fc*(ns/cfg.Fs));
    y  = y .* exp(-1j*2*pi*cfg.fc*(ny/cfg.Fs));
    cfg.signal_domain = 'analytic-baseband';
else
    cfg.signal_domain = 'analytic-passband';
end

Zpad_sim   = 0.1*cfg.Timp_min;



[Nsamp, Nwm] = size(impresp);

df = 1 / (Nsamp*1/cfg.Fs);
f  = ((-floor(Nsamp/2)):(ceil(Nsamp/2)-1)).' * df;

% Apply per-waymark *delay* Δ_m via phase ramp exp(-j 2π f Δ_m)
dCompRef = delaycomp(1);
A_true = zeros(size(impresp), 'like', impresp);
for m = 1:Nwm
    dRel = delaycomp(m) - dCompRef;                             % relative reverse compensation (s)
    Hc   = fftshift(fft(impresp(:,m)));                         % center spectrum
    Hc   = Hc .* exp(-1j*2*pi*f*dRel);                         % h(t - dRel)
    A_true(:,m) = ifft(ifftshift(Hc));                        % back to time (complex)
end


if isfield(cfg,'baseband') && cfg.baseband
    if ~isfield(cfg,'fc'), error('cfg.fc must be set when cfg.baseband==true.'); end
    tau_imp = (0:Nsamp-1).' / cfg.Fs;                        % delay samples (s)
    phase_tau = exp(-1j*2*pi*cfg.fc * tau_imp);         % e^{-j 2π f_c τ}
    A_true = phase_tau .* A_true;                        % broadcast over columns
end



%% downsampling
if cfg.baseband
    [p,q]=rat((cfg.Fs/2) /cfg.Fs , 1e-10);
    s = resample(s, p, q);                 % complex-safe, anti-alias included
    y = resample(y, p, q);
    A_true = resample(A_true, p, q);   % resamples each column over rows (delay)
    cfg.Fs     = cfg.Fs * (p/q);           % update to new Fs
    cfg.Nlfm  = floor(cfg.Fs*cfg.Tp);      % # fast-time samples per ping
    cfg.L  = floor(cfg.max_delay*cfg.Fs);% # delay taps
    cfg.N  = cfg.Nlfm+cfg.L;      % # time samples per ping
    cfg.tap_delays = cfg.initial_delay + (0:cfg.L-1)'/cfg.Fs;
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
fbb0   = cfg.f0 - cfg.fc;                % baseband start (Hz); e.g. -2000 if fc=3000, f0=1000
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

            % Complex baseband LFM on warped time
            xseg = w .*exp(1j*2*pi*( fbb0*tp + 0.5*kchirp*tp.^2 ));
            E = sum(abs(xseg).^2);     
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

    Z = (randn(size(Ymat)) + 1j*randn(size(Ymat))) / sqrt(2);
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



