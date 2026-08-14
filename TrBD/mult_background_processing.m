function [Yb, A_true_clip, s] = mult_background_processing(data, s)

impresp = data.impulseRespRec; % [Nsamp x Nwm]
delaycomp = data.delayCompRec(:).'; % [1 x Nwm]
yb = hilbert(data.outputSignal(:));
lfm = hilbert(s.lfm(:));

Fs0 = s.cfg.Fs;
useBB = isfield(s.cfg, 'baseband') && s.cfg.baseband;

% ------------------------------------------------------------
% Optional complex baseband conversion
% ------------------------------------------------------------
if useBB
    if ~isfield(s.cfg, 'fc') || isempty(s.cfg.fc)
        error('s.cfg.fc must be provided when s.cfg.baseband == true.');
    end

    fc = s.cfg.fc;
    ny = (0:numel(yb)-1).';
    nlf = (0:numel(lfm)-1).';

    yb = yb  .* exp(-1j*2*pi*fc*ny /Fs0);
    lfm = lfm .* exp(-1j*2*pi*fc*nlf/Fs0);

    s.cfg.signal_domain = 'analytic-baseband';
else
    s.cfg.signal_domain = 'analytic-passband';
end

% ------------------------------------------------------------
% True CIR from impulse responses, aligned across waymarks
% ------------------------------------------------------------
[Nsamp, Nwm] = size(impresp);
f = (-floor(Nsamp/2):ceil(Nsamp/2)-1).' * (Fs0/Nsamp);

Hc = fftshift(fft(impresp, [], 1), 1);
Hc = Hc .* exp(-1j*2*pi*f*(delaycomp - delaycomp(1)));
A_true = ifft(ifftshift(Hc, 1), [], 1); % [Nsamp x Nwm]

if useBB
    tau_imp = (0:Nsamp-1).' / Fs0;
    A_true = A_true .* exp(-1j*2*pi*s.cfg.fc*tau_imp);
end

% ------------------------------------------------------------
% Resample to processing grid
% ------------------------------------------------------------
lfm_bb = lfm; % minimal fix: valid also when no resampling/basebanding

if useBB
    Fs_proc = 5000; % keep your current choice
    [p, q] = rat(Fs_proc/Fs0, 1e-10);

    if p ~= 1 || q ~= 1
        yb = resample(yb, p, q);
        lfm_bb = resample(lfm, p, q);
        A_true = resample(A_true, p, q);
    end

    s.cfg.Fs_used = Fs0 * (p/q);
else
    s.cfg.Fs_used = Fs0;
end

% ------------------------------------------------------------
% Waveform / dimensions on processing grid
% ------------------------------------------------------------
s.lfm_bb = lfm_bb(:);
s.t_lfm_bb = (0:numel(s.lfm_bb)-1).' / s.cfg.Fs_used;

s.cfg.Nlfm = numel(s.lfm_bb);
s.cfg.L = round(s.cfg.Trec * s.cfg.Fs_used);
s.cfg.N = s.cfg.Nlfm + s.cfg.L - 1;

% ------------------------------------------------------------
% Slice continuous background into ping records
% ------------------------------------------------------------
Zpad_extra = 0.001;
startIdx0 = floor((s.cfg.initial_delay - Zpad_extra) * s.cfg.Fs_used) + 1;
step_samp = round(s.cfg.PRI * s.cfg.Fs_used);

Yb = zeros(s.cfg.N, s.cfg.Np, 'like', yb);
k_filled = 0;

for k = 1:s.cfg.Np
    idxStart = startIdx0 + (k-1)*step_samp;
    idxEnd = idxStart + s.cfg.N - 1;

    if idxStart < 1 || idxEnd > numel(yb)
        break;
    end

    k_filled = k_filled + 1;
    Yb(:, k_filled) = yb(idxStart:idxEnd);
end

Yb = Yb(:, 1:k_filled);

s.cfg.Np = k_filled;
s.cfg.t_ping = (0:k_filled-1).' * s.cfg.PRI;

% ------------------------------------------------------------
% Delay axis of sliced record
% ------------------------------------------------------------
s.cfg.tau_axis = (0:s.cfg.L-1).' / s.cfg.Fs_used + ...
                         (s.cfg.initial_delay - Zpad_extra);
s.cfg.tau_excess_delay = s.cfg.tau_axis - s.cfg.tau_axis(1);

% ------------------------------------------------------------
% Interpolate true CIR at ping times
% ------------------------------------------------------------
s.cfg.t_wm = (0:Nwm-1) * s.cfg.waymarkPeriod;
A_true_interp = interp1(s.cfg.t_wm, A_true.', s.cfg.t_ping, 'linear', 'extrap').';

% ------------------------------------------------------------
% Clip true CIR to record window
% ------------------------------------------------------------
Zpad_sim = 0.1 * s.cfg.Timp_min;
row0_true = floor(Zpad_sim * s.cfg.Fs_used) + 1;
pad_samp = max(0, floor(Zpad_extra * s.cfg.Fs_used));

A_true_pad = [zeros(pad_samp, s.cfg.Np, 'like', A_true_interp); A_true_interp];

r0 = row0_true;
r1 = r0 + s.cfg.L - 1;

if r1 > size(A_true_pad, 1)
    A_true_pad(end+1:r1, :) = 0;
end

A_true_clip = A_true_pad(r0:r1, :);

end
