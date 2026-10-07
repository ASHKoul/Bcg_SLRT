function [Yb, A_true_clip, s] = mult_background_processing(data, s)

impresp = data.impulseRespRec; % [Nsamp x Nwm]
delaycomp = data.delayCompRec(:).'; % [1 x Nwm]
yb = data.outputSignal(:);
lfm = s.lfm(:);

Fs0 = s.cfg.Fs;
if ~isreal(yb) || ~isreal(lfm) || ~isreal(impresp)
    error('Real-passband processing requires real outputSignal, lfm, and impulseRespRec data.');
end
s.cfg.signal_domain = 'real-passband';

% ------------------------------------------------------------
% True CIR from impulse responses, aligned across waymarks
% ------------------------------------------------------------
[~, Nwm] = size(impresp);
A_true = zeros(size(impresp), 'like', impresp);
for m = 1:Nwm
    dRel = delaycomp(m) - delaycomp(1);
    A_true(:, m) = fractional_delay_real(impresp(:, m), dRel * Fs0);
end

s.cfg.Fs_used = Fs0;

% ------------------------------------------------------------
% Waveform / dimensions on processing grid
% ------------------------------------------------------------
s.lfm_passband = lfm(:);
s.t_lfm_passband = (0:numel(s.lfm_passband)-1).' / s.cfg.Fs_used;

s.cfg.Nlfm = numel(s.lfm_passband);
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
