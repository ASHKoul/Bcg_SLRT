function s = generate_basis_cov_matrix(s)

% ---------------------------------------------------------
% Read settings
% ---------------------------------------------------------

Ts = 1/ s.cfg.Fs_used;
L = s.cfg.L;
N = s.cfg.N;
Nlfm = s.cfg.Nlfm;
tau_axis = s.cfg.tau_axis(:);

% ---------------------------------------------------------
% Waveform on processing grid
% ---------------------------------------------------------
sig = s.lfm_bb(:);
if numel(sig) < Nlfm
    sig(end+1:Nlfm) = 0;
else
    sig = sig(1:Nlfm);
end

% ---------------------------------------------------------
% Derivative-based waveform u
% u(t) = (t-t0) s_dot(t)
% ---------------------------------------------------------
tp = (0:Nlfm-1).' * Ts;
t0 = (Nlfm-1) * Ts / 2;

ds = gradient(sig, Ts);
u = (tp - t0) .* ds;

% ---------------------------------------------------------
% Convolution matrices
% ---------------------------------------------------------
S = convmtx(sig, L);
U = convmtx(u, L);

s.S = S(1:N, 1:L);
s.U = U(1:N, 1:L);

basis_type = string(s.cfg.basis_type);

% ---------------------------------------------------------
% Fixed Gaussian width
% ---------------------------------------------------------
sigma = s.cfg.width_basis ;
d_ref = s.cfg.width_basis;
% ---------------------------------------------------------
% Basis centers
% ---------------------------------------------------------
switch basis_type

    case "uniform"
        % centers are uniformly spaced in linear delay
        tau_lo = tau_axis(1) + 0.5*Ts;
        tau_hi = tau_axis(end) - 0.5*Ts;

        mu = (tau_lo : d_ref : tau_hi).';
        if isempty(mu)
            mu = mean(tau_axis);
        end

    case "log"

        K = min(350, L);

        tau_min = Ts/2;
        tau_max = max(tau_min, (L - 0.5) * Ts);

        if K == 1
            mu_rel = sqrt(tau_min * tau_max);
        else
            mu_rel = exp(linspace(log(tau_min), log(tau_max), K)).';
        end

        mu = tau_axis(1) - 0.5*Ts + mu_rel;

        if isempty(mu)
            mu = mean(tau_axis);
        end

    otherwise
        error('Unknown basis_type: %s. Use "uniform" or "log".', basis_type);
end

% ---------------------------------------------------------
% Build Gaussian basis with fixed sigma
% B is [delay x basis-index]
% ---------------------------------------------------------
D = tau_axis - mu.';
B = exp(-0.5 * (abs(D) / sigma).^2);

% Robust column normalization
col_norm = sqrt(sum(abs(B).^2, 1));
B = B ./ col_norm;

s.B = B;
s.cfg.basis_centers = mu;
s.cfg.basis_sigma = sigma;
end
