function snr_eff = background_tracking(Y, s, target_prop, irec)

cfg = s.cfg;
N = cfg.N;
Np = size(Y, 2);

S = s.S;
B = s.B;
U = s.U;
H = S * B;

K = size(B, 2);
I_N = eye(N, 'like', H);
I_K = eye(K, 'like', H);

target_signal = target_prop.target_signal;
snr_eff_k = zeros(1, Np);

% -----------------------------
% Dynamics
% -----------------------------
F = cfg.Fphi * speye(K);
Q = cfg.sigma_q(irec)^2 * diag(1./(1:K).^2);

% -----------------------------
% LS init
% -----------------------------
theta_prev = pinv(H) * Y(:, 1);
P_prev = 3e1 * I_K;

for k = 2:Np

    yk = Y(:, k);

    % -----------------------------
    % Predict
    % -----------------------------
    theta_pr = F * theta_prev;
    P_pr = F * P_prev * F' + Q;
    P_pr = 0.5 * (P_pr + P_pr');

    a_pr = B * theta_pr;

    % -----------------------------
    % Measurement covariance
    % -----------------------------
    Uc = U .* abs(a_pr(:)).';
    Ud = U * a_pr;

    Rk = s.s2e * I_N ...
       + cfg.sigma_c(irec)^2 * (Uc * Uc') ...
       + cfg.sigma_d(irec)^2 * (Ud * Ud');

    Rk = 0.5 * (Rk + Rk');

    Syy = H * P_pr * H' + Rk;
    Syy = 0.5 * (Syy + Syy');

    % -----------------------------
    % Cholesky
    % -----------------------------
    m = size(Syy, 1);
    jitter = 1e-12 + 1e-10 * real(trace(Syy)) / m;

    [LS, p] = chol(Syy + jitter * eye(m, 'like', Syy), 'lower');

    while p ~= 0
        jitter = jitter * 10;
        [LS, p] = chol(Syy + jitter * eye(m, 'like', Syy), 'lower');
    end

    % -----------------------------
    % Kalman update
    % -----------------------------
    Kgain = ((P_pr * H') / LS') / LS;
    innov = yk - H * theta_pr;

    theta_po = theta_pr + Kgain * innov;

    A = I_K - Kgain * H;
    P_po = A * P_pr * A' + Kgain * Rk * Kgain';
    P_po = 0.5 * (P_po + P_po');

    theta_prev = theta_po;
    P_prev = P_po;

    % -----------------------------
    % Effective SNR
    % -----------------------------
    sk = target_signal(:, k);

    if all(isfinite(sk)) && any(abs(sk) > 0)
        xk = LS \ sk;
        snr_eff_k(k) = real(xk' * xk);
    end

end

snr_eff = sum(snr_eff_k);

end
