function theta_true = dynamical_model_theta(theta, idx, H, s, irec)

K = length(theta);
Np = s.cfg.Np;

idx = idx(:);

sigma_q_idx = s.cfg.sigma_q(irec) ./ idx;
Eref = sum(abs(H * theta).^2);

theta_true = zeros(K, Np, 'like', theta);

theta_true(idx, 1) = theta(idx);
theta_true(:, 1) = normalize_energy(theta_true(:, 1), H, Eref);

for k = 2:Np
    theta_true(idx, k) = s.cfg.Fphi * theta_true(idx, k-1) ...
        + sigma_q_idx .* (randn(numel(idx), 1) + 1j*randn(numel(idx), 1)) / sqrt(2);

    theta_true(:, k) = normalize_energy(theta_true(:, k), H, Eref);
end

end

function theta = normalize_energy(theta, H, Eref)

E = sum(abs(H * theta).^2);

if E > 0 && Eref > 0
    theta = theta * sqrt(Eref / E);
end

end
