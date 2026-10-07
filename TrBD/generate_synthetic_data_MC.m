function y = generate_synthetic_data_MC(s, irec, theta_true)

H = s.S * s.B;
N = size(H, 1);
Np = size(theta_true, 2);

a = s.B * theta_true; % L x Np
Ua = s.U * a; % N x Np

y_clean = H * theta_true; % N x Np

Ze = randn(N, Np);
Zc = randn(size(a));
Zd = randn(1, Np);

ne = sqrt(s.s2e) * Ze;

nc = s.cfg.sigma_c(irec) * s.U * (abs(a) .* Zc);

nd = s.cfg.sigma_d(irec) * (Ua .* Zd);

y = y_clean + ne + nc + nd;

end
