function y = generate_synthetic_data_MC(s, irec, theta_true)

H = s.S * s.B;
N = size(H, 1);
Np = size(theta_true, 2);

a = s.B * theta_true; % L x Np
Ua = s.U * a; % N x Np

y_clean = H * theta_true; % N x Np

Ze = (randn(N, Np) + 1j*randn(N, Np)) / sqrt(2);
Zc = (randn(size(a)) + 1j*randn(size(a))) / sqrt(2);
Zd = (randn(1, Np) + 1j*randn(1, Np)) / sqrt(2);

ne = sqrt(s.s2e) * Ze;

nc = s.cfg.sigma_c(irec) * s.U * (abs(a) .* Zc);

nd = s.cfg.sigma_d(irec) * (Ua .* Zd);

y = y_clean + ne + nc + nd;

end
