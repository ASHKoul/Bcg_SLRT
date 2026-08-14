function ekf_init = background_intialization(y01, y02, use_bg, settings)

N   = numel(y01);
s2e = settings.cfg.sigma_e2;
I_N = eye(N, 'like', settings.H);

if use_bg

    H   = settings.H;
    Kbg = size(settings.B,2);

    lambda = 1e-8 * norm(H,'fro')^2 / Kbg;
    HH     = H' * H + lambda * eye(Kbg, 'like', H);

    th01 = HH \ (H' * y01);
    th02 = HH \ (H' * y02);

    ekf_init.theta_hat_1 = th01;
    ekf_init.theta_hat_2 = th02;
    ekf_init.P_hat_1     = 5e1 * eye(Kbg, 'like', H);
    ekf_init.P_hat_2     = 5e1 * eye(Kbg, 'like', H);
    ekf_init.H           = H;

    a01 = settings.B * th01;
    a02 = settings.B * th02;

    Uc1 = settings.U .* abs(a01(:)).';
    Uc2 = settings.U .* abs(a02(:)).';

    Ua1 = settings.U * a01;
    Ua2 = settings.U * a02;

    R_01 = s2e * eye(N, 'like', H) ...
         + settings.cfg.sigma_c(1)^2 * (Uc1 * Uc1') ...
         + settings.cfg.sigma_d(1)^2 * (Ua1 * Ua1');

    R_02 = s2e * eye(N, 'like', H) ...
         + settings.cfg.sigma_c(2)^2 * (Uc2 * Uc2') ...
         + settings.cfg.sigma_d(2)^2 * (Ua2 * Ua2');

    R_01 = 0.5 * (R_01 + R_01');
    R_02 = 0.5 * (R_02 + R_02');

    ekf_init.Sigma_1 = H * ekf_init.P_hat_1 * H' + R_01;
    ekf_init.Sigma_2 = H * ekf_init.P_hat_2 * H' + R_02;

    ekf_init.Sigma_1 = 0.5 * (ekf_init.Sigma_1 + ekf_init.Sigma_1');
    ekf_init.Sigma_2 = 0.5 * (ekf_init.Sigma_2 + ekf_init.Sigma_2');

    % ------------------------------------------------------------
    % Cholesky factors and log determinants for likelihood reuse
    % ------------------------------------------------------------

    jitter1 = 1e-12 + 1e-10 * real(trace(ekf_init.Sigma_1)) / max(N,1);
    [LS1,p1] = chol(ekf_init.Sigma_1, 'lower');

    while p1 ~= 0
        ekf_init.Sigma_1 = ekf_init.Sigma_1 + jitter1 * I_N;
        jitter1 = 10 * jitter1;
        [LS1,p1] = chol(ekf_init.Sigma_1, 'lower');
    end

    jitter2 = 1e-12 + 1e-10 * real(trace(ekf_init.Sigma_2)) / max(N,1);
    [LS2,p2] = chol(ekf_init.Sigma_2, 'lower');

    while p2 ~= 0
        ekf_init.Sigma_2 = ekf_init.Sigma_2 + jitter2 * I_N;
        jitter2 = 10 * jitter2;
        [LS2,p2] = chol(ekf_init.Sigma_2, 'lower');
    end

    ekf_init.LS_1 = LS1;
    ekf_init.LS_2 = LS2;

    ekf_init.logdet_1 = 2 * sum(log(real(diag(LS1))));
    ekf_init.logdet_2 = 2 * sum(log(real(diag(LS2))));

else

    ekf_init.theta_hat_1 = [];
    ekf_init.theta_hat_2 = [];
    ekf_init.P_hat_1     = [];
    ekf_init.P_hat_2     = [];
    ekf_init.H           = [];

    ekf_init.Sigma_1 = s2e;
    ekf_init.Sigma_2 = s2e;

    ekf_init.LS_1 = [];
    ekf_init.LS_2 = [];

    ekf_init.logdet_1 = [];
    ekf_init.logdet_2 = [];

end

end