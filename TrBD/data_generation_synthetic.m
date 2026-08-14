function [Ybsim, theta_true, s_new] = data_generation_synthetic(data, s)

N_Rec = numel(data);
Nmc = s.cfg.Nmc;
Np = s.cfg.Np;

Ybsim = cell(Nmc, N_Rec);
theta_true = cell(1, N_Rec);

for irec = 1:N_Rec

    %% ============================================================
    % Receiver-specific preprocessing
    % ============================================================

    [Yb, Atrue, s_new] = mult_background_processing(data{irec}, s);
    [~, Atrue, s_new] = power_calibrations(Yb, Atrue, s_new);
    s_new = generate_basis_cov_matrix(s_new);

    N = s_new.cfg.N;
    sigma_e2 = s_new.s2e;

    Ytmp = cell(Nmc, 1);

    %% ============================================================
    % Generate background and noise
    % ============================================================

    if s.cfg.add_background

        B = s_new.B;
        H = s_new.S * B;
        K = size(B, 2);

        lambda = 1e-8 * norm(B, 'fro')^2 / K;
        BtB = B' * B + lambda * eye(K, 'like', B);

        theta_1 = BtB \ (B' * Atrue(:, 1));

        idx_theta = index_finding_theta(theta_1);

        theta_true{irec} = dynamical_model_theta(theta_1, idx_theta, H, s_new, irec);

        parfor inmc = 1:Nmc
            Ytmp{inmc} = generate_synthetic_data_MC(s_new, irec, theta_true{irec}); %#ok<PFBNS>
        end

    else

        %% ============================================================
        % Generate noise-only data
        % ============================================================

        theta_true{irec} = [];

        parfor inmc = 1:Nmc
            Ytmp{inmc} = sqrt(0.5 * sigma_e2) * (randn(N, Np) + 1j * randn(N, Np));
        end

    end

    Ybsim(:, irec) = Ytmp;

end

end
