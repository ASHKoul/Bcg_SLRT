function cfg = generate_tau_beta(cfg)

if ~isfield(cfg,'add_target') || ~cfg.add_target
    return;
end

switch cfg.target_status
  
    case "fixed"
    % Stationary target (vT = 0) => constant tau/beta after arrival.

        T  = cfg.PRI;
        k0 = cfg.Ntrain + 1;

        cfg.tau_target  = NaN(cfg.Np, 1);
        cfg.beta_target = NaN(cfg.Np, 1);

        % initial state at ping k0
        pT = [1500; 300];   % [m]
        vT = [0; 0];        % [m/s]

        % accel std (m/s^2): set to 0 for constant velocity
        sigma_a = 0;

        pTx = cfg.Tx_pos(:);
        pRx = cfg.Rx_pos(:);

        for k = k0:cfg.Np

            aT = sigma_a * randn(2,1);

            if k > k0
                pT = pT + vT*T + 0.5*aT*T^2;
                vT = vT + aT*T;
            end

            dTx = pT - pTx;  rTx = norm(dTx);
            dRx = pT - pRx;  rRx = norm(dRx);

            tau  = (rTx + rRx) / cfg.c;

            uTx  = dTx / max(rTx, eps);
            uRx  = dRx / max(rRx, eps);
            Rdot = dot(uTx, vT) + dot(uRx, vT);

            beta = 1 - (Rdot / cfg.c);

            cfg.tau_target(k)  = tau;
            cfg.beta_target(k) = beta;
        end

    case "moving"
    % Target moves perpendicular to Tx-Rx baseline so it crosses the baseline
    % at k_cross, then tau/beta are held constant after active window.

        T  = cfg.PRI;
        k0 = cfg.Ntrain + 1;

        cfg.tau_target  = NaN(cfg.Np, 1);
        cfg.beta_target = NaN(cfg.Np, 1);

        % -------- DESIGN CHOICES --------
        Kact = 60;              % active pings starting at k0
        pT0  = [1500; 300];     % position at ping k0

        kend    = min(cfg.Np, k0 + Kact - 1);
        k_cross = k0 + floor((Kact - 1)/2);

        if k_cross < k0 || k_cross > kend
            error('k_cross=%d must be within active window [%d,%d].', k_cross, k0, kend);
        end

        pTx = cfg.Tx_pos(:);
        pRx = cfg.Rx_pos(:);

        % Baseline unit vectors
        dTR = pRx - pTx;
        D   = norm(dTR);
        u   = dTR / max(D, eps);      % along baseline
        v   = [-u(2); u(1)];          % perpendicular in 2D

        % Perpendicular offset we want to drive to zero at k_cross
        r0 = pT0(:) - pTx;
        y0 = dot(r0, v);

        % Velocity purely along v to hit y=0 at k_cross
        t_to_cross = (k_cross - k0) * T;
        if abs(y0) < 1e-12 || t_to_cross <= 0
            vT = [0; 0];
        else
            vy_perp = -y0 / t_to_cross;
            vT      = vy_perp * v;
        end

        pT = pT0(:);

        % ---- Compute tau/beta over active window ----
        for k = k0:kend

            if k > k0
                pT = pT + vT*T;
            end

            dTx = pT - pTx;  rTx = norm(dTx);
            dRx = pT - pRx;  rRx = norm(dRx);

            tau  = (rTx + rRx) / cfg.c;

            uTx  = dTx / max(rTx, eps);
            uRx  = dRx / max(rRx, eps);
            Rdot = dot(uTx, vT) + dot(uRx, vT);

            beta = 1 - (Rdot / cfg.c);

            cfg.tau_target(k)  = tau;
            cfg.beta_target(k) = beta;
        end

        % ---- Hold last value constant after active window ----
        if kend < cfg.Np
            cfg.tau_target(kend+1:cfg.Np)  = cfg.tau_target(kend);
            cfg.beta_target(kend+1:cfg.Np) = cfg.beta_target(kend);
        end

        % NOTE: cfg.tau_target(1:k0-1) remains NaN (by design).

    case "fixed_moving"
    % Fixed position but enforce a nonzero bistatic range-rate by choosing vT
    % aligned with (uTx + uRx). Tau is constant; beta constant.

        k0 = cfg.Ntrain + 1;

        cfg.tau_target  = NaN(cfg.Np, 1);
        cfg.beta_target = NaN(cfg.Np, 1);

        pT = [1500; 300];     % fixed target position

        % desired bistatic range-rate (m/s)
        Rdot_des = 5;

        pTx = cfg.Tx_pos(:);
        pRx = cfg.Rx_pos(:);

        dTx = pT(:) - pTx;  rTx = norm(dTx);
        dRx = pT(:) - pRx;  rRx = norm(dRx);

        tau_const = (rTx + rRx) / cfg.c;

        uTx = dTx / max(rTx, eps);
        uRx = dRx / max(rRx, eps);

        s     = uTx + uRx;          % direction that changes bistatic range
        denom = dot(s, s);

        % -------------------- IMPORTANT POTENTIAL ISSUE --------------------
        % If denom ~ 0 (rare geometry where uTx ≈ -uRx), this divides by ~0.

        if denom < 1e-12
            error('Degenerate geometry: ||uTx+uRx||^2 is ~0, cannot set Rdot_des.');
        end

        vT   = (Rdot_des / denom) * s;
        Rdot = dot(uTx, vT) + dot(uRx, vT);

        beta_const = 1 - (Rdot / cfg.c);

        cfg.tau_target(k0:cfg.Np)  = tau_const;
        cfg.beta_target(k0:cfg.Np) = beta_const;

        % NOTE: cfg.tau_target(1:k0-1) remains NaN (by design).

    case "moving_block"
    % Same as "moving" but quantize tau/beta into piecewise-constant blocks

        T  = cfg.PRI;
        k0 = cfg.Ntrain + 1;

        cfg.tau_target  = NaN(cfg.Np, 1);
        cfg.beta_target = NaN(cfg.Np, 1);

        % -------- DESIGN CHOICES --------
        Kact = 60;              % active pings
        Lb   = 10;              % block length (pings)
        pT0  = [1500; 300];

        kend    = min(cfg.Np, k0 + Kact - 1);
        k_cross = k0 + floor((Kact - 1)/2);

        if k_cross < k0 || k_cross > kend
            error('k_cross=%d must be within active window [%d,%d].', k_cross, k0, kend);
        end

        pTx = cfg.Tx_pos(:);
        pRx = cfg.Rx_pos(:);

        dTR = pRx - pTx;
        D   = norm(dTR);
        u   = dTR / max(D, eps);
        v   = [-u(2); u(1)];

        r0 = pT0(:) - pTx;
        y0 = dot(r0, v);

        t_to_cross = (k_cross - k0) * T;
        if abs(y0) < 1e-12 || t_to_cross <= 0
            vT = [0; 0];
        else
            vy_perp = -y0 / t_to_cross;
            vT      = vy_perp * v;
        end

        pT = pT0(:);

        % ---- Continuous tau/beta over active window ----
        for k = k0:kend

            if k > k0
                pT = pT + vT*T;
            end

            dTx = pT - pTx;  rTx = norm(dTx);
            dRx = pT - pRx;  rRx = norm(dRx);

            tau  = (rTx + rRx) / cfg.c;

            uTx  = dTx / max(rTx, eps);
            uRx  = dRx / max(rRx, eps);
            Rdot = dot(uTx, vT) + dot(uRx, vT);

            beta = 1 - (Rdot / cfg.c);

            cfg.tau_target(k)  = tau;
            cfg.beta_target(k) = beta;
        end

        % ---- Quantize tau/beta in blocks within the active window ----
        idxAct   = k0:kend;
        Kact_eff = numel(idxAct);

        if Kact_eff < 1
            error('Active window is empty: k0=%d, kend=%d', k0, kend);
        end

        for b0 = 1:Lb:Kact_eff
            b1 = min(b0 + Lb - 1, Kact_eff);
            kk = idxAct(b0:b1);

            % Representative = first ping in block
            tau_rep  = cfg.tau_target(kk(1));
            beta_rep = cfg.beta_target(kk(1));

            cfg.tau_target(kk)  = tau_rep;
            cfg.beta_target(kk) = beta_rep;
        end

        % ---- Hold last value constant after active window ----
        if kend < cfg.Np
            cfg.tau_target(kend+1:cfg.Np)  = cfg.tau_target(kend);
            cfg.beta_target(kend+1:cfg.Np) = cfg.beta_target(kend);
        end

        % NOTE: cfg.tau_target(1:k0-1) remains NaN (by design).

    %% ====================================================================
    otherwise
        error('Unknown cfg.target_status = "%s".', string(cfg.target_status));
end
end

