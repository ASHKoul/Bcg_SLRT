function cfg =generate_tau_beta(cfg)
if cfg.add_target
    switch cfg.target_status
        case "fixed"
            T  = cfg.PRI;
            k0 = cfg.Ntrain + 1;

            % allocate ping-wise outputs
            cfg.tau_target  = NaN(cfg.Np,1);
            cfg.beta_target = NaN(cfg.Np,1);

            % initial state at ping k0
            pT = [1500; 300];     % [2x1] position (m)
            vT = [0; 0];          % [2x1] velocity (m/s)

            % accel std (m/s^2): set to 0 for constant velocity
            sigma_a = 0;          % <-- change this as you like

            pTx = cfg.Tx_pos(:);
            pRx = cfg.Rx_pos(:);

            for k = k0:cfg.Np

                % random acceleration (2D)
                aT = sigma_a * randn(2,1);

                % propagate
                if k > k0
                    pT = pT + vT*T + 0.5*aT*T^2;
                    vT = vT + aT*T;
                end

                % geometry
                dTx = pT - pTx;
                rTx = norm(dTx);
                dRx = pT - pRx;
                rRx = norm(dRx);

                tau = (rTx + rRx) / cfg.c;

                uTx  = dTx / max(rTx, eps);
                uRx  = dRx / max(rRx, eps);
                Rdot = dot(uTx, vT) + dot(uRx, vT);

                beta = 1 - (Rdot / cfg.c);

                cfg.tau_target(k)  = tau;
                cfg.beta_target(k) = beta;
            end

        case "moving"
            T  = cfg.PRI;
            k0 = cfg.Ntrain + 1;

            cfg.tau_target  = NaN(cfg.Np,1);
            cfg.beta_target = NaN(cfg.Np,1);

            % -------- YOUR DESIGN CHOICES --------
            Kact = 60;                         % only first 60 pings after k0 are "active"
            pT0  = [1500; 300];                % target position at ping k0 (arrival)

            kend = min(cfg.Np, k0 + Kact - 1);

            % crossing at half of the active pings
            % (for Kact=60 => half is 30 pings after start, so k_cross = k0+29)
            k_cross = k0 + floor((Kact-1)/2);

            if k_cross < k0 || k_cross > kend
                error('k_cross=%d must be within active window [%d,%d].', k_cross, k0, kend);
            end

            pTx = cfg.Tx_pos(:);
            pRx = cfg.Rx_pos(:);

            % Baseline unit vectors
            dTR = pRx - pTx;
            D   = norm(dTR);
            u   = dTR / max(D, eps);       % along baseline
            v   = [-u(2); u(1)];           % perpendicular (2D)

            % Express pT0 in baseline coordinates
            r0 = pT0(:) - pTx;
            % x0 = dot(r0, u);               % along-baseline coordinate (constant in this model)
            y0 = dot(r0, v);               % signed perpendicular offset (we drive this to 0)

            % Velocity: purely perpendicular so that y hits 0 at k_cross
            t_to_cross = (k_cross - k0) * T;
            if abs(y0) < 1e-12 || t_to_cross <= 0
                vT = [0; 0];               % already on the line (or degenerate)
            else
                vy_perp = -y0 / t_to_cross; % m/s along v direction
                vT      = vy_perp * v;      % 2x1 velocity vector
            end

            % Initialize at arrival ping
            pT = pT0(:);

            % ---- Compute tau/beta ONLY over active window ----
            for k = k0:kend
                if k > k0
                    pT = pT + vT*T;
                end

                dTx = pT - pTx;  rTx = norm(dTx);
                dRx = pT - pRx;  rRx = norm(dRx);

                tau = (rTx + rRx) / cfg.c;

                uTx  = dTx / max(rTx, eps);
                uRx  = dRx / max(rRx, eps);
                Rdot = dot(uTx, vT) + dot(uRx, vT);

                beta = 1 - (Rdot / cfg.c);

                cfg.tau_target(k)  = tau;
                cfg.beta_target(k) = beta;
            end

            % ---- Hold last tau/beta constant after active window ----
            cfg.tau_target(kend+1:cfg.Np)  = cfg.tau_target(kend);
            cfg.beta_target(kend+1:cfg.Np) = cfg.beta_target(kend);


        case "fixed_moving"
            T  = cfg.PRI; %#ok<NASGU>
            k0 = cfg.Ntrain + 1;

            % allocate ping-wise outputs
            cfg.tau_target  = NaN(cfg.Np,1);
            cfg.beta_target = NaN(cfg.Np,1);

            % --- fixed target position (same as your moving start) ---
            pT = [1500; 300];        % [2x1] position (m), fixed for all k>=k0

            % --- desired bistatic range-rate (m/s) ---
            Rdot_des = 5;            % bistatic velocity / range-rate

            pTx = cfg.Tx_pos(:);
            pRx = cfg.Rx_pos(:);

            % geometry at the fixed position
            dTx = pT(:) - pTx;  rTx = norm(dTx);
            dRx = pT(:) - pRx;  rRx = norm(dRx);

            tau_const = (rTx + rRx) / cfg.c;

            uTx = dTx / max(rTx, eps);
            uRx = dRx / max(rRx, eps);

            % Choose a velocity vector vT such that:
            % Rdot = dot(uTx,vT) + dot(uRx,vT) = Rdot_des
            s = uTx + uRx;                 % 2x1
            denom = dot(s, s);             % ||s||^2

            vT = (Rdot_des / denom) * s;
            Rdot = dot(uTx, vT) + dot(uRx, vT);   % ~ Rdot_des


            beta_const = 1 - (Rdot / cfg.c);

            % fill from arrival ping onward
            cfg.tau_target(k0:cfg.Np)  = tau_const;
            cfg.beta_target(k0:cfg.Np) = beta_const;

        case "moving_block"
            T  = cfg.PRI;
            k0 = cfg.Ntrain + 1;

            cfg.tau_target  = NaN(cfg.Np,1);
            cfg.beta_target = NaN(cfg.Np,1);

            % -------- DESIGN CHOICES --------
            Kact = 60;                 % only first 60 pings after k0 are "active"
            Lb   = 15;                 % block length for piecewise-stationary tau/beta
            pT0  = [1500; 300];        % target position at ping k0 (arrival)

            kend = min(cfg.Np, k0 + Kact - 1);

            % crossing at half of the active pings (e.g., Ntrain=40 => k0=41 => k_cross=70)
            k_cross = k0 + floor((Kact - 1)/2);
            if k_cross < k0 || k_cross > kend
                error('k_cross=%d must be within active window [%d,%d].', k_cross, k0, kend);
            end

            % geometry anchors
            pTx = cfg.Tx_pos(:);
            pRx = cfg.Rx_pos(:);

            % Baseline unit vectors (2D)
            dTR = pRx - pTx;
            D   = norm(dTR);
            u   = dTR / max(D, eps);       % along baseline
            v   = [-u(2); u(1)];           % perpendicular (2D)

            % Express pT0 in baseline coordinates
            r0 = pT0(:) - pTx;
            y0 = dot(r0, v);               % signed perpendicular offset (drive to 0)

            % Velocity: purely perpendicular so that y hits 0 at k_cross
            t_to_cross = (k_cross - k0) * T;
            if abs(y0) < 1e-12 || t_to_cross <= 0
                vT = [0; 0];
            else
                vy_perp = -y0 / t_to_cross;  % m/s along v
                vT      = vy_perp * v;       % 2x1 velocity vector
            end

            % Initialize at arrival ping
            pT = pT0(:);

            % ---- Compute tau/beta ONLY over active window (continuous) ----
            for k = k0:kend
                if k > k0
                    pT = pT + vT*T;   % constant velocity
                end

                dTx = pT - pTx;  rTx = norm(dTx);
                dRx = pT - pRx;  rRx = norm(dRx);

                tau = (rTx + rRx) / cfg.c;

                uTx  = dTx / max(rTx, eps);
                uRx  = dRx / max(rRx, eps);
                Rdot = dot(uTx, vT) + dot(uRx, vT);

                beta = 1 - (Rdot / cfg.c);

                cfg.tau_target(k)  = tau;
                cfg.beta_target(k) = beta;
            end

            % ---- Quantize tau/beta in 10-ping blocks within the 60 active pings ----
            idxAct = k0:kend;
            Kact_eff = numel(idxAct);
            if Kact_eff < 1
                error('Active window is empty: k0=%d, kend=%d', k0, kend);
            end

            for b0 = 1:Lb:Kact_eff
                b1 = min(b0 + Lb - 1, Kact_eff);
                kk = idxAct(b0:b1);

                % Use FIRST ping in each block as representative (piecewise-stationary)
                tau_rep  = cfg.tau_target(kk(1));
                beta_rep = cfg.beta_target(kk(1));

                cfg.tau_target(kk)  = tau_rep;
                cfg.beta_target(kk) = beta_rep;
            end

            % ---- Hold last tau/beta constant after the active window ----
            cfg.tau_target(kend+1:cfg.Np)  = cfg.tau_target(kend);
            cfg.beta_target(kend+1:cfg.Np) = cfg.beta_target(kend);





    end
end