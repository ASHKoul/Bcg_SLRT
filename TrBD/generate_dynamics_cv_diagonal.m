function target_prop = generate_dynamics_cv_diagonal(Yb, model, s, irec)

target_prop = s.cfg;

k0 = target_prop.N_target_appearance;
Np = target_prop.Np;
N = target_prop.N;

target_prop.x_target = NaN(Np, 1);
target_prop.y_target = NaN(Np, 1);
target_prop.vx_target = NaN(Np, 1);
target_prop.vy_target = NaN(Np, 1);
target_prop.v_target = NaN(2, Np);

target_prop.tau_target = NaN(Np, 1);
target_prop.beta_target = NaN(Np, 1);
target_prop.rTx_target = NaN(Np, 1);
target_prop.rRx_target = NaN(Np, 1);
target_prop.rBi_target = NaN(Np, 1);
target_prop.Es_target = zeros(Np, 1);
target_prop.amp_target_true = NaN(Np, 1);
target_prop.eta_target = NaN(Np, 1);

target_prop.TL1_target = NaN(Np, 1); % Tx -> target [dB]
target_prop.TL2_target = NaN(Np, 1); % target -> Rx [dB]
target_prop.eta_eff_target = NaN(Np, 1); % eta after TL1 folding [dB]

target_prop.target_signal = zeros(N, Np);

% ------------------------------------------------------------
% Initial target state: [x; y; vx; vy; eta]
% ------------------------------------------------------------
vT0 = target_prop.speedT * [cosd(-45); sind(-45)];

xT = [target_prop.pT0(:); ...
      vT0(:); ...
      target_prop.trgt_eta];

if k0 <= Np
    for k = k0:Np

        if k > k0
            xT = model.F * xT;
        end

        pT = xT(1:2);
        vT = xT(3:4);
        eta = xT(5);

        target_prop.x_target(k) = pT(1);
        target_prop.y_target(k) = pT(2);
        target_prop.vx_target(k) = vT(1);
        target_prop.vy_target(k) = vT(2);
        target_prop.v_target(:, k) = vT;

        [tau, beta, geom] = cart_to_tau_beta(pT, vT, target_prop, irec);

        target_prop.rTx_target(k) = geom.rTx;
        target_prop.rRx_target(k) = geom.rRx;
        target_prop.rBi_target(k) = geom.rBi;
        target_prop.tau_target(k) = tau;
        target_prop.beta_target(k) = beta;
        target_prop.eta_target(k) = eta;

        % --------------------------------------------------------
        % Split propagation loss:
        % TL1 = Tx -> target  : fold into eta
        % TL2 = target -> Rx  : apply separately in received amplitude
        % --------------------------------------------------------
        if target_prop.enable_TL
            TL1_k = 17 * log10(geom.rTx);
            TL2_k = 17 * log10(geom.rRx);
        else
            TL1_k = target_prop.TL1_fixed;
            TL2_k = 0;
        end

        eta_eff = eta - TL1_k; % eta after Tx->target loss, before target->Rx loss
        Gsig_k = 10^(-TL2_k/10); % propagation from target to receiver

        target_prop.TL1_target(k) = TL1_k;
        target_prop.TL2_target(k) = TL2_k;
        target_prop.eta_eff_target(k) = eta_eff;

        s.amplitude = sqrt(10^(eta_eff/10) * Gsig_k);

        [yk, target_prop.amp_target_true(k)] = target_signal_genration(s, tau, beta);

        target_prop.target_signal(:, k) = yk(:);
        target_prop.Es_target(k) = sum(abs(yk).^2);
    end
end

if target_prop.add_background
    target_prop.SNR_eff = background_tracking(Yb, s, target_prop, irec);
else
    target_prop.SNR_eff = sum(target_prop.Es_target, 'omitnan') / (s.s2e);
end

end
