function data_mat = signal_waveform(Yb_syn, model, s)

N_rec = numel(Yb_syn);
Nmc = s.cfg.Nmc;

data_mat.Yt_clean = cell(Nmc, N_rec);
data_mat.target_true_amplitude = cell(Nmc, N_rec);
data_mat.SNR_eff = zeros(1, N_rec);

for irec = 1:N_rec

    s.t_rec = (0:s.cfg.N-1).' / s.cfg.Fs_used;

    if s.cfg.add_target

        target_prop = generate_dynamics_cv_diagonal(Yb_syn{irec}, model, s, irec);

    else

        Np = s.cfg.Np;
        N = s.cfg.N;
        target_prop = s.cfg;
        target_prop.x_target = NaN(Np, 1);
        target_prop.y_target = NaN(Np, 1);
        target_prop.vx_target = NaN(Np, 1);
        target_prop.vy_target = NaN(Np, 1);
        target_prop.v_target = NaN(2, Np);
        target_prop.a_target_true = NaN(Np, 1);
        target_prop.amp_target_true = NaN(Np, 1);
        target_prop.tau_target = NaN(Np, 1);
        target_prop.beta_target = NaN(Np, 1);
        target_prop.rTx_target = NaN(Np, 1);
        target_prop.rRx_target = NaN(Np, 1);
        target_prop.rBi_target = NaN(Np, 1);
        target_prop.Es_target = zeros(Np, 1);
        target_prop.eta_target = NaN(Np, 1);

        target_prop.target_signal = zeros(N, Np);
        target_prop.SNR_eff = 0;

    end

    data_mat.Yt_clean(:, irec) = repmat({target_prop.target_signal}, Nmc, 1);

    data_mat.target_true_amplitude(:, irec) = repmat({target_prop.amp_target_true}, Nmc, 1);

    data_mat.SNR_eff(irec) = target_prop.SNR_eff;
    s.eta_target{irec} = target_prop.eta_target;
    s.target_prop{irec} = target_prop;
    data_mat.s = s;

end

end
