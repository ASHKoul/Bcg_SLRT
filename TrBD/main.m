clc;
clearvars;
% close all;

%% ============================================================
% Load settings and data
% ============================================================

data_background = cell(1, 2);

data_background{1} = load("Bellhop_iso50_lfm_0030_14_20260216T075948.mat");
data_background{2} = load("Bellhop_R2_lfm_0030_14_20260306T160306.mat");

cfgTxt = fileread("config_iso50_lfm_0030_14_20260216T075948.txt");
[sig, Fs] = audioread("lfm_0030_14.wav");

s = get_settings(cfgTxt, Fs);
s.lfm = sig(1:s.cfg.Nlfm);

N_Rec = numel(s.cfg.Rx_pos);
Nmc = s.cfg.Nmc;

%% ============================================================
% Target dynamic model
% ============================================================

TPRI = s.cfg.PRI;

model.F = [1 0 TPRI 0 0;
           0 1 0 TPRI 0;
           0 0 1 0 0;
           0 0 0 1 0;
           0 0 0 0 1];

model.G = [TPRI^2/2 0 0;
           0 TPRI^2/2 0;
           TPRI 0 0;
           0 TPRI 0;
           0 0 1];

model.Q = diag([s.s2ax, s.s2ay, s.s2_trgt_eta]);

%% ============================================================
% Load generation hyperparameters
% ============================================================

hp_1 = load("learn_hyperparams_uniform_receiver_1_INR_2p00dB_20260408_213113.mat");
hp_2 = load("learn_hyperparams_uniform_receiver_2_INR_2p00dB_20260408_205310.mat");

r1 = hp_1.results_hp.results_models.Mcd;
r2 = hp_2.results_hp.results_models.Mcd;

s.cfg.sigma_q = zeros(1, 2);
s.cfg.sigma_c = zeros(1, 2);
s.cfg.sigma_d = zeros(1, 2);

s.cfg.sigma_q(1) = r1.sigma_q_hat;
s.cfg.sigma_c(1) = r1.sigma_c_hat;
s.cfg.sigma_d(1) = r1.sigma_d_hat;

s.cfg.sigma_q(2) = r2.sigma_q_hat;
s.cfg.sigma_c(2) = r2.sigma_c_hat;
s.cfg.sigma_d(2) = r2.sigma_d_hat;

%% ============================================================
% Generate synthetic background
% ============================================================

[Yb_syn, theta_true, s] = data_generation_synthetic(data_background, s);

%% ============================================================
% Generate target waveform
% ============================================================

data_target = signal_waveform(Yb_syn(1, :), model, s);

%% ============================================================
% Generate combined received data
% ============================================================

[Ymat, data_target.s] = generate_received_signal_with_background(Yb_syn, data_target.Yt_clean, data_target.s);

%% ============================================================
% Diagnostic plots
% ============================================================

plot_data_diagonsitics(Yb_syn{1, 1}, data_target, Ymat{1, 1}, Ymat{1, 2});

%% ============================================================
% Configure tracking filter
% ============================================================

s_flt = settings_for_filter(data_target.s, model);

s_flt.cfg.sigma_q(1) = r1.sigma_q_hat;
s_flt.cfg.sigma_c(1) = 0.9 * r1.sigma_c_hat;
s_flt.cfg.sigma_d(1) = 0.9 * r1.sigma_d_hat;

s_flt.cfg.sigma_q(2) = r2.sigma_q_hat;
s_flt.cfg.sigma_c(2) = 0.9 * r2.sigma_c_hat;
s_flt.cfg.sigma_d(2) = 0.9 * r2.sigma_d_hat;

qshape = (1:size(s_flt.B, 2)).'.^(-2);

s_flt.Qbg_1 = s_flt.cfg.sigma_q(1)^2 * qshape;
s_flt.Qbg_2 = s_flt.cfg.sigma_q(2)^2 * qshape;

s_flt.I_N = eye(size(s_flt.H, 1), 'like', s_flt.H);

%% ============================================================
% Run tracker
% ============================================================

x_out = cell(Nmc, 1);
q_out = cell(Nmc, 1);

for inmc = 1:Nmc

    data = tbd_bis_active(Ymat{inmc, 1}, Ymat{inmc, 2}, model, s_flt);

    x_out{inmc} = data.xh;
    q_out{inmc} = data.q;

end

%% ============================================================
% Store output
% ============================================================

trbd_filter_out.x = x_out;
trbd_filter_out.q = q_out;
