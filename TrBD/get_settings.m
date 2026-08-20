function s = get_settings(cfgTxt, Fs)

%% ================= CONFIGURATION ============================
cfg = struct();

cfg.Np = 60;
cfg.N_target_appearance = 20; % ping index

cfg.f0 = 1000;
cfg.f1 = 5000;
cfg.Tp = 30e-3;
cfg.c = 1500;
cfg.max_delay_spread = 0.1; %used for power calibration of interference

cfg.PRI = 1.4;

cfg.Fphi = 1-eps;
cfg.PFC = 1e-3;
cfg.Nmc = 128; % Monte Carlo Simulations

cfg.Fs = Fs;
cfg.baseband = false;
cfg.add_target = true;
cfg.add_background = true;

cfg.INRdB = 2; % multipath background power in [dB]

% Ocean/environment parameters
cfg.Ws_mps = 3; % wind speed [m/s]
cfg.TS_dB = -16; % Target strength
cfg.speedT = 5; % Speed for the AUV
cfg.pT0 = [50;90]; % initial target position
cfg.enable_TL = false; % true if target amplitude depends on propagation distance

cfg.overlap_frac = 0.2;

cfg.Tx_pos = [0; 0];
cfg.Rx_pos{1} = [1000; 0];
cfg.Rx_pos{2} = [500;-866.0254];
cfg.TL1_fixed = 0;

cfg.trgt_eta = -10; % signal power at the target point

% ---- Parse additional text parameters ---
cfg.fc = str2double(regexp(cfgTxt, 'fCarrier:\s*([\d.]+)', 'tokens', 'once'));
cfg.waymarkPeriod = str2double(regexp(cfgTxt, 'T_WAYMARK:\s*([\d.]+)', 'tokens', 'once'));
cfg.Timp_min = str2double(regexp(cfgTxt, 'T_IMPULSE_MIN:\s*([\d.]+)', 'tokens', 'once'));
tmp = regexp(cfgTxt, 'T_IMPULSE_NEG_START:\s*([\d.]+)', 'tokens', 'once');
cfg.Tneg_cfg = str2double(tmp);
if isempty(cfg.Tneg_cfg)
    cfg.Tneg_cfg = NaN;
end

cfg.initial_delay = norm(cfg.Rx_pos{1}-cfg.Tx_pos)/cfg.c;

% ---- Derived quantities ----
cfg.Nlfm = cfg.Fs * cfg.Tp;
cfg.Trec = 0.12; % processing-record duration [s]
cfg.L = floor( cfg.Trec* cfg.Fs);
cfg.N = cfg.Nlfm + cfg.L-1;
cfg.BW = abs(cfg.f1 - cfg.f0);
cfg.mu = cfg.BW / cfg.Tp; % chirp rate

cfg.t_ping = (0:cfg.Np-1).' * cfg.PRI;
cfg.tau_axis = (0:cfg.L-1)'/cfg.Fs + cfg.initial_delay;
s.t = (0:cfg.N-1).*1/cfg.Fs;

cfg.res_tau = 1/cfg.BW;
cfg.width_basis = cfg.res_tau;
cfg.basis_type = "uniform";

%% ================= TrBD FILTER =================
s.M = 1; % array elements (unused in this block)
s.Ns = 15000; % surviving particles
s.Nb = 15000; % birth particles (can be smaller than Ns)

% --- Existence model ---
s.pb = 1e-3; % birth probability per ping (tune 1e-4..1e-2)
s.ps = 1-eps; % survival probability per ping (0.98..0.999)
s.q0 = 1e-3; % initial existence probability
s.q_threshold = 0.96; % start with 0.9; use 0.95 once stable

% --- Process noise (per ping) for state [x; y; vx; vy; eta] ---
sigma_ax = 0.1; % m/s^2
sigma_ay = 0.1; % m/s^2
sigma_trgt_eta = 1; % eta-state process standard deviation
s.s2ax = sigma_ax^2; % (m/s^2)^2
s.s2ay = sigma_ay^2; % (m/s^2)^2
s.s2_trgt_eta = sigma_trgt_eta^2; % ^2

% --- Cartesian state regions (example; tune to your geometry) ---
s.x_region = [0, 500];
s.y_region = [-300, 100];
s.vx_region = [-5, 5]; % target is mostly perpendicular
s.vy_region = [-5, 5]; % known downward motion

s.tgt_eta_region = [-25, 5]; % allowed target amplitude for surviving particles
s.tgt_eta_region_new_targets = [-25, 5]; % births (can widen a bit if needed)
% --- Noise power (must match generator/likelihood) ---
s.s2e = 1e-2; % E{|n|^2} per complex sample

s.birth_preselect = true; % resample newly generated particles using the likelihood

if cfg.add_background
    s.bg_tracking = true; % enable/disable EKF background tracking
    if s.bg_tracking
        s.flt_nois_inf_const = 1; % measurement-noise inflation factor
    else
       s.flt_nois_inf_const = 1.5e2; % measurement-noise inflation factor when background tracking is disabled
    end
else
    s.bg_tracking = false;
    s.flt_nois_inf_const = 1; % measurement-noise inflation factor
end

%% Plotting
s.do_plot = true;
s.cfg = cfg;

end
