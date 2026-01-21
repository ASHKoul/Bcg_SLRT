function cfg = settings_data(cfgTxt, Fs)

%% ========================================================================
%  DEFAULTS (your experiment settings)
%  ========================================================================

cfg = struct();

% --- experiment length ---
cfg.Np     = 100;    % number of pings
cfg.Ntrain = 40;     % training pings (pre-change), so k0 = Ntrain+1

% --- waveform / propagation ---
cfg.f0        = 1000;      % [Hz]
cfg.f1        = 5000;      % [Hz]
cfg.Tp        = 25e-3;     % [s]
cfg.c         = 1500;      % [m/s]
cfg.max_range = 2000;      % [m]
cfg.max_delay = 90e-3;     % [s]
cfg.PRI       = 0.12;      % [s]
cfg.Fphi      = 1 - eps;   % state transition scale (close to 1)

% --- detection / Monte-Carlo ---
cfg.Pfa = 0.05;      % false alarm probability
cfg.Nmc = 10;        % Monte Carlo repetitions

% --- signal handling ---
cfg.Fs        = Fs;        % [Hz]
cfg.baseband  = true;
cfg.add_noise = true;
cfg.add_target= true;
cfg.target_status = "fixed";

% --- scenario settings ---
cfg.CNRdB      = 30;
cfg.SNRdB_list = 0:5:20;

% --- basis design ---
cfg.overlap_frac = 0.2;

% --- geometry (2D) ---
cfg.Tx_pos = [0; 0];        % [m] transmitter position (2x1)
cfg.Rx_pos = [2000; 0];     % [m] receiver position (2x1)

%% ========================================================================
%  PARSE FIELDS FROM TEXT (cfgTxt)
%  ========================================================================

cfg.fc = str2double(regexp(cfgTxt, 'fCarrier:\s*([\d.]+)', 'tokens', 'once'));

cfg.waymarkPeriod = str2double(regexp(cfgTxt, 'T_WAYMARK:\s*([\d.]+)', 'tokens', 'once'));
cfg.Timp_min      = str2double(regexp(cfgTxt, 'T_IMPULSE_MIN:\s*([\d.]+)', 'tokens', 'once'));

tmp = regexp(cfgTxt, 'T_IMPULSE_NEG_START:\s*([\d.]+)', 'tokens', 'once');
cfg.Tneg_cfg = str2double(tmp);


%% ========================================================================
%  DERIVED PARAMETERS
%  ========================================================================

cfg.initial_delay = cfg.max_range / cfg.c;   % [s] (direct-path reference delay)

cfg.Nlfm = floor(cfg.Fs * cfg.Tp);           % samples in one LFM pulse
cfg.L    = floor(cfg.max_delay * cfg.Fs);    % samples in delay spread window
cfg.N    = cfg.Nlfm + cfg.L;                 % total sample length per ping
cfg.BW   = abs(cfg.f1 - cfg.f0);             % [Hz] bandwidth

% time axes
cfg.t_ping      = (0:cfg.Np-1).' * cfg.PRI;                     % [s] ping times
cfg.tau_axis    = (0:cfg.L-1).' / cfg.Fs + cfg.initial_delay;   % [s] delay axis
cfg.res_tau     = 1 / cfg.BW;                                   % [s] delay resolution
cfg.tap_delays  = cfg.initial_delay + (0:cfg.L-1).' / cfg.Fs;   % [s] tap delays

end
