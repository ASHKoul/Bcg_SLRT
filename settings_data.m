function cfg = settings_data(cfgTxt, Fs)

cfg = struct();
cfg.Np        = 100;
cfg.Ntrain    = 40;

cfg.f0        = 1000;
cfg.f1        = 5000;
cfg.Tp        = 25e-3;
cfg.c         = 1500;
cfg.max_range = 2000;
cfg.max_delay = 90e-3;
cfg.PRI       = 0.12;
cfg.Fphi      = 1-eps;
cfg.Pfa = 0.05;      % P_FA
cfg.Nmc = 10;         % increase for stable results

cfg.Fs        = Fs;
cfg.baseband  = true;
cfg.add_noise = true;
cfg.add_target= true;
cfg.target_status="fixed";
cfg.CNRdB      = 30;
cfg.SNRdB_list = 0:5:20;
cfg.overlap_frac = 0.2;

cfg.Tx_pos = [0; 0];       % [2x1]
cfg.Rx_pos = [2000; 0];    % [2x1]


% parse text
cfg.fc            = str2double(regexp(cfgTxt,'fCarrier:\s*([\d.]+)','tokens','once'));
cfg.waymarkPeriod = str2double(regexp(cfgTxt,'T_WAYMARK:\s*([\d.]+)','tokens','once'));
cfg.Timp_min      = str2double(regexp(cfgTxt,'T_IMPULSE_MIN:\s*([\d.]+)','tokens','once'));
tmp               = regexp(cfgTxt,'T_IMPULSE_NEG_START:\s*([\d.]+)','tokens','once');
cfg.Tneg_cfg      = str2double(tmp); if isempty(cfg.Tneg_cfg), cfg.Tneg_cfg = NaN; end

cfg.initial_delay         = cfg.max_range/cfg.c;

% derived
cfg.Nlfm    = floor(cfg.Fs * cfg.Tp);
cfg.L       = floor(cfg.max_delay * cfg.Fs);
cfg.N       = cfg.Nlfm + cfg.L;
cfg.BW      = abs(cfg.f1 - cfg.f0);


cfg.t_ping        = (0:cfg.Np-1).' * cfg.PRI;
cfg.tau_axis      = (0:cfg.L-1)'/cfg.Fs + cfg.initial_delay;
cfg.res_tau       = 1/cfg.BW;
cfg.tap_delays      = cfg.initial_delay + (0:cfg.L-1).'/cfg.Fs;



end
