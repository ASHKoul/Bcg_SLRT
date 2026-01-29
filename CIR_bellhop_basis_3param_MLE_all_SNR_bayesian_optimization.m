clc; clearvars;

%% ============================ CONFIG =======================================
cfg = struct();
cfg.Np        = 40;                   % # pings
cfg.f0        = 1000;                  % LFM start [Hz]
cfg.f1        = 5000;                  % LFM end [Hz]
cfg.Tp        = 25e-3;                 % pulse duration [s]
cfg.c         = 1500;                  % sound speed [m/s]
cfg.max_range = 2000;                  % [m]
cfg.max_delay = 90e-3;                 % [s]
cfg.PRI       = 0.12;                  % [s]
cfg.Fphi      = 0.999;                 % AR(1) factor
cfg.baseband  = true;
cfg.add_noise = true;

SNR_list = 0:5:30;                     % SNR sweep [dB]

%% ========================== DATA IMPORT ====================================
data   = load("");
cfgTxt = fileread('');

cfg.fc            = str2double(regexp(cfgTxt,'fCarrier:\s*([\d.]+)','tokens','once'));
cfg.waymarkPeriod = str2double(regexp(cfgTxt,'T_WAYMARK:\s*([\d.]+)','tokens','once'));
cfg.Timp_min      = str2double(regexp(cfgTxt,'T_IMPULSE_MIN:\s*([\d.]+)','tokens','once'));
tmp               = regexp(cfgTxt,'T_IMPULSE_NEG_START:\s*([\d.]+)','tokens','once');
cfg.Tneg_cfg      = str2double(tmp); if isempty(cfg.Tneg_cfg), cfg.Tneg_cfg = NaN; end

if isfield(data,'fieldMinDelay') && ~isempty(data.fieldMinDelay)
    cfg.tMin = data.fieldMinDelay;
else
    cfg.tMin = cfg.max_range/cfg.c;
end

[s, cfg.Fs] = audioread("");

cfg.Nlfm          = floor(cfg.Fs * cfg.Tp);
cfg.L             = floor(cfg.max_delay * cfg.Fs);
cfg.N             = cfg.Nlfm + cfg.L;
cfg.BW            = abs(cfg.f1 - cfg.f0);
cfg.initial_delay = cfg.tMin;
cfg.tap_delays    = cfg.initial_delay + (0:cfg.L-1)'/cfg.Fs;
cfg.t_ping        = (0:cfg.Np-1).' * cfg.PRI;
cfg.tau_axis      = (0:cfg.L-1)'/cfg.Fs + cfg.initial_delay;
cfg.res_tau       = 1/cfg.BW;

%% ===================== BASIS BUILDER CONFIG =================================
cfg.overlap_frac = 0.2;     % typical value ~0.2


%% ========================== MAIN LOOP ======================================
for ii = 1:numel(SNR_list)

    fprintf('\n================= SNR = %d dB =================\n', SNR_list(ii));
    cfg.SNRdB = SNR_list(ii);

    % --------- Preprocess data (sets cfg_new.sigma_e, etc.) ----------
    [Y, Atrue, S, U, cfg_new] = data_preprocessing(data, s, cfg); 

    % --------- Build basis ----------
    tau = (0:cfg_new.L-1)'/cfg_new.Fs;
    [B, K] = rbf_basis(tau, cfg_new.res_tau, cfg_new.overlap_frac, cfg_new.Fs); 
    % --------- Fix sigma_e from preprocessing; optimize sigma_q, sigma_c, sigma_d ----------
    kf = cfg_new;
    kf.sigma_e = cfg_new.sigma_e;     % FIXED 

    % Initial guesses (you can tune)
    sigma_q0 = 1e-4;
    %sigma_c0 = 1e-5;
       % sigma_d0 = 1e-5;

    % Bounds (log-domain)
    % LB = log([1e-20; 1e-20; 1e-20]);  % [sigma_q; sigma_c; sigma_d]
    % UB = log([1e+1;  1e+1;  1e+1]);

    %  LB = log([1e-20; 1e-20]);  % [sigma_q; sigma_c; sigma_d]
    % UB = log([1e+1;  1e+1]);

    LB = log(1e-20);  % [sigma_q; sigma_c; sigma_d]
    UB = log(1e+1);



    % --------- bayesopt variables ----------
    vars = [
        optimizableVariable('log_sigma_q',[LB(1), UB(1)], 'Type','real')
        %optimizableVariable('log_sigma_d',[LB(2), UB(2)], 'Type','real')
        %optimizableVariable('log_sigma_d',[LB(3), UB(3)], 'Type','real')
    ];

    % Seed first evaluation at x0 (recommended)
    % x0 = log([sigma_q0; sigma_c0; sigma_d0]);
    % initX = table(x0(1), x0(2), x0(3), ...
    %     'VariableNames', {'log_sigma_q','log_sigma_c','log_sigma_d'});


    %   x0 = log([sigma_q0; sigma_d0]);
    % initX = table(x0(1), x0(2), ...
    %     'VariableNames', {'log_sigma_q','log_sigma_d'});

     x0 = log(sigma_q0);
     initX = table(x0(1), ...
        'VariableNames', {'log_sigma_q' });


    % Objective handle
    objFcn = @(T) nll_bayes(T, Y, S, U, B, kf);

    % Run optimization
 results_bayes = bayesopt(objFcn, vars, ...
    'InitialX', initX, ...
    'IsObjectiveDeterministic', true, ...
    'AcquisitionFunctionName', 'expected-improvement-plus', ...
    'MaxObjectiveEvaluations', 300     , ...
    'UseParallel', true, ...
    'PlotFcn', [], ...          % <-- disables all plots
    'Verbose', 1);              % 0=quiet, 1=summary, 2=more text

    % Extract best point
    bestT = bestPoint(results_bayes);
    % xhat_log = [bestT.log_sigma_q; bestT.log_sigma_c; bestT.log_sigma_d];
     % xhat_log = [bestT.log_sigma_q; bestT.log_sigma_d];
    xhat_log = bestT.log_sigma_q;


    xhat     = exp(xhat_log);   % [sigma_q; sigma_c; sigma_d]

    % --------- Save traces (all evaluated points and their NLLs) ----------
    Xtrace = results_bayes.XTrace;             % table with columns log_sigma_q/c/d
    Ftrace = results_bayes.ObjectiveTrace;     % vector of NLL values (same order)

    % Pack results
    results = struct();
    results.SNRdB        = cfg_new.SNRdB;
    results.xhat_log     = xhat_log;
    results.xhat         = xhat;              % [sigma_q; sigma_c; sigma_d]
    results.minNLL       = results_bayes.MinObjective;
    results.Xtrace       = Xtrace;            % table of log-params
    results.NLLtrace     = Ftrace;            % objective values
    results.cfg          = cfg_new;

    %Save
    stamp    = char(datetime('now','Format','yyyyMMdd_HHmmss'));
    savefile = sprintf('bayesopt_modelH1_reverb_nojitter_SNR%02ddB_%s.mat', ...
        SNR_list(ii), stamp);

    save(savefile, 'results', '-v7.3');
    fprintf('Saved: %s\n', savefile);

    % Display quick summary
    % fprintf('Best: sigma_q=%.3e, sigma_c=%.3e, sigma_d=%.3e | minNLL=%.6g\n', ...
    %     results.xhat(1), results.xhat(2), results.xhat(3), results.minNLL);

   % fprintf('Best: sigma_q=%.3e, sigma_d=%.3e| minNLL=%.6g\n', ...
   %      results.xhat(1), results.xhat(2), results.minNLL);

   fprintf('Best: sigma_q=%.3e | minNLL=%.6g\n', ...
        results.xhat(1), results.minNLL);


end

%% ===================== OBJECTIVE FUNCTIONS =================================
function NLL = nll_bayes(T, Y, S, U, B, kf)
% T is a table row with fields: log_sigma_q, log_sigma_c, log_sigma_d

kf_use         = kf;
kf_use.sigma_q = exp(T.log_sigma_q);
kf_use.sigma_c = 0;%exp(T.log_sigma_c);
kf_use.sigma_d = 0;%exp(T.log_sigma_d);

out = kalman_filter_bellhop_basis(Y, S, U, B, kf_use);

% Must return a scalar objective
NLL = out.NLL;
end

 
