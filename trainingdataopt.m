function results = trainingdataopt(Y, S, U, cfg)

%% ========================================================================
%  OPTIMIZER SETTINGS
%  ========================================================================

opts = optimoptions('fminunc', ...
    'Algorithm','quasi-newton', ...
    'Display','iter', ...
    'UseParallel', true, ...
    'MaxFunctionEvaluations', 200, ...
    'FiniteDifferenceType','central', ...
    'StepTolerance', 1e-10, ...
    'OptimalityTolerance', 1e-10);

disp(cfg.Fs);

%% ========================================================================
%  INITIALIZATION (log-domain)
%  ========================================================================

sigma_q0 = 1e-4;
sigma_c0 = 1e-5;
sigma_d0 = 1e-5;

x0 = log([sigma_q0; sigma_c0; sigma_d0]);

% objective handle
objfun = @(x) nll_local(x, Y, S, U, cfg);

%% ========================================================================
%  OPTIMIZATION
%  ========================================================================

[xhat, ~, ~, ~] = fminunc(objfun, x0, opts);

%% ========================================================================
%  PACKAGE + SAVE
%  ========================================================================

results.xhat = exp(xhat);
results.cfg  = cfg;

stamp    = char(datetime('now','Format','yyyyMMdd_HHmmss'));
savefile = sprintf('H3_target_bellhop_reverb_nojitter_noise_SNR%ddB_%s.mat', cfg.SNRdB, stamp);

save(savefile, 'results', '-v7.3');
fprintf('Saved results to %s\n', savefile);

%% ========================================================================
%  NESTED: NLL OBJECTIVE
%  ========================================================================

function NLL = nll_local(x, Y, S, U, cfg)

    % parameters (positive)
    sigma_q = exp(x(1));
    sigma_c = exp(x(2));
    sigma_d = exp(x(3));

    % basis width is fixed here (not optimized)
    width_bf = cfg.res_tau;

    % build RBF basis
    tau = (0:cfg.L-1)' / cfg.Fs;
    [B, ~] = rbf_basis(tau, width_bf, cfg.overlap_frac, cfg.Fs);

    % KF config
    kf         = cfg;
    kf.sigma_q = sigma_q;
    kf.sigma_c = sigma_c;
    kf.sigma_d = sigma_d;
    kf.sigma_e = cfg.sigma_e;

    % run filter and return scalar NLL
    out = kalman_filter_bellhop_basis(Y, S, U, B, kf);
    NLL = out.NLL;
end

end
