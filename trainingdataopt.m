function results =trainingdataopt(Y,S,U,cfg)



opts = optimoptions('fminunc', ...
    'Algorithm','quasi-newton', ...
    'Display','iter', ...
    'UseParallel',true, ...
    'MaxFunctionEvaluations',200, ...
    'FiniteDifferenceType','central', ...
    'StepTolerance',1e-10,'OptimalityTolerance',1e-10);

disp(cfg.Fs)

sigma_q0= 1e-4;
sigma_d0 = 1e-5;
sigma_c0 = 1e-5;

%  x0=log(sigma_q0);
 x0 = log([sigma_q0; sigma_c0; sigma_d0]);
 %x0 = log([sigma_q0; sigma_d0]);
% LB = log([1e-8; 1e-8; 1e-8; cfg.res_tau]);
% UB = log([1e-0; 1e-0; 1e-0; 5*cfg.res_tau]);

objfun = @(x) nll_local(x, Y, S, U, cfg);


% [xhat, ~, ~, ~] = fmincon(objfun, x0,[], [], [], [], LB, UB, [], opts);
[xhat, ~, ~, ~] = fminunc(objfun, x0, opts);

results.xhat      = exp(xhat);
results.cfg       = cfg;

stamp    = char(datetime('now','Format','yyyyMMdd_HHmmss'));
    savefile = sprintf('H3_target_bellhop_noreverb_jitter_noise_SNR%ddB_%s.mat', cfg.SNRdB, stamp);
    save(savefile, 'results', '-v7.3');
    fprintf('Saved results to %s\n', savefile);



    function NLL = nll_local(x, Y, S, U, cfg)
        sigma_q = exp(x(1));
        sigma_c = exp(x(2));
        sigma_d = exp(x(3));
        width_bf= cfg.res_tau;
        tau = (0:cfg.L-1)'/cfg.Fs;
        [B, ~] = rbf_basis(tau, width_bf, cfg.overlap_frac, cfg.Fs);  % unit-norm columns

        kf          = cfg;
        kf.sigma_q  = sigma_q;
        kf.sigma_c  = sigma_c;
        kf.sigma_d  = sigma_d;
        kf.sigma_e  = cfg.sigma_e;


        out = kalman_filter_bellhop_basis(Y, S, U, B, kf);

        NLL = out.NLL;
    end
end


   


