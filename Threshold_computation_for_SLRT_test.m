clearvars; clc;

%% ----------------------- data load -----------------------
data  = load("Bellhop_reverb_nojitter005_lfm_train_0025_012_20251111T225300.mat");
cfgTxt = fileread('config_reverb_nojitter005_lfm_train_0025_012_20251111T225300.txt');
[s, Fs] = audioread("lfm_train_0025_012.wav");

%% ------------------ Load configurations ------------------
cfg = settings_data(cfgTxt,Fs);

Ns = numel(cfg.SNRdB_list);

%% -------- Load MLE cells (5 x Ns): {sigma_q,sigma_c,sigma_d,width,sigma_e} --------
tmpH1 = load('xhat_mle_H1_reverb_nojitter.mat');
tmpH3 = load('xhat_mle_H3_reverb_nojitter.mat');
tmpH2c = load('xhat_mle_H2c_reverb_nojitter.mat');
tmpH2d = load('xhat_mle_H2d_reverb_nojitter.mat');


xhat_M1 = tmpH1.xhat_mle_H1_reverb_nojitter;   % 5 x Ns cell
xhat_M3 = tmpH3.xhat_mle_H3_reverb_nojitter;   % 5 x Ns cell
xhat_M2c = tmpH2c.xhat_mle_H2c_reverb_nojitter;   % 5 x Ns cell
xhat_M2d = tmpH2d.xhat_mle_H2d_reverb_nojitter;   % 5 x Ns cell





models(1).name = "Noise Only";
models(1).xhat = xhat_M1;


models(2).name = "Noise+Reverb+Jitter";
models(2).xhat = xhat_M3;


models(3).name = "Noise+Reverb";
models(3).xhat = xhat_M2c;


models(4).name = "Noise+Jitter";
models(4).xhat = xhat_M2d;

%% Detector behavior to match online code
mode     = 2;         % continue (do not stop after first detection)
hHuge    = realmax;   % never triggers "detected" branch; T_k evolution unchanged in your code
Pfa      = cfg.Pfa;
pQuant   = 100*(1 - Pfa);

%% ---------------- Allocate outputs ----------------
h_models   = nan(numel(models), Ns);
Gmax_all   = cell(numel(models), Ns);

%% ===================== LOOP OVER assumed SNR/CNR =====================
for isnr = 1:Ns

    % ----- assumed CNR/SNR for TEMPLATE (and for your MLE params) -----
    cfgS = cfg;
    cfgS.SNRdB = cfg.SNRdB_list(isnr);

    % (1) Build H0 clean baseline (NO target, NO added noise)
    cfgH0 = cfgS;
    cfgH0.add_target = false;
    cfgH0.add_noise  = false;
    cfgH0.Ntrain     = 0;
    cfgH0 = generate_tau_beta(cfgH0);
    
    [Y0_clean_full, ~, S, U, ~, cfgH0_out] = data_preprocessing_withtarget(data, s, cfgH0);

    % true noise used to *simulate* H0 realizations
    sigma_true = cfgH0_out.sigma_e;

    % (2) Build template that starts at ping 1 (Ntrain=0)
    cfgTpl = cfgS;
    cfgTpl.add_target = true;
    cfgTpl.add_noise  = false;
    cfgTpl.Ntrain     = 0;
    cfgTpl = generate_tau_beta(cfgTpl);


    % figure;
    % subplot(121)
    % plot(cfgTpl.tau_target)
    % subplot(122)
    % plot(cfgTpl.beta_target)
    [~, ~, ~, ~, S_tgt_full, cfgTpl_out] = data_preprocessing_withtarget(data, s, cfgTpl);

    % (3) Align horizons
    Kuse = min([size(Y0_clean_full,2), size(S_tgt_full,2), cfgH0_out.Np, cfgTpl_out.Np]);

    Y0_clean = Y0_clean_full(:,1:Kuse);
    Stpl     = S_tgt_full(:,1:Kuse);

    % prepare true sigma vector (for noise simulation)
    if isscalar(sigma_true)
        sigma_true_use = sigma_true;
    else
        sigma_true_use = sigma_true(:).';
        if numel(sigma_true_use) < Kuse
            error('True sigma_e length %d < Kuse=%d', numel(sigma_true_use), Kuse);
        end
        sigma_true_use = sigma_true_use(1:Kuse);
    end

    N   = size(Y0_clean,1);
    Np  = size(Y0_clean,2);
    Nmc = cfg.Nmc;

    fprintf('SNR=%g dB: calibrating thresholds with Kuse=%d, Nmc=%d, Pfa=%g\n', ...
        cfgS.SNRdB, Kuse, Nmc, Pfa);

    % ===================== LOOP OVER MODELS =====================
    for m = 1:numel(models)

        xhat = models(m).xhat(:,isnr);

        % detector (assumed) params
        kf = cfgH0_out;
        kf.Nmc     = Nmc;
        kf.sigma_q = xhat{1};
        kf.sigma_c = xhat{2};
        kf.sigma_d = xhat{3};
        width_est  = xhat{4};
        kf.sigma_e = xhat{5};     % assumed sigma_e (used inside KF likelihood)

        % basis
        tau = (0:kf.L-1)'/kf.Fs;
        [B, ~] = rbf_basis(tau, width_est, kf.overlap_frac, kf.Fs);

        Gmax = zeros(Nmc,1);

        parfor nmc = 1:Nmc
            Z = (randn(N,Np) + 1j*randn(N,Np)) / sqrt(2);

            % simulate H0 noise realization
            if isscalar(sigma_true_use)
                E = sigma_true_use * Z;
            else
                E = Z .* sigma_true_use;   % (N×Np).*(1×Np)
            end
            Y0 = Y0_clean + E;

            % run the SAME online detector logic (resets included)
            run = kalman_filter_pagetest(Y0, S, U, B, kf, Stpl, hHuge, mode);

            % maximum Page statistic under H0
            Tk = run.T;
            if ~isvector(Tk) || numel(Tk) < Np
                error('kalman_filter_pagetest: run.T has wrong size.');
            end
            Gmax(nmc) = max(Tk(2:end));   
        end

        % threshold at (1-Pfa) quantile
        h_models(m,isnr) = prctile(Gmax, pQuant);
        Gmax_all{m,isnr} = Gmax;

        fprintf('  [%s] h=%.6g\n', models(m).name, h_models(m,isnr));
    end
end

%% ===================== SAVE =====================
fname = sprintf('page_thresholds_fixed_rnj_CNR%d_pings%d_NMC%d_Pfa%g.mat', ...
     cfg.CNRdB, Kuse,cfg.Nmc, cfg.Pfa);

 save(fname, ...
     'h_models','Gmax_all','models','cfg','Pfa');

fprintf('\nSaved thresholds to: %s\n', fname);
