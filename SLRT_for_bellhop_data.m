clc;
clearvars;

%% ========================================================================
%  FILES / INPUTS
%  ========================================================================

% Data
data   = load("Bellhop_nojitter_reverb005_lfm_train_0025_012_20251111T180648.mat");
N_op   = length(data.outputSignal); %#ok<NASGU>  % kept (not used later)

% Configuration + waveform
cfgTxt = fileread('config_nojitter_reverb005_lfm_train_0025_012_20251111T180648.txt');
[s, Fs] = audioread("lfm_train_0025_012.wav");

%% ========================================================================
%  CONFIG
%  ========================================================================

cfg = settings_data(cfgTxt, Fs);
cfg = generate_tau_beta(cfg);

Ns  = numel(cfg.SNRdB_list);

%% ========================================================================
%  (OPTIONAL) OPTIMIZATION / TRAINING
%  NOTE: This loop computes MLE training results per SNR via trainingdataopt.
%        (Your original behavior is kept: `results` is overwritten each loop.)
%  ========================================================================

for isnr = 1:length(cfg.SNRdB_list)

    cfg.SNRdB      = cfg.SNRdB_list(isnr);
    disp(cfg.Fs);

    cfg.add_target = false;
    cfg.add_noise  = true;

    [Ydata, Atrue, S, U, ~, cfg_new] = data_preprocessing_withtarget(data, s, cfg); %#ok<ASGLU>
    Ktest  = size(Ydata, 2); %#ok<NASGU>

    Ytrain  = Ydata(:, 1:cfg_new.Ntrain);
    results = trainingdataopt(Ytrain, S, U, cfg_new); %#ok<NASGU>
end

%% ========================================================================
%  LOAD PRECOMPUTED MLE PARAMS (xhat) + PAGE THRESHOLDS (h)
%  ========================================================================

tmpH1  = load('');
tmpH3  = load('');
tmpH2c = load('');
tmpH2d = load('');

xhat_M1  = tmpH1.xhat_mle_H1_reverb_nojitter;    % 5 x Ns cell
xhat_M3  = tmpH3.xhat_mle_H3_reverb_nojitter;    % 5 x Ns cell
xhat_M2c = tmpH2c.xhat_mle_H2c_reverb_nojitter;  % 5 x Ns cell
xhat_M2d = tmpH2d.xhat_mle_H2d_reverb_nojitter;  % 5 x Ns cell

thr  = load('');
h_M1  = thr.h_models(1,:).';   % Ns x 1
h_M3  = thr.h_models(2,:).';   % Ns x 1
h_M2c = thr.h_models(3,:).';   % Ns x 1
h_M2d = thr.h_models(4,:).';   % Ns x 1

%% ========================================================================
%  MODEL LIST (names, MLE params, thresholds)
%  ========================================================================

models(1).name = "Noise Only";
models(1).xhat = xhat_M1;
models(1).hvec = h_M1;

models(2).name = "Noise+Reverb+Jitter";
models(2).xhat = xhat_M3;
models(2).hvec = h_M3;

models(3).name = "Noise+Reverb";
models(3).xhat = xhat_M2c;
models(3).hvec = h_M2c;

models(4).name = "Noise+Jitter";
models(4).xhat = xhat_M2d;
models(4).hvec = h_M2d;

%% ========================================================================
%  DETECTION SETTINGS
%  ========================================================================

mode      = 2;    % 1 = stop at first detection
Kwin      = 20;   % Pd within Kwin pings after tau0
isnr_show = 1;    % store one “showcase” run for plots

%% ========================================================================
%  STORAGE ALLOCATION
%  ========================================================================

Res = struct();
for m = 1:numel(models)

    Res(m).name = models(m).name;

    Res(m).tau0 = nan(1, Ns);

    % post-change delay stats
    Res(m).delay_mean = nan(1, Ns);
    Res(m).delay_med  = nan(1, Ns);

    % post-change detection-time stats (tau_det = tau0 + delay_censored)
    Res(m).tau_mean = nan(1, Ns);
    Res(m).tau_med  = nan(1, Ns);

    % per-SNR detailed output
    Res(m).out = cell(1, Ns);

    % detection probabilities
    Res(m).Pd_end = nan(1, Ns);
    Res(m).Pd_win = nan(1, Ns);

    % one showcase run
    Res(m).show = [];
end

%% ========================================================================
%  MAIN SWEEP OVER SNR
%  ========================================================================

for isnr = 1:Ns

    cfgS        = cfg;
    cfgS.SNRdB  = cfg.SNRdB_list(isnr);

    cfgS.add_target = true;
    cfgS.add_noise  = false;        % clean input; detector adds noise internally
    cfgS            = generate_tau_beta(cfgS);

    [Yclean_target, ~, S, U, ~, cfg_new] = data_preprocessing_withtarget(data, s, cfgS);
    Ktest = size(Yclean_target, 2);

    tau0 = cfg_new.Ntrain + 1;

    % --- Build a template that starts at ping 1 (proto)
    cfgProto            = cfgS;
    cfgProto.add_target = true;
    cfgProto.add_noise  = false;
    cfgProto.Ntrain     = 0;
    cfgProto            = generate_tau_beta(cfgProto);

    [~, ~, ~, ~, S_tgt_proto, ~] = data_preprocessing_withtarget(data, s, cfgProto);

    % --- Align lengths
    Kproto = size(S_tgt_proto, 2);
    Kuse   = min(Ktest, Kproto);

    Yuse = Yclean_target(:, 1:Kuse);
    Stpl = S_tgt_proto(:, 1:Kuse);

    %% ---------------------- MODEL LOOP ----------------------
    for m = 1:numel(models)

        xhat = models(m).xhat(:, isnr);
        h    = models(m).hvec(isnr);

        % --- KF params (from cfg_new + MLE xhat)
        kf         = cfg_new;
        kf.Nmc     = cfg.Nmc;

        kf.sigma_q = xhat{1};
        kf.sigma_c = xhat{2};
        kf.sigma_d = xhat{3};

        width_est  = xhat{4};

        kf.sigma_e = xhat{5};

        % --- Build basis B
        tau = (0:kf.L-1)' / kf.Fs;
        [B, ~] = rbf_basis(tau, width_est, kf.overlap_frac, kf.Fs);

        %% ------------------ MC Page detection ------------------
        N   = size(Yuse, 1);
        Np  = size(Yuse, 2);
        Nmc = kf.Nmc;

        ell_all = zeros(Nmc, Np);
        T_all   = zeros(Nmc, Np);
        det_all = zeros(Nmc, Np, 'uint8');

        tau_det = Np * ones(Nmc, 1);
        missed  = true(Nmc, 1);

        % ensure sigma_e is usable in parfor
        sigma_e = kf.sigma_e;
        if ~isscalar(sigma_e)
            sigma_e = sigma_e(:).';   % force row
            if numel(sigma_e) < Np
                error('kf.sigma_e length %d < Np=%d', numel(sigma_e), Np);
            end
            sigma_e = sigma_e(1:Np);
        end

        parfor nmc = 1:Nmc

            Z = (randn(N, Np) + 1j*randn(N, Np)) / sqrt(2);

            if isscalar(sigma_e)
                E = sigma_e * Z;
            else
                E = Z .* sigma_e;     % (N×Np).*(1×Np)
            end

            Yrun = Yuse + E;

            % single-run function: returns .ell, .T, .det, .tau_det, .missed
            run = kalman_filter_pagetest(Yrun, S, U, B, kf, Stpl, h, mode);

            ell_all(nmc,:) = run.ell;
            T_all(nmc,:)   = run.T;
            det_all(nmc,:) = run.det;
            tau_det(nmc)   = run.tau_det;
            missed(nmc)    = run.missed;
        end

        %% ------------------ Package outputs ------------------
        out = struct();

        out.mc.ell_all  = ell_all;
        out.mc.T_all    = T_all;
        out.mc.det_all  = det_all;
        out.mc.tau_det  = tau_det;
        out.mc.missed   = missed;

        out.mc.ell_mean = mean(ell_all, 1);
        out.mc.T_mean   = mean(T_all,   1);
        out.mc.p_detect = mean(det_all == uint8(2), 1);

        out.mc.tau_det_mean = mean(tau_det);
        out.mc.tau_det_med  = median(tau_det);

        % page trace for plotting (use run #1)
        out.page = struct();
        out.page.T_hist = T_all(1,:);
        out.page.det    = det_all(1,:);
        out.page.ell_k  = ell_all(1,:);
        out.page.h      = h;

        %% ------------------ Post-change detection logic ------------------
        det_all_global = out.mc.det_all;     % Nmc x Kuse
        tau_det_global = out.mc.tau_det;     % Nmc x 1
        missed_global  = out.mc.missed;      % Nmc x 1
        Nmc_now        = numel(tau_det_global);

        tau_det_post_change = nan(Nmc_now, 1);

        % If first-ever detection is after tau0, it's post-change
        idx_direct = (~missed_global) & (tau_det_global >= tau0);
        tau_det_post_change(idx_direct) = tau_det_global(idx_direct);

        % If first-ever detection is before tau0, search again after tau0
        idx_need = (~missed_global) & (tau_det_global < tau0);
        rows = find(idx_need).';
        for rr = rows
            idxPost = find(det_all_global(rr, tau0:Kuse) == uint8(2), 1, 'first');
            if ~isempty(idxPost)
                tau_det_post_change(rr) = (tau0 - 1) + idxPost;
            end
        end

        miss_post_change = isnan(tau_det_post_change);

        delay_cens        = max(0, Kuse - tau0 + 1);
        delay_post_change = delay_cens * ones(Nmc_now, 1);

        ok = ~miss_post_change;
        delay_post_change(ok) = tau_det_post_change(ok) - tau0;

        Pd_end = mean(~miss_post_change);
        Pd_win = mean(~miss_post_change & (tau_det_post_change <= (tau0 + Kwin - 1)));

        tau_det_cens = tau0 + delay_post_change;

        %% ------------------ Store summary stats ------------------
        Res(m).tau0(isnr) = tau0;

        Res(m).Pd_end(isnr) = Pd_end;
        Res(m).Pd_win(isnr) = Pd_win;

        Res(m).delay_mean(isnr) = mean(delay_post_change);
        Res(m).delay_med(isnr)  = median(delay_post_change);

        Res(m).tau_mean(isnr) = mean(tau_det_cens);
        Res(m).tau_med(isnr)  = median(tau_det_cens);

        Res(m).out{isnr} = out;

        % store one showcase run for Page trace plots
        if isnr == isnr_show
            Res(m).show.SNRdB = cfgS.SNRdB;
            Res(m).show.h     = h;
            Res(m).show.tau0  = tau0;
            Res(m).show.Kuse  = Kuse;
            Res(m).show.out   = out;
        end

        fprintf('[%s] SNR=%g dB: tau0=%d, Pd_end=%.3f, Pd_win=%.3f, meanDelay=%.2f\n', ...
            Res(m).name, cfgS.SNRdB, tau0, Pd_end, Pd_win, Res(m).delay_mean(isnr));
    end
end

%% ========================================================================
%  PLOTS
%  ========================================================================

SNRv = cfg.SNRdB_list(:).';   % SNR grid

%% (A) Page statistic trace at one SNR (run #1), one figure per model
for m = 1:numel(Res)

    if ~isfield(Res(m), 'show') || isempty(Res(m).show)
        continue;
    end

    sh   = Res(m).show;
    out  = sh.out;
    h    = sh.h;
    tau0 = sh.tau0;

    if ~isfield(out, 'page') || ~isfield(out.page, 'T_hist') || isempty(out.page.T_hist)
        continue;
    end

    T   = out.page.T_hist;
    det = out.page.det;

    figure;
    hold on; grid on; box on;

    plot(T, 'LineWidth', 1.6);
    yline(h, 'r--', 'h');
    xline(tau0, 'k:', '\tau_0');

    kDet = find(det == uint8(2));
    if ~isempty(kDet)
        plot(kDet, T(kDet), 'o', 'MarkerSize', 5);
    end

    xlabel('Ping index k');
    ylabel('T_k');

    % Optional title:
    % title(sprintf('%s: Page statistic (run #1), SNR = %g dB', Res(m).name, sh.SNRdB));
end

%% (B) Mean detection delay vs SNR: E[tau_det - tau0]
figure;
hold on; grid on; box on;

for m = 1:numel(Res)
    delay_mean = Res(m).delay_mean;
    plot(SNRv, delay_mean, 'LineWidth', 2);
end

xlabel('SNR [dB]');
ylabel('Mean delay (pings) to detection');
legend('$\mathcal{M}_0$', '$\mathcal{M}_{cd}$', '$\mathcal{M}_c$', '$\mathcal{M}_d$', ...
    'Interpreter','latex');

%% (D1) P_d vs SNR (by end of horizon)
figure;
hold on; grid on; box on;

for m = 1:numel(Res)
    plot(SNRv, Res(m).Pd_end, 'LineWidth', 2);
end

xlabel('SNR [dB]');
ylabel('P_d');
legend('$\mathcal{M}_0$', '$\mathcal{M}_{cd}$', '$\mathcal{M}_c$', '$\mathcal{M}_d$', ...
    'Interpreter','latex');

%% (D2) P_d vs SNR (within Kwin pings after tau0)
figure;
hold on; grid on; box on;

for m = 1:numel(Res)
    plot(SNRv, Res(m).Pd_win, 'LineWidth', 2);
end

xlabel('SNR [dB]');
ylabel('P_d (windowed)');
legend('$\mathcal{M}_0$', '$\mathcal{M}_{cd}$', '$\mathcal{M}_c$', '$\mathcal{M}_d$', ...
    'Interpreter','latex');
