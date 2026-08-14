function data = tbd_bis_active(Ymat1, Ymat2, model, settings)

% ============================================================
% Setup
% ============================================================
Krun = size(Ymat1, 2) - 1;

data.xh = NaN(size(model.F, 1), Krun);
data.q = NaN(1, Krun);

use_bg = settings.cfg.add_background && settings.bg_tracking;

tping = settings.cfg.t_ping(2:settings.cfg.Np);
do_plot = settings.do_plot;
t_old = -inf;

settings.Qbgmat_1 = diag(settings.Qbg_1);
settings.Qbgmat_2 = diag(settings.Qbg_2);

% ============================================================
% EKF init from first ping
% ============================================================
post_ekf = background_intialization(Ymat1(:, 1), Ymat2(:, 1), use_bg, settings);

% ============================================================
% PF init
% State:
% x = [x; y; vx; vy; eta_dB]
% ============================================================
x = generate_new_particles(Ymat1(:, 1), Ymat2(:, 1), settings, post_ekf, settings.Ns);
w = ones(1, settings.Ns) / settings.Ns;
q = settings.q0;

% ============================================================
% Main PF loop
% ============================================================
for kout = 1:Krun

    k = kout + 1;

    y1k = Ymat1(:, k);
    y2k = Ymat2(:, k);

    % --------------------------------------------------------
    % Background prior
    % --------------------------------------------------------
    prior_ekf = Time_update_ekf(post_ekf, settings, use_bg, q, settings.q_threshold, k, settings.cfg.N_target_appearance);

    % --------------------------------------------------------
    % Bernoulli prediction
    % --------------------------------------------------------
    qp = settings.pb * (1 - q) + settings.ps * q;

    xps = model.F * x + model.G * (settings.modelLQ * randn(size(model.Q, 1), settings.Ns));

    % Clamp surviving particles
    xps(1, :) = min(max(xps(1, :), settings.x_region(1)), settings.x_region(2));
    xps(2, :) = min(max(xps(2, :), settings.y_region(1)), settings.y_region(2));
    xps(3, :) = min(max(xps(3, :), settings.vx_region(1)), settings.vx_region(2));
    xps(4, :) = min(max(xps(4, :), settings.vy_region(1)), settings.vy_region(2));
    xps(5, :) = min(max(xps(5, :), settings.tgt_eta_region(1)), settings.tgt_eta_region(2));

    % Birth particles
    xpb = generate_new_particles(y1k, y2k, settings, prior_ekf, settings.Nb);

    xp = [xps, xpb];

    wp_surv = (settings.ps * q / qp) * w;
    wp_birth = (settings.pb * (1 - q) / qp) * (1/settings.Nb) * ones(1, settings.Nb);
    wp = [wp_surv, wp_birth];

    % --------------------------------------------------------
    % Likelihoods
    % --------------------------------------------------------
    [Lt, Lnt] = get_likelihoods(y1k, y2k, xp, settings, prior_ekf);

    ratio = Lt ./ Lnt;

    I = ratio * wp.';
    q = (I * qp) / ((1 - qp) + qp * I);

    % Particle posterior under H1
    wp = Lt .* wp;
    sw = sum(wp);
    wp = wp / sw;

    % Save estimate
    data.xh(:, kout) = xp * wp.';
    data.q(kout) = q;

    % Resample
    ind = sysresample(wp, settings.Ns);
    x = xp(:, ind);
    w = ones(1, settings.Ns) / settings.Ns;

    % Debug particle diagnostics (kept unchanged; currently always enabled)
    figure(105)
    subplot(311)
    scatter(x(1, :), x(2, :));
    xlabel('x')
    ylabel('y')
    title('posterior p(x/y) after resampling')
    subplot(312)
    scatter(xp(1, :), xp(2, :), 10+1e5*wp);
    xlim([settings.x_region(1), settings.x_region(2)])
    ylim([settings.y_region(1), settings.y_region(2)])
    xlabel('x')
    ylabel('y')
    title('Posterior p(x_k/Y_k) with  their corresponding weights')

    % Plot
    if do_plot && (tping(kout) - t_old) > 1
        t_old = tping(kout);
        plot_results(settings, data, x, kout);
    end

    % --------------------------------------------------------
    % Background measurement update
    % --------------------------------------------------------
    post_ekf = meas_update_ekf(y1k, y2k, prior_ekf, use_bg, q, settings.q_threshold, k, settings.cfg.N_target_appearance);

end

end

% ============================================================
% Birth particles
% ============================================================
function x = generate_new_particles(y1, y2, settings, in_ekf, Np)

dx = diff(settings.x_region);
dy = diff(settings.y_region);
deta = diff(settings.tgt_eta_region_new_targets);

% ------------------------------------------------------------
% Position and eta birth
% ------------------------------------------------------------
px = settings.x_region(1) + dx   * rand(1, Np);
py = settings.y_region(1) + dy   * rand(1, Np);
eta = settings.tgt_eta_region_new_targets(1) + deta * rand(1, Np);

% ------------------------------------------------------------
% Position-dependent velocity birth
% Surveillance center = Tx/Rx centroid
% ------------------------------------------------------------
pc = (settings.cfg.Tx_pos(:) + settings.cfg.Rx_pos{1}(:) + settings.cfg.Rx_pos{2}(:))/3;

Vx = max(abs(settings.vx_region));
Vy = max(abs(settings.vy_region));

vx_abs = Vx * rand(1, Np);
vy_abs = Vy * rand(1, Np);

% left of surveillance center  -> vx positive
% right of surveillance center -> vx negative
sx = 2*(px <= pc(1)) - 1;

% below surveillance center -> vy positive
% above surveillance center -> vy negative
sy = 2*(py <= pc(2)) - 1;

vx = sx .* vx_abs;
vy = sy .* vy_abs;

x = [px;
    py;
    vx;
    vy;
    eta];

if isfield(settings, 'birth_preselect') && settings.birth_preselect

    [Lt, ~] = get_likelihoods(y1, y2, x, settings, in_ekf);
    sLt = sum(Lt);

    if isfinite(sLt) && sLt > 0
        Lt = Lt / sLt;
        % Birth-particle diagnostic plot
        figure(105)
        subplot(313)
        scatter(x(1, :), x(2, :), 10+1e5*Lt);
        xlim([settings.x_region(1), settings.x_region(2)])
        ylim([settings.y_region(1), settings.y_region(2)])
        xlabel('x')
        ylabel('y')
        title('Birth Particles with the Likeihoods as the weights')

        ind = sysresample(Lt, Np);
        x = x(:, ind);
    end
end

end

function [Lt, Lnt] = get_likelihoods(y1, y2, xp, settings, in_ekf)

% ------------------------------------------------------------
% Background-only residuals
% ------------------------------------------------------------
H = in_ekf.H;

if isfield(in_ekf, 'theta_hat_1') && ~isempty(in_ekf.theta_hat_1)
    rbg1 = y1 - H * in_ekf.theta_hat_1;
    rbg2 = y2 - H * in_ekf.theta_hat_2;
else
    rbg1 = y1;
    rbg2 = y2;
end

logLnt = likelihood_gen(rbg1, rbg2, in_ekf);
logLnt = real(logLnt(:));

% ------------------------------------------------------------
% Particle states
% ------------------------------------------------------------
Nparticles = size(xp, 2);
logLt_all = -inf(1, Nparticles);

px = xp(1, :);
py = xp(2, :);
vx = xp(3, :);
vy = xp(4, :);
eta = xp(5, :); % eta_dB

% ------------------------------------------------------------
% State gating
% ------------------------------------------------------------
valid = isfinite(px)  & isfinite(py)  & ...
    isfinite(vx)  & isfinite(vy)  & ...
    isfinite(eta) & ...
    px  >= settings.x_region(1)       & px  <= settings.x_region(2)       & ...
    py  >= settings.y_region(1)       & py  <= settings.y_region(2)       & ...
    vx  >= settings.vx_region(1)      & vx  <= settings.vx_region(2)      & ...
    vy  >= settings.vy_region(1)      & vy  <= settings.vy_region(2)      & ...
    eta >= settings.tgt_eta_region(1) & eta <= settings.tgt_eta_region(2);

idx = find(valid);

if ~isempty(idx)

    px = px(idx);
    py = py(idx);
    vx = vx(idx);
    vy = vy(idx);
    eta = eta(idx);

    t_rel = settings.t_rec(:);

    % --------------------------------------------------------
    % Geometry
    % --------------------------------------------------------
    [tau1_abs, beta1] = cart_to_tau_beta([px; py], [vx; vy], settings.cfg, 1);
    [tau2_abs, beta2] = cart_to_tau_beta([px; py], [vx; vy], settings.cfg, 2);

    tau1_abs = tau1_abs(:).';
    tau2_abs = tau2_abs(:).';
    beta1 = beta1(:).';
    beta2 = beta2(:).';

    tau_direct = settings.cfg.initial_delay;

    tau_rel_1 = tau1_abs - tau_direct;
    tau_rel_2 = tau2_abs - tau_direct;

    valid_geom = isfinite(tau_rel_1) & isfinite(tau_rel_2) & ...
        isfinite(beta1)     & isfinite(beta2)     & ...
        beta1 > 0           & beta2 > 0           & ...
        tau_rel_1 >= 0      & tau_rel_2 >= 0;

    if any(valid_geom)

        idx2 = idx(valid_geom);

        px_g = px(valid_geom);
        py_g = py(valid_geom);
        eta_g = eta(valid_geom);
        tau_rel_1_g = tau_rel_1(valid_geom);
        tau_rel_2_g = tau_rel_2(valid_geom);
        beta1_g = beta1(valid_geom);
        beta2_g = beta2(valid_geom);

        % ----------------------------------------------------
        % Vectorized waveform generation
        % Each column is one particle waveform
        % ----------------------------------------------------
        tp1 = beta1_g .* (t_rel - tau_rel_1_g);
        tp2 = beta2_g .* (t_rel - tau_rel_2_g);

        U1 = sym_chirp_bb_eval(tp1, settings);
        U2 = sym_chirp_bb_eval(tp2, settings);

        E1 = sum(abs(U1).^2, 1);
        E2 = sum(abs(U2).^2, 1);

        % Match the working parfor behavior:
        % both receivers must have valid target support.
        good = isfinite(E1) & isfinite(E2) & E1 > 0 & E2 > 0;
        if any(good)

            idx_good = idx2(good);

            U1g = U1(:, good);
            U2g = U2(:, good);

            E1g = E1(good);
            E2g = E2(good);

            px_good = px_g(good);
            py_good = py_g(good);
            eta_good = eta_g(good);

            % ------------------------------------------------
            % Transmission loss consistent with target generation:
            %   eta_eff = eta - TL1(Tx->target)
            %   receiver gain uses only TL2(target->Rx)
            % ------------------------------------------------
            if isfield(settings.cfg, 'enable_TL') && settings.cfg.enable_TL

                pxy = [px_good; py_good];

                tx = settings.cfg.Tx_pos(:);
                rx1 = settings.cfg.Rx_pos{1}(:);
                rx2 = settings.cfg.Rx_pos{2}(:);

                rTx = sqrt(sum((pxy - tx ).^2, 1));
                rRx1 = sqrt(sum((pxy - rx1).^2, 1));
                rRx2 = sqrt(sum((pxy - rx2).^2, 1));

                TL1 = 17*log10(rTx);
                eta_eff = eta_good - TL1;

                G1 = 10.^(-(17*log10(rRx1))/10);
                G2 = 10.^(-(17*log10(rRx2))/10);

            else
                eta_eff = eta_good - settings.cfg.TL1_fixed;
                G1 = 1;
                G2 = 1;
            end

            eta_eff_lin = 10.^(eta_eff/10);

            a1 = sqrt(eta_eff_lin .* G1 ./ E1g);
            a2 = sqrt(eta_eff_lin .* G2 ./ E2g);

            Yt1 = U1g .* a1;
            Yt2 = U2g .* a2;

            R1 = rbg1 - Yt1;
            R2 = rbg2 - Yt2;

            logLt = likelihood_gen(R1, R2, in_ekf);
            logLt = real(logLt(:)).';

            logLt_all(idx_good) = logLt;
        end
    end
end

% ------------------------------------------------------------
% Stable likelihood scaling
% ------------------------------------------------------------
finite_logs = [logLnt, logLt_all(isfinite(logLt_all))];
c = max(finite_logs);

Lt = exp(logLt_all - c);

Lnt = exp(logLnt - c) * ones(1, Nparticles);

end

% ============================================================
% Plot
% ============================================================
function plot_results(settings, data, x, k)

persistent h ax isInitialized

qth = settings.q_threshold;

kplot = 1:k;
ktrue = 1 + (1:k);

% ------------------------------------------------------------
% Estimated states
% ------------------------------------------------------------
xh = data.xh(1, 1:k);
yh = data.xh(2, 1:k);
qk = data.q(1, 1:k);

eta_hat = data.xh(5, 1:k);

indq = qk >= qth;

xh_plot = xh;
yh_plot = yh;
eta_plot = eta_hat;

xh_plot(~indq) = NaN;
yh_plot(~indq) = NaN;
eta_plot(~indq) = NaN;

if indq(end) && isfinite(xh(end)) && isfinite(yh(end))
    xcur = xh(end);
    ycur = yh(end);
else
    xcur = NaN;
    ycur = NaN;
end

% ------------------------------------------------------------
% Particle cloud
% ------------------------------------------------------------
step = 1;
xs = x(:, 1:step:end);

% ------------------------------------------------------------
% Truth trajectory
% ------------------------------------------------------------
have_truth = isfield(settings, 'target_prop') && ...
    iscell(settings.target_prop) && ...
    numel(settings.target_prop) >= 1 && ...
    isfield(settings.target_prop{1}, 'x_target') && ...
    isfield(settings.target_prop{1}, 'y_target') && ...
    numel(settings.target_prop{1}.x_target) >= max(ktrue) && ...
    numel(settings.target_prop{1}.y_target) >= max(ktrue);

if have_truth
    xt = settings.target_prop{1}.x_target(ktrue);
    yt = settings.target_prop{1}.y_target(ktrue);
else
    xt = NaN(1, k);
    yt = NaN(1, k);
end

% ------------------------------------------------------------
% True eta stored in settings.eta_target{1}
% ------------------------------------------------------------
have_eta_true = isfield(settings, 'eta_target') && ...
    iscell(settings.eta_target) && ...
    numel(settings.eta_target) >= 1 && ...
    ~isempty(settings.eta_target{1});

if have_eta_true
    eta_vec = settings.eta_target{1};

    if isscalar(eta_vec)
        eta_true = eta_vec * ones(1, k);
    elseif numel(eta_vec) >= max(ktrue)
        eta_true = eta_vec(ktrue);
    elseif numel(eta_vec) >= k
        eta_true = eta_vec(1:k);
    else
        eta_true = NaN(1, k);
    end
else
    eta_true = NaN(1, k);
end

xreg = settings.x_region;
yreg = settings.y_region;

% ------------------------------------------------------------
% Initialize figure
% ------------------------------------------------------------
need_init = isempty(isInitialized) || ~isInitialized || ...
    ~isfield(h, 'particles') || ~isgraphics(h.particles) || ...
    ~isfield(h, 'qline')     || ~isgraphics(h.qline)     || ...
    ~isfield(h, 'eta_est')   || ~isgraphics(h.eta_est);

if need_init

    figure(1); clf;
    set(gcf, 'Color', 'w');

    ax.ax1 = subplot(3, 6, [1 2 3 7 8 9 13 14 15]);
    hold(ax.ax1, 'on');
    ax.ax2 = subplot(3, 6, [4 5 6]);
    hold(ax.ax2, 'on');
    ax.ax3 = subplot(3, 6, [10 11 12]);
    hold(ax.ax3, 'on');

    set([ax.ax1 ax.ax2 ax.ax3], 'Units', 'normalized');
    set(ax.ax1, 'Position', [0.06 0.11 0.52 0.80]);
    set(ax.ax2, 'Position', [0.64 0.73 0.31 0.18]);
    set(ax.ax3, 'Position', [0.64 0.42 0.31 0.18]);

    % --------------------------------------------------------
    % Trajectory plot
    % --------------------------------------------------------
    h.particles = plot(ax.ax1, xs(1, :), xs(2, :), 'o', ...
        'LineStyle', 'none', ...
        'MarkerSize', 4, ...
        'MarkerFaceColor', 'none', ...
        'MarkerEdgeColor', [0.65 0.65 0.65]);

    h.true_traj = plot(ax.ax1, xt, yt, 'b-', 'LineWidth', 1.5);
    h.est_traj = plot(ax.ax1, xh_plot, yh_plot, 'm-', 'LineWidth', 1.5);
    h.est_cur = plot(ax.ax1, xcur, ycur, 'mo', 'MarkerSize', 8, 'LineWidth', 1.5);

    h.tx = plot(ax.ax1, NaN, NaN, 'rs', 'MarkerSize', 8, 'LineWidth', 1.5);
    h.rx1 = plot(ax.ax1, NaN, NaN, 'gd', 'MarkerSize', 8, 'LineWidth', 1.5);
    h.rx2 = plot(ax.ax1, NaN, NaN, 'cd', 'MarkerSize', 8, 'LineWidth', 1.5);

    tx = settings.cfg.Tx_pos(:);
    set(h.tx, 'XData', tx(1), 'YData', tx(2));

    rx1 = settings.cfg.Rx_pos{1}(:);
    rx2 = settings.cfg.Rx_pos{2}(:);

    set(h.rx1, 'XData', rx1(1), 'YData', rx1(2));
    set(h.rx2, 'XData', rx2(1), 'YData', rx2(2));

    grid(ax.ax1, 'on');
    axis(ax.ax1, 'manual');
    xlim(ax.ax1, xreg);
    ylim(ax.ax1, yreg);

    xlabel(ax.ax1, 'x (m)');
    ylabel(ax.ax1, 'y (m)');
    h.title1 = title(ax.ax1, sprintf('PF cloud (k = %d)', k));

    set(get(ax.ax1, 'XLabel'), 'Units', 'normalized', 'Position', [0.5 -0.06 0]);
    set(get(ax.ax1, 'YLabel'), 'Units', 'normalized', 'Position', [-0.06 0.5 0]);

    legend(ax.ax1, [h.particles h.est_traj h.est_cur h.true_traj h.tx h.rx1 h.rx2], ...
        {'Particles', 'Estimated traj', 'Estimated current', 'True traj', 'Tx', 'Rx 1', 'Rx 2'}, ...
        'Location', 'northeast');

    % --------------------------------------------------------
    % Existence probability
    % --------------------------------------------------------
    h.qline = plot(ax.ax2, kplot, qk, 'k-', 'LineWidth', 1.5);
    h.qth = yline(ax.ax2, qth, 'r--', 'LineWidth', 1.0);

    grid(ax.ax2, 'on');
    ylim(ax.ax2, [0 1.05]);
    xlim(ax.ax2, [1 max(k, 2)]);
    xlabel(ax.ax2, 'Ping index');
    ylabel(ax.ax2, 'q_{k|k}', 'Interpreter', 'tex');
    title(ax.ax2, 'Existence probability');

    % --------------------------------------------------------
    % Eta plot
    % --------------------------------------------------------
    h.eta_true = plot(ax.ax3, kplot, eta_true, 'b-', 'LineWidth', 1.2);
    h.eta_est = plot(ax.ax3, kplot, eta_plot, 'bo', ...
        'LineStyle', 'none', 'MarkerSize', 4);

    grid(ax.ax3, 'on');
    xlim(ax.ax3, [1 max(k, 2)]);
    xlabel(ax.ax3, 'Ping index');
    ylabel(ax.ax3, '\eta_{dB}', 'Interpreter', 'tex');
    title(ax.ax3, 'Target \eta_{dB}');

    legend(ax.ax3, [h.eta_true h.eta_est], ...
        {'True \eta_{dB}', 'Estimated \eta_{dB}'}, ...
        'Location', 'best', 'Interpreter', 'tex');

    set([ax.ax1 ax.ax2 ax.ax3], 'Box', 'on');

    isInitialized = true;
end

% ------------------------------------------------------------
% Update trajectory
% ------------------------------------------------------------
set(h.particles, 'XData', xs(1, :), 'YData', xs(2, :));
set(h.est_traj, 'XData', xh_plot, 'YData', yh_plot);
set(h.est_cur, 'XData', xcur, 'YData', ycur);
set(h.title1, 'String', sprintf('PF cloud (k = %d)', k));

if have_truth
    set(h.true_traj, 'XData', xt, 'YData', yt);
else
    set(h.true_traj, 'XData', NaN, 'YData', NaN);
end

% ------------------------------------------------------------
% Update q
% ------------------------------------------------------------
set(h.qline, 'XData', kplot, 'YData', qk);
xlim(ax.ax2, [1 max(k, 2)]);

% ------------------------------------------------------------
% Update eta
% ------------------------------------------------------------
set(h.eta_true, 'XData', kplot, 'YData', eta_true);
set(h.eta_est, 'XData', kplot, 'YData', eta_plot);
xlim(ax.ax3, [1 max(k, 2)]);

drawnow limitrate nocallbacks;

end
