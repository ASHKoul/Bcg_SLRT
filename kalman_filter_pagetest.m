function run = kalman_filter_pagetest(Yrun,S,U,B,cfg,Stpl,h,mode)

L   = cfg.L;
N   = cfg.N;
Kb  = size(B,2);
I_N = speye(N);
I_K = speye(Kb);
H   = S*B;
Np  = size(Yrun,2);

if size(Yrun,1) ~= N
    error('Size mismatch: size(Yrun,1)=%d but cfg.N=%d', size(Yrun,1), N);
end
if nargin < 8 || isempty(mode), mode = 2; end

doPage = (nargin >= 6) && ~isempty(Stpl);
if doPage
    if size(Stpl,1) ~= N
        error('Stpl must have %d rows (same as Yrun).', N);
    end
    if isempty(h)
        error('Threshold h must be provided when Stpl is provided.');
    end
    Ktpl = size(Stpl,2);
else
    Ktpl = 0;
end

% Dynamics
F = cfg.Fphi * I_K;
Q = (cfg.sigma_q^2) * I_K;

lam0 = 1e-3;

% ---- sigma_e handling (scalar or vector) ----
if isscalar(cfg.sigma_e)
    sigma_e_vec = [];
else
    sigma_e_vec = cfg.sigma_e(:).';
    if numel(sigma_e_vec) < Np
        error('cfg.sigma_e length %d < Np=%d', numel(sigma_e_vec), Np);
    end
    sigma_e_vec = sigma_e_vec(1:Np);
end

% Allocate outputs (single run)
ell = zeros(1,Np);
T   = zeros(1,Np);
det = zeros(1,Np,'uint8');    % 0=notdetected, 1=monitor, 2=detected

tau_det = Np;
missed  = true;

% ---- INIT at ping 1 ----
y0 = Yrun(:,1);
a0 = (S'*S + lam0*speye(L)) \ (S'*y0);
theta0_prev = (B'*B + 1e-8*speye(Kb)) \ (B'*a0);
P0_prev     = 1e0 * I_K;

Tk = 0;
needReinitBoth = false;

if doPage
    tpl_pos = 1;
    if tpl_pos <= Ktpl
        x1 = Stpl(:,tpl_pos);
    else
        x1 = zeros(N,1,'like',Yrun);
    end
    y1 = Yrun(:,1) - x1;

    a1 = (S'*S + lam0*speye(L)) \ (S'*y1);
    theta1_prev = (B'*B + 1e-8*speye(Kb)) \ (B'*a1);
    P1_prev     = 1e0 * I_K;
    tpl_pos = 2;
end

% ================= MAIN LOOP =================
for k = 2:Np

    if doPage && needReinitBoth
        % restart template and reinit both filters at ping k
        tpl_pos = 1;

        yk0 = Yrun(:,k);
        a0  = (S'*S + lam0*speye(L)) \ (S'*yk0);
        theta0_prev = (B'*B + 1e-8*speye(Kb)) \ (B'*a0);
        P0_prev     = 1e0 * I_K;

        if tpl_pos <= Ktpl
            xk = Stpl(:,tpl_pos);
        else
            xk = zeros(N,1,'like',Yrun);
        end
        yk1 = Yrun(:,k) - xk;

        a1  = (S'*S + lam0*speye(L)) \ (S'*yk1);
        theta1_prev = (B'*B + 1e-8*speye(Kb)) \ (B'*a1);
        P1_prev     = 1e0 * I_K;

        ell(k) = 0;
        T(k)   = Tk;
        det(k) = uint8(0);

        needReinitBoth = false;
        tpl_pos = 2;
        continue
    end

    % -------- choose sigma_e at ping k --------
    if isscalar(cfg.sigma_e)
        sig_e_k = cfg.sigma_e;
    else
        sig_e_k = sigma_e_vec(k);
    end

    % ===================== H0 KF update =====================
    yk = Yrun(:,k);

    theta_pr0 = F*theta0_prev;
    P_pr0     = F*P0_prev*F' + Q;

    a_pr0 = B*theta_pr0;
    a20    = abs(a_pr0).^2;
    Ua0 = U .* sqrt(a20(:)).';
    Rk0 = sig_e_k^2 * I_N ...
        + cfg.sigma_c^2 * (Ua0 * Ua0') ...
        + cfg.sigma_d^2 * ((U*a_pr0)*(U*a_pr0)');
    Syy0 = H*P_pr0*H' + Rk0;
    Syy0 = 0.5*(Syy0 + Syy0');

    [LS0,p0] = chol(Syy0,'lower');
    if p0
        [LS0,~] = chol(Syy0 + 1e-12*eye(size(Syy0,1),'like',Syy0),'lower');
    end

    Kg0    = ((P_pr0*H')/LS0')/LS0;
    innov0 = yk - H*theta_pr0;

    theta0 = theta_pr0 + Kg0*innov0;
    P0     = (I_K - Kg0*H)*P_pr0;

    w0 = LS0\innov0;
    NIS0 = real(w0'*w0);
    logdetS0 = 2*sum(log(diag(LS0)));
    NLL0 = (NIS0 + logdetS0 + cfg.N*log(pi));

    theta0_prev = theta0;
    P0_prev     = P0;

    % ===================== H1 KF update + Page =====================
    if doPage
        if tpl_pos <= Ktpl
            xk = Stpl(:,tpl_pos);
        else
            xk = zeros(N,1,'like',Yrun);
        end
        yk1 = Yrun(:,k) - xk;

        theta_pr1 = F*theta1_prev;
        P_pr1     = F*P1_prev*F' + Q;
        P_pr1 = 0.5*(P_pr1 + P_pr1');

        a_pr1 = B*theta_pr1;
        a21    = abs(a_pr1).^2;
        Ua1 = U .* sqrt(a21(:)).';

        Rk1 = sig_e_k^2 * I_N ...
            + cfg.sigma_c^2 * (Ua1 * Ua1') ...
            + cfg.sigma_d^2 * ((U*a_pr1)*(U*a_pr1)');
        Rk1 = 0.5*(Rk1 + Rk1');

        Syy1 = H*P_pr1*H' + Rk1;
        Syy1 = 0.5*(Syy1 + Syy1');

        [LS1,p1] = chol(Syy1,'lower');
        if p1
            [LS1,~] = chol(Syy1 + 1e-12*eye(size(Syy1,1),'like',Syy1),'lower');
        end

        Kg1    = ((P_pr1*H')/LS1')/LS1;
        innov1 = yk1 - H*theta_pr1;

        theta1 = theta_pr1 + Kg1*innov1;
        P1     = (I_K - Kg1*H)*P_pr1;
        P1 = 0.5*(P1 + P1');

        w1 = LS1\innov1;
        NIS1 = real(w1'*w1);
        logdetS1 = 2*sum(log(diag(LS1)));
        NLL1 = (NIS1 + logdetS1 + cfg.N*log(pi));

        theta1_prev = theta1;
        P1_prev     = P1;

        ell_k  = NLL0 - NLL1;
        Tk_new = Tk + ell_k;
        Tk     = max(0, Tk_new);

        ell(k) = ell_k;
        T(k)   = Tk;

        if Tk == 0
            det(k) = uint8(0);
            needReinitBoth = true;   % if you keep the reinit design
        elseif Tk >= h
            det(k) = uint8(2);
            if missed
                tau_det = k;
                missed  = false;
            end

            if mode == 1
                if k < Np
                    det(k+1:Np) = uint8(2);
                    ell(k+1:Np) = 0;
                    T(k+1:Np)   = Tk;
                end
                break
            end

            tpl_pos = tpl_pos + 1;
        else
            det(k) = uint8(1);
            tpl_pos = tpl_pos + 1;
        end
    end
end

% Package single-run outputs
run.ell     = ell;
run.T       = T;
run.det     = det;
run.tau_det = tau_det;
run.missed  = missed;

% run.page.tpl_idx = tpl_idx_hist;
% run.page.ell_k   = ell;
% run.page.T_hist  = T;
% run.page.det     = det;
end
