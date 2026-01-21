function filter_out = kalman_filter_bellhop_basis(Y,S,U,B,cfg)

% ==== Setup ====
L = cfg.L;
N = cfg.N;
K = size(B,2);
I_N = speye(N);
I_K = speye(K);
H   = S*B;
Np=size(Y,2);





% Outputs
filter_out.Ahat         = complex(zeros(L, Np, 'like', Y));
filter_out.Phat         = complex(zeros(L, L,  Np, 'like', Y));
filter_out.theta_hat    = complex(zeros(K, Np, 'like', Y));
filter_out.P_theta_hat  = complex(zeros(K, K,  Np, 'like', Y));
filter_out.Yhat         = complex(zeros(N, Np, 'like', Y));
filter_out.Syy         = cell(1, Np);
filter_out.NIS         = zeros(1, Np);
filter_out.NLL_k       = zeros(1,Np);
% Dynamics
F = cfg.Fphi * I_K;
Q = (cfg.sigma_q^2) * I_K;


% ---- LS init (ping 1) ----
lam0  = 1e-3;
a0    = (S'*S + lam0*speye(cfg.L)) \ (S'*Y(:,1));
theta_prev = (B'*B + 1e-8*speye(size(B,2))) \ (B'*a0);
P_prev = 1e0 * I_K;

% store ping-1 as initialized (not a proper likelihood update)
% filter_out.theta_hat(:,1)    = theta_prev;
% filter_out.P_theta_hat(:,:,1)= P_prev;
% filter_out.Ahat(:,1)         = B*theta_prev;
% filter_out.Phat(:,:,1)       = B*P_prev*B';
% filter_out.Yhat(:,1)         = H*theta_prev;
% filter_out.NLL_k(1)          = NaN;

kStart = 2;

cumLL  = 0;            % sum of log-likelihood (can be > 0)
cumNLL = 0;            % sum of negative log-likelihood (can be < 0)
% logdet_vec = zeros(1, Np);

for k = kStart:Np
    yk = Y(:,k);

    % ---- Predict in theta-space ----
    theta_pr = F*theta_prev;
    P_pr     = F*P_prev*F' + Q;
    P_pr     = 0.5*(P_pr + P_pr');

    % ---- Induce tap prior moments ----
    a_pr = B*theta_pr;

    % ---- Measurement covariance Rk ----
    a2  = abs(a_pr).^2;
    Ua = U .* sqrt(a2(:)).';           
                       
    % aa2 = a_pr*a_pr';
    if isscalar(cfg.sigma_e)
        Rk = cfg.sigma_e^2 * I_N ...
            + cfg.sigma_c^2 * (Ua * Ua') ...
            + cfg.sigma_d^2 * ((U*a_pr)*(U*a_pr)');
    else

        Rk = cfg.sigma_e(k)^2 * I_N ...
            + cfg.sigma_c^2 * (Ua * Ua') ...
            + cfg.sigma_d^2 * ((U*a_pr)*(U*a_pr)');
    end
   % Rk  = 0.5*(Rk + Rk');

    % ---- Innovation covariance ----
    Syy = H*P_pr*H' + Rk;
    Syy = 0.5*(Syy + Syy');

    m = size(Syy,1);
    jitter = 1e-12 + 1e-10*real(trace(Syy))/m;   % scale-aware start
    [LS,p] = chol(Syy + jitter*eye(m,'like',Syy),'lower');
    while p ~= 0
        jitter = jitter*10;
        [LS,~] = chol(Syy + jitter*eye(m,'like',Syy),'lower');
    end

    Kgain = ((P_pr*H')/LS')/LS;
    % Kgain = P_pr*H'/Syy;
    innov = yk - H*theta_pr;

    % ---- Update ----
    theta_po = theta_pr + Kgain*innov;
    P_po = (I_K - Kgain*H)*P_pr*(I_K - Kgain*H)' + Kgain*Rk*Kgain';
    P_po = 0.5*(P_po + P_po');
    %P_po     = 0.5*(P_po + P_po');

    % ---- NLL_k (complex CN) ----
    w = LS\innov;
    NISk = real(w'*w);
    logdetS = 2*sum(log(real(diag(LS))));
    NLL_k = (NISk + logdetS + cfg.N*log(pi));
    LL_k    = -NLL_k;
    cumLL  = cumLL  + LL_k;
    cumNLL = cumNLL + NLL_k;


    % filter_out.NLL_k(k) = NLL_k;
    % filter_out.NIS(k)   = NISk;
    % filter_out.Syy{k}   = Syy;
    % logdet_vec(k) = logdetS;
    % 
    % % ---- Save ----
    filter_out.theta_hat(:,k)     = theta_po;
    filter_out.P_theta_hat(:,:,k) = P_po;
    % filter_out.Ahat(:,k)          = B*theta_po;
    % filter_out.Phat(:,:,k)        = B*P_po*B';
    % filter_out.Yhat(:,k)          = H*theta_po;

    theta_prev = theta_po;
    P_prev     = P_po;
end


% ---- Summaries
% filter_out.logdetSyy = logdet_vec(:);
% filter_out.loglik_sum = cumLL;          % can be > 0 for complex data
filter_out.NLL        = cumNLL;         % can be < 0 for complex data



end
