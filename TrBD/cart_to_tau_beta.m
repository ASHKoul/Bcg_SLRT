function [tau, beta, geom] = cart_to_tau_beta(pT, vT, cfg, irec)

pTx = cfg.Tx_pos(:);
pRx = cfg.Rx_pos{irec}(:);

% Geometry vectors
dTx = pT  - pTx; % Tx -> target
dRx = pT -  pRx; % Rx -> target

% Ranges
rTx = vecnorm(dTx, 2, 1);
rRx = vecnorm(dRx, 2, 1);

% LOS unit vectors
uTx = dTx ./ rTx;
uRx = dRx ./ rRx;

% Bistatic range and delay
rBi = rTx + rRx;
tau = rBi / cfg.c;

% Bistatic range-rate
Rdot = sum(uTx .* vT, 1) + sum(uRx .* vT, 1);

% Wideband scaling (first-order)
beta = 1 - (Rdot / cfg.c);

% Optional diagnostics
if nargout > 2
    geom = struct();
    geom.dTx = dTx;
    geom.dRx = dRx;
    geom.uTx = uTx;
    geom.uRx = uRx;
    geom.rTx = rTx;
    geom.rRx = rRx;
    geom.rBi = rBi;
    geom.Rdot = Rdot;
end
end
