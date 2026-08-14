function [target_signal, a_target_true] = target_signal_genration(s, tau_abs, beta)

cfg = s.cfg;
N = cfg.N;

target_signal = zeros(N, 1);
a_target_true = 0;

t_rel = s.t_rec;

tau_direct = cfg.initial_delay;
if ~(isfinite(tau_abs) && isfinite(beta) && beta > 0)
    return;
end

tau = tau_abs - tau_direct;
if ~(isfinite(tau) && tau >= 0)
    return;
end

tp = beta * (t_rel - tau);
mask = (tp >= 0) & (tp <= cfg.Tp);

u = sym_chirp_bb_eval(tp, s);
u = u(:);
u(~mask) = 0;

Eu = sum(abs(u).^2);
if ~(isfinite(Eu) && Eu > 0)
    return;
end

u = u / sqrt(Eu);

target_signal = s.amplitude * u;
a_target_true = s.amplitude;

end
