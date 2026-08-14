function y = sym_chirp_bb_eval(t, s)

if s.cfg.baseband
    y = exp(1j*2*pi*((-(s.cfg.BW/2)).*t + 0.5*s.cfg.mu.*t.^2)) .* ...
        (t >= 0 & t <= s.cfg.Tp);
else
    y = exp(1j*2*pi*((s.cfg.fc-(s.cfg.BW/2)).*t + 0.5*s.cfg.mu.*t.^2)) .* ...
        (t >= 0 & t <= s.cfg.Tp);
end
end
