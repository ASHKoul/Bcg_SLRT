function y = sym_chirp_passband(t, s)

phase = 2*pi*((s.cfg.fc - s.cfg.BW/2).*t + 0.5*s.cfg.mu.*t.^2);
y = cos(phase) .* (t >= 0 & t <= s.cfg.Tp);
end
