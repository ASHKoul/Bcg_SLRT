function [Yb, Atrue, s] = power_calibrations(Yb, Atrue, s)

INR_lin = 10^(s.cfg.INRdB / 10);

Neff = round(s.cfg.max_delay_spread * s.cfg.Fs_used);

Yseg_before = Yb(1:Neff, :);

Eb = sum(abs(Yseg_before).^2, 1);
Eb_des = Neff* s.s2e * INR_lin;

scale_ping = sqrt(Eb_des ./ Eb);
scale_ping(~isfinite(scale_ping)) = 1;
Yb = Yb .* scale_ping;
Atrue = Atrue .* scale_ping;

end
