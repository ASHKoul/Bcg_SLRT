function [Ymat, s] = generate_received_signal_with_background(Yb, Yt, s)

Ymat = cellfun(@(yb, yt) yb + fillmissing(yt, 'constant', 0), Yb, Yt, 'UniformOutput', false);

Yb_mc1 = cat(3, Yb{1, :});
Yt_mc1 = cat(3, Yt{1, :});
Yt_mc1 = fillmissing(Yt_mc1, 'constant', 0);

mask = abs(Yt_mc1) > 0;

Et = squeeze(sum(abs(Yt_mc1).^2, 1)); % Np x N_Rec
Ebn = squeeze(sum(abs(Yb_mc1).^2 .* mask, 1)); % Np x N_Rec

SINRdB_measured = NaN(size(Et));

active = squeeze(any(mask, 1)) & Ebn > 0;

SINRdB_measured(active) = 10*log10(Et(active) ./ Ebn(active));

s.SINRdB_measured = SINRdB_measured;

end
