function idx_keep = index_finding_theta(theta)

p = abs(theta(:)).^2;
pmax = max(p);

if pmax <= 0
    idx_keep = 1;
    return
end

mask = p >= 0.04 * pmax;

if ~any(mask)
    [~, imax] = max(p);
    mask(imax) = true;
end

d = diff([false; mask; false]);
i1 = find(d == 1);
i2 = find(d == -1) - 1;

idx_keep = zeros(numel(i1), 1);

for ic = 1:numel(i1)
    idx_seg = i1(ic):i2(ic);
    [~, loc] = max(abs(theta(idx_seg)).^2);
    idx_keep(ic) = idx_seg(loc);
end

end
