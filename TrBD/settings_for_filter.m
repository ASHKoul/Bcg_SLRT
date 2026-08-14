function s_new = settings_for_filter(s, model)

s_new = s;
s_new.H = s.S*s.B;
s_new.modelLQ = chol(model.Q);

s_new.cfg.sigma_e2 = s.flt_nois_inf_const*s.s2e;
s_new.cfg.sigma_e = sqrt(s_new.cfg.sigma_e2);

s_new = rmfield(s_new, 'S');
s_new = rmfield(s_new, 'lfm');
s_new = rmfield(s_new, 'lfm_bb');
s_new = rmfield(s_new, 't_lfm_bb');

end
