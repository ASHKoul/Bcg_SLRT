function prior_out_ekf= Time_update_ekf(in_ekf,settings,use_bg,q,qth,k,kth)

if use_bg
   % if q < qth
        prior_out_ekf = prior_background(settings.H, settings.U, settings.B, settings.cfg, in_ekf,settings);
   % else
       % prior_out_ekf = in_ekf;
   % end
else
    prior_out_ekf.theta_hat_1 = [];
    prior_out_ekf.theta_hat_2 = [];
    prior_out_ekf.P_hat_1     = [];
    prior_out_ekf.P_hat_2     = [];
    prior_out_ekf.H           = [];
    prior_out_ekf.Sigma_1     = settings.cfg.sigma_e2;
    prior_out_ekf.Sigma_2     = settings.cfg.sigma_e2;
end

end

