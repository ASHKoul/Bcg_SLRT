function post_ekf=meas_update_ekf(y1k, y2k,prior_ekf,use_bg,q, qth,k,kth)
if use_bg && q < qth
    post_ekf = posterior_out_ekf(y1k, y2k,prior_ekf);
else
   post_ekf = prior_ekf;
end
end