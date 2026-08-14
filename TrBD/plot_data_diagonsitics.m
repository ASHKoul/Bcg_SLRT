function plot_data_diagonsitics(Yb_syn, data, Ymat1, Ymat2)

Np = data.s.cfg.Np;

figure(100); clf;
tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Background-only synthetic data
nexttile;
imagesc(local_db_norm(Yb_syn));
set(gca, 'YDir', 'normal');
xlabel('Delay sample');
ylabel('Ping index');
title('Synthetic Background Only');
colorbar;
clim([-30 0]);

% Clean target signal Rx-1
nexttile;
imagesc(local_db_norm(data.Yt_clean{1}));
set(gca, 'YDir', 'normal');
xlabel('Delay sample');
ylabel('Ping index');
title('Clean Target Signal Rx-1');
colorbar;
clim([-30 0]);

% Clean target signal Rx-2
nexttile;
imagesc(local_db_norm(data.Yt_clean{2}));
set(gca, 'YDir', 'normal');
xlabel('Delay sample');
ylabel('Ping index');
title('Clean Target Signal Rx-2');
colorbar;
clim([-30 0]);

% True target amplitude
nexttile;
plot(1:Np, data.target_true_amplitude{1}, 'LineWidth', 1.5); hold on;
plot(1:Np, data.target_true_amplitude{2}, '--', 'LineWidth', 1.5);
grid on;
xlabel('Ping index');
ylabel('Amplitude');
title('True Target Amplitude');
legend('Receiver-1', 'Receiver-2', 'Location', 'best');

% Received signal Rx-1
nexttile;
imagesc(local_db_norm(Ymat1));
set(gca, 'YDir', 'normal');
xlabel('Delay sample');
ylabel('Ping index');
title('Received Signal Rx-1');
colorbar;
clim([-30 0]);

% Received signal Rx-2
nexttile;
imagesc(local_db_norm(Ymat2));
set(gca, 'YDir', 'normal');
xlabel('Delay sample');
ylabel('Ping index');
title('Received Signal Rx-2');
colorbar;
clim([-30 0]);

% Measured SINR
nexttile;
plot(1:Np, data.s.SINRdB_measured(:, 1), 'LineWidth', 1.5); hold on;
plot(1:Np, data.s.SINRdB_measured(:, 2), '--', 'LineWidth', 1.5);
grid on;
xlabel('Ping index');
ylabel('SINR [dB]');
title('Measured SINR');
legend('Receiver-1', 'Receiver-2', 'Location', 'best');

% Effective target SNR
nexttile;
plot(1:2, 10*log10(data.SNR_eff(:)), 'o-', 'LineWidth', 1.5);
grid on;
xlim([1 2]);
xticks([1 2]);
xticklabels({'Rx-1', 'Rx-2'});
ylabel('Effective SNR [dB]');
title('Effective Target SNR');

sgtitle('Data Diagnostics');

fprintf('Effective target SNR receiver 1: %.2f dB\n', 10*log10(data.SNR_eff(1)));
fprintf('Effective target SNR receiver 2: %.2f dB\n', 10*log10(data.SNR_eff(2)));

% Received signal Rx-1
Fs = 5000; % or use s.cfg.Fs / settings.cfg.Fs if available

YdB = local_db_norm(Ymat1);

Ndelay = size(Ymat1, 1);
Np = size(Ymat1, 2);

delay_axis_ms = 0.667+(0:Ndelay-1)/Fs; % delay in ms
ping_axis = 1:Np;

figure;
imagesc(delay_axis_ms, ping_axis, YdB);
set(gca, 'YDir', 'normal');

xlabel('Delay time [ms]');
ylabel('Ping index');
title('Received Signal Rx-1');

cb = colorbar;
ylabel(cb, 'Normalized power [dB]');

clim([-30 0]);

end

function YdB = local_db_norm(Y)
YdB = 20*log10(abs(Y).' / max(abs(Y(:))) + eps);
end
