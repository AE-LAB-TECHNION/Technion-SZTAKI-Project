function plot_dryden_frequency_psd(cfg,response)
%PLOT_DRYDEN_FREQUENCY_PSD Plot MATLAB modal PSDs without ZAERO overlays.

positive = response.frequency > 0;
colors = lines(response.Nel);
fig = figure('Name','Dryden modal-response PSD');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile; hold on;
for k = 1:response.Nel
    loglog(response.frequency(positive),response.modalXPsd(positive,k), ...
        'Color',colors(k,:),'LineWidth',cfg.plotLineWidth);
end
grid on; box on; xlim([0.1 100]); xlabel('Frequency [Hz]');
ylabel('PSD of \xi_i [1/Hz]'); title('Modal displacement');
displacementLabels = arrayfun(@(k) sprintf('\\xi_%d',k), ...
    1:response.Nel,'UniformOutput',false);
legend(displacementLabels,'Location','northeast');

nexttile; hold on;
for k = 1:response.Nel
    loglog(response.frequency(positive),response.modalVPsd(positive,k), ...
        'Color',colors(k,:),'LineWidth',cfg.plotLineWidth);
end
grid on; box on; xlim([0.1 100]); xlabel('Frequency [Hz]');
ylabel('PSD of d\xi_i/dt [1/s^2/Hz]'); title('Modal velocity');
velocityLabels = arrayfun(@(k) sprintf('d\\xi_%d/dt',k), ...
    1:response.Nel,'UniformOutput',false);
legend(velocityLabels,'Location','northeast');
apply_plot_style(fig,cfg);
end
