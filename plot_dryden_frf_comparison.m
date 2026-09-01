function plot_dryden_frf_comparison(cfg,comparison)
%PLOT_DRYDEN_FRF_COMPARISON Compare MATLAB and ZAERO complex gust FRFs.

positive = comparison.frequency > 0;
nModes = comparison.nModes;
colors = lines(2*nModes);
fig = figure('Name','Dryden FRF comparison');
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

labels = strings(1,2*nModes);
nexttile; hold on;
for k = 1:nModes
    semilogy(comparison.frequency(positive), ...
        abs(comparison.matlab(positive,k)),'-', ...
        'Color',colors(2*k-1,:),'LineWidth',cfg.plotLineWidth);
    semilogy(comparison.frequency(positive), ...
        abs(comparison.zaero(positive,k)),'--', ...
        'Color',colors(2*k,:),'LineWidth',cfg.plotLineWidth);
    labels(2*k-1) = sprintf('MATLAB \\xi_%d',k);
    labels(2*k) = sprintf('ZAERO \\xi_%d',k);
end
grid on; box on; xlim([0 cfg.drydenFrequencyEnd]);
ylabel('|FRF|'); title('Dryden gust-to-modal-displacement FRF');
legend(labels,'Location','eastoutside','NumColumns',1);

nexttile; hold on;
for k = 1:nModes
    plot(comparison.frequency(positive), ...
        unwrap(angle(comparison.matlab(positive,k)))*180/pi,'-', ...
        'Color',colors(2*k-1,:),'LineWidth',cfg.plotLineWidth);
    plot(comparison.frequency(positive), ...
        unwrap(angle(comparison.zaero(positive,k)))*180/pi,'--', ...
        'Color',colors(2*k,:),'LineWidth',cfg.plotLineWidth);
end
grid on; box on; xlim([0 cfg.drydenFrequencyEnd]);
xlabel('Frequency [Hz]'); ylabel('Phase [deg]');
legend(labels,'Location','eastoutside','NumColumns',1);
apply_plot_style(fig,cfg);
end
