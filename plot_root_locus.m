function plot_root_locus(cfg, eigTraj, plotModelName)
%PLOT_ROOT_LOCUS Plot eigenvalue trajectories colored by airspeed.

figure(cfg.rootLocusFigure); clf; hold on; grid on; box on;

for j = 1:numel(cfg.airspeed)
    Vj = cfg.airspeed(j);

    scatter(real(eigTraj(:,j)), imag(eigTraj(:,j)), ...
        25, Vj*ones(size(eigTraj,1),1), 'filled');
end

colormap(jet)
cb = colorbar;
cb.Label.String = 'Airspeed [m/s]';

xlabel('Real$(\lambda)$','Interpreter','latex');
ylabel('Imag$(\lambda)$','Interpreter','latex');

title(sprintf('Root Locus of %s Eigenvalues', plotModelName), ...
      'Interpreter','latex');

xlim(cfg.rootLocusXLim)
ylim(cfg.rootLocusYLim)

set(gca, 'FontSize', 13);

end
