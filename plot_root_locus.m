function plot_root_locus(cfg, eigTraj, plotModelName)
%PLOT_ROOT_LOCUS Plot eigenvalue trajectories colored by airspeed.

fig = figure(cfg.rootLocusFigure);
clf(fig);
ax = axes(fig);
hold(ax, 'on');
grid(ax, 'on');
box(ax, 'on');

for j = 1:numel(cfg.airspeed)
    Vj = cfg.airspeed(j);

    scatter(ax, real(eigTraj(:,j)), imag(eigTraj(:,j)), ...
        25, Vj*ones(size(eigTraj,1),1), 'filled');
end

colormap(ax, jet)
cb = colorbar(ax);
cb.Label.String = 'Airspeed [m/s]';

xlabel(ax, 'Real(\lambda)');
ylabel(ax, 'Imag(\lambda)');

title(ax, sprintf('Root Locus of %s Eigenvalues', plotModelName));

xlim(ax, cfg.rootLocusXLim)
ylim(ax, cfg.rootLocusYLim)

apply_plot_style(fig, cfg);

end
