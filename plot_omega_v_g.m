function [Vf, ff, data] = plot_omega_v_g(cfg, ~, ~)
%PLOT_OMEGA_V_G Plot physical modal frequency and damping versus airspeed.

data = read_zaero_vgf_data(cfg.ZAERO_filename, cfg.Nel);
fprintf('\nomega-V-g plot uses ZAERO V-G-F data from %s.\n', cfg.ZAERO_filename);

Vf = data.flutter.Vf;
ff = data.flutter.ff;

figNumber = cfg_value(cfg, 'omegaVGFigure', cfg_value(cfg, 'rootLocusFigure', 101) + 1);
fig = figure(figNumber);
clf(fig);
set(fig, 'Name', 'omega-V-g');

figure(fig);
axFreq = subplot(2, 1, 1);
hold(axFreq, 'on');
grid(axFreq, 'on');
box(axFreq, 'on');

axDamp = subplot(2, 1, 2);
hold(axDamp, 'on');
grid(axDamp, 'on');
box(axDamp, 'on');

airspeed = data.airspeed;
numModes = size(data.frequencyHz, 1);
colors = lines(numModes);
labels = mode_labels(cfg, numModes);

for modeIdx = 1:numModes
    plot(axFreq, airspeed, data.frequencyHz(modeIdx, :), ...
        'Color', colors(modeIdx, :), ...
        'Marker', 'o', ...
        'DisplayName', labels{modeIdx});

    plot(axDamp, airspeed, data.dampingG(modeIdx, :), ...
        'Color', colors(modeIdx, :), ...
        'Marker', 'o', ...
        'DisplayName', labels{modeIdx});
end

plot(axDamp, [min(airspeed), max(airspeed)], [0, 0], ...
    'k:', 'HandleVisibility', 'off');

ylabel(axFreq, 'Frequency [Hz]');
ylabel(axDamp, 'Damping g [-]');
xlabel(axDamp, 'Airspeed [m/s]');
ylim(axDamp, cfg_value(cfg, 'omegaVGDampingYLim', [-0.5, 0.5]));

if data.flutter.hasFlutter
    add_flutter_marker(axFreq, Vf, ff);
    add_flutter_marker(axDamp, Vf, 0);
    add_flutter_label(axFreq, Vf, ff, airspeed);
else
    fprintf('\nomega-V-g plot: no damping zero crossing found.\n');
end

linkaxes([axFreq, axDamp], 'x');
xlim(axFreq, [min(airspeed), max(airspeed)]);

apply_plot_style(fig, cfg);

end

% =========================================================================
function add_flutter_marker(ax, Vf, yValue)

yLimits = ylim(ax);
plot(ax, [Vf, Vf], yLimits, 'k--', ...
    'LineWidth', 1.2, ...
    'HandleVisibility', 'off');
plot(ax, Vf, yValue, 'ko', ...
    'MarkerFaceColor', 'k', ...
    'MarkerSize', 4, ...
    'HandleVisibility', 'off');
ylim(ax, yLimits);

end

% =========================================================================
function add_flutter_label(ax, Vf, ff, airspeed)

xLimits = [min(airspeed), max(airspeed)];
yLimits = ylim(ax);
xRange = max(xLimits(2) - xLimits(1), eps);
yRange = max(yLimits(2) - yLimits(1), eps);

if Vf > xLimits(1) + 0.65*xRange
    xText = Vf - 0.02*xRange;
    horizontalAlignment = 'right';
else
    xText = Vf + 0.02*xRange;
    horizontalAlignment = 'left';
end

yText = min(max(ff, yLimits(1) + 0.12*yRange), yLimits(2) - 0.12*yRange);

text(ax, xText, yText, sprintf('Flutter: V = %.3f m/s, f = %.3f Hz', Vf, ff), ...
    'Color', 'k', ...
    'FontWeight', 'normal', ...
    'HorizontalAlignment', horizontalAlignment, ...
    'VerticalAlignment', 'middle');

end

% =========================================================================
function labels = mode_labels(cfg, numModes)

labels = cell(1, numModes);

if isfield(cfg, 'selected_modes') && numel(cfg.selected_modes) >= numModes
    for modeIdx = 1:numModes
        labels{modeIdx} = sprintf('F06 mode %d', cfg.selected_modes(modeIdx));
    end
else
    for modeIdx = 1:numModes
        labels{modeIdx} = sprintf('Mode %d', modeIdx);
    end
end

end

% =========================================================================
function value = cfg_value(cfg, fieldName, defaultValue)

if isfield(cfg, fieldName) && ~isempty(cfg.(fieldName))
    value = cfg.(fieldName);
else
    value = defaultValue;
end

end
