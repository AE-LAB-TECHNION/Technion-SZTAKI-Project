function apply_plot_style(fig, cfg)
%APPLY_PLOT_STYLE Apply shared typography to project figures.

if nargin < 1 || isempty(fig)
    fig = gcf;
end
if nargin < 2 || isempty(cfg)
    cfg = struct();
end

fontName = plot_style_value(cfg, 'plotFontName', 'Helvetica');
fontSize = plot_style_value(cfg, 'plotFontSize', 13);
titleFontSize = plot_style_value(cfg, 'plotTitleFontSize', fontSize + 1);
titleFontWeight = plot_style_value(cfg, 'plotTitleFontWeight', 'normal');
interpreter = plot_style_value(cfg, 'plotTextInterpreter', 'tex');

set(findall(fig, '-property', 'FontName'), 'FontName', fontName);
set(findall(fig, '-property', 'FontSize'), 'FontSize', fontSize);
set(findall(fig, '-property', 'Interpreter'), 'Interpreter', interpreter);
set(findall(fig, '-property', 'TickLabelInterpreter'), ...
    'TickLabelInterpreter', interpreter);

axesHandles = findall(fig, 'Type', 'axes');
for ax = reshape(axesHandles, 1, [])
    set(ax.Title, ...
        'FontName', fontName, ...
        'FontSize', titleFontSize, ...
        'FontWeight', titleFontWeight, ...
        'Interpreter', interpreter);
end

titleHandles = findall(fig, 'Type', 'text');
for textHandle = reshape(titleHandles, 1, [])
    tag = get(textHandle, 'Tag');
    if ~isempty(strfind(lower(tag), 'title'))
        set(textHandle, ...
            'FontName', fontName, ...
            'FontSize', titleFontSize, ...
            'FontWeight', titleFontWeight, ...
            'Interpreter', interpreter);
    end
end

end

function value = plot_style_value(cfg, fieldName, defaultValue)

if isfield(cfg, fieldName)
    value = cfg.(fieldName);
else
    value = defaultValue;
end

end
