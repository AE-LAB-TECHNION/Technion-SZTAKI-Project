function save_analysis_figures(cfg)
%SAVE_ANALYSIS_FIGURES Save all currently open named figures as PNG files.

if ~cfg.output.saveFigures
    return
end

outputDirectory = cfg.output.figureDirectory;
if ~isfolder(outputDirectory)
    mkdir(outputDirectory);
end

figures = findall(groot, 'Type', 'figure');
for fig = reshape(figures, 1, [])
    name = strtrim(get(fig, 'Name'));
    if isempty(name)
        name = sprintf('figure_%d', fig.Number);
    end
    filename = regexprep(lower(name), '[^a-z0-9]+', '_');
    filename = regexprep(filename, '^_|_$', '');
    exportgraphics(fig, fullfile(outputDirectory, [filename '.png']), ...
        'Resolution', cfg.output.figureResolution);
end

end
