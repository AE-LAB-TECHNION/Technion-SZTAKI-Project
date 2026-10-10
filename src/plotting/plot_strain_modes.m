function plot_strain_modes(cfg, inputs)
%PLOT_STRAIN_MODES Plot diagnostic strain mode shapes from DAT centroids and PSI.

datFile = cfg_value(cfg, 'datFile', 'MODEL-0012.dat');
strainElemIDs = cfg.sensor.strainElemIDs;
PSI = inputs.PSI;
selectedModes = cfg.selected_modes;
normalizeModes = cfg_value(cfg, 'normalizeStrainModes', true);

fprintf('\n============================================================\n');
fprintf('Diagnostic strain mode plots\n');
fprintf('============================================================\n');
fprintf('DAT file: %s\n', datFile);
fprintf('exist(datFile,''file'') = %d\n', exist(datFile, 'file'));
fprintf('PSI size = %d x %d\n', size(PSI, 1), size(PSI, 2));
fprintf('Number of strainElemIDs = %d\n', numel(strainElemIDs));

if ~exist(datFile, 'file')
    warning('DAT file was not found: %s. Strain mode plots skipped.', datFile);
    return
end

if isempty(PSI)
    warning('PSI is empty. Strain mode plots skipped.');
    return
end

[strainX, strainY, ~, strainElemInfo] = extract_element_centroids_from_dat(datFile, strainElemIDs);
nFound = sum(~isnan(strainX) & ~isnan(strainY));

fprintf('Requested strain elements : %d\n', numel(strainElemIDs));
fprintf('Found element coordinates : %d\n', nFound);

if nFound == 0
    warning('No strain element coordinates were found. Check the DAT element IDs.');
    return
end

if numel(strainX) ~= size(PSI, 1)
    warning('Number of strain coordinates does not match number of PSI rows.');
end

for im = 1:min(numel(selectedModes), size(PSI, 2))
    modeNumber = selectedModes(im);
    strainMode = PSI(:, im);
    valid = ~isnan(strainX) & ~isnan(strainY) & ~isnan(strainMode);

    fprintf('Mode %d: valid strain points = %d\n', modeNumber, sum(valid));

    if ~any(valid)
        warning('No valid strain points for mode %d. Skipping.', modeNumber);
        continue
    end

    x = strainX(valid);
    y = strainY(valid);
    z = strainMode(valid);

    if normalizeModes
        maxAbsZ = max(abs(z));
        if maxAbsZ > 0
            zPlot = z / maxAbsZ;
            zLabelText = 'Normalized modal strain';
        else
            zPlot = z;
            zLabelText = 'Modal strain';
        end
    else
        zPlot = z;
        zLabelText = 'Modal strain';
    end

    fig3d = figure('Name', sprintf('3D Strain Mode Shape - Mode %d', modeNumber));
    apply_plot_style(fig3d, cfg);
    scatter3(x, y, zPlot, 35, zPlot, 'filled');
    grid on; box on;
    view(45, 25);
    xlabel('X [model units]');
    ylabel('Y [model units]');
    zlabel(zLabelText);
    title(sprintf('3D Strain Mode Shape - Mode %d', modeNumber));
    cb = colorbar;
    cb.Label.String = zLabelText;
    apply_plot_style(fig3d, cfg);

    figTop = figure('Name', sprintf('Top View Strain Mode - Mode %d', modeNumber));
    apply_plot_style(figTop, cfg);
    scatter(x, y, 40, zPlot, 'filled');
    grid on; box on;
    axis equal;
    xlabel('X [model units]');
    ylabel('Y [model units]');
    title(sprintf('Top View of Strain Mode - Mode %d', modeNumber));
    cb = colorbar;
    cb.Label.String = zLabelText;
    apply_plot_style(figTop, cfg);
end

fprintf('\nFirst extracted strain element coordinates:\n');
disp(strainElemInfo(1:min(10, height(strainElemInfo)), :));

end

function [elemX, elemY, elemZ, elemInfo] = extract_element_centroids_from_dat(datFile, elemIDs)

fid = fopen(datFile, 'r');
if fid < 0
    error('Could not open DAT file: %s', datFile);
end
cleaner = onCleanup(@() fclose(fid));

lines = {};
while true
    line = fgetl(fid);
    if ~ischar(line)
        break
    end
    lines{end+1, 1} = line; %#ok<AGROW>
end

gridMap = containers.Map('KeyType', 'double', 'ValueType', 'any');
elemMap = containers.Map('KeyType', 'double', 'ValueType', 'any');

for i = 1:numel(lines)
    line = strtrim(lines{i});
    if isempty(line) || startsWith(line, '$')
        continue
    end

    [fields, ok] = parse_nastran_fields(line);
    if ~ok || isempty(fields)
        continue
    end

    card = upper(fields{1});

    if startsWith(card, 'GRID') && numel(fields) >= 6
        gid = nastran_num(fields{2});
        x = nastran_num(fields{4});
        y = nastran_num(fields{5});
        z = nastran_num(fields{6});
        if ~isnan(gid) && ~isnan(x) && ~isnan(y)
            if isnan(z)
                z = 0;
            end
            gridMap(gid) = [x y z];
        end
    elseif any(strcmp(card, {'CROD', 'CQUAD4', 'CTRIA3'}))
        [eid, gids] = element_connectivity(card, fields);
        if ~isnan(eid) && all(~isnan(gids))
            elemMap(eid) = gids;
        end
    end
end

elemIDs = elemIDs(:);
nElem = numel(elemIDs);
elemX = nan(nElem, 1);
elemY = nan(nElem, 1);
elemZ = nan(nElem, 1);
found = false(nElem, 1);
numNodes = nan(nElem, 1);

for k = 1:nElem
    eid = elemIDs(k);
    if ~isKey(elemMap, eid)
        continue
    end

    gids = elemMap(eid);
    numNodes(k) = numel(gids);
    coords = nan(numel(gids), 3);

    for j = 1:numel(gids)
        gid = gids(j);
        if isKey(gridMap, gid)
            coords(j, :) = gridMap(gid);
        end
    end

    validXY = ~isnan(coords(:, 1)) & ~isnan(coords(:, 2));
    if any(validXY)
        elemX(k) = mean(coords(validXY, 1));
        elemY(k) = mean(coords(validXY, 2));
        elemZ(k) = mean(coords(validXY, 3));
        found(k) = true;
    end
end

elemInfo = table(elemIDs, found, numNodes, elemX, elemY, elemZ, ...
    'VariableNames', {'ElementID', 'Found', 'NumNodes', 'X', 'Y', 'Z'});

end

function [eid, gids] = element_connectivity(card, fields)

eid = NaN;
gids = NaN;

switch card
    case 'CROD'
        if numel(fields) >= 5
            eid = nastran_num(fields{2});
            gids = [nastran_num(fields{4}), nastran_num(fields{5})];
        end
    case 'CQUAD4'
        if numel(fields) >= 7
            eid = nastran_num(fields{2});
            gids = [nastran_num(fields{4}), nastran_num(fields{5}), ...
                nastran_num(fields{6}), nastran_num(fields{7})];
        end
    case 'CTRIA3'
        if numel(fields) >= 6
            eid = nastran_num(fields{2});
            gids = [nastran_num(fields{4}), nastran_num(fields{5}), nastran_num(fields{6})];
        end
end

end

function [fields, ok] = parse_nastran_fields(line)

ok = true;
line = char(line);
idxComment = strfind(line, '$');
if ~isempty(idxComment)
    line = line(1:idxComment(1)-1);
end
line = strtrim(line);

if isempty(line)
    fields = {};
    ok = false;
    return
end

if contains(line, ',')
    fields = strtrim(strsplit(line, ','));
    return
end

rawSpace = regexp(line, '\s+', 'split');
if numel(rawSpace) > 1
    fields = strtrim(rawSpace);
    fields = fields(~cellfun(@isempty, fields));
    return
end

n = length(line);
nFields = ceil(n / 8);
fields = cell(1, nFields);
for k = 1:nFields
    i1 = (k - 1) * 8 + 1;
    i2 = min(k * 8, n);
    fields{k} = strtrim(line(i1:i2));
end
fields = fields(~cellfun(@isempty, fields));

if isempty(fields)
    ok = false;
end

end

function val = nastran_num(s)

if isnumeric(s)
    val = s;
    return
end

s = strtrim(char(s));
if isempty(s)
    val = NaN;
    return
end

val = str2double(s);
if ~isnan(val)
    return
end

s = regexprep(s, '[Dd]', 'e');
val = str2double(s);
if ~isnan(val)
    return
end

if isempty(regexpi(s,'[Ee]','once'))
    tok = regexp(s, '^([+-]?\d*\.?\d*)([+-]\d+)$', 'tokens', 'once');
    if ~isempty(tok)
        val = str2double([tok{1} 'e' tok{2}]);
        return
    end
end

val = NaN;

end

function value = cfg_value(cfg, fieldName, defaultValue)

if isfield(cfg, fieldName)
    value = cfg.(fieldName);
else
    value = defaultValue;
end

end
