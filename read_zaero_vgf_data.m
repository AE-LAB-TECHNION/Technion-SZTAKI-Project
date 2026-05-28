function data = read_zaero_vgf_data(filename, numPhysicalModes)
%READ_ZAERO_VGF_DATA Read compact V-G-F data from a ZAERO output file.

if isstring(filename)
    filename = char(filename);
end

if nargin < 2 || isempty(numPhysicalModes)
    numPhysicalModes = [];
end

lines = read_text_lines(filename);
tables = parse_vgf_tables(lines, numPhysicalModes);

if isempty(tables)
    error('No V-G-F summary table was found in ZAERO output file: %s', filename);
end

tableIdx = choose_vgf_table(tables);
data = tables(tableIdx);
data.filename = filename;
data.flutter = find_vgf_flutter(data.airspeed, data.frequencyHz, data.dampingG);

end

% =========================================================================
function lines = read_text_lines(filename)

fid = fopen(filename, 'r');
if fid == -1
    error('Could not open file: %s', filename);
end
cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

lines = {};
while true
    line = fgetl(fid);
    if ~ischar(line)
        break
    end
    lines{end+1, 1} = line; %#ok<AGROW>
end

end

% =========================================================================
function tables = parse_vgf_tables(lines, numPhysicalModes)

tables = empty_vgf_table(0);
tables = tables([]);
i = 1;

while i <= numel(lines)
    tokens = regexp(lines{i}, ...
        ['THE FOLLOWING V-G-F TABLE LISTS\s+(\d+)\s+NUMBER OF STRUCTURAL MODES', ...
         '.*AND\s+(\d+)\s+NUMBER OF AERODYNAMIC LAG ROOTS'], ...
        'tokens', 'once');

    if isempty(tokens)
        i = i + 1;
        continue
    end

    numStructuralModes = str2double(tokens{1});
    numAeroLagRoots = str2double(tokens{2});
    [table, nextIdx] = parse_vgf_table(lines, i, numStructuralModes, ...
        numAeroLagRoots, numPhysicalModes);

    if ~isempty(table.airspeed)
        tables(end+1) = table; %#ok<AGROW>
    end

    i = max(i + 1, nextIdx);
end

end

% =========================================================================
function [table, nextIdx] = parse_vgf_table(lines, startIdx, ...
    numStructuralModes, numAeroLagRoots, numPhysicalModes)

if isempty(numPhysicalModes)
    numModesToRead = numStructuralModes;
else
    numModesToRead = min(numPhysicalModes, numStructuralModes);
end

table = empty_vgf_table(numModesToRead);
table.numStructuralModes = numStructuralModes;
table.numAeroLagRoots = numAeroLagRoots;
table.startLine = startIdx;
table.analysisType = vgf_analysis_type(lines, startIdx);

j = startIdx + 1;
while j <= numel(lines)
    if j > startIdx + 1 && contains(lines{j}, 'THE FOLLOWING V-G-F TABLE LISTS')
        break
    end

    if has_vgf_rows(table) && contains(lines{j}, 'NON-MATCHED POINT FLUTTER ANALYSIS')
        break
    end

    if isempty(regexp(lines{j}, 'MODE NO\.', 'once'))
        j = j + 1;
        continue
    end

    blockModeNumbers = parse_mode_numbers(lines{j});
    [table, j] = parse_vgf_mode_block(lines, j, blockModeNumbers, table);

    if has_all_requested_modes(table)
        break
    end
end

table = sort_vgf_rows(table);
nextIdx = j;

end

% =========================================================================
function table = empty_vgf_table(numModes)

table = struct();
table.filename = '';
table.analysisType = '';
table.startLine = NaN;
table.numStructuralModes = NaN;
table.numAeroLagRoots = NaN;
table.modeNumbers = 1:numModes;
table.airspeed = [];
table.dynPressure = [];
table.frequencyHz = nan(numModes, 0);
table.dampingG = nan(numModes, 0);
table.flutter = empty_flutter_result();

end

% =========================================================================
function [table, nextIdx] = parse_vgf_mode_block(lines, modeHeaderIdx, ...
    blockModeNumbers, table)

headerIdx = [];
for i = modeHeaderIdx+1:numel(lines)
    if contains(lines{i}, 'V/VREF') && contains(lines{i}, 'F(HZ)')
        headerIdx = i;
        break
    end

    if contains(lines{i}, 'THE FOLLOWING V-G-F TABLE LISTS') || ...
            contains(lines{i}, 'NON-MATCHED POINT FLUTTER ANALYSIS')
        nextIdx = i;
        return
    end
end

if isempty(headerIdx)
    nextIdx = modeHeaderIdx + 1;
    return
end

rowIdx = headerIdx + 1;
while rowIdx <= numel(lines)
    values = parse_numeric_tokens(lines{rowIdx});
    numModeGroups = min(blockModeNumbers_per_row(values), numel(blockModeNumbers));

    if numModeGroups == 0 || numel(values) < 3 || isnan(values(1)) || isnan(values(2))
        break
    end

    airspeed = values(2);
    dynPressure = values(3);

    if airspeed > 0
        [table, airspeedIdx] = ensure_airspeed_column(table, airspeed, dynPressure);

        for k = 1:numModeGroups
            modeNumber = blockModeNumbers(k);
            if modeNumber < 1 || modeNumber > numel(table.modeNumbers)
                continue
            end

            valueIdx = 3 + 3*(k - 1) + 1;
            table.dampingG(modeNumber, airspeedIdx) = values(valueIdx);
            table.frequencyHz(modeNumber, airspeedIdx) = values(valueIdx + 1);
        end
    end

    rowIdx = rowIdx + 1;
end

nextIdx = rowIdx;

end

% =========================================================================
function modeNumbers = parse_mode_numbers(line)

tokens = regexp(line, 'MODE NO\.\s*(\d+)', 'tokens');
modeNumbers = zeros(1, numel(tokens));

for i = 1:numel(tokens)
    modeNumbers(i) = str2double(tokens{i}{1});
end

modeNumbers = modeNumbers(isfinite(modeNumbers));

end

% =========================================================================
function values = parse_numeric_tokens(line)

tokens = regexp(strtrim(line), '\S+', 'match');
values = nan(1, numel(tokens));

for i = 1:numel(tokens)
    values(i) = parse_zaero_number(tokens{i});
end

end

% =========================================================================
function value = parse_zaero_number(token)

token = strtrim(token);

if isempty(token) || contains(token, '*')
    value = NaN;
    return
end

if strcmpi(token, 'INFINT')
    value = Inf;
    return
end

token = strrep(token, 'D', 'E');
token = strrep(token, 'd', 'E');
value = str2double(token);

if ~isnan(value)
    return
end

tokens = regexp(token, '^([+-]?[0-9]*\.?[0-9]+)([+-][0-9]+)$', 'tokens', 'once');
if isempty(tokens)
    value = NaN;
else
    value = str2double([tokens{1}, 'E', tokens{2}]);
end

end

% =========================================================================
function n = blockModeNumbers_per_row(values)

n = floor((numel(values) - 3) / 3);

if n < 0
    n = 0;
end

end

% =========================================================================
function [table, airspeedIdx] = ensure_airspeed_column(table, airspeed, dynPressure)

tol = 1e-8 * max(abs(airspeed), 1);

if isempty(table.airspeed)
    airspeedIdx = 1;
else
    [delta, airspeedIdx] = min(abs(table.airspeed - airspeed));
    if delta <= tol
        if isnan(table.dynPressure(airspeedIdx))
            table.dynPressure(airspeedIdx) = dynPressure;
        end
        return
    end

    airspeedIdx = numel(table.airspeed) + 1;
end

table.airspeed(airspeedIdx) = airspeed;
table.dynPressure(airspeedIdx) = dynPressure;
table.dampingG(:, airspeedIdx) = NaN;
table.frequencyHz(:, airspeedIdx) = NaN;

end

% =========================================================================
function tf = has_vgf_rows(table)

tf = ~isempty(table.airspeed);

end

% =========================================================================
function tf = has_all_requested_modes(table)

tf = has_vgf_rows(table) && all(any(~isnan(table.frequencyHz), 2));

end

% =========================================================================
function table = sort_vgf_rows(table)

if isempty(table.airspeed)
    return
end

[table.airspeed, order] = sort(table.airspeed);
table.dynPressure = table.dynPressure(order);
table.frequencyHz = table.frequencyHz(:, order);
table.dampingG = table.dampingG(:, order);

end

% =========================================================================
function analysisType = vgf_analysis_type(lines, startIdx)

analysisType = 'ZAERO';
contextStart = max(1, startIdx - 30);

for i = contextStart:startIdx
    if contains(lines{i}, 'ASE ANALYSIS')
        analysisType = 'ASE';
        return
    end
end

for i = contextStart:startIdx
    if contains(lines{i}, 'FLUTTER ANALYSIS')
        analysisType = 'FLUTTER';
        return
    end
end

end

% =========================================================================
function tableIdx = choose_vgf_table(tables)

tableIdx = 1;
bestScore = table_score(tables(1));

for i = 2:numel(tables)
    score = table_score(tables(i));
    if score > bestScore
        tableIdx = i;
        bestScore = score;
    end
end

end

% =========================================================================
function score = table_score(table)

score = table.numAeroLagRoots;

if strcmp(table.analysisType, 'ASE')
    score = score + 1000;
end

if has_all_requested_modes(table)
    score = score + 100;
end

end

function flutter = find_vgf_flutter(airspeed, frequencyHz, dampingG)

flutter = empty_flutter_result();
bestV = Inf;

for modeIdx = 1:size(dampingG, 1)
    gBranch = dampingG(modeIdx, :);
    fBranch = frequencyHz(modeIdx, :);

    for iV = 1:numel(airspeed)-1
        g1 = gBranch(iV);
        g2 = gBranch(iV+1);
        f1 = fBranch(iV);
        f2 = fBranch(iV+1);

        if any(isnan([g1, g2, f1, f2]))
            continue
        end

        if g1 <= 0 && g2 >= 0 && (g1 < 0 || g2 > 0)
            V1 = airspeed(iV);
            V2 = airspeed(iV+1);

            if g1 == 0
                Vf = V1;
                ff = f1;
            elseif g2 == 0
                Vf = V2;
                ff = f2;
            else
                Vf = V1 - g1*(V2 - V1)/(g2 - g1);
                ff = f1 + (Vf - V1)*(f2 - f1)/(V2 - V1);
            end
        else
            continue
        end

        if Vf < bestV
            bestV = Vf;
            flutter.hasFlutter = true;
            flutter.Vf = Vf;
            flutter.ff = ff;
            flutter.modeIndex = modeIdx;
            flutter.crossingIndex = iV;
        end
    end
end

end

% =========================================================================
function flutter = empty_flutter_result()

flutter = struct();
flutter.hasFlutter = false;
flutter.Vf = NaN;
flutter.ff = NaN;
flutter.modeIndex = NaN;
flutter.crossingIndex = NaN;

end
