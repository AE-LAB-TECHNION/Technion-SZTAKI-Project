function info = read_zaero_modal_info(filename)
%READ_ZAERO_MODAL_INFO Read FEM modal metadata from a ZAERO output file.

if isstring(filename)
    filename = char(filename);
end

lines = read_text_lines(filename);

modalTable = parse_fem_modal_table(lines, filename);
omittedModes = parse_omitted_modes(lines);
numFEMModes = parse_num_fem_modes(lines);
numControlSurfaceModes = parse_named_mode_count(lines, 'CONTROL SURFACE');
numGustModes = parse_named_mode_count(lines, 'GUST');
numLoadModes = parse_named_mode_count(lines, 'LOADMOD');
hasGustInput = parse_has_gust_input(lines);
hasOmitmodCard = parse_has_omitmod_card(lines);

femModeNumbers = modalTable(:, 1).';
eigenvalues = modalTable(:, 3);
freqHz = modalTable(:, 5);

info = struct();
info.filename = filename;
info.numFEMModes = numFEMModes;
info.numControlSurfaceModes = numControlSurfaceModes;
info.numGustModes = numGustModes;
info.numLoadModes = numLoadModes;
info.hasGustInput = hasGustInput;
info.hasOmitmodCard = hasOmitmodCard;
info.omittedModes = omittedModes;
info.modalTable = modalTable;
info.femModeNumbers = femModeNumbers;
info.eigenvalues = eigenvalues;
info.freqHz = freqHz;
info.eigenvaluesByMode = values_by_mode(femModeNumbers, eigenvalues);
info.freqHzByMode = values_by_mode(femModeNumbers, freqHz);

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
function modalTable = parse_fem_modal_table(lines, filename)

targetRB = 'RIGID BODY DEGREES OF FREEDOM (DEFINED IN THE FEM BASIC COORDINATE SYSTEM)';
rbIdx = find_first_line(lines, targetRB, 1);
if isempty(rbIdx)
    error('Rigid-body DOF line not found in file %s.', filename);
end

headerIdx = [];
for i = rbIdx+1:numel(lines)
    if contains(lines{i}, 'EIGENVALUE') && contains(lines{i}, 'FREQUENCY')
        headerIdx = i;
        break
    end
end

if isempty(headerIdx)
    error('FEM modal eigenvalue table header not found in file %s.', filename);
end

modalTable = [];
tableStarted = false;

for i = headerIdx+1:numel(lines)
    line = normalize_exponents(lines{i});
    nums = sscanf(line, '%f');

    if numel(nums) >= 7 && is_integer_value(nums(1)) && is_integer_value(nums(2))
        modalTable(end+1, :) = nums(1:7).'; %#ok<AGROW>
        tableStarted = true;
        continue
    end

    if tableStarted && isempty(nums)
        break
    end
end

if isempty(modalTable)
    error('No FEM modal eigenvalue rows were found in file %s.', filename);
end

end

% =========================================================================
function numFEMModes = parse_num_fem_modes(lines)

numFEMModes = NaN;

patterns = { ...
    'NUMBER OF FEM MODES=\s*(\d+)', ...
    'THE FOLLOWING V-G-F TABLE LISTS\s+(\d+)\s+NUMBER OF STRUCTURAL MODES' ...
};

for p = 1:numel(patterns)
    for i = 1:numel(lines)
        tokens = regexp(lines{i}, patterns{p}, 'tokens', 'once');
        if ~isempty(tokens)
            numFEMModes = str2double(tokens{1});
            return
        end
    end
end

end

% =========================================================================
function count = parse_named_mode_count(lines, modeName)

count = NaN;
pattern = ['NUMBER OF\s+', modeName, '\s+MODES=\s*(\d+)'];

for i = 1:numel(lines)
    tokens = regexp(lines{i}, pattern, 'tokens', 'once');
    if ~isempty(tokens)
        count = str2double(tokens{1});
        return
    end
end

end

% =========================================================================
function hasGustInput = parse_has_gust_input(lines)

hasGustInput = false;
gustCards = {'GENGUST', 'DGUST', 'CGUST'};

for i = 1:numel(lines)
    for k = 1:numel(gustCards)
        if contains(lines{i}, gustCards{k})
            hasGustInput = true;
            return
        end
    end
end

end

% =========================================================================
function omittedModes = parse_omitted_modes(lines)

reportModes = parse_omitted_modes_report(lines);
cardModes = parse_omitmod_card(lines);

omittedModes = unique([reportModes, cardModes], 'stable');

end

% =========================================================================
function omittedModes = parse_omitted_modes_report(lines)

omittedModes = [];

for i = 1:numel(lines)
    if ~isempty(regexp(lines{i}, 'OMITTED\s+MODES\s*=', 'once'))
        tokens = regexp(lines{i}, 'OMITTED\s+MODES\s*=\s*(.*)$', 'tokens', 'once');
        if isempty(tokens)
            omittedText = text_after_equals(lines{i});
        else
            omittedText = tokens{1};
        end
        j = i + 1;

        while j <= numel(lines) && is_numeric_continuation_line(lines{j})
            omittedText = [omittedText, ' ', lines{j}]; %#ok<AGROW>
            j = j + 1;
        end

        omittedModes = [omittedModes, parse_integer_list(omittedText)]; %#ok<AGROW>
    end
end

omittedModes = unique(omittedModes, 'stable');

end

% =========================================================================
function omittedModes = parse_omitmod_card(lines)

omittedModes = [];

for i = 1:numel(lines)
    if isempty(regexp(lines{i}, '\bOMITMOD\b', 'once'))
        continue
    end

    omittedModes = [omittedModes, parse_integer_list(text_after_marker(lines{i}, 'OMITMOD'))]; %#ok<AGROW>

    j = i + 1;
    while j <= numel(lines)
        line = lines{j};
        if contains(line, '+OMT')
            omittedModes = [omittedModes, parse_integer_list(text_after_marker(line, '+OMT'))]; %#ok<AGROW>
            j = j + 1;
        elseif isempty(strtrim(line))
            j = j + 1;
        else
            break
        end
    end
end

omittedModes = unique(omittedModes, 'stable');

end

% =========================================================================
function hasOmitmodCard = parse_has_omitmod_card(lines)

hasOmitmodCard = false;

for i = 1:numel(lines)
    if ~isempty(regexp(lines{i}, '\bOMITMOD\b', 'once'))
        hasOmitmodCard = true;
        return
    end
end

end

% =========================================================================
function idx = find_first_line(lines, pattern, startIdx)

idx = [];
for i = startIdx:numel(lines)
    if contains(lines{i}, pattern)
        idx = i;
        return
    end
end

end

% =========================================================================
function text = text_after_marker(line, marker)

idx = strfind(line, marker);
if isempty(idx)
    text = line;
else
    text = line(idx(end)+length(marker):end);
end

end

% =========================================================================
function text = text_after_equals(line)

eqIdx = strfind(line, '=');
if isempty(eqIdx)
    text = line;
else
    text = line(eqIdx(end)+1:end);
end

end

% =========================================================================
function nums = parse_integer_list(text)

matches = regexp(text, '[-+]?\d+', 'match');
if isempty(matches)
    nums = [];
else
    nums = str2double(matches);
    nums = nums(isfinite(nums));
    nums = round(nums);
end

end

% =========================================================================
function tf = is_numeric_continuation_line(line)

trimmed = strtrim(line);
tf = ~isempty(trimmed) && ~isempty(regexp(trimmed, '^[0-9,\s]+$', 'once'));

end

% =========================================================================
function out = normalize_exponents(line)

out = strrep(line, 'D', 'E');
out = strrep(out, 'd', 'E');

end

% =========================================================================
function tf = is_integer_value(value)

tf = isfinite(value) && abs(value - round(value)) < 1e-9;

end

% =========================================================================
function valuesByMode = values_by_mode(modeNumbers, values)

valuesByMode = nan(max(modeNumbers), 1);
valuesByMode(modeNumbers) = values;

end
