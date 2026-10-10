function [inputs, cfg] = load_analysis_inputs(cfg)
%LOAD_ANALYSIS_INPUTS Read RFA, ZAERO modal, and sensor/modal data.

validate_input_files(cfg);

RFA_mat = extractRFAmatrices(cfg.RFA_filename, cfg.method);
zaeroModalInfo = read_zaero_modal_info(cfg.ZAERO_filename);
validate_omitted_modes_from_zaero(zaeroModalInfo, cfg.ZAERO_filename);

cfg.Nel = derive_retained_mode_count(RFA_mat, zaeroModalInfo);
cfg.numControlSurfaces = derive_count_from_zaero( ...
    zaeroModalInfo.numControlSurfaceModes, cfg, 'numControlSurfaces', 0, ...
    'control-surface modes');
cfg.numGustInputs = derive_gust_count_from_zaero(zaeroModalInfo, cfg, RFA_mat);
cfg.isgust = double(cfg.numGustInputs > 0);
cfg.selected_modes = derive_selected_modes(zaeroModalInfo.femModeNumbers, ...
    zaeroModalInfo.omittedModes, cfg.Nel);

eigvals = modal_values_for_modes(zaeroModalInfo.eigenvaluesByMode, ...
    cfg.selected_modes, 'eigenvalues');

RBM = load_zaero_rbm_measurement( ...
    cfg.rbmModelFile, cfg.rbmOutputRow, cfg.Nel);

fprintf('ZAERO modal source file: %s\n', cfg.ZAERO_filename);
fprintf('Number of FEM modes retained by ZAERO: %d\n', cfg.Nel);
fprintf('F06 mode numbers omitted by ZAERO: ');
if isempty(zaeroModalInfo.omittedModes)
    fprintf('none');
else
    fprintf('%d ', zaeroModalInfo.omittedModes);
end
fprintf('\n');

fprintf('F06 mode numbers used in aeroelastic model: ');
fprintf('%d ', cfg.selected_modes);
fprintf('\n');
fprintf('Number of control-surface modes in ZAERO RFA: %d\n', cfg.numControlSurfaces);
fprintf('Number of gust vectors used for RFA partition: %d\n', cfg.numGustInputs);
if ~isnan(zaeroModalInfo.numLoadModes)
    fprintf('Number of LOADMOD rows in ZAERO RFA: %d\n', zaeroModalInfo.numLoadModes);
end

freqHz = sqrt(eigvals(:)) / (2*pi);

fprintf('Frequencies used [Hz]:\n');
fprintf('  %-8s %16s\n', 'F06 Mode', 'Frequency');
for i = 1:numel(eigvals)
    fprintf('  %-8d %16.6f\n', cfg.selected_modes(i), freqHz(i));
end

sensor = cfg.sensor;

[PSI, PHI, PHI_ROT, sensorInfo] = load_sensor_modal_matrices(cfg, sensor);

fprintf('\nSensor/modal matrix sizes:\n');
fprintf('  PSI     : %d x %d\n', size(PSI,1),     size(PSI,2));
fprintf('  PHI     : %d x %d\n', size(PHI,1),     size(PHI,2));
fprintf('  PHI_ROT : %d x %d\n', size(PHI_ROT,1), size(PHI_ROT,2));
fprintf('  RBM     : %d x %d\n', ...
    size(RBM.modalCoefficients,1), size(RBM.modalCoefficients,2));

if any(isnan(PSI(:))) || any(isnan(PHI(:))) || any(isnan(PHI_ROT(:)))
    warning(['At least one of PSI/PHI/PHI_ROT contains NaN. ', ...
             'Check that the requested modes, GRID IDs, and element IDs exist in the F06.']);
end

inputs = struct();
inputs.RFA_mat = RFA_mat;
inputs.zaeroModalInfo = zaeroModalInfo;
inputs.all_eigs = zaeroModalInfo.eigenvaluesByMode;
inputs.eigvals = eigvals;
inputs.selected_modes = cfg.selected_modes;
inputs.omittedModes = zaeroModalInfo.omittedModes;
inputs.PSI = PSI;
inputs.PHI = PHI;
inputs.PHI_ROT = PHI_ROT;
inputs.RBM = RBM;
inputs.sensorInfo = sensorInfo;

end

% =========================================================================
function validate_input_files(cfg)

validate_input_file(cfg.RFA_filename, 'cfg.RFA_filename');
validate_input_file(cfg.ZAERO_filename, 'cfg.ZAERO_filename');
validate_input_file(cfg.f06file, 'cfg.f06file');
validate_input_file(cfg.rbmModelFile, 'cfg.rbmModelFile');

if cfg.useSensorCache && ~cfg.rebuildSensorCache && ...
        isfield(cfg, 'sensorCacheFile') && ~isempty(cfg.sensorCacheFile)
    validate_input_file(cfg.sensorCacheFile, 'cfg.sensorCacheFile');
end

end

% =========================================================================
function validate_input_file(filename, fieldName)

if isstring(filename)
    filename = char(filename);
end

[~, basename, ext] = fileparts(filename);
if isempty(basename) && isempty(ext)
    error('%s is empty.', fieldName);
end

if exist(filename, 'file') ~= 2
    error('Required input file %s was not found:\n%s',fieldName,filename);
end

end

% =========================================================================
function validate_omitted_modes_from_zaero(zaeroModalInfo, filename)

if zaeroModalInfo.hasOmitmodCard && isempty(zaeroModalInfo.omittedModes)
    error(['ZAERO output contains an OMITMOD card, but no omitted mode numbers could be parsed.\n', ...
           'Check the OMITMOD formatting in: %s'], filename);
end

end

% =========================================================================
function Nel = derive_retained_mode_count(RFA_mat, zaeroModalInfo)

rfaModeCount = size(RFA_mat.A0, 1);
zaeroModeCount = zaeroModalInfo.numFEMModes;

if isnan(zaeroModeCount)
    warning(['Could not read NUMBER OF FEM MODES from the ZAERO output. ', ...
             'Using the RFA A0 row count instead.']);
    Nel = rfaModeCount;
    return
end

if size(RFA_mat.A0, 1) < zaeroModeCount || size(RFA_mat.A0, 2) < zaeroModeCount
    error(['ZAERO reports %d retained FEM modes, but RFA matrix A0 is only %d x %d. ', ...
           'Check that cfg.ZAERO_filename and cfg.RFA_filename come from the same run.'], ...
           zaeroModeCount, size(RFA_mat.A0, 1), size(RFA_mat.A0, 2));
end

if ~isnan(zaeroModalInfo.numLoadModes)
    expectedRows = zaeroModeCount + zaeroModalInfo.numLoadModes;
    if rfaModeCount ~= expectedRows
        warning(['ZAERO reports NM=%d and NLM=%d, so A0 is expected to have %d rows; ', ...
                 'RFA_mat.A0 has %d rows. Using the first NM=%d structural rows.'], ...
                 zaeroModeCount, zaeroModalInfo.numLoadModes, expectedRows, ...
                 rfaModeCount, zaeroModeCount);
    end
elseif rfaModeCount ~= zaeroModeCount
    warning(['ZAERO reports NM=%d retained FEM modes, while RFA_mat.A0 has %d rows. ', ...
             'Extra rows are treated as non-structural rows and ignored by the state equations.'], ...
             zaeroModeCount, rfaModeCount);
end

Nel = zaeroModeCount;

end

% =========================================================================
function value = derive_count_from_zaero(zaeroValue, cfg, fieldName, defaultValue, label)

hasCfgValue = isfield(cfg, fieldName) && ~isempty(cfg.(fieldName));

if ~isnan(zaeroValue)
    value = zaeroValue;
    if hasCfgValue && cfg.(fieldName) ~= zaeroValue
        warning(['ZAERO reports %d %s, but cfg.%s is %d. ', ...
                 'Using the ZAERO value.'], ...
                 zaeroValue, label, fieldName, cfg.(fieldName));
    end
    return
end

if hasCfgValue
    value = cfg.(fieldName);
else
    value = defaultValue;
end

end

% =========================================================================
function numGustInputs = derive_gust_count_from_zaero(zaeroModalInfo, cfg, RFA_mat)

numGustModes = zaeroModalInfo.numGustModes;
numColumns = size(RFA_mat.A0, 2);
columnsWithoutGust = cfg.Nel + cfg.numControlSurfaces;
extraColumns = max(numColumns - columnsWithoutGust, 0);
hasCfgValue = isfield(cfg, 'numGustInputs') && ~isempty(cfg.numGustInputs);

% The RFA matrix file is authoritative for the input partition. An
% extra column after the retained structural and control-surface columns is
% a gust vector even when the textual ZAERO report does not label it.
if extraColumns > 0
    numGustInputs = extraColumns;

    if ~zaeroModalInfo.hasGustInput
        warning(['The RFA file contains %d extra input column(s) after ', ...
            'the structural and control columns, although the ZAERO report ', ...
            'does not explicitly identify gust vectors. Treating the extra ', ...
            'column(s) as gust inputs.'],numGustInputs);
    elseif ~isnan(numGustModes) && numGustModes ~= numGustInputs
        warning(['ZAERO reports %d gust mode(s), but the RFA file has %d extra input column(s). ', ...
                 'Using the RFA column count for gust vectors.'], ...
                 numGustModes, numGustInputs);
    end
    if hasCfgValue && cfg.numGustInputs ~= numGustInputs
        warning(['The RFA file indicates %d gust vector(s), but cfg.numGustInputs is %d. ', ...
                 'Using the RFA value.'], numGustInputs, cfg.numGustInputs);
    end
    return
end

if ~isnan(numGustModes)
    numGustInputs = numGustModes;
    if hasCfgValue && cfg.numGustInputs ~= numGustInputs
        warning(['ZAERO indicates %d gust vector(s), but cfg.numGustInputs is %d. ', ...
                 'Using the ZAERO value.'], numGustInputs, cfg.numGustInputs);
    end
    return
end

if hasCfgValue
    numGustInputs = cfg.numGustInputs;
else
    numGustInputs = 0;
end

end
% =========================================================================
function selectedModes = derive_selected_modes(femModeNumbers, omittedModes, numRetainedModes)

availableModes = femModeNumbers(:).';
availableModes = availableModes(~ismember(availableModes, omittedModes));

if numel(availableModes) < numRetainedModes
    error(['Could not derive %d retained FEM modes from the ZAERO modal table. ', ...
           'Only %d non-omitted FEM modes were found.'], ...
           numRetainedModes, numel(availableModes));
end

selectedModes = availableModes(1:numRetainedModes);

end

% =========================================================================
function values = modal_values_for_modes(valuesByMode, modeNumbers, valueName)

if max(modeNumbers) > numel(valuesByMode) || any(isnan(valuesByMode(modeNumbers)))
    error('Could not find %s for all derived FEM modes.', valueName);
end

values = valuesByMode(modeNumbers);

end

% =========================================================================
function [PSI, PHI, PHI_ROT, sensorInfo] = load_sensor_modal_matrices(cfg, sensor)

cacheInfo = make_sensor_cache_info(cfg, sensor);

if cfg.useSensorCache && ~cfg.rebuildSensorCache && exist(cfg.sensorCacheFile, 'file') == 2
    cached = load(cfg.sensorCacheFile);

    if is_sensor_cache_match(cached, cacheInfo)
        PSI = cached.PSI;
        PHI = cached.PHI;
        PHI_ROT = cached.PHI_ROT;

        if isfield(cached, 'sensorInfo')
            sensorInfo = cached.sensorInfo;
        else
            sensorInfo = cached.info;
        end

        fprintf('\nLoaded sensor/modal matrices from cache: %s\n', cfg.sensorCacheFile);
        return
    end
end

[PSI, PHI, PHI_ROT, sensorInfo] = extract_PSI_PHI_PHIROT_from_F06( ...
    cfg.f06file, cfg.selected_modes, ...
    sensor.dispGridIDs, sensor.dispComp, ...
    sensor.rotGridIDs, sensor.rotComp, ...
    sensor.strainElemIDs);

if cfg.useSensorCache
    info = sensorInfo;
    save(cfg.sensorCacheFile, 'PSI', 'PHI', 'PHI_ROT', 'sensorInfo', 'info', 'cacheInfo');
    fprintf('\nSaved sensor/modal matrices to cache: %s\n', cfg.sensorCacheFile);
end

end

% =========================================================================
function tf = is_sensor_cache_match(cached, expected)

hasMatrices = isfield(cached, 'PSI') && isfield(cached, 'PHI') && isfield(cached, 'PHI_ROT');
hasInfo = isfield(cached, 'sensorInfo') || isfield(cached, 'info');

if ~hasMatrices || ~hasInfo
    tf = false;
    return
end

if isfield(cached, 'cacheInfo')
    tf = isequal(cached.cacheInfo, expected);
    return
end

if isfield(cached, 'sensorInfo')
    info = cached.sensorInfo;
else
    info = cached.info;
end

tf = isfield(info, 'f06file') && strcmp(info.f06file, expected.f06file) && ...
     isfield(info, 'modesToUse') && isequal(info.modesToUse, expected.selected_modes) && ...
     isfield(info, 'dispGridIDs') && isequal(info.dispGridIDs, expected.dispGridIDs) && ...
     isfield(info, 'dispComp') && strcmp(info.dispComp, expected.dispComp) && ...
     isfield(info, 'rotGridIDs') && isequal(info.rotGridIDs, expected.rotGridIDs) && ...
     isfield(info, 'rotComp') && strcmp(info.rotComp, expected.rotComp) && ...
     isfield(info, 'strainElemIDs') && isequal(info.strainElemIDs, expected.strainElemIDs);

end

% =========================================================================
function cacheInfo = make_sensor_cache_info(cfg, sensor)

cacheInfo = struct();
cacheInfo.f06file = cfg.f06file;
cacheInfo.selected_modes = cfg.selected_modes;
cacheInfo.dispGridIDs = sensor.dispGridIDs;
cacheInfo.dispComp = sensor.dispComp;
cacheInfo.rotGridIDs = sensor.rotGridIDs;
cacheInfo.rotComp = sensor.rotComp;
cacheInfo.strainElemIDs = sensor.strainElemIDs;

fileInfo = dir(cfg.f06file);
if ~isempty(fileInfo)
    cacheInfo.f06Bytes = fileInfo.bytes;
    cacheInfo.f06Datenum = fileInfo.datenum;
else
    cacheInfo.f06Bytes = NaN;
    cacheInfo.f06Datenum = NaN;
end

end
