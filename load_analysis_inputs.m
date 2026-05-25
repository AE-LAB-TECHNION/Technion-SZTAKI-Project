function inputs = load_analysis_inputs(cfg)
%LOAD_ANALYSIS_INPUTS Read RFA, eigenvalue, and sensor/modal data.

RFA_mat = extractRFAmatrices(cfg.RFA_filename, cfg.method);

all_eigs = read_zaero_eigs(cfg.ZAERO_filename, cfg.nmodes);

if max(cfg.selected_modes) > numel(all_eigs)
    error('selected_modes contains an index larger than the number of eigenvalues read.');
end

eigvals = all_eigs(cfg.selected_modes);

fprintf('Selected modes: ');
fprintf('%d ', cfg.selected_modes);
fprintf('\n');

freqHz = sqrt(eigvals(:)) / (2*pi);

fprintf('Frequencies used [Hz]:\n');
fprintf('  %-6s %16s\n', 'Mode', 'Frequency');
for i = 1:numel(eigvals)
    fprintf('  %-6d %16.6f\n', cfg.selected_modes(i), freqHz(i));
end

sensor = cfg.sensor;

[PSI, PHI, PHI_ROT, sensorInfo] = load_sensor_modal_matrices(cfg, sensor);

fprintf('\nSensor/modal matrix sizes:\n');
fprintf('  PSI     : %d x %d\n', size(PSI,1),     size(PSI,2));
fprintf('  PHI     : %d x %d\n', size(PHI,1),     size(PHI,2));
fprintf('  PHI_ROT : %d x %d\n', size(PHI_ROT,1), size(PHI_ROT,2));

if any(isnan(PSI(:))) || any(isnan(PHI(:))) || any(isnan(PHI_ROT(:)))
    warning(['At least one of PSI/PHI/PHI_ROT contains NaN. ', ...
             'Check that the requested modes, GRID IDs, and element IDs exist in the F06.']);
end

inputs = struct();
inputs.RFA_mat = RFA_mat;
inputs.all_eigs = all_eigs;
inputs.eigvals = eigvals;
inputs.PSI = PSI;
inputs.PHI = PHI;
inputs.PHI_ROT = PHI_ROT;
inputs.sensorInfo = sensorInfo;

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
    info = sensorInfo; %#ok<NASGU>
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
