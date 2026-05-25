function inputs = load_analysis_inputs(cfg)
%LOAD_ANALYSIS_INPUTS Read eigenvalues and sensor/modal matrices.

all_eigs = read_zaero_eigs(cfg.ZAERO_filename, cfg.nmodes);

if max(cfg.selected_modes) > numel(all_eigs)
    error('selected_modes contains an index larger than the number of eigenvalues read.');
end

eigvals = all_eigs(cfg.selected_modes);

fprintf('Selected modes: ');
fprintf('%d ', cfg.selected_modes);
fprintf('\n');

fprintf('Eigenvalues used:\n');
disp(eigvals(:).');

sensor = cfg.sensor;

[PSI, PHI, PHI_ROT, sensorInfo] = extract_PSI_PHI_PHIROT_from_F06( ...
    cfg.f06file, cfg.selected_modes, ...
    sensor.dispGridIDs, sensor.dispComp, ...
    sensor.rotGridIDs, sensor.rotComp, ...
    sensor.strainElemIDs);

fprintf('\nSensor/modal matrix sizes:\n');
fprintf('  PSI     : %d x %d\n', size(PSI,1),     size(PSI,2));
fprintf('  PHI     : %d x %d\n', size(PHI,1),     size(PHI,2));
fprintf('  PHI_ROT : %d x %d\n', size(PHI_ROT,1), size(PHI_ROT,2));

if any(isnan(PSI(:))) || any(isnan(PHI(:))) || any(isnan(PHI_ROT(:)))
    warning(['At least one of PSI/PHI/PHI_ROT contains NaN. ', ...
             'Check that the requested modes, GRID IDs, and element IDs exist in the F06.']);
end

inputs = struct();
inputs.all_eigs = all_eigs;
inputs.eigvals = eigvals;
inputs.PSI = PSI;
inputs.PHI = PHI;
inputs.PHI_ROT = PHI_ROT;
inputs.sensorInfo = sensorInfo;

end
