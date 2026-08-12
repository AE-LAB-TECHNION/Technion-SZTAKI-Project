function peakTable = plot_step_response_checks(cfg, inputs, models)
%PLOT_STEP_RESPONSE_CHECKS Plot plant step responses for each control surface.

peakTable = table();

[model, iV] = select_model_at_speed(cfg, models, cfg_value(cfg, 'responseCheckSpeed', 35));

if ~strcmpi(model.modelType, 'plant') || isempty(model.Ap) || isempty(model.Bp)
    fprintf('\nStep response checks skipped because the selected model is not a plant model.\n');
    return
end

[sys, outputIndex] = build_plant_output_system(cfg, inputs, model);
Nctrl = size(model.Bp, 2);

time = cfg_value(cfg, 'responseStepTime', 0:0.001:2);
deltaCommand = deg2rad(cfg_value(cfg, 'responseStepCommandDeg', 1));
nStrongStrain = cfg_value(cfg, 'responseStrongestStrainCount', 5);

peakRate = zeros(Nctrl, 1);
peakDisp = zeros(Nctrl, 1);
peakStrain = zeros(Nctrl, 1);
normBpCol = zeros(Nctrl, 1);

fprintf('\n============================================================\n');
fprintf('State-space model for response checks\n');
fprintf('============================================================\n');
fprintf('Selected V_check = %.2f m/s (cfg.airspeed index %d)\n', model.V, iV);
fprintf('Number of states  = %d\n', size(model.Ap, 1));
fprintf('Number of inputs  = %d\n', size(model.Bp, 2));
fprintf('Number of outputs = %d\n', size(sys.C, 1));

for csID = 1:Nctrl
    fprintf('\n============================================================\n');
    fprintf('Step response for CS %d\n', csID);
    fprintf('============================================================\n');

    inputCommand = zeros(numel(time), Nctrl);
    inputCommand(:, csID) = deltaCommand;

    [y, tout] = lsim(sys, inputCommand, time);

    peakRate(csID) = max(abs(y(:, outputIndex.rate)), [], 'all');
    peakDisp(csID) = max(abs(y(:, outputIndex.disp)), [], 'all');
    peakStrain(csID) = max(abs(y(:, outputIndex.strain)), [], 'all');
    normBpCol(csID) = norm(model.Bp(:, csID));

    plot_output_group(cfg, tout, y(:, outputIndex.rate), ...
        sprintf('Step response: CS %d command to rate outputs', csID), ...
        'Time [s]', 'Rate outputs', 'Rate output');

    plot_output_group(cfg, tout, y(:, outputIndex.disp), ...
        sprintf('Step response: CS %d command to disp/acc outputs', csID), ...
        'Time [s]', 'Disp/acc outputs', 'Disp/acc output');

    strainResponse = y(:, outputIndex.strain);
    strainPeakEachOutput = max(abs(strainResponse), [], 1);
    [~, sortedStrain] = sort(strainPeakEachOutput, 'descend');
    selectedLocal = sortedStrain(1:min(nStrongStrain, numel(sortedStrain)));
    selectedGlobal = outputIndex.strain(selectedLocal);

    strainElemIDs = cfg.sensor.strainElemIDs;
    selectedStrainElemIDs = strainElemIDs(selectedLocal);

    plot_output_group(cfg, tout, 1e6*y(:, selectedGlobal), ...
    sprintf('Step response: CS %d command to strongest strain outputs', csID), ...
    'Time [s]', 'Strain outputs x 10^6', 'Strain element ID', selectedStrainElemIDs);

    fprintf('Max abs rate response      = %.3e\n', peakRate(csID));
    fprintf('Max abs disp/acc response  = %.3e\n', peakDisp(csID));
    fprintf('Max abs strain response    = %.3e\n', peakStrain(csID));
    fprintf('Norm of Bp column          = %.3e\n', normBpCol(csID));
end

CS = (1:Nctrl).';
peakTable = table(CS, normBpCol, peakRate, peakDisp, peakStrain, ...
    'VariableNames', {'CS', 'Norm_Bp_Column', 'Peak_Rate', 'Peak_DispAcc', 'Peak_Strain'});

fprintf('\n============================================================\n');
fprintf('Peak response comparison table\n');
fprintf('============================================================\n');
disp(peakTable);

end

function [sys, outputIndex] = build_plant_output_system(cfg, inputs, model)

nStrain = size(inputs.PSI, 1);
nDisp = size(inputs.PHI, 1);
nRate = size(inputs.PHI_ROT, 1);

outputIndex.strain = 1:nStrain;
outputIndex.disp = nStrain + (1:nDisp);
outputIndex.rate = nStrain + nDisp + (1:nRate);

nCS = cfg_value(cfg, 'numControlSurfaces', size(model.Bp, 2));
DaeSelected = select_control_surface_triplets(model.Dae, nCS);

Cp = [model.Cae, DaeSelected];
Dp = zeros(size(Cp, 1), size(model.Bp, 2));

sys = ss(model.Ap, model.Bp, Cp, Dp);

end

function selected = select_control_surface_triplets(matrix, nCS)

if isempty(matrix)
    selected = matrix;
    return
end

if mod(size(matrix, 2), 3) ~= 0
    error('Control-surface columns must be ordered as [delta, delta_dot, delta_ddot].');
end

nAvailable = size(matrix, 2) / 3;
if nAvailable < nCS
    error('Requested %d control surfaces, but only %d are available.', nCS, nAvailable);
end

selected = [ ...
    matrix(:, 1:nCS), ...
    matrix(:, nAvailable+1:nAvailable+nCS), ...
    matrix(:, 2*nAvailable+1:2*nAvailable+nCS)];

end

function plot_output_group(cfg, time, y, titleText, xLabelText, yLabelText, legendPrefix, legendIDs)

if nargin < 8 || isempty(legendIDs)
    legendIDs = 1:size(y, 2);
end

fig = figure('Name', titleText);
apply_plot_style(fig, cfg);

plot(time, y, 'LineWidth', cfg_value(cfg, 'plotLineWidth', 1.5));
grid on; box on;
xlabel(xLabelText);
ylabel(yLabelText);
title(titleText);

legendLabels = arrayfun(@(k) sprintf('%s %d', legendPrefix, k), ...
    legendIDs, 'UniformOutput', false);
legend(legendLabels, 'Location', 'best');
apply_plot_style(fig, cfg);

end

function [model, iV] = select_model_at_speed(cfg, models, requestedSpeed)

[~, iV] = min(abs(cfg.airspeed - requestedSpeed));
model = models{iV};

end

function value = cfg_value(cfg, fieldName, defaultValue)

if isfield(cfg, fieldName)
    value = cfg.(fieldName);
else
    value = defaultValue;
end

end