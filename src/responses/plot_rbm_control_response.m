function response = plot_rbm_control_response(cfg, models)
%PLOT_RBM_CONTROL_RESPONSE Plot actuator-command influence on the RBM output.

response = struct();

modelSpeeds = cellfun(@(model) model.V, models);
[~, modelIndex] = min(abs(modelSpeeds - cfg.bodeSpeed));
model = models{modelIndex};

if ~strcmpi(model.modelType, 'plant') || isempty(model.Ap) || isempty(model.Bp)
    fprintf('\nRBM control-response plot skipped because it requires the Plant model.\n');
    return
end

if ~isfield(model, 'outputIndex') || ~isfield(model.outputIndex, 'rbm')
    error('The selected model does not define an RBM output row.');
end

rbmRow = model.outputIndex.rbm;
if numel(rbmRow) ~= 1
    error('The RBM control-response plot requires exactly one RBM output.');
end

nControls = min(cfg.numControlSurfaces, size(model.Bp, 2));
if nControls < 1
    fprintf('\nRBM control-response plot skipped because no actuator commands are available.\n');
    return
end

Bcontrol = model.Bp(:, 1:nControls);
Crbm = model.Cp(rbmRow, :);
Drbm = zeros(1, nControls);
sysControlToRbm = ss(model.Ap, Bcontrol, Crbm, Drbm);

frequencyRad = cfg.bodeFrequency(:);
frequencyHz = frequencyRad / (2*pi);
frequencyResponse = freqresp(sysControlToRbm, frequencyRad);
magnitude = permute(abs(frequencyResponse), [3, 2, 1]);
magnitude = reshape(magnitude, numel(frequencyRad), nControls);
magnitudeDb = 20*log10(max(magnitude, realmin));

bodeFig = figure('Name', 'Control commands to RBM - Bode magnitude');
apply_plot_style(bodeFig, cfg);

semilogx(frequencyHz, magnitudeDb, 'LineWidth', cfg.plotLineWidth);
grid on;
box on;
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
xlim([frequencyHz(1), frequencyHz(end)]);

legendLabels = arrayfun(@(controlIndex) sprintf('CS %d command', controlIndex), ...
    1:nControls, 'UniformOutput', false);
legend(legendLabels, 'Location', 'best');

apply_plot_style(bodeFig, cfg);

[peakMagnitudeDb, peakIndices] = max(magnitudeDb, [], 1);
peakFrequencyHz = frequencyHz(peakIndices);
controlSurface = (1:nControls).';
bodePeakTable = table(controlSurface, peakMagnitudeDb.', peakFrequencyHz, ...
    'VariableNames', {'ControlSurface', 'PeakMagnitude_dB', 'PeakFrequency_Hz'});

fprintf('\n============================================================\n');
fprintf('Actuator-command influence on RBM at V = %.2f m/s\n', model.V);
fprintf('============================================================\n');
disp(bodePeakTable);

response.modelIndex = modelIndex;
response.speed = model.V;
response.frequencyHz = frequencyHz;
response.magnitude = magnitude;
response.magnitudeDb = magnitudeDb;
response.bodePeakTable = bodePeakTable;
response.sys = sysControlToRbm;
response.bodeFigure = bodeFig;

end
