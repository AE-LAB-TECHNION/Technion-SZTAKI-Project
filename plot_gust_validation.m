function result = plot_gust_validation(cfg, inputs, models)
%PLOT_GUST_VALIDATION Compare ZAERO GLOADS OMCOS gust response with MATLAB model.
%
% This function compares ZAERO MLDPRNT output with the MATLAB aeroelastic model.
%
% ZAERO expected MLDPRNT order in GUSTXI.DAT:
%   TIME, MODALV(1:Nel), MODALX(1:Nel)
%
% MATLAB model:
%   Uses Aae and Baw.
%   Baw expects normalized gust inputs:
%       u = [w_g/V ; d(w_g/V)/dt]
%
% DGUST OMCOSP interpretation used here:
%   DGUST  SID  OMCOSP  LG  WGV  X0  DIR  DT
%
%   LG  = one-minus-cosine gust duration [s] if lengthMode = 'time'
%   WGV = normalized gust amplitude w_g/V

result = struct();

% -------------------------------------------------------------------------
% Select model at validation speed
% -------------------------------------------------------------------------
[model, iV] = select_model_at_speed(cfg, models, cfg_value(cfg, 'gustValidationSpeed', 35));
Nel = cfg.Nel;

if isempty(model.Aae) || isempty(model.Baw)
    fprintf('\nGust validation skipped because the selected model has no Baw gust input matrix.\n');
    return
end

% -------------------------------------------------------------------------
% Read ZAERO response
% -------------------------------------------------------------------------
zaeroFile = cfg.gustResponseFile;
zaero = read_zaero_gust_modal_response(zaeroFile, Nel);

tCompare = zaero.time;

% Use a finer internal MATLAB time step, then sample MATLAB at ZAERO times
internalDt = cfg_value(cfg, 'gustInternalDt', 0.001);
tCode = make_fine_time(tCompare(1), tCompare(end), internalDt);

% -------------------------------------------------------------------------
% Build OMCOSP gust input
% -------------------------------------------------------------------------
[ugustCode, gust] = build_omcosp_input(cfg, tCode, model.V, size(model.Baw, 2));

% -------------------------------------------------------------------------
% Build modal output system: y = [xi ; xidot]
% -------------------------------------------------------------------------
nStates = size(model.Aae, 1);

Cmodal = zeros(2*Nel, nStates);
Cmodal(1:Nel, 1:Nel) = eye(Nel);                 % xi
Cmodal(Nel+(1:Nel), Nel+(1:Nel)) = eye(Nel);     % xidot

Dmodal = zeros(2*Nel, size(model.Baw, 2));

sys = ss(model.Aae, model.Baw, Cmodal, Dmodal);

% -------------------------------------------------------------------------
% Simulate MATLAB response on fine time grid
% -------------------------------------------------------------------------
yCode = lsim(sys, ugustCode, tCode);

code.time = tCode;
code.modalX = yCode(:, 1:Nel);
code.modalV = yCode(:, Nel+(1:Nel));

% Sample MATLAB response at ZAERO output times
codeCompare.time = tCompare;
codeCompare.modalX = interp1(code.time, code.modalX, tCompare, 'linear', 'extrap');
codeCompare.modalV = interp1(code.time, code.modalV, tCompare, 'linear', 'extrap');

% -------------------------------------------------------------------------
% Errors
% -------------------------------------------------------------------------
errX = codeCompare.modalX - zaero.modalX;
errV = codeCompare.modalV - zaero.modalV;

result.time = tCompare;
result.gust = gust;

result.matlabFine.time = code.time;
result.matlabFine.modalX = code.modalX;
result.matlabFine.modalV = code.modalV;

result.matlabSampled.time = codeCompare.time;
result.matlabSampled.modalX = codeCompare.modalX;
result.matlabSampled.modalV = codeCompare.modalV;

result.zaero.modalX = zaero.modalX;
result.zaero.modalV = zaero.modalV;

result.error.modalX = errX;
result.error.modalV = errV;

result.rms.modalX = sqrt(mean(errX.^2, 1));
result.rms.modalV = sqrt(mean(errV.^2, 1));

result.peakAbs.modalX = max(abs(errX), [], 1);
result.peakAbs.modalV = max(abs(errV), [], 1);

% -------------------------------------------------------------------------
% Print summary
% -------------------------------------------------------------------------
fprintf('\n============================================================\n');
fprintf('Gust validation: ZAERO MLDPRNT vs MATLAB OMCOS\n');
fprintf('============================================================\n');
fprintf('Selected V = %.2f m/s, index %d\n', model.V, iV);
fprintf('ZAERO response file: %s\n', zaeroFile);
fprintf('Aae size: %d x %d\n', size(model.Aae, 1), size(model.Aae, 2));
fprintf('Baw size: %d x %d\n', size(model.Baw, 1), size(model.Baw, 2));
fprintf('DGUST interpretation: LG = %.6g s, WGV = %.6g, ramp time = %.6g s\n', ...
    gust.LG, gust.WGV, gust.Tg);

Mode = (1:Nel).';

rmsModalX = result.rms.modalX(:);
peakModalX = result.peakAbs.modalX(:);
peakZaeroX = max(abs(zaero.modalX), [], 1).';
peakMatlabX = max(abs(codeCompare.modalX), [], 1).';

rmsModalV = result.rms.modalV(:);
peakModalV = result.peakAbs.modalV(:);
peakZaeroV = max(abs(zaero.modalV), [], 1).';
peakMatlabV = max(abs(codeCompare.modalV), [], 1).';

summaryTable = table( ...
    Mode, ...
    rmsModalX, peakModalX, peakZaeroX, peakMatlabX, ...
    rmsModalV, peakModalV, peakZaeroV, peakMatlabV, ...
    'VariableNames', { ...
        'Mode', ...
        'RMS_X', 'PeakErr_X', 'PeakZAERO_X', 'PeakMATLAB_X', ...
        'RMS_V', 'PeakErr_V', 'PeakZAERO_V', 'PeakMATLAB_V'});

disp(summaryTable);
result.summaryTable = summaryTable;

% -------------------------------------------------------------------------
% Plots
% -------------------------------------------------------------------------
plot_gust_input(cfg, code.time, ugustCode, model.V, gust);
plot_modal_overlay(cfg, code, codeCompare, zaero, Nel);

end

% =========================================================================
function t = make_fine_time(t0, tf, dt)

t = (t0:dt:tf).';

if t(end) < tf
    t(end+1, 1) = tf;
end

end

% =========================================================================
function response = read_zaero_gust_modal_response(filename, Nel)

fid = fopen(filename, 'r');

if fid < 0
    error('Could not open ZAERO gust response file: %s', filename);
end

cleanup = onCleanup(@() fclose(fid));

rows = zeros(0, 1 + 2*Nel);

while true
    line = fgetl(fid);

    if ~ischar(line)
        break
    end

    values = sscanf(line, '%f').';

    if numel(values) >= 1 + 2*Nel
        rows(end+1, :) = values(1:1+2*Nel); %#ok<AGROW>
    end
end

if isempty(rows)
    error('No numeric modal time-response rows were found in %s.', filename);
end

% Actual ZAERO GUSTXI.DAT order observed:
% TIME, MODALV(1:Nel), MODALX(1:Nel)
response.time   = rows(:, 1);
response.modalV = rows(:, 1 + (1:Nel));
response.modalX = rows(:, 1 + Nel + (1:Nel));

end

% =========================================================================
function [u, gust] = build_omcosp_input(cfg, time, V, nInputs)

if nInputs < 2 || mod(nInputs, 2) ~= 0
    error('Expected Baw inputs ordered as [w_g/V, wdot_g/V] for each gust vector, but Baw has %d columns.', nInputs);
end

Ngust = nInputs / 2;

% Match ZAERO DGUST OMCOSP:
% DGUST SID OMCOSP LG WGV X0 DIR DT
gust.LG = cfg_value(cfg, 'gustLG', 1.0);                 % [s], if lengthMode='time'
gust.WGV = cfg_value(cfg, 'gustWGV', 0.1);               % normalized amplitude w_g/V
gust.X0 = cfg_value(cfg, 'gustX0', 0.0);
gust.xResponse = cfg_value(cfg, 'gustXResponse', 0.0);
gust.lengthMode = cfg_value(cfg, 'gustLengthMode', 'time');

eta = time(:) - (gust.xResponse - gust.X0)/V;

switch lower(gust.lengthMode)
    case 'time'
        Tg = gust.LG;
    case 'length'
        Tg = gust.LG / V;
    otherwise
        error('cfg.gustLengthMode must be ''time'' or ''length''.');
end

if Tg <= 0
    error('The one-minus-cosine gust duration must be positive.');
end

wbar = zeros(numel(time), 1);      % w_g/V
wdotbar = zeros(numel(time), 1);   % d(w_g/V)/dt = wdot_g/V

active = eta >= 0 & eta <= Tg;

wbar(active) = 0.5*gust.WGV*(1 - cos(2*pi*eta(active)/Tg));
wdotbar(active) = (pi*gust.WGV/Tg)*sin(2*pi*eta(active)/Tg);

% Retain physical values for reporting, while Baw receives normalized inputs.
wg = V*wbar;
wdot = V*wdotbar;

u = zeros(numel(time), nInputs);

for ig = 1:Ngust
    u(:, ig) = wbar;
    u(:, Ngust + ig) = wdotbar;
end

gust.Tg = Tg;
gust.wbar = wbar;
gust.wdotbar = wdotbar;
gust.wg = wg;
gust.wdot = wdot;

end

% =========================================================================
function plot_modal_overlay(cfg, code, codeCompare, zaero, Nel)

% -------------------------------------------------------------------------
% Modal displacement
% -------------------------------------------------------------------------
fig = figure('Name', 'Gust validation: modal displacement');
apply_plot_style(fig, cfg);

tiledlayout(Nel, 2, 'TileSpacing', 'compact');

for k = 1:Nel

    nexttile(2*k - 1);

    markerIdx = marker_indices(numel(zaero.time));

    plot(code.time, code.modalX(:, k), '-', ...
        'Color', [0.75 0.75 0.75], ...
        'LineWidth', 0.8);
    hold on;

    plot(zaero.time, zaero.modalX(:, k), 'k-', ...
        'LineWidth', cfg_value(cfg, 'plotLineWidth', 1.5));

    plot(codeCompare.time, codeCompare.modalX(:, k), 'r--', ...
        'LineWidth', cfg_value(cfg, 'plotLineWidth', 1.5));

    plot(codeCompare.time(markerIdx), codeCompare.modalX(markerIdx, k), 'ro', ...
        'MarkerSize', 3.0, ...
        'MarkerFaceColor', 'w');

    grid on; box on;
    ylabel(sprintf('\\xi_%d', k));

    if k == 1
        legend({'MATLAB fine', 'ZAERO', 'MATLAB sampled', 'Sample points'}, ...
            'Location', 'best');
    end

    nexttile(2*k);

    errX = codeCompare.modalX(:, k) - zaero.modalX(:, k);

    plot(codeCompare.time, errX, 'k-', ...
        'LineWidth', 1.1);
    hold on;

    yline(0, ':', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.0);

    grid on; box on;
    ylabel(sprintf('\\Delta \\xi_%d', k));

end

xlabel('Time [s]');
apply_plot_style(fig, cfg);

% -------------------------------------------------------------------------
% Modal velocity
% -------------------------------------------------------------------------
fig = figure('Name', 'Gust validation: modal velocity');
apply_plot_style(fig, cfg);

tiledlayout(Nel, 2, 'TileSpacing', 'compact');

for k = 1:Nel

    nexttile(2*k - 1);

    markerIdx = marker_indices(numel(zaero.time));

    plot(code.time, code.modalV(:, k), '-', ...
        'Color', [0.75 0.75 0.75], ...
        'LineWidth', 0.8);
    hold on;

    plot(zaero.time, zaero.modalV(:, k), 'k-', ...
        'LineWidth', cfg_value(cfg, 'plotLineWidth', 1.5));

    plot(codeCompare.time, codeCompare.modalV(:, k), 'r--', ...
        'LineWidth', cfg_value(cfg, 'plotLineWidth', 1.5));

    plot(codeCompare.time(markerIdx), codeCompare.modalV(markerIdx, k), 'ro', ...
        'MarkerSize', 3.0, ...
        'MarkerFaceColor', 'w');

    grid on; box on;
    ylabel(sprintf('d\\xi_%d/dt', k));

    if k == 1
        legend({'MATLAB fine', 'ZAERO', 'MATLAB sampled', 'Sample points'}, ...
            'Location', 'best');
    end

    nexttile(2*k);

    errV = codeCompare.modalV(:, k) - zaero.modalV(:, k);

    plot(codeCompare.time, errV, 'k-', ...
        'LineWidth', 1.1);
    hold on;

    yline(0, ':', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.0);

    grid on; box on;
    ylabel(sprintf('\\Delta(d\\xi_%d/dt)', k));

end

xlabel('Time [s]');
apply_plot_style(fig, cfg);

end

% =========================================================================
function plot_gust_input(cfg, time, u, V, gust)

fig = figure('Name', 'Gust validation: OMCOS input');
apply_plot_style(fig, cfg);

plot(time, gust.wbar, 'LineWidth', cfg_value(cfg, 'plotLineWidth', 1.5));
hold on;
plot(time, gust.wdotbar, '--', 'LineWidth', cfg_value(cfg, 'plotLineWidth', 1.5));

grid on; box on;
xlabel('Time [s]');
ylabel('Normalized gust input');
legend({'w_g/V', 'd(w_g/V)/dt'}, 'Location', 'best');

% Optional zoom for the gust input
if cfg_value(cfg, 'gustZoomInput', true)
    xlim([0, min(time(end), 1.2*gust.Tg)]);
end

apply_plot_style(fig, cfg);

end

% =========================================================================
function idx = marker_indices(n)

step = max(1, floor(n/80));
idx = 1:step:n;

if idx(end) ~= n
    idx(end+1) = n;
end

end

% =========================================================================
function [model, iV] = select_model_at_speed(cfg, models, requestedSpeed)

[~, iV] = min(abs(cfg.airspeed - requestedSpeed));
model = models{iV};

end

% =========================================================================
function value = cfg_value(cfg, fieldName, defaultValue)

if isfield(cfg, fieldName) && ~isempty(cfg.(fieldName))
    value = cfg.(fieldName);
else
    value = defaultValue;
end

end
