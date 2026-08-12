function peakTable = plot_control_bodes(cfg, models)
%PLOT_CONTROL_BODES Plot control-surface commands to modal coordinates
% and print the Bode magnitude peak for each transfer function.

peakTable = table();

[~, iV] = min(abs(cfg.airspeed - cfg.bodeSpeed));
modelType = lower(char(cfg.modelType));

switch modelType
    case 'ae'
        A = models{iV}.Aae;
        Bae = models{iV}.Bae;

        if isempty(Bae) || size(Bae,2) == 0
            fprintf('\nBode plots skipped because the AE model has no control-surface inputs.\n');
            return
        end

        if mod(size(Bae,2), 3) ~= 0
            error('Bae columns must be ordered as [delta; delta_dot; delta_ddot].');
        end

        NctrlAvailable = size(Bae,2) / 3;
        Nctrl = control_surface_count_for_plot(cfg, NctrlAvailable);

        [Bdelta, Bdeltadot, Bdeltaddot] = select_ae_input_groups( ...
            Bae, Nctrl, NctrlAvailable);

        inputDescription = 'CS deflection inputs';
        modelDescription = 'aeroelastic Aae';
        inputLabel = 'CS deflection';

    case 'plant'
        A = models{iV}.Ap;
        B = models{iV}.Bp;

        if isempty(B) || size(B,2) == 0
            fprintf('\nBode plots skipped because the plant model has no actuator command inputs.\n');
            return
        end

        Nctrl = control_surface_count_for_plot(cfg, size(B,2));

        if size(B,2) > Nctrl
            warning(['Plant B matrix has %d actuator command inputs, but cfg.numControlSurfaces is %d. ', ...
                     'Plotting only the first %d inputs.'], ...
                     size(B,2), Nctrl, Nctrl);
            B = B(:, 1:Nctrl);
        end

        inputDescription = 'actuator command inputs';
        modelDescription = 'plant Ap';
        inputLabel = 'actuator command';

    otherwise
        error('Unknown cfg.modelType "%s". Use ''ae'' or ''plant''.', char(cfg.modelType));
end

% Use physical FEM mode labels if available: for example [1 2 4]
if isfield(cfg, 'selected_modes') && numel(cfg.selected_modes) >= cfg.Nel
    modeLabels = cfg.selected_modes(:).';
else
    modeLabels = 1:cfg.Nel;
end

fprintf('\n============================================================\n');
fprintf('Bode model selected at V = %.2f m/s\n', models{iV}.V);
fprintf('Bode model uses %s\n', modelDescription);
fprintf('Number of %s = %d\n', inputDescription, Nctrl);
fprintf('Number of modal coordinates = %d\n', cfg.Nel);
fprintf('============================================================\n');

% Storage for peak table
rowModeIndex = [];
rowFEMMode = [];
rowCS = [];
rowPeakGainDb = [];
rowPeakFreqRad = [];
rowPeakFreqHz = [];

for mode_id = 1:cfg.Nel

    physicalMode = modeLabels(mode_id);

    % Output: selected modal coordinate xi_mode
    Cmodal = zeros(1, size(A,1));
    Cmodal(1, mode_id) = 1;

    if strcmpi(modelType, 'ae')
        s = tf('s');

        sys_delta = ss(A, Bdelta, Cmodal, zeros(1, Nctrl));
        sys_deltadot = ss(A, Bdeltadot, Cmodal, zeros(1, Nctrl));
        sys_deltaddot = ss(A, Bdeltaddot, Cmodal, zeros(1, Nctrl));

        sys_xi_mode = sys_delta + s*sys_deltadot + s^2*sys_deltaddot;

    else
        sys_xi_mode = ss(A, B, Cmodal, zeros(1, Nctrl));
    end

    for cs_id = 1:Nctrl

        currentSys = sys_xi_mode(:, cs_id);

        fig = figure;
        apply_plot_style(fig, cfg);

        plot_bode_response(fig, currentSys, cfg.bodeFrequency);

        [peakGainDb, peakFreqRad, peakFreqHz] = compute_bode_peak( ...
            currentSys, cfg.bodeFrequency);

        fprintf('FEM Mode %d, CS %d: peak = %.2f dB at %.2f rad/s = %.2f Hz\n', ...
            physicalMode, cs_id, peakGainDb, peakFreqRad, peakFreqHz);

        add_bode_title(sprintf('Bode: %s %d to \\xi_%d, V = %.1f m/s', ...
            inputLabel, cs_id, physicalMode, models{iV}.V));

        apply_plot_style(fig, cfg);

        % Save values for table
        rowModeIndex(end+1,1) = mode_id; %#ok<AGROW>
        rowFEMMode(end+1,1) = physicalMode; %#ok<AGROW>
        rowCS(end+1,1) = cs_id; %#ok<AGROW>
        rowPeakGainDb(end+1,1) = peakGainDb; %#ok<AGROW>
        rowPeakFreqRad(end+1,1) = peakFreqRad; %#ok<AGROW>
        rowPeakFreqHz(end+1,1) = peakFreqHz; %#ok<AGROW>
    end
end

peakTable = table(rowModeIndex, rowFEMMode, rowCS, ...
    rowPeakGainDb, rowPeakFreqRad, rowPeakFreqHz, ...
    'VariableNames', {'RetainedModeIndex', 'FEM_Mode', 'CS', ...
                      'PeakGain_dB', 'PeakFreq_rad_s', 'PeakFreq_Hz'});

fprintf('\n============================================================\n');
fprintf('Bode peak summary\n');
fprintf('============================================================\n');
disp(peakTable);

end

% =========================================================================
function Nctrl = control_surface_count_for_plot(cfg, availableCount)

if isfield(cfg, 'numControlSurfaces') && ~isempty(cfg.numControlSurfaces)
    Nctrl = cfg.numControlSurfaces;
else
    Nctrl = availableCount;
end

if availableCount < Nctrl
    error('Requested %d control-surface inputs, but only %d are available.', ...
        Nctrl, availableCount);
end

end

% =========================================================================
function [Bdelta, Bdeltadot, Bdeltaddot] = select_ae_input_groups(Bae, Nctrl, availableCount)

if availableCount > Nctrl
    warning(['Bae contains %d control-surface input groups, but cfg.numControlSurfaces is %d. ', ...
             'Plotting only the first %d groups.'], ...
             availableCount, Nctrl, Nctrl);
end

Bdelta = Bae(:, 1:Nctrl);
Bdeltadot = Bae(:, availableCount+1:availableCount+Nctrl);
Bdeltaddot = Bae(:, 2*availableCount+1:2*availableCount+Nctrl);

end

% =========================================================================
function plot_bode_response(fig, sys, frequency)

[mag, phase, wout] = bode(sys, frequency);

magDb = 20*log10(squeeze(mag));
phaseDeg = squeeze(phase);
wout = squeeze(wout);

magDb = magDb(:);
phaseDeg = phaseDeg(:);
wout = wout(:);

figure(fig);

axMag = subplot(2, 1, 1);
semilogx(axMag, wout, magDb);
ylabel(axMag, 'Magnitude [dB]');
grid(axMag, 'on');
box(axMag, 'on');

axPhase = subplot(2, 1, 2);
semilogx(axPhase, wout, phaseDeg);
xlabel(axPhase, 'Frequency [rad/s]');
ylabel(axPhase, 'Phase [deg]');
grid(axPhase, 'on');
box(axPhase, 'on');

linkaxes([axMag, axPhase], 'x');

end

% =========================================================================
function [peakGainDb, peakFreqRad, peakFreqHz] = compute_bode_peak(sys, frequency)

[mag, ~, wout] = bode(sys, frequency);

magDb = 20*log10(squeeze(mag));
wout = squeeze(wout);

magDb = magDb(:);
wout = wout(:);

[peakGainDb, idxPeak] = max(magDb);

peakFreqRad = wout(idxPeak);
peakFreqHz = peakFreqRad / (2*pi);

end

% =========================================================================
function add_bode_title(titleText)

try
    titleHandle = sgtitle(titleText);
    set(titleHandle, 'FontWeight', 'normal');
catch
    axesHandles = findall(gcf, 'Type', 'axes');
    if ~isempty(axesHandles)
        title(axesHandles(end), titleText, 'FontWeight', 'normal');
    end
end

end