function plot_control_bodes(cfg, models)
%PLOT_CONTROL_BODES Plot control-surface commands to modal coordinates.

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
        [Bdelta, Bdeltadot, Bdeltaddot] = select_ae_input_groups(Bae, Nctrl, NctrlAvailable);
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
                     'Plotting only the first %d inputs.'], size(B,2), Nctrl, Nctrl);
            B = B(:, 1:Nctrl);
        end
        inputDescription = 'actuator command inputs';
        modelDescription = 'plant Ap';
        inputLabel = 'actuator command';

    otherwise
        error('Unknown cfg.modelType "%s". Use ''ae'' or ''plant''.', char(cfg.modelType));
end

fprintf('\n============================================================\n');
fprintf('Bode model selected at V = %.2f m/s\n', models{iV}.V);
fprintf('Bode model uses %s\n', modelDescription);
fprintf('Number of %s = %d\n', inputDescription, Nctrl);
fprintf('Number of modal coordinates = %d\n', cfg.Nel);
fprintf('============================================================\n');

for mode_id = 1:cfg.Nel
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
        fig = figure;
        apply_plot_style(fig, cfg);
        plot_bode_response(fig, sys_xi_mode(:, cs_id), cfg.bodeFrequency);

        add_bode_title(sprintf('Bode: %s %d to \\xi_%d, V = %.1f m/s', ...
            inputLabel, cs_id, mode_id, models{iV}.V));
        apply_plot_style(fig, cfg);
    end
end

end

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

function [Bdelta, Bdeltadot, Bdeltaddot] = select_ae_input_groups(Bae, Nctrl, availableCount)

if availableCount > Nctrl
    warning(['Bae contains %d control-surface input groups, but cfg.numControlSurfaces is %d. ', ...
             'Plotting only the first %d groups.'], availableCount, Nctrl, Nctrl);
end

Bdelta = Bae(:, 1:Nctrl);
Bdeltadot = Bae(:, availableCount+1:availableCount+Nctrl);
Bdeltaddot = Bae(:, 2*availableCount+1:2*availableCount+Nctrl);

end

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
