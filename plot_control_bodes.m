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

        Nctrl = size(Bae,2) / 3;
        Bdelta = Bae(:, 1:Nctrl);
        Bdeltadot = Bae(:, Nctrl+1:2*Nctrl);
        Bdeltaddot = Bae(:, 2*Nctrl+1:3*Nctrl);
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

        Nctrl = size(B,2);
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
        figure;
        bode(sys_xi_mode(:, cs_id), cfg.bodeFrequency);
        grid on;

        title(sprintf('Bode: %s %d to \\xi_%d, V = %.1f m/s', ...
            inputLabel, cs_id, mode_id, models{iV}.V));
    end
end

end
