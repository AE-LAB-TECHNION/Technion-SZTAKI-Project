function plot_control_bodes(cfg, models)
%PLOT_CONTROL_BODES Plot control-surface commands to modal coordinates.

if cfg.isinp ~= 1
    fprintf('\nBode plots skipped because isinp = 0.\n');
    fprintf('No control-surface command inputs are included, so Ap/Bp are not built.\n');
    return
end

[~, iV] = min(abs(cfg.airspeed - cfg.bodeSpeed));

Ap = models{iV}.Ap;
Bp = models{iV}.Bp;

Nctrl = size(Bp,2);     % number of control-surface command inputs

fprintf('\n============================================================\n');
fprintf('Bode model selected at V = %.2f m/s\n', models{iV}.V);
fprintf('Bode model uses Ap because isinp = 1\n');
fprintf('Number of control-surface command inputs = %d\n', Nctrl);
fprintf('Number of modal coordinates = %d\n', cfg.Nel);
fprintf('============================================================\n');

for mode_id = 1:cfg.Nel
    Cmodal = zeros(1, size(Ap,1));
    Cmodal(1, mode_id) = 1;

    sys_xi_mode = ss(Ap, Bp, Cmodal, zeros(1, Nctrl));

    for cs_id = 1:Nctrl
        figure;
        bode(sys_xi_mode(:, cs_id), cfg.bodeFrequency);
        grid on;

        title(sprintf('Bode: CS %d command to \\xi_%d, V = %.1f m/s', ...
            cs_id, mode_id, models{iV}.V));
    end
end

end
