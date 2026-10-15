function comparison = compare_omcos_solvers(cfg, model, time, gustInput, zaeroTime, zaero)
%COMPARE_OMCOS_SOLVERS Compare all MATLAB OMCOS solvers against ZAERO.

Nel = cfg.Nel;
nStates = size(model.Aae, 1);
nInputs = size(model.Baw, 2);

Cmodal = zeros(2*Nel, nStates);
Cmodal(1:Nel, 1:Nel) = eye(Nel);
Cmodal(Nel+(1:Nel), Nel+(1:Nel)) = eye(Nel);

sys = ss(model.Aae, model.Baw, Cmodal, ...
    zeros(2*Nel, nInputs));
yLsim = lsim(sys, gustInput, time);

simulations = struct();
simulations.lsim.time = time;
simulations.lsim.modalX = yLsim(:, 1:Nel);
simulations.lsim.modalV = yLsim(:, Nel+(1:Nel));
simulations.ode45 = simulate_omcos_ode45(cfg, model, time);
simulations.ode113 = simulate_omcos_ode113(cfg, model, time);
simulations.expm = simulate_omcos_expm(cfg, model, time);

solverNames = {'lsim', 'ode45', 'ode113', 'expm'};
for k = 1:numel(solverNames)
    name = solverNames{k};
    sampledX = interp1(time, simulations.(name).modalX, ...
        zaeroTime, 'linear', 'extrap');
    sampledV = interp1(time, simulations.(name).modalV, ...
        zaeroTime, 'linear', 'extrap');

    comparison.(name) = simulations.(name);
    comparison.(name).modalX = sampledX;
    comparison.(name).modalV = sampledV;
    comparison.(name).errorX = sampledX - zaero.modalX;
    comparison.(name).errorV = sampledV - zaero.modalV;
end

plot_solver_errors(cfg, comparison, zaeroTime, Nel);

end

function plot_solver_errors(cfg, comparison, time, Nel)

solverNames = {'lsim', 'ode45', 'ode113', 'expm'};
solverLabels = {'lsim', 'ode45', 'ode113', 'expm'};
lineStyles = {'-', '--', '-.', ':'};
colors = [ ...
    0.00 0.00 0.00; ...
    0.85 0.20 0.10; ...
    0.10 0.35 0.80; ...
    0.10 0.55 0.25];

nPlotModes = min(3, Nel);
timeLimit = min(5, time(end));
fig = figure('Name', 'OMCOS solver errors versus ZAERO');
apply_plot_style(fig, cfg);
tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

for modeIndex = 1:nPlotModes
    nexttile;
    hold on;
    for solverIndex = 1:numel(solverNames)
        name = solverNames{solverIndex};
        plot(time, comparison.(name).errorX(:, modeIndex), ...
            'Color', colors(solverIndex, :), ...
            'LineStyle', lineStyles{solverIndex}, ...
            'LineWidth', 1.3);
    end
    yline(0, ':', 'Color', [0.45 0.45 0.45]);
    grid on; box on;
    xlim([0, timeLimit]);
    title(sprintf('MODALX error, mode %d', modeIndex));
    xlabel('Time [s]');
    ylabel('\Delta\xi');
    if modeIndex == 1
        legend(solverLabels, 'Location', 'best');
    end
end

for modeIndex = 1:nPlotModes
    nexttile;
    hold on;
    for solverIndex = 1:numel(solverNames)
        name = solverNames{solverIndex};
        plot(time, comparison.(name).errorV(:, modeIndex), ...
            'Color', colors(solverIndex, :), ...
            'LineStyle', lineStyles{solverIndex}, ...
            'LineWidth', 1.3);
    end
    yline(0, ':', 'Color', [0.45 0.45 0.45]);
    grid on; box on;
    xlim([0, timeLimit]);
    title(sprintf('MODALV error, mode %d', modeIndex));
    xlabel('Time [s]');
    ylabel('\Delta(d\xi/dt)');
end
apply_plot_style(fig, cfg);

end
