function [models, eigTraj, plotModelName] = run_airspeed_sweep(cfg, inputs)
%RUN_AIRSPEED_SWEEP Build models and eigenvalue trajectories over airspeed.

nV = numel(cfg.airspeed);
models = cell(nV, 1);
eigTraj = [];
plotModelName = '';

for i = 1:nV
    V = cfg.airspeed(i);

    [model, Aplot, plotModelName] = build_model_at_speed(cfg, inputs, V);

    models{i} = model;
    eigTraj(:, i) = eig(Aplot);

    if i == 1
        print_model_sizes(model, plotModelName);
    end
end

end

function print_model_sizes(model, plotModelName)

fprintf('\nMatrix sizes at V = %.2f m/s:\n', model.V);
fprintf('  Aae : %d x %d\n', size(model.Aae,1), size(model.Aae,2));
fprintf('  Bae : %d x %d\n', size(model.Bae,1), size(model.Bae,2));
fprintf('  Baw : %d x %d\n', size(model.Baw,1), size(model.Baw,2));
fprintf('  Cae : %d x %d\n', size(model.Cae,1), size(model.Cae,2));
fprintf('  Caw : %d x %d\n', size(model.Caw,1), size(model.Caw,2));
fprintf('  Dae : %d x %d\n', size(model.Dae,1), size(model.Dae,2));

if strcmpi(model.modelType, 'plant')
    fprintf('  Ap  : %d x %d\n', size(model.Ap,1), size(model.Ap,2));
    fprintf('  Bp  : %d x %d\n', size(model.Bp,1), size(model.Bp,2));
    fprintf('  Bpw : %d x %d\n', size(model.Bpw,1), size(model.Bpw,2));
end

fprintf('\nRoot-locus plot uses: %s\n', plotModelName);

end
