function [models, eigTraj, plotModelName] = run_airspeed_sweep(cfg, inputs)
%RUN_AIRSPEED_SWEEP Build models and eigenvalue trajectories over airspeed.

nV = numel(cfg.airspeed);
if nV == 0
    error('cfg.airspeed must contain at least one velocity.');
end

models = cell(nV, 1);

[model, Aplot, plotModelName] = build_model_at_speed( ...
    cfg,inputs,cfg.airspeed(1));
models{1} = model;
eigTraj = zeros(size(Aplot,1),nV);
eigTraj(:,1) = eig(Aplot);
print_model_sizes(model,plotModelName);

for i = 2:nV
    V = cfg.airspeed(i);

    [model, Aplot, plotModelName] = build_model_at_speed(cfg, inputs, V);

    models{i} = model;
    eigTraj(:, i) = eig(Aplot);
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
