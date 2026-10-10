function [model, Aplot, plotModelName] = build_model_at_speed(cfg, inputs, V)
%BUILD_MODEL_AT_SPEED Build aeroelastic and optional plant models at one speed.

includeControlInputs = true;
modelType = lower(char(cfg.modelType));

[Aae, Bae, Baw, Cae, Caw, Dae] = buildAESS_state( ...
    cfg.Nel, inputs.eigvals, cfg.zeta, ...
    inputs.RFA_mat, cfg.method, ...
    cfg.L, V, cfg.rho, ...
    includeControlInputs, cfg.isgust, ...
    inputs.PSI, inputs.PHI, inputs.PHI_ROT, ...
    inputs.RBM.modalCoefficients, ...
    cfg.numControlSurfaces, cfg.numGustInputs);

switch modelType
    case 'ae'
        Ap = [];
        Bp = [];
        Bpw = [];
        Aplot = Aae;
        plotModelName = 'Aae';

    case 'plant'
        [Ap, Bp, Bpw] = buildPlant_from_AESS(Aae, Bae, Baw, cfg.numControlSurfaces);
        Aplot = Ap;
        plotModelName = 'Ap';

    otherwise
        error('Unknown cfg.modelType "%s". Use ''ae'' or ''plant''.', char(cfg.modelType));
end

model = struct();
model.V = V;
model.modelType = modelType;

model.Aae = Aae;
model.Bae = Bae;
model.Baw = Baw;
model.Cae = Cae;
model.Caw = Caw;
model.Dae = Dae;

nStrain = size(inputs.PSI,1);
nAcceleration = size(inputs.PHI,1);
nRate = size(inputs.PHI_ROT,1);
nRbm = size(inputs.RBM.modalCoefficients,1);

model.outputIndex = struct( ...
    'strain',1:nStrain, ...
    'acceleration',nStrain+(1:nAcceleration), ...
    'rate',nStrain+nAcceleration+(1:nRate), ...
    'rbm',nStrain+nAcceleration+nRate+(1:nRbm));

model.Ap = Ap;
model.Bp = Bp;
model.Bpw = Bpw;

if strcmp(modelType,'plant')
    model.Cp = [Cae, Dae];
else
    model.Cp = Cae;
end

% Retained for compatibility with older workspace expectations.
model.Aact = [];
model.Bact = [];

end
