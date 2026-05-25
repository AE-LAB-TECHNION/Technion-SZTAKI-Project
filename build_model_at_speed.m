function [model, Aplot, plotModelName] = build_model_at_speed(cfg, inputs, V)
%BUILD_MODEL_AT_SPEED Build aeroelastic and optional plant models at one speed.

[Aae, Bae, Baw, Cae, Caw, Dae] = buildAESS_state( ...
    cfg.Nel, inputs.eigvals, cfg.zeta, ...
    cfg.RFA_filename, cfg.method, ...
    cfg.L, V, cfg.rho, ...
    cfg.isinp, cfg.isgust, ...
    inputs.PSI, inputs.PHI, inputs.PHI_ROT);

if cfg.isinp == 1
    [Ap, Bp, Bpw] = buildPlant_from_AESS(Aae, Bae, Baw);
    Aplot = Ap;
    plotModelName = 'Ap';
else
    Ap = [];
    Bp = [];
    Bpw = [];
    Aplot = Aae;
    plotModelName = 'Aae';
end

model = struct();
model.V = V;

model.Aae = Aae;
model.Bae = Bae;
model.Baw = Baw;
model.Cae = Cae;
model.Caw = Caw;
model.Dae = Dae;

model.Ap = Ap;
model.Bp = Bp;
model.Bpw = Bpw;

% Retained for compatibility with older workspace expectations.
model.Cp = [];
model.Aact = [];
model.Bact = [];

end
