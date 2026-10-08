function augmented = build_dryden_state_space(cfg, models)
%BUILD_DRYDEN_STATE_SPACE Add Dryden wind and wind-rate states to the plant.
%
% The selected aeroelastic plant must contain one gust vector represented
% by two normalized inputs, ordered as
%
%   u_g = [w_g/V; wdot_g/V].
%
% The Dryden shaping filter converts unit-intensity continuous white noise
% into the physical outputs [w_g; wdot_g].  The returned continuous-time
% model has actuator commands followed by the Dryden white-noise input.

arguments
    cfg struct
    models cell
end

zaeroOutputFile = cfg.ZAERO_filename;
if isfield(cfg,'drydenZaeroOutputFile') && ...
        ~isempty(cfg.drydenZaeroOutputFile)
    zaeroOutputFile = cfg.drydenZaeroOutputFile;
end

drydenData = read_zaero_dryden_psd(zaeroOutputFile);
Vrequested = drydenData.V;

modelSpeeds = cellfun(@(model) model.V, models);
[speedError, modelIndex] = min(abs(modelSpeeds-Vrequested));
model = models{modelIndex};

if speedError > 1e-9*max(1,abs(Vrequested))
    error(['No plant model was built at the ZAERO Dryden velocity %.9g m/s. ', ...
        'The nearest model is at %.9g m/s. Add the Dryden velocity to ', ...
        'cfg.airspeed.'],Vrequested,model.V);
end

if ~strcmpi(model.modelType,'plant') || isempty(model.Ap)
    error('The Dryden augmentation requires cfg.modelType = ''plant''.');
end

if isempty(model.Bpw) || size(model.Bpw,2) ~= 2
    error(['The selected RFA model must contain exactly one gust vector, ', ...
        'giving Bpw = [B_w, B_wdot] with two columns. Use the gust-inclusive ', ...
        '10-root OMCOS APPROX file for the complete control model.']);
end

if size(model.Caw,2) ~= 2
    error('Caw must contain direct-output columns for [w_g/V, wdot_g/V].');
end

tau = 1e-4;
if isfield(cfg,'drydenTau') && ~isempty(cfg.drydenTau)
    tau = cfg.drydenTau;
end

dryden = build_dryden_filter( ...
    model.V,drydenData.Lg,drydenData.sigmaW,tau);

% The RFA gust-input convention is normalized by flight velocity.
CgustNormalized = dryden.C/model.V;
DgustNormalized = dryden.D/model.V;

np = size(model.Ap,1);
nw = size(dryden.A,1);
nControl = size(model.Bp,2);
nOutput = size(model.Cae,1);
noiseScale = dryden.inputNoiseScale;

A = [ ...
    model.Ap,              model.Bpw*CgustNormalized; ...
    zeros(nw,np),          dryden.A];

Bu = [model.Bp; zeros(nw,nControl)];
G = [ ...
    model.Bpw*(noiseScale*DgustNormalized); ...
    noiseScale*dryden.B];

% Sensor outputs.  Dae multiplies the actuator states because the latter
% are ordered as [delta; delta_dot; delta_ddot] in the plant model.
Cp = [model.Cae, model.Dae];
C = [Cp, model.Caw*CgustNormalized];
Du = zeros(nOutput,nControl);
H = model.Caw*(noiseScale*DgustNormalized);

B = [Bu, G];
D = [Du, H];

sys = ss(A,B,C,D);

% Explicit maps are retained so controllers and validation scripts can
% access the physical or normalized gust signals without reconstructing
% state indices.
CgustPhysical = [zeros(2,np), dryden.C];
DgustPhysical = [zeros(2,nControl), noiseScale*dryden.D];
CgustInput = CgustPhysical/model.V;
DgustInput = DgustPhysical/model.V;

augmented = struct();
augmented.A = A;
augmented.B = B;
augmented.Bu = Bu;
augmented.G = G;
augmented.C = C;
augmented.D = D;
augmented.Du = Du;
augmented.H = H;
augmented.sys = sys;
augmented.CgustPhysical = CgustPhysical;
augmented.DgustPhysical = DgustPhysical;
augmented.CgustNormalized = CgustInput;
augmented.DgustNormalized = DgustInput;
augmented.dryden = dryden;
augmented.V = model.V;
augmented.Lg = drydenData.Lg;
augmented.sigmaW = drydenData.sigmaW;
augmented.tau = tau;
augmented.modelIndex = modelIndex;
augmented.zaeroOutputFile = zaeroOutputFile;
augmented.statePartition = struct( ...
    'plant',1:np, ...
    'dryden',np+(1:nw));
augmented.inputPartition = struct( ...
    'control',1:nControl, ...
    'drydenWhiteNoise',nControl+1);

fprintf('\n============================================================\n');
fprintf('Dryden-augmented control state space\n');
fprintf('============================================================\n');
fprintf('ZAERO parameter source     = %s\n',zaeroOutputFile);
fprintf('Selected velocity          = %.6f m/s\n',model.V);
fprintf('Lg                         = %.6f m\n',drydenData.Lg);
fprintf('sigma_w                    = %.6f m/s\n',drydenData.sigmaW);
fprintf('gust-rate LPF tau          = %.6g s\n',tau);
fprintf('plant states               = %d\n',np);
fprintf('Dryden states              = %d\n',nw);
fprintf('total states               = %d\n',size(A,1));
fprintf('control inputs             = %d\n',nControl);
fprintf('stochastic inputs          = %d\n',size(G,2));
fprintf('gust channels into plant   = %d ([w_g/V, wdot_g/V])\n', ...
    size(model.Bpw,2));
fprintf('sensor outputs             = %d\n',nOutput);
fprintf('============================================================\n\n');

end
