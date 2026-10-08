function response = compute_dryden_frequency_response(cfg, inputs, models)
%COMPUTE_DRYDEN_FREQUENCY_RESPONSE Compute modal PSD using ZAERO QHG.
% The configured RFA supplies the aeroelastic dynamics through Aae, and the
% configured QHG file supplies the frequency-dependent generalized gust
% forces. The responses are weighted by the Dryden PSD printed by ZAERO.

arguments
    cfg struct
    inputs struct
    models cell
end

%% ------------------------------------------------------------------------
% Read Dryden parameters and ZAERO input-PSD table
% -------------------------------------------------------------------------
zaeroOutputFile = cfg.ZAERO_filename;
if isfield(cfg,'drydenZaeroOutputFile') && ...
        ~isempty(cfg.drydenZaeroOutputFile)
    zaeroOutputFile = cfg.drydenZaeroOutputFile;
end
inputData = read_zaero_dryden_psd(zaeroOutputFile);

f      = inputData.frequency(:);
Lg     = inputData.Lg;
sigmaW = inputData.sigmaW;
Vzaero = inputData.V;

%% ------------------------------------------------------------------------
% Select the aeroelastic model at the requested validation velocity
% -------------------------------------------------------------------------
modelSpeeds = cellfun(@(m) m.V, models);
[~, iV] = min(abs(modelSpeeds - cfg.gustValidationSpeed));
model = models{iV};

V = model.V;

if abs(V - Vzaero) > 1e-6*max(1,abs(Vzaero))
    warning(['Selected MATLAB model speed is %.9g m/s, while ZAERO used ', ...
             '%.9g m/s.'], V, Vzaero);
end

if ~isfield(model,'Aae') || isempty(model.Aae)
    error('Selected model does not contain a nonempty Aae matrix.');
end

Nel = cfg.Nel;

if size(model.Aae,1) < 2*Nel
    error('Aae has %d states, but at least %d modal states are required.', ...
          size(model.Aae,1), 2*Nel);
end

%% ------------------------------------------------------------------------
% ZAERO generalized aerodynamic gust forces
% -------------------------------------------------------------------------
Qhg = Read_out4_to_mat_form(cfg.drydenQhgFile,'QHGSAV');

if size(Qhg,1) ~= Nel && size(Qhg,2) == Nel
    Qhg = Qhg.';
end

if size(Qhg,1) ~= Nel || size(Qhg,2) ~= numel(f)
    error('QHG dimensions are %d x %d; expected %d x %d.', ...
        size(Qhg,1),size(Qhg,2),Nel,numel(f));
end

%% ------------------------------------------------------------------------
% Rebuild Dryden analytically with the SAME parameters as ZAERO
% -------------------------------------------------------------------------
omega = 2*pi*f;
PhiW = build_dryden_psd(f,Lg,sigmaW,Vzaero);

% This is only a diagnostic: compare the analytical MATLAB Dryden input
% against the input PSD table printed by ZAERO.
inputPsdRelativeError = abs(PhiW-inputData.psd) ./ max(abs(inputData.psd),eps);

%% ------------------------------------------------------------------------
% Generalized-force input matrix for the structural RFA
% -------------------------------------------------------------------------
qdyn = 0.5*cfg.rho*V^2;
A2h = inputs.RFA_mat.A2(1:Nel,1:Nel);
Mbar = eye(Nel) - qdyn*(cfg.L/V)^2*A2h;
nLag = size(model.Aae,1) - 2*Nel;

Bq = [ ...
    zeros(Nel,Nel); ...
    Mbar\(qdyn*eye(Nel)); ...
    zeros(nLag,Nel)];

%% ------------------------------------------------------------------------
% Frequency response from physical w_g to aeroelastic states
% -------------------------------------------------------------------------
nFreq = numel(f);
nState = size(model.Aae,1);
I = eye(nState);

HmodalX = complex(zeros(nFreq,Nel));
HmodalV = complex(zeros(nFreq,Nel));

for j = 1:nFreq
    s = 1i*omega(j);

    forcing = Bq*(Qhg(:,j)/V);
    xPerW = (s*I - model.Aae) \ forcing;

    HmodalX(j,:) = xPerW(1:Nel).';
    HmodalV(j,:) = xPerW(Nel+(1:Nel)).';
end

%% ------------------------------------------------------------------------
% Modal response PSD
%
% IMPORTANT:
% Do NOT multiply the PSD itself by 2.
% ZAERO prints the positive-frequency half spectrum. The factor of 2 is
% applied only when integrating the PSD to obtain RMS.
% -------------------------------------------------------------------------
modalXPsd = abs(HmodalX).^2 .* PhiW;
modalVPsd = abs(HmodalV).^2 .* PhiW;

rmsModalX = sqrt(2*trapz(f,modalXPsd,1));
rmsModalV = sqrt(2*trapz(f,modalVPsd,1));

%% ------------------------------------------------------------------------
% Pack output
% -------------------------------------------------------------------------
response = struct();

response.frequency = f;

response.inputPsdMatlab = PhiW;
response.inputPsdZaero = inputData.psd;
response.inputPsdMaxRelativeError = max(inputPsdRelativeError);

response.modalXFrf = HmodalX;
response.modalVFrf = HmodalV;

response.modalXPsd = modalXPsd;
response.modalVPsd = modalVPsd;

response.rmsModalX = rmsModalX;
response.rmsModalV = rmsModalV;

response.V = V;
response.Vzaero = Vzaero;
response.Lg = Lg;
response.sigmaW = sigmaW;
response.Nel = Nel;
response.modelIndex = iV;
response.qhgFile = cfg.drydenQhgFile;
response.zaeroOutFile = zaeroOutputFile;

fprintf('\n============================================================\n');
fprintf('Dryden frequency-domain modal PSD calculation\n');
fprintf('============================================================\n');
fprintf('ZAERO V                  = %.6f m/s\n',Vzaero);
fprintf('MATLAB model V           = %.6f m/s\n',V);
fprintf('Lg                       = %.6f m\n',Lg);
fprintf('sigma_w                  = %.6f m/s\n',sigmaW);
fprintf('Frequency range          = %.6f to %.6f Hz\n',f(1),f(end));
fprintf('Number of frequencies    = %d\n',nFreq);
fprintf('QHG source               = %s\n',cfg.drydenQhgFile);
fprintf('Input PSD max rel. error = %.6e %%\n', ...
        100*response.inputPsdMaxRelativeError);
fprintf('============================================================\n\n');

end