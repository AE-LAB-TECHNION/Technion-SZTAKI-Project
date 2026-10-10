function response = compute_augmented_dryden_frequency_response(cfg, augmented)
%COMPUTE_AUGMENTED_DRYDEN_FREQUENCY_RESPONSE Compute modal PSD from the full model.
%
% The calculation follows the independent augmented-state-space path
%
%   white noise -> Dryden filter -> Bpw -> plant -> modal response
%
% and does not use ZAERO QHG or Dryden output files.

arguments
    cfg struct
    augmented struct
end

frequency = (cfg.drydenFrequencyStart: ...
    cfg.drydenFrequencyStep:cfg.drydenFrequencyEnd).';

if isempty(frequency) || cfg.drydenFrequencyStep <= 0
    error('The configured Dryden frequency grid is invalid.');
end

Nel = cfg.Nel;
nState = size(augmented.A,1);

if nState < 2*Nel
    error('The augmented model has %d states; at least %d are required.', ...
        nState,2*Nel);
end

if size(augmented.G,2) ~= 1
    error('The augmented Dryden model must have one stochastic input.');
end

% The first 2*Nel plant states are [xi; dxi/dt].
Cmodal = zeros(2*Nel,nState);
Cmodal(:,1:2*Nel) = eye(2*Nel);

modalSystem = ss( ...
    augmented.A,augmented.G,Cmodal,zeros(2*Nel,1));

frequencyResponse = freqresp(modalSystem,2*pi*frequency);
frequencyResponse = permute(frequencyResponse,[3,1,2]);
frequencyResponse = reshape( ...
    frequencyResponse,numel(frequency),2*Nel);

HmodalX = frequencyResponse(:,1:Nel);
HmodalV = frequencyResponse(:,Nel+(1:Nel));

modalXPsd = abs(HmodalX).^2;
modalVPsd = abs(HmodalV).^2;

% ZAERO prints the positive-frequency half spectrum.
rmsModalX = sqrt(2*trapz(frequency,modalXPsd,1));
rmsModalV = sqrt(2*trapz(frequency,modalVPsd,1));

response = struct();
response.frequency = frequency;
response.modalXFrf = HmodalX;
response.modalVFrf = HmodalV;
response.modalXPsd = modalXPsd;
response.modalVPsd = modalVPsd;
response.rmsModalX = rmsModalX;
response.rmsModalV = rmsModalV;
response.V = augmented.V;
response.Lg = augmented.Lg;
response.sigmaW = augmented.sigmaW;
response.Nel = Nel;
response.modelIndex = augmented.modelIndex;
response.calculationPath = 'Dryden filter -> Bpw -> augmented plant';

fprintf('\n============================================================\n');
fprintf('Independent augmented-model Dryden modal PSD\n');
fprintf('============================================================\n');
fprintf('Velocity                   = %.6f m/s\n',augmented.V);
fprintf('Lg                         = %.6f m\n',augmented.Lg);
fprintf('sigma_w                    = %.6f m/s\n',augmented.sigmaW);
fprintf('Frequency range            = %.6f to %.6f Hz\n', ...
    frequency(1),frequency(end));
fprintf('Number of frequencies      = %d\n',numel(frequency));
fprintf('ZAERO QHG used             = no\n');
fprintf('ZAERO Dryden OUT used      = no\n');
fprintf('============================================================\n\n');

end
