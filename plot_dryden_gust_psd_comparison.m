function result = plot_dryden_gust_psd_comparison( ...
    cfg, augmented, zaeroOutputFile)
%PLOT_DRYDEN_GUST_PSD_COMPARISON Compare three Dryden gust spectra.
%
% The single figure contains:
%   1. Welch PSD estimate from a stationary Dryden realization;
%   2. the analytical MATLAB Dryden spectrum;
%   3. the Dryden spectrum printed by ZAERO.

arguments
    cfg struct
    augmented struct
    zaeroOutputFile (1,:) char
end

duration = 60;
if isfield(cfg,'drydenPsdDuration') && ~isempty(cfg.drydenPsdDuration)
    duration = cfg.drydenPsdDuration;
end

dt = 1e-3;
if isfield(cfg,'drydenPsdTimeStep') && ~isempty(cfg.drydenPsdTimeStep)
    dt = cfg.drydenPsdTimeStep;
end

seed = 1;
if isfield(cfg,'drydenRandomSeed') && ~isempty(cfg.drydenRandomSeed)
    seed = cfg.drydenRandomSeed;
end

if duration <= 0 || dt <= 0
    error('Dryden PSD duration and time step must be positive.');
end

zaero = read_zaero_dryden_psd(zaeroOutputFile);

dryden = augmented.dryden;
A = dryden.A;
G = dryden.inputNoiseScale*dryden.B;
C = dryden.C;

% Draw an exactly discretized stationary realization of the continuous
% filter. Sampling the initial state from P avoids a startup transient.
P = lyap(A,G*G.');
Ad = expm(A*dt);
Qd = P-Ad*P*Ad.';

P = (P+P.')/2;
Qd = (Qd+Qd.')/2;
Lp = positive_semidefinite_factor(P);
Lq = positive_semidefinite_factor(Qd);

time = (0:dt:duration).';
nSamples = numel(time);
nState = size(A,1);

rng(seed,'twister');
state = zeros(nSamples,nState);
state(1,:) = (Lp*randn(nState,1)).';

for k = 1:nSamples-1
    state(k+1,:) = ( ...
        Ad*state(k,:).' + Lq*randn(nState,1)).';
end

gust = state*C.';
wg = gust(:,1);
wdotg = gust(:,2);
realizedRms = sqrt(mean(wg.^2));

% MATLAB returns a one-sided PSD whose integral equals the variance.
% ZAERO uses a positive-frequency half spectrum, for which the variance is
% twice the positive-frequency integral. Divide the Welch estimate by two
% so that all three curves use the same convention.
samplingFrequency = 1/dt;
defaultWindowLength = min(8192,2^floor(log2(max(256,nSamples/8))));
windowLength = defaultWindowLength;
if isfield(cfg,'drydenWelchWindowLength') && ...
        ~isempty(cfg.drydenWelchWindowLength)
    windowLength = cfg.drydenWelchWindowLength;
end

windowLength = min(floor(windowLength),nSamples);
if windowLength < 16
    error('The Dryden Welch window must contain at least 16 samples.');
end

windowIndex = (0:windowLength-1).';
window = 0.5-0.5*cos(2*pi*windowIndex/(windowLength-1));
overlap = floor(0.5*windowLength);
nfft = max(32768,2^nextpow2(windowLength));

[oneSidedPsd,frequency] = pwelch( ...
    wg,window,overlap,nfft,samplingFrequency,'onesided');
estimatedPsd = oneSidedPsd/2;
matlabPsd = build_dryden_psd( ...
    frequency,augmented.Lg,augmented.sigmaW,augmented.V);

positiveEstimate = frequency > 0;
positiveZaero = zaero.frequency > 0;

% Keep the ZAERO markers readable on the dense 0.01-Hz grid.
nZaeroMarkers = 180;
zaeroIndices = find(positiveZaero);
markerSkip = max(1,floor(numel(zaeroIndices)/nZaeroMarkers));
zaeroMarkerIndices = zaeroIndices(1:markerSkip:end);

fig = figure( ...
    'Name','Dryden gust PSD comparison', ...
    'Color','white');

loglog( ...
    frequency(positiveEstimate),estimatedPsd(positiveEstimate), ...
    'Color',[0.0000,0.4470,0.7410], ...
    'LineWidth',cfg.plotLineWidth);
hold on;
loglog( ...
    frequency(positiveEstimate),matlabPsd(positiveEstimate),'--', ...
    'Color',[0.8500,0.3250,0.0980], ...
    'LineWidth',1.5*cfg.plotLineWidth);
loglog( ...
    zaero.frequency(zaeroMarkerIndices), ...
    zaero.psd(zaeroMarkerIndices), ...
    'o', ...
    'Color',[0.1500,0.1500,0.1500], ...
    'MarkerFaceColor','none', ...
    'MarkerSize',4, ...
    'LineWidth',0.9);

grid on;
box on;
xlabel('Frequency [Hz]');
ylabel('PSD of w_g [(m/s)^2/Hz]');
legend( ...
    'Estimated (Welch)', ...
    'MATLAB Dryden model', ...
    'ZAERO', ...
    'Location','southwest');
xlim([max(frequency(find(positiveEstimate,1)),0.01),zaero.frequency(end)]);
apply_plot_style(fig,cfg);

% Compare all RMS values over the common ZAERO frequency interval.
inZaeroBand = frequency <= zaero.frequency(end);
estimatedBandRms = sqrt(2*trapz( ...
    frequency(inZaeroBand),estimatedPsd(inZaeroBand)));
matlabBandRms = sqrt(2*trapz( ...
    frequency(inZaeroBand),matlabPsd(inZaeroBand)));
zaeroBandRms = sqrt(2*trapz(zaero.frequency,zaero.psd));

matlabAtZaero = build_dryden_psd( ...
    zaero.frequency,augmented.Lg,augmented.sigmaW,augmented.V);
relativeError = abs(matlabAtZaero-zaero.psd) ./ ...
    max(abs(zaero.psd),eps);

result = struct();
result.time = time;
result.wg = wg;
result.wdotg = wdotg;
result.state = state;
result.frequency = frequency;
result.estimatedPsd = estimatedPsd;
result.matlabPsd = matlabPsd;
result.zaeroFrequency = zaero.frequency;
result.zaeroPsd = zaero.psd;
result.realizedRmsWg = realizedRms;
result.estimatedBandRmsWg = estimatedBandRms;
result.matlabBandRmsWg = matlabBandRms;
result.zaeroBandRmsWg = zaeroBandRms;
result.maxMatlabZaeroRelativeError = max(relativeError);
result.theoreticalRmsWg = augmented.sigmaW;
result.dt = dt;
result.duration = duration;
result.seed = seed;
result.windowLength = windowLength;
result.zaeroOutputFile = zaeroOutputFile;

fprintf('\nCombined Dryden gust PSD comparison\n');
fprintf('  realization duration       = %.6f s\n',duration);
fprintf('  time step                   = %.6g s\n',dt);
fprintf('  Welch window length         = %d samples\n',windowLength);
fprintf('  random seed                 = %d\n',seed);
fprintf('  theoretical RMS(w_g)        = %.9g m/s\n',augmented.sigmaW);
fprintf('  realization RMS(w_g)        = %.9g m/s\n',realizedRms);
fprintf('  estimated RMS to ZAERO Fmax = %.9g m/s\n',estimatedBandRms);
fprintf('  MATLAB RMS to ZAERO Fmax    = %.9g m/s\n',matlabBandRms);
fprintf('  ZAERO RMS to Fmax           = %.9g m/s\n',zaeroBandRms);
fprintf('  max MATLAB--ZAERO PSD error = %.6e %%\n', ...
    100*result.maxMatlabZaeroRelativeError);

end

% =========================================================================
function L = positive_semidefinite_factor(P)

[vectors,values] = eig((P+P.')/2,'vector');
values = real(values);

tolerance = 1e3*eps(max(1,max(abs(values))));
if min(values) < -tolerance
    error('Covariance matrix is not positive semidefinite.');
end

values(values < 0) = 0;
L = real(vectors*diag(sqrt(values)));

end
