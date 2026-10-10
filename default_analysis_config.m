function cfg = default_analysis_config()
%DEFAULT_ANALYSIS_CONFIG Central settings for the aeroelastic analysis.

projectRoot = fileparts(mfilename('fullpath'));
modelData = fullfile(projectRoot,'data','model');
zaeroComparison = fullfile(projectRoot,'zaero_model_comparison');

cfg.rho = 1.225;
cfg.L = 0.05;        % length scale [m], usually half-chord
cfg.zeta = 0;        % modal damping ratio

cfg.modelType = 'plant'; % 'ae' or 'plant'
cfg.numControlSurfaces = 4; % fallback if control-surface count is not reported by ZAERO
cfg.numGustInputs = []; % derived from ZAERO/RFA when gust columns exist

cfg.RFA_filename = fullfile(modelData,'APPROX.DAT');
cfg.ZAERO_filename = fullfile(modelData,'ASE_ANALYSIS.OUT');
cfg.f06file = fullfile(modelData,'model-0012.f06');
cfg.sensorCacheFile = fullfile(modelData,'sensor_modal_matrices.mat');
cfg.useSensorCache = true;
cfg.rebuildSensorCache = false;
cfg.rbmModelFile = fullfile(modelData,'AERBM.S');
cfg.rbmOutputRow = 4;

cfg.method = 'MS';   % 'MS' or 'Rg'

% Retained/omitted FEM modes and RFA partitions are derived from ZAERO output.
cfg.sensor = default_sensor_config();

cfg.airspeed = 25:0.25:50;

% Analysis selection
cfg.analysis.rootLocus = true;
cfg.analysis.omegaVG = false;
cfg.analysis.strainModes = false;
cfg.analysis.gust = true;

% System input-output response figures
cfg.responsePlots.modalBode = false;
cfg.responsePlots.rbmBode = true;
cfg.responsePlots.stepRate = false;
cfg.responsePlots.stepDisp = false;
cfg.responsePlots.stepStrain = false;
cfg.responsePlots.stepRbm = false;

cfg.bodeSpeed = 35;  % [m/s]
cfg.bodeFrequency = logspace(-1, 3, 500); % [rad/s]

% Diagnostic strain-mode settings. Requires cfg.datFile to exist.
cfg.normalizeStrainModes = false;
cfg.datFile = fullfile(modelData,'MODEL-0012.dat');

% Time-domain response settings
cfg.responseCheckSpeed = 35;       % [m/s]
cfg.responseStepTime = 0:0.001:2;  % [s]
cfg.responseStepCommandDeg = 1;    % [deg]
cfg.responseStrongestStrainCount = 5;

% Gust analysis
cfg.gustModel = 'dryden';           % 'omcos' or 'dryden'
cfg.gustValidationSpeed = 35;       % [m/s]
cfg.gustResponseFile = fullfile( ...
    zaeroComparison,'omcos','GUSTXI');
cfg.gustLG = 1.0;                  % [s], DGUST OMCOSP LG
cfg.gustWGV = 0.1;                 % DGUST OMCOSP WGV = w_g/V
cfg.gustX0 = 0.0;                  % [m], DGUST gust reference point
cfg.gustXResponse = 0.0;           % [m], response location for gust delay
cfg.gustLengthMode = 'time';       % DGUST LG is a time duration
cfg.gustInternalDt = 0.001;        % [s], MATLAB integration step

% Continuous Dryden comparison against ZAERO.
% The model parameters are explicit; ZAERO files are used only for validation.
cfg.drydenVelocity = 35;                 % [m/s]
cfg.drydenScaleLength = 0.75;            % Lg [m]
cfg.drydenSigmaW = 0.7;                  % GURMS [m/s]
cfg.drydenFrequencyStart = 0;             % [Hz]
cfg.drydenFrequencyStep = 0.01;           % [Hz]
cfg.drydenFrequencyEnd = 100;             % [Hz]
cfg.drydenZaeroOutputFile = fullfile( ...
    zaeroComparison,'dryden','ASE_ANALYSIS_DRY_RBM.out');
cfg.compareDrydenZaero = true;
cfg.drydenZaeroPsdFile = fullfile( ...
    zaeroComparison,'dryden','DRYPSD.P');
cfg.drydenTau = 1e-4;                    % [s], LPF on gust-rate output only
cfg.drydenPsdDuration = 60;              % [s], stationary PSD realization
cfg.drydenPsdTimeStep = 1e-3;            % [s]
cfg.drydenWelchWindowLength = 8192;      % samples
cfg.drydenRandomSeed = 1;                % reproducible displayed realization

% Plotting settings

cfg.rootLocusFigure = 101;
cfg.omegaVGFigure = 102;
cfg.omegaVGDampingYLim = [-0.5, 0.5];
cfg.rootLocusXLim = [-100, 100];
cfg.rootLocusYLim = [0, 200];

cfg.plotFontName = 'Helvetica';
cfg.plotFontSize = 13;
cfg.plotTitleFontSize = 14;
cfg.plotTitleFontWeight = 'normal';
cfg.plotAxesFontWeight = 'normal';
cfg.plotAxesLineWidth = 1.0;
cfg.plotLineWidth = 1.5;
cfg.plotTextColor = 'black';
cfg.plotTextInterpreter = 'tex';
cfg.plotLegendInterpreter = cfg.plotTextInterpreter;
cfg.plotTickLabelInterpreter = cfg.plotTextInterpreter;

% Output behavior
cfg.output.showFigures = true;
cfg.output.saveFigures = false;
cfg.output.figureDirectory = fullfile(projectRoot,'output','figures');
cfg.output.figureResolution = 300;

end