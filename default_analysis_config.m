function cfg = default_analysis_config()
%DEFAULT_ANALYSIS_CONFIG Central settings for the aeroelastic analysis.

cfg.rho = 1.225;
cfg.L = 0.05;        % length scale [m], usually half-chord
cfg.zeta = 0;        % modal damping ratio

cfg.modelType = 'plant'; % 'ae' or 'plant'
cfg.isgust = 0;      % fallback if gust count is not reported by ZAERO
cfg.numControlSurfaces = 4; % fallback if control-surface count is not reported by ZAERO
cfg.numGustInputs = []; % derived from ZAERO/RFA when gust columns exist

cfg.RFA_filename = 'APPROX_omcos.DAT';
cfg.ZAERO_filename = 'ASE_ANALYSIS_omcos.out';
cfg.f06file = 'model-0012.f06';
cfg.sensorCacheFile = 'sensor_modal_matrices.mat';
cfg.useSensorCache = true;
cfg.rebuildSensorCache = false;

cfg.method = 'MS';   % 'MS' or 'Rg'

% Retained/omitted FEM modes and RFA partitions are derived from ZAERO output.
cfg.sensor = default_sensor_config();

cfg.airspeed = 25:0.25:50;


cfg.plotBode = false;
cfg.bodeSpeed = 35;  % [m/s]
cfg.bodeFrequency = logspace(-1, 3, 500); % [rad/s]

% Diagnostic strain-mode plots. Requires cfg.datFile to exist.
cfg.plotStrainModes = false;
cfg.normalizeStrainModes = false;
cfg.datFile = 'model-0012.dat';

% Time-domain response checks
cfg.plotStepResponse = false;
cfg.responseCheckSpeed = 35;       % [m/s]
cfg.responseStepTime = 0:0.001:2;  % [s]
cfg.responseStepCommandDeg = 1;    % [deg]
cfg.responseStrongestStrainCount = 5;

% Gust validation against ZAERO MLDPRNT output
cfg.plotGustValidation = true;
cfg.gustModel = 'omcos';           % 'omcos' or 'dryden'
cfg.gustValidationSpeed = 35;       % [m/s]
cfg.gustResponseFile = 'GUSTXI';
cfg.gustLG = 1.0;                  % [s], DGUST OMCOSP LG
cfg.gustWGV = 0.1;                 % DGUST OMCOSP WGV = w_g/V
cfg.gustX0 = 0.0;                  % [m], DGUST gust reference point
cfg.gustXResponse = 0.0;           % [m], response location for gust delay
cfg.gustLengthMode = 'time';       % DGUST LG is a time duration
cfg.gustInternalDt = 0.001;        % [s], MATLAB integration step

% Continuous Dryden comparison against ZAERO.
% Lg, GURMS, V and the frequency grid are read from the ZAERO output.
cfg.drydenQhgFile = 'QHG.DAT';          % generalized gust forces from ZAERO
cfg.drydenZaeroOutputFile = 'ASE_ANALYSIS_DRYDEN.out';
cfg.compareDrydenZaero = true;
cfg.drydenZaeroPsdFile = 'DRYPSD.P';       % explicit CROSPSD auto-PSD output
cfg.drydenTau = 1e-4;                    % [s], LPF on gust-rate output only
cfg.drydenPsdDuration = 60;              % [s], stationary PSD realization
cfg.drydenPsdTimeStep = 1e-3;            % [s]
cfg.drydenWelchWindowLength = 8192;      % samples
cfg.drydenRandomSeed = 1;                % reproducible displayed realization

% Plotting settings

cfg.rootLocusFigure = 101;
cfg.plotOmegaVG = false;
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

% Analysis selection and output behavior.
cfg.analysis.rootLocus = true;
cfg.analysis.omegaVG = cfg.plotOmegaVG;
cfg.analysis.stepResponse = cfg.plotStepResponse;
cfg.analysis.strainModes = cfg.plotStrainModes;
cfg.analysis.bode = cfg.plotBode;
cfg.analysis.gust = cfg.plotGustValidation;
cfg.analysis.gustModel = cfg.gustModel;

cfg.output.showFigures = true;
cfg.output.saveFigures = false;
cfg.output.figureDirectory = fullfile('output', 'figures');
cfg.output.figureResolution = 300;

end

