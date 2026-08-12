function cfg = default_analysis_config()
%DEFAULT_ANALYSIS_CONFIG Central settings for the aeroelastic analysis.

cfg.rho = 1.225;
cfg.L = 0.05;        % length scale [m], usually half-chord
cfg.zeta = 0;        % modal damping ratio

cfg.modelType = 'plant'; % 'ae' or 'plant'
cfg.isgust = 0;      % fallback if gust count is not reported by ZAERO
cfg.numControlSurfaces = 4; % fallback if control-surface count is not reported by ZAERO
cfg.numGustInputs = []; % derived from ZAERO/RFA when gust columns exist

cfg.RFA_filename = 'APPROX.DAT';
cfg.ZAERO_filename = 'ASE_ANALYSIS_new.out';
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
cfg.gustModel = 'dryden';           % 'omcos' or 'dryden'
cfg.gustValidationSpeed = 35;       % [m/s]
cfg.gustResponseFile = 'GUSTXI.DAT';
cfg.gustLG = 1.0;                  % [s], DGUST OMCOSP LG
cfg.gustWGV = 0.1;                 % DGUST OMCOSP WGV = w_g/V
cfg.gustX0 = 0.0;                  % [m], DGUST gust reference point
cfg.gustXResponse = 0.0;           % [m], response location for gust delay
cfg.gustLengthMode = 'time';       % DGUST LG is a time duration
cfg.gustInternalDt = 0.001;        % [s], MATLAB integration step

% MIL-STD-1797A low-altitude Dryden turbulence conditions.
cfg.drydenAltitude = 1.5;             % [m] above local terrain
cfg.drydenMeanWindAt20ft = 7.0;     % [m/s], selected u_20 from Figure 263
cfg.drydenReferenceHeight = 20*0.3048; % [m], 20 ft
cfg.drydenRoughnessLength = 0.15*0.3048; % [m], Category C z_0 = 0.15 ft
cfg.drydenScaleLength = [];         % [m], empty uses altitude/2
cfg.drydenSigmaW = [];              % [m/s], empty uses sigma_w = 0.1*u_20
cfg.drydenTau = 0.002;              % [s], added high-frequency LPF
cfg.drydenTime = 0:0.001:60;        % [s]
cfg.drydenRandomSeed = 1;
cfg.drydenWindProfileMaxAltitude = 30; % [m], wind-parameter plot limit

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

end
