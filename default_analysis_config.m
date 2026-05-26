function cfg = default_analysis_config()
%DEFAULT_ANALYSIS_CONFIG Central settings for the aeroelastic analysis.

cfg.rho = 1.225;
cfg.L = 0.05;        % length scale [m], usually half-chord
cfg.zeta = 0;        % modal damping ratio

cfg.modelType = 'plant'; % 'ae' or 'plant'
cfg.isgust = 0;      % fallback if gust count is not reported by ZAERO
cfg.numControlSurfaces = 4; % fallback if control-surface count is not reported by ZAERO

cfg.RFA_filename = 'APPROX.DAT';
cfg.ZAERO_filename = 'ASE_ANALYSIS.out';
cfg.f06file = 'model-0012.f06';
cfg.sensorCacheFile = 'sensor_modal_matrices.mat';
cfg.useSensorCache = true;
cfg.rebuildSensorCache = false;

cfg.method = 'MS';   % 'MS' or 'Rg'

% Retained/omitted FEM modes and RFA partitions are derived from ZAERO output.
cfg.sensor = default_sensor_config();

cfg.airspeed = 25:0.25:50;

% Plotting settings
cfg.rootLocusFigure = 101;
cfg.omegaVGFigure = 102;
cfg.rootLocusXLim = [-100, 100];
cfg.rootLocusYLim = [0, 200];

cfg.bodeSpeed = 35;  % [m/s]
cfg.bodeFrequency = logspace(-1, 3, 500); % [rad/s]

cfg.plotFontName = 'Helvetica';
cfg.plotFontSize = 13;
cfg.plotTitleFontSize = 14;
cfg.plotTitleFontWeight = 'normal';
cfg.plotTextInterpreter = 'tex';

end
