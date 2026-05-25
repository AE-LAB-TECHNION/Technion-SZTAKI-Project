function cfg = default_analysis_config()
%DEFAULT_ANALYSIS_CONFIG Central settings for the aeroelastic analysis.

cfg.rho = 1.225;
cfg.L = 0.05;        % length scale [m], usually half-chord
cfg.zeta = 0;        % modal damping ratio

cfg.isinp = 0;       % 1 = plot Ap, 0 = plot Aae
cfg.isgust = 0;      % 1 = include gust column, 0 = no gust

cfg.nmodes = 4;      % number of structural modes to read from ZAERO output

cfg.RFA_filename = 'APPROX.DAT';
cfg.ZAERO_filename = 'ASE_ANALYSIS_new.out';
cfg.f06file = 'model-0012.f06';

cfg.method = 'MS';   % 'MS' or 'Rg'

% ZAERO uses modes 1,2,4 because mode 3 is omitted.
cfg.selected_modes = [1 2 4];
cfg.Nel = numel(cfg.selected_modes);

cfg.sensor = default_sensor_config();

cfg.airspeed = 25:0.25:50;

cfg.rootLocusFigure = 101;
cfg.rootLocusXLim = [-100, 100];
cfg.rootLocusYLim = [0, 200];

cfg.bodeSpeed = 35;  % [m/s]
cfg.bodeFrequency = logspace(-1, 3, 500); % [rad/s]

end
