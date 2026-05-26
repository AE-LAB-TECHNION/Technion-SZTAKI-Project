clc; clear; close all;
clear functions;

cfg = default_analysis_config();

[inputs, cfg] = load_analysis_inputs(cfg);

[models, eigTraj, plotModelName] = run_airspeed_sweep(cfg, inputs);

plot_root_locus(cfg, eigTraj, plotModelName);
plot_omega_v_g(cfg, models, 'Aae');

[Vf, ff] = report_flutter(cfg, models, 'Aae physical modes');

plot_control_bodes(cfg, models);
