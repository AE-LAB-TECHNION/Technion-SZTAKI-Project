clc; clear; close all;

% -------------------------------------------------------------------------
% Defaults
% -------------------------------------------------------------------------
run('./defaults.m')

cfg = default_analysis_config();

inputs = load_analysis_inputs(cfg);

[models, eigTraj, plotModelName] = run_airspeed_sweep(cfg, inputs);

plot_root_locus(cfg, eigTraj, plotModelName);

[Vf, ff] = report_flutter(cfg, eigTraj, plotModelName);

plot_control_bodes(cfg, models);
