clc; clear; close all;

% -------------------------------------------------------------------------
% Defaults
% -------------------------------------------------------------------------
run('./defaults.m')

cfg = default_analysis_config();

timerInputs = tic;
inputs = load_analysis_inputs(cfg);
print_elapsed_time(cfg, 'Input loading', timerInputs);

timerSweep = tic;
[models, eigTraj, plotModelName] = run_airspeed_sweep(cfg, inputs);
print_elapsed_time(cfg, 'Airspeed sweep', timerSweep);

timerRootLocus = tic;
plot_root_locus(cfg, eigTraj, plotModelName);
print_elapsed_time(cfg, 'Root-locus plot', timerRootLocus);

[Vf, ff] = report_flutter(cfg, eigTraj, plotModelName);

timerBode = tic;
plot_control_bodes(cfg, models);
print_elapsed_time(cfg, 'Bode plots', timerBode);

% -------------------------------------------------------------------------
function print_elapsed_time(cfg, label, timerValue)
if isfield(cfg, 'printTimings') && cfg.printTimings
    fprintf('%s completed in %.3f s\n', label, toc(timerValue));
end
end
