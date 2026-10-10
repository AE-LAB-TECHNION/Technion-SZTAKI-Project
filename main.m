clc; clear; close all;
clear functions;

projectRoot = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(projectRoot,'src')));
addpath(genpath(fullfile(projectRoot,'zaero_model_comparison')));

cfg = default_analysis_config();
results = run_analysis(cfg);
results = plot_analysis_results(cfg, results);
save_analysis_figures(cfg);
