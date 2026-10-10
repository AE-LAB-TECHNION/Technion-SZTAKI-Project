clc; clear; close all;
clear functions;

addpath(genpath(fullfile('src')));
addpath(genpath(fullfile('zaero_model_comparison')));

cfg = default_analysis_config();
results = run_analysis(cfg);
results = plot_analysis_results(cfg, results);
save_analysis_figures(cfg);