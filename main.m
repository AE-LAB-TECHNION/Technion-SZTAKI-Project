clc; clear; close all;

addpath(genpath('src'));
addpath(genpath('zaero_model_comparison'));

cfg = default_analysis_config();
results = run_analysis(cfg);
results = plot_analysis_results(cfg, results);
save_analysis_figures(cfg);