function comparison = show_dryden_zaero_comparison()
%SHOW_DRYDEN_ZAERO_COMPARISON Display one focused MATLAB--ZAERO comparison.
% The function does not save the figure automatically.

cfg = default_analysis_config();
cfg.RFA_filename = 'APPROX.DAT';
cfg.compareDrydenFrf = false;
cfg.output.saveFigures = false;

[inputs,cfg] = load_analysis_inputs(cfg);
[models,~,~] = run_airspeed_sweep(cfg,inputs);
dryden = run_dryden_response(cfg,inputs,models);
comparison = plot_dryden_zaero_comparison( ...
    cfg,dryden,cfg.drydenZaeroResponseFile);
end
