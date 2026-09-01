function plot_analysis_results(cfg, results)
%PLOT_ANALYSIS_RESULTS Create figures for enabled analyses.

if ~cfg.output.showFigures
    return
end

if cfg.analysis.rootLocus
    plot_root_locus(cfg, results.eigTraj, results.plotModelName);
end
if cfg.plotOmegaVG
    plot_omega_v_g(cfg, results.models, 'Aae');
end
if cfg.plotBode
    plot_control_bodes(cfg, results.models);
end
if isfield(results, 'dryden')
    plot_dryden_results(cfg, results.dryden);
    if cfg.compareDrydenZaero && isfile(cfg.drydenZaeroResponseFile)
        plot_dryden_zaero_comparison( ...
            cfg,results.dryden,cfg.drydenZaeroResponseFile);
    end
end
if isfield(results, 'drydenFrequency')
    plot_dryden_frequency_psd(cfg,results.drydenFrequency);
end
if isfield(results, 'drydenFrfComparison')
    plot_dryden_frf_comparison(cfg,results.drydenFrfComparison);
end

end
