function results = plot_analysis_results(cfg, results)
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
if isfield(results, 'drydenFrequency')
    if strcmpi(cfg.gustModel,'dryden')
        zaeroDrydenOutput = cfg.ZAERO_filename;
        if isfield(cfg,'drydenZaeroOutputFile') && ...
                ~isempty(cfg.drydenZaeroOutputFile)
            zaeroDrydenOutput = cfg.drydenZaeroOutputFile;
        end

        if isfield(results,'drydenStateSpace')
            results.drydenGustPsdComparison = ...
                plot_dryden_gust_psd_comparison( ...
                    cfg,results.drydenStateSpace,zaeroDrydenOutput);
        end
    end

    plot_dryden_frequency_psd(cfg,results.drydenFrequency);

    if strcmpi(cfg.gustModel,'dryden') && cfg.compareDrydenZaero
        if ~isfile(cfg.drydenZaeroPsdFile)
            warning('ZAERO Dryden PSD file not found: %s', ...
                cfg.drydenZaeroPsdFile);
        else
            results.drydenPsdComparison = plot_dryden_zaero_comparison( ...
                cfg,results.drydenFrequency);
        end
    end
end

end
