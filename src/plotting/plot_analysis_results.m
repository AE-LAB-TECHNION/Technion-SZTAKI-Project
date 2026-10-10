function results = plot_analysis_results(cfg, results)
%PLOT_ANALYSIS_RESULTS Create figures for enabled analyses.

if ~cfg.output.showFigures
    return
end

if cfg.analysis.rootLocus
    plot_root_locus(cfg, results.eigTraj, results.plotModelName);
end
if cfg.analysis.omegaVG
    plot_omega_v_g(cfg, results.models, 'Aae');
end
if cfg.analysis.strainModes
    plot_strain_modes(cfg, results.inputs);
end

results.responsePlots = run_response_checks( ...
    cfg, results.inputs, results.models);

if isfield(results, 'drydenFrequency')
    plot_dryden_frequency_psd(cfg,results.drydenFrequency);

    if strcmpi(cfg.gustModel,'dryden') && cfg.compareDrydenZaero
        zaeroDrydenOutput = cfg.ZAERO_filename;
        if isfield(cfg,'drydenZaeroOutputFile') && ...
                ~isempty(cfg.drydenZaeroOutputFile)
            zaeroDrydenOutput = cfg.drydenZaeroOutputFile;
        end

        if ~isfile(zaeroDrydenOutput)
            warning('ZAERO Dryden output file not found: %s', ...
                zaeroDrydenOutput);
        elseif isfield(results,'drydenStateSpace')
            results.drydenGustPsdComparison = ...
                plot_dryden_gust_psd_comparison( ...
                    cfg,results.drydenStateSpace,zaeroDrydenOutput);
        end

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
