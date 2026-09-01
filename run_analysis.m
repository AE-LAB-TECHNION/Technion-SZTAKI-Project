function results = run_analysis(cfg)
%RUN_ANALYSIS Load inputs, build models, and run enabled analyses.

[inputs, cfg] = load_analysis_inputs(cfg);
[models, eigTraj, plotModelName] = run_airspeed_sweep(cfg, inputs);

results = struct();
results.cfg = cfg;
results.inputs = inputs;
results.models = models;
results.eigTraj = eigTraj;
results.plotModelName = plotModelName;
[results.flutterSpeed, results.flutterFrequency] = ...
    report_flutter(cfg, models, 'Aae physical modes');
results.responsePeaks = run_response_checks(cfg, inputs, models);

if cfg.plotGustValidation
    switch lower(char(cfg.gustModel))
        case 'omcos'
            results.gust = plot_gust_validation(cfg, inputs, models);
        case 'dryden'
            [~,iV] = min(abs(cfg.airspeed-cfg.gustValidationSpeed));
            if ~isempty(models{iV}.Bpw)
                results.dryden = run_dryden_response(cfg, inputs, models);
            elseif isfile(cfg.drydenQhgFile)
                frequency = (cfg.drydenFrequencyStart: ...
                    cfg.drydenFrequencyStep:cfg.drydenFrequencyEnd).';
                results.drydenFrequency = compute_dryden_frequency_response( ...
                    cfg,inputs,models,cfg.drydenQhgFile,frequency);
            else
                error(['Dryden analysis requires either gust columns in ', ...
                    'APPROX.DAT or the generalized-force file "%s".'], ...
                    cfg.drydenQhgFile);
            end
            if cfg.compareDrydenFrf && isfile(cfg.drydenFrfFile)
                if ~isfield(results,'drydenFrequency')
                    error(['The ZAERO FRF comparison currently requires ', ...
                        'the MATLAB QHG-based frequency response.']);
                end
                results.drydenFrfComparison = compare_dryden_frf( ...
                    results.drydenFrequency,cfg.drydenFrfFile);
            end
        otherwise
            error('Unknown cfg.gustModel "%s". Use ''omcos'' or ''dryden''.', ...
                char(cfg.gustModel));
    end
end

end
