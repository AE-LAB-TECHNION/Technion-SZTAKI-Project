function peakTable = run_response_checks(cfg, inputs, models)
%RUN_RESPONSE_CHECKS Run optional time-domain and diagnostic response plots.

peakTable = table();

runStepResponse = cfg_value(cfg, 'plotStepResponse', false);
runStrainModes = cfg_value(cfg, 'plotStrainModes', false);

if ~runStepResponse && ~runStrainModes
    return
end

if runStepResponse
    peakTable = plot_step_response_checks(cfg, inputs, models);
else
    fprintf('\nStep response checks skipped because cfg.plotStepResponse = false.\n');
end

if runStrainModes
    plot_strain_modes(cfg, inputs);
end

end

function value = cfg_value(cfg, fieldName, defaultValue)

if isfield(cfg, fieldName)
    value = cfg.(fieldName);
else
    value = defaultValue;
end

end
