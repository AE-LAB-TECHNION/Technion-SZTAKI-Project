function response = run_response_checks(cfg, inputs, models)
%RUN_RESPONSE_CHECKS Plot selected system input-output responses.

response = struct();
selection = cfg.responsePlots;

if selection.modalBode
    response.modalBodePeaks = plot_control_bodes(cfg, models);
end

if selection.rbmBode
    response.rbmBode = plot_rbm_control_response(cfg, models);
end

runStep = selection.stepRate || selection.stepDisp || ...
    selection.stepStrain || selection.stepRbm;

if runStep
    response.stepPeaks = plot_step_response_checks(cfg, inputs, models);
end

end
