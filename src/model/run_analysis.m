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

drydenStateSpace = [];

% Build the control model independently of ZAERO validation files.
if strcmpi(cfg.gustModel,'dryden')
    results.drydenStateSpace = build_dryden_state_space(cfg,models);
    drydenStateSpace = results.drydenStateSpace;
end

%% ========================================================================
% Gust validation
% ========================================================================

if cfg.plotGustValidation

    switch lower(char(cfg.gustModel))

        % ---------------------------------------------------------------
        % One-minus-cosine gust
        % ---------------------------------------------------------------
        case 'omcos'

            results.gust = plot_gust_validation( ...
                cfg, inputs, models);


        % ---------------------------------------------------------------
        % Continuous Dryden gust
        % ---------------------------------------------------------------
        case 'dryden'
            results.drydenFrequency = ...
                compute_augmented_dryden_frequency_response( ...
                    cfg, ...
                    results.drydenStateSpace);


        % ---------------------------------------------------------------
        % Invalid option
        % ---------------------------------------------------------------
        otherwise

            error( ...
                'Unknown cfg.gustModel "%s". Use ''omcos'' or ''dryden''.', ...
                char(cfg.gustModel));

    end

end

results.rbmValidation = validate_rbm_measurement( ...
    inputs,models,drydenStateSpace);

end
