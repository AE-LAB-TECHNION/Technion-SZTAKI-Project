clc; clear; close all;
clear functions;

cfg = default_analysis_config();

[inputs, cfg] = load_analysis_inputs(cfg);

[models, eigTraj, plotModelName] = run_airspeed_sweep(cfg, inputs);

plot_root_locus(cfg, eigTraj, plotModelName);
if cfg.plotOmegaVG
    plot_omega_v_g(cfg, models, 'Aae');
end

[Vf, ff] = report_flutter(cfg, models, 'Aae physical modes');

run_response_checks(cfg, inputs, models);

if cfg.plotBode
    plot_control_bodes(cfg, models);
end

if cfg.plotGustValidation
    switch lower(char(cfg.gustModel))
        case 'omcos'
            plot_gust_validation(cfg, inputs, models);
        case 'dryden'
            run_dryden_response(cfg, inputs, models);
        otherwise
            error('Unknown cfg.gustModel "%s". Use ''omcos'' or ''dryden''.', ...
                char(cfg.gustModel));
    end
end
