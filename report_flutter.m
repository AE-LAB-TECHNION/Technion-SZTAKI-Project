function [Vf, ff] = report_flutter(cfg, eigSource, plotModelName)
%REPORT_FLUTTER Estimate and print flutter speed/frequency.

[Vf, ff] = find_flutter(cfg.airspeed, eigSource, cfg.Nel);

fprintf('\nFlutter estimate based on %s eigenvalues:\n', plotModelName);
fprintf('Flutter occurs near V = %.3f m/s\n', Vf);
fprintf('Flutter frequency approx %.3f Hz\n', ff);

end
