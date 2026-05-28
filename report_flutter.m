function [Vf, ff] = report_flutter(cfg, ~, ~)
%REPORT_FLUTTER Estimate and print flutter speed/frequency.

data = read_zaero_vgf_data(cfg.ZAERO_filename, cfg.Nel);
Vf = data.flutter.Vf;
ff = data.flutter.ff;

fprintf('\nFlutter estimate based on ZAERO V-G-F data from %s:\n', cfg.ZAERO_filename);
print_flutter_result(Vf, ff);

end

% =========================================================================
function print_flutter_result(Vf, ff)

if isnan(Vf) || isnan(ff)
    fprintf('No damping zero crossing found.\n');
    return
end

fprintf('Flutter occurs near V = %.3f m/s\n', Vf);
fprintf('Flutter frequency approx %.3f Hz\n', ff);

end
