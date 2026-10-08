function comparison = show_dryden_zaero_comparison(filename)
%SHOW_DRYDEN_ZAERO_COMPARISON Compare ZAERO and MATLAB Dryden PSD.
%
% Example:
%
%   comparison = show_dryden_zaero_comparison( ...
%       'ASE_ANALYSIS_DRYDEN.out');
%
% If no filename is supplied:
%
%   show_dryden_zaero_comparison
%
% uses ASE_ANALYSIS_DRYDEN.out.

    if nargin < 1 || isempty(filename)
        filename = 'ASE_ANALYSIS_DRYDEN.out';
    end

    %% ============================================================
    % 1. Read ZAERO Dryden data
    % =============================================================

    zaero = read_zaero_dryden_psd(filename);

    f = zaero.frequency;

    %% ============================================================
    % 2. Build MATLAB Dryden with exactly the same parameters
    % =============================================================

    matlabPsd = build_dryden_psd( ...
        f, ...
        zaero.Lg, ...
        zaero.sigmaW, ...
        zaero.V);

    %% ============================================================
    % 3. Numerical comparison
    % =============================================================

    absoluteError = abs(matlabPsd - zaero.psd);

    relativeError = ...
        absoluteError ./ max(abs(zaero.psd), eps);

    maxAbsoluteError = max(absoluteError);
    maxRelativeError = max(relativeError);
    meanRelativeError = mean(relativeError);

    %% ============================================================
    % 4. Compare finite-band RMS
    %
    % ZAERO's printed convention is recovered with:
    %
    % sigma^2 = 2 * integral_0^Fmax Phi(f) df
    % =============================================================

    rmsZaero = sqrt(2 * trapz(f, zaero.psd));

    rmsMatlab = sqrt(2 * trapz(f, matlabPsd));

    %% ============================================================
    % 5. Print summary
    % =============================================================

    fprintf('\n');
    fprintf('============================================================\n');
    fprintf('        ZAERO -- MATLAB DRYDEN PSD COMPARISON\n');
    fprintf('============================================================\n');

    fprintf('Source file              : %s\n', filename);
    fprintf('Velocity V               : %.6f m/s\n', zaero.V);
    fprintf('Scale length Lg          : %.6f m\n', zaero.Lg);
    fprintf('Input GURMS sigma_w      : %.6f m/s\n', zaero.sigmaW);
    fprintf('ZAERO reported GURMS     : %.6f m/s\n', ...
            zaero.computedGURMS);

    fprintf('Frequency spacing df     : %.6f Hz\n', zaero.df);
    fprintf('Number of frequencies    : %d\n', zaero.nFrequency);
    fprintf('Maximum frequency        : %.6f Hz\n', f(end));

    fprintf('\n');
    fprintf('Integrated RMS from ZAERO  = %.8f m/s\n', rmsZaero);
    fprintf('Integrated RMS from MATLAB = %.8f m/s\n', rmsMatlab);

    fprintf('\n');
    fprintf('Maximum absolute PSD error = %.6e\n', ...
            maxAbsoluteError);
    fprintf('Maximum relative error     = %.6e %%\n', ...
            100*maxRelativeError);
    fprintf('Mean relative error        = %.6e %%\n', ...
            100*meanRelativeError);

    fprintf('============================================================\n\n');

    %% ============================================================
    % 6. Figure
    % =============================================================

    positive = f > 0;

    figure( ...
        'Name', 'Dryden gust PSD comparison', ...
        'Color', 'white');

    % Plot MATLAB as continuous line
    semilogy( ...
        f(positive), ...
        matlabPsd(positive), ...
        '-', ...
        'LineWidth', 1.8);

    hold on;

    % Show a reduced number of ZAERO markers so the figure remains readable
    nPlotPoints = 200;

    skip = max(1, floor(sum(positive)/nPlotPoints));

    positiveIndices = find(positive);

    markerIndices = positiveIndices(1:skip:end);

    semilogy( ...
        f(markerIndices), ...
        zaero.psd(markerIndices), ...
        'o', ...
        'MarkerSize', 4, ...
        'LineWidth', 1.0);

    grid on;
    box on;

    xlabel('Frequency, f [Hz]');
    ylabel('Dryden gust PSD, \Phi_{w_g}');

    title(sprintf( ...
        ['Dryden PSD: ZAERO vs MATLAB\n' ...
         'V = %.1f m/s, L_g = %.3f m, \\sigma_w = %.3f m/s'], ...
        zaero.V, zaero.Lg, zaero.sigmaW));

    legend( ...
        'MATLAB analytical Dryden', ...
        'ZAERO', ...
        'Location', 'best');

    xlim([max(min(f(positive)),0.01), f(end)]);

    %% ============================================================
    % 7. Return results
    % =============================================================

    comparison = struct();

    comparison.frequency = f;

    comparison.zaeroPsd = zaero.psd;
    comparison.matlabPsd = matlabPsd;

    comparison.absoluteError = absoluteError;
    comparison.relativeError = relativeError;

    comparison.maxAbsoluteError = maxAbsoluteError;
    comparison.maxRelativeError = maxRelativeError;
    comparison.meanRelativeError = meanRelativeError;

    comparison.rmsZaero = rmsZaero;
    comparison.rmsMatlab = rmsMatlab;

    comparison.parameters = struct( ...
        'V', zaero.V, ...
        'Lg', zaero.Lg, ...
        'sigmaW', zaero.sigmaW, ...
        'zaeroComputedGURMS', zaero.computedGURMS, ...
        'df', zaero.df);

end
