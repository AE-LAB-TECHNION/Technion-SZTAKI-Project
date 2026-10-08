function comparison = plot_dryden_zaero_comparison(cfg,response)
%PLOT_DRYDEN_ZAERO_COMPARISON Compare MATLAB and ZAERO modal Dryden PSDs.
%
% This plots the same physical quantities on both sides:
%
%   left column  - modal displacement PSD
%   right column - modal velocity PSD
%
% one row per retained mode.
%
% MATLAB PSD:
%   computed analytically in compute_dryden_frequency_response.m
%
% ZAERO PSD:
%   read from CROSPSD 101:106 in DRYPSD.P
%
% Usage:
%
%   response = compute_dryden_frequency_response(cfg,inputs,models);
%
%   comparison = plot_dryden_zaero_comparison(cfg,response);

arguments
    cfg struct
    response struct
end

zaero = load_zaero_dryden_psd(cfg.drydenZaeroPsdFile);

f = zaero.frequency(:);
nModes = min([response.Nel, size(zaero.modalXPsd,2), ...
              size(zaero.modalVPsd,2)]);

%% ------------------------------------------------------------------------
% Interpolate MATLAB calculation to the exact ZAERO frequency grid
% -------------------------------------------------------------------------
matlabXPsd = zeros(numel(f),nModes);
matlabVPsd = zeros(numel(f),nModes);

for k = 1:nModes
    matlabXPsd(:,k) = interp1( ...
        response.frequency,response.modalXPsd(:,k),f,'linear',NaN);

    matlabVPsd(:,k) = interp1( ...
        response.frequency,response.modalVPsd(:,k),f,'linear',NaN);
end

valid = isfinite(f) & ...
        all(isfinite(matlabXPsd),2) & ...
        all(isfinite(matlabVPsd),2);

f = f(valid);

matlabXPsd = matlabXPsd(valid,:);
matlabVPsd = matlabVPsd(valid,:);

zaeroXPsd = zaero.modalXPsd(valid,1:nModes);
zaeroVPsd = zaero.modalVPsd(valid,1:nModes);

%% ------------------------------------------------------------------------
% RMS comparison
%
% ZAERO prints the positive-frequency half spectrum, therefore:
%
%       RMS = sqrt(2*integral PSD(f) df)
% -------------------------------------------------------------------------
rmsMatlabX = sqrt(2*trapz(f,matlabXPsd,1));
rmsZaeroX  = sqrt(2*trapz(f,zaeroXPsd,1));

rmsMatlabV = sqrt(2*trapz(f,matlabVPsd,1));
rmsZaeroV  = sqrt(2*trapz(f,zaeroVPsd,1));

%% ------------------------------------------------------------------------
% Peak-frequency comparison
% -------------------------------------------------------------------------
[~,idxMx] = max(matlabXPsd,[],1);
[~,idxZx] = max(zaeroXPsd,[],1);

[~,idxMv] = max(matlabVPsd,[],1);
[~,idxZv] = max(zaeroVPsd,[],1);

peakMatlabX = f(idxMx);
peakZaeroX = f(idxZx);

peakMatlabV = f(idxMv);
peakZaeroV = f(idxZv);

%% ------------------------------------------------------------------------
% Figure: one row per mode; displacement left, velocity right
% -------------------------------------------------------------------------
positive = f > 0;

fig = figure( ...
    'Name','MATLAB--ZAERO Dryden modal PSD comparison', ...
    'Color','white');

layout = tiledlayout( ...
    nModes,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');

title(layout,sprintf( ...
    ['Dryden modal-response PSD comparison: ', ...
     'V = %.1f m/s, L_g = %.3f m, \\sigma_w = %.3f m/s'], ...
    response.Vzaero,response.Lg,response.sigmaW));

firstAxis = [];

for k = 1:nModes

    % -------------------------------------------------------------
    % Modal displacement PSD
    % -------------------------------------------------------------
    axX = nexttile;
    hold(axX,'on');

    if isempty(firstAxis)
        firstAxis = axX;
    end

    loglog(axX, ...
        f(positive),matlabXPsd(positive,k),'-', ...
        'LineWidth',cfg.plotLineWidth+0.2);

    loglog(axX, ...
        f(positive),zaeroXPsd(positive,k),'--', ...
        'LineWidth',cfg.plotLineWidth+0.2);

    grid(axX,'on');
    box(axX,'on');

    xlim(axX,[0.1 min(100,f(end))]);

    xlabel(axX,'Frequency [Hz]');
    ylabel(axX,sprintf('S_{\\xi_%d}(f)',k));

    title(axX,sprintf( ...
        ['Mode %d displacement: RMS %.4e / %.4e, ', ...
         'f_p %.2f / %.2f Hz'], ...
        k,rmsMatlabX(k),rmsZaeroX(k), ...
        peakMatlabX(k),peakZaeroX(k)));

    % -------------------------------------------------------------
    % Modal velocity PSD
    % -------------------------------------------------------------
    axV = nexttile;
    hold(axV,'on');

    loglog(axV, ...
        f(positive),matlabVPsd(positive,k),'-', ...
        'LineWidth',cfg.plotLineWidth+0.2);

    loglog(axV, ...
        f(positive),zaeroVPsd(positive,k),'--', ...
        'LineWidth',cfg.plotLineWidth+0.2);

    grid(axV,'on');
    box(axV,'on');

    xlim(axV,[0.1 min(100,f(end))]);

    xlabel(axV,'Frequency [Hz]');
    ylabel(axV,sprintf('S_{d\\xi_%d/dt}(f)',k));

    title(axV,sprintf( ...
        ['Mode %d velocity: RMS %.4e / %.4e, ', ...
         'f_p %.2f / %.2f Hz'], ...
        k,rmsMatlabV(k),rmsZaeroV(k), ...
        peakMatlabV(k),peakZaeroV(k)));
end

legend(firstAxis,{'MATLAB','ZAERO'},'Location','best');

if exist('apply_plot_style','file') == 2
    apply_plot_style(fig,cfg);
end

%% ------------------------------------------------------------------------
% Comparison table
% -------------------------------------------------------------------------
responseName = strings(2*nModes,1);

rmsMatlab = zeros(2*nModes,1);
rmsZaero = zeros(2*nModes,1);

peakMatlab = zeros(2*nModes,1);
peakZaero = zeros(2*nModes,1);

for k = 1:nModes

    responseName(k) = sprintf('xi_%d',k);
    rmsMatlab(k) = rmsMatlabX(k);
    rmsZaero(k) = rmsZaeroX(k);
    peakMatlab(k) = peakMatlabX(k);
    peakZaero(k) = peakZaeroX(k);

    j = nModes+k;

    responseName(j) = sprintf('dxi_%d/dt',k);
    rmsMatlab(j) = rmsMatlabV(k);
    rmsZaero(j) = rmsZaeroV(k);
    peakMatlab(j) = peakMatlabV(k);
    peakZaero(j) = peakZaeroV(k);
end

relativeRmsDifferencePercent = ...
    100*(rmsMatlab-rmsZaero)./rmsZaero;

comparison = table( ...
    responseName, ...
    rmsMatlab, ...
    rmsZaero, ...
    relativeRmsDifferencePercent, ...
    peakMatlab, ...
    peakZaero, ...
    'VariableNames',{ ...
        'Response', ...
        'RmsMatlab', ...
        'RmsZaero', ...
        'RelativeRmsDifferencePercent', ...
        'PeakFrequencyMatlabHz', ...
        'PeakFrequencyZaeroHz'});

fprintf('\n============================================================\n');
fprintf('MATLAB--ZAERO DRYDEN MODAL PSD COMPARISON\n');
fprintf('============================================================\n');
disp(comparison);

end
