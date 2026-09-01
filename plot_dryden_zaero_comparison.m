function comparison = plot_dryden_zaero_comparison(cfg,result,filename)
%PLOT_DRYDEN_ZAERO_COMPARISON Compare MATLAB and ZAERO modal Dryden PSDs.

zaero = load_zaero_dryden_response(filename);
frequency = zaero.frequency(:);
positive = frequency > 0;

naug = size(result.A,1);
Cmodal = zeros(2*result.Nel,naug);
Cmodal(1:result.Nel,1:result.Nel) = eye(result.Nel);
Cmodal(result.Nel+(1:result.Nel),result.Nel+(1:result.Nel)) = eye(result.Nel);
response = squeeze(freqresp( ...
    ss(result.A,result.G,Cmodal,zeros(2*result.Nel,1)),2*pi*frequency));
matlabPsd = 2*abs(response).^2;
matlabXPsd = matlabPsd(1:result.Nel,:).';
matlabVPsd = matlabPsd(result.Nel+(1:result.Nel),:).';

rmsMatlab = [sqrt(trapz(frequency,matlabXPsd)), ...
             sqrt(trapz(frequency,matlabVPsd))].';
rmsZaero = [sqrt(trapz(frequency,zaero.modalXPsd)), ...
            sqrt(trapz(frequency,zaero.modalVPsd))].';
allMatlabPsd = [matlabXPsd,matlabVPsd];
allZaeroPsd = [zaero.modalXPsd,zaero.modalVPsd];
[~,matlabPeakIndex] = max(allMatlabPsd,[],1);
[~,zaeroPeakIndex] = max(allZaeroPsd,[],1);
peakFrequencyMatlab = frequency(matlabPeakIndex(:));
peakFrequencyZaero = frequency(zaeroPeakIndex(:));

matlabColor = [0.0000 0.4470 0.7410];
zaeroColor = [0.8500 0.3250 0.0980];
fig = figure('Name','MATLAB--ZAERO Dryden modal PSD comparison', ...
    'Color','white');
layout = tiledlayout(result.Nel,2,'TileSpacing','compact','Padding','compact');
title(layout,sprintf('Dryden modal PSD comparison at V = %.1f m/s', ...
    cfg.gustValidationSpeed));

firstAxis = [];
for k = 1:result.Nel
    axX = nexttile; hold(axX,'on');
    if isempty(firstAxis), firstAxis = axX; end
    loglog(axX,frequency(positive),matlabXPsd(positive,k),'-', ...
        'Color',matlabColor,'LineWidth',cfg.plotLineWidth+0.3);
    loglog(axX,frequency(positive),zaero.modalXPsd(positive,k),'--', ...
        'Color',zaeroColor,'LineWidth',cfg.plotLineWidth+0.3);
    format_axis(axX,k,'Modal displacement PSD [1/Hz]', ...
        peakFrequencyMatlab(k),peakFrequencyZaero(k), ...
        rmsMatlab(k),rmsZaero(k));

    axV = nexttile; hold(axV,'on');
    loglog(axV,frequency(positive),matlabVPsd(positive,k),'-', ...
        'Color',matlabColor,'LineWidth',cfg.plotLineWidth+0.3);
    loglog(axV,frequency(positive),zaero.modalVPsd(positive,k),'--', ...
        'Color',zaeroColor,'LineWidth',cfg.plotLineWidth+0.3);
    j = result.Nel+k;
    format_axis(axV,k,'Modal velocity PSD [1/s^2/Hz]', ...
        peakFrequencyMatlab(j),peakFrequencyZaero(j), ...
        rmsMatlab(j),rmsZaero(j));
end
legend(firstAxis,{'MATLAB','ZAERO'},'Location','southwest');
apply_plot_style(fig,cfg);

labels = [compose('xi_%d',1:result.Nel), ...
          compose('dxi_%d/dt',1:result.Nel)].';
comparison = table(labels,rmsMatlab,rmsZaero, ...
    100*(rmsMatlab-rmsZaero)./rmsZaero, ...
    peakFrequencyMatlab,peakFrequencyZaero, ...
    'VariableNames',{'Response','RmsMatlab','RmsZaero', ...
    'RelativeDifferencePercent','PeakFrequencyMatlabHz', ...
    'PeakFrequencyZaeroHz'});
fprintf('\nMATLAB--ZAERO Dryden PSD comparison at %.2f m/s\n', ...
    cfg.gustValidationSpeed);
disp(comparison);
end

function format_axis(ax,modeNumber,yLabel,peakMatlab,peakZaero,rmsMatlab,rmsZaero)
grid(ax,'on'); box(ax,'on');
xlim(ax,[0.1 100]);
xlabel(ax,'Frequency [Hz]');
ylabel(ax,yLabel);
title(ax,sprintf(['Mode %d: f_p = %.2f / %.2f Hz, ', ...
    'RMS ratio = %.3f'],modeNumber,peakMatlab,peakZaero, ...
    rmsMatlab/rmsZaero));
end
