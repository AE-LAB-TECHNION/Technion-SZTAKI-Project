function plot_dryden_results(cfg, result)
%PLOT_DRYDEN_RESULTS Plot time histories, spectra, and wind parameters.

plot_gust(cfg, result);
plot_modal_response(cfg, result);
plot_modal_psd(cfg, result);
plot_gust_psd(cfg, result);
plot_wind_parameters(cfg, result.windParameters);

end

function plot_gust(cfg, result)

fig = figure('Name', 'Dryden gust input');
tiledlayout(2,1,'TileSpacing','compact');
nexttile;
plot(result.time, result.gustPhysical(:,1), 'LineWidth', cfg.plotLineWidth);
grid on; box on; ylabel('w_g [m/s]');
nexttile;
plot(result.time, result.gustNormalized(:,1), 'LineWidth', cfg.plotLineWidth);
grid on; box on; xlabel('Time [s]'); ylabel('w_g/V');
apply_plot_style(fig, cfg);

end

function plot_modal_response(cfg, result)

fig = figure('Name', 'Dryden modal response');
tiledlayout(result.Nel,2,'TileSpacing','compact');
for k = 1:result.Nel
    nexttile;
    plot(result.time, result.modalX(:,k), 'LineWidth', cfg.plotLineWidth);
    grid on; box on; ylabel(sprintf('\\xi_%d',k));
    nexttile;
    plot(result.time, result.modalV(:,k), 'LineWidth', cfg.plotLineWidth);
    grid on; box on; ylabel(sprintf('d\\xi_%d/dt',k));
end
xlabel('Time [s]');
apply_plot_style(fig, cfg);

end

function plot_modal_psd(cfg, result)

frequency = logspace(-2,log10(1/(2*result.dt)),1200);
naug = size(result.A,1);
Cmodal = zeros(2*result.Nel,naug);
Cmodal(1:result.Nel,1:result.Nel) = eye(result.Nel);
Cmodal(result.Nel+(1:result.Nel),result.Nel+(1:result.Nel)) = eye(result.Nel);
response = squeeze(freqresp( ...
    ss(result.A,result.G,Cmodal,zeros(2*result.Nel,1)), 2*pi*frequency));
modalPsd = 2*abs(response).^2;

fig = figure('Name', 'Dryden modal-response PSD');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; hold on;
for k = 1:result.Nel
    loglog(frequency,modalPsd(k,:),'LineWidth',cfg.plotLineWidth);
end
grid on; box on; xlim([0.1 100]); xlabel('Frequency [Hz]');
ylabel('PSD of \xi_i [1/Hz]');
legend(compose('\\xi_%d',1:result.Nel),'Location','northeast');
nexttile; hold on;
for k = 1:result.Nel
    loglog(frequency,modalPsd(result.Nel+k,:),'LineWidth',cfg.plotLineWidth);
end
grid on; box on; xlim([0.1 100]); xlabel('Frequency [Hz]');
ylabel('PSD of d\xi_i/dt [1/s^2/Hz]');
legend(compose('d\\xi_%d/dt',1:result.Nel),'Location','northeast');
apply_plot_style(fig, cfg);

end

function plot_gust_psd(cfg, result)

sampleFrequency = 1/result.dt;
wg = result.gustPhysical(:,1);
[estimatedPsd, frequency] = pwelch(wg-mean(wg),[],[],[],sampleFrequency);
dryden = result.dryden;
response = squeeze(freqresp(ss(dryden.A,dryden.B,dryden.C(1,:),0), ...
    2*pi*frequency));
theoreticalPsd = dryden.inputNoiseScale^2*2*abs(response).^2;
fig = figure('Name', 'Dryden gust PSD');
loglog(frequency(2:end),estimatedPsd(2:end),'LineWidth',cfg.plotLineWidth);
hold on;
loglog(frequency(2:end),theoreticalPsd(2:end),'--','LineWidth',cfg.plotLineWidth);
grid on; box on; xlabel('Frequency [Hz]');
ylabel('PSD of w_g [(m/s)^2/Hz]');
legend({'Estimated','Dryden model'},'Location','best');
apply_plot_style(fig, cfg);

end

function plot_wind_parameters(cfg, wind)

maxAltitude = max([cfg.drydenWindProfileMaxAltitude, wind.altitude, ...
    wind.referenceHeight]);
altitude = linspace(max(1.01*wind.roughnessLength,0.05),maxAltitude,500);
meanWind = wind.u20*log(altitude/wind.roughnessLength) / ...
    log(wind.referenceHeight/wind.roughnessLength);

fig = figure('Name', 'MIL-STD-1797A wind parameters');
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
nexttile;
plot(meanWind,altitude,'LineWidth',cfg.plotLineWidth); hold on;
plot(wind.meanWindAtAltitude,wind.altitude,'o','MarkerSize',7,'MarkerFaceColor','auto');
grid on; box on; xlabel('u_w(h) [m/s]'); ylabel('Altitude h [m]');
legend({'Logarithmic profile','Operating point'},'Location','southeast');
nexttile;
plot(altitude/2,altitude,'LineWidth',cfg.plotLineWidth); hold on;
plot(wind.scaleLength,wind.altitude,'o','MarkerSize',7,'MarkerFaceColor','auto');
grid on; box on; xlabel('L_w=h/2 [m]'); ylabel('Altitude h [m]');
nexttile;
plot(wind.sigmaW*ones(size(altitude)),altitude,'LineWidth',cfg.plotLineWidth); hold on;
plot(wind.sigmaW,wind.altitude,'o','MarkerSize',7,'MarkerFaceColor','auto');
grid on; box on; xlabel('\sigma_w=0.1u_{20} [m/s]'); ylabel('Altitude h [m]');
apply_plot_style(fig, cfg);

end

