function result = run_dryden_response(cfg, inputs, models)
%RUN_DRYDEN_RESPONSE Simulate the plant driven by Dryden turbulence.

result = struct();

[model, iV] = select_model_at_speed( ...
    cfg, models, cfg_value(cfg, 'gustValidationSpeed', 35));

if ~strcmpi(model.modelType, 'plant') || isempty(model.Ap) || isempty(model.Bpw)
    error(['Dryden simulation currently requires a plant model with gust ', ...
        'inputs (nonempty Ap and Bpw).']);
end

if size(model.Bpw,2) ~= 2
    error('Dryden simulation expects one gust vector and two Bpw columns.');
end

V = model.V;
h = cfg_value(cfg, 'drydenAltitude', 100);
referenceHeight = cfg_value(cfg, 'drydenReferenceHeight', 20*0.3048);
u20 = cfg_value(cfg, 'drydenMeanWindAt20ft', []);
if isempty(u20)
    % Backward compatibility with configurations created before the
    % MIL-STD-1797A 20-ft reference height was made explicit.
    u20 = cfg_value(cfg, 'drydenMeanWindAt6m', 7.0);
end
z0 = cfg_value(cfg, 'drydenRoughnessLength', 0.15*0.3048);
Lw = cfg_value(cfg, 'drydenScaleLength', []);
sigmaW = cfg_value(cfg, 'drydenSigmaW', []);
tau = cfg_value(cfg, 'drydenTau', 0.002);

if h <= z0 || referenceHeight <= z0 || z0 <= 0
    error('Dryden altitude and reference height must be greater than z_0.');
end

if isempty(Lw)
    Lw = h/2;
end

meanWindAtAltitude = u20*log(h/z0)/log(referenceHeight/z0);
if isempty(sigmaW)
    % MIL-STD-1797A low-altitude vertical turbulence intensity.
    sigmaW = 0.1*u20;
end

dryden = build_dryden_filter(V, Lw, sigmaW, tau);
CwNormalized = dryden.C/V;

np = size(model.Ap,1);
nw = size(dryden.A,1);

Aaug = [ ...
    model.Ap, model.Bpw*CwNormalized;
    zeros(nw,np), dryden.A];
Gaug = [zeros(np,1); dryden.B];
Gaug = dryden.inputNoiseScale*Gaug;

time = cfg_value(cfg, 'drydenTime', 0:0.001:60);
time = time(:);
if numel(time) < 2
    error('cfg.drydenTime must contain at least two samples.');
end

dt = diff(time);
if any(dt <= 0) || max(abs(dt-dt(1))) > 1e-9*max(1,dt(1))
    error('cfg.drydenTime must be strictly increasing and uniformly sampled.');
end
dt = dt(1);

seed = cfg_value(cfg, 'drydenRandomSeed', 1);
rng(seed, 'twister');

% Unit-intensity continuous white noise represented on a grid of width dt.
whiteNoise = randn(numel(time),1)/sqrt(dt);

% Simulate the cascade in two numerically scaled stages. This is
% input-output equivalent to the augmented model above, while avoiding
% ill-conditioning caused by mixing plant and companion-filter states.
sysDryden = ss( ...
    dryden.A, dryden.inputNoiseScale*dryden.B, ...
    dryden.C, dryden.D);
sysDrydenDiscrete = c2d(sysDryden, dt, 'zoh');
gustPhysical = lsim(sysDrydenDiscrete, whiteNoise);
gustNormalized = gustPhysical/V;

sysPlant = ss( ...
    model.Ap, model.Bpw, eye(np), zeros(np,size(model.Bpw,2)));
sysPlantDiscrete = c2d(sysPlant, dt, 'zoh');
xp = lsim(sysPlantDiscrete, gustNormalized);

Nel = cfg.Nel;
modalX = xp(:,1:Nel);
modalV = xp(:,Nel+(1:Nel));

Cp = [model.Cae, model.Dae];
if size(Cp,2) ~= np
    error('Plant sensor-output matrix has %d columns; expected %d.', size(Cp,2), np);
end
sensorOutput = xp*Cp.' + gustNormalized*model.Caw.';

steadyIdx = time >= time(1) + 0.5*(time(end)-time(1));
realizedRmsWg = rms(gustPhysical(steadyIdx,1));

fprintf('\n============================================================\n');
fprintf('Dryden turbulence simulation\n');
fprintf('============================================================\n');
fprintf('Selected V = %.2f m/s, index %d\n', V, iV);
fprintf('Altitude h = %.3f m\n', h);
fprintf('Reference height = %.3f m (20 ft)\n', referenceHeight);
fprintf('Selected u_20 = %.6f m/s\n', u20);
fprintf('Category-C roughness z_0 = %.6f m (0.15 ft)\n', z0);
fprintf('Mean wind u_w(h) = %.6f m/s\n', meanWindAtAltitude);
fprintf('Scale length Lw = %.3f m\n', Lw);
fprintf('Configured sigma_w = %.6f m/s\n', sigmaW);
fprintf('Filter theoretical RMS(w_g) = %.6f m/s\n', dryden.expectedRmsWg);
fprintf('Realized second-half RMS(w_g) = %.6f m/s\n', realizedRmsWg);
fprintf('Continuous white-noise scale = %.6g\n', dryden.inputNoiseScale);
fprintf('LPF tau = %.6g s, random seed = %d\n', tau, seed);
fprintf('Augmented state size = %d\n', size(Aaug,1));

plot_dryden_gust(cfg, time, gustPhysical, gustNormalized);
plot_dryden_modal_response(cfg, time, modalX, modalV);
plot_dryden_modal_psd(cfg, Aaug, Gaug, Nel, dt);
plot_dryden_psd(cfg, time, gustPhysical(:,1), dryden);
plot_dryden_wind_parameters( ...
    cfg, h, referenceHeight, u20, z0, meanWindAtAltitude, Lw, sigmaW);

result.time = time;
result.whiteNoise = whiteNoise;
result.gustPhysical = gustPhysical;
result.gustNormalized = gustNormalized;
result.modalX = modalX;
result.modalV = modalV;
result.sensorOutput = sensorOutput;
result.A = Aaug;
result.G = Gaug;
result.dryden = dryden;
result.realizedRmsWg = realizedRmsWg;
result.windParameters = struct( ...
    'altitude', h, ...
    'referenceHeight', referenceHeight, ...
    'u20', u20, ...
    'roughnessLength', z0, ...
    'meanWindAtAltitude', meanWindAtAltitude, ...
    'scaleLength', Lw, ...
    'sigmaW', sigmaW);

end

% =========================================================================
function plot_dryden_gust(cfg, time, gustPhysical, gustNormalized)

fig = figure('Name', 'Dryden gust input');
apply_plot_style(fig, cfg);
tiledlayout(2,1,'TileSpacing','compact');

nexttile;
plot(time, gustPhysical(:,1), 'LineWidth', cfg_value(cfg,'plotLineWidth',1.5));
grid on; box on;
ylabel('w_g [m/s]');

nexttile;
plot(time, gustNormalized(:,1), 'LineWidth', cfg_value(cfg,'plotLineWidth',1.5));
grid on; box on;
xlabel('Time [s]');
ylabel('w_g/V');

apply_plot_style(fig, cfg);
exportgraphics(fig, 'dryden_gust_input.png', 'Resolution', 300);

end

% =========================================================================
function plot_dryden_modal_response(cfg, time, modalX, modalV)

Nel = size(modalX,2);
fig = figure('Name', 'Dryden modal response');
apply_plot_style(fig, cfg);
tiledlayout(Nel,2,'TileSpacing','compact');

for k = 1:Nel
    nexttile;
    plot(time, modalX(:,k), 'LineWidth', cfg_value(cfg,'plotLineWidth',1.5));
    grid on; box on;
    ylabel(sprintf('\\xi_%d',k));

    nexttile;
    plot(time, modalV(:,k), 'LineWidth', cfg_value(cfg,'plotLineWidth',1.5));
    grid on; box on;
    ylabel(sprintf('d\\xi_%d/dt',k));
end

xlabel('Time [s]');
apply_plot_style(fig, cfg);
exportgraphics(fig, 'dryden_modal_response.png', 'Resolution', 300);

end

% =========================================================================
function plot_dryden_modal_psd(cfg, Aaug, Gaug, Nel, dt)

sampleFrequency = 1/dt;
frequency = logspace(-2,log10(sampleFrequency/2),1200);
naug = size(Aaug,1);

% Outputs are the retained modal displacements and velocities in the plant
% portion of the augmented state vector.
Cmodal = zeros(2*Nel,naug);
Cmodal(1:Nel,1:Nel) = eye(Nel);
Cmodal(Nel+(1:Nel),Nel+(1:Nel)) = eye(Nel);

modalResponse = squeeze(freqresp( ...
    ss(Aaug,Gaug,Cmodal,zeros(2*Nel,1)), 2*pi*frequency));
modalPsd = 2*abs(modalResponse).^2;

fig = figure('Name', 'Dryden modal-response PSD');
apply_plot_style(fig, cfg);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile;
hold on;
for k = 1:Nel
    loglog(frequency, modalPsd(k,:), ...
        'LineWidth', cfg_value(cfg,'plotLineWidth',1.5));
end
grid on; box on;
xlim([0.1 100]);
xlabel('Frequency [Hz]');
ylabel('PSD of \xi_i [1/Hz]');
legend(compose('\\xi_%d',1:Nel),'Location','northeast');

nexttile;
hold on;
for k = 1:Nel
    loglog(frequency, modalPsd(Nel+k,:), ...
        'LineWidth', cfg_value(cfg,'plotLineWidth',1.5));
end
grid on; box on;
xlim([0.1 100]);
xlabel('Frequency [Hz]');
ylabel('PSD of d\xi_i/dt [1/s^2/Hz]');
legend(compose('d\\xi_%d/dt',1:Nel),'Location','northeast');

apply_plot_style(fig, cfg);
exportgraphics(fig, 'dryden_modal_psd.png', 'Resolution', 300);

end

% =========================================================================
function plot_dryden_psd(cfg, time, wg, dryden)

dt = time(2)-time(1);
sampleFrequency = 1/dt;
[estimatedPsd, frequency] = pwelch( ...
    wg-mean(wg), [], [], [], sampleFrequency);

response = squeeze(freqresp( ...
    ss(dryden.A,dryden.B,dryden.C(1,:),0), 2*pi*frequency));
theoreticalPsd = 2*abs(response).^2;
theoreticalPsd = dryden.inputNoiseScale^2*theoreticalPsd;

fig = figure('Name', 'Dryden gust PSD');
apply_plot_style(fig, cfg);
loglog(frequency(2:end), estimatedPsd(2:end), ...
    'LineWidth', cfg_value(cfg,'plotLineWidth',1.5));
hold on;
loglog(frequency(2:end), theoreticalPsd(2:end), '--', ...
    'LineWidth', cfg_value(cfg,'plotLineWidth',1.5));
grid on; box on;
xlabel('Frequency [Hz]');
ylabel('PSD of w_g [(m/s)^2/Hz]');
legend({'Estimated','Dryden model'},'Location','best');
apply_plot_style(fig, cfg);
exportgraphics(fig, 'dryden_gust_psd.png', 'Resolution', 300);

end

% =========================================================================
function plot_dryden_wind_parameters( ...
    cfg, h, referenceHeight, u20, z0, meanWindAtAltitude, Lw, sigmaW)

maxAltitude = cfg_value(cfg, 'drydenWindProfileMaxAltitude', max(30,2*h));
maxAltitude = max([maxAltitude, h, referenceHeight]);
altitude = linspace(max(1.01*z0,0.05), maxAltitude, 500);
meanWind = u20*log(altitude/z0)/log(referenceHeight/z0);
verticalScaleLength = altitude/2;
verticalTurbulenceRms = sigmaW*ones(size(altitude));

fig = figure('Name', 'MIL-STD-1797A wind parameters');
apply_plot_style(fig, cfg);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

nexttile;
plot(meanWind, altitude, 'LineWidth', cfg_value(cfg,'plotLineWidth',1.5));
hold on;
plot(meanWindAtAltitude, h, 'o', 'MarkerSize', 7, ...
    'MarkerFaceColor', 'auto');
grid on; box on;
xlabel('u_w(h) [m/s]');
ylabel('Altitude h [m]');
legend({'Logarithmic profile','Operating point'}, ...
    'Location','southeast', ...
    'FontSize', max(8,cfg_value(cfg,'plotFontSize',13)-3));

nexttile;
plot(verticalScaleLength, altitude, ...
    'LineWidth', cfg_value(cfg,'plotLineWidth',1.5));
hold on;
plot(Lw, h, 'o', 'MarkerSize', 7, 'MarkerFaceColor', 'auto');
grid on; box on;
xlabel('L_w=h/2 [m]');
ylabel('Altitude h [m]');

nexttile;
plot(verticalTurbulenceRms, altitude, ...
    'LineWidth', cfg_value(cfg,'plotLineWidth',1.5));
hold on;
plot(sigmaW, h, 'o', 'MarkerSize', 7, 'MarkerFaceColor', 'auto');
grid on; box on;
xlabel('\sigma_w=0.1u_{20} [m/s]');
ylabel('Altitude h [m]');

apply_plot_style(fig, cfg);
exportgraphics(fig, 'dryden_wind_parameters.png', 'Resolution', 300);

end

% =========================================================================
function [model, iV] = select_model_at_speed(cfg, models, requestedSpeed)

[~, iV] = min(abs(cfg.airspeed-requestedSpeed));
model = models{iV};

end

% =========================================================================
function value = cfg_value(cfg, fieldName, defaultValue)

if isfield(cfg,fieldName)
    value = cfg.(fieldName);
else
    value = defaultValue;
end

end
