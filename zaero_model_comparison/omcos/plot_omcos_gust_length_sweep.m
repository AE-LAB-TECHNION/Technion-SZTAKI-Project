function result = plot_omcos_gust_length_sweep( ...
    cfg,models,gustLengths,timeEnd)
%PLOT_OMCOS_GUST_LENGTH_SWEEP Compare modal responses for OMCOS gust lengths.

comparisonDirectory = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(comparisonDirectory));
addpath(projectRoot);
addpath(genpath(fullfile(projectRoot,'src')));

if nargin < 1 || isempty(cfg)
    cfg = default_analysis_config();
end
if nargin < 2
    models = [];
end
if nargin < 3 || isempty(gustLengths)
    gustLengths = [1.18 2.5 7 35];
end
if nargin < 4 || isempty(timeEnd)
    timeEnd = 5;
end

gustLengths = gustLengths(:).';
if any(~isfinite(gustLengths)) || any(gustLengths <= 0)
    error('gustLengths must contain positive finite lengths in meters.');
end
if ~isscalar(timeEnd) || ~isfinite(timeEnd) || timeEnd <= 0
    error('timeEnd must be a positive finite scalar.');
end

if isempty(models)
    [inputs,cfg] = load_analysis_inputs(cfg);
    cfg.airspeed = cfg.gustValidationSpeed;
    [models,~,~] = run_airspeed_sweep(cfg,inputs);
end

[model,speedIndex] = select_model(models,cfg.gustValidationSpeed);
nModes = cfg.Nel;
nStates = size(model.Aae,1);
nInputs = size(model.Baw,2);

if nInputs < 2 || mod(nInputs,2) ~= 0
    error('Baw must contain [w_g/V, dw_g/dt/V] input pairs.');
end

timeStep = cfg.gustInternalDt;
time = (0:timeStep:timeEnd).';
if time(end) < timeEnd
    time(end+1,1) = timeEnd;
end

modalOutput = zeros(2*nModes,nStates);
modalOutput(1:nModes,1:nModes) = eye(nModes);
modalOutput(nModes+(1:nModes),nModes+(1:nModes)) = eye(nModes);
system = ss(model.Aae,model.Baw,modalOutput,zeros(2*nModes,nInputs));

nLengths = numel(gustLengths);
modalX = zeros(numel(time),nModes,nLengths);
modalV = zeros(numel(time),nModes,nLengths);
gustDuration = gustLengths/model.V;

for lengthIndex = 1:nLengths
    gustInput = build_gust_input( ...
        cfg,time,model.V,gustDuration(lengthIndex),nInputs);
    output = lsim(system,gustInput,time);
    modalX(:,:,lengthIndex) = output(:,1:nModes);
    modalV(:,:,lengthIndex) = output(:,nModes+(1:nModes));
end

fig = plot_responses( ...
    cfg,time,modalX,modalV,gustLengths,timeEnd);

result.figure = fig;
result.time = time;
result.speed = model.V;
result.speedIndex = speedIndex;
result.gustLengths = gustLengths;
result.gustDuration = gustDuration;
result.modalX = modalX;
result.modalV = modalV;

end

% =========================================================================
function gustInput = build_gust_input(cfg,time,V,Tg,nInputs)

nGust = nInputs/2;
eta = time-(cfg.gustXResponse-cfg.gustX0)/V;
active = eta >= 0 & eta <= Tg;

wbar = zeros(size(time));
wdotbar = zeros(size(time));
wbar(active) = 0.5*cfg.gustWGV*(1-cos(2*pi*eta(active)/Tg));
wdotbar(active) = pi*cfg.gustWGV/Tg*sin(2*pi*eta(active)/Tg);

gustInput = zeros(numel(time),nInputs);
gustInput(:,1:nGust) = repmat(wbar,1,nGust);
gustInput(:,nGust+1:end) = repmat(wdotbar,1,nGust);

end

% =========================================================================
function fig = plot_responses( ...
    cfg,time,modalX,modalV,gustLengths,timeEnd)

nModes = size(modalX,2);
nLengths = numel(gustLengths);
colors = lines(nLengths);
referenceLengths = [1.18 2.5 7 35];
referenceColors = [ ...
    0.55 0.82 0.95; ... % light blue
    0.95 0.58 0.10; ... % gold
    0.00 0.45 0.32; ... % dark green
    0.55 0.20 0.65];    % purple
for colorIndex = 1:numel(referenceLengths)
    lengthIndex = find(abs(gustLengths-referenceLengths(colorIndex)) < 1e-9,1);
    if ~isempty(lengthIndex)
        colors(lengthIndex,:) = referenceColors(colorIndex,:);
    end
end
lineStyles = {'-','--','-.',':'};
lineWidth = cfg.plotLineWidth;

fig = figure( ...
    'Name','OMCOS modal response for different gust lengths', ...
    'Units','normalized', ...
    'OuterPosition',[0.05 0.08 0.90 0.82]);
layout = tiledlayout(nModes,2,'TileSpacing','compact','Padding','compact');
legendText = arrayfun( ...
    @(Lg) sprintf('%.3g m',Lg),gustLengths,'UniformOutput',false);

for modeIndex = 1:nModes
    nexttile(2*modeIndex-1);
    responseHandles = plot_response_family( ...
        time,squeeze(modalX(:,modeIndex,:)),colors,lineStyles,lineWidth);
    ylabel(sprintf('\\xi_%d',modeIndex));
    xlim([0 timeEnd]);

    if modeIndex == 1
        legend(responseHandles,legendText,'Location','best');
    end

    nexttile(2*modeIndex);
    plot_response_family( ...
        time,squeeze(modalV(:,modeIndex,:)),colors,lineStyles,lineWidth);
    ylabel(sprintf('d\\xi_%d/dt',modeIndex));
    xlim([0 timeEnd]);
end

nexttile(2*nModes-1);
xlabel('Time [s]');
nexttile(2*nModes);
xlabel('Time [s]');

layout.TileSpacing = 'compact';
apply_plot_style(fig,cfg);

end

% =========================================================================
function handles = plot_response_family( ...
    time,response,colors,lineStyles,lineWidth)

nLengths = size(response,2);
handles = gobjects(nLengths,1);
hold on;

for lengthIndex = 1:nLengths
    styleIndex = mod(lengthIndex-1,numel(lineStyles))+1;
    handles(lengthIndex) = plot(time,response(:,lengthIndex), ...
        'Color',colors(lengthIndex,:), ...
        'LineStyle',lineStyles{styleIndex}, ...
        'LineWidth',lineWidth);
end

grid on;
box on;

end

% =========================================================================
function [model,index] = select_model(models,requestedSpeed)

if ~iscell(models)
    models = num2cell(models);
end

speeds = cellfun(@(item) item.V,models);
[~,index] = min(abs(speeds-requestedSpeed));
model = models{index};

end
