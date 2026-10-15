function comparison = plot_omcos_rfa_comparison(cfg,models)
%PLOT_OMCOS_RFA_COMPARISON Compare old and partitioned RFA responses.

if nargin < 1 || isempty(cfg)
    cfg = default_analysis_config();
end

if nargin < 2
    models = [];
end

comparisonDirectory = fileparts(mfilename('fullpath'));
reference = read_modal_response(fullfile(comparisonDirectory,'GUSTFREQ'));
oldRfa = read_modal_response(fullfile(comparisonDirectory,'GUSTXI'));
partitionedRfa = read_modal_response(fullfile(comparisonDirectory,'GUSTP20'));

validate_time_grids(reference,oldRfa,partitionedRfa);
[matlabLsim,cfg] = build_matlab_lsim_response(cfg,reference.time,models);

comparison.errorTable = build_error_table( ...
    reference,oldRfa,partitionedRfa,matlabLsim);
disp(comparison.errorTable);

comparison.responses = struct( ...
    'frequencyReference',reference, ...
    'oldRfa',oldRfa, ...
    'partitionedRfa',partitionedRfa, ...
    'matlabLsim',matlabLsim);

comparison.figures.modalResponses = plot_modal_responses( ...
    cfg,reference,oldRfa,partitionedRfa,matlabLsim);
comparison.figures.representativeErrors = plot_representative_errors( ...
    cfg,reference,oldRfa,partitionedRfa,matlabLsim);
comparison.figures.modalErrors = plot_modal_errors( ...
    cfg,reference,oldRfa,partitionedRfa);

end

% =========================================================================
function fig = plot_modal_errors( ...
    cfg,reference,oldRfa,partitionedRfa)

nModes = min(3,size(reference.modalX,2));
timeLimit = min(5,reference.time(end));
lineWidth = cfg_value(cfg,'plotLineWidth',1.5);

fig = figure( ...
    'Name','OMCOS RFA errors relative to FREQ', ...
    'Units','normalized', ...
    'OuterPosition',[0.03 0.08 0.94 0.82]);
tiledlayout(2,nModes,'TileSpacing','compact','Padding','compact');

for modeIndex = 1:nModes
    nexttile(modeIndex);
    plot_error_set(reference.time,reference.modalX(:,modeIndex), ...
        oldRfa.modalX(:,modeIndex),partitionedRfa.modalX(:,modeIndex), ...
        lineWidth);
    title(sprintf('MODALX error, mode %d',modeIndex));
    xlabel('Time [s]');
    ylabel('\Delta\xi');
    xlim([reference.time(1),timeLimit]);

    if modeIndex == 1
        legend({'FREQ','Old RFA','Partitioned RFA'}, ...
            'Location','best');
    end

    nexttile(nModes+modeIndex);
    plot_error_set(reference.time,reference.modalV(:,modeIndex), ...
        oldRfa.modalV(:,modeIndex),partitionedRfa.modalV(:,modeIndex), ...
        lineWidth);
    title(sprintf('MODALV error, mode %d',modeIndex));
    xlabel('Time [s]');
    ylabel('\Delta(d\xi/dt)');
    xlim([reference.time(1),timeLimit]);
end

apply_plot_style(fig,cfg);

end

% =========================================================================
function plot_error_set( ...
    time,reference,oldRfa,partitionedRfa,lineWidth)

plot(time,zeros(size(reference)),'k-','LineWidth',lineWidth);
hold on;
plot(time,oldRfa-reference,'--', ...
    'Color',[0.85 0.33 0.10],'LineWidth',lineWidth);
plot(time,partitionedRfa-reference,'-', ...
    'Color',[0 0.45 0.74],'LineWidth',lineWidth);
grid on;
box on;

end

% =========================================================================
function fig = plot_modal_responses( ...
    cfg,reference,oldRfa,partitionedRfa,matlabLsim)

nModes = size(reference.modalX,2);
timeLimit = min(4,reference.time(end));
lineWidth = cfg_value(cfg,'plotLineWidth',1.5);

fig = figure( ...
    'Name','OMCOS RFA comparison: modal responses', ...
    'Units','normalized', ...
    'OuterPosition',[0.05 0.08 0.90 0.82]);
tiledlayout(nModes,2,'TileSpacing','compact','Padding','compact');

for modeIndex = 1:nModes
    nexttile(2*modeIndex-1);
    plot_response_set(reference.time,reference.modalX(:,modeIndex), ...
        oldRfa.modalX(:,modeIndex),partitionedRfa.modalX(:,modeIndex), ...
        matlabLsim.modalX(:,modeIndex),lineWidth);
    ylabel(sprintf('\\xi_%d',modeIndex));
    xlim([reference.time(1),timeLimit]);

    if modeIndex == 1
        responseLegend = legend( ...
            {'Frequency-domain reference','Old RFA','Partitioned RFA', ...
            'MATLAB lsim'}, ...
            'Orientation','horizontal');
        responseLegend.Layout.Tile = 'north';
    end

    nexttile(2*modeIndex);
    plot_response_set(reference.time,reference.modalV(:,modeIndex), ...
        oldRfa.modalV(:,modeIndex),partitionedRfa.modalV(:,modeIndex), ...
        matlabLsim.modalV(:,modeIndex),lineWidth);
    ylabel(sprintf('d\\xi_%d/dt',modeIndex));
    xlim([reference.time(1),timeLimit]);
end

nexttile(2*nModes-1);
xlabel('Time [s]');
nexttile(2*nModes);
xlabel('Time [s]');
apply_plot_style(fig,cfg);

end

% =========================================================================
function plot_response_set( ...
    time,reference,oldRfa,partitionedRfa,matlabLsim,lineWidth)

plot(time,reference,'k-','LineWidth',lineWidth);
hold on;
plot(time,oldRfa,'--','Color',[0.85 0.33 0.10],'LineWidth',lineWidth);
plot(time,partitionedRfa,'-','Color',[0 0.45 0.74],'LineWidth',lineWidth);
plot(time,matlabLsim,'-.','Color',[0.49 0.18 0.56],'LineWidth',lineWidth);
grid on;
box on;

end

% =========================================================================
function fig = plot_representative_errors( ...
    cfg,reference,oldRfa,partitionedRfa,matlabLsim)

modeIndex = min(2,size(reference.modalX,2));
timeLimit = min(4,reference.time(end));
lineWidth = cfg_value(cfg,'plotLineWidth',1.5);

fig = figure( ...
    'Name','OMCOS RFA comparison: normalized errors', ...
    'Units','normalized', ...
    'OuterPosition',[0.10 0.10 0.80 0.75]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

nexttile;
plot_normalized_errors(reference.time,reference.modalX(:,modeIndex), ...
    oldRfa.modalX(:,modeIndex),partitionedRfa.modalX(:,modeIndex), ...
    matlabLsim.modalX(:,modeIndex),lineWidth);
ylabel(sprintf('e_{\\xi_%d} [%%]',modeIndex));
xlim([reference.time(1),timeLimit]);
legend({'Old RFA','Partitioned RFA','MATLAB lsim'},'Location','best');

nexttile;
plot_normalized_errors(reference.time,reference.modalV(:,modeIndex), ...
    oldRfa.modalV(:,modeIndex),partitionedRfa.modalV(:,modeIndex), ...
    matlabLsim.modalV(:,modeIndex),lineWidth);
ylabel(sprintf('e_{d\\xi_%d/dt} [%%]',modeIndex));
xlabel('Time [s]');
xlim([reference.time(1),timeLimit]);
apply_plot_style(fig,cfg);

end

% =========================================================================
function plot_normalized_errors( ...
    time,reference,oldRfa,partitionedRfa,matlabLsim,lineWidth)

normalization = max(abs(reference));
oldError = 100*(oldRfa-reference)/normalization;
partitionedError = 100*(partitionedRfa-reference)/normalization;
matlabError = 100*(matlabLsim-reference)/normalization;

plot(time,oldError,'--','Color',[0.85 0.33 0.10],'LineWidth',lineWidth);
hold on;
plot(time,partitionedError,'-','Color',[0 0.45 0.74],'LineWidth',lineWidth);
plot(time,matlabError,'-.','Color',[0.49 0.18 0.56],'LineWidth',lineWidth);
yline(0,':','Color',[0.45 0.45 0.45],'LineWidth',1.0);
grid on;
box on;

end

% =========================================================================
function errorTable = build_error_table( ...
    reference,oldRfa,partitionedRfa,matlabLsim)

referenceData = [reference.modalX,reference.modalV];
oldData = [oldRfa.modalX,oldRfa.modalV];
partitionedData = [partitionedRfa.modalX,partitionedRfa.modalV];
matlabData = [matlabLsim.modalX,matlabLsim.modalV];
nModes = size(reference.modalX,2);

response = [compose('X%d',(1:nModes).'); compose('V%d',(1:nModes).')];
peakReference = max(abs(referenceData),[],1);
peakErrorOld = 100*abs(max(abs(oldData),[],1)./peakReference-1);
peakErrorPartitioned = 100*abs(max(abs(partitionedData),[],1)./peakReference-1);
peakErrorMatlab = 100*abs(max(abs(matlabData),[],1)./peakReference-1);
responseRange = max(referenceData,[],1)-min(referenceData,[],1);
nrmseOld = 100*sqrt(mean((oldData-referenceData).^2,1))./responseRange;
nrmsePartitioned = 100*sqrt(mean((partitionedData-referenceData).^2,1))./responseRange;
nrmseMatlab = 100*sqrt(mean((matlabData-referenceData).^2,1))./responseRange;

errorTable = table(response,peakErrorOld.',peakErrorPartitioned.', ...
    peakErrorMatlab.',nrmseOld.',nrmsePartitioned.',nrmseMatlab.', ...
    'VariableNames',{'Response','PeakErrorOldPct','PeakErrorPartitionedPct', ...
    'PeakErrorMatlabPct','NRMSEOldPct','NRMSEPartitionedPct','NRMSEMatlabPct'});

end

% =========================================================================
function [response,cfg] = build_matlab_lsim_response(cfg,time,models)

if isempty(models)
    [inputs,cfg] = load_analysis_inputs(cfg);
    cfg.airspeed = cfg.gustValidationSpeed;
    [models,~,~] = run_airspeed_sweep(cfg,inputs);
    model = models{1};
else
    [~,speedIndex] = min(abs(cfg.airspeed-cfg.gustValidationSpeed));
    model = models{speedIndex};
end

nModes = cfg.Nel;
nStates = size(model.Aae,1);
nInputs = size(model.Baw,2);

if nInputs < 2 || mod(nInputs,2) ~= 0
    error('Baw must contain [w_g/V, dw_g/dt/V] input pairs.');
end

nGust = nInputs/2;
time = time(:);
eta = time-(cfg.gustXResponse-cfg.gustX0)/model.V;

switch lower(cfg.gustLengthMode)
    case 'time'
        gustDuration = cfg.gustLG;
    case 'length'
        gustDuration = cfg.gustLG/model.V;
    otherwise
        error('gustLengthMode must be ''time'' or ''length''.');
end

active = eta >= 0 & eta <= gustDuration;
wbar = zeros(size(time));
wdotbar = zeros(size(time));
wbar(active) = 0.5*cfg.gustWGV* ...
    (1-cos(2*pi*eta(active)/gustDuration));
wdotbar(active) = pi*cfg.gustWGV/gustDuration* ...
    sin(2*pi*eta(active)/gustDuration);

gustInput = zeros(numel(time),nInputs);
gustInput(:,1:nGust) = repmat(wbar,1,nGust);
gustInput(:,nGust+1:end) = repmat(wdotbar,1,nGust);

modalOutput = zeros(2*nModes,nStates);
modalOutput(1:nModes,1:nModes) = eye(nModes);
modalOutput(nModes+(1:nModes),nModes+(1:nModes)) = eye(nModes);
system = ss(model.Aae,model.Baw,modalOutput,zeros(2*nModes,nInputs));
output = lsim(system,gustInput,time);

response.time = time;
response.modalX = output(:,1:nModes);
response.modalV = output(:,nModes+(1:nModes));

end

% =========================================================================
function response = read_modal_response(filename)

fid = fopen(filename,'r');

if fid < 0
    error('Could not open ZAERO response file: %s',filename);
end

cleanup = onCleanup(@() fclose(fid));
rows = zeros(0,0);

while true
    line = fgetl(fid);

    if ~ischar(line)
        break
    end

    values = sscanf(line,'%f').';

    if numel(values) >= 3 && mod(numel(values)-1,2) == 0
        if isempty(rows)
            rows = zeros(0,numel(values));
        end

        if numel(values) == size(rows,2)
            rows(end+1,:) = values; %#ok<AGROW>
        end
    end
end

if isempty(rows)
    error('No modal response data found in %s.',filename);
end

nModes = (size(rows,2)-1)/2;
response.time = rows(:,1);
response.modalV = rows(:,1+(1:nModes));
response.modalX = rows(:,1+nModes+(1:nModes));

end

% =========================================================================
function validate_time_grids(reference,oldRfa,partitionedRfa)

if ~isequal(reference.time,oldRfa.time) || ...
        ~isequal(reference.time,partitionedRfa.time)
    error('The three ZAERO responses must use the same time grid.');
end

if size(reference.modalX,2) ~= size(oldRfa.modalX,2) || ...
        size(reference.modalX,2) ~= size(partitionedRfa.modalX,2)
    error('The three ZAERO responses must contain the same number of modes.');
end

end

% =========================================================================
function value = cfg_value(cfg,fieldName,defaultValue)

if isfield(cfg,fieldName) && ~isempty(cfg.(fieldName))
    value = cfg.(fieldName);
else
    value = defaultValue;
end

end
