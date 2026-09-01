function comparison = plot_dryden_qhg_comparison(cfg, inputs, models, qhgFile, dryFile)
%PLOT_DRYDEN_QHG_COMPARISON Compare the MATLAB RFA model with ZAERO GLOADS.
%
% QHG.DAT contains the frequency-dependent generalized gust forces used by
% ZAERO.  Driving the MATLAB structural RFA with the same forces makes this
% an apples-to-apples comparison; it does not require gust columns in the
% APPROX file.

arguments
    cfg struct
    inputs struct
    models cell
    qhgFile char = 'QHG.DAT'
    dryFile char = 'DRY.PLT'
end

[~, iV] = min(abs(cfg.airspeed-cfg.gustValidationSpeed));
model = models{iV};
V = model.V;
Nel = cfg.Nel;

zaero = load_zaero_dryden_response(dryFile);
frequency = zaero.frequency(:);
Qhg = Read_out4_to_mat_form(qhgFile, 'QHGSAV');
if size(Qhg,1) ~= Nel && size(Qhg,2) == Nel
    Qhg = Qhg.';
end
if size(Qhg,1) ~= Nel || size(Qhg,2) ~= numel(frequency)
    error(['QHG dimensions are %d x %d; expected %d modal rows and ', ...
        '%d frequency columns.'], size(Qhg,1), size(Qhg,2), ...
        Nel, numel(frequency));
end

% The saved QHG is the generalized gust-force spectrum used by GLOADS.
% ZAERO defines the gust velocity input in normalized form, w_g/V.
qdyn = 0.5*cfg.rho*V^2;
A2h = inputs.RFA_mat.A2(1:Nel,1:Nel);
Mbar = eye(Nel) - qdyn*(cfg.L/V)^2*A2h;
nLag = size(model.Aae,1)-2*Nel;
Bq = [zeros(Nel); Mbar\(qdyn*eye(Nel)); zeros(nLag,Nel)];

modalX = complex(zeros(numel(frequency),Nel));
modalV = complex(zeros(numel(frequency),Nel));
stateIdentity = eye(size(model.Aae));
for j = 1:numel(frequency)
    stateResponse = ((1i*2*pi*frequency(j))*stateIdentity-model.Aae) ...
        \ (Bq*(Qhg(:,j)/V));
    modalX(j,:) = stateResponse(1:Nel).';
    modalV(j,:) = stateResponse(Nel+(1:Nel)).';
end

matlabXPsd = abs(modalX).^2;
matlabVPsd = abs(modalV).^2;
positive = frequency > 0;

fig = figure('Name','MATLAB--ZAERO Dryden comparison using QHG');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
colors = lines(2*Nel);
matlabColors = colors(1:Nel,:);
zaeroColors = colors(Nel+(1:Nel),:);
nexttile; hold on;
displacementHandles = gobjects(1,2*Nel);
displacementLabels = strings(1,2*Nel);
for k = 1:Nel
    displacementHandles(2*k-1) = loglog( ...
        frequency(positive),matlabXPsd(positive,k),'-', ...
        'Color',matlabColors(k,:),'LineWidth',cfg.plotLineWidth);
    displacementHandles(2*k) = loglog( ...
        frequency(positive),zaero.modalXPsd(positive,k),'--', ...
        'Color',zaeroColors(k,:),'LineWidth',1.25*cfg.plotLineWidth);
    displacementLabels(2*k-1) = sprintf('Mode %d - MATLAB',k);
    displacementLabels(2*k) = sprintf('Mode %d - ZAERO',k);
end
grid on; box on; xlim([0.1 100]); xlabel('Frequency [Hz]');
ylabel('PSD of \xi_i [1/Hz]'); title('Modal displacement');
legend(displacementHandles,displacementLabels, ...
    'Location','northeast','NumColumns',2);

nexttile; hold on;
velocityHandles = gobjects(1,2*Nel);
velocityLabels = strings(1,2*Nel);
for k = 1:Nel
    velocityHandles(2*k-1) = loglog( ...
        frequency(positive),matlabVPsd(positive,k),'-', ...
        'Color',matlabColors(k,:),'LineWidth',cfg.plotLineWidth);
    velocityHandles(2*k) = loglog( ...
        frequency(positive),zaero.modalVPsd(positive,k),'--', ...
        'Color',zaeroColors(k,:),'LineWidth',1.25*cfg.plotLineWidth);
    velocityLabels(2*k-1) = sprintf('Mode %d - MATLAB',k);
    velocityLabels(2*k) = sprintf('Mode %d - ZAERO',k);
end
grid on; box on; xlim([0.1 100]); xlabel('Frequency [Hz]');
ylabel('PSD of d\xi_i/dt [1/s^2/Hz]'); title('Modal velocity');
legend(velocityHandles,velocityLabels, ...
    'Location','northeast','NumColumns',2);
apply_plot_style(fig,cfg);

rmsMatlab = [sqrt(trapz(frequency,matlabXPsd)), ...
             sqrt(trapz(frequency,matlabVPsd))].';
rmsZaero = [sqrt(trapz(frequency,zaero.modalXPsd)), ...
            sqrt(trapz(frequency,zaero.modalVPsd))].';
labels = [compose('xi_%d',1:Nel),compose('dxi_%d/dt',1:Nel)].';
comparison = table(labels,rmsMatlab,rmsZaero, ...
    100*(rmsMatlab-rmsZaero)./rmsZaero, ...
    'VariableNames',{'Response','RmsMatlab','RmsZaero', ...
    'RelativeDifferencePercent'});
fprintf('\nMATLAB--ZAERO Dryden comparison using QHG at %.2f m/s\n',V);
disp(comparison);
end
