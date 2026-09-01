function response = compute_dryden_frequency_response(cfg,inputs,models,qhgFile,frequency)
%COMPUTE_DRYDEN_FREQUENCY_RESPONSE MATLAB modal response using ZAERO QHG.

arguments
    cfg struct
    inputs struct
    models cell
    qhgFile (1,:) char = 'QHG.DAT'
    frequency (:,1) double = (0:0.01:100).'
end

[~,iV] = min(abs(cfg.airspeed-cfg.gustValidationSpeed));
model = models{iV};
V = model.V;
Nel = cfg.Nel;
Qhg = Read_out4_to_mat_form(qhgFile,'QHGSAV');
if size(Qhg,1) ~= Nel && size(Qhg,2) == Nel
    Qhg = Qhg.';
end
if size(Qhg,1) ~= Nel || size(Qhg,2) ~= numel(frequency)
    error('QHG dimensions are %d x %d; expected %d x %d.', ...
        size(Qhg,1),size(Qhg,2),Nel,numel(frequency));
end

qdyn = 0.5*cfg.rho*V^2;
A2h = inputs.RFA_mat.A2(1:Nel,1:Nel);
Mbar = eye(Nel)-qdyn*(cfg.L/V)^2*A2h;
nLag = size(model.Aae,1)-2*Nel;
Bq = [zeros(Nel); Mbar\(qdyn*eye(Nel)); zeros(nLag,Nel)];

modalX = complex(zeros(numel(frequency),Nel));
modalV = complex(zeros(numel(frequency),Nel));
stateIdentity = eye(size(model.Aae));
for j = 1:numel(frequency)
    x = ((1i*2*pi*frequency(j))*stateIdentity-model.Aae) ...
        \ (Bq*(Qhg(:,j)/V));
    modalX(j,:) = x(1:Nel).';
    modalV(j,:) = x(Nel+(1:Nel)).';
end

response = struct('frequency',frequency,'modalX',modalX,'modalV',modalV, ...
    'modalXPsd',abs(modalX).^2,'modalVPsd',abs(modalV).^2, ...
    'V',V,'Nel',Nel);
end
