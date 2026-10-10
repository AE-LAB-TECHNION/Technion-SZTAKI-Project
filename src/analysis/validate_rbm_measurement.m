function validation = validate_rbm_measurement(inputs, models, augmented)
%VALIDATE_RBM_MEASUREMENT Verify RBM output mapping through all model levels.

if nargin < 3
    augmented = [];
end

if ~isfield(inputs,'RBM') || isempty(inputs.RBM)
    error('inputs.RBM is missing. The RBM measurement is required.');
end
if isempty(models)
    error('No aeroelastic models were supplied for RBM validation.');
end

modalCoefficients = inputs.RBM.modalCoefficients;
Nel = size(modalCoefficients,2);
nRbm = size(modalCoefficients,1);
tol = 1e-12;

maxCaeError = 0;
maxDaeError = 0;
maxCawError = 0;
maxCpError = 0;

for modelIndex = 1:numel(models)
    model = models{modelIndex};

    if ~isfield(model,'outputIndex') || ~isfield(model.outputIndex,'rbm')
        error('Model at V=%.9g m/s does not define outputIndex.rbm.',model.V);
    end

    rbmRows = model.outputIndex.rbm;
    if numel(rbmRows) ~= nRbm
        error('Expected %d RBM output row(s), but the model defines %d.', ...
            nRbm,numel(rbmRows));
    end

    expectedCae = zeros(nRbm,size(model.Cae,2));
    expectedCae(:,1:Nel) = modalCoefficients;

    maxCaeError = max(maxCaeError, ...
        norm(model.Cae(rbmRows,:)-expectedCae,inf));
    maxDaeError = max(maxDaeError, ...
        norm(model.Dae(rbmRows,:),inf));
    maxCawError = max(maxCawError, ...
        norm(model.Caw(rbmRows,:),inf));

    if strcmpi(model.modelType,'plant')
        expectedCp = [expectedCae, zeros(nRbm,size(model.Dae,2))];
    else
        expectedCp = expectedCae;
    end
    maxCpError = max(maxCpError,norm(model.Cp(rbmRows,:)-expectedCp,inf));
end

if max([maxCaeError,maxDaeError,maxCawError,maxCpError]) > tol
    error(['RBM measurement validation failed: Cae=%.3e, Dae=%.3e, ', ...
        'Caw=%.3e, Cp=%.3e.'], ...
        maxCaeError,maxDaeError,maxCawError,maxCpError);
end

augmentedCError = NaN;
augmentedDError = NaN;
augmentedSize = [0,0];

if ~isempty(augmented)
    model = models{augmented.modelIndex};
    rbmRows = model.outputIndex.rbm;
    nAddedStates = size(augmented.C,2)-size(model.Cp,2);

    if nAddedStates < 0
        error('The augmented model has fewer state columns than the plant output matrix.');
    end

    expectedAugmented = [ ...
        model.Cp(rbmRows,:), ...
        zeros(nRbm,nAddedStates) ];

    augmentedCError = norm( ...
        augmented.C(rbmRows,:)-expectedAugmented,inf);
    augmentedDError = norm(augmented.D(rbmRows,:),inf);
    augmentedSize = size(augmented.C);

    if max(augmentedCError,augmentedDError) > tol
        error(['Augmented RBM measurement validation failed: ', ...
            'C=%.3e, D=%.3e.'],augmentedCError,augmentedDError);
    end
end

modelSpeeds = cellfun(@(model) model.V,models);
[~,sourceModelIndex] = min(abs(modelSpeeds-inputs.RBM.V));
sourceModel = models{sourceModelIndex};
rbmRows = sourceModel.outputIndex.rbm;

validation = struct();
validation.passed = true;
validation.tolerance = tol;
validation.outputRows = rbmRows;
validation.modalCoefficients = modalCoefficients;
validation.maxCaeError = maxCaeError;
validation.maxDaeError = maxDaeError;
validation.maxCawError = maxCawError;
validation.maxCpError = maxCpError;
validation.augmentedCError = augmentedCError;
validation.augmentedDError = augmentedDError;
validation.CaeSize = size(sourceModel.Cae);
validation.CpSize = size(sourceModel.Cp);
validation.augmentedCSize = augmentedSize;

fprintf('\n============================================================\n');
fprintf('RBM MEASUREMENT VALIDATION\n');
fprintf('============================================================\n');
fprintf('ZAERO source              = %s, output row %d\n', ...
    inputs.RBM.sourceFile,inputs.RBM.outputRow);
fprintf('MATLAB measurement row    =');
fprintf(' %d',rbmRows);
fprintf('\n');
fprintf('modal coefficients        =');
fprintf(' %.12e',modalCoefficients);
fprintf('\n');
fprintf('Cae size                  = %d x %d\n',validation.CaeSize);
fprintf('Cp size                   = %d x %d\n',validation.CpSize);
if ~isempty(augmented)
    fprintf('Caug size                 = %d x %d\n',augmentedSize);
end
fprintf('max Cae mapping error     = %.3e\n',maxCaeError);
fprintf('max Dae RBM-row norm      = %.3e\n',maxDaeError);
fprintf('max Caw RBM-row norm      = %.3e\n',maxCawError);
fprintf('max Cp mapping error      = %.3e\n',maxCpError);
if ~isempty(augmented)
    fprintf('Caug mapping error        = %.3e\n',augmentedCError);
    fprintf('Daug RBM-row norm         = %.3e\n',augmentedDError);
end
fprintf('status                    = PASS\n');
fprintf('============================================================\n\n');

end
