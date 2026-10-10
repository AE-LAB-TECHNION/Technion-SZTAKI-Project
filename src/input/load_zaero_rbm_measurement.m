function rbm = load_zaero_rbm_measurement(filename, outputRow, Nel)
%LOAD_ZAERO_RBM_MEASUREMENT Read the modal RBM output row from ZAERO ASEOUT.

if isstring(filename)
    filename = char(filename);
end

if exist(filename,'file') ~= 2
    error('ZAERO RBM model file was not found: %s',filename);
end

[V,rho,nx,nu,ny,~,~,~,C,D] = sysrd(filename);

if outputRow < 1 || outputRow ~= round(outputRow) || outputRow > ny
    error(['Requested RBM output row %g, but %s contains %d output rows.'], ...
        outputRow,filename,ny);
end

if nx < Nel
    error('%s contains %d states, but Nel=%d modal coefficients are required.', ...
        filename,nx,Nel);
end

modalCoefficients = C(outputRow,1:Nel);
nonModalPart = C(outputRow,Nel+1:end);
directFeedthrough = D(outputRow,:);
tol = 1e-12;

if norm(nonModalPart,inf) > tol
    error(['RBM output row %d in %s depends on states beyond the first ', ...
        '%d modal-displacement states.'],outputRow,filename,Nel);
end

if norm(directFeedthrough,inf) > tol
    error('RBM output row %d in %s has nonzero direct feedthrough.', ...
        outputRow,filename);
end

rbm = struct();
rbm.modalCoefficients = modalCoefficients;
rbm.sourceCRow = C(outputRow,:);
rbm.sourceDRow = directFeedthrough;
rbm.outputRow = outputRow;
rbm.sourceFile = filename;
rbm.V = V;
rbm.rho = rho;
rbm.nx = nx;
rbm.nu = nu;
rbm.ny = ny;

fprintf('\nZAERO RBM measurement source: %s\n',filename);
fprintf('  ASEOUT dimensions : nx=%d, nu=%d, ny=%d\n',nx,nu,ny);
fprintf('  RBM output row    : %d\n',outputRow);
fprintf('  modal coefficients:');
fprintf(' %.12e',modalCoefficients);
fprintf('\n');

end
