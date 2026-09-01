function comparison = compare_dryden_frf(matlabResponse,zaeroFile)
%COMPARE_DRYDEN_FRF Align MATLAB and ZAERO gust frequency responses.

% The ASEOUT sensors are assumed to be modal-displacement sensors ordered
% as the retained modes.  FRF2 contains the ZAERO gust-input responses.

[V,rho,nf,f,~,ny,~,nuGust,FRF2] = frfrd(zaeroFile);
if nuGust < 1
    error('ZAERO file "%s" contains no gust-input FRF block.',zaeroFile);
end

nModes = min([matlabResponse.Nel,ny]);
zaeroModal = complex(zeros(nf,nModes));
for k = 1:nf
    if isempty(FRF2{k}) || size(FRF2{k},2) < 1
        error('Missing ZAERO gust FRF at frequency index %d.',k);
    end
    zaeroModal(k,:) = FRF2{k}(1:nModes,1).';
end

matlabModal = complex(zeros(nf,nModes));
for modeIndex = 1:nModes
    matlabModal(:,modeIndex) = interp1(matlabResponse.frequency, ...
        matlabResponse.modalX(:,modeIndex),f,'linear',NaN);
end

valid = all(isfinite(matlabModal),2);
relativeError = NaN(1,nModes);
for modeIndex = 1:nModes
    denominator = norm(zaeroModal(valid,modeIndex));
    if denominator > 0
        relativeError(modeIndex) = norm(matlabModal(valid,modeIndex)- ...
            zaeroModal(valid,modeIndex))/denominator;
    end
end

comparison = struct('frequency',f,'matlab',matlabModal, ...
    'zaero',zaeroModal,'relativeError',relativeError,'V',V,'rho',rho, ...
    'nModes',nModes,'nGustInputs',nuGust,'sourceFile',zaeroFile);
end
