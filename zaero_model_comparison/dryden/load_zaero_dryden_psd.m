function zaero = load_zaero_dryden_psd(filename)
%LOAD_ZAERO_DRYDEN_PSD Read modal-response PSDs from ZAERO DRYPSD.P.
%
% The current ZAERO output is expected to contain six auto-PSD CROSPSD
% results:
%
%   101 -> xi_1
%   102 -> xi_2
%   103 -> xi_3
%   104 -> dxi_1/dt
%   105 -> dxi_2/dt
%   106 -> dxi_3/dt
%
% ZAERO may print several CROSPSD IDs in one table and the remaining IDs in
% subsequent tables. This reader collects all CROSPSD tables in the file.
%
% Output:
%   zaero.frequency
%   zaero.modalXPsd
%   zaero.modalVPsd
%   zaero.displacementIds
%   zaero.velocityIds
%   zaero.rmsModalX
%   zaero.rmsModalV
%
% RMS is evaluated using ZAERO's positive-frequency convention:
%
%   RMS = sqrt(2 * integral_0^Fmax PSD(f) df)

arguments
    filename (1,:) char = 'DRYPSD.P'
end

if ~isfile(filename)
    error('ZAERO PSD file not found: %s',filename);
end

txt = fileread(filename);
lines = regexp(txt,'\r\n|\n|\r','split');

allIds = [];
allFrequency = {};
allRealPsd = {};

i = 1;

while i <= numel(lines)

    if contains(lines{i},'FREQ,HZ') && contains(lines{i},'CROSPSD')

        if i == numel(lines)
            break
        end

        % The next line contains repeated IDs, e.g.
        % 101) 101) 102) 102) ...
        idTokens = regexp(lines{i+1},'(\d+)\)','tokens');

        if isempty(idTokens)
            i = i + 1;
            continue
        end

        repeatedIds = cellfun(@(c) str2double(c{1}),idTokens);
        ids = unique(repeatedIds,'stable');
        nIds = numel(ids);

        j = i + 2;
        maxBlockRows = numel(lines)-j+1;
        fBlock = zeros(maxBlockRows,1);
        realBlock = zeros(maxBlockRows,nIds);
        blockRowCount = 0;

        while j <= numel(lines)

            line = strtrim(lines{j});

            if isempty(line)
                if blockRowCount > 0
                    break
                end
                j = j + 1;
                continue
            end

            if startsWith(line,'MAX/MIN') || ...
               startsWith(line,'AT FREQ') || ...
               startsWith(line,'R M S') || ...
               startsWith(line,'CYCLES')
                break
            end

            vals = sscanf(lines{j},'%f').';

            % Each auto/cross PSD has real and imaginary columns:
            % f, Re1, Im1, Re2, Im2, ...
            if numel(vals) == 1 + 2*nIds
                blockRowCount = blockRowCount + 1;
                fBlock(blockRowCount) = vals(1);

                row = zeros(1,nIds);
                for k = 1:nIds
                    row(k) = vals(2*k);   % real part
                end

                realBlock(blockRowCount,:) = row;
            elseif blockRowCount > 0
                break
            end

            j = j + 1;
        end

        fBlock = fBlock(1:blockRowCount);
        realBlock = realBlock(1:blockRowCount,:);

        if blockRowCount > 0
            for k = 1:nIds
                allIds(end+1) = ids(k); %#ok<AGROW>
                allFrequency{end+1} = fBlock; %#ok<AGROW>
                allRealPsd{end+1} = realBlock(:,k); %#ok<AGROW>
            end
        end

        i = max(i+1,j);

    else
        i = i + 1;
    end
end

if isempty(allIds)
    error('No CROSPSD tables were found in "%s".',filename);
end

%% ------------------------------------------------------------------------
% Map the six CROSPSD IDs used in the current ZAERO Dryden run
% -------------------------------------------------------------------------
displacementIds = 101:103;
velocityIds = 104:106;

requiredIds = [displacementIds velocityIds];

missing = requiredIds(~ismember(requiredIds,allIds));

if ~isempty(missing)
    error('Missing CROSPSD IDs in "%s": %s', ...
          filename,mat2str(missing));
end

% Use CROSPSD 101 as the master frequency grid.
masterIndex = find(allIds == displacementIds(1),1,'first');
frequency = allFrequency{masterIndex};

nFreq = numel(frequency);
nModes = numel(displacementIds);

modalXPsd = zeros(nFreq,nModes);
modalVPsd = zeros(nFreq,nModes);

for k = 1:nModes

    idxX = find(allIds == displacementIds(k),1,'first');
    idxV = find(allIds == velocityIds(k),1,'first');

    fX = allFrequency{idxX};
    fV = allFrequency{idxV};

    psdX = allRealPsd{idxX};
    psdV = allRealPsd{idxV};

    if numel(fX) ~= nFreq || max(abs(fX-frequency)) > 1e-10
        psdX = interp1(fX,psdX,frequency,'linear',NaN);
    end

    if numel(fV) ~= nFreq || max(abs(fV-frequency)) > 1e-10
        psdV = interp1(fV,psdV,frequency,'linear',NaN);
    end

    modalXPsd(:,k) = psdX;
    modalVPsd(:,k) = psdV;
end

if any(~isfinite(modalXPsd),'all') || any(~isfinite(modalVPsd),'all')
    error('Interpolation of ZAERO modal PSDs produced non-finite values.');
end

% Auto-PSD should be real and nonnegative. Tiny negative roundoff, if any,
% is clipped to zero.
modalXPsd = max(modalXPsd,0);
modalVPsd = max(modalVPsd,0);

rmsModalX = sqrt(2*trapz(frequency,modalXPsd,1));
rmsModalV = sqrt(2*trapz(frequency,modalVPsd,1));

zaero = struct();

zaero.frequency = frequency;
zaero.modalXPsd = modalXPsd;
zaero.modalVPsd = modalVPsd;

zaero.displacementIds = displacementIds;
zaero.velocityIds = velocityIds;

zaero.rmsModalX = rmsModalX;
zaero.rmsModalV = rmsModalV;

zaero.filename = filename;

fprintf('\nRead ZAERO modal PSD file: %s\n',filename);
for k = 1:nModes
    fprintf(['Mode %d: CROSPSD %d RMS(xi)=%.6e, ', ...
             'CROSPSD %d RMS(dxi/dt)=%.6e\n'], ...
            k,displacementIds(k),rmsModalX(k), ...
            velocityIds(k),rmsModalV(k));
end

end
