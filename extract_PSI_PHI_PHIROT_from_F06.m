function [PSI, PHI, PHI_ROT, info] = extract_PSI_PHI_PHIROT_from_F06(f06file, modesToUse, dispGridIDs, dispComp, rotGridIDs, rotComp, strainElemIDs)
%EXTRACT_PSI_PHI_PHIROT_FROM_F06
% Extract modal measurement matrices from a NASTRAN .f06 modal-analysis file.
%
% Outputs:
%   PSI      = strain mode matrix at strain-gauge / rod element locations
%              size: [numStrainElements x numModes]
%   PHI      = displacement mode matrix at selected GRID locations
%              size: [numDispSensors x numModes]
%   PHI_ROT  = rotational mode matrix at selected GRID locations
%              size: [numRotSensors x numModes]
%
% The .f06 is expected to contain:
%   1) REAL EIGEN VECTOR tables with columns:
%          POINT ID. TYPE T1 T2 T3 R1 R2 R3
%   2) CROD strain tables:
%          STRAINS IN ROD ELEMENTS (CROD)
%      with AXIAL STRAIN values.
%
% Important:
%   - For ZAERO after OMITMOD of mode 3, the first 3 modal coordinates usually
%     correspond to structural modes [1 2 4]. Therefore modesToUse=[1 2 4].
%   - PHI uses displacement components T1/T2/T3.
%   - PHI_ROT uses rotational components R1/R2/R3.
%   - PSI uses axial CROD strain.
%
% Example:
%   f06file = 'model-0012.f06';
%   modesToUse = [1 2 4];
%   dispGridIDs = [44 260 2691];
%   dispComp = 'T3';
%   rotGridIDs = [44 260 2691];
%   rotComp = 'R2';   % change according to the gyro axis you need
%   strainElemIDs = [13915:13924 13935:13944 17293:17296 17530:17897];
%   [PSI,PHI,PHI_ROT,info] = extract_PSI_PHI_PHIROT_from_F06( ...
%       f06file,modesToUse,dispGridIDs,dispComp,rotGridIDs,rotComp,strainElemIDs);
%   save('sensor_modal_matrices.mat','PSI','PHI','PHI_ROT','info');

    if nargin < 2 || isempty(modesToUse)
        modesToUse = [1 2 4];
    end
    if nargin < 3
        dispGridIDs = [];
    end
    if nargin < 4 || isempty(dispComp)
        dispComp = 'T3';
    end
    if nargin < 5
        rotGridIDs = dispGridIDs;
    end
    if nargin < 6 || isempty(rotComp)
        rotComp = 'R2';
    end
    if nargin < 7
        strainElemIDs = [];
    end

    compMap = struct('T1',1,'T2',2,'T3',3,'R1',4,'R2',5,'R3',6);
    dispCol = getComponentColumn(dispComp, compMap);
    rotCol  = getComponentColumn(rotComp,  compMap);

    fid = fopen(f06file,'r');
    if fid < 0
        error('Could not open file: %s', f06file);
    end
    cleaner = onCleanup(@() fclose(fid));

    % Raw parsed rows:
    % dispRows   = [mode, gridID, T1, T2, T3, R1, R2, R3]
    % strainRows = [mode, elemID, axialStrain, torsionalStrain]
    dispRows = [];
    strainRows = [];

    currentEig = NaN;
    currentMode = NaN;
    eigValsByMode = containers.Map('KeyType','double','ValueType','double');
    eigToMode = containers.Map('KeyType','char','ValueType','double');

    inDispTable = false;
    inStrainTable = false;

    while true
        line = fgetl(fid);
        if ~ischar(line)
            break;
        end

        % --- Eigenvalue line ---
        tokEig = regexp(line,'EIGENVALUE\s*=\s*([-+0-9.Ee]+)','tokens','once');
        if ~isempty(tokEig)
            currentEig = str2double(tokEig{1});
            inDispTable = false;
            inStrainTable = false;
            % If this eigenvalue was already associated with a mode, recover it.
            k = eigKey(currentEig);
            if isKey(eigToMode,k)
                currentMode = eigToMode(k);
            end
            continue;
        end

        % --- Real eigenvector number line ---
        tokMode = regexp(line,'R E A L\s+E I G E N\s+V E C T O R\s+N O \.\s*(\d+)','tokens','once');
        if ~isempty(tokMode)
            currentMode = str2double(tokMode{1});
            if ~isnan(currentEig)
                eigValsByMode(currentMode) = currentEig;
                eigToMode(eigKey(currentEig)) = currentMode;
            end
            inDispTable = false;
            inStrainTable = false;
            continue;
        end

        % --- Start of displacement table ---
        if contains(line,'POINT ID.') && contains(line,'T1') && contains(line,'R3')
            inDispTable = true;
            inStrainTable = false;
            continue;
        end

        % --- Start of CROD strain table ---
        if contains(line,'S T R A I N S') && contains(line,'C R O D')
            k = eigKey(currentEig);
            if isKey(eigToMode,k)
                currentMode = eigToMode(k);
            end
            inStrainTable = true;
            inDispTable = false;
            continue;
        end

        % --- Parse displacement row: GRID row with 6 modal components ---
        if inDispTable && ~isnan(currentMode)
            tok = regexp(line, '^\s*(\d+)\s+G\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)', 'tokens','once');
            if ~isempty(tok)
                gridID = str2double(tok{1});
                vals = cellfun(@str2double, tok(2:7));
                dispRows(end+1,:) = [currentMode, gridID, vals]; %#ok<AGROW>
            end
            continue;
        end

        % --- Parse CROD strain row.
        % In this file each element is printed as: elemID axialStrain torsionalStrain.
        % There may be two elements per line: [id axial torsion id axial torsion].
        if inStrainTable && ~isnan(currentMode)
            % Require a line that starts with an element ID, to avoid page headers.
            if isempty(regexp(line,'^\s*\d{4,}\s+[-+0-9.Ee]+','once'))
                continue;
            end
            nums = sscanf(line,'%f').';
            if numel(nums) >= 3
                nTriples = floor(numel(nums)/3);
                for j = 1:nTriples
                    base = 3*(j-1);
                    elemID = nums(base+1);
                    axial  = nums(base+2);
                    torsion = nums(base+3);
                    strainRows(end+1,:) = [currentMode, elemID, axial, torsion]; %#ok<AGROW>
                end
            end
            continue;
        end
    end

    if isempty(dispRows)
        error('No REAL EIGEN VECTOR displacement rows were parsed from the file.');
    end
    if isempty(strainRows)
        warning('No CROD strain rows were parsed from the file. PSI will be empty unless strain output exists.');
    end

    % Remove duplicates caused by page repetitions, keeping the last occurrence.
    dispRows = uniqueRowsByModeID(dispRows, 1, 2);
    strainRows = uniqueRowsByModeID(strainRows, 1, 2);

    if isempty(dispGridIDs)
        dispGridIDs = unique(dispRows(:,2)).';
    end
    if isempty(rotGridIDs)
        rotGridIDs = dispGridIDs;
    end
    if isempty(strainElemIDs) && ~isempty(strainRows)
        strainElemIDs = unique(strainRows(:,2)).';
    end

    nm = numel(modesToUse);
    PHI = nan(numel(dispGridIDs), nm);
    PHI_ROT = nan(numel(rotGridIDs), nm);
    PSI = nan(numel(strainElemIDs), nm);

    for im = 1:nm
        m = modesToUse(im);

        for i = 1:numel(dispGridIDs)
            idx = dispRows(:,1)==m & dispRows(:,2)==dispGridIDs(i);
            if any(idx)
                PHI(i,im) = dispRows(find(idx,1,'last'), 2 + dispCol);
            end
        end

        for i = 1:numel(rotGridIDs)
            idx = dispRows(:,1)==m & dispRows(:,2)==rotGridIDs(i);
            if any(idx)
                PHI_ROT(i,im) = dispRows(find(idx,1,'last'), 2 + rotCol);
            end
        end

        for i = 1:numel(strainElemIDs)
            idx = strainRows(:,1)==m & strainRows(:,2)==strainElemIDs(i);
            if any(idx)
                PSI(i,im) = strainRows(find(idx,1,'last'), 3); % axial strain
            end
        end
    end

    info = struct();
    info.f06file = f06file;
    info.modesToUse = modesToUse;
    info.dispGridIDs = dispGridIDs;
    info.dispComp = dispComp;
    info.rotGridIDs = rotGridIDs;
    info.rotComp = rotComp;
    info.strainElemIDs = strainElemIDs;
    info.availableModes_displacement = unique(dispRows(:,1)).';
    if ~isempty(strainRows)
        info.availableModes_strain = unique(strainRows(:,1)).';
    else
        info.availableModes_strain = [];
    end
    info.eigValsByMode = eigValsByMode;

    if any(isnan(PHI(:)))
        warning('PHI contains NaN values. Some requested GRID/mode/component entries were not found.');
    end
    if any(isnan(PHI_ROT(:)))
        warning('PHI_ROT contains NaN values. Some requested GRID/mode/component entries were not found.');
    end
    if ~isempty(PSI) && any(isnan(PSI(:)))
        warning('PSI contains NaN values. Some requested element/mode entries were not found.');
    end
end

% ========================================================================
function col = getComponentColumn(comp, compMap)
    if iscell(comp)
        error('This version expects one component string applied to all sensors, e.g. ''T3'' or ''R2''.');
    end
    comp = upper(strtrim(comp));
    if ~isfield(compMap, comp)
        error('Unknown component %s. Use T1,T2,T3,R1,R2,R3.', comp);
    end
    col = compMap.(comp);
end

% ========================================================================
function key = eigKey(eigVal)
    if isnan(eigVal)
        key = 'NaN';
    else
        key = sprintf('%.10E', eigVal);
    end
end

% ========================================================================
function out = uniqueRowsByModeID(rows, modeCol, idCol)
    if isempty(rows)
        out = rows;
        return;
    end
    keys = string(rows(:,modeCol)) + "_" + string(rows(:,idCol));
    [~, ia] = unique(keys, 'last');
    out = rows(sort(ia),:);
end
