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

    if nargin < 2 || isempty(modesToUse)
        modesToUse = [1 2 4];
    end
    if nargin < 3
        dispGridIDs = [];
    end
    if nargin < 4 || isempty(dispComp)
        dispComp = 'T3';
    end
    if nargin < 5 || isempty(rotGridIDs)
        rotGridIDs = dispGridIDs;
    end
    if nargin < 6 || isempty(rotComp)
        rotComp = 'R2';
    end
    if nargin < 7
        strainElemIDs = [];
    end

    modesToUse = modesToUse(:).';
    dispGridIDs = dispGridIDs(:).';
    rotGridIDs = rotGridIDs(:).';
    strainElemIDs = strainElemIDs(:).';

    compMap = struct('T1',1,'T2',2,'T3',3,'R1',4,'R2',5,'R3',6);
    dispCol = getComponentColumn(dispComp, compMap);
    rotCol  = getComponentColumn(rotComp,  compMap);

    fid = fopen(f06file,'r');
    if fid < 0
        error('Could not open file: %s', f06file);
    end
    cleaner = onCleanup(@() fclose(fid));

    nm = numel(modesToUse);
    collectDispRows = isempty(dispGridIDs) || isempty(rotGridIDs);
    collectStrainRows = isempty(strainElemIDs);

    if collectDispRows
        dispRows = [];
        PHI = [];
        PHI_ROT = [];
    else
        PHI = nan(numel(dispGridIDs), nm);
        PHI_ROT = nan(numel(rotGridIDs), nm);
    end

    if collectStrainRows
        strainRows = [];
        PSI = [];
    else
        PSI = nan(numel(strainElemIDs), nm);
    end

    currentEig = NaN;
    currentMode = NaN;
    eigValsByMode = containers.Map('KeyType','double','ValueType','double');
    eigToMode = containers.Map('KeyType','char','ValueType','double');

    inDispTable = false;
    inStrainTable = false;
    foundDispRows = false;
    foundStrainRows = false;
    availableModesDisp = [];
    availableModesStrain = [];

    while true
        line = fgetl(fid);
        if ~ischar(line)
            break;
        end

        if contains(line,'EIGENVALUE')
            tokEig = regexp(line,'EIGENVALUE\s*=\s*([-+0-9.Ee]+)','tokens','once');
            if ~isempty(tokEig)
                currentEig = str2double(tokEig{1});
                currentMode = NaN;
                inDispTable = false;
                inStrainTable = false;

                k = eigKey(currentEig);
                if isKey(eigToMode,k)
                    currentMode = eigToMode(k);
                end
                continue;
            end
        end

        if contains(line,'R E A L') && contains(line,'E I G E N')
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
        end

        if contains(line,'POINT ID.') && contains(line,'T1') && contains(line,'R3')
            inDispTable = true;
            inStrainTable = false;
            if ~isnan(currentMode)
                availableModesDisp(end+1) = currentMode; %#ok<AGROW>
            end
            continue;
        end

        if contains(line,'S T R A I N S') && contains(line,'C R O D')
            k = eigKey(currentEig);
            if isKey(eigToMode,k)
                currentMode = eigToMode(k);
            end
            inStrainTable = true;
            inDispTable = false;
            if ~isnan(currentMode)
                availableModesStrain(end+1) = currentMode; %#ok<AGROW>
            end
            continue;
        end

        if inDispTable && ~isnan(currentMode)
            [modeIsUsed, modeCol] = ismember(currentMode, modesToUse);
            if ~modeIsUsed
                continue;
            end

            if ~contains(line,' G')
                continue;
            end

            nums = sscanf(line, '%f G %f %f %f %f %f %f').';
            if numel(nums) < 7
                continue;
            end

            foundDispRows = true;
            gridID = nums(1);
            vals = nums(2:7);

            if collectDispRows
                dispRows(end+1,:) = [currentMode, gridID, vals]; %#ok<AGROW>
            else
                iDisp = find(dispGridIDs == gridID, 1);
                if ~isempty(iDisp)
                    PHI(iDisp, modeCol) = vals(dispCol);
                end

                iRot = find(rotGridIDs == gridID, 1);
                if ~isempty(iRot)
                    PHI_ROT(iRot, modeCol) = vals(rotCol);
                end
            end

            continue;
        end

        if inStrainTable && ~isnan(currentMode)
            [modeIsUsed, modeCol] = ismember(currentMode, modesToUse);
            if ~modeIsUsed || ~lineStartsWithDigit(line)
                continue;
            end

            nums = sscanf(line,'%f').';
            if numel(nums) < 3
                continue;
            end

            nTriples = floor(numel(nums)/3);
            for j = 1:nTriples
                base = 3*(j-1);
                elemID = nums(base+1);
                axial = nums(base+2);

                foundStrainRows = true;

                if collectStrainRows
                    torsion = nums(base+3);
                    strainRows(end+1,:) = [currentMode, elemID, axial, torsion]; %#ok<AGROW>
                else
                    iStrain = find(strainElemIDs == elemID, 1);
                    if ~isempty(iStrain)
                        PSI(iStrain, modeCol) = axial;
                    end
                end
            end
        end
    end

    if isempty(availableModesDisp)
        error('No REAL EIGEN VECTOR displacement tables were parsed from the file.');
    end
    if ~foundDispRows
        warning('No displacement rows were parsed for the requested modes.');
    end
    if isempty(availableModesStrain) || ~foundStrainRows
        warning('No CROD strain rows were parsed from the file. PSI will be empty or NaN unless strain output exists.');
    end

    if collectDispRows
        dispRows = uniqueRowsByModeID(dispRows, 1, 2);

        if isempty(dispGridIDs) && ~isempty(dispRows)
            dispGridIDs = unique(dispRows(:,2)).';
        end
        if isempty(rotGridIDs)
            rotGridIDs = dispGridIDs;
        end

        PHI = nan(numel(dispGridIDs), nm);
        PHI_ROT = nan(numel(rotGridIDs), nm);

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
        end
    end

    if collectStrainRows
        strainRows = uniqueRowsByModeID(strainRows, 1, 2);

        if isempty(strainElemIDs) && ~isempty(strainRows)
            strainElemIDs = unique(strainRows(:,2)).';
        end

        PSI = nan(numel(strainElemIDs), nm);

        for im = 1:nm
            m = modesToUse(im);

            for i = 1:numel(strainElemIDs)
                idx = strainRows(:,1)==m & strainRows(:,2)==strainElemIDs(i);
                if any(idx)
                    PSI(i,im) = strainRows(find(idx,1,'last'), 3);
                end
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
    info.availableModes_displacement = unique(availableModesDisp);
    info.availableModes_strain = unique(availableModesStrain);
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
function tf = lineStartsWithDigit(line)
    idx = find(~isspace(line), 1);
    tf = ~isempty(idx) && line(idx) >= '0' && line(idx) <= '9';
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
