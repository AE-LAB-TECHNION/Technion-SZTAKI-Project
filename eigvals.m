function eigvals = read_zaero_eigs(filename, nmodes)
%READ_ZAERO_EIGS  Read first nmodes eigenvalues from a ZAERO output file.
%
%   eigvals = READ_ZAERO_EIGS(filename, nmodes)
%
%   The function looks for the block that starts with the line
%
%     'RIGID BODY DEGREES OF FREEDOM (DEFINED IN THE FEM BASIC COORDINATE SYSTEM) =       0'
%
%   Then it finds the following header line that contains the word 'EIGENVALUE'
%   and reads the next nmodes rows of the table. Each of these rows is assumed
%   to contain, in order:
%       MODE   EXTRACTION ORDER   EIGENVALUE   FREQ(rad/s)   FREQ(Hz)   GMASS   K
%   The function returns the third numeric entry on each row (the EIGENVALUE
%   column, (rad/s)^2) as eigvals(k).

    if nargin < 2
        error('Usage: eigvals = read_zaero_eigs(filename, nmodes)');
    end

    fid = fopen(filename, 'r');
    if fid == -1
        error('Could not open file: %s', filename);
    end
    cleaner = onCleanup(@() fclose(fid));

    %----------------------------------------------------------------------
    % 1. Find the "rigid body DOF = 0" line
    %----------------------------------------------------------------------
    targetRB = 'RIGID BODY DEGREES OF FREEDOM (DEFINED IN THE FEM BASIC COORDINATE SYSTEM)';
    foundRB  = false;

    while true
        line = fgetl(fid);
        if ~ischar(line)
            break
        end
        if contains(line, targetRB) && contains(line, '=       0')
            foundRB = true;
            break
        end
    end

    if ~foundRB
        error('Rigid-body DOF line with "= 0" not found in file %s.', filename);
    end

    %----------------------------------------------------------------------
    % 2. From here, find the header line containing "EIGENVALUE"
    %----------------------------------------------------------------------
    foundHeader = false;
    while true
        line = fgetl(fid);
        if ~ischar(line)
            break
        end
        if contains(line, 'EIGENVALUE')
            foundHeader = true;
            break
        end
    end

    if ~foundHeader
        error('Header line containing "EIGENVALUE" not found after RB DOF line.');
    end

    %----------------------------------------------------------------------
    % 3. Read the next nmodes table lines and extract the 3rd numeric column
    %----------------------------------------------------------------------
    eigvals = nan(nmodes,1);

    for k = 1:nmodes
        line = fgetl(fid);
        if ~ischar(line)
            error('End of file reached before reading %d eigenvalues.', nmodes);
        end

        nums = sscanf(line, '%f');
        % Expect at least: MODE, EXTRACTION ORDER, EIGENVALUE, ...
        if numel(nums) < 3
            error('Could not parse eigenvalue on line %d after header.', k);
        end

        eigvals(k) = nums(3);  % EIGENVALUE column (RAD/S)^2
    end
end
