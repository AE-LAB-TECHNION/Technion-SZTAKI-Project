function zaero = read_zaero_dryden_psd(filename)
%READ_ZAERO_DRYDEN_PSD Read the Dryden gust PSD printed by ZAERO.
%
% Output:
%   zaero.frequency        - frequency vector [Hz]
%   zaero.psd              - Dryden PSD printed by ZAERO
%   zaero.Lg               - turbulence scale length
%   zaero.sigmaW           - input GURMS
%   zaero.computedGURMS    - GURMS computed by ZAERO over [0,FMAX]
%   zaero.V                - flight velocity
%   zaero.df               - ZAERO frequency spacing
%   zaero.nFrequency       - number of PSD frequency points
%   zaero.filename         - source file

    if ~isfile(filename)
        error('ZAERO output file not found: %s', filename);
    end

    txt = fileread(filename);

    %% ------------------------------------------------------------
    % 1. Read Dryden parameters
    %
    % ZAERO line example:
    %
    % LOG10(PSD) OF DRYDEN GUST SPECTRUM,
    % LG=7.500E-01,
    % INPUT GURMS=7.000E-01,
    % COMPUTED GURMS=6.748E-01,
    % VEL=3.500E+01
    % -------------------------------------------------------------

    numberPattern = ...
        '([+\-]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][+\-]?\d+)?)';

    headerPattern = [ ...
        'LOG10\(PSD\)\s+OF\s+DRYDEN\s+GUST\s+SPECTRUM,\s*' ...
        'LG\s*=\s*' numberPattern ',\s*' ...
        'INPUT\s+GURMS\s*=\s*' numberPattern ',\s*' ...
        'COMPUTED\s+GURMS\s*=\s*' numberPattern ',\s*' ...
        'VEL\s*=\s*' numberPattern];

    tokens = regexpi(txt, headerPattern, 'tokens', 'once');

    if isempty(tokens)
        error(['Could not find the ZAERO Dryden spectrum header in ', ...
               '"%s". Make sure CGUST TYPE = DRYP.'], filename);
    end

    Lg            = zaeroNumber(tokens{1});
    sigmaW        = zaeroNumber(tokens{2});
    computedGURMS = zaeroNumber(tokens{3});
    V             = zaeroNumber(tokens{4});

    %% ------------------------------------------------------------
    % 2. Find numerical PSD table
    %
    % Example:
    %
    % GUST PSD AT 10001 FREQUENCY STEPS WITH DF= 0.10000E-01 HZ
    % -------------------------------------------------------------

    tablePattern = [ ...
        'GUST\s+PSD\s+AT\s+(\d+)\s+' ...
        'FREQUENCY\s+STEPS\s+WITH\s+DF\s*=\s*' ...
        numberPattern '\s*HZ'];

    [~, tableEnd, tableTokens] = ...
        regexpi(txt, tablePattern, 'start', 'end', 'tokens', 'once');

    if isempty(tableTokens)
        error('Could not find the numerical GUST PSD table in "%s".', ...
              filename);
    end

    nFrequency = str2double(tableTokens{1});
    df = zaeroNumber(tableTokens{2});

    %% ------------------------------------------------------------
    % 3. Read exactly nFrequency numerical rows
    % -------------------------------------------------------------

    remainingText = txt(tableEnd+1:end);

    frequency = zeros(nFrequency,1);
    psd       = zeros(nFrequency,1);

    dataPattern = [ ...
        '^\s*' numberPattern ...
        '\s+' numberPattern '\s*$'];

    lines = regexp(remainingText, '\r\n|\n|\r', 'split');

    count = 0;

    for i = 1:numel(lines)

        row = regexp(lines{i}, dataPattern, ...
                     'tokens', 'once', 'lineanchors');

        if isempty(row)
            continue
        end

        count = count + 1;

        frequency(count) = zaeroNumber(row{1});
        psd(count)       = zaeroNumber(row{2});

        if count == nFrequency
            break
        end
    end

    if count ~= nFrequency
        error(['ZAERO reports %d Dryden PSD points, but only %d ', ...
               'numerical rows were read.'], ...
              nFrequency, count);
    end

    %% ------------------------------------------------------------
    % 4. Pack output
    % -------------------------------------------------------------

    zaero = struct();

    zaero.frequency     = frequency;
    zaero.psd           = psd;
    zaero.Lg            = Lg;
    zaero.sigmaW        = sigmaW;
    zaero.computedGURMS = computedGURMS;
    zaero.V             = V;
    zaero.df            = df;
    zaero.nFrequency    = nFrequency;
    zaero.filename      = filename;

end


% ========================================================================
function value = zaeroNumber(str)
%ZAERONUMBER Convert ZAERO E/D scientific notation to MATLAB double.

    str = strrep(str, 'D', 'E');
    str = strrep(str, 'd', 'e');

    value = str2double(str);

    if ~isfinite(value)
        error('Could not convert ZAERO number "%s".', str);
    end
end