function data = load_zaero_dryden_response(filename)
%LOAD_ZAERO_DRYDEN_RESPONSE Read ZAERO continuous-gust modal frequency responses.

lines = readlines(filename);
headerIndex = find(contains(lines,"FREQ,HZ"));
if numel(headerIndex) < 2
    error('Expected two frequency-response blocks in %s.',filename);
end

first = read_numeric_block(lines,headerIndex(1)+2,11);
second = read_numeric_block(lines,headerIndex(2)+2,3);
if size(first,1) ~= size(second,1) || any(abs(first(:,1)-second(:,1)) > 1e-10)
    error('The ZAERO frequency grids in %s are inconsistent.',filename);
end

frequency = first(:,1);
modalV = complex(first(:,[2 4 6]),first(:,[3 5 7]));
modalX = [complex(first(:,8),first(:,9)), ...
          complex(first(:,10),first(:,11)), ...
          complex(second(:,2),second(:,3))];

data = struct();
data.frequency = frequency;
data.modalX = modalX;
data.modalV = modalV;
% ZAERO prints the complex spectral-response amplitude for CGUST. Squaring
% its magnitude gives the corresponding auto-PSD on the ZAERO grid.
data.modalXPsd = abs(modalX).^2;
data.modalVPsd = abs(modalV).^2;
end

function values = read_numeric_block(lines,startIndex,numberOfColumns)
rows = cell(0,1);
for i = startIndex:numel(lines)
    line = strtrim(lines(i));
    if strlength(line) == 0
        if ~isempty(rows), break; end
        continue
    end
    numbers = sscanf(line,'%f').';
    if numel(numbers) ~= numberOfColumns
        if ~isempty(rows), break; end
        continue
    end
    rows{end+1,1} = numbers; %#ok<AGROW>
end
if isempty(rows)
    error('Could not read a numeric ZAERO frequency-response block.');
end
values = vertcat(rows{:});
end
