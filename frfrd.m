function [V,rho,nf,f,nu1,ny,FRF1,nu2,FRF2] = frfrd(filename)
%FRFRD Read a frequency-response file written by ZAERO ASEOUT.
%
% The first FRF block is associated with ordinary/control inputs.  The
% optional second block, FRF2, is associated with gust inputs.

arguments
    filename (1,:) char
end

fid = fopen(filename,'rt');
if fid < 0
    error('Could not open ZAERO FRF file "%s".',filename);
end
cleanup = onCleanup(@() fclose(fid));

header = fgetl(fid);
if ~ischar(header)
    error('ZAERO FRF file "%s" is empty.',filename);
end
offset = 3;
if startsWith(header,' ASE=')
    offset = 0;
end
V = fixed_number(header,36+offset,46+offset,'velocity');
rho = fixed_number(header,56+offset,66+offset,'density');

fgetl(fid);
dimensionLine = fgetl(fid);
if ~ischar(dimensionLine) || length(dimensionLine) < 40
    error('Could not read the frequency-count header in "%s".',filename);
end
nf = str2double(strtrim(dimensionLine(35:40)));
if ~isfinite(nf) || nf < 1
    error('Invalid frequency count in "%s".',filename);
end

f = zeros(nf,1);
FRF1 = cell(nf,1);
FRF2 = cell(nf,1);
nu1 = 0;
nu2 = 0;
ny = 0;

for ik = 1:nf
    line = find_record(fid,' IK=');
    index = round(str2double(strtrim(line(5:9))));
    f(index) = str2double(strtrim(line(14:25)));

    require_line(fid,filename);
    require_line(fid,filename);
    [nr,nc,FRF1{index}] = read_output4_matrix(fid,filename);
    ny = nr;
    nu1 = nc;

    require_line(fid,filename);
    secondHeader = require_line(fid,filename);
    if length(secondHeader) >= 50
        [nr2,nc2,FRF2{index}] = read_output4_matrix(fid,filename);
        if nr2 ~= ny
            error('FRF1 and FRF2 output dimensions differ in "%s".',filename);
        end
        nu2 = nc2;
    else
        FRF2{index} = complex(zeros(ny,0));
    end
end

clear cleanup
end

function value = fixed_number(line,first,last,description)
if length(line) < last
    error('ZAERO FRF header is too short to read %s.',description);
end
value = str2double(strtrim(line(first:last)));
if ~isfinite(value)
    error('Could not read %s from the ZAERO FRF header.',description);
end
end

function line = find_record(fid,prefix)
while true
    line = fgetl(fid);
    if ~ischar(line)
        error('Unexpected end of ZAERO FRF file while looking for "%s".',prefix);
    end
    if startsWith(line,prefix)
        return
    end
end
end

function line = require_line(fid,filename)
line = fgetl(fid);
if ~ischar(line)
    error('Unexpected end of ZAERO FRF file "%s".',filename);
end
end

function [nr,nc,matrix] = read_output4_matrix(fid,filename)
header = require_line(fid,filename);
if length(header) < 44
    error('Invalid OUTPUT4 matrix header in "%s".',filename);
end
dims = sscanf(header(1:32),'%d');
if numel(dims) < 2
    error('Could not read OUTPUT4 matrix dimensions in "%s".',filename);
end
nc = dims(1);
nr = dims(2);
valuesPerLine = str2double(header(44));
fieldWidth = 16;
if valuesPerLine == 3
    fieldWidth = 23;
end
matrix = complex(zeros(nr,nc));

for columnRecord = 1:nc+1
    record = require_line(fid,filename);
    location = sscanf(record(1:min(24,length(record))),'%d');
    if numel(location) < 3
        error('Invalid OUTPUT4 column record in "%s".',filename);
    end
    column = location(1);
    firstRow = location(2);
    numberOfRealValues = location(3);
    data = read_fixed_width_values(fid,numberOfRealValues, ...
        valuesPerLine,fieldWidth,filename);
    if column <= nc
        if mod(numberOfRealValues,2) ~= 0
            error('Complex OUTPUT4 record has an odd value count in "%s".',filename);
        end
        numberOfComplexValues = numberOfRealValues/2;
        rows = firstRow:firstRow+numberOfComplexValues-1;
        matrix(rows,column) = complex(data(1:2:end),data(2:2:end));
    end
end
end

function data = read_fixed_width_values(fid,count,perLine,width,filename)
data = zeros(count,1);
position = 0;
while position < count
    line = require_line(fid,filename);
    numberOnLine = min(perLine,count-position);
    for k = 1:numberOnLine
        first = (k-1)*width+1;
        last = min(k*width,length(line));
        position = position+1;
        data(position) = str2double(strrep(strtrim(line(first:last)),'D','E'));
        if ~isfinite(data(position))
            error('Invalid numeric value in ZAERO FRF file "%s".',filename);
        end
    end
end
end
