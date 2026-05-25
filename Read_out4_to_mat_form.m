function A = Read_out4_to_mat_form(filename,Mat_Name)
fid= fopen(filename,'rt');
if (fid == -1)
    error([ filename, ' doesn`t exist']);
end

%% Moving through the file to the requested matrix
tline = fgetl(fid);
while ischar(tline)
    if ~contains(tline,Mat_Name)
        tline = fgetl(fid);
        continue;
    end
    
    fseek(fid,-length(tline),'cof');
    
    break 
    
end

%% Record 1: matrix characteristics
tmp=fscanf(fid,' %d ',4);
COL=tmp(1); ROW=tmp(2); NF=tmp(3); Ntype=tmp(4);
tmp=fgetl(fid);
MatrixName=tmp(1:8);
if ~contains(MatrixName,Mat_Name)
    warning('Matrix names are not right, the function is not working properly')
end
MatrixType=tmp(9:17);
A=zeros(ROW,COL);
%% Record 2: matrix data
switch Ntype
    case {1,2}
        for i=1:COL
            tmp=fscanf(fid,' %d ',3);
            iCOL=tmp(1); iROW=tmp(2); NW=tmp(3);
            if iCOL > COL
                break
            end
            A(iROW:iROW+NW-1,iCOL)=fscanf(fid,' %f ',NW);
        end
    case {3,4}
        for i=1:COL
            tmp=fscanf(fid,' %d ',3);
            iCOL=tmp(1); iROW=tmp(2); NW=tmp(3);
            if iCOL > COL
                break
            end
            tmp=fscanf(fid,' %f ',NW);
            A(iROW:iROW+NW/2-1,iCOL)=tmp(1:2:NW-1)+1i*tmp(2:2:NW);
        end
end
fclose('all');

if sum(sum(isnan(A)+isinf(A)))
    warning('NaN or inf numbers in filename %s',filename)
end
end