function [RFA_mat] = extractRFAmatrices(filename,method)
% This function reads an .mtx file output from ZAERO MINSTAT card in a
% given folder, extracts the data inside, divides it to the different
% matrices saved in the file and outputs them as a structure that contains
% all the matrices
% The method for the RFA (Roger's/Minimum State) should be given as input
% The output is a structure
% The format of the resulting .mtx file, written in OUTPUT4 format, is
% described in the ZAERO user guide in the "ASSIGN MATRIX =" executive
% control input definition

fid= fopen(filename,'rt');
if (fid == -1)
    warndlg({[ filename, ' doesn`t exist !'],'Please choose another file'});
    [file,path] = uigetfile('*.mtx');
    if file==0
        return
    end
    filename = [path,file];
end
fclose(fid);

if ~exist('method','var') || isempty(method)
    
    errordlg({'No RFA method was defined !',...
        'Please define either Rg for Roger''s method or MS for Karpel''s minimus state method'});
    return
    
elseif ~(strcmp(method,'Rg') || strcmp(method,'MS'))
    
    errordlg({'Wrong RFA method definition !',...
        'Please define either Rg for Roger''s method or MS for Karpel''s minimus state method'});
    return
    
end

RFA_mat.A0 = Read_out4_to_mat_form(filename,'AH0RB');
% The size of matrix A0 is [Nhh x Nhh]

RFA_mat.A1 = Read_out4_to_mat_form(filename,'AH1RB');
% The size of matrix A1 is [Nhh x Nhh]

RFA_mat.A2 = Read_out4_to_mat_form(filename,'AH2RB');
% The size of matrix A2 is [Nhh x Nhh]

RFA_mat.D = Read_out4_to_mat_form(filename,'DHRB');
% The size of matrix D depends on the RFA method:
% For Roger's method, the size of D is [Nhh x Nhh*Nlag]
% For Minimum State method, the size of D is [Nhh x Nlag]

RFA_mat.E = Read_out4_to_mat_form(filename,'ERB');
% The size of matrix E depends on the RFA method:
% For Roger's method, the size of E is [Nhh*Nlag x Nhh]
% For Minimum State method, the size of E is [Nlag x Nhh]

% RFA_mat.R = Read_out4_to_mat_form(filename,'RRAAP');
R = Read_out4_to_mat_form(filename,'RRAAP');
RFA_mat.R = diag(R);
% The size of matrix R depends on the RFA method:
% For Roger's method, the size of R is [Nhh*Nlag x Nhh*Nlag]
% For Minimum State method, the size of R is [Nlag x Nlag]

switch method
    
    case 'Rg'
        RFA_mat.Nlag = size(RFA_mat.D,2)/size(RFA_mat.D,1);
        
    case 'MS'
        RFA_mat.Nlag = size(RFA_mat.D,2);
end

% Changing sign to match ZAERO computed QHH, which are positive on the LHS
% of the aeroelastic equation
RFA_mat.A0 = -RFA_mat.A0;
RFA_mat.A1 = -RFA_mat.A1;
RFA_mat.A2 = -RFA_mat.A2;
RFA_mat.D = -RFA_mat.D;

end