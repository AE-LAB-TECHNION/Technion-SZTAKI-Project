%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function SYSRD to read data file created by ZAERO/ASE module           %
%          =====                                                         %
% using ASEOUT bulk data card.                                           %
%                                                                        %
% Input: filename - name of the data file (Character)                    %
%                                                                        %
% Output: V,rho - velocity and density values                            %
%                                                                        %
% nx,nu,ny - number of states, inputs, and outputs of                    %
%            the state-space model                                       %
%                                                                        %
% ABCD - matrix | A B | of state-space model                             %
%               | C D |                                                  %
%                                                                        %
% A,B,C,D - matrices [A], [B], [C], [D] of                               %
%           the state-space model                                        %
%                                                                        %
% G - control gain matrix                                                %
%                                                                        %
% Bw,CG - gust matrices [Bw] and [CG] of                                 %
%         the state-space model                                          %
%                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [V,rho,nx,nu,ny,ABCD,A,B,C,D,G,Bw,CG] = sysrd(filename)

% Open file
fid = fopen(filename);
disp(['Read file ' filename])

% Read and display the header
head = fgetl(fid);
disp(head)
if strcmp(head(1:5), ' ASE=')
    offset=0;
else
    offset=3;
end

V = str2double(head(36+offset:46+offset));
rho = str2double(head(56+offset:66+offset));
% Read dimensions
tmp = fgetl(fid);
tmp = sscanf(tmp,'%i');

if length(tmp) == 1
    nx = tmp(1);
    nu = 0;
    ny = 0;
elseif length(tmp) == 3
    nx = tmp(1);
    nu = tmp(2);
    ny = tmp(3); 
else
    fclose(fid);
    error(['File ' filename ' corrupt.']);
end

% Read state-space matrices
ABCD = fscanf(fid, '%g', [nx+ny, nx+nu]);
A = ABCD(1:nx, 1:nx);
B = ABCD(1:nx, nx+1:nx+nu);
C = ABCD(nx+1:nx+ny, 1:nx);
D = ABCD(nx+1:nx+ny, nx+1:nx+nu);

% Read gain matrix
if strcmp(head(74+offset:80+offset), 'VEHICLE')
    fgetl(fid);
    fgetl(fid);
    G = fscanf(fid, '%g', [nu, ny]);
else
    G = [];
end

% Read gust state-space matrices
fgetl(fid);
tline = fgetl(fid);
if ~ischar(tline)
    Bw = [];
    CG = [];
else
    nG2 = fscanf(fid, '%g', 1);
    Bw = fscanf(fid, '%g', [nx, nG2]);
    
    fgetl(fid);
    tline = fgetl(fid);
    if ~ischar(tline)
        CG = [];
    else
        CG = fscanf(fid, '%g', [ny, nG2]);
    end
end

% Close file
fclose(fid);

if nargout == 1
    V = struct(...
        'V', V, ...
        'rho', rho, ...
        'nx', nx, ...
        'nu', nu, ...
        'ny', ny, ...
        'ABCD', ABCD, ...
        'A', A, ...
        'B', B, ...
        'C', C, ...
        'D', D, ...
        'G', G, ...
        'Bw', Bw, ...
        'CG', CG);
end
end




