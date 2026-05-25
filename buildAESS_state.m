function [Aae, Bae, Baw, Cae, Caw, Dae] = buildAESS_state( ...
    Nel, eigvals, zeta, RFA_source, method, L, V, rho, isinp, isgust, ...
    PSI, PHI, PHI_ROT)
%BUILDAESS_STATE
% Build aeroelastic state-space matrices and output matrices.
%
% Returns:
%   Aae  - aeroelastic state matrix
%   Bae  - input matrix for control-surface inputs u_ae = [delta; delta_dot; delta_ddot]
%   Baw  - input matrix for gust input wbar = [w_g/V; wdot_g/V]
%   Cae  - output matrix for sensor outputs
%   Caw  - direct gust-output matrix
%   Dae  - direct control-output matrix

% -------------------------------------------------------------------------
% Input checks
% -------------------------------------------------------------------------
if ~(strcmp(method,'Rg') || strcmp(method,'MS'))
    error('Wrong RFA method. Use either ''Rg'' or ''MS''.');
end

if nargin < 12 || isempty(PSI)
    PSI = zeros(0, Nel);
end
if nargin < 13 || isempty(PHI)
    PHI = zeros(0, Nel);
end
if nargin < 14 || isempty(PHI_ROT)
    PHI_ROT = zeros(0, Nel);
end

check_modal_matrix(PSI,     Nel, 'PSI');
check_modal_matrix(PHI,     Nel, 'PHI');
check_modal_matrix(PHI_ROT, Nel, 'PHI_ROT');

% -------------------------------------------------------------------------
% RFA matrices
% -------------------------------------------------------------------------
if isstruct(RFA_source)
    RFA_mat = RFA_source;
else
    if isstring(RFA_source)
        RFA_source = char(RFA_source);
    end
    RFA_mat = extractRFAmatrices(RFA_source, method);
end

validate_rfa_matrix_struct(RFA_mat);

Nlag = RFA_mat.Nlag;

mode_idx = 1:Nel;

switch method
    case 'Rg'
        Xae_idx = zeros(1, Nel*Nlag);
        for ii = 1:Nel
            Xae_idx(Nlag*(ii-1)+1:ii*Nlag) = ...
                Nlag*(mode_idx(ii)-1)+1 : Nlag*mode_idx(ii);
        end
        Nxae = length(Xae_idx);

    case 'MS'
        Xae_idx = 1:Nlag;
        Nxae = Nlag;
end

% -------------------------------------------------------------------------
% Structural modal matrices
% -------------------------------------------------------------------------
qdyn  = 0.5*rho*V^2;
Omega = sqrt(eigvals(:));

K = diag(eigvals(:));
Cstruct = 2*zeta*diag(Omega);
M = eye(Nel);

% Initialize outputs
Aae = [];
Bae = [];
Baw = [];
Cae = [];
Caw = [];
Dae = [];

% -------------------------------------------------------------------------
% Case 1: no gust, no control input
% -------------------------------------------------------------------------
if isgust==0 && isinp==0

    A0h = RFA_mat.A0(mode_idx, mode_idx);
    A1h = RFA_mat.A1(mode_idx, mode_idx);
    A2h = RFA_mat.A2(mode_idx, mode_idx);

    Dlag = RFA_mat.D(mode_idx, Xae_idx);
    Rlag = RFA_mat.R(Xae_idx, Xae_idx);
    Eh   = RFA_mat.E(Xae_idx, mode_idx);

    Kbar = K       - qdyn*A0h;
    Cbar = Cstruct - qdyn*(L/V)*A1h;
    Mbar = M       - qdyn*(L/V)^2*A2h;

    invMbar = inv(Mbar);

    Aae = [ zeros(Nel)        eye(Nel)          zeros(Nel,Nxae); ...
           -invMbar*Kbar     -invMbar*Cbar      qdyn*invMbar*Dlag; ...
            zeros(Nxae,Nel)   Eh                (V/L)*Rlag ];

    Bae = zeros(size(Aae,1),0);
    Baw = zeros(size(Aae,1),0);

% -------------------------------------------------------------------------
% Case 2: no gust, with control input
% -------------------------------------------------------------------------
elseif isgust==0 && isinp==1

    Nctrl = size(RFA_mat.A0,2) - Nel;
    if Nctrl < 1
        error('isinp=1, but no control-surface columns were found in RFA_mat.A0.');
    end

    ctrl_idx = Nel+1 : Nel+Nctrl;

    A0h = RFA_mat.A0(mode_idx, mode_idx);
    A1h = RFA_mat.A1(mode_idx, mode_idx);
    A2h = RFA_mat.A2(mode_idx, mode_idx);

    A0c = RFA_mat.A0(mode_idx, ctrl_idx);
    A1c = RFA_mat.A1(mode_idx, ctrl_idx);
    A2c = RFA_mat.A2(mode_idx, ctrl_idx); %#ok<NASGU>
    Ec  = RFA_mat.E(Xae_idx, ctrl_idx);

    Dlag = RFA_mat.D(mode_idx, Xae_idx);
    Rlag = RFA_mat.R(Xae_idx, Xae_idx);
    Eh   = RFA_mat.E(Xae_idx, mode_idx);

    Kbar = K       - qdyn*A0h;
    Cbar = Cstruct - qdyn*(L/V)*A1h;
    Mbar = M       - qdyn*(L/V)^2*A2h;

    invMbar = inv(Mbar);

    Aae = [ zeros(Nel)        eye(Nel)          zeros(Nel,Nxae); ...
           -invMbar*Kbar     -invMbar*Cbar      qdyn*invMbar*Dlag; ...
            zeros(Nxae,Nel)   Eh                (V/L)*Rlag ];

    % u_ae = [delta; delta_dot; delta_ddot]
    %
    % According to your current formulation:
    % xi_ddot receives:
    %   q*invM*A0c*delta + q*(L/V)*invM*A1c*delta_dot + 0*delta_ddot
    % Xa_dot receives:
    %   0*delta + Ec*delta_dot + 0*delta_ddot
    Bae = [ zeros(Nel, Nctrl),              zeros(Nel, Nctrl),                 zeros(Nel, Nctrl); ...
            qdyn*invMbar*A0c,               qdyn*(L/V)*invMbar*A1c,            zeros(Nel, Nctrl); ...
            zeros(Nxae, Nctrl),             Ec,                                zeros(Nxae, Nctrl) ];

    Baw = zeros(size(Aae,1),0);

% -------------------------------------------------------------------------
% Case 3: gust, no control input
% -------------------------------------------------------------------------
elseif isgust==1 && isinp==0

    gust_idx = size(RFA_mat.A0,2);

    A0h = RFA_mat.A0(mode_idx, mode_idx);
    A1h = RFA_mat.A1(mode_idx, mode_idx);
    A2h = RFA_mat.A2(mode_idx, mode_idx);

    A0g = RFA_mat.A0(mode_idx, gust_idx);
    A1g = RFA_mat.A1(mode_idx, gust_idx);

    Dlag = RFA_mat.D(mode_idx, Xae_idx);
    Rlag = RFA_mat.R(Xae_idx, Xae_idx);
    Eh   = RFA_mat.E(Xae_idx, mode_idx);
    Eg   = RFA_mat.E(Xae_idx, gust_idx);

    Kbar = K       - qdyn*A0h;
    Cbar = Cstruct - qdyn*(L/V)*A1h;
    Mbar = M       - qdyn*(L/V)^2*A2h;

    invMbar = inv(Mbar);

    Aae = [ zeros(Nel)        eye(Nel)          zeros(Nel,Nxae); ...
           -invMbar*Kbar     -invMbar*Cbar      qdyn*invMbar*Dlag; ...
            zeros(Nxae,Nel)   Eh                (V/L)*Rlag ];

    Bae = zeros(size(Aae,1),0);

    % wbar = [w_g/V; wdot_g/V]
    Baw = [ zeros(Nel,1),                 zeros(Nel,1); ...
            qdyn/V*invMbar*A0g,           qdyn*L/V^2*invMbar*A1g; ...
            zeros(Nxae,1),                (1/V)*Eg ];

% -------------------------------------------------------------------------
% Case 4: gust + control input
% -------------------------------------------------------------------------
elseif isgust==1 && isinp==1

    Nctrl = size(RFA_mat.A0,2) - Nel - 1;
    if Nctrl < 1
        error('isinp=1 and isgust=1, but no control-surface columns were found.');
    end

    ctrl_idx = Nel+1 : Nel+Nctrl;
    gust_idx = Nel+Nctrl+1;

    A0h = RFA_mat.A0(mode_idx, mode_idx);
    A1h = RFA_mat.A1(mode_idx, mode_idx);
    A2h = RFA_mat.A2(mode_idx, mode_idx);

    A0c = RFA_mat.A0(mode_idx, ctrl_idx);
    A1c = RFA_mat.A1(mode_idx, ctrl_idx);
    A2c = RFA_mat.A2(mode_idx, ctrl_idx); %#ok<NASGU>

    A0g = RFA_mat.A0(mode_idx, gust_idx);
    A1g = RFA_mat.A1(mode_idx, gust_idx);

    Dlag = RFA_mat.D(mode_idx, Xae_idx);
    Rlag = RFA_mat.R(Xae_idx, Xae_idx);

    Eh = RFA_mat.E(Xae_idx, mode_idx);
    Ec = RFA_mat.E(Xae_idx, ctrl_idx);
    Eg = RFA_mat.E(Xae_idx, gust_idx);

    Kbar = K       - qdyn*A0h;
    Cbar = Cstruct - qdyn*(L/V)*A1h;
    Mbar = M       - qdyn*(L/V)^2*A2h;

    invMbar = inv(Mbar);

    Aae = [ zeros(Nel)        eye(Nel)          zeros(Nel,Nxae); ...
           -invMbar*Kbar     -invMbar*Cbar      qdyn*invMbar*Dlag; ...
            zeros(Nxae,Nel)   Eh                (V/L)*Rlag ];

    % u_ae = [delta; delta_dot; delta_ddot]
    Bae = [ zeros(Nel, Nctrl),              zeros(Nel, Nctrl),                 zeros(Nel, Nctrl); ...
            qdyn*invMbar*A0c,               qdyn*(L/V)*invMbar*A1c,            zeros(Nel, Nctrl); ...
            zeros(Nxae, Nctrl),             Ec,                                zeros(Nxae, Nctrl) ];

    % wbar = [w_g/V; wdot_g/V]
    Baw = [ zeros(Nel,1),                 zeros(Nel,1); ...
            qdyn/V*invMbar*A0g,           qdyn*L/V^2*invMbar*A1g; ...
            zeros(Nxae,1),                (1/V)*Eg ];

else
    error('Invalid combination of isinp=%d and isgust=%d.', isinp, isgust);
end

% -------------------------------------------------------------------------
% Build Cae, Caw, Dae
% -------------------------------------------------------------------------
[Cae, Caw, Dae] = build_output_matrices( ...
    PSI, PHI, PHI_ROT, ...
    Nel, Nxae, invMbar, Kbar, Cbar, Dlag, qdyn, ...
    Bae, Baw);

end

% =========================================================================
function [Cae, Caw, Dae] = build_output_matrices( ...
    PSI, PHI, PHI_ROT, ...
    Nel, Nxae, invMbar, Kbar, Cbar, Dlag, qdyn, ...
    Bae, Baw)

nStrain = size(PSI,1);
nAcc    = size(PHI,1);
nGyro   = size(PHI_ROT,1);

% strain = PSI * xi
C_strain = [ ...
    PSI, ...
    zeros(nStrain, Nel), ...
    zeros(nStrain, Nxae) ];

% acceleration = PHI * xi_ddot
%
% xi_ddot =
%   -invMbar*Kbar*xi
%   -invMbar*Cbar*xi_dot
%   + qdyn*invMbar*Dlag*Xa
C_acc = [ ...
   -PHI*invMbar*Kbar, ...
   -PHI*invMbar*Cbar, ...
    qdyn*PHI*invMbar*Dlag ];

% angular_rate = PHI_ROT * xi_dot
C_rate = [ ...
    zeros(nGyro, Nel), ...
    PHI_ROT, ...
    zeros(nGyro, Nxae) ];

Cae = [ ...
    C_strain;
    C_acc;
    C_rate ];

% Direct control-output matrix:
% Only acceleration rows have direct dependence on u through xi_ddot.
if isempty(Bae) || size(Bae,2)==0
    Dae = zeros(nStrain+nAcc+nGyro, 0);
else
    Bae_acc = Bae(Nel+1:2*Nel, :);
    Dae = [ ...
        zeros(nStrain, size(Bae,2));
        PHI*Bae_acc;
        zeros(nGyro, size(Bae,2)) ];
end

% Direct gust-output matrix:
% Only acceleration rows have direct dependence on gust through xi_ddot.
if isempty(Baw) || size(Baw,2)==0
    Caw = zeros(nStrain+nAcc+nGyro, 0);
else
    Baw_acc = Baw(Nel+1:2*Nel, :);
    Caw = [ ...
        zeros(nStrain, size(Baw,2));
        PHI*Baw_acc;
        zeros(nGyro, size(Baw,2)) ];
end

end

% =========================================================================
function check_modal_matrix(X, Nel, name)
if size(X,2) ~= Nel
    error('%s must have Nel=%d columns, but size(%s,2)=%d.', ...
        name, Nel, name, size(X,2));
end
end

% =========================================================================
function validate_rfa_matrix_struct(RFA_mat)
requiredFields = {'A0', 'A1', 'A2', 'D', 'E', 'R', 'Nlag'};
for k = 1:numel(requiredFields)
    if ~isfield(RFA_mat, requiredFields{k})
        error('RFA_mat is missing required field "%s".', requiredFields{k});
    end
end
end
