function [Ap, Bp, Bpw] = buildPlant_from_AESS(Aae, Bae, Baw, nCS)

% Build plant-level state-space from aeroelastic core + actuators
%
% x_p = [ x_ae ; x_act ]
% x_act = [delta_1 ... delta_N, delta_dot_1 ... delta_dot_N, delta_ddot_1 ... delta_ddot_N]
%
% dot{x}_p = Ap x_p + Bp u_act + Bpw w_g

% --- actuator block diagonal ---
act.A0 = 2.151e6;
act.A1 = 4.772e4;
act.A2 = 586.8;

nae = size(Aae,1);
nActStatesPerCS = 3;

if nargin < 4 || isempty(nCS)
    assert(mod(size(Bae,2), nActStatesPerCS) == 0, ...
        'Bae columns inconsistent with actuator order');
    nCS = size(Bae,2) / nActStatesPerCS;
end

Bae = select_control_surface_columns(Bae, nCS);

A_act_big = [ ...
    zeros(nCS)      eye(nCS)        zeros(nCS); ...
    zeros(nCS)      zeros(nCS)      eye(nCS); ...
   -act.A0*eye(nCS) -act.A1*eye(nCS) -act.A2*eye(nCS) ];

B_act_big = [ ...
    zeros(nCS); ...
    zeros(nCS); ...
    act.A0*eye(nCS) ];

% --- PLANT matrices ---
Ap = [ ...
    Aae                         Bae ;
    zeros(nActStatesPerCS*nCS, nae)   A_act_big ];

Bp = [ ...
    zeros(nae, nCS) ;
    B_act_big ];

if isempty(Baw)
    Bpw = [];
else
    Bpw = [ ...
        Baw ;
        zeros(nActStatesPerCS*nCS, size(Baw,2)) ];
end

end

% =========================================================================
function BaeSelected = select_control_surface_columns(Bae, nCS)

nActStatesPerCS = 3;

assert(mod(size(Bae,2), nActStatesPerCS) == 0, ...
    'Bae columns inconsistent with actuator order');

nBaeCS = size(Bae,2) / nActStatesPerCS;

if nBaeCS < nCS
    error('Expected %d control surfaces, but Bae only contains %d.', nCS, nBaeCS);
end

if nBaeCS > nCS
    warning(['Bae contains %d control-surface input groups, but the ZAERO RFA reports %d. ', ...
             'Using the first %d groups and ignoring the rest.'], nBaeCS, nCS, nCS);
end

BaeSelected = [ ...
    Bae(:, 1:nCS), ...
    Bae(:, nBaeCS+1:nBaeCS+nCS), ...
    Bae(:, 2*nBaeCS+1:2*nBaeCS+nCS) ];

end
