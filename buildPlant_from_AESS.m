function [Ap, Bp, Bpw] = buildPlant_from_AESS(Aae, Bae, Baw)

% Build plant-level state-space from aeroelastic core + actuators
%
% x_p = [ x_ae ; x_act ]
%
% dot{x}_p = Ap x_p + Bp u_act + Bpw w_g

% --- actuator block diagonal ---
act.A0 = 2.151e6;
act.A1 = 4.772e4;
act.A2 = 586.8;

A_act = [0 1 0 ; 0 0 1 ; -act.A0  -act.A1  -act.A2 ];
B_act = [0; 0; act.A0];

nae = size(Aae,1);
nActStates = size(A_act,1);
% --- number of control surfaces ---
assert(mod(size(Bae,2), nActStates) == 0, ...
    'Bae columns inconsistent with actuator order');
nCS = size(Bae,2) / nActStates;

A_act_big = kron(eye(nCS), A_act);
B_act_big = kron(eye(nCS), B_act);

% --- PLANT matrices ---
Ap = [ ...
    Aae                    Bae ;
    zeros(nActStates*nCS, nae)   A_act_big ];

Bp = [ ...
    zeros(nae, size(B_act_big,2)) ;
    B_act_big ];

if isempty(Baw)
    Bpw = [];
else
    Bpw = [ ...
        Baw ;
        zeros(nActStates*nCS, size(Baw,2)) ];
end

end
