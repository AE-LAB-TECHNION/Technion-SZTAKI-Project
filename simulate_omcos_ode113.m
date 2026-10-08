function simulation = simulate_omcos_ode113(cfg, model, time)
%SIMULATE_OMCOS_ODE113 Simulate the OMCOS gust with ode113.

time = time(:);
Nel = cfg.Nel;
nInputs = size(model.Baw, 2);

if nInputs < 2 || mod(nInputs, 2) ~= 0
    error('Baw must contain [w_g/V, dw_g/dt/V] input pairs.');
end

Ngust = nInputs / 2;
V = model.V;
LG = cfg.gustLG;
WGV = cfg.gustWGV;
X0 = cfg.gustX0;
XResponse = cfg.gustXResponse;

switch lower(cfg.gustLengthMode)
    case 'time'
        Tg = LG;
    case 'length'
        Tg = LG / V;
    otherwise
        error('gustLengthMode must be ''time'' or ''length''.');
end

if Tg <= 0
    error('OMCOS gust duration must be positive.');
end

options = odeset( ...
    'RelTol', 1e-12, ...
    'AbsTol', 1e-14, ...
    'MaxStep', min(0.001, Tg/100));

x0 = zeros(size(model.Aae, 1), 1);

[timeOut, stateOut] = ode113( ...
    @(t, x) omcos_rhs(t, x, model.Aae, model.Baw, ...
                      V, WGV, Tg, X0, XResponse, Ngust), ...
    time, x0, options);

simulation.time = timeOut;
simulation.states = stateOut;
simulation.modalX = stateOut(:, 1:Nel);
simulation.modalV = stateOut(:, Nel+(1:Nel));
simulation.Tg = Tg;
simulation.WGV = WGV;

end

function dx = omcos_rhs(t, x, A, B, V, WGV, Tg, X0, XResponse, Ngust)

nInputs = size(B, 2);
eta = t - (XResponse - X0) / V;

wbar = 0;
wdotbar = 0;

if eta >= 0 && eta <= Tg
    wbar = 0.5 * WGV * (1 - cos(2*pi*eta/Tg));
    wdotbar = (pi*WGV/Tg) * sin(2*pi*eta/Tg);
end

gustInput = zeros(nInputs, 1);
gustInput(1:Ngust) = wbar;
gustInput(Ngust+1:end) = wdotbar;

dx = A*x + B*gustInput;

end