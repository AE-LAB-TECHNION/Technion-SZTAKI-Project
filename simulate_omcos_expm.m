function simulation = simulate_omcos_expm(cfg, model, time)
%SIMULATE_OMCOS_EXPM Simulate OMCOS using exact LTI propagation.

time = time(:);
Nel = cfg.Nel;
A = model.Aae;
B = model.Baw;
nStates = size(A, 1);
nInputs = size(B, 2);

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

gustStart = (XResponse - X0) / V;
gustEnd = gustStart + Tg;
omega = 2*pi/Tg;

% The same gust is applied to every gust-vector input pair.
bW = sum(B(:, 1:Ngust), 2);
bWdot = sum(B(:, Ngust+1:end), 2);

% Augmented LTI system for the active OMCOS interval:
% wbar     = 0.5*WGV*(1-cos(omega*eta))
% wdotbar = 0.5*WGV*omega*sin(omega*eta)
Factive = zeros(nStates+3);
Factive(1:nStates, 1:nStates) = A;
Factive(1:nStates, nStates+1) = -0.5*WGV*bW;
Factive(1:nStates, nStates+2) = 0.5*WGV*omega*bWdot;
Factive(1:nStates, nStates+3) = 0.5*WGV*bW;
Factive(nStates+1, nStates+2) = -omega;
Factive(nStates+2, nStates+1) = omega;

state = zeros(nStates, 1);
stateOut = zeros(numel(time), nStates);
stateOut(1, :) = state.';

for k = 1:numel(time)-1
    t0 = time(k);
    t1 = time(k+1);
    state = propagate_interval(state, t0, t1, A, Factive, ...
        gustStart, gustEnd, omega, V, X0, XResponse, nStates);
    stateOut(k+1, :) = state.';
end

simulation.time = time;
simulation.states = stateOut;
simulation.modalX = stateOut(:, 1:Nel);
simulation.modalV = stateOut(:, Nel+(1:Nel));
simulation.Tg = Tg;
simulation.WGV = WGV;

end

function state = propagate_interval(state, t0, t1, A, Factive, ...
    gustStart, gustEnd, omega, V, X0, XResponse, nStates)

breaks = [t0, max(t0, min(t1, gustStart)), ...
    max(t0, min(t1, gustEnd)), t1];
breaks = unique(breaks);

for j = 1:numel(breaks)-1
    ta = breaks(j);
    tb = breaks(j+1);
    dt = tb - ta;

    if dt <= 0
        continue
    end

    midpoint = 0.5*(ta + tb);
    eta = midpoint - (XResponse - X0) / V;

    if eta > 0 && eta < (gustEnd - gustStart)
        phase = ta - gustStart;
        augmentedState = [state; cos(omega*phase); sin(omega*phase); 1];
        augmentedState = expm(Factive*dt) * augmentedState;
        state = augmentedState(1:nStates);
    else
        state = expm(A*dt) * state;
    end
end

end