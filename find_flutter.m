function [Vf,ff] = find_flutter(airspeed,eigSource,numPhysicalModes)

if nargin >= 3 && ~isempty(numPhysicalModes)
    data = omega_v_g_data(airspeed, eigSource, numPhysicalModes);

    if ~data.flutter.hasFlutter
        disp('No damping zero crossing found for the physical modes.')
        Vf = NaN;
        ff = NaN;
        return
    end

    Vf = data.flutter.Vf;
    ff = data.flutter.ff;

    fprintf('Flutter occurs near V = %.3f m/s\n', Vf)
    fprintf('Flutter frequency approx %.3f Hz\n', ff)
    return
end

R = real(eigSource);
S = sign(R);

crossMat = (S(:,1:end-1) < 0) & (S(:,2:end) > 0);
[state_cross, j_cross] = find(crossMat, 1, 'first');

% ---- No flutter ----
if isempty(state_cross)
    disp('No zero crossing found for any airspeed.')
    Vf = NaN;
    ff = NaN;
    return
end

% ---- Linear interpolation for better Vf ----
R1 = R(state_cross, j_cross);
R2 = R(state_cross, j_cross+1);

V1 = airspeed(j_cross);
V2 = airspeed(j_cross+1);

Vf = V1 - R1*(V2 - V1)/(R2 - R1);   % Linear interpolation

% ---- Frequency at crossing ----
lambda1 = eigSource(state_cross,j_cross);
lambda2 = eigSource(state_cross,j_cross+1);

omega1 = imag(lambda1);
omega2 = imag(lambda2);

omega_f = omega1 + (Vf - V1)*(omega2 - omega1)/(V2 - V1);
ff = abs(omega_f)/(2*pi);

fprintf('Flutter occurs near V = %.3f m/s\n', Vf)
fprintf('Flutter frequency approx %.3f Hz\n', ff)

end
