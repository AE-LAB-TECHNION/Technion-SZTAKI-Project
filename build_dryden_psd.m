function Phi = build_dryden_psd(frequency, Lg, sigmaW, V)
%BUILD_DRYDEN_PSD Build Dryden vertical-gust PSD matching ZAERO output.
%
% Inputs:
%   frequency - frequency [Hz]
%   Lg        - turbulence scale length
%   sigmaW    - RMS gust velocity
%   V         - flight velocity
%
% Output:
%   Phi       - Dryden PSD using the convention printed by ZAERO
%
% The ZAERO output examined here is reproduced by
%
% Phi(f) = sigmaW^2 * (Lg/V) *
%          [1 + 3*(omega*Lg/V)^2] /
%          [1 + (omega*Lg/V)^2]^2
%
% with
%   omega = 2*pi*f

    arguments
        frequency (:,1) double
        Lg        (1,1) double {mustBePositive}
        sigmaW    (1,1) double {mustBeNonnegative}
        V         (1,1) double {mustBePositive}
    end

    omega = 2*pi*frequency;

    r = omega .* Lg ./ V;

    Phi = sigmaW^2 .* (Lg/V) .* ...
          (1 + 3*r.^2) ./ ...
          (1 + r.^2).^2;

end