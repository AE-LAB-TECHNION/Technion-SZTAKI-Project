clc; clear; close all;

% -------------------------------------------------------------------------
% Defaults
% -------------------------------------------------------------------------
run('./defaults.m')

% -------------------------------------------------------------------------
% User settings
% -------------------------------------------------------------------------
rho = 1.225;
L   = 0.05;        % length scale [m], usually half-chord
zeta = 0;          % modal damping ratio

isinp  = 0;        % 1 = plot Ap, 0 = plot Aae
isgust = 0;        % 1 = include gust column, 0 = no gust

nmodes = 4;        % number of structural modes to read from ZAERO output

RFA_filename   = 'APPROX.DAT';
ZAERO_filename = 'ASE_ANALYSIS_new.out';
f06file        = 'model-0012.f06';

method = 'MS';     % 'MS' or 'Rg'

% ZAERO uses modes 1,2,4 because mode 3 is omitted.
selected_modes = [1 2 4];
Nel = length(selected_modes);

% -------------------------------------------------------------------------
% Sensor definitions for PSI, PHI, PHI_ROT
% -------------------------------------------------------------------------
dispGridIDs = [44 260 2691];
dispComp    = 'T3';      % T1/T2/T3

rotGridIDs = [44 260 2691];
rotComp    = 'R2';       % R1/R2/R3

strainElemIDs = [13915:13924, 13935:13944, 17293:17296, 17530:17897];

% -------------------------------------------------------------------------
% Read structural eigenvalues
% -------------------------------------------------------------------------
all_eigs = read_zaero_eigs(ZAERO_filename, nmodes);
eigvals = all_eigs(selected_modes);

fprintf('Selected modes: ');
fprintf('%d ', selected_modes);
fprintf('\n');

fprintf('Eigenvalues used:\n');
disp(eigvals(:).');

% -------------------------------------------------------------------------
% Extract PSI, PHI, PHI_ROT once.
% -------------------------------------------------------------------------
[PSI, PHI, PHI_ROT, sensorInfo] = extract_PSI_PHI_PHIROT_from_F06( ...
    f06file, selected_modes, ...
    dispGridIDs, dispComp, ...
    rotGridIDs, rotComp, ...
    strainElemIDs);

fprintf('\nSensor/modal matrix sizes:\n');
fprintf('  PSI     : %d x %d\n', size(PSI,1),     size(PSI,2));
fprintf('  PHI     : %d x %d\n', size(PHI,1),     size(PHI,2));
fprintf('  PHI_ROT : %d x %d\n', size(PHI_ROT,1), size(PHI_ROT,2));

if any(isnan(PSI(:))) || any(isnan(PHI(:))) || any(isnan(PHI_ROT(:)))
    warning(['At least one of PSI/PHI/PHI_ROT contains NaN. ', ...
             'Check that the requested modes, GRID IDs, and element IDs exist in the F06.']);
end

% -------------------------------------------------------------------------
% Airspeed loop
% -------------------------------------------------------------------------
airspeed = 25:0.25:50;
nV = length(airspeed);

eigTraj = [];
models = cell(nV,1);

for i = 1:nV

    V = airspeed(i);

    [Aae, Bae, Baw, Cae, Caw, Dae] = buildAESS_state( ...
        Nel, eigvals, zeta, ...
        RFA_filename, method, ...
        L, V, rho, ...
        isinp, isgust, ...
        PSI, PHI, PHI_ROT);

    % ---------------------------------------------------------------------
    % Build plant and choose which A matrix to plot
    % ---------------------------------------------------------------------
    if isinp == 1

        [Ap, Bp, Bpw] = buildPlant_from_AESS(Aae, Bae, Baw);
        Cp = [];
        Aact = [];
        Bact = [];

        Aplot = Ap;
        plotModelName = 'Ap';

    else

        Ap = [];
        Bp = [];
        Bpw = [];
        Cp = [];
        Aact = [];
        Bact = [];

        Aplot = Aae;
        plotModelName = 'Aae';

    end

    % ---------------------------------------------------------------------
    % Eigenvalues for root-locus plot
    % ---------------------------------------------------------------------
    lam = eig(Aplot);
    eigTraj(:,i) = lam;

    % ---------------------------------------------------------------------
    % Store models
    % ---------------------------------------------------------------------
    models{i}.V    = V;

    models{i}.Aae  = Aae;
    models{i}.Bae  = Bae;
    models{i}.Baw  = Baw;
    models{i}.Cae  = Cae;
    models{i}.Caw  = Caw;
    models{i}.Dae  = Dae;

    models{i}.Ap   = Ap;
    models{i}.Bp   = Bp;
    models{i}.Bpw  = Bpw;
    models{i}.Cp   = Cp;
    models{i}.Aact = Aact;
    models{i}.Bact = Bact;

    if i == 1
        fprintf('\nMatrix sizes at V = %.2f m/s:\n', V);
        fprintf('  Aae : %d x %d\n', size(Aae,1), size(Aae,2));
        fprintf('  Bae : %d x %d\n', size(Bae,1), size(Bae,2));
        fprintf('  Baw : %d x %d\n', size(Baw,1), size(Baw,2));
        fprintf('  Cae : %d x %d\n', size(Cae,1), size(Cae,2));
        fprintf('  Caw : %d x %d\n', size(Caw,1), size(Caw,2));
        fprintf('  Dae : %d x %d\n', size(Dae,1), size(Dae,2));

        if isinp == 1
            fprintf('  Ap  : %d x %d\n', size(Ap,1), size(Ap,2));
            fprintf('  Bp  : %d x %d\n', size(Bp,1), size(Bp,2));
            fprintf('  Bpw : %d x %d\n', size(Bpw,1), size(Bpw,2));
        end

        fprintf('\nRoot-locus plot uses: %s\n', plotModelName);
    end
end

% -------------------------------------------------------------------------
% Plot root locus of selected matrix: Ap or Aae
% -------------------------------------------------------------------------
figure(101); clf; hold on; grid on; box on;

for j = 1:nV

    Vj = airspeed(j);

    scatter(real(eigTraj(:,j)), imag(eigTraj(:,j)), ...
        25, Vj*ones(size(eigTraj,1),1), 'filled');

end

colormap(jet)
cb = colorbar;
cb.Label.String = 'Airspeed [m/s]';

xlabel('Real$(\lambda)$','Interpreter','latex');
ylabel('Imag$(\lambda)$','Interpreter','latex');

title(sprintf('Root Locus of %s Eigenvalues', plotModelName), ...
      'Interpreter','latex');

xlim([-100, 100])
ylim([0, 200])

set(gca, 'FontSize', 13);

% -------------------------------------------------------------------------
% Flutter estimate
% -------------------------------------------------------------------------
[Vf, ff] = find_flutter(airspeed, eigTraj);

fprintf('\nFlutter estimate based on %s eigenvalues:\n', plotModelName);
fprintf('Flutter occurs near V = %.3f m/s\n', Vf);
fprintf('Flutter frequency ≈ %.3f Hz\n', ff);

%% ------------------------------------------------------------------------
% Bode plots: control-surface commands -> all modal coordinates xi_m
% Only relevant when isinp = 1
% -------------------------------------------------------------------------

if isinp == 1

    % Choose airspeed for Bode plot
    V_bode = 35;   % [m/s]
    [~, iV] = min(abs(airspeed - V_bode));

    Ap = models{iV}.Ap;
    Bp = models{iV}.Bp;

    Nctrl = size(Bp,2);     % number of control-surface command inputs

    fprintf('\n============================================================\n');
    fprintf('Bode model selected at V = %.2f m/s\n', models{iV}.V);
    fprintf('Bode model uses Ap because isinp = 1\n');
    fprintf('Number of control-surface command inputs = %d\n', Nctrl);
    fprintf('Number of modal coordinates = %d\n', Nel);
    fprintf('============================================================\n');

    % Frequency range [rad/s]
    w = logspace(-1, 3, 500);

    for mode_id = 1:Nel

        Cmodal = zeros(1, size(Ap,1));
        Cmodal(1, mode_id) = 1;

        sys_xi_mode = ss(Ap, Bp, Cmodal, zeros(1, Nctrl));

        for cs_id = 1:Nctrl

            figure;
            bode(sys_xi_mode(:, cs_id), w);
            grid on;

            title(sprintf('Bode: CS %d command to \\xi_%d, V = %.1f m/s', ...
                cs_id, mode_id, models{iV}.V));

        end
    end

else

    fprintf('\nBode plots skipped because isinp = 0.\n');
    fprintf('No control-surface command inputs are included, so Ap/Bp are not built.\n');

end
