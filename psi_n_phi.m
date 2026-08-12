clear; clc;

cfg = default_analysis_config();
f06file = cfg.f06file;

% structural modes 1,2,4
modesToUse = [1 2 4];

% displacement modes
dispGridIDs = [44 260 2691];

% T1 = x, T2 = y, T3 = z
dispComp = 'T3';

% rotational modes
rotGridIDs = [44 260 2691];

% R1 = rotation about x
% R2 = rotation about y
% R3 = rotation about z
rotComp = 'R2';

% strain gauges / rods
strainElemIDs = [13940 17535 17590 17648 17666 17705 17723 17781 17820 17878 17562 17548 17605 17627 17689 17770 17752 17814 17855 17877];

[PSI, PHI, PHI_ROT, info] = extract_PSI_PHI_PHIROT_from_F06( ...
    f06file, modesToUse, ...
    dispGridIDs, dispComp, ...
    rotGridIDs, rotComp, ...
    strainElemIDs);

save('sensor_modal_matrices.mat','PSI','PHI','PHI_ROT','info');

disp('PSI size:')
disp(size(PSI))

disp('PHI size:')
disp(size(PHI))

disp('PHI_ROT size:')
disp(size(PHI_ROT))
