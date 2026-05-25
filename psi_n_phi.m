clear; clc;

f06file = 'model-0012.f06';

% לפי מה שעשית ב-ZAERO:
% mode 3 הוסר, לכן 3 המודים של ZAERO הם כנראה structural modes 1,2,4
modesToUse = [1 2 4];

% נקודות שבהן את רוצה displacement modes
dispGridIDs = [44 260 2691];

% רכיב ההזזה הרצוי:
% T1 = x, T2 = y, T3 = z
dispComp = 'T3';

% נקודות שבהן את רוצה rotational modes
rotGridIDs = [44 260 2691];

% רכיב הסיבוב הרצוי:
% R1 = rotation about x
% R2 = rotation about y
% R3 = rotation about z
rotComp = 'R2';

% אלמנטים של strain gauges / rods
strainElemIDs = [13915:13924 13935:13944 17293:17296 17530:17897];

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