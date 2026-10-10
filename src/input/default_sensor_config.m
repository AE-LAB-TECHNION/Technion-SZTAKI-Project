function sensor = default_sensor_config()
%DEFAULT_SENSOR_CONFIG Sensor and strain-gauge definitions.

sensor.dispGridIDs = [44 260 2691];
sensor.dispComp = 'T3';      % T1/T2/T3

sensor.rotGridIDs = [44 260 2691];
sensor.rotComp = 'R2';       % R1/R2/R3

% Use 20 strain gauges: 10 on the front spar and 10 on the rear spar.
% The IDs were selected from the original 392 strain elements using the DAT
% element centroids, with approximately uniform spacing along the span.
sensor.frontSparStrainElemIDs = [ ...
    13940 17535 17590 17648 17666 ...
    17705 17723 17781 17820 17878];

sensor.rearSparStrainElemIDs = [ ...
    17562 17548 17605 17627 17689 ...
    17770 17752 17814 17855 17877];

sensor.strainElemIDs = [ ...
    sensor.frontSparStrainElemIDs, ...
    sensor.rearSparStrainElemIDs];

end
