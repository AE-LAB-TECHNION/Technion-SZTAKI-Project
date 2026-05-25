function sensor = default_sensor_config()
%DEFAULT_SENSOR_CONFIG Sensor and strain-gauge definitions.

sensor.dispGridIDs = [44 260 2691];
sensor.dispComp = 'T3';      % T1/T2/T3

sensor.rotGridIDs = [44 260 2691];
sensor.rotComp = 'R2';       % R1/R2/R3

sensor.strainElemIDs = [13915:13924, 13935:13944, 17293:17296, 17530:17897];

end
