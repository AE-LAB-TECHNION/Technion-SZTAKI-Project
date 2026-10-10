# Aeroelastic RFA model with Dryden gust input

This repository builds a MATLAB aeroelastic plant from a ZAERO rational
function approximation (RFA). The plant includes actuator dynamics,
strain, displacement/acceleration, rate and RBM measurements, and a
three-state Dryden shaping filter.

The current model uses 20 structural aerodynamic lag states and 20 gust
aerodynamic lag states. With the 6 structural states, 12 actuator states
and 3 Dryden states, the augmented model contains 61 states.

## Run

Open MATLAB in the repository root and run:

```matlab
main
```

`main.m` adds the source folders to the MATLAB path, builds the model,
runs the selected analyses and displays the requested figures.

The workflow requires MATLAB with Control System Toolbox and Signal
Processing Toolbox.

## Repository structure

- `src/input` reads the RFA, modal and sensor data.
- `src/model` builds the aeroelastic, actuator and Dryden-augmented models.
- `src/gust` defines the MATLAB Dryden model and its frequency response.
- `src/responses` contains selectable Bode and step-response figures.
- `src/analysis` contains model checks and flutter reporting.
- `src/plotting` contains general plotting and figure-export functions.
- `data/model` contains the files required to build the MATLAB model.
- `zaero_model_comparison` contains optional MATLAB--ZAERO validation.

The files required for the standard MATLAB run are stored in `data/model`:

- the 20-lag OMCOS RFA file;
- the matching ZAERO modal output;
- the RBM state-space output;
- the F06 structural data and its sensor-matrix cache.

## Response figures

Select the required input-output plots in `default_analysis_config.m`:

```matlab
cfg.responsePlots.modalBode = false;
cfg.responsePlots.rbmBode = true;
cfg.responsePlots.stepRate = false;
cfg.responsePlots.stepDisp = false;
cfg.responsePlots.stepStrain = false;
cfg.responsePlots.stepRbm = false;
```

## Figures

Generated images are not stored in the repository. To save the figures
created by `main`, set:

```matlab
cfg.output.saveFigures = true;
```

The files will be written to `output/figures`.

## ZAERO comparison

The Dryden filter used by the augmented MATLAB model is defined entirely
from the parameters in `default_analysis_config.m`. ZAERO Dryden result
files are used only as validation references. See
`zaero_model_comparison/README.md` for the comparison workflow.

Select `cfg.gustModel = 'dryden'` or `'omcos'` in
`default_analysis_config.m`. ZAERO comparison figures can be disabled
without changing the MATLAB model.
