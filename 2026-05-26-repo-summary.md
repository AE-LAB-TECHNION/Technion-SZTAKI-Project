# Repository Summary - 2026-05-26

This repository is a compact MATLAB aeroelastic and flutter analysis workflow built around ZAERO and F06 output files.

## Main Workflow

The main entry point is `main.m`. It:

1. Builds the analysis configuration with `default_analysis_config.m`.
2. Loads and derives analysis inputs with `load_analysis_inputs.m`.
3. Sweeps airspeed with `run_airspeed_sweep.m`.
4. Plots root-locus and omega-V-g results.
5. Reports the estimated flutter point.
6. Generates control Bode plots.

The default model type is currently the plant model, not the bare aeroelastic model:

```matlab
cfg.modelType = 'plant';
```

## Core Numerical Path

- `extractRFAmatrices.m` reads ZAERO RFA matrices from `APPROX.DAT`.
- `read_zaero_modal_info.m` parses modal metadata from `ASE_ANALYSIS.out`, including retained/omitted modes, control-surface count, gust information, and eigenvalues.
- `extract_PSI_PHI_PHIROT_from_F06.m` extracts displacement, acceleration, rotation, and strain-related modal matrices from `model-0012.f06`.
- `buildAESS_state.m` builds the aeroelastic state-space matrices: `Aae`, `Bae`, `Baw`, `Cae`, `Caw`, and `Dae`.
- `buildPlant_from_AESS.m` augments the aeroelastic model with third-order actuator dynamics to form the plant model.

## Analysis And Plotting Layer

- `plot_root_locus.m` plots eigenvalue trajectories over airspeed.
- `omega_v_g_data.m` tracks physical modal branches and computes modal frequency and damping.
- `plot_omega_v_g.m` creates omega-V-g plots.
- `find_flutter.m` and `report_flutter.m` estimate and print flutter speed and flutter frequency.
- `plot_control_bodes.m` plots control-surface or actuator-command Bode responses.
- `apply_plot_style.m` applies shared per-figure plotting style.

## Important Data Files

- `APPROX.DAT`: RFA matrix data.
- `ASE_ANALYSIS.out`: current ZAERO output used as the modal metadata source.
- `ASE_ANALYSIS_new.out`: older or alternate ZAERO output.
- `model-0012.f06`: large F06 modal/source data file.
- `sensor_modal_matrices.mat`: cached sensor/modal matrices.

## Current Design Direction

The repo has recently moved toward deriving modal, control, and gust dimensions from ZAERO output instead of hard-coding selected modes. The current path is more self-configuring:

- ZAERO metadata determines retained FEM modes.
- Omitted modes are parsed from the ZAERO output.
- Control-surface and gust partitions are inferred from ZAERO/RFA data.
- Sensor/modal matrices can be loaded from cache when compatible.

This keeps the numerical workflow tied to the current ZAERO run files while reducing manual mode-selection assumptions.

## Validation Note

This summary is based on static repository inspection and earlier workflow review. The full `main.m` aeroelastic analysis was not run as part of generating this summary.
