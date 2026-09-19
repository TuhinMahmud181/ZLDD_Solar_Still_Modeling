# ZLDD Cascade Solar Still — MATLAB Model

Transient multi-physics model of a zero-liquid-discharge desalination (ZLDD) system built around a
multi-plate cascade solar still coupled to an atmospheric water generator and an air-to-air
recuperator. The model solves the coupled film, humid-air, glass-cover and absorber-plate balances
as a stiff system with `ode15s` and reports distillate yield, cascade recovery and specific energy
consumption (SEC). Local (OAT) and global (Morris) sensitivity harnesses are included.

## What the model includes

- Falling-film mass, energy and momentum balances on each cascade plate
- Humid-air transport in the inter-plate vapour gaps, with tapered gap geometry
- Glass cover, absorber plate, floor and side-wall energy balances
- Three-surface radiosity exchange per gap, with view factors by Gauss–Legendre quadrature
- Axial wall conduction between gap segments
- Segmented-duct fan pressure drop, AWG condenser coil and air-to-air recuperator
- Closed-loop air circulation with salt mass-balance closure

## Requirements

MATLAB R2021b or newer. Base MATLAB only.

## Files

| File | Purpose |
|---|---|
| `main.m` | Driver. Run the code from here; it declares run controls only. |
| `Baseline_Input_Variable.m` | Design case, single point of declaration for all inputs |
| `ZLDD_complete_modeling.m` | The model |
| `plot_baseline.m` | Results-section figures |
| `Sensitivity.m` | Local (OAT) sensitivity sweep and its post-processing |
| `Morris_Global.m` | Global (Morris) screening |
| `pvgis_irradiance_data.m` | Plane-of-array irradiance series |
| `OAT.mat` | Stored results of the full OAT sweep |
| `Global.mat` | Stored Morris screening results, r = 20 trajectories |

The two `.mat` files are both results and checkpoints: `load` them to reproduce the figures without
re-solving, or leave them in place and the harnesses will resume from where they stopped.

## Usage

Baseline solve on the Nx = 40 grid verified by the GCI study:

```matlab
FC      = Baseline_Input_Variable();
results = ZLDD_complete_modeling(FC);
h = plot_baseline(results);                    % display only
Sensitivity('export', 'svg', h);               % svg format save
```

The model draws the Results-section figures itself when `FC.quiet` is false.

OAT sweep. The sweep overrides the grid to Nx = 20 for run economy, so its results are read
relative to run 1 and are not directly comparable with the baseline figures.

```matlab
Sensitivity('dryrun', true);              % list runs, solve nothing
T = Sensitivity();                        % full sweep -> OAT.mat
Sensitivity('figs', T);                   % figure suite
Sensitivity('tornado', T, 'R_still_pct'); % tornado on any response
Sensitivity('export', 'svg');             % open figures -> vector
```

Morris screening. Gives what OAT cannot: sigma, the dependence of a factor's effect on where the
other factors sit. Cost is (k+1)·r solves.

```matlab
Morris_Global('dryrun', true);              % factor set and cost
M = Morris_Global('r', 20, 'timeout', 900); % ~12 h at r = 20 -> Global.mat
h = Morris_Global('figs', M);
```

Both harnesses checkpoint after every solve and resume from the existing `.mat` file, so an
interrupted run can be restarted without losing completed samples.

## Status

Research code supporting a manuscript in preparation. Interfaces may change.

## Citation

The underlying system is covered by US Patent 12,115,465 B1. If you use this code, please cite the
patent and the accompanying manuscript (details to follow on publication).
Link: https://patents.google.com/patent/US12115465B1/en

## Contact

Tuhin Mahmud
Department of Chemical Engineering, BUET
Email:1024022108@che.buet.ac.bd


