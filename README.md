# ZLDD Cascade Solar Still — MATLAB Model

Transient multi-physics model of a zero-liquid-discharge desalination (ZLDD) system built around a
multi-plate cascade solar still coupled to a dehumidification/condensation unit (the atmospheric
water generator, AWG, of the patented design) and an air-to-air recuperator. The model solves the
coupled film, humid-air, glass-cover and absorber-plate balances
as a stiff system with `ode15s` and reports distillate yield, cascade recovery and specific energy
consumption (SEC). Local (OAT) and global (Morris) sensitivity harnesses are included.

## What the model includes

- Falling-film mass, energy and momentum balances on each cascade plate
- Humid-air transport in the inter-plate vapour gaps, with tapered gap geometry
- Glass cover, absorber plate, floor and side-wall energy balances
- Three-surface radiosity exchange per gap, with tapered-gap view factors from Hottel's flat-plate
  result and a crossed-strings correction
- Axial wall conduction between gap segments
- Segmented-duct fan pressure drop, dehumidification/condensation coil and air-to-air recuperator
- Closed-loop air circulation with salt mass-balance closure
- Optional entrance and turn factors on the gap heat- and mass-transfer coefficients
  (`FC.gap_factor`), used only for the CFD comparison case

## Requirements

MATLAB R2021b or newer. Base MATLAB only.

## Repository structure

| Path | Purpose |
|---|---|
| `code/main.m` | Driver. Run the code from here; it declares run controls only. |
| `code/Baseline_Input_Variable.m` | Design case, single point of declaration for all inputs |
| `code/CFD_Comparison_Input_Variable.m` | Input set of the CFD comparison case (Section 6 of the article) |
| `code/ZLDD_complete_modeling.m` | The model |
| `code/plot_baseline.m` | Results-section figures |
| `code/Sensitivity.m` | Local (OAT) sensitivity sweep and its post-processing |
| `code/Morris_Global.m` | Global (Morris) screening |
| `code/pvgis_irradiance_data.m` | Plane-of-array irradiance series |
| `data/OAT.mat` | Stored results of the full OAT sweep |
| `data/Global.mat` | Stored Morris screening results, r = 20 trajectories |
| `cfd/still2d_day_results.csv` | ANSYS Fluent model, comparison case: results over the operating window |
| `cfd/still2d_flow_results.csv` | ANSYS Fluent model, comparison case: air-flow solution (mass balance, stack pressure drop, fan air power) |
| `appendices/` | Supplementary document for the manuscript |

The two `.mat` files in `data/` are both results and checkpoints: `load` them to reproduce the
figures without re-solving, or leave them in place and the harnesses will resume from where they
stopped.

## Usage

Running `code/main.m` adds all repository folders to the MATLAB path. To use the functions directly
from the Command Window instead, first run this from the repository root folder:

```matlab
addpath(genpath(pwd));
```

Baseline solve on the Nx = 40 grid verified by the GCI study:

```matlab
FC      = Baseline_Input_Variable();
results = ZLDD_complete_modeling(FC);
h = plot_baseline(results);                    % display only
Sensitivity('export', 'svg', h);               % svg format save
```

The model draws the Results-section figures itself when `FC.quiet` is false.

CFD comparison case (Section 6 of the article). Same design case with a feed of 115 L per operating
window and entrance and turn factors on the gap transfer coefficients; the baseline input carries no
factors, so every result in Sections 4 and 5 uses the correlation unchanged. The Fluent results for
the same case are in `cfd/`.

```matlab
FC      = CFD_Comparison_Input_Variable();  % FC.gap_factor = [5.6 1.7 1.2]
results = ZLDD_complete_modeling(FC);
```

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

Research code supporting a manuscript under review. Interfaces may change.

## Citation

The underlying system is covered by US Patent 12,115,465 B1. If you use this code, please cite the
patent and the accompanying manuscript (details to follow on publication).
Link: https://patents.google.com/patent/US12115465B1/en


## Contact

- **Tuhin Mahmud**, MSc Student<br>
  Department of Chemical Engineering, Bangladesh University of Engineering and Technology (BUET), Dhaka, Bangladesh<br>
  Email: 1024022108@che.buet.ac.bd

- **Dr. Mohidus Samad Khan**, Professor<br>
  Department of Chemical Engineering, Bangladesh University of Engineering and Technology (BUET), Dhaka, Bangladesh<br>
  Email: mohid@che.buet.ac.bd

- **Dr. Sayeed Rushd**, Associate Professor<br>
  Department of Chemical Engineering, King Faisal University, Al-Ahsa, Saudi Arabia<br>
  Email: mrushd@kfu.edu.sa
