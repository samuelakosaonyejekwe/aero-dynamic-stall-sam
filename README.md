# Dynamic-Stall Prediction — Industrial Case Study
### UNISTALL™ Universal Unsteady-Aerodynamics & Dynamic-Stall Solver (UIBS core)

**Author:** Akosa Samuel Onyejekwe (independent)
**Date:** June 2026

---

## Overview

This repository contains a complete, reproducible engineering study predicting
**dynamic stall** on a helicopter main-rotor **retreating blade** (NACA 0012
section). The prediction is produced by a novel reduced-order universal solver
and validated against published NACA 0012 static and dynamic-stall experiments.

Retreating-blade dynamic stall sets the maximum forward speed and rotor thrust
of every conventional helicopter. In fast forward flight the retreating blade
(azimuth ψ ≈ 270°) is forced through a rapid once-per-revolution pitch-up at low
dynamic pressure and stalls dynamically: a leading-edge vortex sheds, lift
overshoots the static maximum, and a large nose-down moment break drives
vibration, control loads, fatigue and stall flutter. Accurate, low-cost
prediction of this cycle is essential for blade design and flight-envelope
expansion.

The full technical write-up is included as
[`aero_dynamic_stall_report.pdf`](aero_dynamic_stall_report.pdf).

---

## The solver in one line

UNISTALL marches a **Unified Indicial–Beddoes State-Space (UIBS)** model in
semichord time: attached-flow indicial loads + Kirchhoff trailing-edge
separation (two lags) + a leading-edge dynamic-stall vortex, with
compressibility corrections, an integrated potential-flow field reconstruction,
and a compressible thermal module. It is calibrated per case against a static
polar and validated against McAlister / Carr / McCroskey NACA 0012 data.

---

## Configurations solved

| Case | Description | Chord | Mach | Reduced freq. k | Incidence |
|---|---|---|---|---|---|
| **A — validation rig** | NACA 0012 oscillating aerofoil (matches McAlister/McCroskey deep dynamic-stall test point) | 0.30 m | 0.30 | 0.10 | 10° ± 10° |
| **B — application** | Retreating-blade section r/R = 0.75 (1/rev feathering) | 0.527 m | 0.28 | 0.074 | 12° ± 8° |

**Reference aircraft (generic medium utility helicopter):** 4-blade main rotor,
R = 8.18 m, blade chord 0.527 m, NACA 0012 section, tip speed ΩR ≈ 221 m/s,
advance ratio μ = 0.32, analysis station r/R = 0.75.

---

## Predicted outputs

Unsteady C_L / C_D / C_M hysteresis loops; dynamic-stall onset, lift overshoot
and moment break; trailing-edge separation history and dynamic-stall-vortex
trajectory; surface Cp and 2-D pressure / velocity / vorticity fields;
compressible static & recovery (skin) temperature fields; aerodynamic damping
(stall-flutter indicator); and sensitivity to mean incidence and reduced
frequency.

### Headline results

- Static-polar errors **< 1 %**.
- **All five** integral dynamic-stall metrics fall inside the published
  experimental envelope.
- Matched validation point (M = 0.30, k = 0.10, α = 10° ± 10°): dynamic
  C_L,max ≈ 2.07, C_M,c/4 break ≈ −0.31, C_D,max ≈ 0.31.

---

## Repository structure (pipeline order)

| Folder | Contents |
|---|---|
| `01_geometry/` | Airfoil geometry generation, coordinate CSVs, profile/thickness plots |
| `02_mesh/` | Body-fitted C-grid generation, mesh-quality metrics, mesh plots |
| `03_model_setup/` | Flow conditions, kinematics, thermo properties, solver config, static reference polar |
| `04_solver/` | `unistall_solver.py` (UIBS core + field reconstruction + thermal) and `run_case.py` |
| `05_solution/` | Time histories, Cp distributions, reconstructed fields, integral metrics, convergence residuals |
| `06_postprocessing/` | All plots (`plots/`) plus validation & calibration against experiment (`validation/`) |
| `08_engineering_drawings/` | Dimensioned 3-view, isometric, blade and section A-A drawings |
| `aero_dynamic_stall_report.pdf` | Consolidated technical report |

---

## Validation & calibration data sources

Sheldahl & Klimas SAND80-2114; Abbott & von Doenhoff (1959); McCroskey
NASA TM-100019 (static); McAlister / Carr / McCroskey NASA TP-1100 and
McCroskey et al. NASA TM-84245 (dynamic); Leishman (2006) and Leishman & Beddoes
(1989) for the model formulation. Full citations are given in the technical
report.

---

## Reproducing the pipeline

Each numbered stage is a self-contained Python script that consumes the outputs
of the previous stage. Run them in order, e.g.:

```bash
python3 01_geometry/generate_geometry.py
python3 02_mesh/generate_mesh.py
python3 03_model_setup/generate_setup.py
python3 04_solver/run_case.py
python3 06_postprocessing/make_all_plots.py
python3 06_postprocessing/make_3d_plots.py
python3 06_postprocessing/validation/validate.py
```

**Requirements:** Python 3.12+ with `numpy`, `scipy`, `matplotlib`, and
`pandas`.

---

## License & attribution

© 2026 Akosa Samuel Onyejekwe. Independent research and engineering work.

This repository — code, datasets, figures, engineering drawings, and technical
report — is licensed under the
[Creative Commons Attribution 4.0 International License (CC BY 4.0)](LICENSE).
You may share and adapt the material for any purpose, provided you give
appropriate credit to **Akosa Samuel Onyejekwe**, link to the license, and
indicate any changes.
