# Dynamic-Stall Prediction — Industrial Case Study
### UNISTALL™ Universal Unsteady-Aerodynamics & Dynamic-Stall Solver (UIBS core)

**Author:** Akosa Samuel Onyejekwe (independent)
**Date:** 27 June 2026

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
and a compressible thermal module. The static separation law is recovered in
closed form by inverting the Kirchhoff relation against a single measured static
polar, so no loop-by-loop curve fitting is required; the model is then validated
against McAlister / Carr / McCroskey NACA 0012 data.

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

- Static-polar errors **< 1 %** (lift-curve slope 0.11 %, C_L,max 0.52 %, stall
  angle within the 1° resolution of the reference table — see
  `06_postprocessing/validation/validation_static.csv`). This is a *calibration
  check*: the separation law is fitted to that polar. The predictive claim is
  the held-out dynamic validation below.
- Case A, the matched validation point (M = 0.30, k = 0.10, α = 10° ± 10°):
  dynamic C_L,max = 1.912 at α = 17.5°,
  C_M,c/4 break = −0.236, C_D,max = 0.257,
  a 34 % overshoot above the static maximum
  (`05_solution/metrics_A_validation.csv`).
- Case B, the retreating-blade station (M = 0.28, k = 0.074, α = 12° ± 8°):
  C_L,max = 1.742 at α = 15.8°,
  C_M,c/4 break = −0.209
  (`05_solution/metrics_B_application.csv`).
- Aerodynamic damping: both cases give a *figure-of-eight* C_M loop whose two
  lobes very nearly cancel. The normalised damping Ξ̂ = Ξ / (ΔC_M · Δα) is
  -0.008 (Case A) and -0.0013 (Case B) — negative,
  but well inside the ±0.02 band where the residual is no larger than the
  time-step discretisation error. Both are therefore reported as **neutrally
  damped**, not as a positive stall-flutter finding.

### Held-out validation

The dynamic constants were calibrated against **one** measured oscillating-aerofoil
loop (frame 9302), frozen, and then used to predict four entirely held-out loops
spanning light to deep stall. Across those four blind predictions the frozen model
returns a **mean peak-lift error of 1.7 %** and a **mean moment-break error of
0.023** (mean RMS over the whole C_L loop is 0.195, i.e. the
integral peaks are matched considerably better than the full loop shape —
the residual is dominated by the downstroke/reattachment branch).
Applying the same frozen constants to a *different* aerofoil section
(frame 25104) degrades peak-lift agreement to about 10 %, confirming that the
calibration encodes section-specific physics rather than a generic loop shape.
Per-frame figures are in `06_postprocessing/validation/validation_nasa_real.csv`.

---

## Repository structure (pipeline order)

| Folder | Contents |
|---|---|
| `01_geometry/` | Airfoil geometry generation, coordinate CSVs, profile/thickness plots |
| `02_mesh/` | Body-fitted O-grid generation, mesh-quality metrics, mesh plots |
| `03_model_setup/` | Flow conditions, kinematics, thermo properties, solver config, static reference polar |
| `04_solver/` | `unistall_solver.py` (UIBS core + field reconstruction + thermal) and `run_case.py` |
| `05_solution/` | Time histories, Cp distributions, reconstructed fields, integral metrics, convergence residuals |
| `06_postprocessing/` | All plots (`plots/`) plus validation & calibration against experiment (`validation/`) |
| `08_engineering_drawings/` | Dimensioned 3-view, isometric, blade and section A-A drawings |
| `aero_dynamic_stall_report.pdf` | Consolidated technical report — built by the pipeline (report body + data dossier + plots album), not exported by hand |

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
of the previous stage. Every stage imports the shared plotting module
`aero_style.py` at the repository root, so run them from a full checkout:

```bash
pip install -r requirements.txt

python3 01_geometry/generate_geometry.py
python3 02_mesh/generate_mesh.py
python3 03_model_setup/generate_setup.py
python3 04_solver/run_case.py
python3 06_postprocessing/make_all_plots.py
python3 06_postprocessing/make_3d_plots.py
python3 06_postprocessing/validation/validate.py            # static calibration check
python3 06_postprocessing/validation/validate_nasa_real.py  # held-out dynamic validation
python3 06_postprocessing/validation/validate_digitized.py  # certification harness (optional)
```

Every number, figure and table above is produced by these stages — nothing is
transcribed by hand. `aero_dynamic_stall_report.pdf` is likewise a generated
artefact (report body + data dossier + plots album, assembled from the same
outputs), so it cannot fall out of step with the solver. The report-assembly
tooling itself is not distributed here; the report is included as the finished
PDF.

`03_model_setup/` is the single source of truth for the case conditions: the
solver and every plotting script read `flow_conditions.csv`, `kinematics.csv`
and `solver_config.json` rather than restating any value.

**Requirements:** Python 3.9+ with `numpy`, `scipy`, `matplotlib` and `pandas`
(see `requirements.txt`). Developed and regenerated on Python 3.12.

---

## Related publication

A manuscript based on this study, *"A Rotor Blade Dynamic Stall Model with Closed Form
Separation Calibration,"* is in preparation for submission to the
*Journal of the American Helicopter Society*. The manuscript itself is not
distributed in this repository. The consolidated technical report
(`aero_dynamic_stall_report.pdf`) remains the full write-up of the study.

---

## License & attribution

© 2026 Akosa Samuel Onyejekwe. Independent research and engineering work.

This repository — code, datasets, figures, engineering drawings, and technical
report — is licensed under the
[Creative Commons Attribution 4.0 International License (CC BY 4.0)](LICENSE).
You may share and adapt the material for any purpose, provided you give
appropriate credit to **Akosa Samuel Onyejekwe**, link to the license, and
indicate any changes.
