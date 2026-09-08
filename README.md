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
| **B — application** | Retreating-blade section r/R = 0.75 (1/rev feathering) | 0.527 m | 0.279 | 0.074 | 12° ± 8° |

**Reference aircraft (generic medium utility helicopter):** 4-blade main rotor,
R = 8.18 m, blade chord 0.527 m, NACA 0012
section, Ω = 27.0 rad/s so tip speed
ΩR = 220.9 m/s, advance ratio μ = 0.32,
analysis station r/R = 0.75.

The free-stream state is over-determined — (p, T, ρ), (M, U, a) and, for Case B,
the rotor kinematics all describe the same flow — so only the independent
quantities are specified and the rest are derived: a = √(γRT), ρ = p/(RT),
Re_c = ρUc/µ, with U = M·a for the rig and U = ΩR(r/R − μ) = 94.97 m/s
for the blade station (which is why Case B's Mach is 0.2794, not a
round 0.28). Nothing is quoted twice.

---

## Predicted outputs

Unsteady C_L / C_D / C_M hysteresis loops; dynamic-stall onset, lift overshoot
and moment break; trailing-edge separation history and dynamic-stall-vortex
trajectory; surface Cp and 2-D pressure / velocity / vorticity fields;
compressible static & recovery (skin) temperature fields; aerodynamic damping
(stall-flutter indicator); and sensitivity to mean incidence and reduced
frequency.

The loads come from the UIBS core. The 2-D fields and surface Cp come from a
separate potential-flow reconstruction driven by the UIBS circulation.
Integrating the reconstructed surface Cp returns the C_L it was given to within
**0.9 %** — 0.9 % (Case A) and 0.9 % (Case B) — a number the pipeline measures on every run and
publishes as `Cp_closure_error_pct` in `05_solution/metrics_*.csv`, rather than
a claim in a comment. It previously read -12.4 %, and the explanation recorded
for that deficit was itself wrong; the three real causes (a missing vortex-sheet
self-term, evaluation off the wall instead of at the panel control points, and a
double-counted compressibility factor) are documented in `unistall_solver.py`.

Two limitations of the reconstruction are published as numbers rather than
described in prose. `Cp_TE_jump_max_over_phases` is the trailing-edge pressure
jump; it cannot be driven to zero, because the reconstruction is handed the
indicial C_L, which during dynamic stall departs deliberately from the inviscid
attached circulation that the Kutta condition selects. That circulation is
published beside it as `CL_kutta_inviscid`, and imposing it drives the jump to
~0.001. `Cp_DSV_core_min` is the suction at the centre of the reconstructed
dynamic-stall vortex; it reads about -0.4 where a measured deep-stall core is
usually nearer -3 to -6, and it is reported rather than tuned because nothing
this study ships (integrated cl/cd/cm only) could calibrate a core size.

### Headline results

- Static-polar errors **< 1 %** (lift-curve slope 0.11 %, C_L,max 0.52 %, stall
  angle within the 1° resolution of the reference table — see
  `06_postprocessing/validation/validation_static.csv`). This is a *calibration
  check*: the separation law is fitted to that polar. The predictive claim is
  the held-out dynamic validation below.
- Case A, the matched validation point (M = 0.30, k = 0.10, α = 10° ± 10°):
  dynamic C_L,max = 1.912 at α = 17.5°,
  C_M,c/4 break = -0.236, C_D,max = 0.257,
  dynamic-stall onset at α = 12.50° (the incidence at which the
  model's own vortex-shedding switch fires, C_N′ ≥ C_N1 on the upstroke),
  a 34 % overshoot above the static maximum
  (`05_solution/metrics_A_validation.csv`).
- Case B, the retreating-blade station (M = 0.279, k = 0.074, α = 12° ± 8°):
  C_L,max = 1.742 at α = 15.8°,
  C_M,c/4 break = -0.209, onset at α = 12.49°
  (`05_solution/metrics_B_application.csv`).
- Aerodynamic damping: both cases give a *figure-of-eight* C_M loop whose two
  lobes very nearly cancel. The normalised damping Ξ̂ = Ξ / (ΔC_M · Δα) is
  -0.008 (Case A) and -0.0007 (Case B) — negative, but far inside the
  ±0.08 band within which the model cannot resolve the sign. Both are
  therefore reported as **neutrally damped**, not as a stall-flutter finding.
  That band is *measured*, not assumed: recomputing Ξ̂ from the real NACA 0012
  C_M loops and from the model at the same conditions gives a mean discrepancy
  of 0.072 (max 0.144) — see
  `06_postprocessing/validation/validation_nasa_real.csv`. It is *not* the
  time-step error, which is about 45× smaller (refining 720 → 5760 steps per
  cycle moves Ξ̂ by 0.0004).

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

## Repository structure

Folder numbers are sections, not execution order: `03_model_setup/` runs first
because it defines the case conditions that `01_geometry/` and `02_mesh/`
consume. Execution order is given under *Reproducing the pipeline* below.

| Folder | Contents |
|---|---|
| `03_model_setup/` | Flow conditions, kinematics, thermo properties, solver config, static reference polar — **runs first** |
| `01_geometry/` | Airfoil geometry generation, coordinate CSVs, profile/thickness plots |
| `02_mesh/` | Body-fitted O-grid generation, mesh-quality metrics, mesh plots |
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

Each stage is a self-contained Python script that consumes the outputs of the
stages before it, in the order shown. The plotting stages import the shared
style module `aero_style.py` at the repository root, so run them from a full
checkout (`04_solver/` and `03_model_setup/` do not need it; every other stage
does):

```bash
pip install -r requirements.txt

# NOTE the order: 03_model_setup runs FIRST. It is the single source of truth for
# the case conditions, and both 01_geometry (chord) and 02_mesh (chord, rho, U,
# mu) consume them, so it cannot run after them.
python3 03_model_setup/generate_setup.py
python3 01_geometry/generate_geometry.py
python3 02_mesh/generate_mesh.py
python3 04_solver/run_case.py
python3 06_postprocessing/make_all_plots.py
python3 06_postprocessing/make_3d_plots.py
python3 06_postprocessing/validation/validate.py            # static calibration check
python3 06_postprocessing/validation/validate_nasa_real.py  # held-out dynamic validation
python3 06_postprocessing/validation/validate_digitized.py  # certification harness (optional)

python3 08_engineering_drawings/draw_engineering.py         # the four drawing sheets

# Report build. These three need python-docx, pillow, reportlab and PyMuPDF,
# which requirements.txt installs but the solver pipeline itself does not need.
python3 07_report/build_docx.py           # assembles case.docx
python3 07_report/build_pdfs.py           # plots album + data dossier
python3 07_report/build_report_pdf.py     # -> aero_dynamic_stall_report.pdf

python3 check_claims.py    # asserts every number quoted below still matches the CSVs
```

Or run all fourteen stages in the correct order with a single command:

```bash
python3 run_all.py
```

Every number, figure and table above is produced by these stages. The handful
that are necessarily transcribed — the results quoted in this README and in
`00_overview/case_definition.md`, because Markdown cannot compute — are checked
against the generated CSVs by `check_claims.py`, which the pipeline runs last
and which fails the build if any of them has drifted. `aero_dynamic_stall_report.pdf` is likewise a generated
artefact (report body + data dossier + plots album, assembled from the same
outputs), so it cannot fall out of step with the solver. The report-assembly
tooling in `07_report/` is distributed too, so a clean checkout can rebuild the
report from the solver outputs rather than having to trust the shipped PDF. The
intermediate products of that build -- the plots album, the data dossier,
`case.docx` and the rendered equation images -- are not committed, because they
are regenerated from the same outputs on every run.

`03_model_setup/` is the single source of truth for the case conditions.
`generate_geometry.py`, `generate_mesh.py`, `run_case.py`, `make_all_plots.py`,
`make_3d_plots.py` and the validation scripts all read `flow_conditions.csv`,
`kinematics.csv`, `material_thermo_properties.csv` and `solver_config.json`
rather than restating any value — which is why the setup stage runs first.

**Requirements:** Python 3.9+ with `numpy`, `scipy`, `matplotlib` and `pandas`
(see `requirements.txt`). Developed and regenerated on Python 3.12.

The report stages in `07_report/` additionally need `python-docx`, `pillow`,
`reportlab` and `PyMuPDF`, plus the DejaVu Sans Condensed fonts for the Greek
and symbol glyphs (`fonts-dejavu-core` on Debian/Ubuntu, `dejavu-sans-fonts` on
Fedora, `ttf-dejavu` on Arch). The build searches the usual font directories and
exits with an install hint if it cannot find them. The solver pipeline itself
needs none of this.

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
