# Case definition

**Author:** Akosa Samuel Onyejekwe (independent)

This file is written by `03_model_setup/generate_setup.py`; every number in it
is read from the tables of `03_model_setup`.

## Title
Prediction of dynamic stall on a pitching NACA 0012 section with UNISTALL, a
Leishman–Beddoes dynamic-stall load model with tabulated separation (Python
package `unistall`).

## The problem
On a helicopter in fast forward flight the retreating blade is pitched up
through stall once per revolution at low dynamic pressure. A leading-edge
vortex forms and is shed, lift overshoots its static maximum, and a large
nose-down pitching moment follows. That cycle limits rotor thrust and speed
and drives vibration and control loads, so a cheap prediction of the unsteady
lift, drag and moment on the blade section is wanted in design.

## What is solved
The unsteady loads on a NACA 0012 section in a steady stream with prescribed
sinusoidal pitch, α(t) = α_mean + α_amp sin(ωt), by the load model in
`unistall/dsmodel.py`. Two conditions:

| | Case A: tunnel condition | Case B: rotor-blade-station condition |
|---|---|---|
| Purpose | a condition of the NASA TM-84245 experiment, so the result can be set beside a measured loop (frame 9302, a calibration loop) | an illustration at the condition of a blade section at 0.75 of the radius |
| Mach number | 0.30 | 0.28 |
| Reduced frequency k | 0.1000 | 0.0747 |
| Pitch frequency | 5.327 Hz | 4.300 Hz |
| Incidence | 10° ± 10° | 12° ± 8° |
| Chord | 0.61 m | 0.527 m |
| Chord Reynolds number in sea-level standard air (not the tunnel value; not used by the model) | 4.26e+06 | 3.44e+06 |

Both lie inside the range the model has been compared with measurement
(reduced frequency up to 0.2, peak incidence up to 25°, at the Mach
numbers of the static data).

**Case A is set beside a calibration loop.** The measured loop drawn with it
(frame 9302) was used to fit the model, so the comparison is an
illustration and not an independent test.

**Case B is not a rotor calculation.** Its chord and its once-per-revolution
pitch frequency are those of the UH-60A main rotor (4 blades, radius
26.83 ft, chord 20.76 in, 258 rpm; Bousman, NASA/TP-2003-212265).
The station at 0.75 of the radius turns at 165.7 m/s. Case B
takes the Mach number 0.28 for it, which is what the station sees at
azimuth 270° (the retreating side) at a forward speed of
70.4 m/s (137 kt, advance ratio 0.32); the reduced
frequency 0.0747 follows from the rotor speed, the chord and that speed
(`03_model_setup/station_condition.csv`). Round the azimuth the Mach number at
that station runs from 0.28 to 0.69; Case B holds the lowest value
steady, and its incidence history (12° ± 8°, sinusoidal) is assumed: it
comes from no trim calculation. The real blade uses the SC1095 and SC1094 R8
sections with a swept tip, not a NACA 0012. There is no time-varying
velocity, sweep, inflow, trim or blade motion, and no rotor measurement is
compared with.

**The static inputs carry the offset of the measured static sweeps.** The
section is symmetric, yet the static inputs read from the experiment give
Case A: zero-lift incidence +0.23°, static normal force at zero incidence -0.027; Case B: zero-lift incidence +0.18°, static normal force at zero incidence -0.021
(`03_model_setup/solver_config.json`). The offset is that of the measured
sweeps (a flow angle of the tunnel or an offset of the balance), not a
property of the section, and the model keeps it so that it reads the static
data as measured. The incidence uncertainty the experiment states is not held
in this repository's data files, so the offset is not set against it here.

## What each part of the case study is, and is not

| Folder | Contents | Status |
|---|---|---|
| `01_geometry` | the section's coordinates and thickness | the NACA four-digit thickness form with a closed trailing edge; the same outline in the grid, the panels and the drawings |
| `02_mesh` | a structured O-grid round the section | an illustration: the load model has no spatial mesh and does not read it; no equation of the flow is solved on it; the reconstructed field is evaluated at its nodes for one picture |
| `03_model_setup` | flow conditions, kinematics, air properties, solver configuration, static inputs | derived by script; the constants are the calibrated ones |
| `04_solver` | the driver for the two cases | the solver itself is the `unistall` package: a set of ordinary differential equations marched in time |
| `05_solution` | time histories, surface pressures, reconstructed fields, metrics, convergence | loads are the model's prediction; fields and surface pressures are a reconstruction (below), not a solution of the flow equations |
| `06_postprocessing` | load loops, maps of the reconstructed field | drawn from `05_solution`; the field maps show the reconstruction |
| `07_report` | the case-study report, a plots album and a data dossier | built from the files above |
| `08_engineering_drawings` | four illustrative sheets | every dimension from one table with its source; not for manufacture |

The folder names follow the layout of a flow-solver case. The calculation is
not one: it has no mesh and solves no flow equations. In a paper on the load
model the grid, the field reconstruction and drawing sheets 1 and 2 would be
left out; they are kept in this repository as illustrations, each with its
standing stated where it appears.

**The flow fields are a reconstruction, not a flow solution.** Potential-flow reconstruction drawn round the predicted lift, not a flow solution: incompressible, with a linearised (Karman-Tsien) correction of the pressure that holds only where the local flow is subcritical, and quasi-steady (no shed wake, no pitch-rate boundary condition, no unsteady pressure term, vortex not force-free); no boundary layer, no separated shear layer. A transpiration through the suction surface, of either sign, stands for the separated region, which is not drawn: no reversed flow and no wake appear. The vortex is a marker of the model's vortex clock and does not reproduce a suction footprint. Compared with no measurement; no load depends on it.
Its surface pressure is given the lift only and does not return the model's
moment or drag.

## Accuracy
How well the load model agrees with measurement is established separately, on
held-out loops, and is reported in the repository's `README.md` and
`results/tables.md`, with the targets it meets and the ones it misses. In
particular the model does not predict cycle damping through stall, so the
damping figures of the two cases here are reported as model output and not as
a statement about stall flutter.

## Sources
McCroskey, McAlister, Carr & Pucci (1982), NASA TM-84245 Vol. 1; McAlister, Pucci, McCroskey & Carr (1982), NASA TM-84245 Vol. 2 (measured loops and static data); Leishman & Beddoes
(1989), J. Am. Helicopter Soc. 34(3); Damiani & Hayman (2019),
NREL/TP-5000-66347 (formulation); Bousman, NASA/TP-2003-212265 (rotor parameters of
Case B; not an entry of `docs/references.bib`, and the four parameters have
not been checked against its pages in this repository). The other sources are
listed in full in `docs/references.bib`.
