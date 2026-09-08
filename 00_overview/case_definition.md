# Case Definition

**Author:** Akosa Samuel Onyejekwe (independent)

## Title
Prediction of dynamic stall on a helicopter main-rotor retreating blade using
the UNISTALL™ universal solver (UIBS core).

## Industrial problem
Retreating-blade dynamic stall sets the maximum forward speed and rotor thrust
of every conventional helicopter. In fast forward flight the retreating blade
(azimuth ψ ≈ 270°) is forced through a rapid once-per-revolution pitch-up at low
dynamic pressure and stalls dynamically: a leading-edge vortex sheds, lift
overshoots the static maximum, and a large nose-down moment break drives
vibration, control loads, fatigue and stall flutter. Accurate, cheap prediction
of this cycle is required for blade design and flight-envelope expansion.

## Reference aircraft (CS-MUH)
Generic medium utility helicopter: 4-blade main rotor, R = 8.18 m, blade chord
0.527 m, NACA 0012 section, Ω = 27.0 rad/s so tip speed
ΩR = 220.9 m/s, advance ratio μ = 0.32.
Analysis station r/R = 0.75.

## Configurations solved
Only independent quantities are specified; a, ρ, U (case B) and Re_c are derived
in `03_model_setup/generate_setup.py`, which is the single source of truth.

- **Case A — validation rig:** NACA 0012 oscillating aerofoil, c = 0.30 m,
  M = 0.30, k = 0.10, α = 10° ± 10° (matches McAlister/McCroskey deep dynamic
  stall test point). U = M·a = 102.09 m/s.
- **Case B — application:** retreating-blade section r/R = 0.75, c = 0.527 m,
  U = ΩR(r/R − μ) = 94.97 m/s so M = 0.2794, k = 0.074,
  α = 12° ± 8° (1/rev feathering).

## Required predictions (outputs)
Unsteady C_L/C_D/C_M hysteresis loops; dynamic-stall onset, lift overshoot and
moment break; trailing-edge separation history and dynamic-stall-vortex
trajectory; surface Cp and 2-D pressure/velocity/vorticity fields; compressible
static & recovery (skin) temperature fields; aerodynamic damping (stall-flutter
indicator); and sensitivity to mean incidence and reduced frequency.

## Validation & calibration data sources
Sheldahl & Klimas SAND80-2114; Abbott & von Doenhoff (1959); McCroskey NASA
TM-100019 (static); McAlister/Carr/McCroskey NASA TP-1100 and McCroskey et al.
NASA TM-84245 (dynamic); Leishman (2006) and Leishman & Beddoes (1989) for the
model formulation. Full citations in `aero_dynamic_stall_report.pdf` §16.
