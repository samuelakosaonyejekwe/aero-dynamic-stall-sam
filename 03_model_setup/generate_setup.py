# -*- coding: utf-8 -*-
# Run from the repository root:  PYTHONPATH=. python3 03_model_setup/generate_setup.py
"""
03_model_setup / generate_setup.py
----------------------------------
The conditions of the two cases, written out once. Every later stage (the
solution, the plots, the drawings, the report, the checks) reads them from
the files written here, through project_meta.read_setup.

Author: Akosa Samuel Onyejekwe (independent)

Case A is a condition of the NASA TM-84245 experiment. Case B is an aerofoil
at a rotor-blade-station condition (steady stream, prescribed pitch); it is
not a rotor calculation. Its chord, Mach number and reduced frequency are
derived from the rotor parameters in project_meta (station_condition), and
the derivation is written to station_condition.csv.

Writes, beside this file:
  flow_conditions.csv     chord [m], Mach number, stream speed [m/s], density
                          [kg/m^3], pressure [Pa], temperature [K], viscosity
                          [Pa s], chord Reynolds number, reduced frequency,
                          each with its source
  kinematics.csv          mean incidence and amplitude [deg], reduced frequency,
                          circular frequency [rad/s], frequency [Hz], period [s]
  station_condition.csv   how the Mach number and reduced frequency of Case B
                          follow from the rotor, and the range of Mach number
                          round the azimuth that the steady stream leaves out
  air_properties.csv      the gas properties the cases are run with
  solver_config.json      solver name, calibrated constants, resolution of the
                          march and critical pressure coefficient for each
                          case, settings of the reconstruction
  static_inputs.csv       static normal force, moment and separation point the
                          solver uses at each case's Mach number
and ../00_overview/case_definition.md, the definition of the study in words,
with every number in it taken from the tables above.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

import project_meta as pm
from unistall import dsmodel, flowfield
from unistall.static_model import ALPHA_MAX, StaticModel

HERE = Path(__file__).resolve().parent
OVERVIEW = HERE.parent/"00_overview"
AIR = flowfield.Air(gamma=pm.GAMMA, R_gas=pm.R_GAS, T_inf=pm.T_INF, p_inf=pm.P_INF)
STATIC_STEP_DEG = 0.25               # spacing of static_inputs.csv
MS_TO_KT = 3600.0/1852.0
CASE_B_NOTE = "an aerofoil at a rotor-blade-station condition (steady stream, prescribed pitch); not a rotor calculation"
STANDARD_AIR = "sea-level standard air"
FIELD_PHASES = {"rise": "upstroke, mean incidence plus half the amplitude",
                "peak": "instant of greatest lift",
                "dsv": "instant of greatest vortex normal force",
                "fall": "downstroke, mean incidence plus half the amplitude"}


def description(name: str, case: dict) -> str:
    """One line saying what a case is."""
    if case["measured_frame"]:
        return "tunnel condition of NASA TM-84245, " + case["measured_frame"].replace("_", " ")
    return CASE_B_NOTE


def flow_conditions() -> pd.DataFrame:
    """One column per case, SI units, and the source of each row."""
    derived_b = "station_condition.csv" if pm.CASE_B in pm.CASES else ""
    spec = [("description", "-", lambda n, c, U: description(n, c), "project_meta.CASES"),
            ("airfoil", "-", lambda n, c, U: "NACA " + pm.SECTION, "project_meta.SECTION"),
            ("chord_c", "m", lambda n, c, U: c["chord_m"],
             f"Case A: NASA TM-84245; Case B: {pm.ROTOR_SOURCE}"),
            ("freestream_mach_M", "-", lambda n, c, U: c["M"],
             f"Case A: nominal condition of the frame; Case B: {derived_b}"),
            ("speed_of_sound_a", "m/s", lambda n, c, U: repr(AIR.a_sound), STANDARD_AIR + ", sqrt(gamma R T)"),
            ("freestream_velocity_U", "m/s", lambda n, c, U: repr(U), "M a"),
            ("air_density_rho", "kg/m^3", lambda n, c, U: round(AIR.rho, 5), STANDARD_AIR),
            ("static_pressure_p_inf", "Pa", lambda n, c, U: AIR.p_inf, STANDARD_AIR),
            ("static_temperature_T_inf", "K", lambda n, c, U: AIR.T_inf, STANDARD_AIR),
            ("dynamic_viscosity_mu", "Pa.s", lambda n, c, U: float("%.5g" % AIR.mu), STANDARD_AIR + ", Sutherland"),
            ("reynolds_number_Re_c", "-", lambda n, c, U: float("%.4g" % (AIR.rho*U*c["chord_m"]/AIR.mu)),
             STANDARD_AIR + "; not the tunnel value of the frame and not a flight value; the load model does "
             "not use it"),
            ("reduced_frequency_k", "-", lambda n, c, U: repr(c["k"]),
             f"Case A: nominal condition of the frame; Case B: {derived_b}")]
    table = {"parameter": [row[0] for row in spec]}
    for name, case in pm.CASES.items():
        U = case["M"]*AIR.a_sound
        table[name] = [row[2](name, case, U) for row in spec]
    table["units"] = [row[1] for row in spec]
    table["source"] = [row[3] for row in spec]
    return pd.DataFrame(table)


def kinematics() -> pd.DataFrame:
    """The prescribed pitching motion of each case."""
    rows = []
    for name, c in pm.CASES.items():
        omega = 2.0*c["k"]*c["M"]*AIR.a_sound/c["chord_m"]
        rows.append(dict(case=name, alpha_mean_deg=c["alpha_mean_deg"], alpha_amp_deg=c["alpha_amp_deg"],
                         reduced_freq_k=repr(c["k"]), omega_rad_s=round(omega, 6),
                         freq_Hz=round(omega/(2.0*np.pi), 6), period_s=round(2.0*np.pi/omega, 6),
                         motion="alpha(t) = mean + amp sin(omega t), about the quarter chord"))
    return pd.DataFrame(rows)


def station_condition() -> pd.DataFrame:
    """How the condition of Case B follows from the rotor, row by row."""
    s, r = pm.station_condition(), pm.ROTOR
    src = pm.ROTOR_SOURCE
    psi = pm.STATION_AZIMUTH_DEG
    rows = [
        ("rotor_radius", r["radius_ft"], "ft", src),
        ("rotor_radius_m", round(s["radius_m"], 4), "m", "radius in feet x 0.3048"),
        ("blade_chord", r["chord_in"], "in", src),
        ("blade_chord_m", s["chord_m"], "m", "chord in inches x 0.0254, to three decimals"),
        ("number_of_blades", r["blades"], "-", src),
        ("rotor_speed", r["speed_rpm"], "rpm", src),
        ("rotor_speed_rev_per_s", round(s["rev_per_s"], 6), "1/s", "rpm / 60"),
        ("rotor_speed_Omega", round(s["omega_rad_s"], 5), "rad/s", "2 pi rpm / 60"),
        ("tip_speed", round(s["tip_speed_ms"], 3), "m/s", "Omega R"),
        ("station_r_over_R", pm.STATION_OVER_R, "-", "chosen for this study"),
        ("station_rotational_speed", round(s["rotational_speed_ms"], 3), "m/s", "Omega r"),
        ("station_azimuth", psi, "deg", "chosen for this study: the retreating side"),
        ("station_mach_M", pm.STATION_MACH, "-", "chosen for this study: the Mach number of Case B"),
        ("speed_of_sound", round(s["speed_of_sound_ms"], 3), "m/s", STANDARD_AIR + ", sqrt(gamma R T)"),
        ("station_speed_U", round(s["station_speed_ms"], 3), "m/s", "M a"),
        ("forward_speed_implied", round(s["forward_speed_ms"], 2), "m/s",
         "the forward speed V at which Omega r + V sin(azimuth) = U"),
        ("forward_speed_implied_kt", round(s["forward_speed_ms"]*MS_TO_KT, 1), "kt", "V in knots"),
        ("advance_ratio_implied", round(s["advance_ratio"], 4), "-", "V / (Omega R)"),
        ("mach_at_station_lowest", round(s["mach_min"], 4), "-",
         "(Omega r - V) / a: the value Case B holds steady"),
        ("mach_at_station_highest", round(s["mach_max"], 4), "-",
         "(Omega r + V) / a, on the advancing side: not represented in Case B"),
        ("pitch_frequency", round(s["rev_per_s"], 6), "Hz", "once per revolution: rpm / 60"),
        ("reduced_frequency_k", round(s["reduced_frequency"], 6), "-", "Omega c / (2 U)"),
        ("incidence_history", "%g deg +/- %g deg, sinusoidal" % (pm.CASES[pm.CASE_B]["alpha_mean_deg"],
                                                                pm.CASES[pm.CASE_B]["alpha_amp_deg"]), "-",
         "assumed for this study; it comes from no trim calculation"),
    ]
    return pd.DataFrame(rows, columns=["quantity", "value", "units", "source"])


def air_properties() -> pd.DataFrame:
    """The properties of air the cases are run with (the first five build
    unistall.flowfield.Air; the rest follow from them)."""
    rows = [("gamma", AIR.gamma, "-", "ratio of specific heats"),
            ("R_gas", AIR.R_gas, "J/kg/K", "specific gas constant"),
            ("prandtl", AIR.prandtl, "-", "Prandtl number of air"),
            ("T_inf", AIR.T_inf, "K", "sea-level standard temperature"),
            ("p_inf", AIR.p_inf, "Pa", "sea-level standard pressure"),
            ("a_sound", repr(AIR.a_sound), "m/s", "derived: sqrt(gamma R T)"),
            ("rho", round(AIR.rho, 5), "kg/m^3", "derived: p / (R T)"),
            ("mu", float("%.5g" % AIR.mu), "Pa.s", "derived: Sutherland's law"),
            ("cp", round(AIR.cp, 3), "J/kg/K", "derived: gamma R / (gamma - 1)"),
            ("recovery_factor", round(AIR.recovery, 4), "-", "derived: " + flowfield.RECOVERY_SOURCE)]
    return pd.DataFrame(rows, columns=["property", "value", "units", "note"])


def static_inputs() -> tuple:
    """(table, scalars): the static curves the solver uses at each case's Mach
    number, incidence in degrees, and the scalar static inputs per case."""
    alpha = np.arange(0.0, ALPHA_MAX + STATIC_STEP_DEG/2.0, STATIC_STEP_DEG)
    frames, scalars = [], {}
    for name, c in pm.CASES.items():
        sm = StaticModel(c["M"])
        frames.append(pd.DataFrame({"case": name, "mach_M": c["M"], "alpha_deg": alpha,
                                    "CN_static": np.round(sm.CN_static(alpha), 5),
                                    "CM_static_c4": np.round(sm.cm_static(alpha), 5),
                                    "f_separation": np.round(sm.f(alpha), 5)}))
        scalars[name] = dict(CN_alpha_per_rad=round(float(sm.CN_alpha), 5),
                             alpha0_deg=round(float(sm.alpha0_deg), 4),
                             CN_static_at_zero_incidence=round(float(sm.CN_static(np.array([0.0]))[0]), 5),
                             alpha_stall_deg=round(float(sm.alpha_stall_deg), 4),
                             CN1_onset=round(float(sm.CN1), 5),
                             inside_static_mach_range=bool(sm.in_range))
    return pd.concat(frames, ignore_index=True), scalars


def solver_config(static_scalars: dict) -> dict:
    """Everything the solution stage is run with, as one record."""
    cases = {}
    for name, c in pm.CASES.items():
        steps, cycles = dsmodel.march_resolution(c["k"])
        cases[name] = dict(label=c["label"], description=description(name, c),
                           alpha_mean_deg=c["alpha_mean_deg"], alpha_amp_deg=c["alpha_amp_deg"],
                           reduced_frequency_k=c["k"], mach_M=c["M"], chord_m=c["chord_m"],
                           measured_frame=c["measured_frame"],
                           march=dict(steps_per_cycle=steps, cycles=cycles,
                                      step_semichords=round(2.0*np.pi/c["k"]/steps, 6)),
                           critical=dict(zip(("sonic_speed_ms", "Cp_critical"),
                                             (round(v, 4) for v in flowfield.sonic_values(c["M"], AIR)),
                                             strict=True), Cp_vacuum=round(flowfield.vacuum_cp(c["M"], AIR), 4)),
                           static_inputs=static_scalars[name])
    ff = flowfield
    return dict(
        author=pm.AUTHOR_FULL, title=pm.TITLE, solver_name=pm.SOLVER, solver=pm.SOLVER_LONG,
        method=pm.METHOD, package="unistall", load_model="unistall.dsmodel.solve", section=pm.SECTION,
        constants=dsmodel.load_constants(),
        constants_source="results/calibrated_constants.json, read by unistall.dsmodel.load_constants",
        march_rule=dict(largest_step_semichords=dsmodel.DS_MAX, settle_semichords=dsmodel.SETTLE_SEMICHORDS,
                        steps_min=dsmodel.STEPS_MIN, cycles_min=dsmodel.CYCLES_MIN,
                        cycles_max=dsmodel.CYCLES_MAX),
        compared_range=dict(k_max=dsmodel.K_MAX_COMPARED, peak_alpha_max_deg=dsmodel.PEAK_ALPHA_MAX_COMPARED),
        cases=cases,
        field_reconstruction=dict(
            module="unistall.flowfield",
            what_it_is=ff.STANDING_NOTE,
            quasi_steady_omissions=list(ff.QUASI_STEADY_OMISSIONS),
            read_by_the_driver=dict(
                note="04_solver/run_case.py passes these four to unistall.flowfield; changing one here changes "
                     "the output",
                n_panels=ff.N_PANELS, domain_chords=list(ff.DOMAIN_CHORDS), grid_nx=ff.GRID_NX, grid_ny=ff.GRID_NY),
            record_of_module_constants=dict(
                note="a record of the constants of unistall.flowfield at the time of writing; the driver does "
                     "not read them, and a check asserts that they still equal the module's",
                panels="straight panels, each with a source of constant strength and a vortex sheet of one "
                       "uniform strength, velocities from the closed-form panel integrals",
                near_wall_panel_lengths=ff.NEAR_WALL_PANEL_LENGTHS,
                near_wall_rings_blanked=ff.NEAR_WALL_RINGS_BLANKED,
                attached_outflow_fraction=ff.ATTACHED_OUTFLOW_FRACTION,
                vortex_x0=ff.VORTEX_X0, vortex_dx=ff.VORTEX_DX, vortex_y0=ff.VORTEX_Y0, vortex_dy=ff.VORTEX_DY,
                vortex_core_over_c=ff.VORTEX_CORE_OVER_C,
                lamb_oseen_peak_coefficient=round(ff.LAMB_OSEEN_PEAK, 6)),
            transpiration=dict(
                what="an equivalent transpiration through the suction surface standing for the separated "
                     "region: an outflow where positive, an inflow where negative; its strength is solved so "
                     "that the Kutta condition holds at the trailing edge while the circulation stays that of "
                     "the predicted lift",
                wall_condition="v_n = v_t w on the suction-surface panels aft of x/c = min(f_sep, 1 - "
                               "attached_outflow_fraction); w a smooth step from 0 there to 1 at the trailing edge",
                kutta_condition="equal pressure at the two control points next to the trailing edge"),
            compressibility=dict(
                what="a linearised correction applied to an incompressible reconstruction, valid only where the "
                     "local flow is subcritical; not a compressible solution",
                pressure="Karman-Tsien correction of the incompressible Cp at the stream Mach number: "
                         "Cp_c = Cp / (beta + M^2 Cp / (2 (1 + beta))), beta = sqrt(1 - M^2)",
                vacuum_cap=ff.VACUUM_CAP,
                vacuum_cap_note="the corrected Cp is not allowed below this fraction of the vacuum value "
                                "-2/(gamma M^2) (cases.<case>.critical.Cp_vacuum)",
                beyond_critical="where the corrected Cp is below cases.<case>.critical.Cp_critical the local Mach "
                                "number is above one, the correction does not hold and the values are not "
                                "physical; those points are flagged, hatched in the figures and counted",
                local_state="local Mach number, static temperature and recovery temperature from the corrected "
                            "Cp by the isentropic relations with the stagnation state of the stream",
                recovery_factor=round(AIR.recovery, 4), recovery_factor_source=ff.RECOVERY_SOURCE,
                circulation="the sheet carries the fraction of the circulation of the lift at which the "
                            "corrected surface pressure returns the lift given; the fraction is solved at each "
                            "instant, starting from the Prandtl-Glauert factor sqrt(1 - M^2)"),
            field_files=dict(format="gzip-compressed CSV, one row per grid node",
                             columns=["x_m", "y_m", "u_ms", "v_ms", "Cp", "Cp_incompressible", "Mach", "T_static_K",
                                      "T_recovery_K", "vorticity_1s", "beyond_critical", "inside"],
                             velocity="u_ms, v_ms: the incompressible velocity",
                             beyond_critical="1 where the corrected Cp is below the critical value: Cp, Mach and "
                                             "the temperatures are not physical there; 0 elsewhere",
                             inside="1 at a node inside the section, where the values are those of the "
                                    "singularities and have no physical meaning; 0 elsewhere"),
            bound_circulation="Gamma = 0.5 CL U c, uniform vortex sheet on the surface",
            vortex=dict(what="a marker of the load model's vortex clock; it does not reproduce a suction footprint",
                        circulation="0.5 CN_vortex U c", profile="Lamb-Oseen",
                        position="x/c = x0 + dx tau, y/c = y0 + dy tau, tau = tau_v / Tvl; over the trailing "
                                 "edge at tau = 1",
                        path_and_core="assumed: chosen for illustration, neither derived nor measured"),
            field_phases=FIELD_PHASES))


def case_definition(flow: pd.DataFrame, kin: pd.DataFrame, station: pd.DataFrame, scalars: dict,
                    cfg: dict) -> str:
    """The text of 00_overview/case_definition.md. Every number in it is read
    from the tables handed in."""
    f, k, st = flow.set_index("parameter"), kin.set_index("case"), station.set_index("quantity")["value"]
    a, b = pm.CASE_A, pm.CASE_B
    rng, fr = cfg["compared_range"], cfg["field_reconstruction"]

    def row(label: str, fmt: str, get) -> str:
        return f"| {label} | {fmt.format(get(a))} | {fmt.format(get(b))} |"

    frame = pm.CASES[a]["measured_frame"].split("_")[1]
    table = "\n".join([
        f"| | {pm.CASES[a]['label']} | {pm.CASES[b]['label']} |", "|---|---|---|",
        f"| Purpose | a condition of the NASA TM-84245 experiment, so the result can be set beside a measured loop "
        f"(frame {frame}, a calibration loop) | an illustration at the condition of a blade section at "
        f"{st['station_r_over_R']} of the radius |",
        row("Mach number", "{:.2f}", lambda n: float(f.loc["freestream_mach_M", n])),
        row("Reduced frequency k", "{:.4f}", lambda n: float(k.loc[n, "reduced_freq_k"])),
        row("Pitch frequency", "{:.3f} Hz", lambda n: float(k.loc[n, "freq_Hz"])),
        f"| Incidence | {k.loc[a, 'alpha_mean_deg']:g}° ± {k.loc[a, 'alpha_amp_deg']:g}° | "
        f"{k.loc[b, 'alpha_mean_deg']:g}° ± {k.loc[b, 'alpha_amp_deg']:g}° |",
        row("Chord", "{:g} m", lambda n: float(f.loc["chord_c", n])),
        row("Chord Reynolds number in sea-level standard air (not the tunnel value; not used by the model)",
            "{:.3g}", lambda n: float(f.loc["reynolds_number_Re_c", n]))])
    offsets = "; ".join(
        f"{pm.short_label(n)}: zero-lift incidence {scalars[n]['alpha0_deg']:+.2f}°, static normal force at zero "
        f"incidence {scalars[n]['CN_static_at_zero_incidence']:+.3f}" for n in (a, b))
    return f"""# Case definition

**Author:** {pm.AUTHOR_FULL}

This file is written by `03_model_setup/generate_setup.py`; every number in it
is read from the tables of `03_model_setup`.

## Title
Prediction of dynamic stall on a pitching NACA {pm.SECTION} section with UNISTALL, a
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
The unsteady loads on a NACA {pm.SECTION} section in a steady stream with prescribed
sinusoidal pitch, α(t) = α_mean + α_amp sin(ωt), by the load model in
`unistall/dsmodel.py`. Two conditions:

{table}

Both lie inside the range the model has been compared with measurement
(reduced frequency up to {rng['k_max']:g}, peak incidence up to {rng['peak_alpha_max_deg']:g}°, at the Mach
numbers of the static data).

**Case A is set beside a calibration loop.** The measured loop drawn with it
(frame {frame}) was used to fit the model, so the comparison is an
illustration and not an independent test.

**Case B is not a rotor calculation.** Its chord and its once-per-revolution
pitch frequency are those of the UH-60A main rotor ({st['number_of_blades']} blades, radius
{st['rotor_radius']} ft, chord {st['blade_chord']} in, {st['rotor_speed']} rpm; {pm.ROTOR_SOURCE.split(' (')[0]}).
The station at {st['station_r_over_R']} of the radius turns at {float(st['station_rotational_speed']):.1f} m/s. Case B
takes the Mach number {float(st['station_mach_M']):.2f} for it, which is what the station sees at
azimuth {float(st['station_azimuth']):.0f}° (the retreating side) at a forward speed of
{float(st['forward_speed_implied']):.1f} m/s ({float(st['forward_speed_implied_kt']):.0f} kt, advance ratio {float(st['advance_ratio_implied']):.2f}); the reduced
frequency {float(st['reduced_frequency_k']):.4f} follows from the rotor speed, the chord and that speed
(`03_model_setup/station_condition.csv`). Round the azimuth the Mach number at
that station runs from {float(st['mach_at_station_lowest']):.2f} to {float(st['mach_at_station_highest']):.2f}; Case B holds the lowest value
steady, and its incidence history ({st['incidence_history'].replace(" deg +/- ", "° ± ").replace(" deg", "°")}) is assumed: it
comes from no trim calculation. The real blade uses the SC1095 and SC1094 R8
sections with a swept tip, not a NACA {pm.SECTION}. There is no time-varying
velocity, sweep, inflow, trim or blade motion, and no rotor measurement is
compared with.

**The static inputs carry the offset of the measured static sweeps.** The
section is symmetric, yet the static inputs read from the experiment give
{offsets}
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

**The flow fields are a reconstruction, not a flow solution.** {fr['what_it_is']}
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
{EXPERIMENT_CITATION} (measured loops and static data); Leishman & Beddoes
(1989), J. Am. Helicopter Soc. 34(3); Damiani & Hayman (2019),
NREL/TP-5000-66347 (formulation); {pm.ROTOR_SOURCE.split(' (')[0]} (rotor parameters of
Case B; not an entry of `docs/references.bib`, and the four parameters have
not been checked against its pages in this repository). The other sources are
listed in full in `docs/references.bib`.
"""


EXPERIMENT_CITATION = pm.experiment_citation()


def main() -> None:
    """Write the six tables and the case definition."""
    flow, kin, station = flow_conditions(), kinematics(), station_condition()
    flow.to_csv(HERE/"flow_conditions.csv", index=False)
    kin.to_csv(HERE/"kinematics.csv", index=False)
    station.to_csv(HERE/"station_condition.csv", index=False)
    air_properties().to_csv(HERE/"air_properties.csv", index=False)
    table, scalars = static_inputs()
    table.to_csv(HERE/"static_inputs.csv", index=False)
    cfg = solver_config(scalars)
    with open(HERE/"solver_config.json", "w", encoding="utf-8") as fh:
        json.dump(cfg, fh, indent=2, ensure_ascii=False)
    (OVERVIEW/"case_definition.md").write_text(case_definition(flow, kin, station, scalars, cfg), encoding="utf-8")
    for name, c in pm.CASES.items():
        U = c["M"]*AIR.a_sound
        print("[setup] %s: U = %.2f m/s, k = %.4f, Re_c = %.3g (standard air), march %s steps x %s cycles"
              % ((name, U, c["k"], AIR.rho*U*c["chord_m"]/AIR.mu) + dsmodel.march_resolution(c["k"])))


if __name__ == "__main__":
    main()
