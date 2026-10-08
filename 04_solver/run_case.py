# -*- coding: utf-8 -*-
# Run from the repository root:  PYTHONPATH=. python3 04_solver/run_case.py
"""
04_solver / run_case.py
-----------------------
Solves the two cases with the load model (unistall.dsmodel.solve) and writes the
solution data to 05_solution/. The flow fields and surface pressures written
here are the potential-flow reconstruction of unistall.flowfield, drawn around
the predicted lift: a reconstruction, not a flow solution, incompressible with
a linearised (Karman-Tsien) correction of the pressure that holds only where
the local flow is subcritical, quasi-steady, compared with no measurement. It is given the lift only; the
moment and the pressure drag of its surface pressure are written beside the
model's and do not agree with them. No load depends on it.

Author: Akosa Samuel Onyejekwe (independent)

Reads 03_model_setup/ (through project_meta.read_setup),
01_geometry/naca0012_coordinates.csv and results/validation_table.csv. Writes:
  time_history_<case>.csv           loads and states at every step of the last cycle
  cp_distribution_<case>.csv        surface Cp (corrected and incompressible),
                                    local Mach number, static and recovery
                                    temperature, speed, transpiration and the
                                    flag of the critical pressure at six instants
  field_<case>_<phase>_a<deg>.csv.gz   reconstructed field at four instants
                                    (gzip-compressed CSV): rise, peak (greatest
                                    lift), dsv (greatest vortex normal force),
                                    fall; <deg> is the incidence in whole
                                    degrees. Columns: x_m, y_m, u_ms, v_ms (the
                                    incompressible velocity), Cp (corrected),
                                    Cp_incompressible, Mach, T_static_K,
                                    T_recovery_K, vorticity_1s, beyond_critical
                                    (1 where the corrected Cp is below the
                                    critical value: not physical there) and
                                    inside (1 at a node inside the section,
                                    where the values have no meaning).
  reconstruction_<case>.csv         one row per instant: the transpiration, the
                                    circulations, the vortex and the pressure
                                    beneath it, the extremes of Cp, where Cp is
                                    below the critical value, and the closure of
                                    lift, moment and drag
  cycle_scan_<case>.csv             the transpiration, the closure and the vortex
                                    at instants spread over the whole cycle
  metrics_<case>.csv                scalar results of the case
  summary_all_cases.csv             one row per case
  held_out_error_bands.csv          error of the model in each headline quantity
                                    on the held-out loops (mean and spread)
  static_station_spread.csv         Case B run with the static inputs of each
                                    Mach station in place of their blend
  model_static_polar_<case>.csv     lift, drag, moment and separation point of
                                    the model in a slow sweep at the case's
                                    Mach number: its quasi-steady limit, not a
                                    measurement
  convergence/residuals_<case>.csv  change of the last cycle when one more
                                    cycle is marched
  convergence/timestep_refinement.csv   loads against the number of steps
  convergence/timestep_order.csv    observed order, Richardson estimate and
                                    discretisation uncertainty of each load
  convergence/panel_convergence_<case>.csv   the reconstruction at the six
                                    instants against the number of panels
  convergence/panel_convergence_summary.csv  largest change on doubling the
                                    panels, against the stated tolerances
  convergence/reconstruction_timing.csv      CPU time of the reconstruction
  response_surface.csv              peak lift, minimum moment and damping over
                                    mean incidence and reduced frequency at
                                    Case B's Mach number and amplitude
  runtime_environment.csv           machine, library versions, CPU time
"""
import importlib.metadata as md
import platform
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd

import project_meta as pm
from unistall import dsmodel, flowfield
from unistall.metrics import cycle_damping
from unistall.static_model import StaticModel, stations

ROOT = Path(__file__).resolve().parent.parent
SOL = ROOT/"05_solution"
GEO = ROOT/"01_geometry"/"naca0012_coordinates.csv"
VTABLE = ROOT/"results"/"validation_table.csv"

CYCLE_SAMPLES = 361                                # instants of the cycle at which the reconstruction is scanned
REFINEMENT_FACTORS = (0.125, 0.25, 0.5, 1.0, 2.0, 4.0)   # multiples of the reported steps per cycle
ORDER_MIN, ORDER_MAX = 0.5, 3.0                    # an observed order outside this band is not used to extrapolate
PANEL_FACTORS = (0.25, 0.5, 1.0, 2.0)              # multiples of the configured number of panels
# Largest change tolerated in a quoted quantity when the panels are doubled from the configured number.
PANEL_TOLERANCES = {"Cp_surface_min_incompressible": ("percent of the value", 1.0),
                    "transpiration_vt_over_U": ("absolute", 0.005),
                    "closure_error_pct_incompressible": ("percentage points", 0.5),
                    "beyond_critical_surface_fraction": ("absolute", 0.005)}
POLAR_SWEEP = dict(mean_deg=7.5, amp_deg=12.5, k=0.002, cycles=2)   # slow sweep, -5 to 20 deg, for the polar
POLAR_STEP_DEG = 0.25                              # spacing of the incidences of the polar
POLAR_NOTE = ("the model's slow-sweep limit (upstroke of a sinusoidal sweep at reduced frequency %g), "
              "not a measurement" % POLAR_SWEEP["k"])
SURFACE_MEANS_DEG = np.arange(6.0, 16.0 + 1e-9, 2.0)     # response surface: mean incidence
SURFACE_KS = np.arange(0.04, 0.16 + 1e-9, 0.02)          # response surface: reduced frequency
SURFACE_CASE = pm.CASE_B                           # its Mach number and amplitude are used
LIBRARIES = ("numpy", "scipy", "pandas", "matplotlib")
FIELD_INSTANTS = ("rise", "peak", "dsv", "fall")
FIELD_FILE = re.compile(r"field_.+_(%s)_a\d+\.csv(\.gz)?" % "|".join(FIELD_INSTANTS))
GZIP = {"method": "gzip", "compresslevel": 9, "mtime": 0}     # mtime 0: the same bytes on every run
# Decimals kept in the field files: 1e-4 m, 0.01 m/s, 1e-4 in Cp, 0.1 1/s
FIELD_DECIMALS = {"x_m": 4, "y_m": 4, "u_ms": 2, "v_ms": 2, "Cp": 4, "Cp_incompressible": 4, "Mach": 4,
                  "T_static_K": 2, "T_recovery_K": 2, "vorticity_1s": 1, "beyond_critical": 0, "inside": 0}
FIELD_KEYS = {"x_m": "X", "y_m": "Y", "u_ms": "u", "v_ms": "v", "Cp": "Cp", "Cp_incompressible": "Cp_incompressible",
              "Mach": "Mach", "T_static_K": "T_static", "T_recovery_K": "T_recovery", "vorticity_1s": "vorticity",
              "beyond_critical": "beyond_critical", "inside": "inside"}
# Held-out error of the model in each headline quantity: (quantity, model column, measured column, kind)
BAND_QUANTITIES = (("CL_max", "CLmax_model", "CLmax_exp", "percent"),
                   ("CM_min_c4", "CMmin_model", "CMmin_exp", "absolute"),
                   ("CD_max", "CDmax_model", "CDmax_exp", "absolute"),
                   ("alpha_at_CLmax_deg", "alpha_CLmax_model", "alpha_CLmax_exp", "absolute"),
                   ("moment_stall_alpha_deg", "alpha_Mstall_model", "alpha_Mstall_exp", "absolute"))


# --------------------------------------------------------------------------- #
#  inputs
# --------------------------------------------------------------------------- #
def march(case: dict, air: flowfield.Air, steps: int | None = None, cycles: int | None = None,
          static: StaticModel | None = None, **cond: float) -> dict:
    """The load model for one case (last cycle). `steps` per cycle and
    `cycles` default to the case's own; `static` replaces the static inputs
    of the case's Mach number; `cond` overrides a_mean, a_amp [deg] or k. A
    condition outside the compared range raises ValueError."""
    c = {**case, **cond}
    return dsmodel.solve(c["a_mean"], c["a_amp"], c["k"], c["M"], static=static, n_per_cycle=steps or case["steps"],
                         n_cycles=cycles or case["cycles"], chord=c["chord"], a_sound=air.a_sound, strict=True)


class Instants:
    """The reconstruction of one case, one instant at a time. Each instant is
    solved once and kept; the panels themselves are built and inverted once
    by unistall.flowfield. An instant is solved for the Mach number of the
    case unless another is asked for: its sheet then carries the fraction of
    the circulation of the lift at which the pressure corrected for
    compressibility returns the lift."""

    def __init__(self, out: dict, case: dict, air: flowfield.Air, cfg: dict) -> None:
        self.out, self.case, self.air = out, case, air
        self.settings = cfg["field_reconstruction"]["read_by_the_driver"]
        self.Tvl = float(cfg["constants"]["Tvl"])
        self._kept = {}

    def at(self, i: int, vortex: bool = True, outflow: bool = True, n_panel: int | None = None,
           mach: float | None = None, factor: float | None = None) -> flowfield.Reconstruction:
        """The reconstruction at step i; the vortex or the transpiration may
        be switched off, the number of panels changed from the configured
        one, the Mach number changed from the case's (0: incompressible) and
        the fraction of the circulation fixed instead of solved."""
        mach = self.case["M"] if mach is None else mach
        key = (int(i), vortex, outflow, n_panel or self.settings["n_panels"], mach, factor)
        if key not in self._kept:
            self._kept[key] = flowfield.from_solution(GEO, self.case["chord"], self.case["U"], self.out, key[0],
                                                      self.Tvl, vortex=vortex, outflow=outflow, n_panel=key[3],
                                                      mach=mach, factor=factor, air=self.air)
        return self._kept[key]

    def field(self, i: int) -> dict:
        """The field of step i on the configured grid."""
        s = self.settings
        return flowfield.reconstruct_field(self.at(i), self.case["M"], self.air, tuple(s["domain_chords"]),
                                           s["grid_nx"], s["grid_ny"])


# --------------------------------------------------------------------------- #
#  the cycle
# --------------------------------------------------------------------------- #
def time_history(out: dict, Tvl: float) -> pd.DataFrame:
    """Loads and states at every step of the last cycle. Time [s], angles
    [deg], rate [rad/s], vortex clock in semichords."""
    t = out["t"]
    table = pd.DataFrame({
        "time_s": t, "phase_deg": 360.0*t/t[-1], "alpha_deg": out["alpha_deg"],
        "alpha_dot_rad_s": out["alpha_dot"], "alpha_f_deg": out["alpha_f_deg"],
        "CL": out["CL"], "CD": out["CD"], "CM_c4": out["CM"], "CN": out["CN"], "CC": out["CC"],
        "CN_prime": out["CN_prime"], "CN_vortex": out["CN_vortex"],
        "f_separation": out["f_sep"], "f_static": out["f_static"],
        "tau_v_semichords": out["tau_v"],
        "vortex_on_chord": ((out["tau_v"] > 0.0) & (out["tau_v"] <= Tvl)).astype(int)}).round(6)
    table["time_s"] = t.round(9)
    return table


def clock_runs(tau_v: np.ndarray) -> int:
    """How many times the vortex clock starts or restarts in the periodic
    cycle whose steps are `tau_v` (first and last step are the same instant):
    the number of ramps of the clock."""
    tau = np.asarray(tau_v, float)[:-1]
    before = np.roll(tau, 1)
    return int(np.sum((tau > 0.0) & ((before <= 0.0) | (tau < before))))


def stroke_index(out: dict, target_deg: float, upstroke: bool) -> int:
    """Step of the cycle nearest incidence target_deg [deg] on the named stroke."""
    idx = np.where((out["alpha_dot"] > 0.0) == upstroke)[0]
    return int(idx[np.argmin(np.abs(out["alpha_deg"][idx] - target_deg))])


def instants(out: dict, case: dict) -> dict:
    """Step index of each named instant. The four of the fields are rise,
    peak, dsv and fall; the surface pressures add the two mean-incidence
    crossings."""
    half = case["a_mean"] + 0.5*case["a_amp"]
    return {"mean_up": stroke_index(out, case["a_mean"], True),
            "rise": stroke_index(out, half, True),
            "peak": int(np.argmax(out["CL"])),
            "dsv": int(np.argmax(out["CN_vortex"])),
            "fall": stroke_index(out, half, False),
            "mean_down": stroke_index(out, case["a_mean"], False)}


def cp_distribution(R: Instants) -> pd.DataFrame:
    """Surface Cp (corrected for compressibility, and incompressible), local
    Mach number, static and recovery temperature [K], speed [m/s] of the
    incompressible solution, transpiration through the wall over U and the
    flag of the critical pressure at the panel control points, at each
    instant of `instants`."""
    out, frames = R.out, []
    for tag, i in instants(out, R.case).items():
        st = flowfield.surface_state(R.at(i), R.case["M"], R.air)
        frames.append(pd.DataFrame({
            "phase": tag, "alpha_deg": round(float(out["alpha_deg"][i]), 3),
            "stroke": "up" if out["alpha_dot"][i] > 0.0 else "down",
            "CL": round(float(out["CL"][i]), 5), "CN_vortex": round(float(out["CN_vortex"][i]), 5),
            "surface": np.where(st["upper"], "upper", "lower"),
            "x_over_c": st["x_over_c"].round(6), "s_over_c": st["s_over_c"].round(6),
            "Cp": st["Cp"].round(5), "Cp_incompressible": st["Cp_incompressible"].round(5),
            "Mach": st["Mach"].round(5), "T_static_K": st["T_static"].round(3),
            "T_recovery_K": st["T_recovery"].round(3), "speed_ms": st["speed"].round(3),
            "vn_over_U": st["vn_over_U"].round(5), "beyond_critical": st["beyond_critical"].astype(int)}))
    return pd.concat(frames, ignore_index=True)


def write_field(name: str, tag: str, alpha_deg: float, fld: dict) -> str:
    """Write one reconstructed field as gzip-compressed CSV; returns the file name."""
    table = pd.DataFrame({col: np.asarray(fld[key], float).ravel().round(FIELD_DECIMALS[col])
                          for col, key in FIELD_KEYS.items()})
    for flag in ("beyond_critical", "inside"):
        table[flag] = table[flag].astype(int)
    fname = "field_%s_%s_a%.0f.csv.gz" % (name, tag, alpha_deg)
    table.to_csv(SOL/fname, index=False, compression=GZIP)
    return fname


# --------------------------------------------------------------------------- #
#  the reconstruction as numbers
# --------------------------------------------------------------------------- #
def closure_row(R: Instants, i: int, n_panel: int | None = None) -> dict:
    """Lift, moment and drag returned by the surface pressure at step i.

    closure_error_pct: the stored pressure (corrected for compressibility)
    against the lift given, with the circulation solved so that it closes:
    what is left is the tolerance of that solution. circulation_factor is the
    fraction of the circulation of the lift the sheet then carries, beside
    the Prandtl-Glauert factor it starts from, and
    closure_error_pct_prandtl_glauert_factor is the closure with the
    circulation fixed at that factor. closure_error_pct_incompressible: the
    incompressible reconstruction (no correction, the whole circulation)
    against the lift given. closure_error_pct_panels_only: the same with a
    solid wall and no vortex: discretisation error. The moment and drag are
    those of the stored pressure, set beside the model's."""
    out, M, air = R.out, R.case["M"], R.air
    rec = R.at(i, n_panel=n_panel)
    clo = flowfield.surface_load_closure(rec, M, air)
    beta = float(np.sqrt(1.0 - M*M))
    fixed = flowfield.surface_load_closure(R.at(i, n_panel=n_panel, factor=beta), M, air)
    plain = flowfield.surface_load_closure(R.at(i, n_panel=n_panel, mach=0.0))
    bare = flowfield.surface_load_closure(R.at(i, vortex=False, outflow=False, n_panel=n_panel, mach=0.0))
    return dict(closure_error_pct=clo["error_pct"], closure_dCL=clo["dCL"], circulation_factor=rec.factor,
                prandtl_glauert_factor=beta, closure_error_pct_prandtl_glauert_factor=fixed["error_pct"],
                closure_error_pct_incompressible=plain["error_pct"],
                closure_error_pct_panels_only=bare["error_pct"], Cp_trailing_edge_jump=clo["te_jump"],
                CM_model=out["CM"][i], CM_from_Cp=clo["CM_from_Cp"],
                CM_from_Cp_minus_model=clo["CM_from_Cp"] - out["CM"][i],
                CD_model=out["CD"][i], CD_from_Cp=clo["CD_from_Cp"],
                CD_from_Cp_minus_model=clo["CD_from_Cp"] - out["CD"][i])


def state_row(st: dict, cover: dict, fld: dict | None) -> dict:
    """The corrected pressure and the local state as numbers: the extremes
    over the surface and, where a field is stored, over the grid, taken
    where the correction holds (outside the region beyond the critical
    pressure), and the size of that region."""
    ok = ~st["beyond_critical"]
    row = dict(Cp_surface_min=st["Cp"].min(), Cp_surface_max=st["Cp"].max(),
               Cp_surface_min_incompressible=st["Cp_incompressible"].min(),
               beyond_critical_surface_fraction=cover["surface_fraction"],
               capped_surface_fraction=cover["surface_fraction_capped"],
               Mach_surface_max_subcritical=st["Mach"][ok].max(),
               T_static_surface_min_subcritical_K=st["T_static"][ok].min(),
               T_recovery_surface_min_subcritical_K=st["T_recovery"][ok].min(),
               T_recovery_surface_max_K=st["T_recovery"].max())
    if fld is not None:
        good = ~fld["inside"] & ~fld["beyond_critical"]
        row.update(Cp_field_min_subcritical=fld["Cp"][good].min(), Cp_field_max=fld["Cp"][good].max(),
                   Mach_field_max_subcritical=fld["Mach"][good].max(),
                   T_static_field_min_subcritical_K=fld["T_static"][good].min(),
                   T_static_field_max_K=fld["T_static"][good].max(),
                   beyond_critical_field_area_over_c2=cover["field_area_over_c2"],
                   beyond_critical_field_fraction=cover["field_fraction"])
    return row


def instant_row(tag: str, R: Instants, i: int, fld: dict | None) -> dict:
    """One row of reconstruction_<case>.csv: the reconstruction at step i as
    numbers. `fld` is its field on the grid, or None if none is stored."""
    out, M = R.out, R.case["M"]
    rec = R.at(i)
    st = flowfield.surface_state(rec, M, R.air)
    row = dict(phase=tag, alpha_deg=out["alpha_deg"][i], stroke="up" if out["alpha_dot"][i] > 0.0 else "down",
               CL=out["CL"][i], f_separation=out["f_sep"][i])
    row.update({"transpiration_" + k: v for k, v in flowfield.outflow_state(rec).items()})
    row.update({"circulation_%s_over_Uc" % k: v for k, v in flowfield.circulation_state(rec).items()})
    vs, fp = flowfield.vortex_state(rec), flowfield.vortex_footprint(rec)
    row.update(vortex_clock_tau_over_Tvl=out["tau_v"][i]/R.Tvl, vortex_x_over_c=vs["x_over_c"],
               vortex_y_over_c=vs["y_over_c"], vortex_peak_vorticity_c_over_U=vs["peak_vorticity_c_over_U"],
               vortex_over_chord=int(fp["over_chord"]), Cp_under_vortex=fp["Cp_under"],
               Cp_under_vortex_without_vortex=fp["Cp_under_without_vortex"], dCp_under_vortex=fp["dCp_under"])
    row["Cp_critical"] = flowfield.sonic_values(M, R.air)[1]
    row.update(state_row(st, flowfield.critical_coverage(rec, M, R.air, fld), fld))
    row.update(closure_row(R, i))
    return row


def reconstruct_instants(name: str, R: Instants) -> tuple:
    """Write the four fields and return (table of reconstruction_<case>.csv, field file names)."""
    rows, written = [], []
    for tag, i in instants(R.out, R.case).items():
        fld = None
        if tag in FIELD_INSTANTS:
            fld = R.field(i)
            written.append(write_field(name, tag, float(R.out["alpha_deg"][i]), fld))
        rows.append(instant_row(tag, R, i, fld))
    return significant(pd.DataFrame(rows)), written


def significant(table: pd.DataFrame, digits: int = 5) -> pd.DataFrame:
    """The table with its numbers rounded to `digits` significant figures."""
    table = table.copy()
    num = table.select_dtypes("number").columns
    table[num] = table[num].apply(lambda c: c.map(lambda v: float("%.*g" % (digits, v))))
    return table


def cycle_scan(R: Instants) -> pd.DataFrame:
    """The reconstruction at CYCLE_SAMPLES instants spread evenly over the
    cycle: the transpiration, the fraction of the circulation at which the
    corrected pressure returns the lift and what is left of the closure, the
    part of the surface beyond the critical pressure, and the vortex."""
    out, rows = R.out, []
    for i in np.unique(np.linspace(0, len(out["CL"]) - 1, CYCLE_SAMPLES).astype(int)):
        rec = R.at(int(i))
        o, vs = flowfield.outflow_state(rec), flowfield.vortex_state(rec)
        clo = flowfield.surface_load_closure(rec, R.case["M"], R.air)
        cover = flowfield.critical_coverage(rec, R.case["M"], R.air)
        rows.append(dict(phase_deg=360.0*out["t"][i]/out["t"][-1], alpha_deg=out["alpha_deg"][i], CL=out["CL"][i],
                         CN_vortex=out["CN_vortex"][i], vortex_clock_tau_over_Tvl=out["tau_v"][i]/R.Tvl,
                         transpiration_vt_over_U=o["vt_over_U"], transpiration_flux_over_Uc=o["flux_over_Uc"],
                         circulation_factor=rec.factor, closure_dCL=clo["dCL"],
                         beyond_critical_surface_fraction=cover["surface_fraction"],
                         Cp_trailing_edge_jump=clo["te_jump"],
                         vortex_circulation_over_Uc=vs["Gamma_over_Uc"],
                         vortex_peak_vorticity_c_over_U=vs["peak_vorticity_c_over_U"],
                         vortex_x_over_c=vs["x_over_c"]))
    return significant(pd.DataFrame(rows), 6)


def panel_convergence(R: Instants) -> pd.DataFrame:
    """The quoted quantities of the reconstruction at the six instants,
    against the number of panels."""
    rows = []
    base = R.settings["n_panels"]
    for tag, i in instants(R.out, R.case).items():
        for fac in PANEL_FACTORS:
            n = int(round(base*fac))
            rec = R.at(i, n_panel=n)
            st, clo = flowfield.surface_state(rec, R.case["M"], R.air), closure_row(R, i, n)
            cover = flowfield.critical_coverage(rec, R.case["M"], R.air)
            aft = rec.xc >= 0.98*rec.chord
            steep = np.degrees(np.arctan2(np.abs(rec.yp[1:] - rec.yp[:-1]), np.abs(rec.xp[1:] - rec.xp[:-1])))
            rows.append(dict(phase=tag, n_panels=n, configured=n == base,
                             Cp_surface_min_incompressible=st["Cp_incompressible"].min(),
                             transpiration_vt_over_U=rec.vt/rec.U,
                             transpiration_flux_over_Uc=flowfield.outflow_state(rec)["flux_over_Uc"],
                             circulation_factor=rec.factor,
                             beyond_critical_surface_fraction=cover["surface_fraction"],
                             closure_error_pct=clo["closure_error_pct"],
                             closure_error_pct_incompressible=clo["closure_error_pct_incompressible"],
                             closure_error_pct_panels_only=clo["closure_error_pct_panels_only"],
                             Cp_trailing_edge_jump=clo["Cp_trailing_edge_jump"],
                             CM_from_Cp=clo["CM_from_Cp"], CD_from_Cp=clo["CD_from_Cp"],
                             steepest_panel_aft_of_0p98_deg=float(steep[aft].max())))
    return significant(pd.DataFrame(rows), 6)


def panel_summary(name: str, table: pd.DataFrame) -> list:
    """For each quoted quantity: the largest change over the six instants
    when the panels are doubled from half the configured number to it and
    from it to twice it, against the stated tolerance."""
    base = int(table[table["configured"]]["n_panels"].iloc[0])
    rows = []
    for q, (kind, tol) in PANEL_TOLERANCES.items():
        piv = table.pivot(index="phase", columns="n_panels", values=q)
        change = {}
        for lo, hi in ((base//2, base), (base, 2*base)):
            d = (piv[hi] - piv[lo]).abs()
            change[(lo, hi)] = float((100.0*d/piv[hi].abs()).max() if kind.startswith("percent of") else d.max())
        rows.append(dict(case=name, quantity=q, change_measured_as=kind, n_panels_configured=base,
                         largest_change_half_to_configured=float("%.4g" % change[(base//2, base)]),
                         largest_change_configured_to_double=float("%.4g" % change[(base, 2*base)]),
                         tolerance=tol, within_tolerance=change[(base, 2*base)] <= tol))
    return rows


# --------------------------------------------------------------------------- #
#  scalar results
# --------------------------------------------------------------------------- #
def load_rows(out: dict, case: dict) -> list:
    """Scalar results of the load model."""
    a, CL, CM, CD = out["alpha_deg"], out["CL"], out["CM"], out["CD"]
    iL, iM, iD = int(np.argmax(CL)), int(np.argmin(CM)), int(np.argmax(CD))
    return [("CL_max", CL[iL], "-"), ("alpha_at_CLmax_deg", a[iL], "deg"),
            ("CM_min_c4", CM[iM], "-"), ("alpha_at_CMmin_deg", a[iM], "deg"),
            ("CD_max", CD[iD], "-"), ("alpha_at_CDmax_deg", a[iD], "deg"), ("CD_min", CD.min(), "-"),
            ("stall_onset_alpha_deg", out["onset_alpha_deg"], "deg"),
            ("vortex_sheddings_per_cycle", out["n_sheddings"], "-"),
            ("vortex_clock_runs_per_cycle", clock_runs(out["tau_v"]), "-"),
            ("cycle_damping_Xi", cycle_damping(a, CM, case["a_amp"]), "-")]


def reconstruction_rows(R: Instants, table: pd.DataFrame, scan: pd.DataFrame) -> list:
    """Scalar results of the reconstruction: closure at greatest lift and over
    the cycle, the critical pressure, the transpiration over the cycle and
    the per-instant figures (from `table`, the rows of
    reconstruction_<case>.csv, and `scan`, those of cycle_scan_<case>.csv)."""
    t = table.set_index("phase")
    v_sonic, cp_star = flowfield.sonic_values(R.case["M"], R.air)
    vt = scan["transpiration_vt_over_U"]
    rows = [("Cp_closure_error_pct", t.loc["peak", "closure_error_pct"], "%"),
            ("Cp_closure_dCL", t.loc["peak", "closure_dCL"], "-"),
            ("circulation_factor_at_CLmax", t.loc["peak", "circulation_factor"], "-"),
            ("prandtl_glauert_factor", t.loc["peak", "prandtl_glauert_factor"], "-"),
            ("Cp_closure_error_pct_prandtl_glauert_factor",
             t.loc["peak", "closure_error_pct_prandtl_glauert_factor"], "%"),
            ("Cp_closure_error_pct_incompressible", t.loc["peak", "closure_error_pct_incompressible"], "%"),
            ("Cp_closure_error_pct_panels_only", t.loc["peak", "closure_error_pct_panels_only"], "%"),
            ("Cp_closure_worst_abs_pct_stored_instants", t["closure_error_pct"].abs().max(), "%"),
            ("Cp_closure_worst_abs_pct_incompressible_stored_instants",
             t["closure_error_pct_incompressible"].abs().max(), "%"),
            ("Cp_closure_worst_abs_dCL_over_cycle", scan["closure_dCL"].abs().max(), "-"),
            ("circulation_factor_cycle_min", scan["circulation_factor"].min(), "-"),
            ("circulation_factor_cycle_max", scan["circulation_factor"].max(), "-"),
            ("Cp_moment_minus_model_at_CLmax", t.loc["peak", "CM_from_Cp_minus_model"], "-"),
            ("Cp_drag_minus_model_at_CLmax", t.loc["peak", "CD_from_Cp_minus_model"], "-"),
            ("Cp_moment_minus_model_worst_abs_stored_instants", t["CM_from_Cp_minus_model"].abs().max(), "-"),
            ("Cp_drag_minus_model_worst_abs_stored_instants", t["CD_from_Cp_minus_model"].abs().max(), "-"),
            ("Cp_moment_sign_agrees_with_model_instants", int((np.sign(t["CM_from_Cp"]) == np.sign(t["CM_model"])).sum()), "-"),
            ("stored_instants", len(t), "-"),
            ("Cp_trailing_edge_jump_worst_stored_instants", t["Cp_trailing_edge_jump"].max(), "-"),
            ("Cp_critical", cp_star, "-"), ("Cp_vacuum", flowfield.vacuum_cp(R.case["M"], R.air), "-"),
            ("sonic_speed_over_U", v_sonic/R.case["U"], "-"),
            ("total_temperature_K", R.air.T_inf*(1.0 + 0.5*(R.air.gamma - 1.0)*R.case["M"]**2), "K"),
            ("recovery_factor", R.air.recovery, "-"),
            ("Cp_surface_min_stored_instants", t["Cp_surface_min"].min(), "-"),
            ("Cp_surface_min_incompressible_stored_instants", t["Cp_surface_min_incompressible"].min(), "-"),
            ("beyond_critical_surface_fraction_worst_stored_instants", t["beyond_critical_surface_fraction"].max(), "-"),
            ("beyond_critical_surface_fraction_worst_over_cycle", scan["beyond_critical_surface_fraction"].max(), "-"),
            ("beyond_critical_fraction_of_cycle", float((scan["beyond_critical_surface_fraction"] > 0.0).mean()), "-"),
            ("beyond_critical_field_area_over_c2_worst_stored_instants",
             t["beyond_critical_field_area_over_c2"].max(), "-"),
            ("beyond_critical_field_fraction_worst_stored_instants", t["beyond_critical_field_fraction"].max(), "-"),
            ("cycle_instants_scanned", len(scan), "-"),
            ("transpiration_vt_over_U_cycle_min", vt.min(), "-"), ("transpiration_vt_over_U_cycle_max", vt.max(), "-"),
            ("transpiration_inflow_fraction_of_cycle", float((vt < 0.0).mean()), "-"),
            ("transpiration_flux_over_Uc_cycle_min", scan["transpiration_flux_over_Uc"].min(), "-"),
            ("transpiration_flux_over_Uc_cycle_max", scan["transpiration_flux_over_Uc"].max(), "-")]
    per_instant = ("transpiration_vt_over_U", "transpiration_flux_over_Uc", "transpiration_x_sep",
                   "circulation_total_over_Uc", "circulation_factor", "Cp_surface_min", "Cp_surface_max",
                   "Cp_field_min_subcritical", "Cp_field_max", "Mach_surface_max_subcritical",
                   "Mach_field_max_subcritical", "T_static_surface_min_subcritical_K",
                   "T_static_field_min_subcritical_K", "T_recovery_surface_max_K",
                   "beyond_critical_surface_fraction", "beyond_critical_field_area_over_c2",
                   "beyond_critical_field_fraction", "closure_error_pct", "closure_error_pct_incompressible",
                   "dCp_under_vortex")
    for tag in FIELD_INSTANTS:
        rows += [(f"{col}_{tag}", t.loc[tag, col], "K" if col.endswith("_K") else "-") for col in per_instant]
    return rows


def vortex_rows(R: Instants) -> list:
    """The reconstructed vortex at the instant of its greatest normal force.
    Its position and core radius follow from constants chosen for
    illustration and are marked as assumed."""
    out = R.out
    iV = int(np.argmax(out["CN_vortex"]))
    vs = flowfield.vortex_state(R.at(iV))
    return [("CN_vortex_max", out["CN_vortex"][iV], "-"), ("alpha_at_CN_vortex_max_deg", out["alpha_deg"][iV], "deg"),
            ("vortex_clock_tau_over_Tvl", out["tau_v"][iV]/R.Tvl, "-"),
            ("vortex_circulation_over_Uc", vs["Gamma_over_Uc"], "-"),
            ("vortex_peak_swirl_over_U", vs["peak_swirl_over_U"], "-"),
            ("vortex_peak_vorticity_c_over_U", vs["peak_vorticity_c_over_U"], "-"),
            ("vortex_centre_Cp", vs["Cp_centre"], "-"),
            ("vortex_x_over_c", vs["x_over_c"], "- (assumed path)"),
            ("vortex_y_over_c", vs["y_over_c"], "- (assumed path)"),
            ("vortex_core_radius_over_c", vs["rc_over_c"], "- (assumed)")]


def metrics(R: Instants, table: pd.DataFrame, scan: pd.DataFrame) -> pd.DataFrame:
    """Scalar results of one case as (metric, value, units) rows."""
    out, case = R.out, R.case
    rows = load_rows(out, case) + reconstruction_rows(R, table, scan) + vortex_rows(R)
    rows += [("mean_alpha_deg", case["a_mean"], "deg"), ("amp_alpha_deg", case["a_amp"], "deg"),
             ("reduced_frequency_k", case["k"], "-"), ("mach_M", case["M"], "-"),
             ("chord_m", case["chord"], "m"), ("freestream_velocity_U", case["U"], "m/s"),
             ("reynolds_number_Re_c", case["Re"], "- (sea-level standard air)"),
             ("steps_per_cycle", out["meta"]["steps_per_cycle"], "-"), ("cycles_marched", out["meta"]["cycles"], "-")]
    return pd.DataFrame([(m, float("%.6g" % float(v)), u) for m, v, u in rows], columns=["metric", "value", "units"])


def summary_row(name: str, met: pd.DataFrame) -> dict:
    """The headline scalars of one case."""
    v = met.set_index("metric")["value"]
    keys = ("CL_max", "alpha_at_CLmax_deg", "CM_min_c4", "CD_max", "stall_onset_alpha_deg",
            "cycle_damping_Xi", "Cp_closure_error_pct", "Cp_closure_error_pct_incompressible",
            "Cp_closure_error_pct_panels_only", "circulation_factor_at_CLmax", "CN_vortex_max",
            "vortex_circulation_over_Uc")
    return dict(case=name, **{k: v[k] for k in keys})


def held_out_error_bands() -> pd.DataFrame:
    """The error of the model in each headline quantity on the held-out loops
    of results/validation_table.csv (tabulated model, inside the static Mach
    range): mean, standard deviation and largest size of model minus
    measured, as a percentage of the measured value or in the units of the
    quantity."""
    df = pd.read_csv(VTABLE)
    h = df[(df["model"] == "tabulated") & (df["set"] == "held_out") & df["in_mach_range"]]
    rows = []
    for quantity, model, measured, kind in BAND_QUANTITIES:
        pair = h[[model, measured]].dropna()
        err = pair[model] - pair[measured]
        if kind == "percent":
            err = 100.0*err/pair[measured].abs()
        rows.append(dict(quantity=quantity, error_is="model minus measured, " + (
            "percent of measured" if kind == "percent" else "in the units of the quantity"),
            n_loops=len(pair), mean_error=round(float(err.mean()), 4), std_error=round(float(err.std(ddof=1)), 4),
            mean_abs_error=round(float(err.abs().mean()), 4), largest_abs_error=round(float(err.abs().max()), 4),
            source="results/validation_table.csv: %s, %s" % (model, measured)))
    return pd.DataFrame(rows)


def static_station_spread(name: str, case: dict, air: flowfield.Air) -> pd.DataFrame:
    """The case run with the static inputs of each Mach station in place of
    their interpolation to the case's Mach number: what the blend between
    the stations is worth in the results."""
    def row(label: str, mach: float, o: dict) -> dict:
        return dict(case=name, static_inputs=label, static_mach=mach, CL_max=o["CL"].max(), CM_min=o["CM"].min(),
                    CD_max=o["CD"].max(), stall_onset_alpha_deg=o["onset_alpha_deg"])

    rows = [row("interpolated to the case's Mach number (reported)", case["M"], march(case, air))]
    for st in stations():
        rows.append(row("Mach station %.3f" % st.M, st.M, march(case, air, static=StaticModel(st.M))))
    df = pd.DataFrame(rows)
    spread = {q: df[q].iloc[1:].max() - df[q].iloc[1:].min()
              for q in ("CL_max", "CM_min", "CD_max", "stall_onset_alpha_deg")}
    df.loc[len(df)] = dict(case=name, static_inputs="spread between the stations", static_mach=float("nan"), **spread)
    return df.round(5)


def model_static_polar(case: dict, air: flowfield.Air) -> pd.DataFrame:
    """The static polar of the model at the Mach number of a case: lift,
    drag, moment about the quarter chord, normal and chord force and
    separation point on the
    upstroke of a sweep slow enough to be quasi-steady (POLAR_SWEEP), read at
    even steps of incidence. `outside_static_data` is 1 at an incidence
    outside the incidences the static data of the model are tabulated over.
    It is the model's slow-sweep limit, not a measurement."""
    w = POLAR_SWEEP
    o = dsmodel.solve(w["mean_deg"], w["amp_deg"], w["k"], case["M"], n_cycles=w["cycles"], chord=case["chord"],
                      a_sound=air.a_sound)
    up = o["alpha_dot"] > 0.0
    order = np.argsort(o["alpha_deg"][up])
    a = o["alpha_deg"][up][order]
    grid = np.arange(w["mean_deg"] - w["amp_deg"], w["mean_deg"] + w["amp_deg"] + 1e-9, POLAR_STEP_DEG)
    lo = max(float(st.alpha.min()) for st in stations())
    hi = min(float(st.alpha.max()) for st in stations())
    table = pd.DataFrame({"alpha_deg": grid})
    for col, key in (("CL", "CL"), ("CD", "CD"), ("CM_c4", "CM"), ("CN", "CN"), ("CC", "CC"), ("f_separation", "f_sep")):
        table[col] = np.interp(grid, a, o[key][up][order]).round(5)
    table["outside_static_data"] = ((grid < lo) | (grid > hi)).astype(int)
    table["mach_M"] = case["M"]
    table["note"] = POLAR_NOTE
    return table


# --------------------------------------------------------------------------- #
#  convergence
# --------------------------------------------------------------------------- #
def cycle_residuals(case: dict, air: flowfield.Air) -> pd.DataFrame:
    """The last cycle after marching n = 1 .. CYCLES_MAX cycles, and the
    largest change of each load over the cycle from n - 1 to n."""
    rows, before = [], None
    for n in range(1, dsmodel.CYCLES_MAX + 1):
        o = march(case, air, cycles=n)
        row = dict(cycles_marched=n, CL_max=o["CL"].max(), CM_min=o["CM"].min(), CD_max=o["CD"].max())
        for q in ("CL", "CM", "CD"):
            row["max_abs_change_" + q] = float(np.max(np.abs(o[q] - before[q]))) if before else float("nan")
        row["reported"] = n == case["cycles"]
        rows.append(row)
        before = o
    return pd.DataFrame(rows)


def timestep_refinement(cases: dict, air: flowfield.Air) -> pd.DataFrame:
    """Loads of each case against the number of steps per cycle, with the
    difference from the finest march."""
    frames = []
    for name, case in cases.items():
        rows = []
        for fac in REFINEMENT_FACTORS:
            n = int(round(case["steps"]*fac))
            o = march(case, air, steps=n)
            rows.append(dict(case=name, steps_per_cycle=n, step_semichords=2.0*np.pi/case["k"]/n,
                             CL_max=o["CL"].max(), CM_min=o["CM"].min(), CD_max=o["CD"].max(),
                             stall_onset_alpha_deg=o["onset_alpha_deg"],
                             cycle_damping_Xi=cycle_damping(o["alpha_deg"], o["CM"], case["a_amp"]),
                             reported=n == case["steps"]))
        df = pd.DataFrame(rows)
        for q in ("CL_max", "CM_min", "CD_max"):
            df["pct_from_finest_" + q] = 100.0*(df[q] - df[q].iloc[-1]).abs()/abs(df[q].iloc[-1])
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def timestep_order(refine: pd.DataFrame) -> pd.DataFrame:
    """Observed order of convergence, Richardson estimate and discretisation
    uncertainty of each load, from the three finest marches of each case
    (the reported number of steps, twice it and four times it).

    With f1 the finest, f2 the middle and f3 the reported value, the order is
    p = ln((f3 - f2)/(f2 - f1))/ln 2 where the two differences have the same
    sign, and the Richardson estimate is f1 + (f1 - f2)/(2^p - 1). An order
    is called observed only for p between ORDER_MIN and ORDER_MAX; then the
    uncertainty of the reported value is its distance from the estimate.
    Otherwise no order is observed and the uncertainty is the range of the
    three marches."""
    rows = []
    for name, sub in refine.groupby("case", sort=False):
        fine = sub.sort_values("steps_per_cycle").iloc[-3:]
        if not bool(fine["reported"].iloc[0]):
            raise ValueError("the three finest marches do not start at the reported one")
        for q in ("CL_max", "CM_min", "CD_max", "stall_onset_alpha_deg"):
            f3, f2, f1 = (float(v) for v in fine[q])
            d32, d21 = f3 - f2, f2 - f1
            monotone = d32*d21 > 0.0
            p = float(np.log(d32/d21)/np.log(2.0)) if monotone else float("nan")
            observed = bool(monotone and ORDER_MIN <= p <= ORDER_MAX)
            estimate = f1 + (f1 - f2)/(2.0**p - 1.0) if observed else float("nan")
            unc = abs(f3 - estimate) if observed else max(f1, f2, f3) - min(f1, f2, f3)
            rows.append(dict(case=name, quantity=q, steps_reported=int(fine["steps_per_cycle"].iloc[0]),
                             value_reported=f3, value_twice_the_steps=f2, value_four_times_the_steps=f1,
                             monotone=monotone, observed_order=p, order_observed=observed,
                             richardson_estimate=estimate, discretisation_uncertainty=unc,
                             uncertainty_pct_of_reported=100.0*unc/abs(f3),
                             uncertainty_basis="distance of the reported value from the Richardson estimate"
                             if observed else "range of the three finest marches (no order observed)"))
    return pd.DataFrame(rows)


def response_surface(case: dict, air: flowfield.Air) -> pd.DataFrame:
    """Peak lift, minimum moment and cycle damping over mean incidence and
    reduced frequency, at the Mach number and amplitude of `case`. Every
    point is inside the compared range (the march refuses any that is not)."""
    rows = []
    for am in SURFACE_MEANS_DEG:
        for k in SURFACE_KS:
            steps, cycles = dsmodel.march_resolution(float(k))
            o = march(case, air, steps=steps, cycles=cycles, a_mean=float(am), k=float(k))
            rows.append(dict(alpha_mean_deg=am, reduced_freq_k=round(float(k), 4), alpha_amp_deg=case["a_amp"],
                             mach_M=case["M"], peak_alpha_deg=am + case["a_amp"],
                             steps_per_cycle=steps, cycles=cycles,
                             CL_max=round(float(o["CL"].max()), 5), CM_min=round(float(o["CM"].min()), 5),
                             cycle_damping_Xi=round(cycle_damping(o["alpha_deg"], o["CM"], case["a_amp"]), 5),
                             stall_onset_alpha_deg=round(float(o["onset_alpha_deg"]), 3)))
    df = pd.DataFrame(rows)
    if df["peak_alpha_deg"].max() > dsmodel.PEAK_ALPHA_MAX_COMPARED:
        raise ValueError("response surface reaches beyond the compared range")
    return df


def runtime_environment(cases: dict, cpu: dict) -> pd.DataFrame:
    """The machine, the library versions and the CPU time of each march [s]."""
    rows = [("python", platform.python_version()), ("platform", platform.platform()),
            ("processor", platform.processor() or "unknown")]
    rows += [("version_" + lib, md.version(lib)) for lib in LIBRARIES]
    for name, case in cases.items():
        rows += [("steps_per_cycle_" + name, case["steps"]), ("cycles_" + name, case["cycles"]),
                 ("march_cpu_time_s_" + name, round(cpu[name], 3))]
    rows.append(("timing_note", "CPU time of the load-model march on this machine; it varies from run to run"))
    return pd.DataFrame(rows, columns=["property", "value"])


# --------------------------------------------------------------------------- #
def solve_case(name: str, case: dict, air: flowfield.Air, cfg: dict) -> tuple:
    """Solve one case and write its files; returns (summary row, CPU time of
    the march [s], row of the reconstruction timing, rows of the panel
    summary)."""
    t0 = time.process_time()
    out = march(case, air)
    cpu = time.process_time() - t0
    R = Instants(out, case, air, cfg)
    time_history(out, R.Tvl).to_csv(SOL/f"time_history_{name}.csv", index=False)
    cp_distribution(R).to_csv(SOL/f"cp_distribution_{name}.csv", index=False)
    t0 = time.process_time()
    table, fields = reconstruct_instants(name, R)
    scan = cycle_scan(R)
    recon_cpu = time.process_time() - t0
    table.to_csv(SOL/f"reconstruction_{name}.csv", index=False)
    scan.to_csv(SOL/f"cycle_scan_{name}.csv", index=False)
    met = metrics(R, table, scan)
    met.to_csv(SOL/f"metrics_{name}.csv", index=False)
    panels = panel_convergence(R)
    panels.to_csv(SOL/"convergence"/f"panel_convergence_{name}.csv", index=False)
    cycle_residuals(case, air).round(10).to_csv(SOL/"convergence"/f"residuals_{name}.csv", index=False)
    model_static_polar(case, air).to_csv(SOL/f"model_static_polar_{name}.csv", index=False)
    row = summary_row(name, met)
    timing = dict(case=name, what="the six instants, the four fields and the scan of the cycle",
                  instants=len(table), fields=len(fields), cycle_instants_scanned=len(scan),
                  n_panels=R.settings["n_panels"], cpu_time_s=round(recon_cpu, 2))
    print("[run] %s: CL_max %.3f at %.1f deg, CM_min %.3f, CD_max %.3f, onset %.2f deg, damping %.4f; "
          "surface-Cp closure at greatest lift %+.4f %% corrected (circulation factor %.3f), %+.2f %% "
          "incompressible, %+.2f %% panels only; %d fields; reconstruction %.1f s CPU"
          % (name, row["CL_max"], row["alpha_at_CLmax_deg"], row["CM_min_c4"], row["CD_max"],
             row["stall_onset_alpha_deg"], row["cycle_damping_Xi"], row["Cp_closure_error_pct"],
             row["circulation_factor_at_CLmax"], row["Cp_closure_error_pct_incompressible"],
             row["Cp_closure_error_pct_panels_only"], len(fields), recon_cpu))
    return row, cpu, timing, panel_summary(name, panels)


def main() -> None:
    """Solve both cases, then the convergence studies and the response surface."""
    (SOL/"convergence").mkdir(parents=True, exist_ok=True)
    for old in SOL.glob("field_*"):
        if FIELD_FILE.fullmatch(old.name):          # the fields of the four instants only
            old.unlink()
    cases, air, cfg = pm.read_setup()
    summary, cpu, timing, panel_rows = [], {}, [], []
    for name, case in cases.items():
        row, cpu[name], t_row, p_rows = solve_case(name, case, air, cfg)
        summary.append(row)
        timing.append(t_row)
        panel_rows += p_rows
    pd.DataFrame(summary).to_csv(SOL/"summary_all_cases.csv", index=False)
    pd.DataFrame(timing).to_csv(SOL/"convergence"/"reconstruction_timing.csv", index=False)
    pd.DataFrame(panel_rows).to_csv(SOL/"convergence"/"panel_convergence_summary.csv", index=False)
    held_out_error_bands().to_csv(SOL/"held_out_error_bands.csv", index=False)
    static_station_spread(pm.CASE_B, cases[pm.CASE_B], air).to_csv(SOL/"static_station_spread.csv", index=False)
    refine = timestep_refinement(cases, air)
    refine.round(8).to_csv(SOL/"convergence"/"timestep_refinement.csv", index=False)
    timestep_order(refine).round(8).to_csv(SOL/"convergence"/"timestep_order.csv", index=False)
    response_surface(cases[SURFACE_CASE], air).to_csv(SOL/"response_surface.csv", index=False)
    runtime_environment(cases, cpu).to_csv(SOL/"runtime_environment.csv", index=False)
    print("[run] done: solution written to 05_solution/")


if __name__ == "__main__":
    main()
