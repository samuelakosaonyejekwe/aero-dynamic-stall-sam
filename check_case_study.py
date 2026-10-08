# -*- coding: utf-8 -*-
"""
check_case_study.py
-------------------
Checks on the case study's outputs. Run from the repository root:

    PYTHONPATH=. python3 check_case_study.py

Author: Akosa Samuel Onyejekwe (independent)

Each check is a statement about the files in the numbered folders that must
hold if they were produced by the current scripts and are consistent with one
another. A failed check is printed and the script exits with status 1; no
check is skipped silently. Conditions are read from 03_model_setup, as every
stage reads them, and one check holds that folder against project_meta.
"""
import importlib.util
import json
import math
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import project_meta as pm
from unistall import dsmodel, flowfield, metrics, naca4

ROOT = Path(__file__).resolve().parent
SOL = ROOT/"05_solution"
CONV = SOL/"convergence"
SETUP = ROOT/"03_model_setup"
PLOTS = ROOT/"06_postprocessing"/"plots"
REPORT = ROOT/"07_report"
SHEETS = ROOT/"08_engineering_drawings"
GEO = ROOT/"01_geometry"/"naca0012_coordinates.csv"
CASE_FOLDERS = ("00_overview", "01_geometry", "02_mesh", "03_model_setup", "04_solver", "06_postprocessing",
                "06_postprocessing/validation", "07_report", "08_engineering_drawings")
VALIDATION = ROOT/"06_postprocessing"/"validation"
POLAR_STALL_TOL_PCT = 5.0       # slow-sweep lift of the model at static stall against C_N1 resolved to lift
SOLVER_ENTRY_TOL = 1e-6         # the single-condition entry against unistall.dsmodel.solve (it writes six decimals)
LOAD_TOL = 5e-6                 # the time histories are written to six decimals
OUTLINE_TOL = 1e-6              # chords: coordinates and mesh wall nodes against the analytic outline
PANEL_OUTLINE_TOL = 1e-5        # chords: panel end points (a spline through the coordinates) against it
CLOSURE_MAX_PCT = 5.0           # incompressible surface-pressure lift against the lift given, at every stored instant
TE_JUMP_MAX = 0.02              # difference of Cp across the trailing edge, at every stored instant
STEP_MAX_PCT = 0.5              # reported step against the finest march
SOLID_WALL_PCT = 1.0            # lowest Cp of a solid-wall solution: configured panels against twice as many
RECOMPUTE_TOL = 2e-4            # relative: a stored reconstruction row against the same row computed now
LETTERING_MIN_PT = 6.0          # smallest lettering of a figure as printed in the report
CORE_NODES_MIN = 4.0            # grid nodes across the core of the vortex marker
PHASES = ("rise", "peak", "dsv", "fall")
FIELD_QUANTITIES = ("Cp", "Mach", "Tstatic", "Trecovery", "speed_stream", "vorticity", "vectors")
CORRECTED_QUANTITIES = ("Cp", "Mach", "Tstatic", "Trecovery")     # maps that follow from the corrected pressure
FIELD_COLUMNS = ["x_m", "y_m", "u_ms", "v_ms", "Cp", "Cp_incompressible", "Mach", "T_static_K", "T_recovery_K",
                 "vorticity_1s", "beyond_critical", "inside"]
SURFACE_FIGURES = ("fig3d_response_surface.png", "fig3d_cp_phase_surface.png", "fig3d_field_surface_Cp.png",
                   "fig3d_field_surface_speed.png", "fig3d_section_vectors.png")
CORRECTED_CLOSURE_MAX_PCT = 0.01   # lift of the stored (corrected) pressure against the lift given
ISENTROPIC_TOL = 5e-3              # local Mach number against the stored corrected Cp, as files round them
TEMPERATURE_TOL_K = 0.05           # stored temperatures against the isentropic relations
DOCUMENTS = ("case.docx", "UNISTALL_report.pdf", "UNISTALL_plots_album.pdf", "UNISTALL_data_dossier.pdf")


def _metric(name: str) -> pd.Series:
    return pd.read_csv(SOL/f"metrics_{name}.csv").set_index("metric")["value"]


def _script(rel: str) -> object:
    """A script of a numbered folder, loaded as a module."""
    spec = importlib.util.spec_from_file_location(Path(rel).stem, ROOT/rel)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _sources() -> dict:
    """The text of every Python script of the case study, by relative path."""
    files = [p for folder in CASE_FOLDERS for p in sorted((ROOT/folder).glob("*.py"))]
    files += [ROOT/"run_all.py", ROOT/"project_meta.py", ROOT/"assets"/"make_favicon.py"]
    return {str(p.relative_to(ROOT)): p.read_text(encoding="utf-8") for p in files}


def _peak_reconstruction(name: str, case: dict, cfg: dict) -> tuple:
    """(reconstruction, stored row of the time history) at greatest lift, from the stored history."""
    th = pd.read_csv(SOL/f"time_history_{name}.csv")
    row = th.loc[th["CL"].idxmax()]
    n = cfg["field_reconstruction"]["read_by_the_driver"]["n_panels"]
    return flowfield.solve(GEO, case["chord"], case["U"], float(row.alpha_deg), float(row.CL), float(row.CN_vortex),
                           float(row.tau_v_semichords)/case["Tvl"], float(row.f_separation), n_panel=n,
                           mach=case["M"]), row


def setup_agrees_with_the_definition(cases: dict, air: flowfield.Air, cfg: dict) -> list:
    """03_model_setup holds the cases of project_meta, and nothing else reads project_meta for numbers."""
    out = []
    for name, c in pm.CASES.items():
        s = cases[name]
        same = (s["a_mean"] == c["alpha_mean_deg"] and s["a_amp"] == c["alpha_amp_deg"] and s["M"] == c["M"]
                and s["chord"] == c["chord_m"] and math.isclose(s["k"], c["k"], rel_tol=1e-12)
                and math.isclose(s["U"], c["M"]*air.a_sound, rel_tol=1e-12))
        out.append((f"{name}: the conditions in 03_model_setup are those of project_meta, the stream speed "
                    f"without rounding", same))
    number = r'\[\s*"(k|M|chord_m|alpha_mean_deg|alpha_amp_deg)"\s*\]'
    readers = [rel for rel, text in _sources().items()
               if rel not in ("03_model_setup/generate_setup.py", "project_meta.py")
               and (re.search(r"CASES\[[^\]]+\]" + number, text)
                    or (re.search(r"CASES\.items\(\)", text) and re.search(r"\b(c|case)" + number, text)))]
    out.append(("no script but the setup stage takes a condition from project_meta.CASES"
                + (f" (found: {readers})" if readers else ""), not readers))
    text = (ROOT/"00_overview"/"case_definition.md").read_text(encoding="utf-8")
    gs = _script("03_model_setup/generate_setup.py")
    table, scalars = gs.static_inputs()
    fresh = gs.case_definition(gs.flow_conditions(), gs.kinematics(), gs.station_condition(), scalars,
                               gs.solver_config(scalars))
    out.append(("00_overview/case_definition.md is what the setup stage writes from its tables", text == fresh))
    return out


def case_b_follows_from_the_rotor(cases: dict, air: flowfield.Air) -> list:
    """Mach number and reduced frequency of Case B are derived, and the pitch frequency is the rotor speed."""
    kin = pd.read_csv(SETUP/"kinematics.csv").set_index("case")
    st = pd.read_csv(SETUP/"station_condition.csv").set_index("quantity")["value"]
    rev = pm.ROTOR["speed_rpm"]/60.0
    freq = float(kin.loc[pm.CASE_B, "freq_Hz"])
    b = cases[pm.CASE_B]
    omega = 2.0*math.pi*rev
    out = [(f"Case B: pitch frequency {freq:.4f} Hz equals rotor speed / 60 = {rev:.4f} Hz to four figures",
            float("%.4g" % freq) == float("%.4g" % rev)),
           (f"Case B: reduced frequency {b['k']:.5f} equals Omega c / (2 U)",
            math.isclose(b["k"], omega*b["chord"]/(2.0*b["U"]), rel_tol=1e-9)),
           ("Case B: station_condition.csv gives the forward speed and the range of Mach number round the azimuth "
            f"({float(st['mach_at_station_lowest']):.2f} to {float(st['mach_at_station_highest']):.2f}), the lowest "
            "being the Mach number of the case",
            math.isclose(float(st["mach_at_station_lowest"]), b["M"], abs_tol=5e-5)
            and float(st["mach_at_station_highest"]) > b["M"] and float(st["forward_speed_implied"]) > 0.0)]
    dims = pd.read_csv(SHEETS/"dimensions.csv").set_index("name")["value"]
    same = (float(dims["Case B rotor radius, in feet"]) == pm.ROTOR["radius_ft"]
            and float(dims["Case B blade chord, in inches"]) == pm.ROTOR["chord_in"]
            and float(dims["Case B number of blades"]) == pm.ROTOR["blades"]
            and float(dims["Case B rotor speed"]) == pm.ROTOR["speed_rpm"])
    out.append(("drawings: the rotor parameters of dimensions.csv are those of project_meta.ROTOR", same))
    out.append(("drawings: the chords are the cases' chords",
                float(dims["Case A chord of the tunnel model"]) == cases[pm.CASE_A]["chord"]
                and float(dims["Case B blade chord"]) == b["chord"]))
    return out


def one_speed_of_sound(air: flowfield.Air) -> list:
    """The speed of sound has one definition and every march of the case study is given it."""
    prop = pd.read_csv(SETUP/"air_properties.csv").set_index("property")["value"]
    flow = pd.read_csv(SETUP/"flow_conditions.csv").set_index("parameter")
    a = math.sqrt(air.gamma*air.R_gas*air.T_inf)
    out = [(f"speed of sound {a:.3f} m/s is sqrt(gamma R T) in unistall.flowfield.Air, air_properties.csv and "
            f"flow_conditions.csv", air.a_sound == a and float(prop["a_sound"]) == a
            and all(float(flow.loc["speed_of_sound_a", n]) == a for n in pm.CASES))]
    literal, bare = [], []
    for rel, text in {**_sources(), "check_case_study.py": Path(__file__).read_text(encoding="utf-8")}.items():
        code = text.split("def one_speed_of_sound")[0] if rel == "check_case_study.py" else text
        if re.search(r"\b340(\.\d+)?\b", code):
            literal.append(rel)
        for call in re.finditer(r"dsmodel\.solve\((.*?)\n\n", code, re.S):
            head = call.group(1).split(")\n")[0]
            if "a_sound" not in head and "n_cycles=2" not in head:
                bare.append(rel)
    out.append(("no script of the case study types a speed of sound" + (f" (found in {literal})" if literal else ""),
                not literal))
    out.append(("every march whose loads are kept is given the speed of sound of 03_model_setup"
                + (f" (not in {bare})" if bare else ""), not bare))
    return out


def one_outline(cfg: dict) -> list:
    """Geometry file, mesh wall and panels lie on the one closed outline."""
    t = naca4.thickness_ratio(pm.SECTION)
    co = pd.read_csv(GEO)
    summ = pd.read_csv(ROOT/"01_geometry"/"section_geometry_summary.csv").set_index("property")["value"]
    err = float(np.max(np.abs(np.abs(co.y_over_c) - naca4.closed_thickness(co.x_over_c.values, t))))
    out = [(f"geometry: the coordinate file lies on the closed outline of unistall.naca4 (largest difference "
            f"{err:.1e} c) and its trailing-edge gap is {float(summ['TE_gap_percent_chord']):g} %",
            err <= OUTLINE_TOL and float(summ["TE_gap_percent_chord"]) == 0.0
            and co.y_over_c.iloc[0] == 0.0 and co.y_over_c.iloc[-1] == 0.0)]
    nodes = pd.read_csv(ROOT/"02_mesh"/"mesh_nodes.csv")
    wall = nodes[nodes.j == 0]
    err = float(np.max(np.abs(np.abs(wall.y_over_c) - naca4.closed_thickness(wall.x_over_c.values, t))))
    out.append((f"mesh: the wall nodes lie on the same outline (largest difference {err:.1e} c)", err <= OUTLINE_TOL))
    xp, yp = flowfield.panel_geometry(GEO, 1.0, cfg["field_reconstruction"]["read_by_the_driver"]["n_panels"])
    X, Y, _, _ = naca4.outline(pm.SECTION, 20001)
    ax, ay, bx, by = X[:-1], Y[:-1], np.diff(X), np.diff(Y)

    def distance(x: float, y: float) -> float:
        """Distance from a point to the outline, taken as the fine polygon (X, Y)."""
        s = np.clip(((x - ax)*bx + (y - ay)*by)/(bx*bx + by*by), 0.0, 1.0)
        return float(np.min(np.hypot(x - ax - s*bx, y - ay - s*by)))

    far = max(distance(x, y) for x, y in zip(xp, yp, strict=True))
    out.append((f"reconstruction: the panel end points lie on the same outline (largest distance {far:.1e} c)",
                far <= PANEL_OUTLINE_TOL))
    steep = max(float(pd.read_csv(CONV/f"panel_convergence_{n}.csv")["steepest_panel_aft_of_0p98_deg"].max())
                for n in pm.CASES)
    out.append((f"reconstruction: no panel aft of x/c = 0.98 is within 45 deg of the y axis (steepest "
                f"{steep:.1f} deg from the chord line)", steep < 45.0))
    typed = [rel for rel, text in _sources().items() if re.search(r"(?<![\d.])0\.12(?!\d)", text)]
    out.append(("the thickness ratio is typed in no script of the case study: it follows from the designation "
                f"NACA {pm.SECTION}" + (f" (found in {typed})" if typed else ""), not typed))
    return out


def loads_come_from_the_solver(cases: dict, air: flowfield.Air) -> list:
    """The stored time histories are what the solver returns now."""
    out = []
    for name, c in cases.items():
        th = pd.read_csv(SOL/f"time_history_{name}.csv")
        o = dsmodel.solve(c["a_mean"], c["a_amp"], c["k"], c["M"], n_per_cycle=c["steps"], n_cycles=c["cycles"],
                          chord=c["chord"], a_sound=air.a_sound, strict=True)
        worst = max(float(np.max(np.abs(th[col].values - o[key]))) for col, key in
                    (("CL", "CL"), ("CD", "CD"), ("CM_c4", "CM"), ("CN", "CN"), ("CN_vortex", "CN_vortex")))
        out.append((f"{name}: stored loads equal the solver's (largest difference {worst:.1e})", worst <= LOAD_TOL))
        out.append((f"{name}: inside the range the model was compared with measurement",
                    o["meta"]["outside_compared_range"] == []))
    return out


def metrics_follow_from_the_histories(cases: dict) -> list:
    """Peak lift, minimum moment, damping and the count of vortex-clock runs are those of the stored histories."""
    out = []
    run_case = _script("04_solver/run_case.py")
    for name, c in cases.items():
        th, m = pd.read_csv(SOL/f"time_history_{name}.csv"), _metric(name)
        xi = metrics.cycle_damping(th.alpha_deg.values, th.CM_c4.values, c["a_amp"])
        ok = (abs(th.CL.max() - float(m["CL_max"])) < 1e-4 and abs(th.CM_c4.min() - float(m["CM_min_c4"])) < 1e-4
              and abs(xi - float(m["cycle_damping_Xi"])) < 1e-4)
        out.append((f"{name}: metrics agree with the time history", ok))
        summ = pd.read_csv(SOL/"summary_all_cases.csv").set_index("case").loc[name]
        out.append((f"{name}: summary row agrees with the metrics file",
                    abs(float(summ.CL_max) - float(m["CL_max"])) < 1e-4))
        runs = run_case.clock_runs(th.tau_v_semichords.values)
        out.append((f"{name}: the vortex clock runs {runs} times in the stored cycle, which is the count in the "
                    f"metrics ({int(m['vortex_clock_runs_per_cycle'])}) and the sheddings the load model counts "
                    f"({int(m['vortex_sheddings_per_cycle'])})",
                    runs == int(m["vortex_clock_runs_per_cycle"]) == int(m["vortex_sheddings_per_cycle"])))
        out.append((f"{name}: lowest drag {float(m['CD_min']):+.4f} is not negative", float(m["CD_min"]) >= 0.0))
    return out


def reconstruction_is_current_and_consistent(cases: dict, air: flowfield.Air, cfg: dict) -> list:
    """The stored reconstruction is what the module computes now; at every stored instant its surface
    pressure returns the lift given within the limit and is continuous at the trailing edge; its moment
    and drag are stored beside the model's."""
    out = []
    for name, c in cases.items():
        rec = pd.read_csv(SOL/f"reconstruction_{name}.csv")
        stored = rec.set_index("phase").loc["peak"]
        now, _ = _peak_reconstruction(name, c, cfg)
        clo, st = flowfield.surface_load_closure(now, c["M"], air), flowfield.surface_state(now, c["M"], air)
        fresh = dict(transpiration_vt_over_U=now.vt/now.U, Cp_surface_min_incompressible=st["Cp_incompressible"].min(),
                     circulation_factor=now.factor, CM_from_Cp=clo["CM_from_Cp"], CD_from_Cp=clo["CD_from_Cp"],
                     beyond_critical_surface_fraction=flowfield.critical_coverage(now, c["M"], air)["surface_fraction"])
        worst = max(abs(v - float(stored[k]))/max(abs(v), 1e-3) for k, v in fresh.items())
        out.append((f"{name}: the stored reconstruction at greatest lift is what unistall.flowfield computes now "
                    f"(largest relative difference {worst:.1e})", worst <= RECOMPUTE_TOL))
        for r in rec.itertuples():
            out.append((f"{name}, {r.phase}: lift of the stored (corrected) surface pressure {r.closure_error_pct:+.1e} % "
                        f"from the lift given (limit {CORRECTED_CLOSURE_MAX_PCT} %), the sheet carrying "
                        f"{r.circulation_factor:.3f} of the circulation; incompressible pressure with the whole "
                        f"circulation {r.closure_error_pct_incompressible:+.2f} % (limit {CLOSURE_MAX_PCT} %); Cp "
                        f"across the trailing edge differs by {r.Cp_trailing_edge_jump:.1e} (limit {TE_JUMP_MAX})",
                        abs(r.closure_error_pct) <= CORRECTED_CLOSURE_MAX_PCT
                        and abs(r.closure_error_pct_incompressible) <= CLOSURE_MAX_PCT
                        and 0.5 < r.circulation_factor < 1.05
                        and r.Cp_trailing_edge_jump < TE_JUMP_MAX))
        cols = {"CM_model", "CM_from_Cp", "CM_from_Cp_minus_model", "CD_model", "CD_from_Cp", "CD_from_Cp_minus_model",
                "circulation_bound_over_Uc", "circulation_total_over_Uc", "closure_error_pct",
                "closure_error_pct_incompressible", "closure_error_pct_prandtl_glauert_factor", "circulation_factor",
                "beyond_critical_surface_fraction"}
        out.append((f"{name}: moment and drag of the surface pressure are stored beside the model's, and the total "
                    f"circulation beside that of the lift, at all {len(rec)} instants (largest moment and drag "
                    f"differences {rec.CM_from_Cp_minus_model.abs().max():.3f} and "
                    f"{rec.CD_from_Cp_minus_model.abs().max():.3f}; reported, not a pass criterion)"
                    if cols <= set(rec.columns) else f"{name}: columns of the reconstruction are missing",
                    cols <= set(rec.columns) and bool(rec[sorted(cols)].notna().all().all())))
        pc = pd.read_csv(CONV/f"panel_convergence_{name}.csv")
        jump = pc.groupby("n_panels")["Cp_trailing_edge_jump"].max()
        out.append((f"{name}: the difference of Cp across the trailing edge stays below {TE_JUMP_MAX} at every "
                    f"panel count (largest {jump.max():.1e})", bool((jump < TE_JUMP_MAX).all())))
        over = rec[rec.vortex_over_chord == 1]
        out.append((f"{name}: the change of surface Cp beneath the vortex marker is stored at each of the "
                    f"{len(over)} instants at which it is over the chord (it raises the pressure at "
                    f"{int((over.dCp_under_vortex > 0).sum())} of them)",
                    bool(over.dCp_under_vortex.notna().all())))
    return out


def reconstruction_is_converged_in_panels(cases: dict, cfg: dict) -> list:
    """The quoted quantities change by less than their tolerance when the panels are doubled."""
    summ = pd.read_csv(CONV/"panel_convergence_summary.csv")
    out = [(f"{r.case}: {r.quantity} changes by {r.largest_change_configured_to_double:g} ({r.change_measured_as}) "
            f"when the {r.n_panels_configured} panels are doubled (tolerance {r.tolerance:g})",
            bool(r.within_tolerance) and r.largest_change_configured_to_double <= r.tolerance)
           for r in summ.itertuples()]
    n = cfg["field_reconstruction"]["read_by_the_driver"]["n_panels"]
    case = next(iter(cases.values()))
    low = [float(flowfield.surface_state(flowfield.solve(GEO, case["chord"], case["U"], 5.0, 0.6, n_panel=k), 0.0)
                 ["Cp"].min()) for k in (n, 2*n)]
    pct = 100.0*abs(low[0] - low[1])/abs(low[1])
    out.append((f"solid wall at 5 deg: the lowest Cp with {n} panels is within {pct:.2f} % of that with {2*n} "
                f"(limit {SOLID_WALL_PCT} %)", pct <= SOLID_WALL_PCT))
    return out


def vortex_marker_follows_the_clock(cases: dict, cfg: dict) -> list:
    """The marker is over the trailing edge when the clock says so, its core is resolved by the grid and its
    vorticity grows with the vortex normal force."""
    drv = cfg["field_reconstruction"]["read_by_the_driver"]
    x1 = flowfield.vortex_centre(1.0)[0]
    out = [(f"vortex marker: x/c = {x1:.2f} at tau_v / T_VL = 1", abs(x1 - 1.0) <= 0.01)]
    spacing = max((drv["domain_chords"][1] - drv["domain_chords"][0])/(drv["grid_nx"] - 1),
                  (drv["domain_chords"][3] - drv["domain_chords"][2])/(drv["grid_ny"] - 1))
    across = 2.0*flowfield.VORTEX_CORE_OVER_C/spacing
    out.append((f"vortex marker: its core is {across:.1f} grid nodes across (at least {CORE_NODES_MIN:g})",
                across >= CORE_NODES_MIN))
    for name in cases:
        scan = pd.read_csv(SOL/f"cycle_scan_{name}.csv").sort_values("vortex_circulation_over_Uc")
        steps = np.diff(scan["vortex_peak_vorticity_c_over_U"].values)
        out.append((f"{name}: the core vorticity of the marker rises with its circulation over the "
                    f"{len(scan)} instants scanned", bool((steps >= -1e-9).all())))
    return out


def transpiration_is_reported_with_its_sign(cases: dict) -> list:
    """The cycle minimum of the transpiration and the fraction of the cycle with inflow are those of the scan."""
    out = []
    for name in cases:
        scan, m = pd.read_csv(SOL/f"cycle_scan_{name}.csv"), _metric(name)
        vt = scan["transpiration_vt_over_U"]
        frac = float((vt < 0.0).mean())
        out.append((f"{name}: transpiration at the trailing edge from {vt.min():+.3f} U to {vt.max():+.3f} U over "
                    f"the cycle, an inflow over a fraction {frac:.3f} of it, as the metrics give",
                    abs(float(m["transpiration_vt_over_U_cycle_min"]) - vt.min()) < 1e-5
                    and abs(float(m["transpiration_inflow_fraction_of_cycle"]) - frac) < 1e-5))
    return out


def local_state_follows_from_the_corrected_pressure(cases: dict, air: flowfield.Air) -> list:
    """The corrected pressure is the Karman-Tsien correction of the incompressible one; the local Mach number
    and the temperatures stored beside it follow from it by the isentropic relations; the flag of the region
    beyond the critical pressure is where that Mach number is above one; and the size of the region is stored."""
    out = []
    g = air.gamma
    for name, c in cases.items():
        M = c["M"]
        cp_star = flowfield.sonic_values(M, air)[1]
        T0 = air.T_inf*(1.0 + 0.5*(g - 1.0)*M*M)
        tables = {"surface pressures": pd.read_csv(SOL/f"cp_distribution_{name}.csv")}
        tables.update({f.name: pd.read_csv(f) for f in sorted(SOL.glob(f"field_{name}_*_a*.csv.gz"))})
        worst = dict(kt=0.0, mach=0.0, T=0.0, Tr=0.0)
        flag_ok, n = True, 0
        for df in tables.values():
            df = df[df["inside"] == 0] if "inside" in df else df
            sample = df.iloc[::max(1, len(df)//400)]
            kt, _ = flowfield.karman_tsien(sample["Cp_incompressible"].values, M, air)
            p = np.maximum(1.0 + 0.5*g*M*M*sample["Cp"].values, 1e-12)
            mach = np.sqrt(np.maximum(5.0*((T0/air.T_inf)*p**(-(g - 1.0)/g) - 1.0), 0.0))
            T = T0/(1.0 + 0.2*mach**2)
            Tr = T*(1.0 + air.recovery*0.2*mach**2)
            soft = sample["Cp"].values > flowfield.VACUUM_CAP*flowfield.vacuum_cp(M, air) + 1e-3
            worst["kt"] = max(worst["kt"], float(np.max(np.abs(kt - sample["Cp"].values)[soft], initial=0.0)))
            worst["mach"] = max(worst["mach"], float(np.max(np.abs(mach - sample["Mach"].values)/(1.0 + mach))))
            worst["T"] = max(worst["T"], float(np.max(np.abs(T - sample["T_static_K"].values)[soft], initial=0.0)))
            worst["Tr"] = max(worst["Tr"], float(np.max(np.abs(Tr - sample["T_recovery_K"].values)[soft], initial=0.0)))
            clear = np.abs(sample["Cp"].values - cp_star) > 1e-3
            flag_ok = flag_ok and bool(((sample["beyond_critical"].values == 1) == (sample["Cp"].values < cp_star))[clear].all())
            flag_ok = flag_ok and bool(((sample["Mach"].values > 1.0) == (sample["beyond_critical"].values == 1))[clear].all())
            n += len(sample)
        out.append((f"{name}: at {n} sampled points of the surface and field files the stored Cp is the Kármán-Tsien "
                    f"correction of the incompressible one (largest difference {worst['kt']:.1e}), and the local Mach "
                    f"number, static and recovery temperature follow from it by the isentropic relations (largest "
                    f"differences {worst['mach']:.1e}, {worst['T']:.2f} K, {worst['Tr']:.2f} K)",
                    worst["kt"] < 2e-3 and worst["mach"] < ISENTROPIC_TOL and worst["T"] < TEMPERATURE_TOL_K
                    and worst["Tr"] < TEMPERATURE_TOL_K))
        out.append((f"{name}: the flag of the region beyond the critical pressure is set exactly where the corrected "
                    f"Cp is below Cp* = {cp_star:.2f}, which is where the local Mach number is above one", flag_ok))
    one = flowfield.compressible_state(np.array([cp_star]), M, air)
    out.append((f"the critical pressure coefficient and local Mach number one coincide (Mach "
                f"{float(one['Mach'][0]):.6f} at Cp*)", abs(float(one["Mach"][0]) - 1.0) < 1e-9))
    return out


def restored_figures_carry_their_statistic(cases: dict) -> list:
    """The Mach, temperature and surface figures exist, and each figure of a corrected quantity carries the size
    of the region beyond the critical pressure that the reconstruction file gives for its instant."""
    out = []
    masked = pd.read_csv(ROOT/"06_postprocessing"/"masked_region_record.csv")
    for name in cases:
        rec = pd.read_csv(SOL/f"reconstruction_{name}.csv").set_index("phase")
        maps = [p.name for q in CORRECTED_QUANTITIES for p in sorted(PLOTS.glob(f"contour_{q}_{name}_*.png"))]
        out.append((f"{name}: {len(maps)} maps of the corrected pressure, the local Mach number and the two "
                    f"temperatures exist, and the temperature profile along the surface",
                    len(maps) == len(CORRECTED_QUANTITIES)*len(PHASES)
                    and (PLOTS/f"temperature_profile_{name}.png").exists()))
        ok, n = True, 0
        for fname in maps + [f"temperature_profile_{name}.png", f"cp_distribution_{name}.png"]:
            rows = masked[masked["file"] == fname]
            ok = ok and len(rows) > 0
            for r in rows.itertuples():
                n += 1
                ok = ok and math.isclose(r.surface_fraction, rec.loc[r.phase, "beyond_critical_surface_fraction"],
                                         rel_tol=1e-9, abs_tol=1e-12)
                if fname.startswith("contour_"):
                    ok = ok and math.isclose(r.field_fraction, rec.loc[r.phase, "beyond_critical_field_fraction"],
                                             rel_tol=1e-9, abs_tol=1e-12)
        out.append((f"{name}: every one of those figures prints the part of the surface and of the field beyond the "
                    f"critical pressure at its instant ({n} entries of masked_region_record.csv, equal to the "
                    f"reconstruction file; largest part of the surface "
                    f"{100.0*rec.beyond_critical_surface_fraction.max():.1f} %)", ok and n > 0))
    missing = [f for f in SURFACE_FIGURES if not (PLOTS/f).exists()]
    out.append((f"the {len(SURFACE_FIGURES)} surface and vector figures exist" + (f" (missing {missing})" if missing else ""),
                not missing))
    return out


def fields_hold_what_they_say(cases: dict, cfg: dict) -> list:
    """Field files, the field at the mesh nodes and the configuration."""
    out = []
    fr = cfg["field_reconstruction"]
    n_panels = fr["read_by_the_driver"]["n_panels"]
    varying = {col: False for col in FIELD_COLUMNS}
    for name in cases:
        files = sorted(SOL.glob(f"field_{name}_*_a*.csv.gz"))
        out.append((f"{name}: four reconstructed fields", len(files) == len(PHASES)))
        for f in files:
            df = pd.read_csv(f)
            if list(df.columns) != FIELD_COLUMNS:
                out.append((f"{f.name}: columns are {FIELD_COLUMNS}", False))
            for col in FIELD_COLUMNS:
                varying[col] = varying[col] or df[col].nunique() > 1
        cp = pd.read_csv(SOL/f"cp_distribution_{name}.csv")
        out.append((f"{name}: the surface pressures hold {n_panels} control points an instant, the number of "
                    f"panels the configuration gives the driver",
                    bool((cp.groupby("phase").size() == n_panels).all())))
    out.append(("no column of the field files is the same at every node of every file", all(varying.values())))
    rec = fr["record_of_module_constants"]
    same = (rec["vortex_x0"] == flowfield.VORTEX_X0 and rec["vortex_dx"] == flowfield.VORTEX_DX
            and rec["vortex_y0"] == flowfield.VORTEX_Y0 and rec["vortex_dy"] == flowfield.VORTEX_DY
            and rec["vortex_core_over_c"] == flowfield.VORTEX_CORE_OVER_C
            and rec["attached_outflow_fraction"] == flowfield.ATTACHED_OUTFLOW_FRACTION
            and rec["near_wall_panel_lengths"] == flowfield.NEAR_WALL_PANEL_LENGTHS)
    out.append(("configuration: the constants it records are those of unistall.flowfield now", same))
    out.append(("configuration: the four omissions of the quasi-steady reconstruction are listed",
                fr["quasi_steady_omissions"] == list(flowfield.QUASI_STEADY_OMISSIONS)
                and len(flowfield.QUASI_STEADY_OMISSIONS) == 4))
    edge = np.zeros((5, 7), bool)
    edge[:, 0] = True
    grown = flowfield._dilate(edge, 1)
    out.append(("mask dilation does not wrap: a mask on one edge of the grid leaves the opposite edge clear",
                bool(grown[:, 1].all()) and not grown[:, -1].any()))
    return out


def mesh_field_is_what_it_says(cases: dict, air: flowfield.Air, cfg: dict) -> list:
    """Wall values are the surface values, nothing nearer the wall than the near-wall distance carries a
    field value, and the script marches nothing."""
    out = []
    script = _script("02_mesh/field_on_mesh.py")
    text = (ROOT/"02_mesh"/"field_on_mesh.py").read_text(encoding="utf-8")
    out.append(("mesh field: the script calls no solver of the loads", "dsmodel" not in text))
    for name, c in cases.items():
        df = pd.read_csv(SOL/f"field_{name}_peak_on_mesh.csv.gz")
        rec, row = _peak_reconstruction(name, c, cfg)
        near = flowfield.smoothing_radius(rec)/c["chord"]
        inner = df[(df.j > 0) & (df.wall_distance_over_c <= near)]
        wall = df[df.j == 0]
        cp_wall, _, _ = script.wall_values(rec, c, air, wall.x_over_c.values, wall.y_over_c.values)
        err = float(np.max(np.abs(wall.Cp.values - cp_wall)))
        out.append((f"{name}: mesh wall values are the surface values (largest difference {err:.1e}); the "
                    f"{len(inner)} nodes nearer the wall than {near:.4f} c carry no value; incidence "
                    f"{row.alpha_deg:.3f} deg is the stored peak row",
                    err < 1e-4 and bool(inner.Cp.isna().all()) and bool(wall.Cp.notna().all())
                    and bool((wall.source == "surface").all())))
    mesh = pd.read_csv(ROOT/"02_mesh"/"mesh_quality_metrics.csv").set_index("metric")["value"]
    nodes = pd.read_csv(ROOT/"02_mesh"/"mesh_nodes.csv")
    out.append(("mesh: node file has the stated number of nodes", len(nodes) == int(float(mesh["total_nodes"]))))
    return out


def numerics_are_reported_as_observed() -> list:
    """Step refinement, observed order and repeatability."""
    ts = pd.read_csv(CONV/"timestep_refinement.csv")
    rep = ts[ts.reported]
    cols = [c for c in ts.columns if c.startswith("pct_from_finest")]
    worst = float(rep[cols].abs().max().max())
    out = [(f"reported step within {worst:.2f} % of the finest march (limit {STEP_MAX_PCT} %)", worst <= STEP_MAX_PCT)]
    order = pd.read_csv(CONV/"timestep_order.csv")
    fresh = _script("04_solver/run_case.py").timestep_order(ts)
    same = bool(np.allclose(order["discretisation_uncertainty"], fresh["discretisation_uncertainty"], rtol=2e-2,
                            atol=1e-7) and (order["order_observed"].values == fresh["order_observed"].values).all())
    seen = order[order.order_observed]
    out.append((f"observed order of convergence: found in {len(seen)} of {len(order)} quantities"
                + (f" (between {seen.observed_order.min():.2f} and {seen.observed_order.max():.2f})" if len(seen) else "")
                + f"; largest discretisation uncertainty {order.uncertainty_pct_of_reported.max():.3f} % of the "
                f"reported value; the table is what the refinement gives", same))
    texts = "".join((REPORT/f).read_text(encoding="utf-8") for f in ("rpt_text_mid.py", "rpt_text_back.py"))
    out.append(("the report text states no order of convergence of its own (it quotes the observed-order table)",
                not re.search(r"(first|second)[ -]order(?! lag| differential)", texts)))
    return out


def error_bands_and_spread_are_current(cases: dict) -> list:
    """The held-out error bands quoted beside the case results are those of results/validation_table.csv."""
    fresh = _script("04_solver/run_case.py").held_out_error_bands()
    stored = pd.read_csv(SOL/"held_out_error_bands.csv")
    out = [(f"held-out error bands are those of results/validation_table.csv ({int(stored.n_loops.max())} loops; "
            f"peak lift {stored.mean_error.iloc[0]:+.1f} % mean, {stored.std_error.iloc[0]:.1f} % standard deviation)",
            bool(np.allclose(stored[["mean_error", "std_error"]], fresh[["mean_error", "std_error"]], atol=1e-4)))]
    sp = pd.read_csv(SOL/"static_station_spread.csv")
    out.append((f"Case B: the spread of its results between the two Mach stations of the static inputs is stored "
                f"(peak lift {sp.CL_max.iloc[-1]:.3f}, onset {sp.stall_onset_alpha_deg.iloc[-1]:.2f} deg)",
                len(sp) == 4 and abs(sp.CL_max.iloc[0] - float(_metric(pm.CASE_B)["CL_max"])) < 1e-4))
    return out


def figures_are_what_the_files_support() -> list:
    """Figures present, none of a quantity that is not derived, the measured loop drawn where there is one."""
    out = []
    for name in pm.CASES:
        n_plots = sum(len(list(PLOTS.glob(f"contour_{q}_{name}_*.png"))) for q in FIELD_QUANTITIES)
        out.append((f"{name}: {n_plots} field maps of {len(FIELD_QUANTITIES)*len(PHASES)}",
                    n_plots == len(FIELD_QUANTITIES)*len(PHASES)))
    overlay = pd.read_csv(ROOT/"06_postprocessing"/"overlay_record.csv").set_index("case")
    for name, c in pm.CASES.items():
        if c["measured_frame"]:
            out.append((f"{name}: the measured loop {c['measured_frame']} is drawn with the load loops "
                        f"(06_postprocessing/overlay_record.csv)", bool(overlay.loc[name, "drawn"])))
    return out


def documents_are_current() -> list:
    """The documents exist, their register of numbers is whole and current, figures are legible and pages full."""
    out = [(f"report: {doc} exists", (REPORT/doc).exists()) for doc in DOCUMENTS]
    reg = json.loads((REPORT/"report_numbers.json").read_text(encoding="utf-8"))
    out.append((f"report: the register counts {reg['documents']['report_numbers_quoted']} numbers and holds "
                f"{len(reg['numbers'])}", reg["documents"]["report_numbers_quoted"] == len(reg["numbers"])))
    checked, stale, tables = 0, [], {}
    for n in reg["numbers"]:
        m = re.fullmatch(r"(\w+)((?: \[[^\]=]+=[^\]]*\])+)", n["field"])
        if not (m and n["file"].endswith(".csv") and (ROOT/n["file"]).exists()):
            continue
        df = tables.setdefault(n["file"], pd.read_csv(ROOT/n["file"]))
        mask = pd.Series(True, index=df.index)
        for key, wanted in re.findall(r"\[([^\]=]+)=([^\]]*)\]", m.group(2)):
            mask &= df[key].astype(str) == wanted
        if mask.sum() != 1:                     # a figure taken over several rows: not one cell to compare
            continue
        now = df.loc[mask, m.group(1)].iloc[0]
        checked += 1
        if not (str(now) == str(n["value"]) or (isinstance(n["value"], (int, float)) and not isinstance(now, str)
                                               and math.isclose(float(now), float(n["value"]), rel_tol=1e-12))):
            stale.append(f"{n['file']}: {n['field']} is {now}, the report quotes {n['value']}")
    out.append((f"report: {checked} quoted numbers read from CSV files still equal their source values"
                + (f"; stale: {stale[:3]}" if stale else ""), checked > 0 and not stale))
    let = pd.read_csv(REPORT/"figure_lettering.csv")
    rec = let[let.recorded]
    small = rec[rec.smallest_lettering_as_printed_pt < LETTERING_MIN_PT]
    out.append((f"report: the smallest lettering of the {len(rec)} figures of the case study is "
                f"{rec.smallest_lettering_as_printed_pt.min():.1f} pt as printed (at least {LETTERING_MIN_PT:g} pt); "
                f"the {int((~let.recorded).sum())} figures of results/figures keep no record and are not checked"
                + (f"; too small: {list(small.file)}" if len(small) else ""), len(rec) > 0 and small.empty))
    fill = pd.read_csv(REPORT/"page_fill.csv")
    short = fill[fill.short]
    out.append((f"documents: of {len(fill)} pages, {len(short)} are less than half full without being the end of a "
                f"part" + (f" ({[(r.document, r.page) for r in short.itertuples()][:6]})" if len(short) else ""),
                short.empty))
    return out


def drawings_state_what_they_show() -> list:
    """Section A-A is on the side its arrows require, the body is one object, every sheet states its scale."""
    rec = pd.read_csv(SHEETS/"sheet_record.csv")
    val = rec.set_index(["sheet", "item"])["value"]
    drawn, due = val[(3, "section A-A: leading edge drawn on the")], val[(3, "section A-A: leading edge belongs on the")]
    draw = _script("08_engineering_drawings/draw_engineering.py")
    out = [(f"sheet 3: section A-A is seen looking {val[(3, 'section A-A: cutting-plane arrows point')]}, which "
            f"puts the leading edge on the {due}; it is drawn on the {drawn}",
            drawn == due == draw.section_side(val[(3, "section A-A: cutting-plane arrows point")] == "outboard")
            and draw.section_side(True) != draw.section_side(False))]
    body = rec[rec.item.str.startswith("generic body")].pivot(index="item", columns="sheet", values="value")
    out.append((f"sheets 1 and 2: the generic body has the same mast and fin heights in rotor radii "
                f"({dict(body[1].astype(float))})", bool((body[1].astype(float) == body[2].astype(float)).all())))
    scales = rec[rec.item == "scale stated in the title block"]
    out.append((f"every one of the {len(scales)} sheets states its scale or that it has none",
                len(scales) == 4 and bool(scales.value.str.len().gt(0).all())))
    dims = pd.read_csv(SHEETS/"dimensions.csv").set_index("name")["value"]
    out.append((f"drawings: the root end of the blade is one value, {dims['Root end of the blade as drawn, fraction of rotor radius']} R, "
                f"listed in dimensions.csv", float(dims["Root end of the blade as drawn, fraction of rotor radius"])
                == pm.ROOT_END_OVER_R))
    return out


def words_match_the_standing() -> list:
    """Landing page, folder notes and run commands."""
    index = (ROOT/"index.md").read_text(encoding="utf-8")
    definition = (ROOT/"00_overview"/"case_definition.md").read_text(encoding="utf-8")
    cite = _script("03_model_setup/generate_setup.py").EXPERIMENT_CITATION
    out = [("landing page: says that Case A is set beside a calibration loop and that the fields are a "
            "reconstruction", "calibration loop" in index and "reconstruction" in index),
           ("the experiment report is cited with one string, authors as in docs/references.bib, on the landing "
            "page and in the case definition", cite in index and cite in definition)]
    for folder in ("02_mesh", "05_solution", "06_postprocessing"):
        first = (ROOT/folder/"README.md").read_text(encoding="utf-8").splitlines()[0]
        out.append((f"{folder}: the folder note opens with its standing", first.startswith("**Standing:")))
    bare = [rel for rel, text in _sources().items()
            if re.search(r"^# Run from the repository root:\s+python3", text, re.M)]
    bare += ["04_solver/README.md"] if re.search(r"`python3 0", (ROOT/"04_solver"/"README.md").read_text()) else []
    out.append(("every run command printed at the head of a script sets PYTHONPATH" + (f" (not in {bare})" if bare else ""),
                not bare))
    return out


def validation_folder_is_the_results(cases: dict) -> list:
    """06_postprocessing/validation holds the numbers of results/, and the loops its rule gives."""
    mv = _script("06_postprocessing/validation/make_validation.py")
    names = ("validation_static.csv", "fig_validation_static.png", "validation_loops_summary.csv",
             "validation_targets.csv", "fig_validation_loops.png", "calibration_constants.csv", "PROVENANCE.md",
             "README.md")
    missing = [n for n in names if not (VALIDATION/n).exists()]
    out = [("validation: every file of the folder is present" + (f" (missing {missing})" if missing else ""),
            not missing)]
    table = mv.scored_loops()
    summary = pd.read_csv(VALIDATION/"validation_loops_summary.csv")
    same = (len(summary) == len(table) and list(summary["frame"]) == list(table["frame"])
            and bool(np.allclose(summary.select_dtypes("number"), table[summary.columns].select_dtypes("number"),
                                 equal_nan=True)))
    out.append((f"validation: the summary holds the {len(table)} rows of results/validation_table.csv for the "
                f"{mv.MODEL} model", same))
    targets, source = pd.read_csv(VALIDATION/"validation_targets.csv"), pd.read_csv(ROOT/"results"/"validation_targets.csv")
    out.append((f"validation: the targets file is results/validation_targets.csv ({int((~source['met']).sum())} of "
                f"{len(source)} targets not met, as that file gives them)", targets.equals(source)))
    chosen = [frame for frame, _ in mv.chosen_loops(table)]
    have = sorted({p.name[len("loop_"):-len("_CL.csv")] for p in VALIDATION.glob("loop_*_CL.csv")})
    complete = all((VALIDATION/f"loop_{f}_{q}.csv").exists() for f in chosen for q in ("CL", "CM", "CD"))
    out.append((f"validation: the loops drawn are the ones the rule gives ({', '.join(chosen)})",
                have == sorted(set(chosen)) and complete
                and chosen[-1] == pm.CASES[pm.CASE_A]["measured_frame"]))
    cal = json.loads((ROOT/"results"/"calibrated_constants.json").read_text(encoding="utf-8"))
    cc = pd.read_csv(VALIDATION/"calibration_constants.csv").set_index("constant")["value"]
    out.append(("validation: the constants listed are those of results/calibrated_constants.json",
                all(math.isclose(float(cc[k]), float(v)) for k, v in {**cal["constants"], **cal["fixed"]}.items())))
    return out


def solver_entry_and_static_polar(cases: dict, air: flowfield.Air, cfg: dict) -> list:
    """The single-condition entry returns the packaged solver's cycle, and the model's static polar is stored
    and passes through the static-stall point of its inputs."""
    entry = _script("04_solver/unistall_solver.py")
    c = cases[pm.CASE_A]
    table = entry.cycle_table(c["a_mean"], c["a_amp"], c["k"], c["M"])
    o = dsmodel.solve(c["a_mean"], c["a_amp"], c["k"], c["M"], a_sound=air.a_sound)
    worst = max(float(np.max(np.abs(table[col].values - o[key]))) for col, key in entry.COLUMNS)
    out = [(f"04_solver/unistall_solver.py: its cycle for Case A is that of unistall.dsmodel.solve (largest "
            f"difference {worst:.1e})", len(table) == len(o["CL"]) and worst <= SOLVER_ENTRY_TOL)]
    for name in cases:
        path = SOL/f"model_static_polar_{name}.csv"
        if not path.exists():
            out.append((f"{name}: the static polar of the model is stored", False))
            continue
        polar = pd.read_csv(path)
        si = cfg["cases"][name]["static_inputs"]
        a_ss = si["alpha_stall_deg"]
        cl, cc = (float(np.interp(a_ss, polar["alpha_deg"], polar[q])) for q in ("CL", "CC"))
        due = si["CN1_onset"]*math.cos(math.radians(a_ss)) + cc*math.sin(math.radians(a_ss))
        pct = 100.0*(cl - due)/due
        out.append((f"{name}: the slow-sweep lift of the model at the static-stall incidence {a_ss:.2f} deg is "
                    f"{cl:.3f}, {pct:+.1f} % from the static input C_N1 resolved to lift with the model's chord "
                    f"force there, {due:.3f} (tolerance {POLAR_STALL_TOL_PCT} %); the file says it is not a "
                    f"measurement",
                    abs(pct) <= POLAR_STALL_TOL_PCT and bool(polar["note"].str.contains("not a measurement").all())
                    and {"alpha_deg", "CL", "CD", "CM_c4", "CN", "CC", "f_separation", "outside_static_data"}
                    <= set(polar.columns)))
    return out


def folder_list_is_complete() -> list:
    """Every file of the numbered folders is covered by an entry of 00_overview/folder_contents.md."""
    fc = _script("00_overview/write_folder_contents.py")
    loose = fc.uncovered()
    listed = (ROOT/"00_overview"/"folder_contents.md").read_text(encoding="utf-8")
    absent = [e[0] for e in fc.ENTRIES if e[0].split("/", 1)[1] not in listed]
    return [(f"00_overview/folder_contents.md: every one of the {len(fc.files())} files of the numbered folders is "
             f"covered by one of its {len(fc.ENTRIES)} entries" + (f" (not covered: {loose[:5]})" if loose else ""),
             not loose and not absent)]


def run() -> list:
    """Every (description, passed) pair."""
    cases, air, cfg = pm.read_setup()
    groups = (setup_agrees_with_the_definition(cases, air, cfg), case_b_follows_from_the_rotor(cases, air),
              one_speed_of_sound(air), one_outline(cfg), loads_come_from_the_solver(cases, air),
              metrics_follow_from_the_histories(cases), reconstruction_is_current_and_consistent(cases, air, cfg),
              reconstruction_is_converged_in_panels(cases, cfg), vortex_marker_follows_the_clock(cases, cfg),
              transpiration_is_reported_with_its_sign(cases), fields_hold_what_they_say(cases, cfg),
              local_state_follows_from_the_corrected_pressure(cases, air), restored_figures_carry_their_statistic(cases),
              mesh_field_is_what_it_says(cases, air, cfg), numerics_are_reported_as_observed(),
              error_bands_and_spread_are_current(cases), figures_are_what_the_files_support(),
              documents_are_current(), drawings_state_what_they_show(), words_match_the_standing(),
              validation_folder_is_the_results(cases), solver_entry_and_static_polar(cases, air, cfg),
              folder_list_is_complete())
    return [item for group in groups for item in group]


def main() -> int:
    results = run()
    for text, ok in results:
        print(("  ok    " if ok else "  FAIL  ") + text)
    failed = [t for t, ok in results if not ok]
    print(f"[check] {len(results) - len(failed)} of {len(results)} checks passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
