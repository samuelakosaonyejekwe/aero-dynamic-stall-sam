# -*- coding: utf-8 -*-
"""
07_report / rpt_text_back.py
----------------------------
The last part of the report: results of the two cases, the reconstructed flow
field and its limits, the grid, the response over mean incidence and reduced
frequency, the drawings, the limits of the study and the references. Every
number in the text is read from a file through rpt_data.Data, which registers
it.

Author: Akosa Samuel Onyejekwe (independent)

Units: as in the files read; angles in degrees, temperatures in kelvin.
"""
import project_meta as meta
from rpt_data import cell, frame_rows, sibling
from unistall import flowfield

front = sibling("rpt_text_front")
mid = sibling("rpt_text_mid")

CFG = front.CFG
VTABLE = front.VTABLE
DIMS = front.DIMS
PLOTS = mid.PLOTS
MESH = "02_mesh/mesh_quality_metrics.csv"
SURFACE = "05_solution/response_surface.csv"
SHEETS = "08_engineering_drawings"
SHEET_RECORD = "08_engineering_drawings/sheet_record.csv"
OVERLAY = "06_postprocessing/overlay_record.csv"
PANELS = "05_solution/convergence/panel_convergence_summary.csv"
DAMPING = "results/validation_damping.csv"

# metric, label, format, quantity of the held-out error bands (or "")
LOAD_METRICS = (
    ("CL_max", "Peak lift coefficient", "{:.2f}", "CL_max"),
    ("alpha_at_CLmax_deg", "Incidence at peak lift, deg", "{:.1f}", "alpha_at_CLmax_deg"),
    ("CM_min_c4", "Minimum moment coefficient about the quarter chord", "{:.2f}", "CM_min_c4"),
    ("alpha_at_CMmin_deg", "Incidence at minimum moment, deg", "{:.1f}", ""),
    ("CD_max", "Peak drag coefficient", "{:.2f}", "CD_max"),
    ("alpha_at_CDmax_deg", "Incidence at peak drag, deg", "{:.1f}", ""),
    ("CD_min", "Lowest drag coefficient", "{:.3f}", ""),
    ("stall_onset_alpha_deg", "Incidence at which the onset criterion is passed, deg", "{:.1f}", ""),
    ("vortex_sheddings_per_cycle", "Vortex sheddings per cycle, as the load model counts them", "{:d}", ""),
    ("vortex_clock_runs_per_cycle", "Runs of the vortex clock in the cycle (the ramps of the states figure)",
     "{:d}", ""),
    ("CN_vortex_max", "Greatest vortex normal force", "{:.2f}", ""),
    ("alpha_at_CN_vortex_max_deg", "Incidence at greatest vortex normal force, deg", "{:.1f}", ""),
    ("cycle_damping_Xi", "Cycle damping (model output only)", "{:.3f}", ""),
)
FIELD_METRICS = (
    ("Cp_critical", "Critical pressure coefficient at the stream Mach number", "{:.2f}"),
    ("Cp_vacuum", "Pressure coefficient of a vacuum at the stream Mach number", "{:.2f}"),
    ("total_temperature_K", "Total temperature of the stream, K", "{:.2f}"),
    ("recovery_factor", "Recovery factor", "{:.3f}"),
    ("Cp_surface_min_incompressible_stored_instants", "Lowest incompressible surface pressure coefficient over the "
     "stored instants", "{:.1f}"),
    ("Cp_surface_min_stored_instants", "Lowest corrected surface pressure coefficient over the stored instants "
     "(held at the vacuum-side bound)", "{:.1f}"),
    ("beyond_critical_surface_fraction_worst_stored_instants", "Largest part of the surface beyond the critical "
     "pressure over the stored instants, fraction", "{:.3f}"),
    ("beyond_critical_surface_fraction_worst_over_cycle", "The same over the whole cycle", "{:.3f}"),
    ("beyond_critical_fraction_of_cycle", "Fraction of the cycle with any of the surface beyond the critical "
     "pressure", "{:.2f}"),
    ("beyond_critical_field_area_over_c2_worst_stored_instants", "Largest area of the field beyond the critical "
     "pressure, over c²", "{:.4f}"),
    ("beyond_critical_field_fraction_worst_stored_instants", "The same as a fraction of the stored field", "{:.5f}"),
    ("Cp_closure_error_pct", "Lift from the stored (corrected) surface pressure at peak lift, difference from the "
     "lift given, %", "pow10"),
    ("circulation_factor_at_CLmax", "Fraction of the circulation of the lift the sheet carries for that, at peak "
     "lift", "{:.3f}"),
    ("prandtl_glauert_factor", "Prandtl-Glauert factor the fraction starts from", "{:.3f}"),
    ("Cp_closure_error_pct_prandtl_glauert_factor", "Lift difference of the corrected pressure with the circulation "
     "fixed at the Prandtl-Glauert factor, at peak lift, %", "{:+.2f}"),
    ("Cp_closure_error_pct_incompressible", "Lift difference of the incompressible pressure with the whole "
     "circulation, at peak lift, %", "{:+.2f}"),
    ("Cp_closure_error_pct_panels_only", "The same with the panels alone (solid wall, no vortex), %", "{:+.2f}"),
    ("Cp_closure_worst_abs_pct_incompressible_stored_instants", "Largest lift difference of the incompressible "
     "pressure over the stored instants, %", "{:.2f}"),
    ("circulation_factor_cycle_min", "Smallest fraction of the circulation over the cycle", "{:.3f}"),
    ("circulation_factor_cycle_max", "Largest fraction of the circulation over the cycle", "{:.3f}"),
    ("Cp_moment_minus_model_at_CLmax", "Pitching moment of the surface pressure less the model's, at peak lift",
     "{:+.3f}"),
    ("Cp_drag_minus_model_at_CLmax", "Pressure drag of the surface pressure less the model's drag, at peak lift",
     "{:+.3f}"),
    ("Cp_moment_minus_model_worst_abs_stored_instants", "Largest moment difference over the stored instants",
     "{:.3f}"),
    ("Cp_drag_minus_model_worst_abs_stored_instants", "Largest drag difference over the stored instants", "{:.3f}"),
    ("Cp_moment_sign_agrees_with_model_instants", "Stored instants at which the two moments have the same sign",
     "{:.0f}"),
    ("stored_instants", "Stored instants", "{:.0f}"),
    ("Cp_trailing_edge_jump_worst_stored_instants", "Largest difference of pressure coefficient across the "
     "trailing edge over the stored instants", "pow10"),
    ("transpiration_vt_over_U_cycle_max", "Greatest transpiration speed at the trailing edge over the cycle, over U "
     "(outflow)", "{:.3f}"),
    ("transpiration_vt_over_U_cycle_min", "Least transpiration speed over the cycle, over U (negative: inflow)",
     "{:+.3f}"),
    ("transpiration_inflow_fraction_of_cycle", "Fraction of the cycle with inflow", "{:.3f}"),
    ("transpiration_flux_over_Uc_cycle_max", "Greatest flux through the suction surface over the cycle, over U c",
     "{:.3f}"),
    ("vortex_circulation_over_Uc", "Circulation of the vortex marker at its greatest, over U c", "{:.3f}"),
    ("vortex_peak_vorticity_c_over_U", "Vorticity at the centre of the vortex marker at its greatest, times c / U",
     "{:.1f}"),
    ("vortex_centre_Cp", "Incompressible pressure coefficient at the centre of the vortex marker", "{:.2f}"),
)
ASSUMED_METRICS = (
    ("vortex_x_over_c", "Vortex marker position x/c at its greatest (assumed path)", "{:.2f}"),
    ("vortex_y_over_c", "Vortex marker position y/c at its greatest (assumed path)", "{:.2f}"),
    ("vortex_core_radius_over_c", "Core radius of the vortex marker, over c (assumed)", "{:.2f}"),
)
INSTANT_FLOW = (
    ("phase", "Instant"), ("alpha_deg", "Incidence, deg"), ("f_separation", "Separation point, x/c"),
    ("transpiration_vt_over_U", "Transpiration at the edge, over U"),
    ("transpiration_flux_over_Uc", "Flux through the surface, over U c"),
    ("circulation_bound_over_Uc", "Circulation of the sheet, over U c"),
    ("circulation_total_over_Uc", "Total circulation, over U c"),
    ("dCp_under_vortex", "Change of Cp beneath the vortex"),
)
INSTANT_STATE = (
    ("phase", "Instant"), ("Cp_surface_min_incompressible", "Lowest surface Cp, incompressible"),
    ("Cp_surface_min", "Lowest surface Cp, corrected"),
    ("beyond_critical_surface_fraction", "Surface beyond critical, fraction"),
    ("beyond_critical_field_area_over_c2", "Field beyond critical, over c²"),
    ("beyond_critical_field_fraction", "Field beyond critical, fraction"),
    ("Mach_surface_max_subcritical", "Highest surface Mach number outside that region"),
    ("T_static_surface_min_subcritical_K", "Lowest surface static temperature outside it, K"),
    ("T_recovery_surface_max_K", "Highest surface recovery temperature, K"),
)
INSTANT_LOADS = (
    ("phase", "Instant"), ("circulation_factor", "Fraction of the circulation carried"),
    ("closure_error_pct", "Lift difference, corrected pressure, %"),
    ("closure_error_pct_prandtl_glauert_factor", "The same at the Prandtl-Glauert factor, %"),
    ("closure_error_pct_incompressible", "Lift difference, incompressible, %"),
    ("closure_error_pct_panels_only", "Lift difference, panels alone, %"),
    ("CM_model", "Moment, model"), ("CM_from_Cp", "Moment, surface pressure"),
    ("CD_model", "Drag, model"), ("CD_from_Cp", "Pressure drag, surface pressure"),
)
RECON = flowfield.STANDING_NOTE


def reconstruction_file(name: str) -> str:
    """The per-instant file of the reconstruction of a case."""
    return f"05_solution/reconstruction_{name}.csv"


def metrics_file(name: str) -> str:
    """The metrics file of a case."""
    return f"05_solution/metrics_{name}.csv"


def m(D, name: str, metric: str, fmt: str) -> str:
    """One metric of a case, as text."""
    return D.c(metrics_file(name), "value", fmt, metric=metric)


def plot(D, pattern: str):
    """The one plot of 06_postprocessing/plots that matches a pattern."""
    found = sorted((D.root/PLOTS).glob(pattern))
    if len(found) != 1:
        raise FileNotFoundError(f"{PLOTS}/{pattern}: {len(found)} files")
    return found[0]


def measured_loop(R, D, name: str) -> None:
    """For a case with a measured frame: what the comparison is and is not,
    and whether the measured points are in the figures (from the record the
    plotting stage writes)."""
    frame = meta.CASES[name]["measured_frame"]
    if not frame:
        R.p("No measurement is compared with in this case.")
        return

    def v(col: str, fmt: str) -> str:
        return D.c(VTABLE, col, fmt, model="tabulated", frame=frame)

    drawn = str(D.c(OVERLAY, "drawn", "{}", case=name)) == "True"
    if drawn:
        first = (f"The loops are drawn with the measured points of {frame.replace('_', ' ')} of NASA TM-84245. The "
                 f"curve is the model at the nominal condition of the case; the points were measured at Mach "
                 f"{D.c(OVERLAY, 'M', '{:.3f}', case=name)}, reduced frequency "
                 f"{D.c(OVERLAY, 'k', '{:.3f}', case=name)}, incidence "
                 f"{D.c(OVERLAY, 'alpha_mean_deg', '{:.1f}', case=name)}° ± "
                 f"{D.c(OVERLAY, 'alpha_amp_deg', '{:.1f}', case=name)}° ({OVERLAY}); both sets of conditions are "
                 f"printed in the legend of each loop.")
    else:
        first = (f"The loops are drawn WITHOUT the measured points of {frame.replace('_', ' ')}: "
                 f"{D.c(OVERLAY, 'note', '{}', case=name)} ({OVERLAY}).")
    R.p(first + f" That frame is a calibration loop, so agreement with it is not an independent test. In the "
        f"accuracy study the frame is run at its own measured conditions (Mach {v('M', '{:.3f}')}, reduced "
        f"frequency {v('k', '{:.3f}')}, incidence {v('alpha0_deg', '{:.1f}')}° ± {v('amp_deg', '{:.1f}')}°). "
        f"There its loop errors are "
        f"{v('nRMS_CL', '{:.3f}')} in lift, {v('nRMS_CM', '{:.3f}')} in moment and {v('nRMS_CD', '{:.3f}')} in "
        f"drag; the predicted peak lift is {v('CLmax_model', '{:.3f}')} against a measured "
        f"{v('CLmax_exp', '{:.3f}')}, and the predicted minimum moment {v('CMmin_model', '{:.3f}')} against "
        f"{v('CMmin_exp', '{:.3f}')} ({VTABLE}).")


def case_results(R, D, name: str) -> None:
    """Section: the results of one case."""
    case = meta.CASES[name]
    short = mid.short_case(name)
    R.h1(f"Results: {case['label']}", name)
    R.p(f"{short} is the NACA {meta.SECTION} section at {front.case_line(D, name)}, on a chord of "
        f"{D.c(front.FLOW, name, '{:.3f}', parameter='chord_c')} m.")
    if name == meta.CASE_B:
        R.note("<b>Case B is not a rotor calculation.</b> It is an aerofoil at a rotor-blade-station condition "
               "in a steady stream with prescribed pitch, at the lowest Mach number the station sees round the "
               f"azimuth and with an assumed incidence history (Section {R.sec('solved')}). Its static inputs are "
               f"a blend of two Mach stations, and its results are uncertain within the spread given in Section "
               f"{R.sec('static')}. Nothing here stands for a rotor in flight.")
    b = front.band
    R.p(f"Lift rises above its static maximum and reaches {m(D, name, 'CL_max', '{:.2f}')} [{b(D, 'CL_max')}] at "
        f"an incidence of {m(D, name, 'alpha_at_CLmax_deg', '{:.1f}')}° [{b(D, 'alpha_at_CLmax_deg', '°')}]. The "
        f"onset criterion of the vortex is passed at {m(D, name, 'stall_onset_alpha_deg', '{:.1f}')}°, and the "
        f"vortex normal force reaches {m(D, name, 'CN_vortex_max', '{:.2f}')} at "
        f"{m(D, name, 'alpha_at_CN_vortex_max_deg', '{:.1f}')}°. The moment about the quarter chord falls to "
        f"{m(D, name, 'CM_min_c4', '{:.2f}')} [{b(D, 'CM_min_c4')}] at "
        f"{m(D, name, 'alpha_at_CMmin_deg', '{:.1f}')}° and the drag reaches {m(D, name, 'CD_max', '{:.2f}')} "
        f"[{b(D, 'CD_max')}] at {m(D, name, 'alpha_at_CDmax_deg', '{:.1f}')}°; the lowest drag over the cycle is "
        f"{m(D, name, 'CD_min', '{:.3f}')}. {front.band_note(D)} The vortex clock runs "
        f"{m(D, name, 'vortex_clock_runs_per_cycle', 'word')} times in the cycle, which is the number of "
        f"sheddings the load model counts ({m(D, name, 'vortex_sheddings_per_cycle', '{:d}')}). The cycle damping "
        f"the model gives is {m(D, name, 'cycle_damping_Xi', '{:.2f}')}; as Section {R.sec('accuracy')} shows, the "
        f"model does not predict cycle damping through stall, so that figure is model output and says nothing "
        f"about stall flutter. {R.next_tab()} lists the scalar results.")
    D.whole(metrics_file(name), "the load metrics")
    R.table(["Quantity", short, "Held-out error of the model in this quantity (model minus measured)"],
            [[label, m(D, name, key, fmt), b(D, q, "°" if "deg" in q else "") if q else "no measured counterpart "
              "scored"] for key, label, fmt, q in LOAD_METRICS],
            f"Scalar results of {short} ({metrics_file(name)}), with the held-out error of the model where the "
            f"quantity is scored against measurement ({front.BANDS}).")
    measured_loop(R, D, name)
    R.p(f"{R.next_fig()} shows the loads through the reported cycle, {R.next_fig(2)} to {R.next_fig(4)} the "
        f"loops of lift, moment and drag against incidence, and {R.next_fig(5)} the states of the model: the "
        "static and the lagged separation point, the vortex normal force, the lagged normal force as a "
        "fraction of its onset level and the vortex clock as a fraction of the travel time over the chord.")
    R.fig(plot(D, f"timehist_loads_{name}.png"), f"{short}: lift, drag and moment through the reported cycle.",
          width=0.8, max_h_in=6.4)
    for tag, load in (("cl", "lift"), ("cm", "moment about the quarter chord"), ("cd", "drag")):
        R.fig(plot(D, f"hyst_{tag}_{name}.png"), f"{short}: loop of {load} against incidence.", width=0.7,
              max_h_in=5.2)
    R.fig(plot(D, f"states_{name}.png"), f"{short}: states of the load model through the reported cycle.",
          width=0.85, max_h_in=6.6)


def field_description(R, D) -> None:
    """What the reconstruction is."""
    fr = "field_reconstruction"
    drv = f"{fr}.read_by_the_driver"
    nx, ny = D.j(CFG, f"{drv}.grid_nx", "{:d}"), D.j(CFG, f"{drv}.grid_ny", "{:d}")
    box = [D.j(CFG, f"{drv}.domain_chords.{i}", "{:.1f}") for i in range(4)]
    rec = f"{fr}.record_of_module_constants"
    a, b = list(meta.CASES)
    R.h1("The reconstructed flow field", "field")
    R.note("<b>The flow fields are a reconstruction, not a flow solution.</b> " + RECON + " It is given the lift "
           "only: its surface pressure does not return the model's moment or drag.")
    R.p(f"The reconstruction (unistall/flowfield.py) is built from these parts ({CFG}).")
    R.bullets([
        f"<b>The section and its lift.</b> The closed section is {D.j(CFG, f'{drv}.n_panels', '{:d}')} straight "
        "panels, each carrying a source of constant strength and a vortex sheet of one uniform strength that "
        "carries the circulation of the predicted lift. The velocity of every panel comes from the closed-form "
        "integrals of the constant-strength panel, the same in the wall condition, on the surface and in the "
        "field.",
        "<b>A transpiration for the separated region.</b> The predicted lift is not, in general, the lift at "
        "which potential flow leaves the trailing edge smoothly; in stall it is far below it. Flow is therefore "
        "passed through the suction-surface panels aft of the separation point of the load model, rising "
        "smoothly from nothing there to its full strength at the trailing edge. Its strength is the single "
        "number that makes the pressure the same at the two control points next to the trailing edge (the "
        "Kutta condition) while the circulation stays that of the predicted lift. It takes either sign: an "
        "outflow where the predicted lift is below the Kutta value and an inflow, a suction into the wall, "
        "where it is above it. Where the flow of the model is attached the transpiration is confined to the "
        f"last {D.j(CFG, f'{rec}.attached_outflow_fraction', '{:.2f}')} of the suction surface.",
        "<b>The vortex marker.</b> A single Lamb-Oseen vortex marks the vortex clock of the load model. It "
        "carries the circulation of the predicted vortex normal force and moves along a straight path that "
        "puts it over the trailing edge when the clock says the vortex leaves the chord, and downstream of it "
        f"afterwards. Its path and its core radius ({D.j(CFG, f'{rec}.vortex_core_over_c', '{:.2f}')} of the "
        "chord) are assumed.",
        "<b>A linearised correction for compressibility.</b> The incompressible pressure coefficient is "
        "corrected by the Kármán-Tsien rule at the stream Mach number of the case, and the sheet is given the "
        "fraction of the circulation of the lift at which the corrected surface pressure returns the lift "
        "given; the fraction is solved at each instant, starting from the Prandtl-Glauert factor. From the "
        "corrected pressure the local Mach number, the static temperature and the recovery temperature follow "
        "by the isentropic relations with the stagnation state of the stream; the recovery factor is "
        f"{D.j(CFG, f'{fr}.compressibility.recovery_factor', '{:.3f}')} "
        f"({D.j(CFG, f'{fr}.compressibility.recovery_factor_source', '{}')}). This is a correction applied to an "
        "incompressible reconstruction, not a compressible solution, and it holds only where the local flow is "
        f"subcritical. Where the corrected pressure coefficient is below the critical value, "
        f"{m(D, a, 'Cp_critical', '{:.2f}')} in {mid.short_case(a)} and {m(D, b, 'Cp_critical', '{:.2f}')} in "
        f"{mid.short_case(b)}, the local Mach number it gives is above one and the values are not physical: "
        "that region is hatched in every figure and its size is given at each instant. The corrected value is "
        f"not allowed below {D.j(CFG, f'{fr}.compressibility.vacuum_cap', '{:.2f}')} of the vacuum value.",
    ])
    n_inst = D.note(len(D.json(CFG)[fr]["field_phases"]), CFG, f"number of entries of {fr}.field_phases", "word")
    R.p(f"The incompressible velocity, the pressure coefficient (corrected and incompressible), the local Mach "
        f"number, the static and recovery temperatures, the vorticity and the flag of the region beyond the "
        f"critical pressure are stored on a rectangular grid of {nx} by {ny} points over x/c from {box[0]} to "
        f"{box[1]} and y/c from {box[2]} to {box[3]}, at {n_inst} instants of the cycle.")


def footprint(D, name: str) -> str:
    """What the vortex marker does to the surface pressure beneath it in one case."""
    rel = reconstruction_file(name)
    df = D.csv(rel)
    over = df[df["vortex_over_chord"] == 1]
    n_over = D.note(len(over), rel, "number of rows with vortex_over_chord = 1", "word")
    if over.empty:
        return f"in {mid.short_case(name)} the marker is over the chord at none of the stored instants"
    raised = D.note(int((over["dCp_under_vortex"] > 0).sum()), rel,
                    "number of rows with vortex_over_chord = 1 and dCp_under_vortex > 0", "word")
    lo = D.note(over["dCp_under_vortex"].min(), rel, "smallest dCp_under_vortex", "{:+.2f}")
    hi = D.note(over["dCp_under_vortex"].max(), rel, "largest dCp_under_vortex", "{:+.2f}")
    return (f"in {mid.short_case(name)} the marker is over the chord at {n_over} of the stored instants and "
            f"raises the surface pressure coefficient beneath it at {raised} of them (change from {lo} to {hi})")


def field_limits(R, D) -> None:
    """The limits of the reconstruction, with what the files show of each."""
    a, b = list(meta.CASES)

    def peak(name: str, col: str, fmt: str) -> str:
        return D.c(reconstruction_file(name), col, fmt, phase="peak")

    def inflow(name: str) -> str:
        return (f"{m(D, name, 'transpiration_vt_over_U_cycle_min', '{:+.3f}')} U at the trailing edge at its lowest "
                f"and an inflow over a fraction {m(D, name, 'transpiration_inflow_fraction_of_cycle', '{:.2f}')} of "
                f"the cycle in {mid.short_case(name)}")

    omit = D.json(CFG)["field_reconstruction"]["quasi_steady_omissions"]
    D.note(len(omit), CFG, "number of entries of field_reconstruction.quasi_steady_omissions", "{:d}")
    R.h2("Limits of the reconstruction")
    R.bullets([
        "<b>It is potential flow round the predicted lift.</b> There is no boundary layer and no separated "
        "shear layer.",
        "<b>It is quasi-steady.</b> Each instant is a steady flow at that instant's incidence and lift, with "
        "steady Bernoulli pressure: " + "; ".join(omit) + ". The circulation the sheet carries is that of the "
        "lift; with the vortex marker present the total circulation round section and marker is larger, and "
        "the two are tabulated side by side below.",
        f"<b>The compressibility correction is linearised and holds only for subcritical flow.</b> At peak "
        f"lift the lowest incompressible surface pressure coefficient is "
        f"{peak(a, 'Cp_surface_min_incompressible', '{:.1f}')} in {mid.short_case(a)} and "
        f"{peak(b, 'Cp_surface_min_incompressible', '{:.1f}')} in {mid.short_case(b)}; corrected, the pressure is "
        f"below the critical value over a fraction {peak(a, 'beyond_critical_surface_fraction', '{:.3f}')} and "
        f"{peak(b, 'beyond_critical_surface_fraction', '{:.3f}')} of the surface at the leading edge, and over "
        f"{peak(a, 'beyond_critical_field_area_over_c2', '{:.4f}')} c² and "
        f"{peak(b, 'beyond_critical_field_area_over_c2', '{:.4f}')} c² of the field. Inside that region the "
        f"pressure, the local Mach number and the temperatures are not physical: no shock and no supersonic "
        f"region is represented. The region is hatched in the figures and its size tabulated per instant; some "
        f"of the surface is in it over a fraction {m(D, a, 'beyond_critical_fraction_of_cycle', '{:.2f}')} of the "
        f"cycle in {mid.short_case(a)} and {m(D, b, 'beyond_critical_fraction_of_cycle', '{:.2f}')} in "
        f"{mid.short_case(b)}.",
        "<b>The temperatures are isentropic estimates.</b> They follow from the corrected pressure with a "
        "uniform total temperature; there is no heat conduction. The recovery temperature is that of an "
        "adiabatic wall under a turbulent boundary layer, which the reconstruction does not have; in a field "
        "map it is only what such a wall would take at the local state.",
        f"<b>The transpiration is an equivalent one, and it is not always an outflow.</b> It is chosen to "
        f"satisfy the Kutta condition, and its shape is chosen for illustration. Over the cycle it reaches "
        f"{inflow(a)}, and {inflow(b)}. It is not a model of the flow inside a separated region.",
        "<b>The separated region is not drawn.</b> The transpiration produces no reversed flow, no shear layer "
        "and no wake, so the maps at the instants of deep stall look like attached flow. They are not pictures "
        "of a stalled flow; each map marks the separation point of the load model on the section.",
        f"<b>The vortex is a marker of the vortex clock, not a model of the dynamic-stall vortex.</b> A "
        f"clockwise vortex above the surface slows the surface flow beneath it and so raises the pressure "
        f"there, where a dynamic-stall vortex leaves a suction footprint: {footprint(D, a)}; {footprint(D, b)}. "
        f"Its path and core radius are assumed.",
        "<b>The surface pressure does not integrate exactly to the lift it was given, and does not return the "
        "model's moment or drag</b> (next paragraphs).",
        "<b>The vorticity maps show the singularities, not a vorticity field.</b> The vorticity of the "
        "reconstruction is zero everywhere except in the core of the vortex marker; the bound sheet lies on the "
        "surface itself. No boundary-layer or wake vorticity is in them, and no sign information about a real "
        "flow should be read from them.",
        "<b>Not compared with measurement.</b> No measured pressure or velocity is set beside any of it.",
    ])
    R.p(f"<b>Closure of the lift.</b> The model's lift is the lift of a compressible flow, so the pressure "
        f"compared with it is the corrected one. With the whole circulation of the lift on the sheet, the "
        f"incompressible surface pressure returns the lift at peak lift within "
        f"{m(D, a, 'Cp_closure_error_pct_incompressible', '{:+.2f}')} % in {mid.short_case(a)} and "
        f"{m(D, b, 'Cp_closure_error_pct_incompressible', '{:+.2f}')} % in {mid.short_case(b)} (panels alone, a "
        f"solid wall and no vortex: {m(D, a, 'Cp_closure_error_pct_panels_only', '{:+.2f}')} % and "
        f"{m(D, b, 'Cp_closure_error_pct_panels_only', '{:+.2f}')} %); correcting that pressure would raise the "
        f"lift above the one given. With the circulation reduced by the Prandtl-Glauert factor "
        f"({m(D, a, 'prandtl_glauert_factor', '{:.3f}')} and {m(D, b, 'prandtl_glauert_factor', '{:.3f}')}), the "
        f"corrected pressure returns the lift within "
        f"{m(D, a, 'Cp_closure_error_pct_prandtl_glauert_factor', '{:+.2f}')} % and "
        f"{m(D, b, 'Cp_closure_error_pct_prandtl_glauert_factor', '{:+.2f}')} %: the Kármán-Tsien rule raises a "
        f"strong suction by more than that factor. The reconstruction that is stored therefore solves the "
        f"fraction of the circulation at each instant so that the corrected pressure returns the lift given: "
        f"at peak lift the fraction is {m(D, a, 'circulation_factor_at_CLmax', '{:.3f}')} and "
        f"{m(D, b, 'circulation_factor_at_CLmax', '{:.3f}')} (over the cycle from "
        f"{m(D, a, 'circulation_factor_cycle_min', '{:.3f}')} to "
        f"{m(D, a, 'circulation_factor_cycle_max', '{:.3f}')} and from "
        f"{m(D, b, 'circulation_factor_cycle_min', '{:.3f}')} to "
        f"{m(D, b, 'circulation_factor_cycle_max', '{:.3f}')}), and what is left of the difference is the "
        f"tolerance of that solution ({m(D, a, 'Cp_closure_error_pct', 'pow10')} % and "
        f"{m(D, b, 'Cp_closure_error_pct', 'pow10')} %). The lift closes because it is made to; the corrected "
        f"pressure inside the region beyond the critical value, which is not physical, is part of the integral. "
        f"The reconstructed surface pressures are therefore a picture and not a second source of the loads.")
    R.p(f"<b>The reconstruction is given the lift and nothing else.</b> Its surface pressure does not return "
        f"the model's moment or drag. At peak lift the pitching moment about the quarter chord found from the "
        f"surface pressure differs from the model's by {m(D, a, 'Cp_moment_minus_model_at_CLmax', '{:+.3f}')} in "
        f"{mid.short_case(a)} and {m(D, b, 'Cp_moment_minus_model_at_CLmax', '{:+.3f}')} in {mid.short_case(b)}, "
        f"and the two have the same sign at only "
        f"{m(D, a, 'Cp_moment_sign_agrees_with_model_instants', '{:.0f}')} and "
        f"{m(D, b, 'Cp_moment_sign_agrees_with_model_instants', '{:.0f}')} of the "
        f"{m(D, a, 'stored_instants', '{:.0f}')} stored instants: through stall the model's moment is nose "
        f"down and that of the reconstruction is not. The pressure drag of the surface pressure differs from "
        f"the model's drag by {m(D, a, 'Cp_drag_minus_model_at_CLmax', '{:+.3f}')} and "
        f"{m(D, b, 'Cp_drag_minus_model_at_CLmax', '{:+.3f}')} at peak lift and by up to "
        f"{m(D, a, 'Cp_drag_minus_model_worst_abs_stored_instants', '{:.3f}')} and "
        f"{m(D, b, 'Cp_drag_minus_model_worst_abs_stored_instants', '{:.3f}')} over the stored instants. A "
        f"potential flow with a transpiration carries neither the aft loading of a separated flow nor its "
        f"drag. No pressure, moment or drag should be read from these fields.")
    panel_convergence(R, D)


def panel_convergence(R, D) -> None:
    """How far the quoted quantities of the reconstruction depend on the
    number of panels."""
    df = D.csv(PANELS)
    D.whole(PANELS, "all rows")
    base = D.note(int(df["n_panels_configured"].iloc[0]), PANELS, "n_panels_configured", "{:d}")
    ok = D.note(int(df["within_tolerance"].sum()), PANELS, "number of rows with within_tolerance = True", "{:d}")
    n = D.note(len(df), PANELS, "number of rows", "{:d}")
    labels = {"Cp_surface_min_incompressible": "Lowest incompressible surface pressure coefficient",
              "transpiration_vt_over_U": "Transpiration at the trailing edge, over U",
              "closure_error_pct_incompressible": "Lift difference of the incompressible surface pressure",
              "beyond_critical_surface_fraction": "Part of the surface beyond the critical pressure"}
    R.p(f"<b>Dependence on the number of panels.</b> The six stored instants of each case were solved again with "
        f"a quarter, a half and twice the {base} panels used. {R.next_tab()} gives the largest change of each "
        f"quoted quantity over the six instants when the panels are doubled, with the tolerance set for it: "
        f"{ok} of the {n} entries change by less than their tolerance on doubling from {base} panels. The "
        f"lowest pressure coefficient is the value at a control point of an incompressible solution with a "
        f"sharp suction peak; it is quoted to the figures this table supports and no further.")
    rows = [[mid.short_case(r.case).replace(" ", "\u00a0"), labels[r.quantity], r.change_measured_as,
             cell(r.largest_change_half_to_configured, 3), cell(r.largest_change_configured_to_double, 3),
             cell(r.tolerance, 3), "yes" if r.within_tolerance else "no"] for r in df.itertuples()]
    R.table(["Case", "Quantity", "Change measured as", "Half to the panels used", "Panels used to twice as many",
             "Tolerance", "Within tolerance on doubling"], rows,
            f"Largest change of the quoted quantities of the reconstruction over the six instants when the number "
            f"of panels is doubled ({PANELS}).", size=7.5)


def instants_table(R, D, name: str) -> None:
    """The reconstruction of one case, instant by instant, as three tables:
    the flow quantities, the corrected pressure with the region beyond the
    critical value, and the loads the surface pressure returns."""
    rel = reconstruction_file(name)
    for columns, what, extra in (
            (INSTANT_FLOW, "the transpiration and the circulations", " Total circulation is that of the sheet plus "
             "that of the vortex marker. A blank in the last column means the marker is not over the chord at "
             "that instant."),
            (INSTANT_STATE, "the corrected pressure, the region beyond the critical pressure and the local state",
             " The field is stored at four of the six instants; the others are blank in the field columns."),
            (INSTANT_LOADS, "the loads its surface pressure returns beside the model's", "")):
        df = D.csv(rel)[[c for c, _ in columns]]
        D.whole(rel, "columns " + ", ".join(c for c, _ in columns))
        _, rows = frame_rows(df, 3)
        R.table([label for _, label in columns], [[c.replace("_", " ") for c in r] for r in rows],
                f"{mid.short_case(name)}: the reconstruction instant by instant, {what} ({rel}).{extra}", size=7.5)


def field_numbers(R, D) -> None:
    """The scalar figures of the reconstruction and its pictures."""
    a, b = list(meta.CASES)
    R.h2("Figures and pictures of the reconstruction")
    R.p(f"{R.next_tab()} gives the scalar figures of the reconstruction, {R.next_tab(2)} the quantities that "
        f"are assumed and not results, and {R.next_tab(3)} to {R.next_tab(8)} the transpiration, the "
        f"circulations, the corrected pressure with the size of the region beyond the critical value, the local "
        f"state and the closure at each stored instant.")
    sa = mid.short_case(a)
    shown = "The separated region is not drawn; the triangle marks the separation point of the load model. "
    hatched = flowfield.CRITICAL_NOTE + " "
    figures = [
        (f"contour_Cp_{a}_dsv_a*.png", f"{sa}, instant of greatest vortex normal force: corrected pressure "
         f"coefficient. {hatched}{shown}{RECON}", {}),
        (f"contour_Mach_{a}_peak_a*.png", f"{sa}, instant of greatest lift: local Mach number from the corrected "
         f"pressure. {hatched}{shown}{RECON}", {}),
        (f"contour_Tstatic_{a}_peak_a*.png", f"{sa}, instant of greatest lift: static temperature, an isentropic "
         f"estimate from the corrected pressure. {hatched}{RECON}", {}),
        (f"contour_Trecovery_{a}_peak_a*.png", f"{sa}, instant of greatest lift: recovery temperature. It belongs "
         f"to an adiabatic wall under a turbulent boundary layer, which the reconstruction does not have; away "
         f"from the wall the map shows what such a wall would take at the local state. {hatched}{RECON}", {}),
        (f"temperature_profile_{a}.png", f"{sa}: static and recovery temperature along the surface; the curves "
         f"stop where the corrected pressure is below the critical value. {RECON}", {}),
        (f"contour_speed_stream_{a}_dsv_a*.png", f"{sa}, instant of greatest vortex normal force: speed of the "
         f"incompressible solution and streamlines. {shown}{RECON}", {}),
        (f"contour_speed_stream_{a}_fall_a*.png", f"{sa}, downstroke, deep in stall by the load model: speed and "
         f"streamlines. The flow drawn looks attached because the separated region is not represented. "
         f"{shown}{RECON}", {}),
        (f"contour_vorticity_{a}_dsv_a*.png", f"{sa}, instant of greatest vortex normal force: vorticity. A "
         f"picture of the singularities with no physical sign information: it is zero outside the core of the "
         f"vortex marker, and no boundary-layer or wake vorticity is in it. {RECON}", {}),
        (f"cp_distribution_{a}.png", f"{sa}: reconstructed surface pressure at instants of the cycle, with the "
         f"band beyond the critical value hatched. The note under the figure gives the part of the surface in "
         f"that band, the closure of the lift and the change of pressure beneath the vortex marker at each "
         f"instant with its sign. {RECON}", {}),
        (f"cp_phase_map_{a}.png", f"{sa}: reconstructed upper-surface pressure over the chord and the cycle. "
         f"{hatched}{RECON}", dict(width=0.9)),
        ("fig3d_cp_phase_surface.png", f"{sa}: the same as a surface, suction upwards, cut off at the critical "
         f"value. {RECON}", dict(width=0.85)),
        ("fig3d_field_surface_Cp.png", f"{sa}, instant of greatest lift: the corrected pressure field as a "
         f"surface, suction upwards, cut off at the critical value. {RECON}", dict(width=0.85)),
        ("fig3d_field_surface_speed.png", f"{sa}, instant of greatest lift: the speed of the incompressible "
         f"solution as a surface. {RECON}", dict(width=0.85)),
        ("fig3d_section_vectors.png", f"{sa}, instant of greatest lift: the section drawn with span, with the "
         f"reconstructed velocity on its near end plane. The span is for the picture only; the calculation is "
         f"two-dimensional. {RECON}", dict(width=0.85)),
        (f"contour_Cp_{b}_dsv_a*.png", f"{mid.short_case(b)}, instant of greatest vortex normal force: corrected "
         f"pressure coefficient. {hatched}{shown}{RECON}", {}),
        (f"contour_Mach_{b}_peak_a*.png", f"{mid.short_case(b)}, instant of greatest lift: local Mach number "
         f"from the corrected pressure. {hatched}{shown}{RECON}", {}),
    ]
    R.p(f"{R.next_fig()} to {R.next_fig(len(figures))} show the reconstruction; every one of them is subject to "
        "the limits above. In every figure the pressure coefficient is shown with suction upwards, and each "
        "figure of a corrected quantity prints the part of the surface and of the field beyond the critical "
        "pressure at its instant.")
    queue = [(plot(D, pattern), caption, kw) for pattern, caption, kw in figures]

    def place(n: int) -> None:
        """Place the next n figures of the queue, so that they come between the tables."""
        for path, caption, kw in [queue.pop(0) for _ in range(min(n, len(queue)))]:
            R.fig(path, caption, **kw)

    for name in (a, b):
        D.whole(metrics_file(name), "the metrics of the reconstruction")
    place(3)
    R.table(["Quantity"] + [mid.short_case(n) for n in (a, b)],
            [[label, m(D, a, key, fmt), m(D, b, key, fmt)] for key, label, fmt in FIELD_METRICS],
            f"Scalar figures of the reconstruction ({metrics_file(a)}, {metrics_file(b)}).")
    place(3)
    R.table(["Assumed quantity"] + [mid.short_case(n) for n in (a, b)],
            [[label, m(D, a, key, fmt), m(D, b, key, fmt)] for key, label, fmt in ASSUMED_METRICS],
            f"Quantities of the vortex marker that follow from constants chosen for illustration: settings, not "
            f"results, and the same in the two cases at the same clock ({metrics_file(a)}, {metrics_file(b)}).")
    place(3)
    instants_table(R, D, a)
    place(4)
    instants_table(R, D, b)
    place(len(queue))


def grid(R, D) -> None:
    """Section: the grid."""
    def q(metric: str, fmt: str) -> str:
        return D.c(MESH, "value", fmt, metric=metric)

    R.h1("The grid", "grid")
    R.p(f"The folder 02_mesh holds a structured O-grid round the section: {q('nodes_wrap_i', '{:d}')} nodes "
        f"round the section by {q('nodes_normal_j', '{:d}')} outwards, {q('total_cells', '{:d}')} cells, a "
        f"far boundary at {q('farfield_radius', '{:.0f}')} chords, a first layer "
        f"{q('first_layer_height', '{:.3f}')} of the chord high and a growth ratio of "
        f"{q('wall_normal_growth_ratio', '{:.3f}')}. Its status, as the file states it, is: "
        f"\"{q('status', '{}')}\" ({MESH}).")
    R.note("<b>What the grid is for.</b> The load model has no spatial mesh and does not read this grid; no "
           "load in this report depends on it. No equation of the flow is solved on it, and its wall spacing "
           "is not sized for a boundary layer. Its one use is as a set of points at which the reconstructed "
           "flow field is evaluated for one picture. The folders of this case study follow the layout of a "
           "flow-solver case, but the calculation is a set of ordinary differential equations marched in time; "
           "in a paper on the load model the grid and the field reconstruction would be left out.")
    R.p(f"The largest equiangle skewness is {q('max_skewness_equiangle', '{:.3f}')}, in the cells at the "
        f"sharp trailing edge; without those cells it is "
        f"{q('max_skewness_without_trailing_edge_cells', '{:.3f}')}. The smallest orthogonality angle is "
        f"{q('min_orthogonality_angle', '{:.1f}')}°, or "
        f"{q('min_orthogonality_angle_without_trailing_edge_cells', '{:.1f}')}° without the trailing-edge "
        f"cells, and the number of inverted cells is {q('inverted_cells', '{:d}')}.")
    R.fig(D.root/"02_mesh/fig_mesh_full.png", "The O-grid round the section (02_mesh). An illustration: no "
          "equation of the flow is solved on it, and the load model does not read it.")
    R.fig([D.root/"02_mesh/fig_mesh_le_zoom.png", D.root/"02_mesh/fig_mesh_te_zoom.png"],
          "The grid at the leading edge (left) and at the trailing edge (right) (02_mesh).")
    field_on_grid(R, D)


def field_on_grid(R, D) -> None:
    """The reconstructed field at the grid nodes."""
    drv = "field_reconstruction.read_by_the_driver"
    name = next(iter(meta.CASES))
    rel = f"05_solution/field_{name}_peak_on_mesh.csv.gz"
    df = D.csv(rel)
    nodes = D.note(len(df), rel, "number of rows", "{:d}")
    wall = D.note(int((df["source"] == "surface").sum()), rel, "number of rows with source = surface", "{:d}")
    none = D.note(int((df["source"] == "none").sum()), rel, "number of rows with source = none", "{:d}")
    R.p(f"The script 02_mesh/field_on_mesh.py evaluates the reconstructed field at the {nodes} nodes of the "
        f"O-grid, for each case at the instant of peak lift read from the stored time history, and writes it to "
        f"05_solution (field_&lt;case&gt;_peak_on_mesh.csv.gz). The {wall} nodes of the wall line take the values "
        f"on the surface itself; the {none} nodes nearer the wall than the distance within which the field "
        f"shows the junctions of the panels carry no value; the rest take the field. {R.next_fig()} shows the "
        f"pressure coefficient for {mid.short_case(name)}. It is the same reconstruction as in Section "
        f"{R.sec('field')}, with the same limits: the grid supplies the points and nothing else. The maps of "
        f"Section {R.sec('field')} are drawn from the reconstruction evaluated on its own rectangular grid of "
        f"{D.j(CFG, f'{drv}.grid_nx', '{:d}')} by {D.j(CFG, f'{drv}.grid_ny', '{:d}')} points ({CFG}).")
    R.fig(D.root/f"02_mesh/fig_mesh_field_{name}.png",
          f"{mid.short_case(name)}: pressure coefficient of the reconstructed field at the nodes of the "
          f"O-grid, at the instant of peak lift. {RECON}")


def pivot_rows(df, value: str, fmt: str) -> tuple:
    """A response-surface quantity as a table: mean incidence down, k across."""
    piv = df.pivot(index="alpha_mean_deg", columns="reduced_freq_k", values=value)
    header = ["Mean incidence, deg"] + [f"k = {k:.2f}" for k in piv.columns]
    rows = [[f"{a:.0f}"] + [fmt.format(v) for v in piv.loc[a]] for a in piv.index]
    return header, rows


def response(R, D) -> None:
    """Section: the response over mean incidence and reduced frequency."""
    df = D.csv(SURFACE)

    def rng(col: str, fmt: str) -> str:
        return (f"{D.note(df[col].min(), SURFACE, f'smallest {col}', fmt)} to "
                f"{D.note(df[col].max(), SURFACE, f'largest {col}', fmt)}")

    n_neg = D.note(int((df["cycle_damping_Xi"] < 0).sum()), SURFACE, "number of rows with cycle_damping_Xi < 0",
                   "{:d}")
    R.h1("The response over mean incidence and reduced frequency", "response")
    R.p(f"The load model is run at {D.note(len(df), SURFACE, 'number of rows', '{:d}')} conditions: mean "
        f"incidence from {rng('alpha_mean_deg', '{:.0f}')}° and reduced frequency from "
        f"{rng('reduced_freq_k', '{:.2f}')}, at the Mach number "
        f"({D.note(df['mach_M'].iloc[0], SURFACE, 'mach_M', '{:.2f}')}) and pitch amplitude "
        f"({D.note(df['alpha_amp_deg'].iloc[0], SURFACE, 'alpha_amp_deg', '{:.0f}')}°) of Case B. As in Case B, "
        f"this is a section in a steady stream with prescribed pitch, not a rotor calculation. The highest "
        f"peak incidence is {D.note(df['peak_alpha_deg'].max(), SURFACE, 'largest peak_alpha_deg', '{:.0f}')}°, "
        f"inside the {D.j(CFG, 'compared_range.peak_alpha_max_deg', '{:.0f}')}° reached by the measured loops "
        f"the model was compared with.")
    def dmp(col: str, fmt: str) -> str:
        return D.c(DAMPING, col, fmt, group=mid.ALL_LOOPS)

    sign = (f"on the {dmp('n_frames', '{:d}')} held-out loops with a closed measured moment loop the model gets "
            f"the sign of the damping right on a fraction {dmp('sign_agreement', '{:.2f}')}, where a prediction "
            f"of positive damping everywhere would score {dmp('sign_agreement_if_always_positive', '{:.2f}')}, "
            f"and {dmp('measured_negative', '{:d}')} of those loops measure negative damping ({DAMPING})")
    R.p(f"Over these conditions the peak lift runs from {rng('CL_max', '{:.2f}')}, the minimum moment from "
        f"{rng('CM_min', '{:.2f}')} and the onset incidence from {rng('stall_onset_alpha_deg', '{:.1f}')}°. The "
        f"cycle damping the model gives runs from {rng('cycle_damping_Xi', '{:.3f}')} and is negative at "
        f"{n_neg} of the conditions. That is not evidence that the damping is positive there: {sign}. The "
        f"damping surface is model output: it does not locate a stall-flutter boundary, and the sign of the "
        f"damping in it is not evidence of what a measurement would show. It is kept here only with that "
        f"figure beside it, which is also printed on the plot, and it would be left out of a paper. "
        f"{R.next_fig()} shows the surfaces and {R.next_tab()} and {R.next_tab(2)} the values ({SURFACE}).")
    R.fig(plot(D, "fig3d_response_surface.png"),
          f"Peak lift and cycle damping given by the load model over mean incidence and reduced frequency. The "
          f"damping is model output only: {sign}.", max_h_in=4.6)
    D.whole(SURFACE, "CL_max and cycle_damping_Xi against alpha_mean_deg and reduced_freq_k")
    R.table(*pivot_rows(df, "CL_max", "{:.3f}"), "Peak lift coefficient over mean incidence and reduced "
            f"frequency ({SURFACE}).")
    R.table(*pivot_rows(df, "cycle_damping_Xi", "{:.3f}"), "Cycle damping given by the model over mean "
            f"incidence and reduced frequency; model output only, with the sign agreement on the held-out loops "
            f"given in the text above ({SURFACE}).")


def drawings(R, D) -> None:
    """Section: the drawings."""
    dims = D.csv(DIMS)
    D.whole(DIMS, "all rows")
    rec = D.csv(SHEET_RECORD)
    D.whole(SHEET_RECORD, "all rows")
    sheets = sorted((D.root/SHEETS).glob("sheet*.png"))
    n = D.note(len(sheets), SHEETS, "number of files sheet*.png", "word")
    scales = rec[rec["item"] == "scale stated in the title block"].set_index("sheet")["value"]
    R.h1("Drawings", "drawings")
    R.p(f"The folder {SHEETS} holds {n} sheets. They are illustrative: they show the section of Case A, and "
        f"for Case B the rotor whose blade station gives the condition, so that a reader can see where the "
        f"dimensions come from. They are not for manufacture, and the helicopter outline in them is a generic "
        f"sketch drawn from one set of body coordinates on both sheets that show it. Every dimension on them "
        f"is taken from the one table reproduced as {R.next_tab()}, with its source. Each sheet states its "
        f"scale, or that it has none, in its title block. Section A-A of sheet 3 is drawn as seen looking the "
        f"way its cutting-plane arrows point "
        f"({D.c(SHEET_RECORD, 'value', '{}', item='section A-A: cutting-plane arrows point')}), with the leading "
        f"edge on the {D.c(SHEET_RECORD, 'value', '{}', item='section A-A: leading edge drawn on the')} "
        f"({SHEET_RECORD}).")
    R.note("<b>Sheets 1 and 2 show a helicopter, and this study contains no rotor calculation.</b> They are kept "
           "in the repository only to show where the chord, the radius and the rotor speed of Case B come from. "
           "In a paper they would be left out, with sheet 4 and the station sketch of sheet 3 kept.")
    header, rows = frame_rows(dims)
    R.table(["Dimension", "Value", "Unit", "Source"], rows, f"The dimensions of the drawings ({DIMS}).")
    for path in sheets:
        number = int(path.stem.split("_")[0].replace("sheet", ""))
        words = path.stem.split("_", 1)[1].replace("_", " ").replace("AA", "A-A").replace("3view", "three views")
        R.fig(path, f"Drawing sheet {number}: {words}. Scale as the sheet states it: {scales[number]}. Shown "
              "reduced from the sheet format. Illustrative; not for manufacture.", sheet=True)


def limits(R, D) -> None:
    """Section: limits of the study."""
    T = front.TARGETS

    def miss(key: str, fmt: str) -> str:
        return (f"{D.c(T, 'measured', fmt, measure=key)} against a target of {D.c(T, 'kind', '{}', measure=key)} "
                f"{D.c(T, 'target', fmt, measure=key)}")

    R.h1("Limits of the study", "limits")
    R.bullets([
        f"<b>Accuracy.</b> On the held-out loops the loop error in lift is {miss('mean_nRMS_CL', '{:.3f}')}, "
        f"the loop error in moment {miss('mean_nRMS_CM', '{:.3f}')}, and the fraction of loops with the sign "
        f"of cycle damping correct {miss('damping_sign_agreement', '{:.2f}')}. Those targets are not met.",
        "<b>Cycle damping.</b> The model does not predict cycle damping through stall and cannot locate a "
        "stall-flutter boundary.",
        "<b>The unsteady-moment factor</b> is an empirical correction of this model, taken from few loops at "
        "one Mach number. It is not a property of the flow.",
        "<b>The separation law.</b> On this data set the tabulated separation point gives the same accuracy as "
        "the usual fitted law; it is not shown to be better.",
        f"<b>Range.</b> The results are for {front.compared_range(D)}. Most of the measured loops are near the "
        "upper end of that Mach range. Nothing is claimed outside it, and nothing is claimed for another "
        "section.",
        "<b>Case A</b> is set beside a calibration loop, which is not an independent test.",
        "<b>Case B</b> is an aerofoil at a rotor-blade-station condition in a steady stream with prescribed "
        "pitch, held at the lowest Mach number the station sees round the azimuth, with an assumed incidence "
        "history. It is not a rotor calculation and no rotor measurement is compared with.",
        f"<b>The results of the cases carry the error of the model.</b> {front.band_note(D)} The results of "
        f"Case B are further uncertain within the spread between the two Mach stations of the static inputs "
        f"(Section {R.sec('static')}).",
        "<b>The flow fields</b> are a reconstruction, not a flow solution. " + RECON + " Its surface pressure "
        f"does not return the model's moment or drag. Section {R.sec('field')} gives its limits one by one and "
        "how far the surface pressure is from returning the lift it was given.",
        "<b>The grid</b> is not read by the load model, and no equation of the flow is solved on it.",
        "<b>The drawings</b> are illustrative; sheets 1 and 2 show a helicopter although the study contains no "
        "rotor calculation.",
        f"<b>The state-space form</b> agrees with the indicial march to the figures of Section "
        f"{R.sec('statespace')}; where the model switches a time constant the two forms differ by an amount "
        "that depends on the step. The growth or decay it gives for a section free in torsion is model output, with "
        "the same standing as the cycle damping.",
    ])


def references(R, D) -> None:
    """Section: references."""
    src = D.c(DIMS, "source", "{}", name="Case B rotor radius")
    R.h1("References", "references")
    for line in R.reference_lines():
        R.p(line)
    D.whole("docs/references.bib", "entries " + ", ".join(R.cited))
    R.p(f"The rotor dimensions of Case B are those the dimension table attributes to \"{src}\" ({DIMS}). That "
        "report is not an entry of docs/references.bib, and the four parameters taken from it (radius, chord, "
        "number of blades, rotor speed) have not been checked against its pages in this repository.")


def build(R, D) -> None:
    """Add the last part of the report."""
    for name in meta.CASES:
        case_results(R, D, name)
    field_description(R, D)
    field_limits(R, D)
    field_numbers(R, D)
    grid(R, D)
    response(R, D)
    drawings(R, D)
    limits(R, D)
    references(R, D)
