# -*- coding: utf-8 -*-
"""
07_report / rpt_text_front.py
-----------------------------
The first part of the report: title page, summary, the problem, what is
solved, and the account of the model with its equations and constants. Every
number in the text is read from a file through rpt_data.Data, which registers
it.

Author: Akosa Samuel Onyejekwe (independent)

Units: as in the files read; angles in degrees, lengths in metres.
"""
import re

import project_meta as meta
from rpt_blocks import escape
from unistall import flowfield

FLOW = "03_model_setup/flow_conditions.csv"
KIN = "03_model_setup/kinematics.csv"
CFG = "03_model_setup/solver_config.json"
SUMMARY = "05_solution/summary_all_cases.csv"
TARGETS = "results/validation_targets.csv"
BOARD = "results/targets_scoreboard.csv"
VTABLE = "results/validation_table.csv"
DIMS = "08_engineering_drawings/dimensions.csv"
GEOM = "01_geometry/section_geometry_summary.csv"
CASE_MD = "00_overview/case_definition.md"
FORM_MD = "docs/formulation.md"
TABLES_MD = "results/tables.md"
MOMENT = "results/attached_moment_factor.json"
CAL = "results/calibrated_constants.json"
BANDS = "05_solution/held_out_error_bands.csv"
STATION = "03_model_setup/station_condition.csv"

MEASURES = {
    "mean_nRMS_CL": ("Loop error, lift", "{:.3f}", ""),
    "mean_nRMS_CM": ("Loop error, moment", "{:.3f}", ""),
    "mean_nRMS_CD": ("Loop error, drag", "{:.3f}", ""),
    "mean_abs_dalpha_CLmax_deg": ("Incidence of maximum lift", "{:.2f}", "°"),
    "mean_abs_dalpha_Mstall_deg": ("Incidence of moment stall", "{:.2f}", "°"),
    "damping_sign_agreement": ("Fraction of loops with the sign of cycle damping correct", "{:.2f}", ""),
}


def clean_md(text: str) -> str:
    """A Markdown cell without its code ticks, emphasis marks and links."""
    text = re.sub(r"\[([^\]]+)\]\([^)]+\)", r"\1", text)
    return text.replace("`", "").replace("**", "")


def held_out(D) -> dict:
    """Count and range of the held-out loops inside the static Mach range."""
    df = D.csv(VTABLE)
    h = df[(df["model"] == "tabulated") & (df["set"] == "held_out") & df["in_mach_range"]]
    sel = "rows with model=tabulated, set=held_out, in_mach_range=True"
    return dict(n=D.note(len(h), VTABLE, f"number of {sel}", "{:d}"),
                m_lo=D.note(h["M"].min(), VTABLE, f"smallest M over {sel}", "{:.2f}"),
                m_hi=D.note(h["M"].max(), VTABLE, f"largest M over {sel}", "{:.2f}"))


def compared_range(D) -> str:
    """The range the model has been compared with measurement over, in words."""
    h = held_out(D)
    k_max = D.j(CFG, "compared_range.k_max", "{:.2f}")
    a_max = D.j(CFG, "compared_range.peak_alpha_max_deg", "{:.0f}")
    return (f"Mach {h['m_lo']} to {h['m_hi']}, reduced frequency up to {k_max} and peak incidence up to "
            f"{a_max}° on the NACA 0012 section")


def case_line(D, name: str) -> str:
    """One case's condition in words."""
    mach = D.c(FLOW, name, "{:.2f}", parameter="freestream_mach_M")
    k = D.c(KIN, "reduced_freq_k", "{:.4f}", case=name)
    mean = D.c(KIN, "alpha_mean_deg", "{:.0f}", case=name)
    amp = D.c(KIN, "alpha_amp_deg", "{:.0f}", case=name)
    return f"Mach {mach}, reduced frequency {k}, incidence {mean}° ± {amp}°"


def band(D, quantity: str, unit: str = "") -> str:
    """The held-out error of the model in one headline quantity, as text:
    mean and standard deviation of model minus measured over the held-out
    loops."""
    pct = "percent" in str(D.csv(BANDS).set_index("quantity").loc[quantity, "error_is"])
    fmt = "{:+.1f}" if pct or unit == "°" else "{:+.2f}"
    spread = fmt.replace("+", "")
    unit = " %" if pct else unit
    return (f"{D.c(BANDS, 'mean_error', fmt, quantity=quantity)}{unit} mean, "
            f"{D.c(BANDS, 'std_error', spread, quantity=quantity)}{unit} standard deviation")


def band_note(D) -> str:
    """A sentence on what the error bands are."""
    n = D.c(BANDS, "n_loops", "{:d}", quantity="CL_max")
    return (f"The figure in brackets after each load is the error of the model in that quantity on the {n} "
            f"held-out loops, model minus measured, as mean and standard deviation ({BANDS}); the loads are "
            f"rounded to match.")


def case_result(D, name: str) -> str:
    """One case's main results in a sentence, each with the held-out error
    of the model in that quantity."""
    label = meta.CASES[name]["label"]
    return (f"<b>{label}</b> ({case_line(D, name)}): lift reaches "
            f"{D.c(SUMMARY, 'CL_max', '{:.2f}', case=name)} [{band(D, 'CL_max')}] at an incidence of "
            f"{D.c(SUMMARY, 'alpha_at_CLmax_deg', '{:.1f}', case=name)}° [{band(D, 'alpha_at_CLmax_deg', '°')}], "
            f"the moment about the quarter chord falls to {D.c(SUMMARY, 'CM_min_c4', '{:.2f}', case=name)} "
            f"[{band(D, 'CM_min_c4')}], drag reaches {D.c(SUMMARY, 'CD_max', '{:.2f}', case=name)} "
            f"[{band(D, 'CD_max')}], and the onset criterion of the vortex is passed at "
            f"{D.c(SUMMARY, 'stall_onset_alpha_deg', '{:.1f}', case=name)}° (no measured counterpart; the "
            f"incidence of moment stall is in error by {band(D, 'moment_stall_alpha_deg', '°')}).")


def scoreboard_sentence(D) -> str:
    """How many targets are met and not met, by kind of target, in words."""
    sb = D.csv(BOARD)
    left = sb[~sb["met"].isin(["yes", "no"])]
    if len(left):
        raise ValueError(f"{BOARD}: entries in column met other than yes or no: {sorted(set(left['met']))}")
    words = {"physical": "on the model against measurement or theory", "numerical": "on the numerical method",
             "software": "on the software"}
    parts = []
    for kind, label in words.items():
        d = sb[sb["kind_of_target"] == kind]
        met = D.note(int((d["met"] == "yes").sum()), BOARD, f"rows with met=yes [kind_of_target={kind}]", "{:d}")
        n = D.note(len(d), BOARD, f"number of rows [kind_of_target={kind}]", "{:d}")
        parts.append(f"{met} of {n} {label}")
    missed = D.note(int((sb["met"] == "no").sum()), BOARD, "number of rows with met=no", "{:d}")
    return "the targets met are " + ", ".join(parts) + f", and {missed} are not met"


def title_page(R, D) -> None:
    """The cover."""
    R.title(title=meta.TITLE, subtitle=meta.SUBTITLE, solver=meta.SOLVER_LONG, method=meta.METHOD,
            author=meta.AUTHOR_FULL, date=meta.STUDY_DATE, kind_line="Case-study report")
    for key in ("TITLE", "SUBTITLE", "SOLVER_LONG", "METHOD", "AUTHOR_FULL", "STUDY_DATE"):
        D.whole("project_meta.py", key)


def summary(R, D) -> None:
    """Section: summary."""
    h = held_out(D)
    lift = D.c(TARGETS, "measured", "{:.3f}", measure="mean_nRMS_CL")
    lift_t = D.c(TARGETS, "target", "{:.3f}", measure="mean_nRMS_CL")
    mom = D.c(TARGETS, "measured", "{:.3f}", measure="mean_nRMS_CM")
    mom_t = D.c(TARGETS, "target", "{:.3f}", measure="mean_nRMS_CM")
    sign = D.c(TARGETS, "measured", "{:.2f}", measure="damping_sign_agreement")
    sign_t = D.c(TARGETS, "target", "{:.2f}", measure="damping_sign_agreement")
    n_cases = D.note(len(meta.CASES), "project_meta.py", "number of entries of CASES", "word")
    R.h1("Summary")
    R.p(f"This report gives the unsteady lift, drag and pitching moment on a NACA 0012 section in a steady "
        f"stream with prescribed sinusoidal pitch, as predicted by {meta.SOLVER_LONG}. The load model is the "
        f"{meta.METHOD}. It is run at {n_cases} conditions.")
    R.bullets([case_result(D, name) for name in meta.CASES])
    R.p(band_note(D))
    R.p(f"How well the load model agrees with measurement is established separately, on {h['n']} held-out "
        f"oscillating-aerofoil loops of the NASA TM-84245 experiment {R.cite('mccroskey1982v1')}, none of which "
        f"is used in a fit. Of the targets set for the model, {scoreboard_sentence(D)}.")
    R.bullets([
        f"<b>Targets not met.</b> The loop error in lift is {lift} against a target of at most {lift_t}; the "
        f"loop error in moment is {mom} against at most {mom_t}; the sign of cycle damping is correct on a "
        f"fraction {sign} of the loops against a target of at least {sign_t}.",
        "<b>Cycle damping.</b> The model does not predict cycle damping through stall. It cannot locate a "
        "stall-flutter boundary and is not offered for that purpose. The damping figures of the two cases are "
        "reported as model output only.",
        "<b>The flow fields are a reconstruction, not a flow solution.</b> " + flowfield.STANDING_NOTE
        + " Its surface pressure is given the lift only and does not return the model's moment or drag.",
        "<b>Case A is set beside a calibration loop,</b> which is not an independent test.",
        "<b>Case B is not a rotor calculation.</b> It is an aerofoil at a rotor-blade-station condition in a "
        "steady stream with prescribed pitch.",
        "<b>The unsteady-moment factor</b> of the model is an empirical correction of this model, not a "
        "property of the flow.",
        f"<b>Range.</b> The results are for {compared_range(D)}. Nothing is claimed outside it.",
    ])


def problem(R, D) -> None:
    """Section: the problem."""
    R.h1("The problem")
    R.p("On a helicopter in fast forward flight the retreating blade is pitched up through stall once per "
        "revolution, at low dynamic pressure. A vortex forms near the leading edge, travels over the upper "
        "surface and is shed. While it is over the chord the lift exceeds its static maximum; as it leaves, "
        f"the lift collapses and a large nose-down pitching moment follows {R.cite('mccroskey1982arfm', 'carr1988')}"
        ". That cycle limits the thrust and speed of the rotor and drives vibration and control loads.")
    R.p("A design calculation needs the unsteady lift, drag and moment of the blade section many times over, "
        "so a load model that costs a small fraction of a flow solution is wanted. The Leishman-Beddoes "
        f"method {R.cite('leishman1989')} is such a model: a set of first-order lags, marched in time, that "
        "stand for the delay of the pressure, of the boundary layer and of the vortex. This report applies one "
        "implementation of it and states how far its predictions agree with measurement.")


def cases_table(R, D) -> str:
    """The table of the two cases."""
    names = list(meta.CASES)
    rows = [["Purpose"] + [D.c(FLOW, n, "{}", parameter="description") for n in names]]
    for label, par, fmt in (("Mach number", "freestream_mach_M", "{:.2f}"), ("Chord, m", "chord_c", "{:.3f}"),
                            ("Stream speed, m/s", "freestream_velocity_U", "{:.2f}"),
                            ("Reynolds number on chord in sea-level standard air (not the tunnel value of the "
                             "frame, nor a flight value; the load model does not use it)", "reynolds_number_Re_c",
                             "pow10")):
        rows.append([label] + [D.c(FLOW, n, fmt, parameter=par) for n in names])
    for label, col, fmt in (("Reduced frequency k", "reduced_freq_k", "{:.4f}"),
                            ("Mean incidence, deg", "alpha_mean_deg", "{:.1f}"),
                            ("Pitch amplitude, deg", "alpha_amp_deg", "{:.1f}"),
                            ("Pitch frequency, Hz", "freq_Hz", "{:.3f}"), ("Period, s", "period_s", "{:.4f}")):
        rows.append([label] + [D.c(KIN, col, fmt, case=n) for n in names])
    for label, key, fmt in (("Steps per cycle", "steps_per_cycle", "{:d}"), ("Cycles marched", "cycles", "{:d}")):
        rows.append([label] + [D.j(CFG, f"cases.{n}.march.{key}", fmt) for n in names])
    return R.table(["Quantity"] + [meta.CASES[n]["label"] for n in names], rows,
                   f"The cases ({FLOW}, {KIN}, {CFG}).")


def what_is_solved(R, D) -> None:
    """Section: what is solved."""
    R.h1("What is solved", "solved")
    R.p("The unsteady loads on a NACA 0012 section in a steady stream, pitching sinusoidally about the "
        "quarter chord, by the load model in unistall/dsmodel.py. The motion and the reduced frequency are")
    R.eq("motion")
    R.p(f"with U the stream speed and c the chord. {R.next_tab()} gives the conditions. Both lie inside "
        f"the range over which the model has been compared with measurement: {compared_range(D)}.")
    cases_table(R, D)
    frame = meta.CASES[meta.CASE_A]["measured_frame"].split("_")[1]
    D.note(frame, "project_meta.py", f"CASES.{meta.CASE_A}.measured_frame", "{}")
    R.p(f"<b>Case A</b> is a condition of the NASA TM-84245 experiment, so its loops can be set beside a "
        f"measured loop (frame {frame}). That frame is one of the calibration loops, so the comparison in "
        "Case A is an illustration and not an independent test; the independent test is the held-out set of "
        f"Section {R.sec('accuracy')}.")

    def dim(name: str) -> str:
        return D.c(DIMS, "value", "{}", name=name)

    def st(quantity: str, fmt: str) -> str:
        return D.c(STATION, "value", fmt, quantity=quantity)

    R.note(f"<b>Case B is not a rotor calculation.</b> Its chord and its once-per-revolution pitch frequency are "
           f"those of the UH-60A main rotor ({dim('Case B number of blades')} blades, radius "
           f"{dim('Case B rotor radius, in feet')} ft, chord {dim('Case B blade chord, in inches')} in, "
           f"{dim('Case B rotor speed')} rpm; {D.c(DIMS, 'source', '{}', name='Case B rotor radius')}). The "
           f"station at {dim('Case B analysis station, fraction of rotor radius')} of the radius turns at "
           f"{st('station_rotational_speed', '{:.1f}')} m/s. Case B takes the Mach number "
           f"{st('station_mach_M', '{:.2f}')} for it, which the station sees at azimuth "
           f"{st('station_azimuth', '{:.0f}')}° (the retreating side) at a forward speed of "
           f"{st('forward_speed_implied', '{:.1f}')} m/s ({st('forward_speed_implied_kt', '{:.0f}')} kt, advance "
           f"ratio {st('advance_ratio_implied', '{:.2f}')}); the pitch frequency is the rotor speed, "
           f"{st('pitch_frequency', '{:.2f}')} Hz, and the reduced frequency "
           f"{st('reduced_frequency_k', '{:.4f}')} follows from it, the chord and the stream speed ({STATION}). "
           f"Round the azimuth the Mach number at that station runs from "
           f"{st('mach_at_station_lowest', '{:.2f}')} to {st('mach_at_station_highest', '{:.2f}')}: Case B holds "
           f"the lowest value steady. Its incidence history is assumed and comes from no trim calculation. The "
           f"real blade uses the SC1095 and SC1094 R8 sections with a swept tip, not a NACA 0012. There is no "
           f"time-varying velocity, sweep, inflow, trim or blade motion, and no rotor measurement is compared "
           f"with.")
    t = D.c(GEOM, "value", "{:.4f}", property="max_thickness_t_c")
    xt = D.c(GEOM, "value", "{:.3f}", property="x_at_max_thickness_x_c")
    r = D.c(GEOM, "value", "{:.4f}", property="LE_radius_r_c")
    gap = D.c(GEOM, "value", "{:.1f}", property="TE_gap_percent_chord")
    R.p(f"The section is drawn in {R.next_fig()}. Its coordinates come from the NACA four-digit thickness form "
        f"with the closed trailing edge (a gap of {gap} % of the chord): the greatest thickness is {t} of the "
        f"chord at {xt} of the chord from the leading edge, and the leading-edge radius is {r} of the chord "
        f"({GEOM}). The same outline is used by the grid, by the panels of the flow-field reconstruction and by "
        f"the drawings.")
    R.fig(D.root/"01_geometry/fig_geometry_profile.png", "The NACA 0012 section (01_geometry).", width=0.85)
    parts = D.md(CASE_MD)["What each part of the case study is, and is not"]
    D.whole(CASE_MD, "table 'What each part of the case study is, and is not'")
    R.p(f"The case study is kept in numbered folders. {R.next_tab()} says what each holds and what standing "
        "it has.")
    R.table(list(parts.columns), [[clean_md(c) for c in row] for row in parts.itertuples(index=False)],
            f"The parts of the case study and their standing ({CASE_MD}).")


def model_attached(R, D) -> None:
    """Subsection: attached flow."""
    lo, hi = [s.strip() for s in fixed_constant(D, "blend limits in (18a)").split(",")]
    R.h2("Attached flow")
    R.p("Time is counted in semichords of travel, and compressibility enters through the Prandtl-Glauert "
        "factor:")
    R.eq("step")
    R.p("The pitch rate is made dimensionless with the chord and the stream speed, and the incidence is taken "
        "at the three-quarter chord, where x<sub>p</sub> is the pitch axis as a fraction of the chord:")
    R.eq("rate")
    R.p("The circulatory lift lags the incidence. The lag is carried by a pair of deficiency functions, each "
        "advanced by an exact exponential recurrence over the step:")
    R.eq("deficiency")
    R.p("They give the effective incidence and the circulatory normal force:")
    R.eq("circulatory")
    R.p("The impulsive (non-circulatory) force, the pressure wave of a sudden change of incidence, decays with "
        "its own time constant T<sub>α</sub>. With K<sub>α</sub> the rate of change of incidence, a second "
        "recurrence of the same form gives it, and a like term C<sub>N</sub><sup>nc,q</sup> follows from the "
        "pitch rate. The sum is the normal force the section would carry with the flow attached:")
    R.eq("impulsive")
    R.p(f"The indicial constants are those of compressible flow. Below Mach {hi} the attached-flow loads are "
        f"blended with incompressible loads marched by the same recurrences, with a weight w that falls from "
        f"one at Mach {hi} to zero at Mach {lo}; at the Mach numbers of this study w is one and the equations "
        f"above are unchanged. In attached flow the model is checked against Theodorsen's solution "
        f"{R.cite('theodorsen1935')} (Section {R.sec('accuracy')}).")


def model_separated(R, D) -> None:
    """Subsection: separation, the Kirchhoff relation and the two lags."""
    f_min = fixed_constant(D, "F_MIN")
    R.h2("Separated flow: the pressure lag, the boundary-layer lag and the Kirchhoff relation")
    R.p("The leading-edge pressure does not follow the normal force at once. A first-order lag with time "
        "constant T<sub>p</sub> gives a delayed normal force C'<sub>N</sub>:")
    R.eq("pressure_lag")
    R.p("The delayed force is turned into an equivalent incidence, at which the static position of the "
        "separation point f (the fraction of the chord over which the flow is attached) is read:")
    R.eq("lookup")
    R.p("Here the model departs from the usual form. The Leishman-Beddoes method fits f with a "
        "exponential law; this model reads f from a table made from the static measurements "
        f"themselves, by inverting the Kirchhoff relation at each measured incidence, with f kept between "
        f"{f_min} and one:")
    R.eq("static_f")
    R.p("The boundary layer in turn does not follow the pressure at once. A further lag with time constant "
        "T<sub>f</sub> gives the dynamic separation point f''. The multiplier σ<sub>1</sub> changes the time "
        "constant with the state of the flow (separating or reattaching, with or without a vortex over the "
        "chord), by fixed rules taken from the source and not fitted:")
    R.eq("boundary_layer_lag")
    R.p("The Kirchhoff relation then reduces the circulatory normal force for the separated flow:")
    R.eq("kirchhoff")


def model_vortex(R, D) -> None:
    """Subsection: the vortex, the chord force and the moment."""
    R.h2("The vortex lift")
    R.p("The vortex starts when the delayed normal force first exceeds the normal force at static stall, "
        "C<sub>N1</sub>. From then a clock τ<sub>v</sub> counts the semichords travelled. While the vortex "
        "is over the chord (τ<sub>v</sub> up to T<sub>VL</sub>) it is fed by the difference between the "
        "attached and the separated circulatory force, and its lift decays with time constant T<sub>V</sub>:")
    R.eq("vortex_feed")
    R.eq("vortex_lift")
    R.p("Once the vortex has passed the trailing edge the feed stops and its lift only decays. The multiplier "
        "σ<sub>3</sub> again follows fixed rules.")
    R.h2("The chord force, lift and drag")
    R.p("The chord force is the leading-edge suction, reduced by separation and by the recovery factor η:")
    R.eq("chord_force")
    R.p("The normal and chord forces are resolved into lift and drag, with the drag at zero lift "
        "C<sub>D0</sub> added:")
    R.eq("lift_drag")
    R.h2("The moment")
    g_m = D.j(MOMENT, "cm_unsteady", "{:.3f}")
    n_loops = D.note(len(D.json(MOMENT)["frames"]), MOMENT, "number of entries of frames", "word")
    R.p("The vortex lift acts at a point that moves aft as the vortex travels:")
    R.eq("vortex_arm")
    R.p("The moment about the quarter chord is the static moment read at a lagged incidence, the unsteady "
        "moment terms of attached flow, and the moment of the vortex lift:")
    R.eq("moment")
    rng = D.json(MOMENT)["range_over_loops_with_a_closed_moment_loop"]
    g_lo = D.note(rng[0], MOMENT, "range_over_loops_with_a_closed_moment_loop, first", "{:.2f}")
    g_hi = D.note(rng[1], MOMENT, "range_over_loops_with_a_closed_moment_loop, second", "{:.2f}")
    R.note(f"<b>The factor g<sub>M</sub> is empirical.</b> Below static stall the model gives a wider moment "
           f"loop than is measured. g<sub>M</sub> = {g_m} is the one number that minimises the moment loop "
           f"error on the {n_loops} calibration loops whose peak incidence is below static stall ({MOMENT}). "
           f"Part of the excess comes from reading the static moment at the lagged incidence and not from "
           f"the unsteady terms, so the factor corrects the model as a whole below stall; the file gives "
           f"both parts for each loop. Taken alone, the loops with a closed moment loop would choose {g_lo} "
           f"and {g_hi}. It rests on few loops at one Mach number and is not a property of the flow.")


GREEK = (("alpha", "α"), ("eta", "η"))
FILE_EQUATION = re.compile(r"(?<![\w.])\((\d+[a-z]?)\)")
SUBSCRIPT = re.compile(r"([A-Za-zα-ω]+)_([A-Za-z0-9]+(?:,[A-Za-z0-9]+)*)")


def typeset(text: str) -> str:
    """A cell of the constants table of docs/formulation.md for the report:
    symbols with their subscripts set as subscripts, Greek names as letters,
    names of the code left as they are, and each bare equation number marked
    as one of that file (the report numbers its own equations differently)."""
    out = []
    for i, part in enumerate(re.sub(r"\[([^\]]+)\]\([^)]+\)", r"\1", text).replace("**", "").split("`")):
        part = escape(part)
        if i % 2 == 0:
            for name, letter in GREEK:
                part = re.sub(rf"\b{name}(?=_|\b)", letter, part)
            part = SUBSCRIPT.sub(r"\1<sub>\2</sub>", part)
            part = FILE_EQUATION.sub(r"(\1) of formulation.md", part)
        out.append(part)
    return "".join(out)


def fixed_constant(D, name: str) -> str:
    """The value of a fixed constant as printed in section 5 of the formulation."""
    table = D.md(FORM_MD)["5. Constants"]
    row = table[table["Constant"] == name]
    return D.note(row["Value"].iloc[0], FORM_MD, f"section 5, constant '{name}', column Value", "{}")


def model_constants(R, D) -> None:
    """Subsection: the step rule and the tables of constants."""
    cfg_step = D.j(CFG, "march_rule.largest_step_semichords", "{:.2f}")
    steps_min = D.j(CFG, "march_rule.steps_min", "{:d}")
    settle = D.j(CFG, "march_rule.settle_semichords", "{:.0f}")
    c_lo, c_hi = D.j(CFG, "march_rule.cycles_min", "{:d}"), D.j(CFG, "march_rule.cycles_max", "{:d}")
    if "2 pi" not in R.equations["step_rule"]["source_text"]["54"]:
        raise ValueError(f"{FORM_MD}: line (54) no longer reads as the step rule typeset in the report")
    R.h2("The step rule")
    R.p("The model has no spatial mesh: its only discretisation is the time step. Because the lags are "
        "measured in semichords of travel, the step is set in semichords and not as a fraction of the cycle. "
        "The number of steps per cycle is the smallest that satisfies")
    R.eq("step_rule")
    R.p(f"with Δs<sub>max</sub> = {cfg_step}, and is never fewer than {steps_min}. Enough cycles are marched "
        f"for {settle} semichords of travel to precede the reported cycle, between "
        f"{c_lo} and {c_hi} cycles in all ({CFG}). Section {R.sec('numerics')} gives the resolution of the cases and what the "
        "convergence studies show.")
    R.h2("Constants")
    cal = D.md(TABLES_MD)["Table 3. Constants"]
    D.whole(TABLES_MD, "Table 3. Constants, columns Constant, Literature, Tabulated model, Status, Range")
    cost = D.j(CAL, "cost", "{:.3f}")
    cost_lit = D.j(CAL, "cost_at_literature_values", "{:.3f}")
    n_cal = D.j(CAL, "n_frames", "{:d}")
    n_fit = D.note(len(D.json(CAL)["constants"]), CAL, "number of entries of constants", "word")
    R.p(f"{R.next_tab()} lists the constants that are fitted or empirical, with their literature values. The "
        f"{n_fit} stall constants are fitted on {n_cal} calibration loops by least squares; the cost, the sum "
        f"over those loops of the squared loop errors in lift, moment and drag, falls from {cost_lit} at the "
        f"literature values to {cost} ({CAL}). {R.next_tab(2)} lists the constants that are fixed.")
    keep = ["Constant", "Literature", "Tabulated model", "Status", "Range over starts"]
    R.table(["Constant", "Literature value", "Value used", "Status", "Range over the starts of the fit"],
            [[typeset(c) for c in row] for row in cal[keep].itertuples(index=False)],
            f"Fitted and empirical constants; the time constants are in semichords ({TABLES_MD}, Table 3).",
            rich=True)
    fixed = D.md(FORM_MD)["5. Constants"]
    fixed = fixed[fixed["Status"] == "fixed"]
    D.whole(FORM_MD, "section 5, rows with Status = fixed")
    R.table(["Constant", "Value", "Source"],
            [[typeset(c) for c in row] for row in fixed[["Constant", "Value", "Source"]].itertuples(index=False)],
            f"Fixed constants ({FORM_MD}, section 5; [DH19] is Damiani and Hayman). Equation numbers in this "
            f"table are those of that file, not of this report; names in brackets are those of the code.",
            rich=True)


def model(R, D) -> None:
    """Section: the model."""
    R.h1("The model", "model")
    R.p(f"The load model is a chain of first-order lags marched in time. Attached-flow loads come from "
        f"indicial functions; a lagged normal force decides where the flow separates; a lagged separation "
        f"point reduces the circulatory load through the Kirchhoff relation; and a vortex, started when a "
        f"critical normal force is passed, adds lift and nose-down moment while it travels over the chord. "
        f"The equations are those of Leishman and Beddoes {R.cite('leishman1989')} in the form given by "
        f"Damiani and Hayman {R.cite('damiani2019')}. Every equation and constant the code evaluates is listed, "
        f"with the places where the implementation departs from that source, in {FORM_MD}; the principal ones "
        "follow. A subscript n marks the value at step n.")
    model_attached(R, D)
    model_separated(R, D)
    model_vortex(R, D)
    model_constants(R, D)


def build(R, D) -> None:
    """Add the first part of the report."""
    title_page(R, D)
    summary(R, D)
    problem(R, D)
    what_is_solved(R, D)
    model(R, D)
