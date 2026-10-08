# -*- coding: utf-8 -*-
"""
07_report / rpt_text_mid.py
---------------------------
The middle part of the report: static inputs, calibration and the comparison
with measurement, and numerical resolution. Every number in the text is read
from a file through rpt_data.Data, which registers it.

Author: Akosa Samuel Onyejekwe (independent)

Units: as in the files read; angles in degrees, travel in semichords.
"""
import re

import project_meta as meta
from rpt_data import cell, frame_rows, power_of_ten, sibling

front = sibling("rpt_text_front")

CFG = front.CFG
TARGETS = front.TARGETS
BOARD = front.BOARD
TABLES_MD = front.TABLES_MD
CAL = front.CAL
STATIC = "results/static_by_mach_naca0012.csv"
ONSET = "results/onset_level_check.csv"
SS = "results/statespace_summary.csv"
STRUCT = "results/structural_checks.csv"
CHOICE = "results/static_curve_choice.csv"
FLOOR = "results/uncertainty_summary.csv"
DEPTH = "results/validation_by_stall_depth.csv"
SAME = "results/comparison_same_constants.csv"
EQUIV = "results/equivalence.csv"
OWN = "results/comparison.csv"
DAMP = "results/validation_damping.csv"
CONV = "results/convergence_summary.csv"
REFINE = "05_solution/convergence/timestep_refinement.csv"
ORDER = "05_solution/convergence/timestep_order.csv"
SPREAD = "05_solution/static_station_spread.csv"
AIR = "03_model_setup/air_properties.csv"
STEP_TOLERANCE = 1e-4        # relative; the configuration gives the step to six decimals
ORDER_LABELS = {"CL_max": "Peak lift", "CM_min": "Minimum moment", "CD_max": "Peak drag",
                "stall_onset_alpha_deg": "Onset incidence, deg"}
FIGS = "results/figures"
PLOTS = "06_postprocessing/plots"

ALL_LOOPS = "all primary held-out loops with a closed moment loop"
OWN_LABELS = {
    "mean |CLmax_err_pct|": ("peak-lift error", "{:.2f}", " %"),
    "mean |dalpha_CLmax_deg|": ("incidence of maximum lift", "{:.2f}", "°"),
    "mean |dalpha_Mstall_deg|": ("incidence of moment stall", "{:.2f}", "°"),
    "mean |nRMS_CL|": ("loop error in lift", "{:.3f}", ""),
    "mean |nRMS_CM|": ("loop error in moment", "{:.3f}", ""),
    "mean |nRMS_CD|": ("loop error in drag", "{:.3f}", ""),
}


def short_case(name: str) -> str:
    """'Case A' or 'Case B' for a case key."""
    return meta.short_label(name)


def static_inputs(R, D) -> None:
    """Section: static inputs."""
    st = D.csv(STATIC)
    n_st = D.note(len(st), STATIC, "number of rows", "word")
    machs = " and ".join(D.note(m, STATIC, "M", "{:.3f}") for m in st["M"])
    more = D.c(CHOICE, "cost", "{:.2f}", post_stall_curve="more separated")
    less = D.c(CHOICE, "cost", "{:.2f}", post_stall_curve="less separated")
    R.h1("Static inputs", "static")
    R.p(f"The model takes from static measurement the lift-curve slope, the zero-lift incidence, the normal "
        f"force at static stall C<sub>N1</sub>, the separation point f and the static moment, each as a "
        f"function of Mach number. They are tabulated at {n_st} Mach stations, {machs}, from quasi-steady "
        f"sweeps of the same experiment and from the static figures of the NASA report "
        f"{R.cite('mccroskey1982v1', 'mcalister1982v2')}, and interpolated between the stations "
        f"({R.next_tab()}).")
    cols = ["M", "CL_alpha_per_deg", "alpha_zero_lift_deg", "CL_max", "alpha_at_CL_max_deg", "CN_at_CL_max", "frame"]
    D.whole(STATIC, "columns " + ", ".join(cols))
    _, rows = frame_rows(st[cols])
    R.table(["Mach", "Lift-curve slope, per deg", "Zero-lift incidence, deg", "Maximum static lift coefficient",
             "Static-stall incidence, deg", "Normal force at static stall", "Source sweep"], rows,
            f"Static inputs at the Mach stations ({STATIC}).")
    rows = []
    for label, key, fmt in (("Lift-curve slope, per rad", "CN_alpha_per_rad", "{:.3f}"),
                            ("Zero-lift incidence, deg", "alpha0_deg", "{:.2f}"),
                            ("Static normal force at zero incidence", "CN_static_at_zero_incidence", "{:+.3f}"),
                            ("Static-stall incidence, deg", "alpha_stall_deg", "{:.2f}"),
                            ("Normal force at onset", "CN1_onset", "{:.3f}")):
        rows.append([label] + [D.j(CFG, f"cases.{n}.static_inputs.{key}", fmt) for n in meta.CASES])
    R.table(["Quantity"] + [meta.CASES[n]["label"] for n in meta.CASES], rows,
            f"Static inputs interpolated to the Mach numbers of the cases ({CFG}).")
    a, b = list(meta.CASES)
    R.p(f"<b>The static inputs carry the offset of the measured sweeps.</b> The section is symmetric, yet the "
        f"static inputs give it a zero-lift incidence of "
        f"{D.j(CFG, f'cases.{a}.static_inputs.alpha0_deg', '{:+.2f}')}° in {short_case(a)} and "
        f"{D.j(CFG, f'cases.{b}.static_inputs.alpha0_deg', '{:+.2f}')}° in {short_case(b)}, and a static normal "
        f"force at zero incidence of {D.j(CFG, f'cases.{a}.static_inputs.CN_static_at_zero_incidence', '{:+.3f}')} "
        f"and {D.j(CFG, f'cases.{b}.static_inputs.CN_static_at_zero_incidence', '{:+.3f}')}. The offset is that of "
        f"the measured static sweeps (a flow angle of the tunnel or an offset of the balance), not a property "
        f"of the section; the model keeps it so that it reads the static data as measured. The incidence "
        f"uncertainty the experiment states is not held in the data files of this repository, so the offset is "
        f"not set against it here.")
    station_spread(R, D, b)
    R.p(f"Past static stall the measurements are not single-valued: they show a more separated and a less "
        f"separated state. The model follows the attached curve up to static stall and the more separated "
        f"state beyond it. That choice was made on the calibration loops, where it gives a cost of {more} "
        f"against {less} for the less separated state ({CHOICE}). {R.next_fig()} shows the static normal "
        f"force, separation point and moment at the Mach numbers of the cases, and {R.next_fig(2)} the same "
        "inputs at the Mach stations with the measurements they come from.")
    R.fig(D.root/PLOTS/"static_inputs.png",
          "Static inputs of the load model at the Mach numbers of the cases (06_postprocessing).")
    R.fig(D.root/FIGS/"fig01_static_inputs.png",
          "Static normal force, separation point and moment at the Mach stations (results/figures).")


def station_spread(R, D, name: str) -> None:
    """What the blend of the static inputs between the Mach stations is
    worth in the results of a case that lies between them."""
    df = D.csv(SPREAD)
    lo, hi = (D.note(m, SPREAD, "static_mach of a station row", "{:.3f}")
              for m in df["static_mach"].dropna().iloc[1:3])

    def v(col: str, fmt: str, row: str) -> str:
        return D.c(SPREAD, col, fmt, static_inputs=row)

    blend, gap = "interpolated to the case's Mach number (reported)", "spread between the stations"
    R.note(f"<b>The static inputs of {short_case(name)} are a straight-line blend of two Mach stations, and no "
           f"static measurement at its Mach number supports the blend.</b> Run with the inputs of the station at "
           f"Mach {lo} and of the station at Mach {hi} in turn, {short_case(name)} gives a peak lift of "
           f"{v('CL_max', '{:.2f}', 'Mach station ' + lo)} and {v('CL_max', '{:.2f}', 'Mach station ' + hi)} "
           f"(reported {v('CL_max', '{:.2f}', blend)}), a minimum moment of "
           f"{v('CM_min', '{:.2f}', 'Mach station ' + lo)} and {v('CM_min', '{:.2f}', 'Mach station ' + hi)} "
           f"(reported {v('CM_min', '{:.2f}', blend)}) and an onset incidence of "
           f"{v('stall_onset_alpha_deg', '{:.1f}', 'Mach station ' + lo)}° and "
           f"{v('stall_onset_alpha_deg', '{:.1f}', 'Mach station ' + hi)}° (reported "
           f"{v('stall_onset_alpha_deg', '{:.1f}', blend)}°). The spread between the stations is "
           f"{v('CL_max', '{:.2f}', gap)} in peak lift, {v('CM_min', '{:.3f}', gap)} in minimum moment, "
           f"{v('CD_max', '{:.2f}', gap)} in peak drag and {v('stall_onset_alpha_deg', '{:.1f}', gap)}° in onset "
           f"incidence ({SPREAD}); the results of {short_case(name)} are uncertain within it.")


def calibration(R, D) -> None:
    """Subsection: the measured loops and the fit."""
    loops = D.md(TABLES_MD)["Table 1. Measured loops used"]
    D.whole(TABLES_MD, "Table 1. Measured loops used")
    n_cal = D.j(CAL, "n_frames", "{:d}")
    starts = D.j(CAL, "starts", "{:d}")
    cost, lit = D.json(CAL)["cost"], D.json(CAL)["cost_at_literature_values"]
    drop = D.note(100.0*(lit - cost)/lit, CAL, "100 (cost_at_literature_values - cost) / cost_at_literature_values",
                  "{:.0f}")
    costs = sorted(D.json(CAL).get("first_stage_record", D.json(CAL))["polished_costs"])
    polished = (f"from {D.note(costs[0], CAL, 'lowest of polished_costs', '{:.3f}')} to "
                f"{D.note(costs[-1], CAL, 'highest of polished_costs', '{:.3f}')}")
    R.h2("The measured loops and the fit")
    R.p(f"All measured loops are from the NASA TM-84245 experiment on oscillating aerofoils "
        f"{R.cite('mccroskey1982v1')}. {R.next_tab()} lists how they are divided. The loops outside the Mach "
        "range of the static data are scored but support no claim.")
    R.table(list(loops.columns), [list(r) for r in loops.itertuples(index=False)],
            f"Measured loops used ({TABLES_MD}, Table 1).")
    R.p(f"The stall constants are fitted on the {n_cal} calibration loops inside the static Mach range, by "
        f"bounded least squares from {starts} starting points, each refined at the resolution the "
        f"results are reported at. The fit lowers the cost by {drop} % from the literature constants, and the "
        f"starts refined at full resolution reach costs of {polished} ({CAL}). {R.next_fig()} shows the cost "
        "of every start, the sensitivity of the cost to each constant, the cross-validation inside the "
        "calibration set and a section of the cost surface.")
    onset = D.csv(ONSET)
    best = onset.loc[onset["cost"].idxmin()]
    R.p(f"As a check on the onset criterion, the level of normal force at which the vortex starts was "
        f"multiplied by factors from {D.note(onset['onset_level_factor'].min(), ONSET, 'smallest onset_level_factor', '{:.2f}')} "
        f"to {D.note(onset['onset_level_factor'].max(), ONSET, 'largest onset_level_factor', '{:.2f}')}; the "
        f"calibration cost is lowest, {D.note(best['cost'], ONSET, 'smallest cost', '{:.3f}')}, at a factor of "
        f"{D.note(best['onset_level_factor'], ONSET, 'onset_level_factor at the smallest cost', '{:.2f}')} "
        f"({ONSET}).")
    R.fig(D.root/FIGS/"fig14_calibration.png",
          "The calibration: cost of every start, sensitivity, cross-validation and cost surface "
          "(results/figures).", max_h_in=5.6)


def targets_rows(D) -> list:
    """The six held-out measures as rows of the accuracy table."""
    rows = []
    for key, (label, fmt, unit) in front.MEASURES.items():
        met = str(D.csv(TARGETS).set_index("measure").loc[key, "met"]) == "True"
        D.note(met, TARGETS, f"met [measure={key}]", "{}")
        rows.append([label, D.c(TARGETS, "measured", fmt, measure=key) + unit,
                     f"{D.c(TARGETS, 'kind', '{}', measure=key)} {D.c(TARGETS, 'target', fmt, measure=key)}{unit}",
                     "met" if met else "not met"])
    return rows


def accuracy(R, D) -> None:
    """Subsection: accuracy on the held-out loops."""
    h = front.held_out(D)
    n_meas = D.note(len(D.csv(TARGETS)), TARGETS, "number of rows", "word")
    floor = {q: D.c(FLOOR, "measurement_floor", "{:.3f}", metric=q) for q in ("nRMS_CL", "nRMS_CM", "nRMS_CD")}
    R.h2("Accuracy on the held-out loops")
    R.p(f"The model is scored on {h['n']} held-out NACA 0012 loops, each run at its own measured conditions. "
        f"A loop error is the root-mean-square difference between the predicted and the measured loop, "
        f"divided by the measured range of the quantity over the loop. {R.next_tab()} gives the {n_meas} "
        "measures with the targets set for them.")
    R.table(["Measure", "Measured on the held-out loops", "Target", "Result"], targets_rows(D),
            f"Accuracy on the held-out loops against the targets ({TARGETS}).")
    R.note("<b>The lift and moment loop-error targets and the damping-sign target are not met.</b> The drag "
           "loop-error target and the targets on the incidence of maximum lift and of moment stall are met.")
    R.p(f"The error a perfect model would show from the stated measurement uncertainty alone is "
        f"{floor['nRMS_CL']} in lift, {floor['nRMS_CM']} in moment and {floor['nRMS_CD']} in drag ({FLOOR}); "
        f"the errors of the model are several times that. Over every target set for the model, including "
        f"those on attached flow, structure, code and numerics, {front.scoreboard_sentence(D)} ({BOARD}); "
        f"{R.next_tab()} lists those.")
    sb = D.csv(BOARD)
    D.whole(BOARD, "rows with met other than yes")
    cols = ["group", "kind_of_target", "target", "kind", "value", "measured", "met"]
    header, rows = frame_rows(sb[sb["met"] != "yes"][cols])
    R.table(["Group", "Kind of target", "Target", "Kind", "Target value", "Measured", "Met"],
            [[c.replace("_", " ") for c in r] for r in rows],
            f"Targets not met ({BOARD}; names as in the file, with spaces for underscores).")
    R.p(f"{R.next_fig()} sets the held-out errors beside their intervals, the measurement floor and the "
        f"targets, and {R.next_fig(2)} shows measured and predicted loops for some of the held-out loops.")
    R.fig(D.root/FIGS/"fig16_accuracy_summary.png",
          "Held-out errors with their confidence intervals, the measurement floor and the targets (results/figures).",
          width=0.72)
    R.fig(D.root/FIGS/"fig05_held_out_loops.png",
          "Measured and predicted lift, moment and drag loops for held-out loops, with both separation laws "
          "(results/figures).", max_h_in=7.9)
    R.p(f"{R.next_fig()} sets measured and predicted loops side by side for loops chosen by a fixed rule from "
        f"{front.VTABLE}: the held-out loops with the lowest, the median and the highest loop error in lift, and "
        f"the calibration loop that Case A is set beside. The folder 06_postprocessing/validation holds that "
        f"figure, the static comparison, a row for every scored loop and the targets with their verdicts; they "
        f"are the same numbers as those of results/, set out beside the case study, not a second assessment.")
    R.fig(D.root/"06_postprocessing/validation/fig_validation_loops.png",
          "Measured against predicted loops of lift, moment and drag for loops chosen by rule; the last row is a "
          "calibration loop (06_postprocessing/validation).", max_h_in=9.0)
    depth = D.csv(DEPTH)
    D.whole(DEPTH, "all rows")
    _, rows = frame_rows(depth)
    R.p(f"{R.next_tab()} groups the same loops by how far the peak incidence goes beyond static stall.")
    R.table(["Group", "Loops", "Loop error, lift", "Loop error, moment", "Loop error, drag",
             "Damping sign correct, fraction", "Mean size of measured damping"], rows,
            f"Held-out accuracy by stall depth ({DEPTH}).")


def attached_check(R, D) -> None:
    """A paragraph on the check against Theodorsen's solution."""
    t9 = D.md(TABLES_MD)["Table 9. Attached flow against Theodorsen's solution"]
    row = t9[t9["Variant"] == "the model (low-Mach blend)"].iloc[0]
    src = "Table 9, row 'the model (low-Mach blend)'"
    la = D.note(row["Lift amplitude, %"], TABLES_MD, f"{src}, Lift amplitude, %", "{:.2f}")
    lp = D.note(row["Lift phase, deg"], TABLES_MD, f"{src}, Lift phase, deg", "{:.2f}")
    ma = D.note(row["Moment amplitude, %"], TABLES_MD, f"{src}, Moment amplitude, %", "{:.2f}")
    mp = D.note(row["Moment phase, deg"], TABLES_MD, f"{src}, Moment phase, deg", "{:.2f}")
    R.p(f"In attached incompressible flow the model is compared with Theodorsen's solution "
        f"{R.cite('theodorsen1935')}: over the reduced frequencies tested the largest error is {la} % in lift "
        f"amplitude, {lp}° in lift phase, {ma} % in moment amplitude and {mp}° in moment phase ({TABLES_MD}, "
        "Table 9).")


def equivalence(R, D) -> None:
    """Subsection: tabulated against fitted separation law."""
    same = D.csv(SAME).dropna(subset=["verdict"])
    loop = same[same["metric"].str.contains("nRMS")]["paired_difference"].abs().max()
    timing = same[same["metric"].str.contains("dalpha")]["paired_difference"].abs().max()
    n = D.note(len(same), SAME, "number of rows with a verdict", "{:d}")
    d_loop = D.note(loop, SAME, "largest |paired_difference| over the nRMS rows", "{:.4f}")
    d_time = D.note(timing, SAME, "largest |paired_difference| over the dalpha rows", "{:.2f}")
    n_ref = D.note(int((same["verdict"] == "reference better").sum()), SAME, "rows with verdict 'reference better'",
                   "{:d}")
    n_tab = D.note(int((same["verdict"] == "tabulated better").sum()), SAME, "rows with verdict 'tabulated better'",
                   "{:d}")
    level = D.note(re.search(r"(\d+) %", " ".join(same["verdict"])).group(1), SAME, "level named in verdict", "{}")
    n_sets = D.note(same["constants"].nunique(), SAME, "number of distinct values of constants", "word")
    own = D.csv(OWN)
    differ = []
    for _, r in own[own["verdict"] == "reference better"].iterrows():
        label, fmt, unit = OWN_LABELS[r["metric"]]
        differ.append(f"{label} {D.note(r['tabulated'], OWN, 'tabulated [metric=' + r['metric'] + ']', fmt)}{unit} "
                      f"against {D.note(r['reference'], OWN, 'reference [metric=' + r['metric'] + ']', fmt)}{unit}")
    R.h2("Tabulated against fitted separation law")
    R.p(f"The one departure of this model from the usual method is the tabulated separation point. To see "
        f"what it changes, a reference model with the usual fitted exponential law is run on the same "
        f"held-out loops. With the same dynamic constants in both models, over {n_sets} sets of constants, "
        f"the two score alike against measurement: mean loop errors differ by at most {d_loop} and the incidence "
        f"of stall by at most {d_time}°, and the 95 % interval of every loop-error difference lies inside the "
        f"equivalence margin recorded beforehand in data/targets.json (one tenth of each accuracy target). Of "
        f"the {n} paired comparisons, {n_ref} favour the fitted law and {n_tab} the tabulated one at the "
        f"{level} % level ({SAME}).")
    eq = D.csv(EQUIV)
    parts, short = [], []
    for q, name in (("CL", "lift"), ("CM", "moment"), ("CD", "drag")):
        e = eq[eq["coefficient"] == q]
        where = f"[coefficient={q}]"
        parts.append(
            f"{name} {D.note(e['mean_rms_difference_over_range'].max(), EQUIV, 'largest mean_rms_difference_over_range ' + where, '{:.4f}')} "
            f"(upper 95 % limit {D.note(e['ci95_hi'].max(), EQUIV, 'largest ci95_hi ' + where, '{:.4f}')}, margin "
            f"{D.note(e['margin'].iloc[0], EQUIV, 'margin ' + where, '{:.3f}')}, worst loop "
            f"{D.note(e['worst_loop_rms_over_range'].max(), EQUIV, 'largest worst_loop_rms_over_range ' + where, '{:.3f}')})")
        if not e["loads_within_margin"].all():
            short.append(name)
    R.p(f"The loads themselves are compared directly, as the RMS difference between the two predictions over "
        f"the cycle divided by the measured range and averaged over the loops; the largest over the sets of "
        f"constants is: {'; '.join(parts)} ({EQUIV}).")
    R.note("<b>On this data set the tabulated and the fitted separation laws are equally accurate, and the "
           "loads they predict are close but not identical.</b> "
           + (f"At the recorded margin the loads are not shown to be equivalent in {' and '.join(short)}. "
              if short else "At the recorded margin the loads are equivalent in all three coefficients. ")
           + "The tabulated law is not shown to be more accurate than the fitted one.")
    if differ:
        R.p(f"Calibrated separately by the same procedure, the two models differ at the {level} % level in: "
            f"{'; '.join(differ)} (tabulated model first; the fitted-law model is the better in each). With "
            f"identical constants those differences shrink to the size given above, so they come from the "
            f"values each calibration gave the constants and not from the separation law ({OWN}). "
            f"{R.next_fig()} shows the paired differences.")
    R.fig(D.root/FIGS/"fig12_separation_laws.png",
          "Tabulated against fitted separation law: paired differences on the held-out loops (results/figures).")


def damping(R, D) -> None:
    """Subsection: cycle damping."""
    def d(col: str, group: str, fmt: str) -> str:
        return D.c(DAMP, col, fmt, group=group)

    deep = "deep stall (6 deg or more beyond static stall)"
    R.h2("Cycle damping")
    R.p("The cycle damping is the work the moment does on the pitching section over a cycle, made "
        "dimensionless with the pitch amplitude α<sub>1</sub>; a negative value means the air feeds the "
        "motion:")
    R.eq("damping")
    R.p(f"On the {d('n_frames', ALL_LOOPS, '{:d}')} held-out loops with a closed measured moment loop the sign "
        f"of the damping is predicted correctly on a fraction {d('sign_agreement', ALL_LOOPS, '{:.2f}')}; a "
        f"prediction of positive damping everywhere would score "
        f"{d('sign_agreement_if_always_positive', ALL_LOOPS, '{:.2f}')}. The correlation between predicted and "
        f"measured damping is {d('correlation', ALL_LOOPS, '{:.2f}')} and the root-mean-square error "
        f"{d('rms_error', ALL_LOOPS, '{:.3f}')}. {d('measured_negative', ALL_LOOPS, '{:d}')} loops measure "
        f"negative damping; the model predicts negative damping on "
        f"{d('predicted_negative', 'measured damping negative', '{:d}')} of them, and on "
        f"{d('predicted_negative', 'measured damping positive', '{:d}')} loops whose measured damping is "
        f"positive. On the {d('n_frames', deep, '{:d}')} loops in deep stall the sign is correct on a fraction "
        f"{d('sign_agreement', deep, '{:.2f}')} ({DAMP}; {R.next_fig()}).")
    R.note("<b>The model does not predict cycle damping through stall.</b> It cannot locate a stall-flutter "
           "boundary and is not offered for that purpose. Where this report gives a damping figure for a "
           "case, it is model output and not a statement about stall flutter.")
    R.fig(D.root/FIGS/"fig13_cycle_damping.png", "Predicted against measured cycle damping on the held-out loops "
          "(results/figures).", width=0.6, max_h_in=5.4)


def comparison(R, D) -> None:
    """Section: calibration and the comparison with measurement."""
    R.h1("Calibration and the comparison with measurement", "accuracy")
    calibration(R, D)
    accuracy(R, D)
    attached_check(R, D)
    equivalence(R, D)
    damping(R, D)


def check_step(D) -> None:
    """Stop the build if the reduced-time step the load model marches for a
    case is not the step the configuration gives for it."""
    from unistall import dsmodel
    for name, c in D.json(CFG)["cases"].items():
        out = dsmodel.solve(c["alpha_mean_deg"], c["alpha_amp_deg"], c["reduced_frequency_k"], c["mach_M"],
                            n_per_cycle=c["march"]["steps_per_cycle"], n_cycles=1, chord=c["chord_m"],
                            a_sound=float(D.csv(AIR).set_index("property").loc["a_sound", "value"]))
        marched, stated = float(out["meta"]["ds"]), float(c["march"]["step_semichords"])
        if abs(marched - stated) > STEP_TOLERANCE*stated:
            raise ValueError(f"{name}: the model marches a step of {marched} semichords, the configuration gives "
                             f"{stated}")


def order_sentence(D) -> str:
    """What the three finest marches of each case show of an order of
    convergence, in words."""
    df = D.csv(ORDER)
    seen = df[df["order_observed"]]
    n_seen = D.note(len(seen), ORDER, "number of rows with order_observed = True", "word")
    n_all = D.note(len(df), ORDER, "number of rows", "word")
    text = (f"An order of convergence was looked for in the three finest marches of each case (the reported "
            f"number of steps, twice it and four times it). It is observed in {n_seen} of the {n_all} quantities "
            f"examined")
    if len(seen):
        text += (f", where it lies between {D.note(seen['observed_order'].min(), ORDER, 'smallest observed_order', '{:.1f}')} "
                 f"and {D.note(seen['observed_order'].max(), ORDER, 'largest observed_order', '{:.1f}')}")
    rest = df[~df["order_observed"]]
    if len(rest):
        names = sorted({f"{ORDER_LABELS[q].lower().replace(', deg', '')} of {short_case(c)}"
                        for c, q in zip(rest["case"], rest["quantity"], strict=True)})
        text += (f"; in the others ({', '.join(names)}) the three marches do not change monotonically and no order "
                 f"is observed")
    worst = df.loc[df["uncertainty_pct_of_reported"].idxmax()]
    return (text + f". Where an order is observed the reported value is set against the Richardson estimate; "
            f"elsewhere the range of the three marches is taken. The largest discretisation uncertainty of a "
            f"reported value found in this way is "
            f"{D.note(worst['uncertainty_pct_of_reported'], ORDER, 'largest uncertainty_pct_of_reported', '{:.2f}')} % "
            f"({ORDER_LABELS[worst['quantity']].lower().replace(', deg', '')}, {short_case(worst['case'])}; {ORDER}).")


def residual_sentence(D, names: list) -> str:
    """How far the reported cycle repeats, from the residual files."""
    parts = []
    for n in names:
        rel = f"05_solution/convergence/residuals_{n}.csv"
        df = D.csv(rel)
        settled = df[df["cycles_marched"] >= int(df[df["reported"]]["cycles_marched"].iloc[0])]
        floor = settled[["max_abs_change_CL", "max_abs_change_CM", "max_abs_change_CD"]].max().max()
        last = D.note(int(df["cycles_marched"].max()), rel, "largest cycles_marched", "{:d}")
        parts.append(f"for {short_case(n)} the largest change of any load over the cycle, from the reported cycle "
                     f"to the last of {last} marched, is "
                     f"{D.note(floor, rel, 'largest max_abs_change over the rows from the reported cycle on', 'pow10')}")
    return ("The reported cycle repeats to this extent: " + "; ".join(parts) + ". That figure, and not the change "
            "of one particular cycle, is the floor of repeatability of the loop: at a turning point of the motion "
            "the sign of the pitch rate decides a branch of the model, and rounding can place it on either side "
            "of zero, so a change of this size can reappear after the loop has settled.")


def numerics(R, D) -> None:
    """Section: numerical resolution."""
    names = list(meta.CASES)
    check_step(D)
    R.h1("Numerical resolution", "numerics")
    R.p("There is no spatial mesh. The model is a set of first-order lags marched in time, so the only "
        "discretisation is the time step. The linear lags use exact exponential recurrences and the march has "
        "no stability limit; but the switching of the time-constant multipliers, the vortex leaving the chord "
        "and repeated shedding are resolved only to the nearest step, so the error of the march is set by "
        "those events and need not fall smoothly with the step.")
    parts = []
    for n in names:
        parts.append(f"{short_case(n)} is marched with {D.j(CFG, f'cases.{n}.march.steps_per_cycle', '{:d}')} "
                     f"steps per cycle, a step of "
                     f"{D.j(CFG, f'cases.{n}.march.step_semichords', '{:.4f}')} semichords, for "
                     f"{D.j(CFG, f'cases.{n}.march.cycles', '{:d}')} cycles")
    R.p("; ".join(parts) + f" ({CFG}). The last cycle is the one reported. The step the load model marches is "
        "checked against the configured one when this report is built.")
    n_fr = D.c(CONV, "n_frames", "{:d}", study="step size")
    finer = D.c(CONV, "default_compared_with", "{}", study="step size")
    longer = D.c(CONV, "default_compared_with", "{}", study="cycles marched")
    R.p(f"On the {n_fr} calibration loops, a march with {finer} changes the loop error of any loop by at most "
        f"{D.c(CONV, 'worst_change_in_loop_error', '{:.4f}', study='step size')}, and one with {longer} "
        f"changes it by at most {D.c(CONV, 'worst_change_in_loop_error', 'pow10', study='cycles marched')}; "
        f"the target for both is {D.c(CONV, 'target_change_in_loop_error', '{:.3f}', study='step size')} "
        f"({CONV}; {R.next_fig()}).")
    R.fig(D.root/FIGS/"fig17_convergence.png",
          "Loads and loop errors of the calibration loops against step size and cycles marched "
          "(results/figures).")
    refine = D.csv(REFINE)
    parts = []
    for n in names:
        finest = D.note(refine[refine["case"] == n]["steps_per_cycle"].max(), REFINE,
                        f"largest steps_per_cycle [case={n}]", "{:d}")
        pct = [D.c(REFINE, f"pct_from_finest_{q}", "{:.3f}", case=n, reported=True)
               for q in ("CL_max", "CM_min", "CD_max")]
        parts.append(f"For {short_case(n)} the reported march differs from one with {finest} steps per cycle by "
                     f"{pct[0]} % in peak lift, {pct[1]} % in minimum moment and {pct[2]} % in peak drag.")
    R.p(" ".join(parts) + f" {R.next_tab()} gives the refinement of the step for the cases and "
        f"{R.next_fig()} shows it.")
    cols = ["case", "steps_per_cycle", "step_semichords", "CL_max", "CM_min", "CD_max", "pct_from_finest_CL_max",
            "pct_from_finest_CM_min", "pct_from_finest_CD_max", "reported"]
    D.whole(REFINE, "columns " + ", ".join(cols))
    view = refine[cols].copy()
    view["case"] = [short_case(n).replace(" ", "\u00a0") for n in view["case"]]
    _, rows = frame_rows(view)
    R.table(["Case", "Steps per cycle", "Step, semichords", "Peak lift", "Minimum moment", "Peak drag",
             "Peak lift, % from finest", "Minimum moment, % from finest", "Peak drag, % from finest",
             "Reported march"], rows, f"Refinement of the time step ({REFINE}).", size=7.5)
    R.fig(D.root/PLOTS/"timestep_refinement.png",
          "Refinement of the time step for the cases: difference from the finest march (06_postprocessing).")
    R.p(order_sentence(D) + f" {R.next_tab()} gives the figures.")
    order = D.csv(ORDER)
    D.whole(ORDER, "all rows")
    rows = []
    for r in order.itertuples():
        rows.append([short_case(r.case).replace(" ", "\u00a0"), ORDER_LABELS[r.quantity], cell(r.value_reported),
                     power_of_ten(r.value_twice_the_steps - r.value_reported, 2),
                     power_of_ten(r.value_four_times_the_steps - r.value_twice_the_steps, 2),
                     f"{r.observed_order:.2f}" if r.order_observed else "none observed",
                     power_of_ten(r.richardson_estimate - r.value_reported, 2) if r.order_observed else "not formed",
                     power_of_ten(r.discretisation_uncertainty, 2), f"{r.uncertainty_pct_of_reported:.3f}"])
    R.table(["Case", "Quantity", "Value at the reported steps", "Change on doubling the steps",
             "Change on doubling them again", "Observed order", "Richardson estimate less the reported value",
             "Discretisation uncertainty", "Uncertainty, % of reported"], rows,
            f"Observed order of convergence, Richardson estimate and discretisation uncertainty of the reported "
            f"loads, from the three finest marches ({ORDER}).", size=7.5)
    R.p(residual_sentence(D, names) + f" {R.next_fig()} shows the cycle-to-cycle change.")
    R.fig(D.root/PLOTS/"convergence_residuals.png",
          "Cycle-to-cycle change of the loads against the number of cycles marched (06_postprocessing).",
          width=0.85)


def statespace(R, D) -> None:
    """Section: the state-space form."""
    def v(col: str, fmt: str) -> str:
        return D.note(D.csv(SS)[col].iloc[0], SS, col, fmt)

    lo = front.fixed_constant(D, "blend limits in (18a)").split(",")[1].strip()
    R.h1("The state-space form", "statespace")
    R.p(f"The same model is also written as {v('states', '{:d}')} first-order differential equations "
        f"(unistall/statespace.py), integrated by the {v('integrator', '{}').split(',')[0]} method with the "
        f"step of Eq. {R.eq_ref('step_rule')}. With s the travel in semichords, each lag is written on the "
        "quantity it delays. The circulatory lags of the incidence and of the pitch rate become")
    R.eq("ss_circulatory")
    R.p("the impulsive terms of the normal force and of the moment become")
    R.eq("ss_impulsive")
    R.eq("ss_moment")
    R.p("and the pressure lag, the boundary-layer lag, the lag of the incidence at which the static moment is "
        "read, and the vortex lift become")
    R.eq("ss_lags")
    R.eq("ss_vortex")
    R.p("Each recurrence of Section @@model@@ is the exact solution of the matching equation over a step with "
        "a linearly varying input, so the two forms are the same model. The multipliers on the time constants "
        "are evaluated at the start of each step and held over it, and the vortex clock is advanced at the "
        "end of the step; across a switch the two forms therefore differ by an amount that depends on the step.")
    R.p(f"On the {v('n_frames', '{:d}')} calibration loops the largest difference between the two forms in any "
        f"loop error is {v('worst_difference_in_loop_error', '{:.4f}')} (target "
        f"{v('target_difference_in_loop_error', '{:.3f}')}; mean "
        f"{v('mean_difference_in_loop_error', '{:.4f}')}), and in cycle damping "
        f"{v('worst_difference_in_damping', '{:.4f}')}. At any instant the largest difference is "
        f"{v('worst_difference_in_CL_pct_of_range', '{:.1f}')} % of the range in lift, "
        f"{v('worst_difference_in_CM_pct_of_range', '{:.1f}')} % in moment and "
        f"{v('worst_difference_in_CD_pct_of_range', '{:.1f}')} % in drag; those largest differences sit where the "
        f"model switches a time constant. The state-space march costs "
        f"{v('cpu_ratio_state_space_to_indicial', '{:.1f}')} times the CPU time of the indicial one ({SS}; "
        f"{R.next_fig()}).")
    R.p(f"The state-space form covers the whole Mach range of the indicial one, including the low-Mach blend: "
        f"{v('n_frames_below_mach_0_20', '{:d}')} of those loops are below Mach {lo}, down to Mach "
        f"{v('lowest_mach', '{:.2f}')}, where further states carry the incompressible circulatory lift. No "
        "dynamic-stall claim is made for those loops; they test only that the two forms agree.")
    R.fig(D.root/FIGS/"fig18_state_space.png", "The state-space form against the indicial march "
          "(results/figures).", max_h_in=5.6)
    torsion(R, D)


def torsion(R, D) -> None:
    """The coupled torsional-section result."""
    st = D.csv(STRUCT)
    n = D.note(len(st), STRUCT, "number of rows", "{:d}")
    agree = D.note(int(st["state_space_agrees_in_sign"].sum()), STRUCT,
                   "number of rows with state_space_agrees_in_sign = True", "{:d}")
    grows = D.note(int((st["free_log_growth_per_cycle_state_space"] > 0).sum()), STRUCT,
                   "number of rows with free_log_growth_per_cycle_state_space > 0", "{:d}")
    R.p(f"Coupled to a section free in torsion (unistall/structural.py), the pitch angle and its rate join the "
        f"state vector and the whole is integrated as one system. In {agree} of {n} conditions it gives the "
        f"same sign of growth or decay of the free motion as the step-by-step coupling of the indicial march; "
        f"the motion grows in {grows} of them ({STRUCT}; {R.next_tab()}). This shows that the two forms of the "
        "model agree. It is not a prediction of stall flutter: the growth or decay follows from the cycle "
        "damping of the model, which is not predicted through stall.")
    cols = ["alpha_mean_deg", "alpha_amp_deg", "k", "M", "cycle_damping", "free_log_growth_per_cycle",
            "free_log_growth_per_cycle_state_space", "state_space_agrees_in_sign"]
    D.whole(STRUCT, "columns " + ", ".join(cols))
    _, rows = frame_rows(st[cols])
    R.table(["Mean incidence, deg", "Amplitude, deg", "k", "Mach", "Cycle damping of the model",
             "Growth per cycle, indicial march", "Growth per cycle, state-space form", "Same sign"], rows,
            f"A section free in torsion: logarithmic growth of the motion per cycle from the two forms of the "
            f"model; model output only ({STRUCT}).")


def build(R, D) -> None:
    """Add the middle part of the report."""
    static_inputs(R, D)
    comparison(R, D)
    numerics(R, D)
    statespace(R, D)
