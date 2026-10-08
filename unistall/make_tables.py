# -*- coding: utf-8 -*-
"""
unistall / make_tables.py
---------------------
The tables of the study as one readable document, results/tables.md, written
from the result files. Nothing in it is typed by hand; the same numbers are in
the .csv and .json files named under each table.

Author: Akosa Samuel Onyejekwe (independent)

Usage:  python3 -m unistall.make_tables
"""
import json
import numpy as np
import pandas as pd

from unistall import dsmodel as dm
from unistall.paths import RESULTS

LABELS = {"mean_nRMS_CL": "Loop error, lift", "mean_nRMS_CM": "Loop error, moment",
          "mean_nRMS_CD": "Loop error, drag", "mean_abs_dalpha_CLmax_deg": "Incidence of maximum lift, deg",
          "mean_abs_dalpha_Mstall_deg": "Incidence of moment stall, deg",
          "damping_sign_agreement": "Fraction of loops with the damping sign correct"}
NAMES = {"Tp": "pressure lag T_p", "Tf0": "boundary-layer lag T_f", "Tv0": "vortex decay T_v",
         "Tvl": "vortex travel time T_vl", "eta": "chord-force recovery eta"}


HEADERS = {
    "C_N1": "Stall level CN1", "alpha_deg": "Incidence, deg", "alpha0_deg": "Mean incidence, deg",
    "amp_deg": "Amplitude, deg", "mean_deg": "Mean incidence, deg", "n_frames": "Loops", "k": "Reduced frequency k",
    "M": "Mach number", "nRMS_CL": "Lift loop error", "nRMS_CM": "Moment loop error", "nRMS_CD": "Drag loop error",
    "mean_nRMS_CL": "Mean lift loop error", "mean_nRMS_CM": "Mean moment loop error",
    "mean_nRMS_CD": "Mean drag loop error", "post_stall_curve": "Post-stall static curve",
    "kind_of_target": "Kind of target", "measurement_floor": "Error implied by the quoted uncertainty",
    "CD_static_data": "Static drag, measured", "CD_model_slow_sweep": "Drag, model slow sweep",
    "CD_model_without_CD0": "Drag, model slow sweep, without CD0", "error_without_CD0": "Error without CD0",
    "CC_min": "Lowest chord force", "CD_min": "Lowest drag", "CD_range": "Range of drag over the cycle",
    "max_abs_CN_sum": "Largest sum of normal force, loop and mirror image",
    "max_abs_CM_sum": "Largest sum of moment, loop and mirror image",
    "max_abs_CL_sum": "Largest sum of lift, loop and mirror image",
    "max_abs_CC_difference": "Largest difference of chord force, loop and mirror image",
    "max_abs_CD_difference": "Largest difference of drag, loop and mirror image",
    "fraction_of_cycle_CC_negative": "Fraction of cycle with chord force negative",
    "fraction_of_cycle_CD_negative": "Fraction of cycle with drag negative",
    "sign_agreement_if_always_positive": "Sign agreement of a prediction of positive damping everywhere",
    "worst_change_in_loop_error": "Largest change in a loop error", "worst_change_in_CL": "Largest change in lift",
    "worst_change_in_CM": "Largest change in moment", "worst_change_in_CD": "Largest change in drag",
    "worst_error_against_zero_step": "Largest error against zero step (extrapolated)",
    "target_change_in_loop_error": "Target", "default_compared_with": "Default resolution compared with",
    "in_both": "Moment stall in measurement and model", "in_neither": "In neither",
    "measured_only_missed_by_model": "Measured, missed by the model",
    "model_only_not_measured": "Predicted, not measured",
    "mean_abs_incidence_error_deg_where_both": "Mean incidence error where both, deg",
}


def _header(name: str) -> str:
    """A column name as a table header: the readable name if one is given in
    HEADERS, otherwise the name with its underscores read as spaces."""
    text = HEADERS.get(name, name.replace("_", " "))
    return text[:1].upper() + text[1:]


def _cell(value: object, spec: str) -> str:
    """One cell as text: booleans as yes or no, floats by `spec`, and any bar
    in text escaped so that it cannot be read as a column break."""
    if isinstance(value, (bool, np.bool_)):
        return "yes" if value else "no"
    if isinstance(value, float):
        return "" if pd.isna(value) else format(value, spec)
    return str(value).replace("|", "\\|")


def _md(df, fmt=None):
    """A DataFrame as a Markdown table."""
    fmt = fmt or {}
    cols = list(df.columns)
    lines = ["| " + " | ".join(_header(str(c)) for c in cols) + " |", "|" + "---|"*len(cols)]
    lines += ["| " + " | ".join(_cell(v, fmt.get(c, ".3f")) for c, v in zip(cols, row, strict=False)) + " |"
              for row in df.itertuples(index=False)]
    return "\n".join(lines)


def test_matrix():
    t = pd.read_csv(RESULTS/"validation_table.csv")
    t = t[t.model == "tabulated"]
    rows = []
    for (s, inside), d in t.groupby(["set", "in_mach_range"]):
        rows.append({"Set": s.replace("_", " "), "Inside the static Mach range": bool(inside), "Loops": len(d),
                     "Mach": f"{d.M.min():.2f} to {d.M.max():.2f}", "Reduced frequency": f"{d.k.min():.3f} to {d.k.max():.3f}",
                     "Mean incidence, deg": f"{d.alpha0_deg.min():.0f} to {d.alpha0_deg.max():.0f}",
                     "Amplitude, deg": f"{d.amp_deg.min():.0f} to {d.amp_deg.max():.0f}"})
    return pd.DataFrame(rows)


def static_inputs():
    by = pd.read_csv(RESULTS/"static_by_mach_naca0012.csv").sort_values("M")
    return pd.DataFrame({"Mach": by.M, "Lift-curve slope, per deg": by.CL_alpha_per_deg,
                         "Zero-lift incidence, deg": by.alpha_zero_lift_deg, "Maximum static C_L": by.CL_max,
                         "Static-stall incidence, deg": by.alpha_at_CL_max_deg, "C_N1": by.CN_at_CL_max,
                         "Source sweep": by.frame.str.replace("frame_", "")})


def constants():
    tab = json.load(open(RESULTS/"calibrated_constants.json"))
    refc = json.load(open(RESULTS/"calibrated_constants_reference.json"))
    rows = []
    for n, label in NAMES.items():
        spread = tab["first_stage_record"]["spread_over_starts"].get(n, [float("nan")]*2)
        rows.append({"Constant": label, "Literature": dm.DEFAULTS[n],
                     "Tabulated model": tab["constants"].get(n, tab["fixed"].get(n)),
                     "Fitted-law model": refc["constants"].get(n, refc["fixed"].get(n)),
                     "Status": ("fitted; the two refined starts differ by "
                                f"{tab['polished_constants_differ_by'].get(n, float('nan')):.3f}")
                     if n in tab["constants"] else "literature value",
                     "Range over starts": f"{spread[0]:.2f} to {spread[1]:.2f}"})
    rows.append({"Constant": "unsteady-moment factor g_M", "Literature": 1.0,
                 "Tabulated model": tab["fixed"]["cm_unsteady"], "Fitted-law model": refc["fixed"]["cm_unsteady"],
                 "Status": "empirical, from loops below static stall", "Range over starts": ""})
    return pd.DataFrame(rows)


def accuracy():
    tg = pd.read_csv(RESULTS/"validation_targets.csv").set_index("measure")
    summ = pd.read_csv(RESULTS/"validation_summary.csv")
    s = summ[(summ.model == "tabulated") & (summ.group == "primary")].set_index("metric")
    rows = []
    for key, label in LABELS.items():
        m = key if key in s.index else "mean_Xi_sign_agree"
        rows.append({"Measure": label, "Loops": int(s.n_frames[m]), "Measured": tg.measured[key],
                     "95 % interval": f"{s.ci95_lo[m]:.3f} to {s.ci95_hi[m]:.3f}",
                     "Target": f"{tg.kind[key]} {tg.target[key]:.2f}", "Met": bool(tg.met[key])})
    return pd.DataFrame(rows)


GROUPS = {"primary": "all primary held-out loops", "primary_at_calibrated_mach": "Mach 0.285 and above (where the calibration loops are)",
          "primary_below_calibrated_mach": "below Mach 0.285", "primary_not_repeating_calibration": "not repeating a calibration condition",
          "outside_mach": "held out, below the static Mach range (no claim)", "calibration": "calibration loops"}


def by_group():
    summ = pd.read_csv(RESULTS/"validation_summary.csv")
    s = summ[summ.model == "tabulated"]
    rows = []
    for g, label in GROUPS.items():
        d = s[s.group == g].set_index("metric")
        rows.append({"Group": label, "Loops": int(d.n_frames["mean_nRMS_CL"]), "Lift": d.value["mean_nRMS_CL"],
                     "Moment": d.value["mean_nRMS_CM"], "Drag": d.value["mean_nRMS_CD"],
                     "Incidence of maximum lift, deg": d.value["mean_abs_dalpha_CLmax_deg"],
                     "Incidence of moment stall, deg": d.value["mean_abs_dalpha_Mstall_deg"],
                     "Damping sign correct": d.value["mean_Xi_sign_agree"]})
    return pd.DataFrame(rows)


def separation_laws():
    own = pd.read_csv(RESULTS/"comparison.csv")
    own.insert(0, "constants", "each model's own calibration")
    d = pd.concat([own, pd.read_csv(RESULTS/"comparison_same_constants.csv")], ignore_index=True)
    d = d[d.metric.str.startswith("mean")]
    return pd.DataFrame({"Constants": d.constants, "Measure": d.metric.str.replace("mean ", ""), "Tabulated": d.tabulated,
                         "Fitted law": d.reference, "Difference": d.paired_difference,
                         "95 % interval": [f"{lo:+.4f} to {hi:+.4f}" for lo, hi in zip(d.ci95_lo, d.ci95_hi, strict=False)],
                         "Verdict": d.verdict})


def equivalence():
    e = pd.read_csv(RESULTS/"equivalence.csv")
    return pd.DataFrame({
        "Constants": e.constants, "Coefficient": e.coefficient, "Margin": e.margin,
        "Mean RMS difference / range": e.mean_rms_difference_over_range.map("{:.4f}".format),
        "95 % interval": [f"{lo:.4f} to {hi:.4f}" for lo, hi in zip(e.ci95_lo, e.ci95_hi, strict=False)],
        "Worst loop": e.worst_loop_rms_over_range.map("{:.4f}".format),
        "Largest difference at any instant (coefficient)": e.worst_loop_max_abs_raw.map("{:.4f}".format),
        "Loads within margin": e.loads_within_margin.map({True: "yes", False: "no"}),
        "Loop errors within margin": e.loop_error_within_margin.map({True: "yes", False: "no"})})


def build():
    att = json.load(open(RESULTS/"attached_moment_factor.json"))
    parts = [
        "# Tables\n\nAuthor: Akosa Samuel Onyejekwe (independent)\n\nWritten by `python3 -m unistall.make_tables` "
        "from the result files; no number is typed by hand.",
        "## Table 1. Measured loops used\n\nNACA 0012, NASA TM-84245 (`results/validation_table.csv`, `data/data_manifest.csv`).\n\n"
        + _md(test_matrix()),
        "## Table 2. Static inputs at the two Mach stations\n\n(`results/static_by_mach_naca0012.csv`)\n\n"
        + _md(static_inputs(), {"Mach": ".3f", "Lift-curve slope, per deg": ".4f", "Zero-lift incidence, deg": ".2f",
                                "Static-stall incidence, deg": ".1f"}),
        "## Table 3. Constants\n\n(`results/calibrated_constants.json`, `results/calibrated_constants_reference.json`)\n\n"
        + _md(constants()),
        "## Table 4. Unsteady moment below static stall\n\nCalibration loops whose peak incidence is below static stall; the "
        f"empirical factor {att['cm_unsteady']} is the value that minimises their moment loop error, and the loops with a closed "
        f"moment loop would alone choose {att['range_over_loops_with_a_closed_moment_loop'][0]} and "
        f"{att['range_over_loops_with_a_closed_moment_loop'][1]}. The static part is the damping the model gives with the "
        "unsteady moment terms removed, which comes from reading the static table at the lagged incidence "
        "(`results/attached_moment_factor.json`).\n\n"
        + _md(pd.DataFrame(att["frames"]).drop(columns=["moment_loop_closed", "M"]).rename(columns={
            "frame": "Frame", "damping_measured": "Damping, measured", "damping_without_factor": "Damping, full strength",
            "damping_with_factor": "Damping, with factor", "damping_static_part": "Damping, static part",
            "damping_unlagged_static": "Damping, full strength, static table unlagged",
            "nRMS_CM_without_factor": "Moment error, full strength", "nRMS_CM_with_factor": "Moment error, with factor",
            "nRMS_CM_unlagged_static": "Moment error, full strength, static table unlagged",
            "best_factor_this_loop": "Factor this loop would choose", "peak_vortex_normal_force": "Peak vortex normal force",
            "peak_CN_prime_over_CN1": "Peak lagged normal force / onset level"})),
        "## Table 5. Accuracy on the held-out loops\n\n(`results/validation_targets.csv`, `results/validation_summary.csv`)\n\n"
        + _md(accuracy()),
        "## Table 6. Accuracy by stall depth, by Mach number and without repeated conditions\n\n"
        "(`results/validation_by_stall_depth.csv`, `results/validation_summary.csv`)\n\n"
        + _md(pd.read_csv(RESULTS/"validation_by_stall_depth.csv")) + "\n\n" + _md(by_group()),
        "## Table 6a. Moment stall: measured, predicted, missed\n\nThe incidence error of moment stall in Table 5 is defined "
        "only for loops where both the measurement and the model show one on the up-stroke "
        "(`results/validation_moment_stall.csv`).\n\n" + _md(pd.read_csv(RESULTS/"validation_moment_stall.csv")),
        "## Table 6b. Intervals when repeated conditions are not counted as independent\n\nThe held-out means with a 95 % "
        "interval from resampling the loops one by one and from resampling groups of loops that repeat one condition "
        "(`results/validation_intervals_by_condition.csv`).\n\n"
        + (_md(pd.read_csv(RESULTS/"validation_intervals_by_condition.csv"))
           if (RESULTS/"validation_intervals_by_condition.csv").exists() else "(written by the validation)"),
        "## Table 7. Tabulated against fitted separation law\n\nMean absolute error on the held-out loops; difference is "
        "tabulated minus fitted (`results/comparison.csv`, `results/comparison_same_constants.csv`).\n\n"
        + _md(separation_laws(), {"Tabulated": ".4f", "Fitted law": ".4f", "Difference": "+.4f"}),
        "## Table 7a. Loads predicted by the two separation laws with identical constants\n\nRMS difference between the "
        "two predictions over the cycle, divided by the measured range of the coefficient and averaged over the held-out "
        "loops, against the margin recorded beforehand in `data/targets.json` (`results/equivalence.csv`).\n\n"
        + _md(equivalence()),
        "## Table 8. Cycle damping\n\n(`results/validation_damping.csv`)\n\n" + _md(pd.read_csv(RESULTS/"validation_damping.csv")),
        "## Table 9. Attached flow against Theodorsen's solution\n\nLargest error over the reduced frequencies tested "
        "(`results/attached_checks.csv`).\n\n" + _md(_theodorsen(), {c: ".2f" for c in
                                                                    ("Lift amplitude, %", "Lift phase, deg", "Moment amplitude, %", "Moment phase, deg")}),
        "## Table 9a. Compressible attached flow against exact linear theory\n\nStarting values after a step and the exact "
        "short-time lift, at the Mach numbers of the study, where Theodorsen's solution does not apply "
        "(`results/indicial_checks.csv`).\n\n" + _md(pd.read_csv(RESULTS/"indicial_checks.csv"), {
            "model": ".4f", "theory": ".4f", "error_pct": "+.2f"}),
        "## Table 9b. The model beside the published Leishman-Beddoes model\n\nTwo NACA 0012 cases at Mach 0.3 and reduced "
        "frequency 0.1 for which Leishman and Crouse (1989) plot their model against the measurements. Loop errors of the "
        "published curves and of this model against the same traced measured points, and the RMS difference between the two "
        "models over the cycle as a fraction of the measured range (`results/published_reference.csv`; the curves are "
        "tracings of printed figures).\n\n"
        + _md(pd.read_csv(RESULTS/"published_reference.csv").drop(columns=["source", "measured_points", "published_curve_points"]),
              {"loop_error_published_model": ".4f", "loop_error_this_model": ".4f",
               "rms_difference_between_models_over_range": ".4f"}),
        "## Table 10. Second aerofoil (Ames A-01)\n\nDynamic constants unchanged (`results/validation_other_airfoil.csv`).\n\n"
        + _md(_other_airfoil()),
        "## Table 11. Every target against its measured value\n\n(`data/targets.json`, `results/targets_scoreboard.csv`)\n\n"
        + _md(pd.read_csv(RESULTS/"targets_scoreboard.csv").fillna(""), {"value": "g", "measured": "g"}),
        "## Table 12. Step size and cycles marched\n\nThe default resolution against a march with four times as many steps "
        "and one with four more cycles, on the calibration loops (`results/convergence_summary.csv`).\n\n"
        + _md(pd.read_csv(RESULTS/"convergence_summary.csv"), {c: ".4f" for c in
                                                             ("worst_change_in_loop_error", "worst_change_in_CL", "worst_change_in_CM", "worst_change_in_CD")}),
        "## Table 12a. State-space form against the indicial march\n\nThe same model integrated as differential "
        "equations, on the calibration loops (`results/statespace_summary.csv`).\n\n"
        + _md(pd.read_csv(RESULTS/"statespace_summary.csv").T.reset_index().rename(columns={"index": "quantity", 0: "value"})
              .astype(str)),
        "## Table 12b. The model at its limits\n\nWritten by `python3 -m unistall.check_limits`. A loop and its mirror image "
        "about the zero-lift incidence: sums and differences that are zero for a model symmetric about it "
        "(`results/limit_mirror_pair.csv`).\n\n" + _md(pd.read_csv(RESULTS/"limit_mirror_pair.csv"), {
            "max_abs_CN_sum": ".1e", "max_abs_CM_sum": ".4f", "max_abs_CC_difference": ".4f", "max_abs_CL_sum": ".4f",
            "max_abs_CD_difference": ".4f", "CD_range": ".3f"})
        + "\n\nDrag of a slow sweep against the static drag data at Mach 0.30 (`results/limit_static_drag.csv`).\n\n"
        + _md(pd.read_csv(RESULTS/"limit_static_drag.csv").drop(columns=["sweep_k"]))
        + "\n\nLoops that cross zero lift (`results/limit_zero_lift.csv`).\n\n" + _md(pd.read_csv(RESULTS/"limit_zero_lift.csv"))
        + "\n\nSwitching rules on every usable loop (`results/limit_vortex_rules.csv`).\n\n" + _md(_vortex_rules()),
        "## Table 13. Post-stall static curve\n\nCalibration cost with each of the two measured post-stall states, at the "
        "literature stall constants (`results/static_curve_choice.csv`).\n\n" + _md(pd.read_csv(RESULTS/"static_curve_choice.csv")),
    ]
    return "\n\n".join(parts) + "\n"


def _vortex_rules():
    v = pd.read_csv(RESULTS/"limit_vortex_rules.csv")
    rows = []
    for name, d in list(v.groupby("set")) + [("all", v)]:
        rows.append({"Set": name.replace("_", "-"), "Loops": len(d),
                     "Loops with a vortex shed": int((d.sheddings > 0).sum()),
                     "Loops with repeated shedding": int((d.repeated_sheddings > 0).sum()),
                     "Most sheddings in a cycle": int(d.sheddings.max()),
                     "Loops on which the sign rule removes vortex lift": int((d.steps_with_vortex_lift_removed > 0).sum()),
                     "Largest vortex normal force removed": f"{d.largest_vortex_lift_removed.max():.4f}"})
    return pd.DataFrame(rows)


def _theodorsen():
    a = pd.read_csv(RESULTS/"attached_checks.csv")
    a = a[a.check == "Theodorsen"]
    g = a.groupby("variant", sort=False)[["lift_amp_err_pct", "lift_phase_err_deg", "moment_amp_err_pct",
                                          "moment_phase_err_deg"]].agg(lambda v: v.abs().max()).reset_index()
    g.columns = ["Variant", "Lift amplitude, %", "Lift phase, deg", "Moment amplitude, %", "Moment phase, deg"]
    return g


def _other_airfoil():
    oa = pd.read_csv(RESULTS/"validation_other_airfoil.csv")
    g = oa.groupby("static_inputs").agg(Loops=("frame", "nunique"), lift=("nRMS_CL", "mean"), moment=("nRMS_CM", "mean"),
                                        drag=("nRMS_CD", "mean"), peak=("CLmax_err_pct", lambda v: v.abs().mean())).reset_index()
    g.columns = ["Static inputs", "Loops", "Loop error, lift", "Loop error, moment", "Loop error, drag", "Peak-lift error, %"]
    return g


if __name__ == "__main__":
    text = build()
    (RESULTS/"tables.md").write_text(text, encoding="utf-8")
    print(f"[tables] results/tables.md written ({len(text.splitlines())} lines)")
