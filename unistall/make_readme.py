# -*- coding: utf-8 -*-
"""
unistall / make_readme.py
---------------------
Write README.md from the result files, so that no number in it is typed by
hand. tests/test_readme_numbers.py fails if README.md differs from what this
script produces.

Author: Akosa Samuel Onyejekwe (independent)

Usage:  python3 -m unistall.make_readme            write README.md
        python3 -m unistall.make_readme --print    print it instead
"""
import sys, json
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parent
from unistall.paths import DATA, DOCS, RESULTS
ROOT = HERE.parent

LABELS = {"mean_nRMS_CL": ("Loop error, lift", "{:.3f}"),
          "mean_nRMS_CM": ("Loop error, moment", "{:.3f}"),
          "mean_nRMS_CD": ("Loop error, drag", "{:.3f}"),
          "mean_abs_dalpha_CLmax_deg": ("Incidence of maximum lift", "{:.2f}°"),
          "mean_abs_dalpha_Mstall_deg": ("Incidence of moment stall", "{:.2f}°"),
          "damping_sign_agreement": ("Sign of cycle damping correct", "{:.0%}")}
SHORT = {"mean |nRMS_CL|": "lift loop error", "mean |nRMS_CM|": "moment loop error",
         "mean |nRMS_CD|": "drag loop error", "mean |CLmax_err_pct|": "peak-lift error (per cent)",
         "mean |dalpha_CLmax_deg|": "incidence of maximum lift (degrees)",
         "mean |dalpha_Mstall_deg|": "incidence of moment stall (degrees)"}


def _targets_rows(tg):
    rows = []
    for key, (label, fmt) in LABELS.items():
        r = tg.loc[key]
        rows.append(f"| {label} | {fmt.format(r.measured)} | {r.kind} {fmt.format(r.target)} | "
                    f"{'yes' if r.met else 'no'} |")
    return "\n".join(rows)


def _equivalence_paragraph(same, eq, n_prim, tg):
    """Finding 1, written from results/equivalence.csv: what the two
    separation laws share (accuracy) and where their loads differ."""
    verdicts = same[same.metric.str.startswith("mean")].verdict.value_counts()
    loops = same[same.metric.str.contains("nRMS")]
    worst_same = loops.paired_difference.abs().max()
    worst_angle = same[same.metric.str.contains("dalpha")].paired_difference.abs().max()
    verdicts = same[same.metric.str.startswith("mean")].verdict.value_counts()
    n_ref, n_tab = int(verdicts.get("reference better", 0)), int(verdicts.get("tabulated better", 0))
    fitted_lower = int((loops.paired_difference > 0).sum())
    names = {"CL": "lift", "CM": "moment", "CD": "drag"}
    parts, shown, not_shown = [], [], []
    for q, d in eq.groupby("coefficient", sort=False):
        w = d.loc[d.worst_loop_max_abs_raw.idxmax()]
        parts.append(f"{names[q]} {d.mean_rms_difference_over_range.min():.4f} to "
                     f"{d.mean_rms_difference_over_range.max():.4f} (upper 95 % limit {d.ci95_hi.max():.4f}, margin "
                     f"{d.margin.iloc[0]:.3f}; worst loop {d.worst_loop_rms_over_range.max():.3f}; largest difference "
                     f"at any instant {w.worst_loop_max_abs_raw:.3f} in the coefficient, {w.worst_loop_frame.replace('frame_', 'frame ')})")
        (shown if d.loads_within_margin.all() else not_shown).append(names[q])
    verdict = (f"At that margin the loads are equivalent in {' and '.join(shown)}" if shown else
               "At that margin the loads are not shown to be equivalent in any coefficient")
    if not_shown:
        verdict += f" and are not shown to be equivalent in {' and '.join(not_shown)}"
    title = ("equally accurate on this data set" if eq.loop_error_within_margin.all()
             else "not shown to be equally accurate on this data set")
    return (
        f"**1. Tabulated and fitted separation laws are {title}; the loads they predict are close, not identical.** "
        f"With the same dynamic constants in both models, over three sets of constants (literature values, and those "
        f"fitted for each model), a separation point tabulated from the static measurements and the standard "
        f"three-parameter exponential fit of the same measurements score alike against {n_prim} held-out loops: mean "
        f"loop errors differ by at most {worst_same:.4f} and stall timing by at most {worst_angle:.2f}°, and the 95 % "
        f"interval of every loop-error difference lies inside the equivalence margin recorded beforehand in "
        f"`data/targets.json` (one tenth of each accuracy target). Of the {int(verdicts.sum())} paired comparisons "
        f"{n_ref} favour the fitted law and {n_tab} the tabulated one at the 95 % level, and the fitted law has the "
        f"lower mean loop error in {fitted_lower} of {len(loops)} rows, so there is a small bias in its favour. "
        f"The predicted loads themselves are compared directly, as the RMS difference between the two predictions "
        f"over the cycle divided by the measured range, averaged over the loops: {'; '.join(parts)}. {verdict}. "
        f"For scale, the model's own error against measurement is {tg.measured['mean_nRMS_CL']:.2f} to "
        f"{tg.measured['mean_nRMS_CM']:.2f} "
        f"(`results/comparison_same_constants.csv`, `results/equivalence.csv`, `results/load_difference_per_frame.csv`)."
    )


def _split_sentence(T):
    """How the loops were divided, written from data/split.json and
    results/split_departures.csv, naming every frame fixed by name."""
    f = pd.read_csv(RESULTS/"split_fixed_frames.csv")
    m = pd.read_csv(RESULTS/"split_departures.csv").set_index(["group", "metric"]).value
    W = T["held_out_validation"]

    def num(frame):
        return frame.replace("frame_", "")
    fixed = ", ".join(f"{num(r.frame)} ({r.set_as_fixed.replace('_', '-')})" for r in f.itertuples())
    dep = f[f.departs_from_rule]
    departs = "; ".join(f"{num(r.frame)} is in the {r.set_as_fixed.replace('_', '-')} set where the rule gives "
                        f"{r.set_by_rule.replace('_', '-')}" for r in dep.itertuples())
    cols = [("mean_nRMS_CL", "lift"), ("mean_nRMS_CM", "moment"), ("mean_Xi_sign_agree", "damping sign")]

    def row(g):
        return ", ".join(f"{name} {m[(g, k)]:.4f}" for k, name in cols)
    limits = {"mean_nRMS_CL": W["mean_nRMS_CL_max"], "mean_nRMS_CM": W["mean_nRMS_CM_max"]}
    groups = ("split as fixed", "departing frames moved to the set the rule gives", "departing frames left out")
    same_verdict = all(len({m[(g, k)] <= lim for g in groups}) == 1 for k, lim in limits.items()) and \
        len({m[(g, "mean_Xi_sign_agree")] >= W["damping_sign_agreement_min"] for g in groups}) == 1
    return (
        f"- The split into calibration and held-out loops (`data/split.json`) is a rule on the frame number "
        f"(divisible by three: calibration) for every loop but {len(f)}, which were in use before the rule was "
        f"written and were fixed by name: {fixed}. {len(dep)} of them depart from the rule: {departs}. "
        f"Held-out measures with the split as fixed: {row(groups[0])}; with those {len(dep)} frames moved to "
        f"the set the rule gives (constants unchanged, so the frame moved out of calibration has been fitted): "
        f"{row(groups[1])}; with them left out: {row(groups[2])}. "
        + ("No target changes from met to not met or back. " if same_verdict else
           "At least one target changes between met and not met. ")
        + "(`results/split_fixed_frames.csv`, `results/split_departures.csv`)")


def _scoreboard_counts(board):
    """Targets met and not met, counted separately for the model, the
    numerical method and the software, with the missed ones named."""
    words = {"physical": "on the model against measurement or theory", "numerical": "on the numerical method",
             "software": "on the software"}
    parts = [f"{int((d.met == 'yes').sum())} of {len(d)} {words[k]}" for k, d in board.groupby("kind_of_target", sort=False)]
    missed = board[board.met != "yes"]
    names = "; ".join(f"`{r.target}` ({r.measured:g} against {r.kind} {r.value:g})" for r in missed.itertuples())
    return ("Targets met: " + ", ".join(parts) + f". Not met ({len(missed)}): {names}.")


def _attached_moment_paragraph(att_m, unc):
    """What the empirical unsteady-moment factor is, what it corrects and how
    far it is determined, from results/attached_moment_factor.json."""
    fr = att_m["frames"]
    closed = [r for r in fr if r["moment_loop_closed"]]
    lo, hi = att_m["range_over_loops_with_a_closed_moment_loop"]
    rows = "; ".join(
        f"frame {r['frame'].split('_')[1]}, k = {r['k']:.2f}: measured {r['damping_measured']:.3f}, model at full "
        f"strength {r['damping_without_factor']:.3f}, of which {r['damping_static_part']:.3f} comes from reading "
        f"the static moment table at the lagged incidence (with the table read at the incidence itself "
        f"{r['damping_unlagged_static']:.3f}), with the factor {r['damping_with_factor']:.3f}" for r in closed)
    share = [r["damping_static_part"]/(r["damping_without_factor"] - r["damping_measured"]) for r in closed]
    shed = [r for r in fr if r["peak_vortex_normal_force"] and r["peak_vortex_normal_force"] > 0]
    worse = [r for r in fr if r["nRMS_CM_with_factor"] > r["nRMS_CM_without_factor"]]
    b = att_m["basis"]
    text = (
        f"- *Before stall onset.* With its unsteady moment terms at full strength the model gives more pitch "
        f"damping than is measured in attached flow. An empirical factor of {att_m['cm_unsteady']:.3f} on those "
        f"terms is the value that minimises the moment error over the attached part of {b['n_loops']} calibration "
        f"loops ({b['n_points']} measured points that the model places before stall onset); its 95 % bootstrap "
        f"interval over loops is {b['interval_95'][0]:.2f} to {b['interval_95'][1]:.2f}. On the {len(fr)} loops whose "
        f"peak incidence is below static stall ({rows}), {min(share):.0%} to {max(share):.0%} of the excess damping "
        f"comes from the lagged reading of the static table; moving the attached part of that reading to the "
        f"circulatory incidence changes it little, so the lag is largely physical and the factor corrects an excess "
        f"of the inviscid unsteady moment over measurement. Taken alone, the {len(closed)} of those loops with a closed "
        f"moment loop would choose {lo:.2f} and {hi:.2f}. Over the interval of the factor "
        f"the held-out moment error moves between {unc['moment factor_min']['nRMS_CM']:.4f} and "
        f"{unc['moment factor_max']['nRMS_CM']:.4f} and the damping-sign agreement between "
        f"{unc['moment factor_min']['Xi_sign_agree']:.3f} and {unc['moment factor_max']['Xi_sign_agree']:.3f} "
        f"(`results/uncertainty_summary.csv`).")
    if shed:
        text += (" The loops are below static stall but not all free of stall in the model: on "
                 + ", ".join(f"frame {r['frame'].split('_')[1]}" for r in shed)
                 + f" the lagged normal force reaches {max(r['peak_CN_prime_over_CN1'] for r in shed):.2f} of the onset "
                   f"level and a vortex of normal force {max(r['peak_vortex_normal_force'] for r in shed):.3f} is shed.")
    if worse:
        text += (" The factor makes the moment error worse on "
                 + ", ".join(f"frame {r['frame'].split('_')[1]} (k = {r['k']:.2f}, {r['nRMS_CM_without_factor']:.3f} to "
                             f"{r['nRMS_CM_with_factor']:.3f})" for r in worse) + ".")
    return text + (" The factor is a correction of this model, from few loops at one Mach number, not a "
                   "measurement of the flow.")


def _reconstruction_sentence():
    """What the reconstructed surface pressure returns besides lift, from
    05_solution/metrics_<case>.csv."""
    parts = []
    for name, label in (("A_validation", "Case A"), ("B_application", "Case B")):
        v = pd.read_csv(ROOT/"05_solution"/f"metrics_{name}.csv").set_index("metric").value.astype(float)
        parts.append(f"{label}: moment {v['Cp_moment_minus_model_at_CLmax']:+.3f}, drag "
                     f"{v['Cp_drag_minus_model_at_CLmax']:+.3f}, same sign of moment at "
                     f"{int(v['Cp_moment_sign_agrees_with_model_instants'])} of {int(v['stored_instants'])} stored instants")
    return ("the pitching moment and pressure drag of its surface pressure do not agree with the model's "
            "(difference at peak lift, " + "; ".join(parts) + "), so no moment or drag should be read from it.")


def _limit_checks_paragraph():
    """What the model does where no measured loop tests it directly, from the
    files written by unistall/check_limits.py."""
    mp = pd.read_csv(RESULTS/"limit_mirror_pair.csv")
    sd = pd.read_csv(RESULTS/"limit_static_drag.csv")
    zero = pd.read_csv(RESULTS/"limit_zero_lift.csv")
    v = pd.read_csv(RESULTS/"limit_vortex_rules.csv")
    below = sd[sd.alpha_deg.between(0.0, 12.5)]
    worst_zero = zero.loc[zero.fraction_of_cycle_CD_negative.idxmax()]
    return (
        "Behaviour at the limits (`python3 -m unistall.check_limits`, Table 12b). "
        f"*Mirror image:* normal force and moment are odd about the zero-lift incidence to {mp.max_abs_CN_sum.max():.0e} "
        f"in normal force; lift and drag are resolved on the experiment's incidence scale, so a loop and its mirror "
        f"image differ in drag by up to {mp.max_abs_CD_difference.max():.3f} (`results/limit_mirror_pair.csv`). "
        f"*Static limit:* the drag of a slow sweep is not the measured static drag, which is total drag from a wake "
        f"survey where the loops are pressure drag: between 0 and 12° the model is higher by up to "
        f"{below.error.max():.3f}, and through stall the two differ by up to {sd.error.abs().max():.2f} "
        f"(`results/limit_static_drag.csv`). "
        f"*Zero lift:* on the {len(zero)} loops that cross it the lowest chord force is {zero.CC_min.min():.4f} and the "
        f"lowest drag {zero.CD_min.min():+.4f}; drag is below zero for {worst_zero.fraction_of_cycle_CD_negative:.3f} of "
        f"the cycle on frame {worst_zero.frame.split('_')[1]} ({worst_zero.set.replace('_', '-')}) "
        f"(`results/limit_zero_lift.csv`). "
        f"*Switching rules:* the vortex clock restarts within a cycle on {int((v.repeated_sheddings > 0).sum())} of "
        f"{len(v)} loops, up to {int(v.sheddings.max())} sheddings in a cycle on the slowest; the rule that removes "
        f"vortex lift opposing the separated force acts on {int((v.steps_with_vortex_lift_removed > 0).sum())} loops and "
        f"removes at most {v.largest_vortex_lift_removed.max():.4f} in normal force (`results/limit_vortex_rules.csv`).")


def _refined_starts_sentence(cal):
    """The spread of the refined starts of the first-stage fit."""
    costs = sorted(cal.get("first_stage_record", cal)["polished_costs"])
    near = sum(c <= 1.01*costs[0] for c in costs)
    return (f"The {len(costs)} starts refined at full resolution reach costs from {costs[0]:.3f} to {costs[-1]:.3f}; "
            f"{near} of them are within 1 % of the lowest, so the constants are one of several near-equal "
            f"solutions (`results/calibration_runs.csv` lists them all).")


def _other_airfoil_sentence(own, naca):
    """Every measure of the second aerofoil, with its own static data and with
    those of the NACA 0012, whichever way each moves."""
    rows = (("lift loop error", "nRMS_CL", False, ".3f", ""), ("moment loop error", "nRMS_CM", False, ".3f", ""),
            ("drag loop error", "nRMS_CD", False, ".3f", ""), ("peak-lift error", "CLmax_err_pct", True, ".1f", " %"))
    parts = []
    for label, col, absolute, fmt, unit in rows:
        a = own[col].abs().mean() if absolute else own[col].mean()
        b = naca[col].abs().mean() if absolute else naca[col].mean()
        parts.append(f"{label} {a:{fmt}}{unit} against {b:{fmt}}{unit} ({'better' if a < b else 'worse'})")
    return "; ".join(parts)


def _candidates_paragraph():
    """Changes to the model that were tried on the calibration loops, with the
    cross-validated cost of each beside the present form, from
    results/candidates_tried.csv (written by unistall/try_candidates.py)."""
    f = RESULTS/"candidates_tried.csv"
    if not f.exists():
        return ""
    c = pd.read_csv(f)
    base = c[c.candidate == "present form"].iloc[0]
    rows = "; ".join(f"{r.candidate}: {r.left_out_cost:.3f} ({'adopted' if r.adopted else 'not adopted'})"
                     for r in c[c.candidate != "present form"].itertuples())
    return ("\n\n**What was tried.** Changes to the form of the model were judged on the calibration loops alone, by the "
            f"cross-validation rule of `data/targets.json`; the present form has a left-out cost of {base.left_out_cost:.3f}. "
            f"{rows} (`results/candidates_tried.csv`). The held-out loops played no part in these choices.")


def _identifiability_sentence(cal):
    """How far the calibration loops determine each constant: the change of
    cost for a 10 % change, and the spread over the refined starts."""
    sel = cal["selection"]
    first = cal.get("first_stage_record", cal)
    sens = sel["largest_cost_change_pct_for_10pct"]
    spread = first.get("polished_constants_differ_by", {})
    parts = "; ".join(f"{n}: {sens[n]:.2f} %" + (f", refined starts differ by {spread[n]:.2f}" if n in spread else "")
                      for n in sens)
    steps = "; ".join(f"{st['larger']} against {st['smaller']}: mean improvement {st['mean_improvement_per_fold']:.4f} "
                      f"a fold, standard error {st['standard_error']:.4f}, "
                      f"{'kept' if st['larger_kept'] else 'not kept'}" for st in sel.get("steps", []))
    text = f"Change of the calibration cost for a 10 % change of each constant in the fit of all five: {parts}."
    if steps:
        text += f" Cross-validation over {len(sel['sensitive_set_by_fold'])} folds ({steps})."
    return text


def _target_origin(T):
    """How each held-out limit was formed, from data/targets.json."""
    o = T["held_out_validation"]["origin"]
    names = {"mean_nRMS_CL_max": "lift", "mean_nRMS_CM_max": "moment", "mean_nRMS_CD_max": "drag",
             "mean_abs_dalpha_CLmax_deg_max": "incidence of maximum lift",
             "mean_abs_dalpha_Mstall_deg_max": "incidence of moment stall"}
    return "; ".join(f"{label} {o[k]['absolute']:g} and {o[k]['three_quarters_of_baseline']:g} (baseline {o[k]['baseline']:g})"
                     for k, label in names.items())


def _published_reference_paragraph():
    """Where this implementation stands beside the published model, from
    results/published_reference.csv."""
    d = pd.read_csv(RESULTS/"published_reference.csv")
    names = {"CN": "normal force", "CM": "moment"}
    parts = "; ".join(f"{r.case.split(',')[0]}, {names[r.coefficient]}: published model {r.loop_error_published_model:.3f}, "
                      f"this model {r.loop_error_this_model:.3f}" for r in d.itertuples())
    better = int((d.loop_error_this_model < d.loop_error_published_model).sum())
    return ("Beside the published model (`python3 -m unistall.published_reference`, Table 9b). For two NACA 0012 cases "
            "at Mach 0.3 and reduced frequency 0.1, Leishman and Crouse (1989) plot the Leishman–Beddoes model against "
            "the measurements. Scored against the same traced points with the same loop error: " + parts
            + f". This model has the lower error in {better} of {len(d)} comparisons; the two models differ from each "
            f"other by {d.rms_difference_between_models_over_range.min():.3f} to "
            f"{d.rms_difference_between_models_over_range.max():.3f} of the measured range. The published curves are "
            "for the constants their authors chose for these conditions; this model's constants are fitted to 20 loops "
            "at once. Neither case is a new test of this model (`results/published_reference.csv`).")


def _own_calibration_sentence(own, same):
    """What the comparison shows when each model has its own calibration, and
    whether those differences survive identical constants: for each measure
    that differs with own constants, the largest difference with identical
    constants over the three sets, and whether any of those excludes zero."""
    d = own[own.metric.isin(SHORT)]
    diff = d[d.verdict != "no difference at 95 %"]
    if not len(diff):
        return ("Calibrated separately by the same script, the two models also score the same: no "
                "measure differs at the 95 % level.")
    names = {"tabulated better": "tabulated law better", "reference better": "fitted law better"}
    parts, survive = [], []
    for r in diff.itertuples():
        ident = same[same.metric == r.metric]
        biggest = ident.loc[ident.paired_difference.abs().idxmax()]
        parts.append(f"{SHORT[r.metric]} {r.tabulated:.3f} against {r.reference:.3f} ({names[r.verdict]}; with identical "
                     f"constants the difference is at most {abs(biggest.paired_difference):.3f})")
        if (ident.verdict == r.verdict).any() and abs(biggest.paired_difference) >= 0.5*abs(r.paired_difference):
            survive.append(SHORT[r.metric])
    text = ("Calibrated separately by the same script, the two models differ at the 95 % level in: "
            + "; ".join(parts) + ". ")
    if survive:
        return text + ("For " + ", ".join(survive) + " at least half of the difference remains with identical "
                       "constants, so it belongs to the separation law; the rest come from the values each "
                       "calibration gave the constants.")
    return text + ("With identical constants each of those differences falls to less than half its size, so they "
                   "come from the values each calibration gave the constants and not from the separation law.")


def _constants_sentence(cal):
    fitted = ", ".join(f"{k} = {v:.3f}" for k, v in cal["constants"].items()) or "none"
    dropped = ", ".join(cal["selection"]["returned_to_literature_values"]) or "none"
    gain = 100*(cal["cost_at_literature_values"] - cal["cost"])/cal["cost_at_literature_values"]
    return fitted, dropped, gain


def _depth_rows(depth):
    return "\n".join(f"| {r.group} | {r.n_frames} | {r.nRMS_CL:.3f} | {r.nRMS_CM:.3f} | {r.nRMS_CD:.3f} | "
                     f"{r.damping_sign_agreement:.0%} |" for r in depth.itertuples())


def _not_met(tg):
    """One bullet per target that is not met, with the measured value."""
    out = []
    for key, (label, fmt) in LABELS.items():
        r = tg.loc[key]
        if not r.met:
            out.append(f"- **{label}:** {fmt.format(r.measured)} against a target of {r.kind} "
                       f"{fmt.format(r.target)}.")
    return "\n".join(out) if out else "- Every target on the held-out loops is met."


def _case_rows():
    """The two cases' headline numbers, from 05_solution/summary_all_cases.csv."""
    f = ROOT/"05_solution"/"summary_all_cases.csv"
    if not f.exists():
        raise FileNotFoundError(f"{f} is missing: run python3 run_all.py before building the README")
    names = {"A_validation": "A: tunnel condition", "B_application": "B: rotor-blade-station condition"}
    rows = ["| Case | Peak lift | at incidence | Minimum moment | Peak drag | Stall onset |",
            "|---|---|---|---|---|---|"]
    for r in pd.read_csv(f).itertuples():
        rows.append(f"| {names.get(r.case, r.case)} | {r.CL_max:.3f} | {r.alpha_at_CLmax_deg:.1f}° | {r.CM_min_c4:.3f} | "
                    f"{r.CD_max:.3f} | {r.stall_onset_alpha_deg:.1f}° |")
    return "\n".join(rows)


def _first_plot(prefix):
    """Path of the first plot whose name starts with `prefix`."""
    hits = sorted((ROOT/"06_postprocessing"/"plots").glob(prefix + "*.png"))
    if not hits:
        raise FileNotFoundError(f"no plot named {prefix}*.png in 06_postprocessing/plots: run python3 run_all.py")
    return f"06_postprocessing/plots/{hits[0].name}"


def _onset_sentence():
    """How the calibration cost moves with the stall-onset level."""
    d = pd.read_csv(RESULTS/"onset_level_check.csv")
    best = d.loc[d.cost.idxmin()]
    others = d[d.onset_level_factor != best.onset_level_factor]
    return (f"{others.onset_level_factor.min():.2f} to {others.onset_level_factor.max():.2f} raises the calibration "
            f"cost from {best.cost:.3f} at a factor of {best.onset_level_factor:.2f} to between "
            f"{others.cost.min():.3f} and {others.cost.max():.3f}")


def _numerics_numbers():
    """The step of the march, the state-space check and the count of
    torsional conditions on which the two couplings agree."""
    from unistall import dsmodel
    stc = pd.read_csv(RESULTS/"structural_checks.csv")
    return (dsmodel.DS_MAX, pd.read_csv(RESULTS/"statespace_summary.csv").iloc[0],
            stc.state_space_agrees_in_sign.sum(), len(stc))


def build():
    import project_meta as meta
    tg = pd.read_csv(RESULTS/"validation_targets.csv").set_index("measure")
    same = pd.read_csv(RESULTS/"comparison_same_constants.csv")
    own = pd.read_csv(RESULTS/"comparison.csv")
    eq = pd.read_csv(RESULTS/"equivalence.csv")
    depth = pd.read_csv(RESULTS/"validation_by_stall_depth.csv")
    by = pd.read_csv(RESULTS/"static_by_mach_naca0012.csv").sort_values("M")
    man = pd.read_csv(DATA/"data_manifest.csv")
    tab = pd.read_csv(RESULTS/"validation_table.csv")
    unc = pd.read_csv(RESULTS/"uncertainty_summary.csv").set_index("metric")
    att = pd.read_csv(RESULTS/"attached_checks.csv")
    oa = pd.read_csv(RESULTS/"validation_other_airfoil.csv")
    cal = json.load(open(RESULTS/"calibrated_constants.json"))
    refs = pd.read_csv(DOCS/"bibliography_resolution.csv")
    quality = pd.read_csv(RESULTS/"code_quality.csv").set_index("measure")["value"]
    T = json.load(open(DATA/"targets.json"))

    t = tab[tab.model == "tabulated"]
    n_prim = int(((t["set"] == "held_out") & t.in_mach_range).sum())
    n_verified = int((man.verification == "verified").sum())
    n_in = int(t.in_mach_range.sum())
    from unistall.validate import MACH_CALIBRATED as MACH_CAL
    n_cal = int(((t["set"] == "calibration") & t.in_mach_range).sum())
    th = att[att.check == "Theodorsen"]
    model_th = th[th.variant.str.startswith("the model")]
    alone = th[th.variant.str.startswith("compressible indicial constants alone")]
    own_static = oa[oa.static_inputs == "own static sweep"]
    naca_static = oa[oa.static_inputs != "own static sweep"]
    lo, hi = by.iloc[0], by.iloc[-1]
    fitted, dropped, gain = _constants_sentence(cal)
    damp = pd.read_csv(RESULTS/"validation_damping.csv").set_index("group")
    d_all = damp[damp.index.str.startswith("all primary")].iloc[0]
    d_pos, d_neg = damp.loc["measured damping positive"], damp.loc["measured damping negative"]
    d_deep = damp[damp.index.str.startswith("deep stall")].iloc[0]
    d_rand = damp[damp.index.str.contains("random")].iloc[0]
    d_sure = damp[damp.index.str.contains("follows the stroke")].iloc[0]
    att_m = json.load(open(RESULTS/"attached_moment_factor.json"))
    stall = pd.read_csv(RESULTS/"validation_moment_stall.csv").set_index("model").loc["tabulated"]
    summ = pd.read_csv(RESULTS/"validation_summary.csv")
    grp = {g: d.set_index("metric") for g, d in summ[summ.model == "tabulated"].groupby("group")}
    hi_m, lo_m, fresh = (grp[g] for g in ("primary_at_calibrated_mach", "primary_below_calibrated_mach",
                                          "primary_not_repeating_calibration"))
    prim = t[(t["set"] == "held_out") & t.in_mach_range]
    conv = pd.read_csv(RESULTS/"convergence_summary.csv").set_index("study")
    board = pd.read_csv(RESULTS/"targets_scoreboard.csv")
    curve = pd.read_csv(RESULTS/"static_curve_choice.csv").set_index("post_stall_curve").cost
    dm_ds, ssf, struct_ok, struct_n = _numerics_numbers()
    S = T["separated_flow"]
    t_in = t[t.in_mach_range]
    drag_all_ok = bool((t_in.CD_neg_frac.max() <= S["CD_negative_fraction_of_cycle_max"])
                       and (t_in.CD_min_model.min() >= S["CD_min"]))

    return f"""# {meta.SOFTWARE_TITLE}
### {meta.SOFTWARE_SUBTITLE}

**Author:** {meta.AUTHOR} ({meta.AFFILIATION})

---

## What this repository contains

UNISTALL (Python package `unistall`) computes the unsteady lift, drag and pitching
moment on an aerofoil pitching through stall, by the Leishman–Beddoes method
in indicial and in state-space form, with the trailing-edge separation point
and the static moment read from measured static data instead of a fitted
curve.

| Where | What |
|---|---|
| [`unistall/dsmodel.py`](unistall/dsmodel.py) | **The solver.** It marches the loads in time; `attached_flow.py` and `static_model.py` supply the attached-flow loads and the static inputs |
| [`results/figures/`](results/figures) | **Every figure**: loops, contour maps, error plots, calibration, model states |
| [`results/`](results) | **Every table of results** (`.csv`, `.json`), each written by one script |
| [`data/`](data) | Inputs: verified conditions of the measured loops, digitised static data, calibration and held-out sets, targets |
| [`docs/formulation.md`](docs/formulation.md) | Every equation and constant |
| [`tests/`](tests) | Test suite |

The measurements are those of McCroskey, McAlister, Carr & Pucci (1982), *An
Experimental Study of Dynamic Stall on Advanced Airfoil Sections*, NASA
TM-84245: NACA 0012 and Ames A-01 sections, chord 0.61 m, Mach numbers up to
0.30. Of their NACA 0012 oscillating-aerofoil loops, {n_in} lie inside the Mach
range of the static data ({n_cal} used for calibration, {n_prim} held out) and
carry the claims made here; {len(t) - n_in} more, below that range, are scored
and reported but support no claim. {oa.frame.nunique()} Ames A-01 loops are scored
as a second section.

It is a load model, not a flow solver: it returns forces and moments on the
section. The model itself has no mesh and computes no flow field. The case
study adds a reconstructed field, drawn round the predicted lift for
illustration (see Case study); no load depends on it.

---

## Results at a glance

Every figure and table below is written by a script that runs the solver or
reads a file the solver wrote; a test fails if one is present that no script
draws. What is not computed by the solver is measured: the points of the
oscillating loops it is compared with, and the static data it reads as inputs
(the digitised static lift, moment and drag, and the slope, zero-lift
incidence and stall level taken from the quasi-steady sweeps).

![Measured and predicted loops](results/figures/fig05_held_out_loops.png)

*Fig. 5. Measured and predicted lift, moment and drag loops for three held-out
frames chosen by rule from the model's own lift error: the lowest, the median
and the highest. Each panel gives the loop error and the largest difference
between the tabulated and the fitted separation law; the Mach number of each
loop is in its title.*

![Dynamic loops against the static curves](results/figures/fig07_dynamic_and_static.png)

*Fig. 7. Predicted lift, moment and drag loops, deep stall and light stall,
against the static curves the model reads (lines) and the measured static
points (symbols). Static drag is measured and is not an input of the model's
drag. Three surface-pressure drag points between 14.0 and 15.0 degrees lie
above their neighbours; they are plotted as measured.*

![Design-space contours](results/figures/fig08_design_space.png)

*Fig. 8. Predicted peak lift, minimum moment and peak drag over mean incidence
and reduced frequency, as filled maps with a scale bar. The roughness above
about 15 degrees mean incidence is the model's switching and repeated
shedding, not noise in the plot.*

| Figure | Shows |
|---|---|
| [Fig. 1](results/figures/fig01_static_inputs.png) | Static normal force, separation point and moment at the two Mach stations |
| [Fig. 2](results/figures/fig02_attached_flow.png) | Attached-flow loads against Theodorsen's solution |
| [Fig. 3](results/figures/fig03_attached_moment.png) | Moment loops that peak below static stall: unsteady terms at full strength and with the empirical factor (literature stall constants; frame 7110 reaches the model's onset level, which is the kink near its peak) |
| [Fig. 4](results/figures/fig04_calibration_loops.png) | Three calibration loops with their RMS errors. The secondary bumps and slope breaks in the predicted loops come from repeated shedding and the switching of time constants; the measured loops have no counterpart to them |
| [Fig. 5](results/figures/fig05_held_out_loops.png) | Three held-out loops chosen by rule, both separation laws |
| [Fig. 6](results/figures/fig06_effect_of_constants.png) | Effect of the constants on one held-out deep-stall loop |
| [Fig. 7](results/figures/fig07_dynamic_and_static.png) | Dynamic loops against the static curves |
| [Fig. 8](results/figures/fig08_design_space.png) | Contours of the predicted peaks over mean incidence and reduced frequency |
| [Fig. 9](results/figures/fig09_damping_map.png) | Cycle damping over the same plane, with the measured loops |
| [Fig. 10](results/figures/fig10_mach_trends.png) | Contours of the predicted peaks over reduced frequency and Mach number |
| [Fig. 11](results/figures/fig11_error_by_condition.png) | Loop error against stall depth and frequency |
| [Fig. 12](results/figures/fig12_separation_laws.png) | Tabulated against fitted separation law, paired differences |
| [Fig. 13](results/figures/fig13_cycle_damping.png) | Predicted against measured cycle damping |
| [Fig. 14](results/figures/fig14_calibration.png) | Cost of every start, sensitivity, cross-validation, cost surface |
| [Fig. 15](results/figures/fig15_model_states.png) | The solver's internal states over one cycle, with the number of times the vortex clock starts in the cycle and how many of those are repeated sheddings |
| [Fig. 16](results/figures/fig16_accuracy_summary.png) | Held-out errors, 95 % intervals, measurement floor and targets |
| [Fig. 17](results/figures/fig17_convergence.png) | Loads and loop errors against step size and cycles marched |
| [Fig. 18](results/figures/fig18_state_space.png) | The state-space form against the indicial march |

The tables are in [`results/tables.md`](results/tables.md): the measured loops
used, static inputs, constants, the attached-flow moment, accuracy against the
targets, accuracy by stall depth, the two separation laws, cycle damping, the
attached-flow check, the second aerofoil, every target against its measured
value, numerical convergence and the post-stall curve.

---

## Findings

{_equivalence_paragraph(same, eq, n_prim, tg)}

**2. What separate calibration shows.** {_own_calibration_sentence(own, same)}
(`results/comparison.csv`)

**3. Accuracy on held-out loops.** {n_prim} held-out NACA 0012 loops, each
run at its own measured conditions (`results/validation_targets.csv`). Loop
errors are RMS differences divided by the measured range of the quantity over
the loop. The targets are in `data/targets.json` and were set before the
held-out loops were scored with this form of the model
([`docs/provenance.md`](docs/provenance.md) gives the dated record). Each
limit is the stricter of an absolute limit and three quarters of what the
earlier form of the model had scored on the same loops:
{_target_origin(T)}. The absolute limits are the author's working limits, not
a standard, and the held-out loops had been looked at with the earlier form,
so they are a consulted test set and not an untouched one. They are never
used in a fit.
`results/targets_scoreboard.csv` lists every target with its measured value:
{_scoreboard_counts(board)}

| Measure | Measured | Target | Met |
|---|---|---|---|
{_targets_rows(tg)}

The error a perfect model would show from the stated measurement uncertainty
alone is {unc.measurement_floor['nRMS_CL']:.3f} in lift,
{unc.measurement_floor['nRMS_CM']:.3f} in moment and
{unc.measurement_floor['nRMS_CD']:.3f} in drag (`results/uncertainty_summary.csv`);
the errors above are {tg.measured['mean_nRMS_CL']/unc.measurement_floor['nRMS_CL']:.1f},
{tg.measured['mean_nRMS_CM']/unc.measurement_floor['nRMS_CM']:.1f} and
{tg.measured['mean_nRMS_CD']/unc.measurement_floor['nRMS_CD']:.1f} times that. The
figure for drag rests on an assumed uncertainty of the pressure drag, for
which the report gives no number, and all three treat the quoted uncertainty
as a random error of each point; if it is a bias instead, it does not average
out and the same figures are its size.

The incidence error of moment stall is defined only where both the
measurement and the model show a moment stall on the up-stroke:
{int(stall.in_both)} of the {int(stall.n_frames)} loops. The model misses a
measured moment stall on {int(stall.measured_only_missed_by_model)} loops and
predicts one that is not measured on {int(stall.model_only_not_measured)}
(`results/validation_moment_stall.csv`).

All {cal['n_frames']} calibration loops are at Mach 0.285 or above.
{int(hi_m.n_frames['mean_nRMS_CL'])} of the held-out loops are too (lift
{hi_m.value['mean_nRMS_CL']:.3f}, moment {hi_m.value['mean_nRMS_CM']:.3f},
incidence of maximum lift {hi_m.value['mean_abs_dalpha_CLmax_deg']:.2f}°); the
{int(lo_m.n_frames['mean_nRMS_CL'])} below it, at Mach
{prim[prim.M < 0.285].M.min():.2f} to {prim[prim.M < 0.285].M.max():.2f}, give
lift {lo_m.value['mean_nRMS_CL']:.3f}, moment {lo_m.value['mean_nRMS_CM']:.3f}
and {lo_m.value['mean_abs_dalpha_CLmax_deg']:.2f}°. Leaving out the
{n_prim - int(fresh.n_frames['mean_nRMS_CL'])} held-out loops whose conditions
repeat a calibration loop changes the lift and moment errors to
{fresh.value['mean_nRMS_CL']:.3f} and {fresh.value['mean_nRMS_CM']:.3f}.

**4. Where the error sits.** The same held-out loops grouped by how far the
peak incidence goes beyond static stall (`results/validation_by_stall_depth.csv`):

| Group | Loops | Lift | Moment | Drag | Damping sign correct |
|---|---|---|---|---|---|
{_depth_rows(depth)}

**5. Cycle damping is not predicted through stall**
(`results/validation_damping.csv`, `results/attached_moment_factor.json`).

{_attached_moment_paragraph(att_m, unc)}
- *Stalled loops.* On the {int(d_all.n_frames)} held-out loops with a closed
  measured moment loop the sign of the damping is predicted correctly on
  {d_all.sign_agreement:.0%}; a prediction of positive damping everywhere
  would score {d_all.sign_agreement_if_always_positive:.0%}. The correlation
  between predicted and measured damping is {d_all.correlation:.2f} and the
  RMS error {d_all.rms_error:.3f}. On the {int(d_deep.n_frames)} loops that
  pass static stall by 6° or more the sign is correct on
  {d_deep.sign_agreement:.0%}.
- *How certain the measured sign is.* The stated uncertainty of the measured
  moment gives the measured damping an uncertainty of
  {d_all.mean_uncertainty_random:.3f} if the error is random from point to
  point and up to {d_all.mean_uncertainty_stroke:.3f} if it follows the
  stroke. The measured sign is certain in the first sense on
  {int(d_rand.n_frames)} loops (model sign correct on
  {d_rand.sign_agreement:.0%}, always-positive
  {d_rand.sign_agreement_if_always_positive:.0%}) and in the second on
  {int(d_sure.n_frames)} (model {d_sure.sign_agreement:.0%}, always-positive
  {d_sure.sign_agreement_if_always_positive:.0%}).
- *Negative damping.* {int(d_neg.n_frames)} held-out loops measure negative
  damping; the model predicts negative damping on
  {int(d_neg.predicted_negative)} of them, and on
  {int(d_pos.predicted_negative)} loops whose measured damping is positive.
  **The model does not predict cycle damping through stall with useful
  accuracy. It cannot locate a stall-flutter boundary and is not offered for
  that purpose.**

**6. The model is provided in indicial and state-space form.** The indicial
march of `unistall/dsmodel.py` and the same equations written as
{int(ssf.states)} first-order differential equations and integrated by the
classical fourth-order Runge–Kutta method (`unistall/statespace.py`) give the
same loads to the differences that follow: on the {int(ssf.n_frames)} loops of the calibration set
({int(ssf.n_frames) - int(ssf.n_frames_below_mach_0_20)} inside the Mach range of the static data, which are the ones
used in the fit, and {int(ssf.n_frames_below_mach_0_20)} below Mach 0.20, used here only to exercise the low-Mach blend) the largest
difference in any loop error is {ssf.worst_difference_in_loop_error:.4f}
(target {ssf.target_difference_in_loop_error}), in cycle damping
{ssf.worst_difference_in_damping:.4f}, and at any instant
{ssf.worst_difference_in_CL_pct_of_range:.1f} % of the range in lift and
{ssf.worst_difference_in_CM_pct_of_range:.1f} % in moment; those largest
instantaneous differences sit where the model switches a time constant, which
the two forms resolve to one step (`results/statespace_summary.csv`, Fig. 18).
Coupled to a torsional section, the pitch angle, its rate and the aerodynamic
states are integrated as one system and give the same growth or decay as the
step-by-step coupling in all {int(struct_ok)} of {int(struct_n)} conditions
(`results/structural_checks.csv`). The state-space form covers the same Mach
range as the indicial one: {int(ssf.n_frames_below_mach_0_20)} of those loops
are below Mach 0.20, down to Mach {ssf.lowest_mach:.2f}, where the
incompressible loads are blended in.

---

## Limits of the model

Targets not met on the held-out loops:

{_not_met(tg)}

The calibration lowers the cost by {gain:.0f} % from the literature
constants. {_refined_starts_sentence(cal)}

{_limit_checks_paragraph()}

{_published_reference_paragraph()}

The stall-onset level is taken from the static-stall normal force and is not
fitted. Multiplying it by {_onset_sentence()}
(`results/onset_level_check.csv`). That check moves the level with the other
constants held, on the calibration cost only, so it does not show that the
misses above are independent of where onset is placed.{_candidates_paragraph()}

Other limits:

- **Range.** Static data are at Mach {lo.M:.3f} and {hi.M:.3f}. Of the {n_prim}
  held-out loops, {int((prim.M >= MACH_CAL).sum())} are at Mach {prim.M[prim.M >= MACH_CAL].min():.3f} to
  {prim.M.max():.3f}, where every calibration loop also lies; the other
  {int((prim.M < MACH_CAL).sum())}, at Mach {prim.M.min():.3f} to {prim.M[prim.M < MACH_CAL].max():.3f}, are probes of
  the interpolation in Mach number and are predicted worse (finding 3), so
  the comparison is in effect at Mach {prim.M[prim.M >= MACH_CAL].min():.2f} to {prim.M.max():.2f}. Reduced frequency runs up
  to {prim.k.max():.3f} and peak incidence up to {(prim.alpha0_deg + prim.amp_deg).max():.0f}° on the NACA 0012
  section, in sinusoidal pitch only. The experiment stops at Mach
  0.30. The model reports any condition outside that range and refuses it
  when called with `strict=True`. Nothing is claimed outside it.
- **Numerical resolution.** There is no spatial mesh; the discretisation is
  the time step, set to at most {dm_ds:.2f} semichords of travel. Against a march
  with four times as many steps the loop error of any calibration loop
  changes by at most {conv.worst_change_in_loop_error['step size']:.4f}, and
  against four more cycles by
  {conv.worst_change_in_loop_error['cycles marched']:.4f}
  (`results/convergence_summary.csv`, Fig. 17).
- **Static inputs.** Past static stall the measurements have two states; the
  model follows the more separated one, which gives a calibration cost of
  {curve['more separated']:.2f} against {curve['less separated']:.2f}
  (`results/static_curve_choice.csv`). The static data above 17° at the lower
  Mach station come from a sweep at Mach 0.204.
- **Low Mach number.** The compressible indicial constants alone differ from
  Theodorsen's solution at Mach 0.05 by {alone.lift_amp_err_pct.abs().max():.1f} %
  in lift amplitude and {alone.moment_amp_err_pct.abs().max():.1f} % in moment
  amplitude. The model blends them with the incompressible loads below Mach
  0.20 and then agrees to {model_th.lift_amp_err_pct.abs().max():.1f} % and
  {model_th.moment_amp_err_pct.abs().max():.1f} %
  (`results/attached_checks.csv`). No dynamic-stall loop below Mach 0.20 is used
  for a claim.
- **Drag.** Over every scored loop inside the Mach range, calibration and
  held-out, drag is negative for at most {100*t_in.CD_neg_frac.max():.1f} % of a
  cycle (limit {100*S["CD_negative_fraction_of_cycle_max"]:.0f} %, on
  frame {t_in.loc[t_in.CD_neg_frac.idxmax(), "frame"].split("_")[1]}) and its
  lowest value is {t_in.CD_min_model.min():+.4f} (limit {S["CD_min"]:+.2f});
  {"both limits are met on every loop" if drag_all_ok else "a limit is exceeded on at least one loop"}
  (`results/validation_table.csv`, `results/separated_checks.csv`).
- **Other sections.** On a second aerofoil (Ames A-01,
  {own_static.frame.nunique()} loops) with the dynamic constants unchanged,
  that aerofoil's own static data against the NACA 0012 static data give:
  {_other_airfoil_sentence(own_static, naca_static)}. Transfer to another
  section without re-calibration is not shown.
- **Scope.** This is an aerofoil model in a steady stream with prescribed
  pitch. It is not a rotor calculation and no rotor data are compared with.
  The dynamic measurements are not corrected for wall interference.

---

## Case study

The numbered folders hold a case study built on the solver: the section at a
tunnel condition (Case A) and at a rotor-blade-station condition (Case B),
with its geometry, a grid, the solution files, reconstructed flow fields,
plots, a report and illustrative drawings. `python3 run_all.py` reproduces it.
[`00_overview/case_definition.md`](00_overview/case_definition.md) states what
each part is and is not.

{_case_rows()}

| Folder | Contents |
|---|---|
| [`00_overview/`](00_overview) | The case definition, and a list of what every numbered folder holds and which script writes it |
| [`01_geometry/`](01_geometry) | Section coordinates and thickness |
| [`02_mesh/`](02_mesh) | An O-grid round the section; the solver does not read it, the reconstructed field is evaluated at its nodes |
| [`03_model_setup/`](03_model_setup) | Flow conditions, kinematics, air properties, solver configuration, static inputs |
| [`04_solver/`](04_solver) | The driver for the two cases and a command-line entry for a single condition (the solver itself is the `unistall` package) |
| [`05_solution/`](05_solution) | Time histories, surface pressures, reconstructed fields, the model's slow-sweep polar, metrics, convergence |
| [`06_postprocessing/`](06_postprocessing) | Hysteresis loops, contour maps of the reconstructed field (pressure, Mach number, temperatures, speed, vorticity), surface plots, and `validation/`, the comparison with measurement set out from the result files |
| [`07_report/`](07_report) | The case-study report, a plots album and a data dossier |
| [`08_engineering_drawings/`](08_engineering_drawings) | Four illustrative sheets, every dimension from one sourced table |

**The flow fields are a reconstruction, not a flow solution**: a potential
flow drawn round the lift the model predicts, with an outflow through the
suction surface standing for the separated region and chosen so that the flow
leaves the trailing edge smoothly. The field is incompressible; its pressure
is corrected for compressibility by the Kármán–Tsien rule, and local Mach
number and temperatures follow from the corrected pressure by the isentropic
relations. That correction holds only where the flow is subcritical: where
the corrected pressure is below the critical value the points are flagged,
hatched in every figure and counted, and are not physical. There is no
boundary layer or separated shear layer. It is given the lift only: {_reconstruction_sentence()} None of it
is compared with a measurement and no load depends on it. **Case B is not a
rotor calculation**: it is the section in a steady stream at a blade station's
Mach number and reduced frequency. The evidence for the model's accuracy is in
the findings above, not in the case study.

![Case A hysteresis loops](06_postprocessing/plots/hyst_cl_A_validation.png)

The measured loop drawn with Case A is frame 9302, a calibration loop: the
constants were fitted with it, so this figure shows a fit and not a prediction.

![Reconstructed pressure field, Case A at peak lift]({_first_plot("contour_Cp_A_validation_peak")})

---

## The model in brief

- **Attached flow.** Circulatory lift from two deficiency functions on the
  three-quarter-chord incidence; impulsive loads for incidence and pitch rate
  with the time constant c/a; unsteady moment terms for both.
- **Separation.** Pressure lag on the normal force, boundary-layer lag on the
  separation point, Kirchhoff's relation for the separated normal force.
- **Dynamic stall.** A leading-edge vortex is shed when the lagged normal
  force exceeds the static-stall level; its lift accumulates while it is over
  the chord and decays afterwards; repeated shedding is allowed.
- **Static inputs vary with Mach number.** Maximum static lift in this
  experiment is {lo.CL_max:.2f} at Mach {lo.M:.3f} and {hi.CL_max:.2f} at Mach
  {hi.M:.3f}, so slope, stall level and separation point are interpolated
  between two Mach stations.
- **Constants.** Fitted on {cal['n_frames']} calibration loops: {fitted}.
  Returned to literature values because the data do not identify them:
  {dropped}. {_identifiability_sentence(cal)}

Every equation and constant is in [`docs/formulation.md`](docs/formulation.md).

---

## Data

- {n_verified} of {len(man)} measured loops have their Mach number, reduced
  frequency, mean angle, amplitude and aerofoil checked against the tables of
  NASA TM-84245 Volume 1 (`data/data_manifest.csv`).
- Static lift, moment and drag at Mach 0.30 are digitised from the report's
  own figures and tables, each point with its source page and uncertainty
  (`data/static_naca0012_M030_*.csv`).
{_split_sentence(T)}
- **The measured loops are third-party data and are not distributed here.**
  `python3 -m unistall.fetch_frames` downloads them from one fixed commit of
  the source repository and checks each file against the SHA-256 recorded in
  `data/frames_inventory.csv`; it names any file that is missing or has
  changed. If the source is ever withdrawn or altered, the results cannot be
  reproduced from it: the same loops are printed in Volume 2 of the NASA
  report and would have to be digitised again, and the recorded checksums
  show whether a copy obtained elsewhere is the one used here.
- `data/data_sources.csv` lists each source with the tables, figures and
  pages that the data files use, and `python3 -m unistall.build_manifest
  --check` rebuilds the list of usable loops from the report's tables.

---

## Reproducing the results

**A quick check.** `python3 -m pytest -q` runs in minutes. It recomputes the
held-out measures from the stored constants and the fetched loops and
requires them to equal the recorded ones, re-derives the data manifest and
the source log, checks the model against theory and against its own
invariances, and fails if a table, figure or document is not what its script
writes. That is the path for a reader who wants to confirm the results
without repeating the calibration.

**The whole chain** takes many hours, most of it the calibration of the two
models and its cross-validation; the commands below are in the order they
are run. The state-space form is slower than the indicial march by the
factor given in `results/statespace_summary.csv` and is used for the checks
only; every reported result comes from the indicial march.

```bash
pip install -r requirements-lock.txt     # run everything from the repository root
pip install -e .                         # or `pip install .` to use the model outside the clone

python3 -m unistall.fetch_frames              # measured loops (needs network)
python3 -m unistall.static_stations           # slope, zero-lift incidence and stall level from the two static sweeps
python3 -m unistall.static_polar              # static polar at Mach 0.30 from the digitised points
python3 -m unistall.static_model              # static inputs by Mach number
python3 -m unistall.check_attached            # attached-flow checks against Theodorsen
python3 -m unistall.check_indicial            # compressible branch against exact starting values and short-time lift
python3 -m unistall.calibrate --attached      # unsteady-moment factor from loops below static stall
python3 -m unistall.calibrate --curve         # cost with each post-stall static curve
python3 -m unistall.calibrate --fit && python3 -m unistall.calibrate --sensitivity && python3 -m unistall.calibrate --select
python3 -m unistall.calibrate --fit --model=reference && python3 -m unistall.calibrate --sensitivity --model=reference && python3 -m unistall.calibrate --select --model=reference
python3 -m unistall.calibrate --repeat && python3 -m unistall.calibrate --repeat --model=reference
python3 -m unistall.calibrate --onset         # cost against the stall-onset level
python3 -m unistall.convergence               # step size and cycles marched
python3 -m unistall.check_statespace          # the state-space form against the indicial march
python3 -m unistall.check_separated           # drag sign, onset step-independence
python3 -m unistall.check_limits              # mirror pair, static-limit drag, zero lift, switching rules
python3 -m unistall.published_reference       # beside the published Leishman-Beddoes curves for two cases
python3 -m unistall.validate                  # held-out scores, comparisons, second aerofoil
python3 -m unistall.uncertainty               # bands and measurement floor
python3 -m unistall.trend_sweep               # trends inside the range the data cover
python3 -m unistall.structural                # torsional section, dimensional loads
python3 -m unistall.code_quality              # regression and code measures
python3 -m unistall.scoreboard                # every target against its measured value
python3 -m unistall.make_figures --refresh    # every figure in results/figures/
python3 -m unistall.make_tables               # results/tables.md
python3 -m unistall.make_readme               # this file

python3 -m pytest                           # {len(list((ROOT/'tests').glob('test_*.py')))} test files
```

Using the model directly:

```python
from unistall import dsmodel
out = dsmodel.solve(alpha_mean_deg=10.0, alpha_amp_deg=8.0, k=0.10, M=0.30, strict=True)
out["alpha_deg"], out["CL"], out["CD"], out["CM"]      # the last cycle
```

---

## Repository structure

| File or folder | Contents |
|---|---|
| `unistall/dsmodel.py` | The solver: the dynamic-stall load model, indicial form |
| `unistall/statespace.py` | The same model in state-space form |
| `unistall/attached_flow.py` | Attached-flow indicial loads, the incompressible limit, Theodorsen's solution |
| `unistall/static_model.py`, `static_polar.py` | Static inputs by Mach number from the measured static data |
| `unistall/reference_lb.py` | The reference model: the same equations with the fitted exponential separation law |
| `unistall/calibrate.py` | Scripted calibration: attached-flow moment factor, fit, sensitivity, selection by cross-validation |
| `unistall/metrics.py`, `validate.py`, `uncertainty.py`, `acceptance.py` | Error measures, held-out scores, comparisons, uncertainty, acceptance thresholds |
| `unistall/check_attached.py`, `check_separated.py`, `check_limits.py`, `structural.py`, `trend_sweep.py` | Checks against theory, drag sign and step independence, behaviour at the limits, a torsional section, trends |
| `unistall/make_figures.py`, `make_tables.py`, `make_readme.py` | Every figure, the tables and this file, from the results |
| `unistall/fetch_frames.py` | Downloads the measured loops and checks them against `data/frames_inventory.csv` |
| `data/` | Inputs: `data_manifest.csv` (verified conditions), `split.json` (calibration and held-out sets), `targets.json` (every target), digitised static data |
| `results/` | Result files and `tables.md`; `results/figures/` holds the figures |
| `docs/` | `formulation.md` (every equation and constant), `provenance.md` (when the targets and the split were set), `references.bib` |
| `tests/` | Test suite |
| `00_overview/` to `08_engineering_drawings/` | The case study (see above) |
| `unistall/flowfield.py` | The flow-field reconstruction used by the case study |
| `run_all.py`, `check_case_study.py`, `project_meta.py` | One-command run of the case study, checks on its outputs, its identity |
| `index.md`, `_config.yml`, `_layouts/`, `_includes/`, `assets/` | The project site |
| `pyproject.toml`, `requirements.txt`, `requirements-lock.txt` | Package definition (`pip install -e .`), dependency bounds, exact versions |
| `.github/` | Continuous-integration workflow |
| `CITATION.cff` | How to cite the work |
| `LICENSE`, `NOTICE` | Copyright terms and the notice on third-party data |

Code measures (`results/code_quality.csv`): {int(quality['lint_findings_whole_repository'])} lint
findings, {int(quality['unused_parameters_in_package'])} unused parameters,
{quality['public_package_functions_typed_and_documented_pct']:.0f} % of public functions typed
and documented, longest function
{int(quality['longest_package_function_lines_excluding_docstring'])} lines.

---

## References

`docs/references.bib` holds {len(refs)} references:
{int((refs.registry != 'other').sum())} with an identifier resolved through
Crossref or the NASA Technical Reports Server and
{int((refs.registry == 'other').sum())} entered by hand from the record named
for each in `docs/bibliography_resolution.csv`, being works those registries
do not hold.

**Prior work.** Reading the separation point from the static data instead of
fitting it is established practice and is not proposed here: Hansen, Gaunaa &
Madsen (2004), *A Beddoes-Leishman Type Dynamic Stall Model in State-Space and
Indicial Formulations*, Risø-R-1354(EN), obtain it by inverting Kirchhoff's
relation against the static lift curve, in state-space and indicial form, and
Larsen, Nielsen & Krenk (2007) and Damiani & Hayman (2019) do likewise. What
this repository adds is a measurement of how much that choice changes the
predicted loads when everything else is held fixed, and an assessment on
held-out loops.

The formulation follows Damiani & Hayman (2019), *The Unsteady Aerodynamics
Module for FAST 8*, NREL/TP-5000-66347, doi:10.2172/1576488, and Leishman &
Beddoes (1989), *A Semi-Empirical Model for Dynamic Stall*, J. Am. Helicopter
Soc. 34(3), doi:10.4050/jahs.34.3.3.

---

## Copyright and citation

© 2026 Akosa Samuel Onyejekwe. All rights reserved.

This work is published so that it can be read and its results checked. You
may view it, run the unmodified code to verify the reported results, and cite
it. Copying, modifying, redistributing or reusing any part of it, and any
commercial use, need the author's prior written permission. The full terms are
in [`LICENSE`](LICENSE).

**Third-party data are not distributed here.** The measured loops are not the
author's. The original measurements are NASA TM-84245 (a work of the U.S.
Government); the digitised files come from
[BL-DSM-JFS-2021](https://github.com/luizpancini/BL-DSM-JFS-2021), which states
no licence of its own. Details in [`NOTICE`](NOTICE).

To cite this work, use [`CITATION.cff`](CITATION.cff).
"""


if __name__ == "__main__":
    text = build()
    if "--print" in sys.argv:
        print(text)
    else:
        (ROOT/"README.md").write_text(text, encoding="utf-8")
        print(f"[readme] README.md written ({len(text.splitlines())} lines)")
