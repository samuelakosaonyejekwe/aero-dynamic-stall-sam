# -*- coding: utf-8 -*-
"""
unistall / scoreboard.py
--------------------
Every target in data/targets.json against the value measured for it, met or
not, in one table: results/targets_scoreboard.csv. A target with no measured
value is listed as "not evaluated", so none can be passed over. A target is
met or not met; nothing else. Each is marked physical (the model against
measurement or theory), numerical (the method) or software (the code).

Author: Akosa Samuel Onyejekwe (independent)

Usage:  python3 -m unistall.scoreboard
"""
import json
import pandas as pd

from unistall.paths import DATA, RESULTS

SKIP = ("theodorsen_M", "theodorsen_k", "group", "note", "lint_rule_set", "second_pitch_axis_x_c",
        "sensitivity_perturbation_pct", "drop_constant_if_cost_change_below_pct", "state_space_note",
        "selection_rule", "selection_rule_note", "origin", "added_on", "margin_fraction_of_accuracy_target", "quantity", "test", "worst_single_loop")

# what a target is a target on: the model against measurement or theory, the
# numerical method, or the software. Counts are reported by kind.
KIND_OF_GROUP = {"attached_flow": "physical", "separated_flow": "physical", "held_out_validation": "physical",
                 "equivalence_of_separation_laws": "physical", "structural": "numerical", "calibration": "numerical",
                 "numerical": "numerical", "code_quality": "software"}


def _measured() -> dict:
    """target name -> (measured value, where it comes from)."""
    att = pd.read_csv(RESULTS/"attached_checks.csv")
    th = att[(att.check == "Theodorsen") & att.variant.str.startswith("the model")]
    ns = att[att.check == "no-stall frame, the model"]
    ind = att[att.check == "chord and speed independence"]
    sep = pd.read_csv(RESULTS/"separated_checks.csv")
    tg_all = pd.read_csv(RESULTS/"validation_targets.csv").set_index("measure")
    tg = tg_all.measured
    st = pd.read_csv(RESULTS/"structural_checks.csv")
    cq = pd.read_csv(RESULTS/"code_quality.csv").set_index("measure")["value"]
    conv = pd.read_csv(RESULTS/"convergence_summary.csv")
    rep = json.load(open(RESULTS/"calibration_repeatability.json"))
    out = {
        "lift_and_moment_amplitude_error_max_pct": (max(th.lift_amp_err_pct.abs().max(), th.moment_amp_err_pct.abs().max()), "attached_checks.csv"),
        "lift_and_moment_phase_error_max_deg": (max(th.lift_phase_err_deg.abs().max(), th.moment_phase_err_deg.abs().max()), "attached_checks.csv"),
        "chord_speed_independence_max_abs_dCL": (float(ind.max_abs_difference.max()), "attached_checks.csv"),
        "no_stall_frame_nRMS_CL_max": (float(ns.nRMS_CL.max()), "attached_checks.csv"),
        "no_stall_frame_nRMS_CM_max": (float(ns.nRMS_CM.max()), "attached_checks.csv"),
        "CD_negative_fraction_of_cycle_max": (float(sep.CD_neg_frac.max()), "separated_checks.csv"),
        "CD_min": (float(sep.CD_min.min()), "separated_checks.csv"),
        "onset_shift_with_4x_steps_max_deg": (float(sep.onset_shift_deg.max()), "separated_checks.csv"),
        "repeatability_max_abs_difference": (float(rep["max_abs_difference_in_constants"]), "calibration_repeatability.json"),
        "mean_nRMS_CL_max": (tg["mean_nRMS_CL"], "validation_targets.csv"),
        "mean_nRMS_CM_max": (tg["mean_nRMS_CM"], "validation_targets.csv"),
        "mean_nRMS_CD_max": (tg["mean_nRMS_CD"], "validation_targets.csv"),
        "mean_abs_dalpha_CLmax_deg_max": (tg["mean_abs_dalpha_CLmax_deg"], "validation_targets.csv"),
        "mean_abs_dalpha_Mstall_deg_max": (tg["mean_abs_dalpha_Mstall_deg"], "validation_targets.csv: "
                                           + str(tg_all.get("note", pd.Series(dtype=str)).get("mean_abs_dalpha_Mstall_deg", ""))),
        "damping_sign_agreement_min": (tg["damping_sign_agreement"], "validation_targets.csv"),
        "cycle_work_identity_max_error_pct": (float(st.work_identity_error_pct.abs().max()), "structural_checks.csv"),
        "free_response_sign_agreement_min_conditions": (int(st.agree.sum()), "structural_checks.csv"),
        "dimensional_load_identity_max_abs": (float(st.dimensional_identity_max_abs.max()), "structural_checks.csv"),
        "regression_max_abs_difference": (cq["regression_max_abs_difference"], "code_quality.csv"),
        "import_path_edits_load_model_path_max": (cq["import_path_edits_load_model_path"], "code_quality.csv"),
        "module_level_executable_lines_per_load_model_file_max": (cq["module_level_executable_lines_worst_load_model_file"], "code_quality.csv"),
        "longest_package_function_lines_max": (cq["longest_package_function_lines_excluding_docstring"], "code_quality.csv"),
        "unused_parameters_max": (cq["unused_parameters_in_package"], "code_quality.csv"),
        "public_functions_typed_and_unit_documented_pct_min": (cq["public_package_functions_typed_and_documented_pct"], "code_quality.csv"),
        "lint_findings_max": (cq["lint_findings_whole_repository"], "code_quality.csv"),
        "silently_skipped_checks_max": (cq["silently_skipped_checks"], "code_quality.csv"),
        "revision_history_phrases_in_package_max": (cq["revision_history_phrases_in_package"], "code_quality.csv"),
        "loop_error_change_max": (float(max(conv.worst_error_against_zero_step.max(),
                                            conv.worst_change_in_loop_error[conv.study != "step size"].max())),
                                  "convergence_summary.csv"),
        "state_space_loop_error_difference_max": (
            float(pd.read_csv(RESULTS/"statespace_summary.csv").worst_difference_in_loop_error.iloc[0]),
            "statespace_summary.csv"),
    }
    eq = pd.read_csv(RESULTS/"equivalence.csv")
    for q in ("CL", "CM", "CD"):
        out[f"margin_{q}"] = (float(eq[eq.coefficient == q].ci95_hi.max()), "equivalence.csv")
    th2 = att[(att.check == "Theodorsen second axis")] if "Theodorsen second axis" in set(att.check) else None
    if th2 is not None and len(th2):
        out["second_axis_theodorsen_amplitude_max_pct"] = (max(th2.lift_amp_err_pct.abs().max(), th2.moment_amp_err_pct.abs().max()), "attached_checks.csv")
        out["second_axis_theodorsen_phase_max_deg"] = (max(th2.lift_phase_err_deg.abs().max(), th2.moment_phase_err_deg.abs().max()), "attached_checks.csv")
    return out


def _uncertainty_levels() -> dict:
    """For the held-out loop errors, the error a perfect model would show from
    the stated measurement uncertainty alone (results/uncertainty_summary.csv).
    Given for scale; it plays no part in whether a target is met."""
    unc = pd.read_csv(RESULTS/"uncertainty_summary.csv").set_index("metric").measurement_floor
    return {"mean_nRMS_CL_max": float(unc["nRMS_CL"]), "mean_nRMS_CM_max": float(unc["nRMS_CM"]),
            "mean_nRMS_CD_max": float(unc["nRMS_CD"])}


def build() -> pd.DataFrame:
    """One row per target: its value, the measured value, and whether it is met."""
    T = json.load(open(DATA/"targets.json"))
    measured = _measured()
    levels = _uncertainty_levels()
    rows = []
    for group, d in T.items():
        if not isinstance(d, dict):
            continue
        for name, target in d.items():
            if name in SKIP:
                continue
            at_least = name.endswith("_min") or name == "CD_min" or "min_conditions" in name
            got = measured.get(name)
            if got is None:
                rows.append(dict(group=group, kind_of_target=KIND_OF_GROUP[group], target=name,
                                 kind="at least" if at_least else "at most", value=target,
                                 measured=None, measurement_floor=None, met="not evaluated", source=""))
                continue
            v = float(got[0])
            met = v >= target if at_least else v <= target
            floor = levels.get(name)
            rows.append(dict(group=group, kind_of_target=KIND_OF_GROUP[group], target=name,
                             kind="at least" if at_least else "at most", value=target,
                             measured=float(f"{v:.4g}"), measurement_floor=None if floor is None else float(f"{floor:.4g}"),
                             met="yes" if met else "no", source=got[1]))
    return pd.DataFrame(rows)


if __name__ == "__main__":
    df = build()
    df.to_csv(RESULTS/"targets_scoreboard.csv", index=False)
    print(df.to_string(index=False))
    for kind, d in df.groupby("kind_of_target", sort=False):
        print(f"[scoreboard] {kind}: {int((d.met == 'yes').sum())} met, {int((d.met == 'no').sum())} not met, "
              f"{int((d.met == 'not evaluated').sum())} not evaluated")
