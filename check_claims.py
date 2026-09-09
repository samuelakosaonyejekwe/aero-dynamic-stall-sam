"""
check_claims.py — assert that every number quoted in README.md and
00_overview/case_definition.md is still the number the pipeline produced.

Why this exists
---------------
The README states roughly two dozen results: Mach numbers, C_L,max, the
moment break, stall onset, the damping band, the held-out validation errors.
Markdown cannot compute, so each of those is a transcribed copy of a value that
lives in a generated CSV. Copies drift. This script is the thing that stops them
drifting: it re-reads the artifacts, formats each claim exactly as the prose
states it, and fails if the prose no longer contains it.

Run it as the last pipeline stage (run_all.py does), or on its own:

    python3 check_claims.py

Exit status 0 = every claim traces to an artifact; 1 = at least one has drifted,
and the offending claim is printed with the value that should replace it.
"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from project_meta import STUDY_DATE           # the one date of record


def _series(rel, key, val):
    return pd.read_csv(ROOT / rel).set_index(key)[val]


flow = pd.read_csv(ROOT / "03_model_setup/flow_conditions.csv").set_index("parameter")
A, B = flow["case_A_validation"], flow["case_B_application"]
mA = _series("05_solution/metrics_A_validation.csv", "metric", "value")
mB = _series("05_solution/metrics_B_application.csv", "metric", "value")
vs = _series("06_postprocessing/validation/validation_realdata_summary.csv",
             "metric", "value")
vst = pd.read_csv(ROOT / "06_postprocessing/validation/validation_static.csv")
nr = pd.read_csv(ROOT / "06_postprocessing/validation/validation_nasa_real.csv"
                 ).set_index("frame")
_x = nr.loc["frame_25104"]
xcheck_pct = 100.0*abs(_x.CLmax_model - _x.CLmax_exp)/_x.CLmax_exp

kin = pd.read_csv(ROOT / "03_model_setup/kinematics.csv").set_index("case_id")
kA, kB = kin.loc["A_validation_rig"], kin.loc["B_application_rotor"]

f = float
sl = vst[vst.metric.str.startswith("lift-curve")].iloc[0]
cm = vst[vst.metric.str.startswith("CL_max")].iloc[0]


def _row(flow_col, kin_row, dp, k_dp):
    """The trailing cells of one row of the README's "Configurations solved"
    table. Case B's chord, Mach and k were already guarded one by one while
    Case A's were not guarded at all and neither case's incidence was, so half
    the table could drift. Guarding the rendered row covers all four at once.

    dp is the decimal places the table gives the chord and the Mach (the two
    happen to share one in both rows: 2 for the rig, 3 for the blade station);
    k_dp is the k column's."""
    return (f"| {f(flow_col['chord_c']):.{dp}f} m "
            f"| {f(flow_col['freestream_mach_M']):.{dp}f} "
            f"| {f(kin_row['reduced_freq_k']):.{k_dp}f} "
            f"| {f(kin_row['alpha_mean_deg']):.0f}° ± {f(kin_row['alpha_amp_deg']):.0f}° |")

# (document, exact substring that must be present, what it came from)
CLAIMS = [
    # ---- identity ---------------------------------------------------------
    # The date of record lives in project_meta.py and is rendered onto every
    # report cover and every drawing title block from there; the README states
    # it in prose, which is a copy, so it is guarded like every other copy.
    ("README.md", f"**Date:** {STUDY_DATE}", "project_meta.STUDY_DATE"),
    # ---- case conditions --------------------------------------------------
    ("README.md", _row(A, kA, 2, 2), "flow_conditions + kinematics, case A table row"),
    ("README.md", _row(B, kB, 3, 3), "flow_conditions + kinematics, case B table row"),
    ("README.md", f"analysis station r/R = {f(B['radial_station_r_R']):.2f}",
     "flow_conditions radial_station_r_R"),
    ("README.md", f"{f(B['chord_c']):.3f} m", "flow_conditions chord_c (B)"),
    ("README.md", f"{f(B['freestream_mach_M']):.3f}", "flow_conditions M (B)"),
    ("README.md", f"{f(B['reduced_frequency_k']):.3f}", "flow_conditions k (B)"),
    ("README.md", f"R = {f(B['rotor_radius_R']):.2f} m", "flow_conditions rotor_radius_R"),
    ("README.md", f"{f(B['rotor_speed_Omega']):.1f} rad/s", "flow_conditions rotor_speed_Omega"),
    ("README.md", f"ΩR = {f(B['rotor_tip_speed_OmegaR']):.1f} m/s", "flow_conditions OmegaR"),
    ("README.md", f"μ = {f(B['advance_ratio_mu']):.2f}", "flow_conditions advance_ratio_mu"),
    ("README.md", f"= {f(B['freestream_velocity_U']):.2f} m/s", "flow_conditions U (B)"),
    ("README.md", f"Mach is {f(B['freestream_mach_M']):.4f}", "flow_conditions M (B), 4 dp"),
    # The near-trailing-edge panel oscillation and the leading-edge peak it is
    # measured against: both are metrics, so both are guarded like every other
    # transcribed number.
    ("README.md", f"{f(mA['Cp_TE_panel_oscillation_max_abs']):.3f}",
     "metrics_A Cp_TE_panel_oscillation_max_abs"),
    ("README.md", f"{f(mA['Cp_max_abs_outside_TE_zone']):.3f}",
     "metrics_A Cp_max_abs_outside_TE_zone"),
    ("README.md", f"x/c > {f(mA['Cp_TE_panel_oscillation_zone_x_c']):.3f}",
     "metrics_A Cp_TE_panel_oscillation_zone_x_c"),
    # ---- headline results -------------------------------------------------
    ("README.md", f"C_L,max = {mA['CL_max_dynamic']} at α = {f(mA['alpha_at_CLmax_deg']):.1f}°",
     "metrics_A CL_max_dynamic / alpha_at_CLmax_deg"),
    ("README.md", f"C_M,c/4 break = {mA['CM_min(c/4)']}, C_D,max = {mA['CD_max']}",
     "metrics_A CM_min / CD_max"),
    ("README.md", f"onset at α = {f(mA['stall_onset_alpha_deg']):.2f}°",
     "metrics_A stall_onset_alpha_deg"),
    ("README.md", f"a {100*(f(mA['dynamic_overshoot_ratio'])-1):.0f} % overshoot",
     "metrics_A dynamic_overshoot_ratio"),
    ("README.md", f"C_L,max = {mB['CL_max_dynamic']} at α = {f(mB['alpha_at_CLmax_deg']):.1f}°",
     "metrics_B CL_max_dynamic / alpha_at_CLmax_deg"),
    ("README.md", f"C_M,c/4 break = {mB['CM_min(c/4)']}, onset at α = {f(mB['stall_onset_alpha_deg']):.2f}°",
     "metrics_B CM_min / stall onset"),
    ("README.md", f"{mA['Cp_closure_error_pct']} % (Case A)) and {mB['Cp_closure_error_pct']} % (Case B)".replace("))", ")"),
     "metrics_A/B Cp_closure_error_pct"),
    # the cycle-wide closure is an ABSOLUTE residual: the worst instantaneous
    # percentage is not quotable, because the loop passes through C_L = 0.09
    ("README.md", f"{mA['Cp_closure_worst_dCL_cycle']} (Case A) and "
                  f"{mB['Cp_closure_worst_dCL_cycle']} (Case B)",
     "metrics_A/B Cp_closure_worst_dCL_cycle"),
    ("README.md", f"{mA['Cp_closure_worst_dCL_pct_of_CLmax']} % and "
                  f"{mB['Cp_closure_worst_dCL_pct_of_CLmax']} % of",
     "metrics_A/B Cp_closure_worst_dCL_pct_of_CLmax"),
    ("README.md", f"{mA['Cp_DSV_core_min']} (Case A) and {mB['Cp_DSV_core_min']} (Case B)",
     "metrics_A/B Cp_DSV_core_min"),
    # the vortex's derived properties, quoted in the same paragraph
    ("README.md", f"({mA['DSV_circulation_over_Uc']} U c for Case A)",
     "metrics_A DSV_circulation_over_Uc"),
    ("README.md", f"({mA['DSV_core_radius_chords']} c, swirling at "
                  f"{mA['DSV_peak_swirl_over_U']} U)",
     "metrics_A DSV_core_radius_chords / DSV_peak_swirl_over_U"),
    ("README.md", f"{mA['DSV_core_radius_cells']} of a cell",
     "metrics_A DSV_core_radius_cells"),
    ("README.md", f"({mA['DSV_induced_lift_dCL']} for Case A)",
     "metrics_A DSV_induced_lift_dCL"),
    # ---- damping ----------------------------------------------------------
    ("README.md", f"{mA['aero_damping_Xi_normalised']} (Case A) and {mB['aero_damping_Xi_normalised']} (Case B)",
     "metrics_A/B aero_damping_Xi_normalised"),
    ("README.md", f"±{vs['solver damping neutral band (us.DAMPING_TOL)']} band",
     "validation_realdata_summary damping band"),
    ("README.md", f"of {vs['mean |Xi_hat| model-exp discrepancy']} (max {vs['max |Xi_hat| model-exp discrepancy']})",
     "validation_realdata_summary Xi_hat spread"),
    # ---- validation -------------------------------------------------------
    ("README.md", f"mean peak-lift error of {vs['mean |CLmax| error [%]']} %",
     "validation_realdata_summary mean |CLmax| error"),
    ("README.md", f"mean moment-break error of\n{vs['mean |CMmin| error [abs]']}",
     "validation_realdata_summary mean |CMmin| error"),
    ("README.md", f"C_L loop is {vs['mean RMS_CL']}", "validation_realdata_summary mean RMS_CL"),
    # two of the four held-out frames peak past the static polar's calibration
    # limit, and are about twice as inaccurate; the split is quoted, so guard it
    ("README.md", f"{vs['mean RMS_CL, inside the calibration range']} over the two held-out frames",
     "validation_realdata_summary mean RMS_CL inside calibration"),
    ("README.md", f"and {vs['mean RMS_CL, peaking past it']} over the two that do not",
     "validation_realdata_summary mean RMS_CL past calibration"),
    ("README.md", f"past the {f(vs['static polar calibrated to [deg]']):.0f}° the static polar",
     "validation_realdata_summary static polar calibration limit"),
    ("README.md", f"slope {sl.pct_error} %, C_L,max {cm.pct_error} %",
     "validation_static pct_error"),
    # the AMES-01 cross-check is quoted loosely ("about 10 %"); assert the
    # rounded value the prose implies still holds
    ("README.md", f"agreement to about {xcheck_pct:.0f} %",
     "validation_nasa_real frame_25104 CLmax model vs exp"),
    # ---- overview ---------------------------------------------------------
    ("00_overview/case_definition.md", f"R = {f(B['rotor_radius_R']):.2f} m",
     "flow_conditions rotor_radius_R"),
    ("00_overview/case_definition.md", f"ΩR = {f(B['rotor_tip_speed_OmegaR']):.1f} m/s",
     "flow_conditions OmegaR"),
    ("00_overview/case_definition.md", f"U = M·a = {f(A['freestream_velocity_U']):.2f} m/s",
     "flow_conditions U (A)"),
    ("00_overview/case_definition.md",
     f"= {f(B['freestream_velocity_U']):.2f} m/s so M = {f(B['freestream_mach_M']):.4f}",
     "flow_conditions U and M (B)"),
]

_cache = {}
def _text(rel):
    """Document text, or None if the document is not part of this checkout.
    Both documents ship today (00_overview/case_definition.md is tracked; it was
    excluded from the repository in an earlier revision and this note still said
    so). The skip is kept anyway, so that running the checker against a partial
    export reports what it could not check instead of crashing on it -- and
    main() names anything skipped rather than counting it as passed."""
    if rel not in _cache:
        p = ROOT / rel
        _cache[rel] = p.read_text(encoding="utf-8") if p.exists() else None
    return _cache[rel]


def main():
    bad, skipped = [], set()
    for doc, needle, source in CLAIMS:
        t = _text(doc)
        if t is None:
            skipped.add(doc); continue
        if needle not in t:
            bad.append((doc, needle, source))
    if bad:
        print("[claims] %d claim(s) no longer match the artifacts:\n" % len(bad))
        for doc, needle, source in bad:
            print(f"  {doc}: expected to find {needle!r}")
            print(f"      (derived from {source})")
        print("\nUpdate the prose to the values above, or explain the difference.")
        return 1
    n = len([c for c in CLAIMS if c[0] not in skipped])
    docs = sorted({c[0] for c in CLAIMS} - skipped)
    print(f"[claims] {n} quoted numbers in {', '.join(docs)} all trace to "
          "generated artifacts"
          + (f" (skipped, not in this checkout: {', '.join(sorted(skipped))})"
             if skipped else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
