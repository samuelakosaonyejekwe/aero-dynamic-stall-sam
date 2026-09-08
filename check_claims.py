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

f = float
sl = vst[vst.metric.str.startswith("lift-curve")].iloc[0]
cm = vst[vst.metric.str.startswith("CL_max")].iloc[0]

# (document, exact substring that must be present, what it came from)
CLAIMS = [
    # ---- case conditions --------------------------------------------------
    ("README.md", f"{f(B['chord_c']):.3f} m", "flow_conditions chord_c (B)"),
    ("README.md", f"{f(B['freestream_mach_M']):.3f}", "flow_conditions M (B)"),
    ("README.md", f"{f(B['reduced_frequency_k']):.3f}", "flow_conditions k (B)"),
    ("README.md", f"R = {f(B['rotor_radius_R']):.2f} m", "flow_conditions rotor_radius_R"),
    ("README.md", f"{f(B['rotor_speed_Omega']):.1f} rad/s", "flow_conditions rotor_speed_Omega"),
    ("README.md", f"ΩR = {f(B['rotor_tip_speed_OmegaR']):.1f} m/s", "flow_conditions OmegaR"),
    ("README.md", f"μ = {f(B['advance_ratio_mu']):.2f}", "flow_conditions advance_ratio_mu"),
    ("README.md", f"= {f(B['freestream_velocity_U']):.2f} m/s", "flow_conditions U (B)"),
    ("README.md", f"Mach is {f(B['freestream_mach_M']):.4f}", "flow_conditions M (B), 4 dp"),
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
    ("README.md", f"{mA['Cp_closure_error_pct']} % (Case A)", "metrics_A Cp_closure_error_pct"),
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
    00_overview/ is excluded from the public repository, so a clean checkout has
    the README and nothing else -- the check must skip what is not there rather
    than crash on it."""
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
