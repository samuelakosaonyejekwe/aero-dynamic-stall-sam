# -*- coding: utf-8 -*-
"""
06_postprocessing / validation / validate_digitized.py   (OPTION 2)
------------------------------------------------------------------
Certification-grade harness. Drop a CSV of EXPERIMENTAL points digitised from a
specific figure of a cited report (e.g. McAlister TP-1100 / McCroskey TM-84245)
into  06_postprocessing/validation/experimental/  and this script will:

  * run the UNISTALL solver at the matching condition (constants frozen),
  * match each experimental point to the model loop on the correct stroke,
  * compute point-by-point error metrics: RMS C_L / C_M, max |ΔC_L|,
    moment-break error, and lift-loop-area error, and
  * produce an experiment-vs-model overlay figure.

If no experimental data are present, it writes a TEMPLATE + instructions and
exits WITHOUT fabricating any data (honesty by construction).

Experimental CSV schema (one file per test condition):
    alpha_deg, CL [, CM] [, CD] [, stroke]     stroke in {up,down} (optional)
Conditions are read from experimental/conditions.csv:
    file, mean, amp, k, M, c, U, source
"""
import sys, json
import numpy as np
import pandas as pd
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
EXP = HERE/"experimental"; EXP.mkdir(exist_ok=True)
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT/"04_solver"))
from aero_style import apply_style, PALETTE, INK, INK_SOFT
import matplotlib.pyplot as plt
import unistall_solver as us
apply_style(); plt.rcParams["figure.constrained_layout.use"] = True

SETUP = ROOT/"03_model_setup"; SOL = ROOT/"05_solution"
cfg = json.load(open(SETUP/"solver_config.json"))
CNALPHA = cfg["lift_curve_slope_CNalpha_per_rad"]
FROZEN = dict(**cfg["indicial_circulatory"], **cfg["time_constants_semichords"])
FROZEN.update({k: v for k, v in cfg["calibrated_constants"].items() if k != "comment"})
stat = pd.read_csv(SETUP/"static_polar_reference.csv")
f_static = us.calibrate_separation(stat["alpha_deg"], stat["Cl"], stat["Cd"], CNALPHA)

def write_template():
    tmpl = pd.DataFrame({
        "alpha_deg": [0, 5, 10, 14, 17, 19, 20, 18, 14, 10, 5, 0],
        "CL": ["" for _ in range(12)], "CM": ["" for _ in range(12)],
        "CD": ["" for _ in range(12)],
        "stroke": ["up"]*6 + ["down"]*6})
    tmpl.to_csv(EXP/"TEMPLATE_experiment.csv", index=False)
    pd.DataFrame({"file": ["TEMPLATE_experiment.csv"], "mean": [10], "amp": [10],
                  "k": [0.10], "M": [0.30], "c": [0.30], "U": [102.0],
                  "source": ["McAlister TP-1100 Fig. X / McCroskey TM-84245 frame Y"]}
                 ).to_csv(EXP/"conditions.csv", index=False)
    (EXP/"README.txt").write_text(
        "Digitised-experiment validation harness\n"
        "=======================================\n"
        "1. Open the target figure in McAlister NASA TP-1100 or McCroskey NASA\n"
        "   TM-84245 (e.g. a C_L vs alpha or C_M vs alpha dynamic-stall loop).\n"
        "2. Digitise it (e.g. WebPlotDigitizer https://automeris.io) into a CSV\n"
        "   with columns: alpha_deg, CL [, CM] [, CD] [, stroke=up/down].\n"
        "3. Add a matching row to conditions.csv (file, mean, amp, k, M, c, U,\n"
        "   source) describing the exact test point.\n"
        "4. Re-run:  python3 validate_digitized.py\n"
        "   -> writes validation_digitized_<file>.csv + overlay figure with\n"
        "      true point-by-point RMS / peak / loop-area errors.\n"
        "No data are fabricated; results appear only for files you provide.\n")

def branch(a, y, dadt):
    up = dadt > 0
    out = {}
    for name, m in [("up", up), ("down", ~up)]:
        order = np.argsort(a[m])
        out[name] = (a[m][order], y[m][order])
    return out

def interp_branch(branches, stroke, aq):
    xa, ya = branches[stroke if stroke in branches else "up"]
    return np.interp(aq, xa, ya)

cond_path = EXP/"conditions.csv"
exp_files = [f for f in EXP.glob("*.csv") if f.name not in ("conditions.csv",)
             and not f.name.startswith("validation_")]
real = [f for f in exp_files if not f.name.startswith("TEMPLATE")]

if not cond_path.exists() or not real:
    write_template()
    print("[digitized] no experimental data found — wrote TEMPLATE + README to "
          f"{EXP.relative_to(ROOT)} . Add digitised CSV(s) + conditions.csv, then re-run.")
    sys.exit(0)

cond = pd.read_csv(cond_path).set_index("file")
summary = []
for f in real:
    if f.name not in cond.index:
        print(f"[digitized] skip {f.name}: no row in conditions.csv"); continue
    c = cond.loc[f.name]
    exp = pd.read_csv(f)
    o = us.solve_dynamic_stall(c["mean"], c["amp"], c["k"], c["M"], c["c"], c["U"],
                               f_static, CNalpha=CNALPHA, consts=FROZEN,
                               n_per_cycle=720, n_cycles=6)
    a = o["alpha_deg"]; dadt = o["alpha_dot"]
    brCL = branch(a, o["CL"], dadt); brCM = branch(a, o["CM"], dadt)
    strokes = exp["stroke"] if "stroke" in exp else pd.Series(["up"]*len(exp))
    mCL = np.array([interp_branch(brCL, s, av) for s, av in zip(strokes, exp["alpha_deg"])])
    rms_cl = float(np.sqrt(np.mean((mCL-exp["CL"])**2)))
    maxe_cl = float(np.max(np.abs(mCL-exp["CL"])))
    rec = {"file": f.name, "n_points": len(exp), "RMS_CL": round(rms_cl, 4),
           "maxAbs_CL": round(maxe_cl, 4)}
    if "CM" in exp and exp["CM"].notna().any():
        mCM = np.array([interp_branch(brCM, s, av) for s, av in zip(strokes, exp["alpha_deg"])])
        rec["RMS_CM"] = round(float(np.sqrt(np.mean((mCM-exp["CM"])**2))), 4)
        rec["CMbreak_err"] = round(float(abs(o["CM"].min()-exp["CM"].min())), 4)
    rec["CLmax_err"] = round(float(abs(o["CL"].max()-exp["CL"].max())), 4)
    rec["source"] = str(c["source"])
    summary.append(rec)
    pd.DataFrame([rec]).to_csv(HERE/f"validation_digitized_{f.stem}.csv", index=False)

    # overlay figure
    fig, axs = plt.subplots(1, 2 if "CM" in exp else 1,
                            figsize=(11 if "CM" in exp else 6.5, 4.8), squeeze=False)
    axs[0][0].plot(a, o["CL"], color=PALETTE[0], lw=2, label="UNISTALL (frozen)")
    axs[0][0].plot(exp["alpha_deg"], exp["CL"], "o", color=PALETTE[1], ms=5,
                   label="experiment (digitised)")
    axs[0][0].set_xlabel("α [deg]"); axs[0][0].set_ylabel("$C_L$")
    axs[0][0].set_title(f"Lift loop — {f.stem}  (RMS={rms_cl:.3f})", pad=10)
    axs[0][0].legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2)
    if "CM" in exp:
        axs[0][1].plot(a, o["CM"], color=PALETTE[0], lw=2, label="UNISTALL (frozen)")
        axs[0][1].plot(exp["alpha_deg"], exp["CM"], "o", color=PALETTE[1], ms=5,
                       label="experiment (digitised)")
        axs[0][1].set_xlabel("α [deg]"); axs[0][1].set_ylabel("$C_{M,c/4}$")
        axs[0][1].set_title("Moment loop", pad=10)
        axs[0][1].legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2)
    fig.savefig(HERE/f"fig_validation_digitized_{f.stem}.png"); plt.close(fig)

if summary:
    pd.DataFrame(summary).to_csv(HERE/"validation_digitized_summary.csv", index=False)
    print(f"[digitized] validated {len(summary)} experimental file(s); "
          "see validation_digitized_summary.csv")
