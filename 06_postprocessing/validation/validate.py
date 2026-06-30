# -*- coding: utf-8 -*-
"""
06_postprocessing / validation / validate.py
--------------------------------------------
Layer 1 of validation: STATIC polar (calibration check) + the calibration
record. The dynamic validation is performed against REAL digitised experimental
loops in validate_nasa_real.py (NASA TM-84245), and the certification harness is
validate_digitized.py.

STATIC SOURCES (recorded):
  [S1] Sheldahl & Klimas (1981) SAND80-2114, Sandia National Laboratories.
  [S2] Abbott & von Doenhoff (1959) "Theory of Wing Sections", Dover.
  [S3] McCroskey (1987) NASA TM-100019.
DYNAMIC DATA SOURCE:
  [S4] McAlister, Pucci, McCroskey & Carr (1982) NASA TM-84245 (digitised via
       the Pancini BL-DSM-JFS-2021 repository); airfoil identity confirmed from
       that repository's load_frame.m (frames 7019-14220 = NACA 0012).
  [S6] Leishman (2006) "Principles of Helicopter Aerodynamics" (model basis).
"""
import sys
import numpy as np
import pandas as pd
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
from aero_style import apply_style, PALETTE, INK, INK_SOFT
import matplotlib.pyplot as plt
apply_style(); plt.rcParams["figure.constrained_layout.use"] = True
SOL = ROOT/"05_solution"; SETUP = ROOT/"03_model_setup"

# ---------------- STATIC VALIDATION (calibration check) ----------------
mp = pd.read_csv(SOL/"model_static_polar.csv")
ref = pd.read_csv(SETUP/"static_polar_reference.csv")
Cl_model = np.interp(ref["alpha_deg"], mp["alpha_deg"], mp["Cl_model"])
def slope(a, cl):
    m = a <= 8; return np.polyfit(a[m], cl[m], 1)[0]
sm = slope(ref["alpha_deg"].values, Cl_model); sr = slope(ref["alpha_deg"].values, ref["Cl"].values)
clmax_m, clmax_r = Cl_model.max(), ref["Cl"].max()
ast_m = ref["alpha_deg"].values[np.argmax(Cl_model)]; ast_r = ref["alpha_deg"].values[np.argmax(ref["Cl"].values)]
lin = ref["alpha_deg"] <= 12; rmse = float(np.sqrt(np.mean((Cl_model[lin]-ref["Cl"][lin])**2)))
pd.DataFrame([
 ["lift-curve slope a0 [1/deg]", round(sm,4), round(sr,4), round(abs(sm-sr),4),
  round(100*abs(sm-sr)/sr,2), "S1,S2"],
 ["CL_max (static)", round(clmax_m,3), round(clmax_r,3), round(abs(clmax_m-clmax_r),3),
  round(100*abs(clmax_m-clmax_r)/clmax_r,2), "S1,S3"],
 ["stall angle [deg]", round(ast_m,2), round(ast_r,2), round(abs(ast_m-ast_r),2),
  round(100*abs(ast_m-ast_r)/ast_r,2), "S1,S3"],
 ["RMSE C_L (α≤12°)", round(rmse,4), 0.0, round(rmse,4), "-", "S1,S2"],
], columns=["metric","model","reference","abs_error","pct_error","source"]
).to_csv(HERE/"validation_static.csv", index=False)

fig, ax = plt.subplots(figsize=(7.6, 5))
ax.plot(mp["alpha_deg"], mp["Cl_model"], color=PALETTE[0], lw=2.4, label="UNISTALL model")
ax.plot(ref["alpha_deg"], ref["Cl"], "s", color=PALETTE[1], ms=6, label="published static [S1,S2,S3]")
ax.set_xlabel("angle of attack  α [deg]"); ax.set_ylabel("$C_L$")
ax.set_title("Static validation — NACA 0012 lift polar", pad=10)
ax.text(0.03, 0.97, f"a₀ err {100*abs(sm-sr)/sr:.1f}%\nCLmax err {100*abs(clmax_m-clmax_r)/clmax_r:.1f}%"
        f"\nRMSE {rmse:.3f}", transform=ax.transAxes, va="top", ha="left", fontsize=10,
        bbox=dict(boxstyle="round", fc="white", ec=INK_SOFT, alpha=0.9))
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2)
fig.savefig(HERE/"fig_validation_static_polar.png"); plt.close(fig)

# ---------------- CALIBRATION RECORD ----------------
import json
cfg = json.load(open(SETUP/"solver_config.json")); cc = cfg["calibrated_constants"]
pd.DataFrame([
 ["A1, A2 (indicial circulatory)", "0.30, 0.70", "literature [S6]", "fixed"],
 ["b1, b2 (indicial circulatory)", "0.14, 0.53", "literature [S6]", "fixed"],
 ["Tp (pressure lag)", "1.7", "literature [S6]", "fixed"],
 ["Tf (boundary-layer lag)", str(cc["Tf"]), "calibrated to real frame 9302", "tuned"],
 ["Tv, Tvl (vortex)", f"{cc['Tv']}, {cc['Tvl']}", "calibrated to real frame 9302", "tuned"],
 ["CN1 (DSV onset)", str(cc["CN1"]), "near static CLmax", "tuned"],
 ["k1 (CP aft-travel)", str(cc["k1"]), "calibrated to real frame 9302 C_M", "tuned"],
 ["cpv_amp (vortex CP travel)", str(cc["cpv_amp"]), "calibrated to real frame 9302 C_M", "tuned"],
 ["Bv (vortex-feed gain)", str(cc["Bv"]), "calibrated to real frame 9302 C_L overshoot", "tuned"],
 ["f(α) Kirchhoff fit", "from static polar", "calibrated to [S1,S2,S3]", "calibrated"],
], columns=["constant","value","basis","status"]).to_csv(HERE/"calibration_constants.csv", index=False)

print("[validate] static: slope err %.1f%%, CLmax err %.1f%%, stall err %.1f°, RMSE %.3f"
      % (100*abs(sm-sr)/sr, 100*abs(clmax_m-clmax_r)/clmax_r, abs(ast_m-ast_r), rmse))
