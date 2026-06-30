# -*- coding: utf-8 -*-
"""
06_postprocessing / validation / validate_nasa_real.py
------------------------------------------------------
PRIMARY dynamic validation against REAL digitised experimental loops.

DATA PROVENANCE
  Digitised C_L(α), C_M(α) oscillating-airfoil loops from McAlister, Pucci,
  McCroskey & Carr (1982), NASA TM-84245, via the open repository
  L. Pancini, "BL-DSM-JFS-2021" (NASA Data/frame_*.mat),
  https://github.com/luizpancini/BL-DSM-JFS-2021 .

AIRFOIL IDENTITY — CONFIRMED (no longer a guess).
  The repository's src/functions/load_frame.m maps frame number -> airfoil:
      frame in [ 7019 , 14220 ]  ->  NACA 0012
      frame in [24022 , 31310 ]  ->  AMES-01
      frame >= 67000             ->  NLR-7301
  All frames used below as NACA 0012 lie in [7019, 14220]; frame 25104 (AMES-01)
  is included only as a labelled cross-check.

PROTOCOL (calibrate-once / predict-the-rest):
  * the dynamic constants are calibrated ONLY on frame 9302 (= Case A,
    10deg +/- 10deg, M0.30, k0.10) and stored in solver_config.json;
  * they are then FROZEN and used to PREDICT the other (held-out) NACA 0012
    frames spanning light->deep stall and reduced frequency.
"""
import sys, json
import numpy as np
import pandas as pd
import scipy.io as sio
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FR = HERE/"experimental"/"nasa_frames"
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT/"04_solver"))
from aero_style import apply_style, PALETTE, INK, INK_SOFT
import matplotlib.pyplot as plt
import unistall_solver as us
apply_style(); plt.rcParams["figure.constrained_layout.use"] = True

cfg = json.load(open(ROOT/"03_model_setup"/"solver_config.json"))
CNALPHA = cfg["lift_curve_slope_CNalpha_per_rad"]
FROZEN = dict(**cfg["indicial_circulatory"], **cfg["time_constants_semichords"])
FROZEN.update({k: v for k, v in cfg["calibrated_constants"].items() if k != "comment"})
stat = pd.read_csv(ROOT/"03_model_setup"/"static_polar_reference.csv")
f_static = us.calibrate_separation(stat["alpha_deg"], stat["Cl"], stat["Cd"], CNALPHA)

# role: 'calibration' (9302=Case A) or 'held-out prediction'; all NACA0012 except the cross-check
FRAMES = [
 ("frame_9302.mat",  "NACA 0012", "calibration (10°±10°, k0.10 = Case A)"),
 ("frame_9217.mat",  "NACA 0012", "held-out: deep stall 15°±10°, k0.10"),
 ("frame_9214.mat",  "NACA 0012", "held-out: deep stall 15°±10°, k0.05"),
 ("frame_7113.mat",  "NACA 0012", "held-out: light stall 10°±5°, k0.10"),
 ("frame_10118.mat", "NACA 0012", "held-out: 15°±5°, k0.10"),
 ("frame_25104.mat", "AMES-01",   "cross-check (NOT NACA0012): 10°±10°, k0.10"),
]
CHORD = 0.30

def loadframe(fn):
    d = sio.loadmat(str(FR/fn)); g = lambda k: float(d[k].ravel()[0])
    return dict(M=g("M"), k=g("k"), a0=g("alpha_0")*180/np.pi, da=g("delta_alpha")*180/np.pi,
                acl=d["alpha_exp_cl"].ravel(), cl=d["cl_exp"].ravel(),       # already in deg
                acm=d["alpha_exp_cm"].ravel(), cm=d["cm_exp"].ravel())

def stroke_split(a):
    ip = int(np.argmax(a)); s = np.full(len(a), "up", dtype="<U4"); s[ip+1:] = "down"; return s

def model_branches(o):
    a = o["alpha_deg"]; up = o["alpha_dot"] > 0; br = {}
    for nm, m in [("up", up), ("down", ~up)]:
        order = np.argsort(a[m]); br[nm] = (a[m][order], o["CL"][m][order], o["CM"][m][order])
    return br

rows = []
fig, axs = plt.subplots(len(FRAMES), 2, figsize=(11, 3.0*len(FRAMES)))
for i, (fn, airfoil, role) in enumerate(FRAMES):
    fr = loadframe(fn); U = fr["M"]*340.0
    o = us.solve_dynamic_stall(fr["a0"], fr["da"], fr["k"], fr["M"], CHORD, U, f_static,
                               CNalpha=CNALPHA, consts=FROZEN, n_per_cycle=720, n_cycles=6)
    br = model_branches(o)
    pd.DataFrame({"alpha_deg": np.round(fr["acl"],3), "CL_exp": np.round(fr["cl"],4)}
                 ).to_csv(HERE/f"exp_{fn.replace('.mat','')}_CL.csv", index=False)
    pd.DataFrame({"alpha_deg": np.round(fr["acm"],3), "CM_exp": np.round(fr["cm"],4)}
                 ).to_csv(HERE/f"exp_{fn.replace('.mat','')}_CM.csv", index=False)
    scl = stroke_split(fr["acl"]); scm = stroke_split(fr["acm"])
    mcl = np.array([np.interp(av, br[s][0], br[s][1]) for s, av in zip(scl, fr["acl"])])
    mcm = np.array([np.interp(av, br[s][0], br[s][2]) for s, av in zip(scm, fr["acm"])])
    rms_cl = float(np.sqrt(np.mean((mcl-fr["cl"])**2)))
    rms_cm = float(np.sqrt(np.mean((mcm-fr["cm"])**2)))
    rows.append([fn.replace(".mat",""), airfoil, role, round(fr["M"],3), round(fr["k"],3),
                 round(fr["a0"],1), round(fr["da"],1), round(rms_cl,4), round(rms_cm,4),
                 round(float(o["CL"].max()),3), round(float(fr["cl"].max()),3),
                 round(float(o["CM"].min()),3), round(float(fr["cm"].min()),3)])
    axs[i,0].plot(o["alpha_deg"], o["CL"], color=PALETTE[0], lw=2, label="UNISTALL")
    axs[i,0].plot(fr["acl"], fr["cl"], "o", color=PALETTE[1], ms=3.2, label="experiment")
    axs[i,0].set_ylabel("$C_L$"); axs[i,0].set_xlabel("α [deg]")
    axs[i,0].set_title(f"{fn.replace('.mat','')} [{airfoil}] — {role}   RMS$_{{CL}}$={rms_cl:.3f}",
                       pad=6, fontsize=8.5)
    axs[i,1].plot(o["alpha_deg"], o["CM"], color=PALETTE[0], lw=2, label="UNISTALL")
    axs[i,1].plot(fr["acm"], fr["cm"], "o", color=PALETTE[1], ms=3.2, label="experiment")
    axs[i,1].set_ylabel("$C_{M,c/4}$"); axs[i,1].set_xlabel("α [deg]")
    axs[i,1].set_title(f"moment   RMS$_{{CM}}$={rms_cm:.3f}", pad=6, fontsize=8.5)
    if i == 0:
        for ax in axs[i]:
            ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.75), ncol=2, fontsize=8)
fig.savefig(HERE/"fig_validation_nasa_real.png"); plt.close(fig)

res = pd.DataFrame(rows, columns=["frame","airfoil","role","M","k","alpha0_deg","amp_deg",
        "RMS_CL","RMS_CM","CLmax_model","CLmax_exp","CMmin_model","CMmin_exp"])
res["airfoil_source"] = "load_frame.m mapping (Pancini repo); data NASA TM-84245"
res.to_csv(HERE/"validation_nasa_real.csv", index=False)

# headline metrics over held-out NACA0012 frames only
ho = res[(res.airfoil == "NACA 0012") & (res.role.str.startswith("held-out"))]
clpe = (100*(ho.CLmax_model-ho.CLmax_exp).abs()/ho.CLmax_exp).mean()
cmpe = (ho.CMmin_model-ho.CMmin_exp).abs().mean()
pd.DataFrame({"metric": ["NACA0012 held-out frames", "mean RMS_CL", "mean RMS_CM",
                         "mean |CLmax| error [%]", "mean |CMmin| error [abs]",
                         "airfoil identity"],
             "value": [len(ho), round(ho.RMS_CL.mean(),3), round(ho.RMS_CM.mean(),3),
                       round(clpe,1), round(cmpe,3),
                       "CONFIRMED via load_frame.m (frames 7019-14220 = NACA0012)"]}
            ).to_csv(HERE/"validation_realdata_summary.csv", index=False)
print(f"[nasa-real] {len(ho)} held-out NACA0012 frames: meanRMS_CL={ho.RMS_CL.mean():.3f}, "
      f"mean|CLmax|err={clpe:.1f}%, mean|CMmin|err={cmpe:.3f}; airfoil CONFIRMED")
