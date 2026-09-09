# -*- coding: utf-8 -*-
"""
06_postprocessing / validation / validate_nasa_real.py
------------------------------------------------------
PRIMARY dynamic validation against REAL digitised experimental loops.

Outputs (all to 06_postprocessing/validation/):
  validation_nasa_real.csv         one row per frame: conditions, role, errors,
                                   and whether the frame's peak incidence stays
                                   inside the static polar's calibration range
  validation_realdata_summary.csv  headline metrics over the held-out frames
  exp_<frame>_CL.csv               digitised lift points extracted from each .mat
  exp_<frame>_CM.csv               digitised moment points, likewise
  fig_validation_nasa_real.png     model-vs-experiment overlay, all frames

DATA PROVENANCE
  Digitised C_L(α), C_M(α) oscillating-airfoil loops from McCroskey, McAlister,
  Carr & Pucci (1982), NASA TM-84245 [S5], via the open repository
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
from aero_style import apply_style, PALETTE
import matplotlib.pyplot as plt
import unistall_solver as us
apply_style(); plt.rcParams["figure.constrained_layout.use"] = True

cfg = json.load(open(ROOT/"03_model_setup"/"solver_config.json"))
CNALPHA = cfg["lift_curve_slope_CNalpha_per_rad"]
FROZEN = dict(**cfg["indicial_circulatory"], **cfg["time_constants_semichords"])
FROZEN.update({k: v for k, v in cfg["calibrated_constants"].items() if k != "comment"})
stat = pd.read_csv(ROOT/"03_model_setup"/"static_polar_reference.csv")
f_static = us.calibrate_separation(stat["alpha_deg"], stat["Cl"], stat["Cd"], CNALPHA)
# The separation law is fitted by inverse Kirchhoff to a polar tabulated only to
# ALPHA_CAL_MAX; past that, calibrate_separation() decays f towards full
# separation on an assumption. Everything else in this study flags that boundary
# -- the response surface is bounded by it, the model polar carries a
# within_calibration column, the sampling incidences are capped just inside it,
# and both calibration figures draw the tail dashed. The dynamic validation did
# not, and two of its four held-out frames peak at 25 deg, five degrees past it.
# It is published per frame below rather than left to a reader to notice.
ALPHA_CAL_MAX = float(stat["alpha_deg"].max())

# role: 'calibration' (9302=Case A) or 'held-out prediction'; all NACA0012 except the cross-check
FRAMES = [
 ("frame_9302.mat",  "NACA 0012", "calibration (10°±10°, k0.10 = Case A)"),
 ("frame_9217.mat",  "NACA 0012", "held-out: deep stall 15°±10°, k0.10"),
 ("frame_9214.mat",  "NACA 0012", "held-out: deep stall 15°±10°, k0.05"),
 ("frame_7113.mat",  "NACA 0012", "held-out: light stall 10°±5°, k0.10"),
 ("frame_10118.mat", "NACA 0012", "held-out: 15°±5°, k0.10"),
 ("frame_25104.mat", "AMES-01",   "cross-check (NOT NACA0012): 10°±10°, k0.10"),
]
# The UIBS march is chord- AND speed-independent: omega = 2kU/c, ds = 2U dt/c and
# dt/T_I = pi/(k n_per_cycle K_alpha) all scale so that every recurrence depends
# on (k, M) alone. CHORD and the speed of sound below are therefore nominal
# values that do not affect any load reported here; they are read from the case
# setup rather than restated so that nothing in this file can silently disagree
# with it. (test_frame_independence() at the bottom asserts the independence.)
_fl = pd.read_csv(ROOT/"03_model_setup"/"flow_conditions.csv"
                  ).set_index("parameter")["case_A_validation"]
CHORD = float(_fl["chord_c"])
A_SND = float(_fl["speed_of_sound_a"])

def loadframe(fn):
    d = sio.loadmat(str(FR/fn)); g = lambda k: float(d[k].ravel()[0])
    return dict(M=g("M"), k=g("k"), a0=g("alpha_0")*180/np.pi, da=g("delta_alpha")*180/np.pi,
                acl=d["alpha_exp_cl"].ravel(), cl=d["cl_exp"].ravel(),       # already in deg
                acm=d["alpha_exp_cm"].ravel(), cm=d["cm_exp"].ravel())

def stroke_split(a):
    """Label each experimental point up/down. The record is a closed loop that may
    start anywhere, so split on BOTH turning points rather than assuming it opens
    on the upstroke."""
    a = np.asarray(a, float)
    imax, imin = int(np.argmax(a)), int(np.argmin(a))
    s = np.empty(len(a), dtype="<U4")
    if imin <= imax:                       # ... min ... max ...  -> rising between them
        s[:] = "down"; s[imin:imax+1] = "up"
    else:                                  # ... max ... min ...  -> falling between them
        s[:] = "up";   s[imax+1:imin+1] = "down"
    return s

def model_branches(o):
    a = o["alpha_deg"]; up = o["alpha_dot"] > 0; br = {}
    for nm, m in [("up", up), ("down", ~up)]:
        order = np.argsort(a[m]); br[nm] = (a[m][order], o["CL"][m][order], o["CM"][m][order])
    return br

rows = []
# The chord-independence invariant is checked HERE, before any result is
# written. It used to run at the END of the script, after every output file
# had already been published -- so a march that violated the assumption the
# whole validation rests on would have left those results standing as if
# they were sound.
# ---- invariant: the march really is chord/speed independent, as claimed above
def _frame_independence_check():
    o1 = us.solve_dynamic_stall(10.0, 10.0, 0.10, 0.30, 0.30, 0.30*A_SND, f_static,
                                CNalpha=CNALPHA, consts=FROZEN, n_per_cycle=360, n_cycles=3)
    o2 = us.solve_dynamic_stall(10.0, 10.0, 0.10, 0.30, 1.70, 0.30*A_SND, f_static,
                                CNalpha=CNALPHA, consts=FROZEN, n_per_cycle=360, n_cycles=3)
    d = float(np.max(np.abs(o1["CL"] - o2["CL"])))
    assert d < 1e-10, f"march is NOT chord-independent: max |dCL| = {d:.3e}"
_frame_independence_check()

fig, axs = plt.subplots(len(FRAMES), 2, figsize=(11, 3.0*len(FRAMES)))
for i, (fn, airfoil, role) in enumerate(FRAMES):
    fr = loadframe(fn); U = fr["M"]*A_SND
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
    # ---- normalised aerodynamic damping, model vs EXPERIMENT, at the same
    #      condition. This is what sets the neutral band in the solver
    #      (us.DAMPING_TOL): a damping residual smaller than the discrepancy
    #      below cannot be claimed as a finding. The band used to be justified
    #      as "no larger than the time-step discretisation error", which is
    #      false by a factor of ~45 -- refining n_per_cycle 720 -> 5760 moves
    #      Xi_hat by 0.0004. ----
    _, xh_mod = us.aerodynamic_damping(o["alpha_deg"], o["CM"], normalise=True)
    _ae = np.append(fr["acm"], fr["acm"][0])          # close the measured loop
    _ce = np.append(fr["cm"], fr["cm"][0])
    _ar = np.radians(_ae)
    _box = (_ce.max()-_ce.min())*(_ar.max()-_ar.min())
    xh_exp = float(-us._trapz(_ce, _ar)/_box) if _box > 0 else np.nan
    _peak = fr["a0"] + fr["da"]
    rows.append([fn.replace(".mat",""), airfoil, role, round(fr["M"],3), round(fr["k"],3),
                 round(fr["a0"],1), round(fr["da"],1),
                 round(_peak,1), bool(_peak <= ALPHA_CAL_MAX + 1e-9),
                 round(rms_cl,4), round(rms_cm,4),
                 round(float(o["CL"].max()),3), round(float(fr["cl"].max()),3),
                 round(float(o["CM"].min()),3), round(float(fr["cm"].min()),3),
                 round(xh_mod,4), round(xh_exp,4)])
    axs[i,0].plot(o["alpha_deg"], o["CL"], color=PALETTE[0], lw=2, label="UNISTALL")
    axs[i,0].plot(fr["acl"], fr["cl"], "o", color=PALETTE[1], ms=3.2, label="experiment")
    axs[i,0].set_ylabel("$C_L$"); axs[i,0].set_xlabel("α [deg]")
    # on its own line: appended inline it pushed the frame_9217 title off the
    # left edge of the figure, since the title is centred on the left axes
    _ext = "" if _peak <= ALPHA_CAL_MAX + 1e-9 else \
           f"\npeaks at {_peak:.0f}° — past the {ALPHA_CAL_MAX:.0f}° f(α) is fitted to"
    axs[i,0].set_title(f"{fn.replace('.mat','')} [{airfoil}] — {role}   "
                       f"RMS$_{{CL}}$={rms_cl:.3f}{_ext}", pad=6, fontsize=8.5)
    if _ext:
        axs[i,0].axvspan(ALPHA_CAL_MAX, max(fr["acl"].max(), _peak), color=PALETTE[3],
                         alpha=0.10, zorder=0)
        axs[i,1].axvspan(ALPHA_CAL_MAX, max(fr["acm"].max(), _peak), color=PALETTE[3],
                         alpha=0.10, zorder=0)
    axs[i,1].plot(o["alpha_deg"], o["CM"], color=PALETTE[0], lw=2, label="UNISTALL")
    axs[i,1].plot(fr["acm"], fr["cm"], "o", color=PALETTE[1], ms=3.2, label="experiment")
    axs[i,1].set_ylabel("$C_{M,c/4}$"); axs[i,1].set_xlabel("α [deg]")
    axs[i,1].set_title(f"moment   RMS$_{{CM}}$={rms_cm:.3f}", pad=6, fontsize=8.5)
# one figure-level legend (the two per-axes legends were identical and orphaned).
# constrained_layout does not reserve space for a figure legend, so shrink the
# layout rect first, otherwise the legend lands on top of the first row's title.
fig.get_layout_engine().set(rect=(0, 0, 1, 0.972))
h, l = axs[0, 0].get_legend_handles_labels()
fig.legend(h, l, loc="upper center", bbox_to_anchor=(0.5, 1.0), ncol=2,
           fontsize=9, frameon=True)
fig.savefig(HERE/"fig_validation_nasa_real.png"); plt.close(fig)

res = pd.DataFrame(rows, columns=["frame","airfoil","role","M","k","alpha0_deg","amp_deg",
        "peak_alpha_deg","within_static_calibration",
        "RMS_CL","RMS_CM","CLmax_model","CLmax_exp","CMmin_model","CMmin_exp",
        "Xihat_model","Xihat_exp"])
res["airfoil_source"] = "load_frame.m mapping (Pancini repo); data NASA TM-84245"
res.to_csv(HERE/"validation_nasa_real.csv", index=False)

# headline metrics over held-out NACA0012 frames only
ho = res[(res.airfoil == "NACA 0012") & (res.role.str.startswith("held-out"))]
clpe = (100*(ho.CLmax_model-ho.CLmax_exp).abs()/ho.CLmax_exp).mean()
cmpe = (ho.CMmin_model-ho.CMmin_exp).abs().mean()
# resolution of the damping metric: the model-vs-experiment spread in Xi_hat over
# every real NACA 0012 frame. This is the evidence behind us.DAMPING_TOL.
n0012 = res[res.airfoil == "NACA 0012"]
dxi = (n0012.Xihat_model - n0012.Xihat_exp).abs()
# Two of the four held-out frames peak at 25 deg, i.e. 5 deg past the incidence
# the separation law is fitted to, and they are twice as inaccurate as the two
# that stay inside it. Reporting only the pooled mean hid that, in a study that
# flags the same boundary in three other places.
_in  = ho[ho.within_static_calibration]
_out = ho[~ho.within_static_calibration]
pd.DataFrame({"metric": ["NACA0012 held-out frames", "mean RMS_CL", "mean RMS_CM",
                         "mean |CLmax| error [%]", "mean |CMmin| error [abs]",
                         "held-out frames inside the static calibration range",
                         "held-out frames peaking past it",
                         "static polar calibrated to [deg]",
                         "mean RMS_CL, inside the calibration range",
                         "mean RMS_CL, peaking past it",
                         "mean |Xi_hat| model-exp discrepancy",
                         "max |Xi_hat| model-exp discrepancy",
                         "solver damping neutral band (us.DAMPING_TOL)",
                         "airfoil identity"],
             "value": [len(ho), round(ho.RMS_CL.mean(),3), round(ho.RMS_CM.mean(),3),
                       round(clpe,1), round(cmpe,3),
                       len(_in), len(_out), ALPHA_CAL_MAX,
                       round(float(_in.RMS_CL.mean()),3) if len(_in) else float("nan"),
                       round(float(_out.RMS_CL.mean()),3) if len(_out) else float("nan"),
                       round(float(dxi.mean()),3), round(float(dxi.max()),3),
                       us.DAMPING_TOL,
                       "CONFIRMED via load_frame.m (frames 7019-14220 = NACA0012)"]}
            ).to_csv(HERE/"validation_realdata_summary.csv", index=False)
if us.DAMPING_TOL < float(dxi.mean()) - 1e-9:
    print(f"[nasa-real] WARNING: solver DAMPING_TOL={us.DAMPING_TOL} is tighter than the "
          f"measured model-vs-experiment spread {dxi.mean():.3f}; damping verdicts "
          "outside the band are not supported by the validation data.")
print(f"[nasa-real] {len(_out)} of {len(ho)} held-out frames peak past the "
      f"{ALPHA_CAL_MAX:.0f}deg calibration limit: meanRMS_CL {_in.RMS_CL.mean():.3f} "
      f"inside vs {_out.RMS_CL.mean():.3f} outside")
print(f"[nasa-real] {len(ho)} held-out NACA0012 frames: meanRMS_CL={ho.RMS_CL.mean():.3f}, "
      f"mean|CLmax|err={clpe:.1f}%, mean|CMmin|err={cmpe:.3f}; airfoil CONFIRMED; "
      f"Xi_hat model-exp spread mean={dxi.mean():.3f} max={dxi.max():.3f} "
      f"(solver band {us.DAMPING_TOL})")

