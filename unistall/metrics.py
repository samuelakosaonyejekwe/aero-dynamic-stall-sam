# -*- coding: utf-8 -*-
"""
unistall / metrics.py
-----------------
ONE definition of every error measure by which a dynamic-stall model is judged
against a measured loop. Both models (tabulated and fitted separation law)
are scored by this file and by nothing else, on every frame, so that two numbers in the manuscript are always the same
kind of number.

Author: Akosa Samuel Onyejekwe (independent)

Per frame, against the measured loop:
  nRMS_CL, nRMS_CM, nRMS_CD   RMS error divided by the MEASURED range of that
                              quantity over the loop
  CLmax_err_pct               100*(model - measured)/measured, signed
  dalpha_CLmax_deg            incidence of C_L,max, model - measured
  dalpha_Mstall_deg           incidence of moment stall, model - measured
  Xi_model, Xi_exp            standard cycle damping coefficient
  Xi_sign_agree               whether the two have the same sign
  CD_neg_frac                 fraction of the cycle with model C_D < 0
The raw RMS values and the measured ranges are written beside them so the
normalisation can be undone.

Aggregates are means over the calibration frames and over the held-out frames
SEPARATELY, each with a bootstrap 95 % interval over frames.

Definitions that need stating because they are choices:

  MOMENT STALL. The incidence, on the upstroke, at which C_M first falls
  MSTALL_DROP below its attached-flow level, the attached-flow level being the
  median C_M over the lowest third of the upstroke's incidence range. The same
  rule is applied to the model and to the measurement. A loop that never drops
  that far has no moment stall and reports NaN.

  CYCLE DAMPING. Xi = -(1/(pi*alpha_1^2)) * closed integral of C_M d(alpha),
  alpha_1 the pitch amplitude in radians. Positive is positive damping. This is
  the customary normalisation.

"""
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.io as sio

HERE = Path(__file__).resolve().parent
from unistall.strokes import stroke_split
from unistall.paths import CACHE, DATA

_trapz = getattr(np, "trapezoid", None) or np.trapz

MSTALL_DROP = 0.05          # C_M drop below the attached level that marks moment stall
CLOSED_FRACTION = 0.8         # each stroke of a measured loop must span this much of the incidence range
BOOT_N, BOOT_SEED = 10000, 20261006

FRAME_COLUMNS = ["frame", "airfoil", "set", "M", "k", "alpha0_deg", "amp_deg",
                 "RMS_CL", "RMS_CM", "RMS_CD", "range_CL", "range_CM", "range_CD",
                 "nRMS_CL", "nRMS_CM", "nRMS_CD",
                 "CLmax_model", "CLmax_exp", "CLmax_err_pct",
                 "alpha_CLmax_model", "alpha_CLmax_exp", "dalpha_CLmax_deg",
                 "alpha_Mstall_model", "alpha_Mstall_exp", "dalpha_Mstall_deg",
                 "Mstall_model_found", "Mstall_exp_found",
                 "CMmin_model", "CMmin_exp", "CDmax_model", "CDmax_exp",
                 "Xi_model", "Xi_exp", "moment_loop_closed", "Xi_sign_agree", "CD_neg_frac", "CD_min_model"]


# --------------------------------------------------------------------------- #
#  measured frames
# --------------------------------------------------------------------------- #
def load_frame(path: object) -> dict:
    """Read one frame_*.mat. Incidences are stored in degrees, the mean angle
    and amplitude in radians."""
    d = sio.loadmat(str(path))
    def g(k: str) -> float:
        """Scalar field k of the frame file."""
        return float(d[k].ravel()[0])

    def r(k: str) -> np.ndarray:
        """Array field k of the frame file."""
        return d[k].ravel().astype(float)
    return dict(frame=Path(path).stem, M=g("M"), k=g("k"),
                a0=np.degrees(g("alpha_0")), da=np.degrees(g("delta_alpha")),
                acl=r("alpha_exp_cl"), cl=r("cl_exp"),
                acm=r("alpha_exp_cm"), cm=r("cm_exp"),
                acd=r("alpha_exp_cd"), cd=r("cd_exp"))


# --------------------------------------------------------------------------- #
#  pieces
# --------------------------------------------------------------------------- #
def _branches(alpha, alpha_dot, y):
    """Model curve split into up- and down-stroke, each sorted on incidence."""
    up = np.asarray(alpha_dot) > 0
    out = {}
    for name, m in (("up", up), ("down", ~up)):
        o = np.argsort(alpha[m])
        out[name] = (alpha[m][o], y[m][o])
    return out


def loop_rms(a_exp: np.ndarray, y_exp: np.ndarray, alpha: np.ndarray, alpha_dot: np.ndarray,
             y_model: np.ndarray) -> float:
    """RMS of (model - measured) at the measured points, each compared with the
    model on its own stroke."""
    br = _branches(alpha, alpha_dot, y_model)
    s = stroke_split(a_exp)
    ym = np.array([np.interp(av, *br[st]) for st, av in zip(s, a_exp, strict=False)])
    return float(np.sqrt(np.mean((ym - y_exp)**2)))


def moment_stall_alpha(alpha_up: np.ndarray, cm_up: np.ndarray) -> float:
    """Upstroke incidence at which C_M first falls MSTALL_DROP below its
    attached level. `alpha_up` must be increasing. NaN if it never does."""
    a = np.asarray(alpha_up, float)
    c = np.asarray(cm_up, float)
    if len(a) < 4:
        return float("nan")
    lo = a <= a[0] + (a[-1] - a[0])/3.0
    thr = float(np.median(c[lo])) - MSTALL_DROP
    below = np.where(c <= thr)[0]
    below = below[below > 0]
    if len(below) == 0:
        return float("nan")
    i = int(below[0])
    if c[i-1] == c[i]:
        return float(a[i])
    return float(a[i-1] + (thr - c[i-1])*(a[i] - a[i-1])/(c[i] - c[i-1]))


def cycle_damping(alpha_deg: np.ndarray, cm: np.ndarray, amp_deg: float) -> float:
    """Xi = -(1/(pi*alpha_1^2)) * closed integral of C_M d(alpha)."""
    a = np.radians(np.asarray(alpha_deg, float))
    c = np.asarray(cm, float)
    if a[0] != a[-1] or c[0] != c[-1]:              # close the loop
        a = np.append(a, a[0])
        c = np.append(c, c[0])
    a1 = np.radians(amp_deg)
    return float(-_trapz(c, a)/(np.pi*a1*a1))


def loop_is_closed(alpha_deg: np.ndarray) -> bool:
    """True if the traced points cover both strokes: each stroke spans at
    least CLOSED_FRACTION of the incidence range of the whole record. A record
    of one stroke has no enclosed area, so no cycle damping can be taken from
    it."""
    a = np.asarray(alpha_deg, float)
    s = stroke_split(a)
    full = float(np.ptp(a))
    return all((s == st).sum() >= 3 and float(np.ptp(a[s == st])) >= CLOSED_FRACTION*full for st in ("up", "down"))


def _exp_upstroke(a, y):
    s = stroke_split(a)
    m = s == "up"
    o = np.argsort(a[m])
    return a[m][o], y[m][o]


# --------------------------------------------------------------------------- #
#  one frame
# --------------------------------------------------------------------------- #
def score_frame(fr: dict, out: dict, airfoil: str = "NACA 0012", set_name: str = "") -> dict:
    """Score one model solution `out` (the dict dsmodel.solve returns, or
    any dict with alpha_deg, alpha_dot, CL, CM, CD over one closed cycle)
    against one measured frame `fr` (from load_frame)."""
    a = np.asarray(out["alpha_deg"], float)
    ad = np.asarray(out["alpha_dot"], float)
    CL, CM, CD = (np.asarray(out[k], float) for k in ("CL", "CM", "CD"))
    rms = {q: loop_rms(fr["a"+q], fr[q], a, ad, y)
           for q, y in (("cl", CL), ("cm", CM), ("cd", CD))}
    rng = {q: float(np.ptp(fr[q])) for q in ("cl", "cm", "cd")}
    up = ad > 0
    o = np.argsort(a[up])
    ms_mod = moment_stall_alpha(a[up][o], CM[up][o])
    ms_exp = moment_stall_alpha(*_exp_upstroke(fr["acm"], fr["cm"]))
    xi_mod = cycle_damping(a, CM, fr["da"])
    closed = loop_is_closed(fr["acm"])
    xi_exp = cycle_damping(fr["acm"], fr["cm"], fr["da"]) if closed else float("nan")
    iM, iE = int(np.argmax(CL)), int(np.argmax(fr["cl"]))
    # Scores are returned unrounded: the calibration differentiates them.
    return dict(
        frame=fr["frame"], airfoil=airfoil, set=set_name,
        M=round(fr["M"], 3), k=round(fr["k"], 3),
        alpha0_deg=round(fr["a0"], 2), amp_deg=round(fr["da"], 2),
        RMS_CL=rms["cl"], RMS_CM=rms["cm"], RMS_CD=rms["cd"],
        range_CL=rng["cl"], range_CM=rng["cm"], range_CD=rng["cd"],
        nRMS_CL=rms["cl"]/rng["cl"], nRMS_CM=rms["cm"]/rng["cm"], nRMS_CD=rms["cd"]/rng["cd"],
        CLmax_model=float(CL[iM]), CLmax_exp=float(fr["cl"][iE]),
        CLmax_err_pct=float(100.0*(CL[iM] - fr["cl"][iE])/fr["cl"][iE]),
        alpha_CLmax_model=float(a[iM]), alpha_CLmax_exp=float(fr["acl"][iE]),
        dalpha_CLmax_deg=float(a[iM] - fr["acl"][iE]),
        alpha_Mstall_model=ms_mod, alpha_Mstall_exp=ms_exp, dalpha_Mstall_deg=ms_mod - ms_exp,
        Mstall_model_found=bool(np.isfinite(ms_mod)), Mstall_exp_found=bool(np.isfinite(ms_exp)),
        CMmin_model=float(CM.min()), CMmin_exp=float(fr["cm"].min()),
        CDmax_model=float(CD.max()), CDmax_exp=float(fr["cd"].max()),
        Xi_model=xi_mod, Xi_exp=xi_exp, moment_loop_closed=closed,
        Xi_sign_agree=float(np.sign(xi_mod) == np.sign(xi_exp)) if closed else float("nan"),
        CD_neg_frac=float(np.mean(CD[:-1] < 0.0)), CD_min_model=float(CD.min()))


# --------------------------------------------------------------------------- #
#  aggregates
# --------------------------------------------------------------------------- #
#  what is averaged, and how: "abs" means the mean of the absolute value
AGGREGATES = [("nRMS_CL", "mean"), ("nRMS_CM", "mean"), ("nRMS_CD", "mean"),
              ("RMS_CL", "mean"), ("RMS_CM", "mean"), ("RMS_CD", "mean"),
              ("CLmax_err_pct", "abs"), ("dalpha_CLmax_deg", "abs"),
              ("dalpha_Mstall_deg", "abs"), ("Xi_sign_agree", "mean"),
              ("CD_neg_frac", "mean")]


def _boot(x, rng):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return float("nan"), float("nan"), float("nan"), 0
    if len(x) == 1:
        return float(x[0]), float("nan"), float("nan"), 1
    m = x[rng.integers(0, len(x), size=(BOOT_N, len(x)))].mean(axis=1)
    lo, hi = np.percentile(m, [2.5, 97.5])
    return float(x.mean()), float(lo), float(hi), len(x)


def summarise(frames_df: pd.DataFrame, groups: dict) -> pd.DataFrame:
    """`groups` maps a label to a boolean mask over frames_df. Returns one row
    per (group, metric): mean, bootstrap 95 % interval over frames, n."""
    rng = np.random.default_rng(BOOT_SEED)
    rows = []
    for label, mask in groups.items():
        sub = frames_df[mask]
        for col, how in AGGREGATES:
            v = sub[col].astype(float)
            v = v.abs() if how == "abs" else v
            mean, lo, hi, n = _boot(v.values, rng)
            rows.append([label, ("mean_abs_" if how == "abs" else "mean_") + col,
                         round(mean, 4), round(lo, 4), round(hi, 4), n])
    return pd.DataFrame(rows, columns=["group", "metric", "value", "ci95_lo", "ci95_hi",
                                       "n_frames"])


# --------------------------------------------------------------------------- #
#  the verified data set
# --------------------------------------------------------------------------- #
MANIFEST = DATA/"data_manifest.csv"
FRAME_CACHE = CACHE/"frames"       # filled by python3 -m unistall.fetch_frames
QUASI_STEADY_K = 0.005                    # below this a frame is a static sweep


def manifest_frames(airfoil: str = "NACA 0012", sets: tuple = ("calibration", "held_out"),
                    dynamic: bool = True) -> pd.DataFrame:
    """Rows of data_manifest.csv that this study may use: verified against the
    NASA report's tables, not boundary-layer-trip cases, peaking inside the
    static lift range. Each frame is then run at ITS OWN measured conditions
    (M, k, mean angle and amplitude as stored in the frame), never at nominal
    ones."""
    m = pd.read_csv(MANIFEST)
    keep = (m.airfoil == airfoil) & m.use_in_study & m["set"].isin(sets)
    keep &= (m.k >= QUASI_STEADY_K) if dynamic else (m.k < QUASI_STEADY_K)
    return m[keep].reset_index(drop=True)


def score_manifest(run: object, airfoil: str = "NACA 0012", sets: tuple = ("calibration", "held_out"),
                   dynamic: bool = True) -> pd.DataFrame:
    """Score a model on every usable frame of the manifest. `run(fr)` returns
    the model solution for one loaded frame."""
    rows = []
    for r in manifest_frames(airfoil, sets, dynamic).itertuples():
        p = FRAME_CACHE/f"{r.frame}.mat"
        if not p.exists():
            raise FileNotFoundError(f"{p} missing: run python3 -m unistall.fetch_frames first")
        fr = load_frame(p)
        row = score_frame(fr, run(fr), airfoil=r.airfoil, set_name=r.set)
        row["beyond_static_moment_drag_range"] = bool(r.beyond_static_moment_drag_range)
        rows.append(row)
    return pd.DataFrame(rows, columns=FRAME_COLUMNS + ["beyond_static_moment_drag_range"])
