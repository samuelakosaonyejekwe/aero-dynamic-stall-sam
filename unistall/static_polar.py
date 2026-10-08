# -*- coding: utf-8 -*-
"""
unistall / static_polar.py
----------------------
Assemble the static polar the load model reads, from the measured points on
file and nothing else:

  static_naca0012_M030_CL.csv   NASA TM-84245 Vol. 1 Fig. 16(a)   -5 to 25 deg
  static_naca0012_M030_CM.csv   NASA TM-84245 Vol. 1 Fig. 9       -5 to 20 deg
  static_naca0012_M030_CD.csv   NASA TM-84245 Vol. 1 Table 7 (wake survey,
                                -5 to 14 deg) and Vol. 2 Fig. 4 (pressure
                                drag, 10 to 20 deg)

Author: Akosa Samuel Onyejekwe (independent)

Output: static_polar_naca0012_M030.csv, one row per half degree from 0 to
ALPHA_MAX, with C_L, C_D, C_M and C_N on three curves.

WHY THREE CURVES. Past static stall the measurements do not lie on one curve.
Between about 14.5 and 20 deg the lift points fall on two branches, near 1.1-1.25
and near 0.8-1.0, and the moment points do the same (near zero, and near -0.08).
That is the bistable post-stall behaviour of this section, not scatter, and a
single line through it is a choice. Three are written so the choice is visible
and its effect can be measured:
    lower   the more separated branch (lower lift, more nose-down moment)
    upper   the less separated branch
    mean    the average of the two
The model follows the lower one (unistall/static_model.py, curve "more"); the
cost with each is in results/static_curve_choice.csv.

WHAT IS MEASURED AND WHAT IS NOT. Every row carries flags:
    CL_measured   True up to the last lift point (25 deg)
    CM_measured   True up to 20 deg; held at the 20-deg value beyond
    CD_measured   True up to 20 deg
    CN_basis      "CL,CD" where both are measured; "CL/cos(alpha)" beyond 20
                  deg, i.e. zero chord force, the fully separated limit. At
                  20 deg, where both forms can be evaluated, they are compared
                  and the difference is printed.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator

HERE = Path(__file__).resolve().parent
from unistall.paths import DATA, RESULTS
ALPHA_MAX = 25.0
ALPHA_STEP = 0.5
STALL_SPLIT_FROM = 14.3     # post-stall branches exist beyond this incidence
CD_CM_MEASURED_TO = 20.0
MERGE_WITHIN_DEG = 0.3      # points closer than this in incidence are one station


def _stations(a, y):
    """Average points that sit at the same nominal incidence."""
    o = np.argsort(a)
    a, y = np.asarray(a)[o], np.asarray(y)[o]
    ga, gy, i = [], [], 0
    while i < len(a):
        j = i
        while j + 1 < len(a) and a[j+1] - a[i] <= MERGE_WITHIN_DEG:
            j += 1
        ga.append(a[i:j+1].mean())
        gy.append(y[i:j+1].mean())
        i = j + 1
    return np.array(ga), np.array(gy)


def _curve(a, y):
    a, y = _stations(a, y)
    return PchipInterpolator(a, y, extrapolate=False), float(a.min()), float(a.max())


def build():
    cl = pd.read_csv(DATA/"static_naca0012_M030_CL.csv")
    cm = pd.read_csv(DATA/"static_naca0012_M030_CM.csv")
    cd = pd.read_csv(DATA/"static_naca0012_M030_CD.csv")

    # ---- lift: one curve to stall, two branches after it, rejoined at 25 deg
    pre = cl[cl.alpha_deg <= STALL_SPLIT_FROM]
    post = cl[cl.alpha_deg > STALL_SPLIT_FROM]
    last = post[post.alpha_deg > 24.0]                  # the 25-deg pair: one station
    mid = post[post.alpha_deg <= 24.0]
    # the two branches are separated by the widest gap in lift at each station;
    # 1.05 lies in that gap at every station between 14.5 and 20 deg
    CL_GAP = 1.05
    up, lo = mid[mid.Cl > CL_GAP], mid[mid.Cl <= CL_GAP]
    anchor_a, anchor_y = [pre.alpha_deg.max()], [float(pre.Cl[pre.alpha_deg.idxmax()])]
    end_a, end_y = [last.alpha_deg.mean()], [last.Cl.mean()]
    f_pre, _, _ = _curve(pre.alpha_deg, pre.Cl)
    f_lo, _, _ = _curve(np.r_[anchor_a, lo.alpha_deg, end_a], np.r_[anchor_y, lo.Cl, end_y])
    f_up, _, _ = _curve(np.r_[anchor_a, up.alpha_deg, end_a], np.r_[anchor_y, up.Cl, end_y])

    # ---- moment: same construction; the less separated branch is the one near zero
    mpre = cm[cm.alpha_deg <= STALL_SPLIT_FROM]
    mpost = cm[cm.alpha_deg > STALL_SPLIT_FROM]
    CM_GAP = -0.05
    mup, mlo = mpost[mpost.Cm_c4 > CM_GAP], mpost[mpost.Cm_c4 <= CM_GAP]
    m_anchor_a, m_anchor_y = [mpre.alpha_deg.max()], [float(mpre.Cm_c4[mpre.alpha_deg.idxmax()])]
    g_pre, _, _ = _curve(mpre.alpha_deg, mpre.Cm_c4)
    g_lo, _, lo_end = _curve(np.r_[m_anchor_a, mlo.alpha_deg], np.r_[m_anchor_y, mlo.Cm_c4])
    g_up, _, up_end = _curve(np.r_[m_anchor_a, mup.alpha_deg], np.r_[m_anchor_y, mup.Cm_c4])

    # ---- drag: tabulated wake-survey values to 14 deg, pressure drag beyond.
    #      Where both exist (10-14 deg) the wake survey is used: it is total
    #      drag and the more accurate of the two. The pressure-drag trace is
    #      taken up only beyond the last wake station by more than the merging
    #      distance, so the unstalled wake value at 14 deg and the stalled
    #      pressure value at 14.03 deg are never averaged into one station.
    #      Past stall there is one drag trace, with one excursion (14.5-15 deg)
    #      that belongs to the more separated state; it is used as measured
    #      for all three curves.
    wake = cd[cd.kind.str.startswith("total")]
    press = cd[~cd.kind.str.startswith("total")
               & (cd.alpha_deg > wake.alpha_deg.max() + MERGE_WITHIN_DEG)]
    h, _, _ = _curve(np.r_[wake.alpha_deg, press.alpha_deg], np.r_[wake.Cd, press.Cd])

    a = np.arange(0.0, ALPHA_MAX + 1e-9, ALPHA_STEP)
    out = pd.DataFrame({"alpha_deg": a})
    a_split = anchor_a[0]
    m_split = m_anchor_a[0]
    Cd = h(np.minimum(a, CD_CM_MEASURED_TO))
    cols = {}
    for name, fb, gb, gend in (("lower", f_lo, g_lo, lo_end), ("upper", f_up, g_up, up_end)):
        Cl = np.where(a <= a_split, f_pre(np.minimum(a, a_split)), fb(np.clip(a, a_split, end_a[0])))
        Cm = np.where(a <= m_split, g_pre(np.minimum(a, m_split)), gb(np.clip(a, m_split, gend)))
        cols[name] = (Cl, Cm)
    cols["mean"] = (0.5*(cols["lower"][0] + cols["upper"][0]), 0.5*(cols["lower"][1] + cols["upper"][1]))
    ar = np.radians(a)
    for name, (Cl, Cm) in cols.items():
        Cn_full = Cl*np.cos(ar) + Cd*np.sin(ar)
        Cn_sep = Cl/np.cos(ar)
        out[f"Cl_{name}"] = Cl.round(4)
        out[f"Cm_{name}"] = Cm.round(4)
        out[f"Cn_{name}"] = np.where(a <= CD_CM_MEASURED_TO, Cn_full, Cn_sep).round(4)
    out["Cd"] = Cd.round(5)
    out["CL_measured"] = a <= cl.alpha_deg.max() + 0.05
    out["CM_measured"] = a <= CD_CM_MEASURED_TO
    out["CD_measured"] = a <= CD_CM_MEASURED_TO
    out["CN_basis"] = np.where(a <= CD_CM_MEASURED_TO, "CL,CD", "CL/cos(alpha)")
    out["source"] = "NASA TM-84245: Vol. 1 Fig. 16(a), Fig. 9, Table 7; Vol. 2 Fig. 4 (see static_naca0012_M030_*.csv)"
    out.to_csv(RESULTS/"static_polar_naca0012_M030.csv", index=False)

    # the zero-chord-force form against the measured one, where both exist
    i20 = int(np.argmin(np.abs(a - CD_CM_MEASURED_TO)))
    chk = {n: float(cols[n][0][i20]/np.cos(ar[i20])
                    - (cols[n][0][i20]*np.cos(ar[i20]) + Cd[i20]*np.sin(ar[i20]))) for n in cols}
    return out, chk


if __name__ == "__main__":
    out, chk = build()
    summ = pd.read_csv(DATA/"static_naca0012_M030_summary.csv").set_index("quantity")["value"]
    i = int(out.Cl_mean.idxmax())
    lin = out[out.alpha_deg <= 8.0]
    slope = float(np.polyfit(lin.alpha_deg, lin.Cl_mean, 1)[0])
    print(f"[static] {len(out)} rows, 0 to {out.alpha_deg.max():.0f} deg; "
          f"slope {slope:.4f}/deg (report {summ['CL_alpha_per_deg']}); "
          f"C_L,max {out.Cl_mean[i]:.3f} at {out.alpha_deg[i]:.1f} deg "
          f"(report {summ['CL_max']} at {summ['alpha_static_stall_deg']})")
    print("[static] C_N from C_L/cos(alpha) minus C_N from C_L and C_D, at 20 deg: "
          + ", ".join(f"{k} {v:+.3f}" for k, v in chk.items()))
