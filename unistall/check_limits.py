# -*- coding: utf-8 -*-
"""
unistall / check_limits.py
----------------------
Checks of the load model at its limits, each a number in a file: what the
model does where no measured loop tests it directly.

Author: Akosa Samuel Onyejekwe (independent)

  MIRROR PAIR       a loop and its mirror image about the zero-lift incidence:
                    normal force and moment must change sign exactly; lift and
                    drag are resolved on the experiment's incidence scale and
                    differ by the rotation through the zero-lift incidence.
  STATIC LIMIT      the drag of a slow sweep against the static drag data.
  ZERO LIFT         the chord force and drag on the loops that cross zero lift.
  VORTEX SIGN RULE  how often the vortex lift is set to zero because it
                    opposes the separated force, and the size of the step.
  REPEATED SHEDDING loops on which the vortex clock restarts within a cycle.

Outputs (results/): limit_mirror_pair.csv, limit_static_drag.csv,
limit_zero_lift.csv, limit_vortex_rules.csv.

Usage:  python3 -m unistall.check_limits
"""
import numpy as np
import pandas as pd

from unistall import dsmodel as dm
from unistall import metrics as mt
from unistall.paths import DATA, RESULTS
from unistall.static_model import StaticModel

MIRROR_CASES = ((10.0, 10.0, 0.10, 0.30), (15.0, 10.0, 0.10, 0.29), (9.0, 5.0, 0.10, 0.30), (12.0, 8.0, 0.05, 0.22))
SWEEP = dict(mean=7.5, amp=12.5, k=0.002, M=0.30)      # slow sweep, -5 to 20 deg


def mirror_pair() -> pd.DataFrame:
    """Each case and its mirror image about the zero-lift incidence alpha_0
    (mean 2 alpha_0 - mean, amplitude reversed). Largest absolute values over
    the cycle of C_N + C_N', C_M + C_M' (zero if the model is symmetric about
    alpha_0), C_L + C_L' and C_D - C_D' (not zero: see module text)."""
    consts, rows = dm.load_constants(), []
    for mean, amp, k, M in MIRROR_CASES:
        sm = StaticModel(M)
        a = dm.solve(mean, amp, k, M, consts=consts, static=sm)
        b = dm.solve(2.0*sm.alpha0_deg - mean, -amp, k, M, consts=consts, static=sm)
        rows.append(dict(mean_deg=mean, amp_deg=amp, k=k, M=M, alpha0_deg=round(sm.alpha0_deg, 4),
                         max_abs_CN_sum=float(np.abs(a["CN"] + b["CN"]).max()),
                         max_abs_CM_sum=float(np.abs(a["CM"] + b["CM"]).max()),
                         max_abs_CC_difference=float(np.abs(a["CC"] - b["CC"]).max()),
                         max_abs_CL_sum=float(np.abs(a["CL"] + b["CL"]).max()),
                         max_abs_CD_difference=float(np.abs(a["CD"] - b["CD"]).max()),
                         CD_range=float(np.ptp(a["CD"]))))
    return pd.DataFrame(rows)


def static_drag() -> pd.DataFrame:
    """Drag of a slow sweep (up-stroke) at the incidences of the static drag
    data at Mach 0.30, with and without the constant C_D0. The data are total
    drag from a wake survey; the oscillating loops are pressure drag."""
    d = pd.read_csv(DATA/"static_naca0012_M030_CD.csv")
    o = dm.solve(SWEEP["mean"], SWEEP["amp"], SWEEP["k"], SWEEP["M"], consts=dm.load_constants(), n_cycles=2)
    up = o["alpha_dot"] > 0.0
    order = np.argsort(o["alpha_deg"][up])
    cd = np.interp(d.alpha_deg, o["alpha_deg"][up][order], o["CD"][up][order])
    cd0 = o["meta"]["constants"]["CD0"]
    return pd.DataFrame(dict(alpha_deg=d.alpha_deg, CD_static_data=d.Cd, kind=d.kind, CD_model_slow_sweep=cd.round(5),
                             CD_model_without_CD0=(cd - cd0).round(5), error=(cd - d.Cd).round(5),
                             error_without_CD0=(cd - cd0 - d.Cd).round(5), sweep_k=SWEEP["k"]))


def _frames() -> list:
    m = mt.manifest_frames()
    m = m[m.within_static_mach_range]
    return [(r, mt.load_frame(mt.FRAME_CACHE/f"{r.frame}.mat")) for r in m.itertuples()]


def loop_rules() -> tuple:
    """(zero-lift table, vortex-rule table) over every usable loop inside the
    static Mach range, with the calibrated constants."""
    consts, zero, vort = dm.load_constants(), [], []
    for r, fr in _frames():
        sm = StaticModel(fr["M"])
        o = dm.solve(fr["a0"], fr["da"], fr["k"], fr["M"], consts=consts, static=sm)
        if fr["a0"] - abs(fr["da"]) < sm.alpha0_deg:
            zero.append(dict(frame=r.frame, set=r.set, M=r.M, k=r.k, alpha0_deg=r.alpha0_deg, amp_deg=r.amp_deg,
                             CC_min=float(o["CC"].min()), CD_min=float(o["CD"].min()),
                             fraction_of_cycle_CC_negative=float(np.mean(o["CC"][:-1] < 0.0)),
                             fraction_of_cycle_CD_negative=float(np.mean(o["CD"][:-1] < 0.0))))
        rem = np.abs(o["CN_vortex_removed"][1:])
        vort.append(dict(frame=r.frame, set=r.set, M=r.M, k=r.k, alpha0_deg=r.alpha0_deg, amp_deg=r.amp_deg,
                         steps_with_vortex_lift_removed=int((rem > 0.0).sum()),
                         largest_vortex_lift_removed=float(rem.max()),
                         sheddings=o["n_sheddings"], repeated_sheddings=o["n_repeated_sheddings"]))
    return pd.DataFrame(zero), pd.DataFrame(vort)


def main() -> None:
    mp = mirror_pair()
    mp.to_csv(RESULTS/"limit_mirror_pair.csv", index=False)
    sd = static_drag()
    sd.to_csv(RESULTS/"limit_static_drag.csv", index=False)
    zero, vort = loop_rules()
    zero.round(5).to_csv(RESULTS/"limit_zero_lift.csv", index=False)
    vort.round(6).to_csv(RESULTS/"limit_vortex_rules.csv", index=False)
    pd.set_option("display.width", 220)
    print(mp.to_string(index=False))
    print(sd.drop(columns=["kind"]).to_string(index=False))
    print(zero.round(4).to_string(index=False))
    print(f"[limits] vortex lift removed by the sign rule on {int((vort.steps_with_vortex_lift_removed > 0).sum())} of "
          f"{len(vort)} loops, largest {vort.largest_vortex_lift_removed.max():.5f}; repeated shedding on "
          f"{int((vort.repeated_sheddings > 0).sum())} loops")


if __name__ == "__main__":
    main()
