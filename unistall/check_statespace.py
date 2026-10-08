# -*- coding: utf-8 -*-
"""
unistall / check_statespace.py
----------------------------
The state-space form (statespace.py) against the indicial march (dsmodel.py)
on every calibration loop, those at low Mach number included: the largest
difference in each load over the cycle, and the difference in each normalised
loop error against the measurement.

Author: Akosa Samuel Onyejekwe (independent)

Output: results/statespace_check.csv (one row per loop) and
        results/statespace_summary.csv (worst case against the target in
        data/targets.json, group "numerical").
Usage:  python3 -m unistall.check_statespace
"""
import json
import time

import numpy as np
import pandas as pd

from unistall import attached_flow as af
from unistall import dsmodel as dm
from unistall import metrics as mt
from unistall import statespace as ss
from unistall.paths import DATA, RESULTS

LOADS = ("CL", "CM", "CD")


def compare(fr: dict, consts: dict) -> dict:
    """Both forms on one loaded frame, at the frame's own conditions."""
    t0 = time.process_time()
    a = dm.solve(fr["a0"], fr["da"], fr["k"], fr["M"], consts=consts)
    t1 = time.process_time()
    b = ss.solve(fr["a0"], fr["da"], fr["k"], fr["M"], consts=consts)
    t2 = time.process_time()
    sa, sb = mt.score_frame(fr, a), mt.score_frame(fr, b)
    return {"frame": fr["frame"], "M": fr["M"], "k": fr["k"], "compressible_weight": af.compressible_weight(fr["M"]),
            "alpha_mean_deg": float(fr["a0"]),
            "alpha_amp_deg": float(fr["da"]), "steps_per_cycle": a["meta"]["steps_per_cycle"],
            **{f"max_abs_d{q}": float(np.max(np.abs(a[q] - b[q]))) for q in LOADS},
            **{f"range_{q}": float(np.ptp(a[q])) for q in LOADS},
            **{f"nRMS_{q}_indicial": sa[f"nRMS_{q}"] for q in LOADS},
            **{f"nRMS_{q}_state_space": sb[f"nRMS_{q}"] for q in LOADS},
            "damping_indicial": sa["Xi_model"], "damping_state_space": sb["Xi_model"],
            "cpu_s_indicial": t1 - t0, "cpu_s_state_space": t2 - t1}


def main() -> pd.DataFrame:
    """Compare on every calibration loop and write the two result files."""
    target = json.load(open(DATA/"targets.json"))["numerical"]["state_space_loop_error_difference_max"]
    m = mt.manifest_frames(sets=("calibration",))          # every calibration loop, low Mach number included
    consts = dm.load_constants()
    df = pd.DataFrame([compare(mt.load_frame(mt.FRAME_CACHE/f"{f}.mat"), consts) for f in m.frame])
    df.round(7).to_csv(RESULTS/"statespace_check.csv", index=False)
    diff = pd.concat([(df[f"nRMS_{q}_indicial"] - df[f"nRMS_{q}_state_space"]).abs() for q in LOADS], axis=1)
    worst = float(diff.max().max())
    summ = pd.DataFrame([dict(
        n_frames=len(df), n_frames_below_mach_0_20=int((df.compressible_weight < 1.0).sum()),
        lowest_mach=round(float(df.M.min()), 3), states=ss.N_STATES, integrator="classical fourth-order Runge-Kutta, fixed step",
        worst_difference_in_loop_error=round(worst, 6), mean_difference_in_loop_error=round(float(diff.mean().mean()), 6),
        **{f"worst_difference_in_{q}": round(float(df[f"max_abs_d{q}"].max()), 6) for q in LOADS},
        **{f"worst_difference_in_{q}_pct_of_range": round(float((100*df[f"max_abs_d{q}"]/df[f"range_{q}"]).max()), 4)
           for q in LOADS},
        worst_difference_in_damping=round(float((df.damping_indicial - df.damping_state_space).abs().max()), 6),
        cpu_ratio_state_space_to_indicial=round(float(df.cpu_s_state_space.sum()/df.cpu_s_indicial.sum()), 2),
        target_difference_in_loop_error=target, met=bool(worst <= target))])
    summ.to_csv(RESULTS/"statespace_summary.csv", index=False)
    print(summ.T.to_string(header=False))
    return summ


if __name__ == "__main__":
    main()
