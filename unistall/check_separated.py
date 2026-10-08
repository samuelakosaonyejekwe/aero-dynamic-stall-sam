# -*- coding: utf-8 -*-
"""
unistall / check_separated.py
-------------------------
The separated-flow acceptance measures of work package 3, on the CALIBRATION
frames inside the static Mach range, written to results/separated_checks.csv:

  * fraction of the cycle with model C_D < 0, and the minimum C_D
  * shift in stall-onset incidence between the default step and one four times smaller
  * number of vortex sheddings per cycle

Author: Akosa Samuel Onyejekwe (independent)
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
from unistall import dsmodel as dm
from unistall import metrics as mt
from unistall.paths import DATA, RESULTS

T = json.load(open(DATA/"targets.json"))["separated_flow"]


def main():
    m = mt.manifest_frames(sets=("calibration",))
    m = m[m.within_static_mach_range]
    consts = dm.load_constants()
    rows = []
    for r in m.itertuples():
        fr = mt.load_frame(mt.FRAME_CACHE/f"{r.frame}.mat")
        o = dm.solve(fr["a0"], fr["da"], fr["k"], fr["M"], consts=consts)
        fine = dm.solve(fr["a0"], fr["da"], fr["k"], fr["M"], consts=consts,
                        n_per_cycle=4*dm.march_resolution(fr["k"])[0])
        rows.append(dict(frame=r.frame, M=r.M, k=r.k, alpha0_deg=r.alpha0_deg, amp_deg=r.amp_deg,
                         CD_neg_frac=round(float(np.mean(o["CD"][:-1] < 0.0)), 4),
                         CD_min=round(float(o["CD"].min()), 4),
                         onset_default_deg=round(o["onset_alpha_deg"], 4),
                         onset_fine_deg=round(fine["onset_alpha_deg"], 4),
                         onset_shift_deg=round(abs(o["onset_alpha_deg"] - fine["onset_alpha_deg"]), 4),
                         sheddings_per_cycle=o["n_sheddings"],
                         repeated_sheddings_per_cycle=o["n_repeated_sheddings"]))
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS/"separated_checks.csv", index=False)
    print(f"[separated] {len(df)} calibration frames: worst C_D<0 fraction {df.CD_neg_frac.max():.3f} "
          f"(target {T['CD_negative_fraction_of_cycle_max']}), lowest C_D {df.CD_min.min():+.4f} "
          f"(target {T['CD_min']}), largest onset shift {df.onset_shift_deg.max():.4f} deg "
          f"(target {T['onset_shift_with_4x_steps_max_deg']}), "
          f"frames with more than one shedding {int((df.sheddings_per_cycle > 1).sum())}, "
          f"most sheddings in a cycle {int(df.sheddings_per_cycle.max())}")


if __name__ == "__main__":
    main()
