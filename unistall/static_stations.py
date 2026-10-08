# -*- coding: utf-8 -*-
"""
unistall / static_stations.py
-------------------------
The static description of the section at the two Mach numbers where a
quasi-steady sweep was measured in the calibration set: lift-curve slope,
zero-lift incidence, maximum static lift, the incidence where it occurs and
the normal force there. unistall/static_model.py reads the table this writes.

Author: Akosa Samuel Onyejekwe (independent)

RULE, the same for both sweeps (up-stroke only):
  slope, zero-lift incidence   least-squares line through C_L between
                               FIT_RANGE_DEG[0] and FIT_RANGE_DEG[1] degrees
  maximum static lift          the largest C_L of the up-stroke; its incidence
                               is the static-stall incidence
  C_N1                         the normal force at that incidence, from the
                               sweep's own C_L and pressure C_D

Output: results/static_by_mach_naca0012.csv
Usage:  python3 -m unistall.static_stations
"""
import numpy as np
import pandas as pd

from unistall import metrics as mt
from unistall.strokes import stroke_split
from unistall.paths import RESULTS

SWEEPS = ("frame_13308", "frame_12102")          # quasi-steady sweeps (k = 0.001) of the calibration set
FIT_RANGE_DEG = (-5.0, 8.0)


def station(name: str) -> dict:
    """The static description from one quasi-steady sweep frame `name`."""
    fr = mt.load_frame(mt.FRAME_CACHE/f"{name}.mat")
    up = stroke_split(fr["acl"]) == "up"
    a, cl = fr["acl"][up], fr["cl"][up]
    lin = (a >= FIT_RANGE_DEG[0]) & (a <= FIT_RANGE_DEG[1])
    slope, intercept = np.polyfit(a[lin], cl[lin], 1)
    i = int(np.argmax(cl))
    upd = stroke_split(fr["acd"]) == "up"                   # drag of the up sweep, as the lift
    o = np.argsort(fr["acd"][upd])
    cd = float(np.interp(a[i], fr["acd"][upd][o], fr["cd"][upd][o]))
    cn = cl[i]*np.cos(np.radians(a[i])) + cd*np.sin(np.radians(a[i]))
    return dict(M=round(fr["M"], 3), frame=name, set="calibration", CL_alpha_per_deg=round(float(slope), 4),
                CN_alpha_per_rad=round(float(np.degrees(slope)), 3), alpha_zero_lift_deg=round(float(-intercept/slope), 2),
                CL_max=round(float(cl[i]), 3), alpha_at_CL_max_deg=round(float(a[i]), 1), CN_at_CL_max=round(float(cn), 3),
                fit_range_deg=f"{FIT_RANGE_DEG[0]:.0f} to {FIT_RANGE_DEG[1]:.0f} (up sweep)",
                source=f"NASA TM-84245 quasi-steady sweep {name} (k = {fr['k']:.3f}), not wall-corrected, "
                       "same measurement chain as the dynamic frames")


if __name__ == "__main__":
    df = pd.DataFrame([station(n) for n in SWEEPS]).sort_values("M")
    print(df.drop(columns=["source"]).to_string(index=False))
    df.to_csv(RESULTS/"static_by_mach_naca0012.csv", index=False)
