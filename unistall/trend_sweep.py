# -*- coding: utf-8 -*-
"""
unistall / trend_sweep.py
---------------------
Trends of the calibrated model with reduced frequency and Mach number, at the
reporting resolution (720 steps per cycle, sixth cycle), written to
results/trend_sweep.csv.

Author: Akosa Samuel Onyejekwe (independent)

Every row carries `inside_span_of_measured_loops`, which is True only if BOTH the
reduced frequency and the Mach number lie inside the ranges spanned by the
verified frames the model was calibrated and scored on (read from
data/data_manifest.csv, not typed here) and the peak incidence is inside the
static lift range. Points outside are computed and flagged, never hidden.
"""
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
from unistall import dsmodel as dm
from unistall import metrics as mt
from unistall.static_model import ALPHA_MAX
from unistall.paths import RESULTS

MEANS_AMPS = ((10.0, 5.0), (10.0, 10.0), (15.0, 5.0))
K_VALUES = (0.01, 0.025, 0.05, 0.075, 0.10, 0.125, 0.15, 0.20, 0.25, 0.30)
MACH_VALUES = (0.18, 0.22, 0.25, 0.28, 0.30)


def compared_ranges():
    m = mt.manifest_frames()
    m = m[m.within_static_mach_range]
    return (float(m.k.min()), float(m.k.max())), (float(m.M.min()), float(m.M.max()))


def main():
    (k0, k1), (M0, M1) = compared_ranges()
    rows = []
    for a0, a1 in MEANS_AMPS:
        for M in MACH_VALUES:
            for k in K_VALUES:
                o = dm.solve(a0, a1, k, M)
                i = int(np.argmax(o["CL"]))
                rows.append(dict(alpha_mean_deg=a0, alpha_amp_deg=a1, M=M, k=k, steps_per_cycle=o["meta"]["steps_per_cycle"], n_cycles=o["meta"]["cycles"],
                                 CL_max=round(float(o["CL"][i]), 4), alpha_at_CL_max_deg=round(float(o["alpha_deg"][i]), 2),
                                 CM_min=round(float(o["CM"].min()), 4), CD_max=round(float(o["CD"].max()), 4),
                                 onset_alpha_deg=round(o["onset_alpha_deg"], 2),
                                 cycle_damping=round(mt.cycle_damping(o["alpha_deg"], o["CM"], a1), 4),
                                 k_in_range=bool(k0 <= k <= k1), M_in_range=bool(M0 - 1e-9 <= M <= M1 + 1e-9),
                                 peak_in_static_range=bool(a0 + a1 <= ALPHA_MAX),
                                 inside_span_of_measured_loops=bool(k0 <= k <= k1 and M0 - 1e-9 <= M <= M1 + 1e-9
                                                             and a0 + a1 <= ALPHA_MAX)))
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS/"trend_sweep.csv", index=False)
    print(f"[trend] {len(df)} points; span of the measured loops: k {k0}-{k1}, Mach {M0}-{M1}; "
          f"{int(df.inside_span_of_measured_loops.sum())} inside, {int((~df.inside_span_of_measured_loops).sum())} flagged outside")


if __name__ == "__main__":
    main()
