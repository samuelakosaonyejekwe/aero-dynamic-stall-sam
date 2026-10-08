# -*- coding: utf-8 -*-
"""
unistall / convergence.py
---------------------
Independence of the results from the two numerical choices of the march: the
number of time steps per cycle and the number of cycles marched before the
reported one. The model has no spatial mesh; these are its discretisation.

Author: Akosa Samuel Onyejekwe (independent)

Every calibration loop inside the static Mach range is marched at each
multiple in STEP_FACTORS of the default step count, and with each number in
EXTRA_CYCLES of cycles added to the default (dsmodel.march_resolution). Loads are compared with the finest march on a common
grid of phase; scores are the normalised loop errors against the measurement.

The step study is judged against the limit of zero step, by Richardson
extrapolation (results/convergence_extrapolated.csv), not against the finest
march run.

Outputs: results/convergence_steps.csv, results/convergence_cycles.csv,
         results/convergence_summary.csv (worst case against the targets in
         data/targets.json, group "numerical").
"""
import json
import numpy as np
import pandas as pd

from unistall import dsmodel as dm
from unistall import metrics as mt
from unistall.paths import DATA, RESULTS

STEP_FACTORS = (0.125, 0.25, 0.5, 1.0, 2.0, 4.0)     # multiples of the default steps per cycle
EXTRA_CYCLES = (0, 2, 4)                             # cycles added to the default
PHASE = np.linspace(0.0, 1.0, 361)
LOADS = ("CL", "CM", "CD")


def _on_phase(o: dict) -> dict:
    """The last cycle of a solution on the common phase grid."""
    ph = np.linspace(0.0, 1.0, len(o["CL"]))
    return {q: np.interp(PHASE, ph, o[q]) for q in LOADS}


def _row(fr: dict, o: dict, fine: dict, fine_score: dict, **tags: object) -> dict:
    cur, s = _on_phase(o), mt.score_frame(fr, o)
    return {"frame": fr["frame"], "k": fr["k"], **tags,
            **{f"max_abs_d{q}": float(np.max(np.abs(cur[q] - fine[q]))) for q in LOADS},
            **{f"d_nRMS_{q}": abs(s[f"nRMS_{q}"] - fine_score[f"nRMS_{q}"]) for q in LOADS},
            **{f"nRMS_{q}": s[f"nRMS_{q}"] for q in LOADS}}


def _studies(frames: list) -> tuple:
    """March every frame at each multiple of the default step count and with
    each number of extra cycles; differences are from the finest march."""
    consts = dm.load_constants()
    steps, cycles = [], []
    for fr in frames:
        n0, c0 = dm.march_resolution(fr["k"])

        def run(n: int, c: int, fr: dict = fr) -> dict:
            return dm.solve(fr["a0"], fr["da"], fr["k"], fr["M"], consts=consts, n_per_cycle=n, n_cycles=c)
        runs = [(f, run(int(round(f*n0)), c0)) for f in STEP_FACTORS]
        fine = runs[-1][1]
        steps += [_row(fr, o, _on_phase(fine), mt.score_frame(fr, fine), step_factor=f, steps_per_cycle=int(round(f*n0)),
                       step_semichords=round(2.0*np.pi/(fr["k"]*f*n0), 5)) for f, o in runs]
        runs = [(e, run(n0, c0 + e)) for e in EXTRA_CYCLES]
        fine = runs[-1][1]
        cycles += [_row(fr, o, _on_phase(fine), mt.score_frame(fr, fine), extra_cycles=e, cycles=c0 + e) for e, o in runs]
    return pd.DataFrame(steps), pd.DataFrame(cycles)


def extrapolated_step_error(steps: pd.DataFrame) -> pd.DataFrame:
    """Error of each loop error at the default step against its value at zero
    step, by Richardson extrapolation from the default step and the two finer
    marches (2 and 4 times as many steps).

    With s1, s2, s4 the loop error at those three resolutions, the observed
    order is p = log2((s1 - s2)/(s2 - s4)) and the limit is
    s4 + (s4 - s2)/(2^p - 1). Where the three values are not monotone or p is
    outside 0.5 to 3 (a switching event moved by a step), first order is
    assumed, which gives the limit s4 + (s4 - s2). One row per loop and
    coefficient."""
    rows = []
    for frame, d in steps.groupby("frame", sort=False):
        d = d.set_index("step_factor")
        for q in LOADS:
            s1, s2, s4 = (float(d.loc[f, f"nRMS_{q}"]) for f in (1.0, 2.0, 4.0))
            a, b = s1 - s2, s2 - s4
            p = float(np.log2(a/b)) if a*b > 0.0 else float("nan")
            observed = bool(np.isfinite(p) and 0.5 <= p <= 3.0)
            limit = s4 + b/(2.0**p - 1.0) if observed else s4 + b
            rows.append(dict(frame=frame, k=float(d.k.iloc[0]), coefficient=q,
                             step_semichords=float(d.loc[1.0, "step_semichords"]),
                             loop_error_default=s1, loop_error_2x=s2, loop_error_4x=s4,
                             observed_order=p if observed else float("nan"),
                             order_used="observed" if observed else "first order assumed",
                             loop_error_at_zero_step=limit, error_of_default_step=abs(s1 - limit)))
    return pd.DataFrame(rows)


def main() -> pd.DataFrame:
    """Run both studies and compare the default resolution with the targets."""
    T = json.load(open(DATA/"targets.json"))["numerical"]
    m = mt.manifest_frames(sets=("calibration",))
    m = m[m.within_static_mach_range]
    frames = [mt.load_frame(mt.FRAME_CACHE/f"{f}.mat") for f in m.frame]
    steps, cycles = _studies(frames)
    steps.round(6).to_csv(RESULTS/"convergence_steps.csv", index=False)
    cycles.round(6).to_csv(RESULTS/"convergence_cycles.csv", index=False)
    extra = extrapolated_step_error(steps)
    extra.round(6).to_csv(RESULTS/"convergence_extrapolated.csv", index=False)
    cols = [f"d_nRMS_{q}" for q in LOADS]
    rows = []
    for name, df, var, finest in (("step size", steps, "step_factor", f"{STEP_FACTORS[-1]:g} times as many steps"),
                                  ("cycles marched", cycles, "extra_cycles", f"{EXTRA_CYCLES[-1]} more cycles")):
        worst = df[df[var] == (1.0 if var == "step_factor" else 0)]
        rows.append(dict(study=name, default_compared_with=finest, n_frames=worst.frame.nunique(),
                         worst_change_in_loop_error=round(float(worst[cols].max().max()), 6),
                         worst_change_in_CL=round(float(worst.max_abs_dCL.max()), 6),
                         worst_change_in_CM=round(float(worst.max_abs_dCM.max()), 6),
                         worst_change_in_CD=round(float(worst.max_abs_dCD.max()), 6),
                         worst_error_against_zero_step=(round(float(extra.error_of_default_step.max()), 6)
                                                        if var == "step_factor" else float("nan")),
                         target_change_in_loop_error=T["loop_error_change_max"],
                         met=bool((extra.error_of_default_step.max() if var == "step_factor"
                                   else worst[cols].max().max()) <= T["loop_error_change_max"])))
    summ = pd.DataFrame(rows)
    summ.to_csv(RESULTS/"convergence_summary.csv", index=False)
    print(summ.to_string(index=False))
    return summ


if __name__ == "__main__":
    main()
