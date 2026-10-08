# -*- coding: utf-8 -*-
# Run from the repository root:  PYTHONPATH=. python3 04_solver/unistall_solver.py 10 10 0.10 0.30 --out cycle.csv
"""
04_solver / unistall_solver.py
------------------------------
A command-line entry to the packaged solver for ONE condition. It holds no
physics: the load model is the package unistall (unistall.dsmodel for the
indicial form, unistall.statespace for the state-space form), and this
script only calls it and writes the last cycle to a CSV file.

Author: Akosa Samuel Onyejekwe (independent)

    PYTHONPATH=. python3 04_solver/unistall_solver.py MEAN AMP K MACH [--form indicial|state-space] [--out FILE]

MEAN and AMP are the mean incidence and the pitch amplitude in degrees, K
the reduced frequency omega c / (2 U) and MACH the stream Mach number. The
constants are the calibrated ones, the chord that of the tunnel model and
the speed of sound that of the case study's air; none of the three changes
the coefficients. Columns written: phase_deg, alpha_deg, CL, CD, CM_c4,
f_separation, CN_vortex.
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

import project_meta as pm
from unistall import dsmodel, statespace

FORMS = {"indicial": dsmodel.solve, "state-space": statespace.solve}
COLUMNS = (("alpha_deg", "alpha_deg"), ("CL", "CL"), ("CD", "CD"), ("CM_c4", "CM"), ("f_separation", "f_sep"),
           ("CN_vortex", "CN_vortex"))
DECIMALS = 6


def cycle_table(mean_deg: float, amp_deg: float, k: float, mach: float, form: str = "indicial") -> pd.DataFrame:
    """The last cycle of the packaged solver for one condition, as a table:
    cycle phase [deg], incidence [deg], lift, drag, moment about the quarter
    chord, separation point and vortex normal force. `form` is "indicial" or
    "state-space"."""
    out = FORMS[form](mean_deg, amp_deg, k, mach, a_sound=pm.speed_of_sound())
    table = pd.DataFrame({"phase_deg": 360.0*out["t"]/out["t"][-1]})
    for column, key in COLUMNS:
        table[column] = out[key]
    return table.round(DECIMALS)


def main(argv: list | None = None) -> int:
    """Read the arguments, run the solver and write the cycle."""
    parser = argparse.ArgumentParser(description="Run the packaged UNISTALL load model for one condition.")
    parser.add_argument("mean", type=float, help="mean incidence [deg]")
    parser.add_argument("amp", type=float, help="pitch amplitude [deg]")
    parser.add_argument("k", type=float, help="reduced frequency omega c / (2 U)")
    parser.add_argument("mach", type=float, help="stream Mach number")
    parser.add_argument("--form", choices=sorted(FORMS), default="indicial", help="form of the model")
    parser.add_argument("--out", type=Path, default=Path("unistall_cycle.csv"), help="CSV file to write")
    args = parser.parse_args(argv)
    table = cycle_table(args.mean, args.amp, args.k, args.mach, args.form)
    table.to_csv(args.out, index=False)
    print("[solver] %s form: CL_max %.4f, CM_min %.4f, CD_max %.4f; %d steps written to %s"
          % (args.form, table["CL"].max(), table["CM_c4"].min(), table["CD"].max(), len(table), args.out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
