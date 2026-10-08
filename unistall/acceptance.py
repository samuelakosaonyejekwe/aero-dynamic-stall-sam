# -*- coding: utf-8 -*-
"""
unistall / acceptance.py
--------------------
A record of the held-out measures as the validation gave them, so that a test
can confirm the stored results are what the code and data reproduce.

Author: Akosa Samuel Onyejekwe (independent)

`python3 -m unistall.acceptance --record` writes results/acceptance_thresholds.json
from the current results/validation_targets.csv. tests/test_acceptance.py
recomputes the measures from the model and the fetched frames and requires
each to equal its recorded value within TOLERANCE. It is a check that the
result files are current and reproducible, not a limit on how the results may
move: the held-out scores are not used to steer the model.
"""
import sys, json
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parent
from unistall.paths import RESULTS
TOLERANCE = 2.0e-4          # recomputed against recorded value; the record is rounded to four decimals


def record():
    """Write the record of the held-out measures and return it."""
    from unistall.paths import DATA
    t = pd.read_csv(RESULTS/"validation_targets.csv").set_index("measure")
    group = json.load(open(DATA/"targets.json"))["held_out_validation"]["group"]
    out = dict(tolerance=TOLERANCE, group=group,
               use="reproduction check only: each recomputed measure must equal its recorded value within the tolerance",
               measures={m: dict(recorded=float(r.measured), target=float(r.target), kind=r.kind, met=bool(r.met))
                         for m, r in t.iterrows()})
    json.dump(out, open(RESULTS/"acceptance_thresholds.json", "w"), indent=2)
    open(RESULTS/"acceptance_thresholds.json", "a").write("\n")
    return out


if __name__ == "__main__":
    if "--record" in sys.argv:
        print(json.dumps(record()["measures"], indent=1))
    else:
        sys.exit(__doc__)
