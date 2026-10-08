# -*- coding: utf-8 -*-
"""
unistall / published_reference.py
-----------------------------
An outside reference for the stall model: the Leishman-Beddoes model as its
authors published it, on two NACA 0012 cases of the same experiment.

Author: Akosa Samuel Onyejekwe (independent)

Leishman and Crouse (1989) plot their model's normal force and moment against
the measurements of McAlister et al. for Mach 0.3 and reduced frequency 0.1,
at 10 +/- 10 deg and at 15 +/- 10 deg. Digitised tracings of those curves,
model and measurement, are in the third-party repository the measured loops
come from (folder "Other Data"), which also states the conditions. They are
fetched at the same fixed commit into the git-ignored cache and checksummed
(data/published_curves_inventory.csv); none is redistributed here.

For each case this script scores, with the loop error of unistall/metrics.py
(RMS difference at the measured points, each on its own stroke, over the
measured range):
  the published model against the published measured points;
  this model, with its calibrated constants, against the same points;
and the RMS difference between this model and the published model over the
cycle. The tracings are of printed figures, so the comparison is as good as
the tracing; the conditions are those the source repository gives.

Neither case is new to this model: both lie at conditions of loops in the
calibration or held-out sets. The comparison shows where this implementation
stands beside the published one; it is not a further held-out test.

Output: results/published_reference.csv
Usage:  python3 -m unistall.published_reference
"""
import hashlib
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd
import scipy.io as sio

from unistall import dsmodel as dm
from unistall import metrics as mt
from unistall import paths
from unistall.fetch_frames import SOURCE_COMMIT
from unistall.paths import DATA, RESULTS
from unistall.strokes import stroke_split

CACHE = paths.CACHE/"published"
RAW = f"https://raw.githubusercontent.com/luizpancini/BL-DSM-JFS-2021/{SOURCE_COMMIT}/Other%20Data/"
INVENTORY = DATA/"published_curves_inventory.csv"
# conditions as stated in the source repository (src/functions/load_other.m)
CASES = {"A0_10_A1_10": dict(M=0.3, k=0.1, mean_deg=9.9, amp_deg=9.9),
         "A0_10_A1_15": dict(M=0.3, k=0.1, mean_deg=15.0, amp_deg=10.0)}   # the file name has the two swapped
SOURCE = "Leishman and Crouse (1989), as traced in the source repository"


def _name(quantity: str, case: str, kind: str) -> str:
    return f"{quantity}_alpha_LeishmanCrouse1989_{case}_{kind}.mat"


def fetch() -> pd.DataFrame:
    """Download the eight curve files if absent and return their inventory
    (file, bytes, SHA-256); a file whose checksum differs from the committed
    inventory raises."""
    CACHE.mkdir(parents=True, exist_ok=True)
    known = pd.read_csv(INVENTORY).set_index("file").sha256.to_dict() if INVENTORY.exists() else {}
    rows = []
    for case in CASES:
        for q in ("CN", "CM"):
            for kind in ("EXP", "BL"):
                name = _name(q, case, kind)
                dst = CACHE/name
                if not dst.exists():
                    with urllib.request.urlopen(RAW + urllib.parse.quote(name), timeout=60) as r:
                        dst.write_bytes(r.read())
                sha = hashlib.sha256(dst.read_bytes()).hexdigest()
                if name in known and known[name] != sha:
                    raise RuntimeError(f"{name} differs from the checksum recorded in {INVENTORY.name}")
                rows.append(dict(file=name, bytes=dst.stat().st_size, sha256=sha, commit=SOURCE_COMMIT))
    return pd.DataFrame(rows)


def _curve(q: str, case: str, kind: str) -> tuple:
    """(incidence [deg], coefficient) of one traced curve, in the order traced."""
    d = sio.loadmat(str(CACHE/_name(q, case, kind)))
    a = next(v for k, v in d.items() if not k.startswith("__"))
    return a[:, 0].astype(float), a[:, 1].astype(float)


def _loop_error(a_exp: np.ndarray, y_exp: np.ndarray, a_mod: np.ndarray, y_mod: np.ndarray,
                rate: np.ndarray) -> float:
    return mt.loop_rms(a_exp, y_exp, a_mod, rate, y_mod)/float(np.ptp(y_exp))


def compare() -> pd.DataFrame:
    """One row per case and coefficient: the two models against the measured
    points, and against each other."""
    consts = dm.load_constants()
    rows = []
    for case, c in CASES.items():
        o = dm.solve(c["mean_deg"], c["amp_deg"], c["k"], c["M"], consts=consts)
        for q in ("CN", "CM"):
            a_e, y_e = _curve(q, case, "EXP")
            a_p, y_p = _curve(q, case, "BL")
            rate_p = np.where(stroke_split(a_p) == "up", 1.0, -1.0)      # the traced loop, split at its turning points
            # this model read at the incidences of the published curve, each on its own stroke
            br = mt._branches(o["alpha_deg"], o["alpha_dot"], o[q])
            mine = np.array([np.interp(a, *br[s]) for a, s in zip(a_p, stroke_split(a_p), strict=True)])
            rows.append(dict(
                case=f"{c['mean_deg']:g} +/- {c['amp_deg']:g} deg, k {c['k']:g}, M {c['M']:g}", coefficient=q,
                measured_points=len(a_e), published_curve_points=len(a_p), measured_range=float(np.ptp(y_e)),
                loop_error_published_model=_loop_error(a_e, y_e, a_p, y_p, rate_p),
                loop_error_this_model=_loop_error(a_e, y_e, o["alpha_deg"], o[q], o["alpha_dot"]),
                rms_difference_between_models_over_range=float(np.sqrt(np.mean((mine - y_p)**2))/np.ptp(y_e)),
                extreme_published_model=float(y_p.max() if q == "CN" else y_p.min()),
                extreme_this_model=float(o[q].max() if q == "CN" else o[q].min()),
                extreme_measured=float(y_e.max() if q == "CN" else y_e.min()), source=SOURCE))
    return pd.DataFrame(rows)


if __name__ == "__main__":
    fetch().to_csv(INVENTORY, index=False)
    df = compare()
    df.round(5).to_csv(RESULTS/"published_reference.csv", index=False)
    pd.set_option("display.width", 220)
    print(df.drop(columns=["source"]).round(4).to_string(index=False))
