# -*- coding: utf-8 -*-
"""
unistall / build_manifest.py
------------------------
Write data/data_manifest.csv, the list that decides which measured loops the
study may use, from the files it rests on:

  data/frames_inventory.csv          every fetched frame, its conditions as
                                     stored in the frame, its set and checksum
                                     (written by unistall/fetch_frames.py)
  data/report_tables_naca0012.csv    the conditions the NASA report's tables
  data/report_tables_other_airfoils.csv   give for each frame number
  data/report_table9_zero_drift.csv  frames the report marks for transducer
                                     zero drift

Author: Akosa Samuel Onyejekwe (independent)

A frame is "verified" if the report lists its number and the Mach number,
reduced frequency, mean angle and amplitude stored in the frame agree with a
listed row inside TOLERANCE; otherwise the differing quantities are named, or
the frame is "not listed". A frame is used in the study if it is verified, is
not a boundary-layer-trip case and peaks inside the static lift range.

Usage:  python3 -m unistall.build_manifest          write the manifest
        python3 -m unistall.build_manifest --check  compare with the committed file
"""
import sys
import numpy as np
import pandas as pd

from unistall.paths import DATA

CHORD_M = 0.61
RE_PER_MACH = 14.0e6
RE_BASIS = "Re = 14e6 x M (Vol. 1, p.9 of text / PDF p.18)"
TOLERANCE = dict(M=0.02, k=0.011, alpha0_deg=0.5, amp_deg=0.5)   # stored against tabulated conditions
QUASI_STEADY_K = 0.005                    # below this a frame is a static sweep
STATIC_LIFT_MAX_DEG, STATIC_MOMENT_DRAG_MAX_DEG = 25.0, 20.0     # extent of the static data in incidence
STATIC_MACH_RANGE = (0.199, 0.318)        # the static stations 0.215 and 0.302 with a margin of 0.016
COLUMNS = ["frame", "airfoil", "set", "M", "k", "alpha0_deg", "amp_deg", "peak_alpha_deg", "chord_m", "Re_c",
           "Re_basis", "report_location", "verification", "differences", "flags", "use_in_study", "sha256",
           "beyond_static_lift_range", "beyond_static_moment_drag_range", "within_static_mach_range"]


def _report_rows() -> pd.DataFrame:
    naca = pd.read_csv(DATA/"report_tables_naca0012.csv").assign(airfoil="NACA 0012")
    return pd.concat([naca, pd.read_csv(DATA/"report_tables_other_airfoils.csv")], ignore_index=True)


def _verify(frame: object, rows: pd.DataFrame) -> tuple:
    """(report location, verification, differences, remarks) for one inventory
    row against the report rows that carry its frame number."""
    if not len(rows):
        return "", "not listed", "", []
    where = "; ".join(dict.fromkeys(f"Table {int(r.table)} p.{int(r.report_page)}" for r in rows.itertuples()))
    best = []
    for r in rows.itertuples():
        off = [q for q, tol in TOLERANCE.items() if abs(getattr(frame, q) - getattr(r, q)) > tol]
        if not off:
            best = []
            break
        best = off if not best or len(off) < len(best) else best
    remarks = [str(v) for v in rows.remark.dropna().unique()]
    return where, "verified" if not best else "differs", "; ".join(best), remarks


def build() -> pd.DataFrame:
    """The manifest as a table, one row per fetched frame."""
    inv = pd.read_csv(DATA/"frames_inventory.csv")
    report = _report_rows()
    drift = pd.read_csv(DATA/"report_table9_zero_drift.csv").set_index("frame").flag.to_dict()
    out = []
    for f in inv.itertuples():
        where, verification, differences, remarks = _verify(f, report[report.frame_number == f.frame_number])
        flags = []
        if any("low Reynolds number" in r for r in remarks):
            flags.append(next(r for r in remarks if "low Reynolds number" in r))
        if f.k < QUASI_STEADY_K:
            flags.append("quasi-steady")
        if "trip" in remarks:
            flags.append("boundary-layer trip")
        if f.frame in drift:
            flags.append(drift[f.frame])
        beyond_lift = bool(f.peak_alpha_deg > STATIC_LIFT_MAX_DEG)
        out.append(dict(
            frame=f.frame, airfoil=f.airfoil, set=f.set, M=f.M, k=f.k, alpha0_deg=f.alpha0_deg, amp_deg=f.amp_deg,
            peak_alpha_deg=f.peak_alpha_deg, chord_m=CHORD_M, Re_c=float(round(RE_PER_MACH*f.M, -4)), Re_basis=RE_BASIS,
            report_location=where or np.nan, verification=verification, differences=differences or np.nan,
            flags="; ".join(flags) or np.nan,
            use_in_study=bool(verification == "verified" and "boundary-layer trip" not in flags and not beyond_lift),
            sha256=f.sha256, beyond_static_lift_range=beyond_lift,
            beyond_static_moment_drag_range=bool(f.peak_alpha_deg > STATIC_MOMENT_DRAG_MAX_DEG),
            within_static_mach_range=bool(STATIC_MACH_RANGE[0] <= f.M <= STATIC_MACH_RANGE[1])))
    return pd.DataFrame(out, columns=COLUMNS)


def differences_from_committed() -> list:
    """Cells in which the committed manifest differs from what `build` gives,
    as (frame, column, committed, built)."""
    new, old = build(), pd.read_csv(DATA/"data_manifest.csv")
    if list(new.frame) != list(old.frame) or list(new.columns) != list(old.columns):
        return [("(table)", "frames or columns", len(old), len(new))]
    diff = []
    for c in COLUMNS:
        a, b = old[c], new[c]
        same = (a == b) | (a.isna() & b.isna())
        if a.dtype.kind == "f" and b.dtype.kind == "f":
            same |= np.isclose(a, b, rtol=1e-12, atol=0.0)
        diff += [(old.frame[i], c, a[i], b[i]) for i in np.where(~same)[0]]
    return diff


if __name__ == "__main__":
    if "--check" in sys.argv:
        d = differences_from_committed()
        for row in d:
            print(row)
        print(f"[manifest] {len(d)} cells differ from the committed file")
        sys.exit(1 if d else 0)
    df = build()
    df.to_csv(DATA/"data_manifest.csv", index=False)
    print(f"[manifest] {len(df)} frames; {int((df.verification == 'verified').sum())} verified; "
          f"{int(df.use_in_study.sum())} used in the study")
