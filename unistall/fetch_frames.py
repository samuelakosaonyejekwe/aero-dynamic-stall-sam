# -*- coding: utf-8 -*-
"""
unistall / fetch_frames.py
----------------------
Fetch the digitised NASA TM-84245 frame files at run time and write an
inventory of them. The files are third-party data with no stated licence
(see NOTICE), so they are downloaded into the git-ignored unistall/cache/frames/
and are NOT committed; what is committed is the inventory, which pins each
file by SHA-256 so a later fetch can be shown to be the same data.

Author: Akosa Samuel Onyejekwe (independent)

Source: L. Pancini, BL-DSM-JFS-2021, directory "NASA Data",
        https://github.com/luizpancini/BL-DSM-JFS-2021  (branch master)

Aerofoil identity follows that repository's src/functions/load_frame.m:
    7019-14220 NACA 0012 ; 24022-31310 AMES-01 ; >= 67000 NLR-7301.
Frames outside those three bands are listed with aerofoil "unassigned" and are
not used until identified from the NASA report itself.

Usage:
  python3 -m unistall.fetch_frames            download what is missing, verify
                                          checksums against the committed
                                          inventory (or create it), write
                                          data/frames_inventory.csv
"""
import sys, json, hashlib, urllib.request, urllib.parse
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.io as sio

HERE = Path(__file__).resolve().parent
from unistall import paths
from unistall.paths import DATA
CACHE = paths.CACHE/"frames"
INV = DATA/"frames_inventory.csv"
# The source repository is read at one fixed commit (its latest, of 17 May
# 2022), not at the tip of a branch, so a later change there cannot alter
# what is fetched; each file is also checked against the SHA-256 recorded in
# data/frames_inventory.csv.
SOURCE_COMMIT = "84e7945c5b16c6da4075c7fb467c01125dace36b"
API = f"https://api.github.com/repos/luizpancini/BL-DSM-JFS-2021/git/trees/{SOURCE_COMMIT}?recursive=1"
RAW = f"https://raw.githubusercontent.com/luizpancini/BL-DSM-JFS-2021/{SOURCE_COMMIT}/"


def airfoil_of(n):
    if 7019 <= n <= 14220:
        return "NACA 0012"
    if 24022 <= n <= 31310:
        return "AMES-01"
    if n >= 67000:
        return "NLR-7301"
    return "unassigned"


def split_of(n, airfoil, fixed):
    """The rule frozen in data/split.json."""
    name = f"frame_{n}"
    if name in fixed:
        return fixed[name]
    if airfoil == "NACA 0012":
        return "calibration" if n % 3 == 0 else "held_out"
    return "other_airfoil" if airfoil != "unassigned" else "unassigned"


def main():
    CACHE.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(API, timeout=60) as r:
        tree = json.load(r)["tree"]
    paths = sorted(p["path"] for p in tree
                   if p["path"].startswith("NASA Data/frame_") and p["path"].endswith(".mat"))
    split = json.load(open(DATA/"split.json"))
    fixed = {e["frame"]: s for s, lst in split["frames"].items() for e in lst}
    known = pd.read_csv(INV).set_index("frame")["sha256"].to_dict() if INV.exists() else {}
    rows, changed = [], []
    for p in paths:
        name = Path(p).stem
        dst = CACHE/f"{name}.mat"
        if not dst.exists():
            with urllib.request.urlopen(RAW + urllib.parse.quote(p), timeout=60) as r:
                dst.write_bytes(r.read())
        sha = hashlib.sha256(dst.read_bytes()).hexdigest()
        if name in known and known[name] != sha:
            changed.append(name)
        d = sio.loadmat(str(dst))
        def g(k, d=d):
            return float(d[k].ravel()[0])
        n = int(name.split("_")[1])
        af = airfoil_of(n)
        a0, da = np.degrees(g("alpha_0")), np.degrees(g("delta_alpha"))
        def npt(k, d=d):
            return int(d[k].size) if k in d else 0
        rows.append(dict(frame=name, frame_number=n, airfoil=af, set=split_of(n, af, fixed),
                         M=round(g("M"), 3), k=round(g("k"), 4),
                         alpha0_deg=round(a0, 2), amp_deg=round(da, 2),
                         peak_alpha_deg=round(a0 + da, 2),
                         n_cl=npt("cl_exp"), n_cm=npt("cm_exp"), n_cd=npt("cd_exp"),
                         cl_max=round(float(d["cl_exp"].max()), 3) if npt("cl_exp") else np.nan,
                         cm_min=round(float(d["cm_exp"].min()), 3) if npt("cm_exp") else np.nan,
                         sha256=sha))
    if changed:
        sys.exit(f"[fetch] checksum changed for {changed}: the upstream data are not the "
                 "data this study recorded. Nothing written.")
    df = pd.DataFrame(rows).sort_values("frame_number")
    df.to_csv(INV, index=False)
    print(f"[fetch] {len(df)} frames in {CACHE}; "
          f"{dict(df.airfoil.value_counts())}")


if __name__ == "__main__":
    main()
