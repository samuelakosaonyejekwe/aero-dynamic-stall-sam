# -*- coding: utf-8 -*-
"""
unistall / source_log.py
--------------------
Write, into data/data_sources.csv, which tables, figures and pages of each
source the data files of this repository actually use. The entries are
collected from the source, location and basis columns of the files in data/,
so the log cannot say more, or less, than the data files do.

Author: Akosa Samuel Onyejekwe (independent)

Usage:  python3 -m unistall.source_log            rewrite the column parts_used
        python3 -m unistall.source_log --check    compare with the committed file
"""
import re
import sys
import pandas as pd

from unistall.paths import DATA

LOG = DATA/"data_sources.csv"
VOLUME_KEYS = {"19820024438": ("Vol. 1", "NTRS 19820024438"), "19830003778": ("Vol. 2", "NTRS 19830003778"),
               "19830009234": ("Vol. 3", "NTRS 19830009234"), "19880002254": ("TM-100019", "NTRS 19880002254")}
TEXT_COLUMNS = ("source", "report_location", "Re_basis", "flag")
VOLUME_1_FILES = ("data_manifest.csv", "report_tables_naca0012.csv", "report_tables_other_airfoils.csv",
                  "report_table9_zero_drift.csv")          # their locations are tables of Volume 1


def _strings() -> list:
    """(file name, text) for every citation string in the data files."""
    out = []
    for f in sorted(DATA.glob("*.csv")):
        if f.name == LOG.name:
            continue
        d = pd.read_csv(f)
        for c in TEXT_COLUMNS:
            if c in d.columns:
                out += [(f.name, str(v)) for v in d[c].dropna().unique()]
        if {"table", "report_page"} <= set(d.columns):
            out += [(f.name, f"Table {int(r.table)} p.{int(r.report_page)}")
                    for r in d[["table", "report_page"]].drop_duplicates().itertuples()]
    return out


def parts_used(ntrs_id: str) -> str:
    """The tables, figures and pages of one NASA document that the data files cite."""
    keys = VOLUME_KEYS.get(str(ntrs_id))
    if keys is None:
        return ""
    tables, figures, pages = set(), set(), set()
    for name, text in _strings():
        mine = any(k in text for k in keys) or (keys[0] == "Vol. 1" and name in VOLUME_1_FILES
                                                and not any(k in text for ks in VOLUME_KEYS.values() for k in ks))
        if not mine:
            continue
        tables |= {int(n) for n in re.findall(r"Table (\d+)", text)}
        figures |= set(re.findall(r"Fig\. (\d+(?:\([a-z]\))?)", text))
        pages |= {int(n) for n in re.findall(r"p\.\s?(\d+)", text)}
    parts = []
    if tables:
        parts.append("Tables " + ", ".join(str(t) for t in sorted(tables)))
    if figures:
        parts.append("Figs " + ", ".join(sorted(figures)))
    if pages:
        parts.append("printed pages " + ", ".join(str(p) for p in sorted(pages)))
    return "; ".join(parts) if parts else "no table, figure or page of it is cited by a data file"


def build() -> pd.DataFrame:
    """The source log with its column parts_used filled from the data files."""
    log = pd.read_csv(LOG, dtype=str)
    log["parts_used"] = [parts_used(i) or "see status" for i in log.ntrs_id]
    return log


if __name__ == "__main__":
    new = build()
    if "--check" in sys.argv:
        old = pd.read_csv(LOG, dtype=str)
        same = "parts_used" in old.columns and list(old.parts_used.fillna("")) == list(new.parts_used.fillna(""))
        print("[sources] parts_used", "agrees with" if same else "differs from", "the data files")
        sys.exit(0 if same else 1)
    new.to_csv(LOG, index=False)
    print(new[["source", "parts_used"]].to_string(index=False))
