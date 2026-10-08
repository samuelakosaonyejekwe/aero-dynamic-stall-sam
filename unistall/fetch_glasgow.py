# -*- coding: utf-8 -*-
"""
unistall / fetch_glasgow.py
---------------------------
Fetch the University of Glasgow NACA 0012 oscillating-aerofoil runs, and the
static runs of the same wind-tunnel model, into the git-ignored cache; write
the inventories that pin every file by SHA-256; and read one run.

Author: Akosa Samuel Onyejekwe (independent)

TWO SOURCES

  (1) the sinusoidal runs, as redistributed in
      L. Pancini, BL-DSM-JFS-2021, directory "Glasgow Data",
      https://github.com/luizpancini/BL-DSM-JFS-2021 at the fixed commit
      SOURCE_COMMIT (the commit the NASA frames are fetched at). Per run:
        NNNNNNNN.dat          the original data file of the Glasgow database
        NNNNNNNN_coeffs.dat   force coefficients exported from it
        GUD_NNNNNNNN.mat      the same coefficients repacked by that
                              repository's generate_GUD.m
      That repository states no licence; nothing from it is committed here.

  (2) the original deposit, for what (1) does not carry:
      Green, R. B. & Giuni, M. (2017), "Dynamic stall database R and D
      1570-AM-01: Final Report", University of Glasgow,
      doi:10.5525/gla.researchdata.464 (CC BY 4.0), files
      dynamic_stall_database_final_report.pdf and GU_dynamic_stall.zip.
      From the archive only these members are read, by byte range, so the
      131 MB file is not downloaded whole: the 28 static runs of model 11
      (CorrectedFiles/1100*.dat), TransducerLocations.dat, the two run lists
      of model 11, and the two MATLAB sources that define how the coefficient
      files were computed (DSplot.m, aerofoil_taps.m). The archive's listing
      also gives the CRC-32 of its own copy of every sinusoidal data file of
      model 11; the inventory records whether each redistributed
      NNNNNNNN.dat has that CRC-32, that is, whether it is the deposited
      file unchanged.

WHAT A RUN FILE CONTAINS (final report, sections 2 and 4, Tables 7 and 12;
DSplot.m, function export_coeffs_Callback)

  model       NACA 0012, "model 11", chord 0.55 m, spanning the 5 ft
              dimension of the 7 ft x 5 ft closed-return low-speed tunnel,
              pitched about its QUARTER CHORD
  file name   abcdefgh: ab model, c type of experiment (0 = clean, no trip),
              d motion (0 static, 1 sinusoidal), efg test number, h attempt
  NNNNNNNN.dat
              first row, the run information block (RIB), 32 numbers; the
              ones used here, numbered from 1 as in Table 7:
                1 run number, 2-4 day, month, year, 5 temperature [deg C],
                6 barometric pressure [mm Hg], 7 motion type,
                8 mean incidence [deg] (static: starting incidence),
                9 amplitude [deg] (static: arc), 10 frequency [Hz],
                11 samples per channel, 13 number of cycles,
                18 sampling frequency [Hz], 19 dynamic pressure [Pa],
                20 Reynolds number, 21 Mach number, 22 reduced frequency
                k = omega c / (2 U), 23 wind speed [m/s],
                26 averaged (1) or not (0), 30 model number.
              then one row per sample: dynamic pressure [Pa], the pressure
              coefficient at each of 30 surface transducers (transducer 1 at
              the upper-surface trailing edge, 30 at the lower-surface
              trailing edge; chordwise positions in Table 12), incidence
              [deg]. For a sinusoidal run the rows are one cycle, the average
              of the sampled cycles.
  NNNNNNNN_coeffs.dat
              columns: phase omega t [rad] from the first sample, incidence
              [deg], C_n, C_t, C_m. They are sums of pressure times panel
              length over the 30 transducers:
                C_n  normal to the chord, positive towards the upper surface
                C_t  along the chord, positive TOWARDS THE LEADING EDGE
                     (leading-edge suction is positive)
                C_m  about the quarter chord, positive nose-up
              They are PRESSURE forces only: no skin friction is in them.

WHAT IS DERIVED HERE, and nothing else:
    C_l = C_n cos(alpha) + C_t sin(alpha)
    C_d = C_n sin(alpha) - C_t cos(alpha)        (pressure drag)
    measured mean incidence = (max + min)/2, amplitude = (max - min)/2 of the
    incidence column; each run is solved at these, not at the nominal values.

ASSUMPTIONS, each because the sources do not say otherwise:
  * no wall or blockage correction has been applied to any number in the
    files, and none is applied here (the final report names none; the
    directory of the deposit is called "CorrectedFiles", which the report
    explains as corrections of the test records against the paper reports);
  * the incidence column is the geometric pitch angle of the model;
  * the Mach number, reduced frequency and Reynolds number of the RIB are
    taken as stored.
data/glasgow_notes.md lists what this makes weaker.

HELD-OUT RUNS. data/glasgow_split.json assigns every sinusoidal run to
"calibration" or "held_out" by its run number. `load_run` refuses a held-out
run unless unistall/glasgow_blind.py has opened the single scoring; the
inventory reads only the conditions and the number of points of each run.

Usage:
  python3 -m unistall.fetch_glasgow     download what is missing, verify the
                                        checksums against the committed
                                        inventories (or create them), write
                                        data/glasgow_inventory.csv,
                                        data/glasgow_static_inventory.csv and
                                        data/glasgow_sources.csv
"""
import hashlib
import io
import json
import sys
import urllib.parse
import urllib.request
import zipfile
import zlib
from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io as sio

from unistall import paths
from unistall.fetch_frames import SOURCE_COMMIT
from unistall.paths import DATA

CACHE = paths.CACHE/"glasgow"
RUNS = CACHE/"runs"                 # files of source (1)
DEPOSIT = CACHE/"deposit"           # files of source (2)
INV = DATA/"glasgow_inventory.csv"
STATIC_INV = DATA/"glasgow_static_inventory.csv"
SOURCES = DATA/"glasgow_sources.csv"
SPLIT = DATA/"glasgow_split.json"

API = f"https://api.github.com/repos/luizpancini/BL-DSM-JFS-2021/git/trees/{SOURCE_COMMIT}?recursive=1"
RAW = f"https://raw.githubusercontent.com/luizpancini/BL-DSM-JFS-2021/{SOURCE_COMMIT}/"
FOLDER = "Glasgow Data/"
DEPOSIT_DOI = "10.5525/gla.researchdata.464"
DEPOSIT_ZIP = "https://researchdata.gla.ac.uk/464/1/GU_dynamic_stall.zip"
DEPOSIT_REPORT = "https://researchdata.gla.ac.uk/464/2/dynamic_stall_database_final_report.pdf"
REPORT_NAME = "dynamic_stall_database_final_report.pdf"
DEPOSIT_MEMBERS = ("TransducerLocations.dat", "ListOfRuns/model11_static.txt", "ListOfRuns/model11_osc_sine.txt",
                   "MatlabR2013a/DSplot.m", "MatlabR2013a/aerofoil_taps.m")
STATIC_PREFIX = "CorrectedFiles/1100"       # model 11, clean, static
SINE_PREFIX = "CorrectedFiles/1101"         # model 11, clean, sinusoidal: only the checksums are read
CRC_NAME = "model11_sinusoidal_crc32.csv"

MODEL_NUMBER = 11                   # NACA 0012 in the Glasgow database
CHORD_M = 0.55
PITCH_AXIS_X_C = 0.25
N_TRANSDUCERS = 30
RIB_LENGTH = N_TRANSDUCERS + 2
DOWNLOAD_THREADS = 8
# RIB entries, numbered from 1 as in Table 7 of the final report
RIB = dict(run=1, day=2, month=3, year=4, temperature_C=5, pressure_mmHg=6, motion=7, mean_deg=8, amp_deg=9,
           frequency_Hz=10, samples=11, cycles=13, sampling_Hz=18, q_Pa=19, Re=20, M=21, k=22, U_m_s=23,
           averaged=26, model=30)
MOTION_STATIC, MOTION_SINUSOIDAL = 0, 1


class HeldOutAccess(RuntimeError):
    """A held-out run was asked for outside the single blind scoring."""


_HELD_OUT_OPEN = False       # set only by unistall.glasgow_blind, inside its scoring


# --------------------------------------------------------------------------- #
#  the split
# --------------------------------------------------------------------------- #
def split_of(run: int) -> str:
    """The set of sinusoidal run number `run` by the rule frozen in
    data/glasgow_split.json: "calibration" or "held_out"."""
    with open(SPLIT) as fh:
        rule = json.load(fh)["rule"]
    return "calibration" if int(run) % rule["modulus"] == rule["calibration_remainder"] else "held_out"


# --------------------------------------------------------------------------- #
#  reading the files
# --------------------------------------------------------------------------- #
def read_dat(path: Path) -> tuple:
    """(rib, rows) of one database file: the 32 numbers of the run information
    block and the sample rows (dynamic pressure [Pa], 30 pressure
    coefficients, incidence [deg]). Files that carry three calibration rows
    after the RIB (DSplot.m allows for them) have those rows dropped."""
    with open(path) as fh:
        table = [[float(v) for v in line.split()] for line in fh if line.strip()]
    rib = np.array(table[0])
    if len(rib) != RIB_LENGTH:
        raise ValueError(f"{path.name}: run information block of {len(rib)} numbers, expected {RIB_LENGTH}")
    rows = np.array(table[1:])
    if len(rows) == int(rib[RIB["samples"] - 1]) + 3:
        rows = rows[3:]
    if rows.shape[1] != RIB_LENGTH:
        raise ValueError(f"{path.name}: sample rows of {rows.shape[1]} numbers, expected {RIB_LENGTH}")
    return rib, rows


def rib_value(rib: np.ndarray, name: str) -> float:
    """Entry `name` (a key of RIB) of a run information block."""
    return float(rib[RIB[name] - 1])


def read_coeffs(path: Path) -> np.ndarray:
    """The table of a _coeffs.dat file: columns phase [rad], incidence [deg],
    C_n, C_t, C_m."""
    return np.loadtxt(path, comments="%", ndmin=2)


def _conditions(run: int) -> dict:
    """Conditions and counts of one sinusoidal run, for the inventory. Only
    the run information block, the incidence column and the number of rows of
    each file are used; no coefficient value is returned."""
    rib, rows = read_dat(RUNS/f"{run}.dat")
    co = read_coeffs(RUNS/f"{run}_coeffs.dat")
    mat = sio.loadmat(str(RUNS/f"GUD_{run}.mat"))
    a = co[:, 1]
    same = all(len(mat[k].ravel()) == len(co) and np.allclose(mat[k].ravel(), co[:, j], rtol=0, atol=1e-12)
               for k, j in (("alpha", 1), ("cn", 2), ("cc", 3), ("cm", 4)))
    finite = [bool(np.isfinite(co[:, j]).all()) for j in (2, 3, 4)]
    return dict(
        run=run, test_number=(run // 10) % 1000, attempt=run % 10, set=split_of(run),
        model=int(rib_value(rib, "model")), experiment_type=(run // 100000) % 10,
        motion_type=int(rib_value(rib, "motion")),
        date=f"19{int(rib_value(rib, 'year')):02d}-{int(rib_value(rib, 'month')):02d}-{int(rib_value(rib, 'day')):02d}",
        M=rib_value(rib, "M"), k=rib_value(rib, "k"), Re=rib_value(rib, "Re"), U_m_s=rib_value(rib, "U_m_s"),
        frequency_Hz=rib_value(rib, "frequency_Hz"), mean_nominal_deg=rib_value(rib, "mean_deg"),
        amp_nominal_deg=rib_value(rib, "amp_deg"),
        alpha0_deg=round(float((a.max() + a.min())/2.0), 4), amp_deg=round(float((a.max() - a.min())/2.0), 4),
        peak_alpha_deg=round(float(a.max()), 4), min_alpha_deg=round(float(a.min()), 4),
        cycles_averaged=int(rib_value(rib, "cycles")), averaged=int(rib_value(rib, "averaged")),
        n_pressure_rows=len(rows), n_cn=int(len(co)) if finite[0] else 0, n_ct=int(len(co)) if finite[1] else 0,
        n_cm=int(len(co)) if finite[2] else 0,
        coefficients_stored="C_n C_t C_m (pressure forces; C_l and pressure C_d derived)",
        incidence_same_in_both_files=bool(np.allclose(rows[:, -1], a, rtol=0, atol=1e-9)),
        mat_equals_coeffs=bool(same),
        dat_crc32_equals_deposit=_same_as_deposit(run),
        sha256_dat=sha256_of(RUNS/f"{run}.dat"), sha256_coeffs=sha256_of(RUNS/f"{run}_coeffs.dat"),
        sha256_mat=sha256_of(RUNS/f"GUD_{run}.mat"))


@lru_cache(maxsize=1)
def _deposit_crc() -> dict:
    """CRC-32 of each sinusoidal data file of model 11 as the archive of the
    deposit lists it (written by fetch_deposit)."""
    t = pd.read_csv(DEPOSIT/CRC_NAME)
    return dict(zip(t.file, t.crc32, strict=True))


def _same_as_deposit(run: int) -> bool:
    """Whether the redistributed NNNNNNNN.dat has the CRC-32 the deposit's
    archive records for its own file of that name."""
    return bool(_deposit_crc().get(f"{run}.dat") == zlib.crc32((RUNS/f"{run}.dat").read_bytes()))


def load_run(run: int) -> dict:
    """One sinusoidal run in the form unistall.metrics.score_frame reads:
    frame, M, k, a0 and da [deg] (measured mean and amplitude), and the
    incidence [deg] and value of C_l, C_m, C_d point by point in the order
    sampled; also cn and ct as stored, U [m/s] and Re. Raises HeldOutAccess
    for a held-out run unless the blind scoring has been opened."""
    run = int(run)
    if split_of(run) != "calibration" and not _HELD_OUT_OPEN:
        raise HeldOutAccess(f"run {run} is held out: it is read only by unistall.glasgow_blind score")
    rib, _ = read_dat(RUNS/f"{run}.dat")
    co = read_coeffs(RUNS/f"{run}_coeffs.dat")
    a, cn, ct, cm = co[:, 1], co[:, 2], co[:, 3], co[:, 4]
    ar = np.radians(a)
    return dict(frame=str(run), M=rib_value(rib, "M"), k=rib_value(rib, "k"), Re=rib_value(rib, "Re"),
                U=rib_value(rib, "U_m_s"), a0=float((a.max() + a.min())/2.0), da=float((a.max() - a.min())/2.0),
                acl=a, cl=cn*np.cos(ar) + ct*np.sin(ar), acm=a, cm=cm,
                acd=a, cd=cn*np.sin(ar) - ct*np.cos(ar), cn=cn, ct=ct)


# --------------------------------------------------------------------------- #
#  inventories, built from the cached files only
# --------------------------------------------------------------------------- #
def sha256_of(path: Path) -> str:
    """SHA-256 of the file at `path`, as hexadecimal text."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def cached_runs() -> list:
    """Run numbers of the sinusoidal runs whose three files are in the cache."""
    runs = sorted(int(p.stem) for p in RUNS.glob("[0-9]" * 8 + ".dat"))
    return [r for r in runs if (RUNS/f"{r}_coeffs.dat").exists() and (RUNS/f"GUD_{r}.mat").exists()]


def inventory_from_cache() -> pd.DataFrame:
    """data/glasgow_inventory.csv as the cached files give it: one row per
    sinusoidal run, in order of run number."""
    return pd.DataFrame([_conditions(r) for r in cached_runs()])


def static_inventory_from_cache() -> pd.DataFrame:
    """data/glasgow_static_inventory.csv as the cached files give it: one row
    per static run of model 11."""
    rows = []
    for p in sorted((DEPOSIT/"static").glob("*.dat")):
        rib, d = read_dat(p)
        rows.append(dict(
            run=int(p.stem), model=int(rib_value(rib, "model")), experiment_type=(int(p.stem) // 100000) % 10,
            motion_type=int(rib_value(rib, "motion")),
            date=f"19{int(rib_value(rib, 'year')):02d}-{int(rib_value(rib, 'month')):02d}-{int(rib_value(rib, 'day')):02d}",
            M=rib_value(rib, "M"), Re=rib_value(rib, "Re"), U_m_s=rib_value(rib, "U_m_s"),
            start_nominal_deg=rib_value(rib, "mean_deg"), arc_nominal_deg=rib_value(rib, "amp_deg"),
            alpha_first_deg=round(float(d[0, -1]), 4), alpha_last_deg=round(float(d[-1, -1]), 4),
            alpha_min_deg=round(float(d[:, -1].min()), 4), alpha_max_deg=round(float(d[:, -1].max()), 4),
            n_pressure_rows=len(d), sha256_dat=sha256_of(p)))
    return pd.DataFrame(rows)


def sources_from_cache() -> pd.DataFrame:
    """data/glasgow_sources.csv: the documents and code read from the deposit
    and from the redistributing repository, each with its size and SHA-256."""
    rows = []
    for p in sorted(q for q in DEPOSIT.iterdir() if q.is_file()):
        where = DEPOSIT_REPORT if p.name == REPORT_NAME else f"{DEPOSIT_ZIP} (member {p.name})"
        if p.name == CRC_NAME:
            where = f"{DEPOSIT_ZIP} (CRC-32 of members {SINE_PREFIX}*.dat, from the archive listing)"
        rows.append(dict(file=p.name, source=where, doi=DEPOSIT_DOI, bytes=p.stat().st_size, sha256=sha256_of(p)))
    for p in sorted((CACHE/"repository").glob("*")):
        rows.append(dict(file=p.name, source=f"https://github.com/luizpancini/BL-DSM-JFS-2021 at {SOURCE_COMMIT}",
                         doi="", bytes=p.stat().st_size, sha256=sha256_of(p)))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
#  download
# --------------------------------------------------------------------------- #
class _HttpRange(io.RawIOBase):
    """A remote file read by HTTP byte ranges, so single members of a large
    zip archive can be taken without downloading the archive."""

    def __init__(self, url: str) -> None:
        self.url, self.pos = url, 0
        with urllib.request.urlopen(urllib.request.Request(url, method="HEAD"), timeout=60) as r:
            self.size = int(r.headers["Content-Length"])

    def seekable(self) -> bool:
        """Always true: any byte range can be asked for."""
        return True

    def readable(self) -> bool:
        """Always true."""
        return True

    def tell(self) -> int:
        """Current position [bytes]."""
        return self.pos

    def seek(self, offset: int, whence: int = 0) -> int:
        """Move to `offset` [bytes] from the start (0), here (1) or the end (2)."""
        self.pos = offset if whence == 0 else self.pos + offset if whence == 1 else self.size + offset
        return self.pos

    def readinto(self, buffer: bytearray) -> int:
        """Fill `buffer` from the current position; returns the bytes read."""
        n = min(len(buffer), self.size - self.pos)
        if n <= 0:
            return 0
        rq = urllib.request.Request(self.url, headers={"Range": f"bytes={self.pos}-{self.pos + n - 1}"})
        with urllib.request.urlopen(rq, timeout=120) as r:
            got = r.read()
        buffer[:len(got)] = got
        self.pos += len(got)
        return len(got)


def _get(url: str, dst: Path) -> None:
    if not dst.exists():
        with urllib.request.urlopen(url, timeout=120) as r:
            dst.write_bytes(r.read())


def fetch_runs() -> int:
    """Download the files of every sinusoidal run, and the three source files
    of the redistributing repository that document them, at SOURCE_COMMIT.
    Returns the number of runs."""
    RUNS.mkdir(parents=True, exist_ok=True)
    (CACHE/"repository").mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(API, timeout=60) as r:
        tree = [p["path"] for p in json.load(r)["tree"] if p["type"] == "blob"]
    names = sorted(Path(p).name for p in tree if p.startswith(FOLDER) and p.endswith((".dat", ".mat")))
    jobs = [(RAW + urllib.parse.quote(FOLDER + n), RUNS/n) for n in names]
    jobs += [(RAW + urllib.parse.quote(p), CACHE/"repository"/Path(p).name)
             for p in ("README.md", "Glasgow Data/generate_GUD.m", "src/functions/load_GUD.m")]
    with ThreadPoolExecutor(DOWNLOAD_THREADS) as pool:
        list(pool.map(lambda job: _get(*job), jobs))
    return len([n for n in names if n.endswith("_coeffs.dat")])


def fetch_deposit() -> int:
    """Download the final report, and from the archive of the deposit the
    members named in the module text. Returns the number of static runs."""
    (DEPOSIT/"static").mkdir(parents=True, exist_ok=True)
    _get(DEPOSIT_REPORT, DEPOSIT/REPORT_NAME)
    wanted = {m: DEPOSIT/Path(m).name for m in DEPOSIT_MEMBERS}
    if all(p.exists() for p in [*wanted.values(), DEPOSIT/CRC_NAME]) and any((DEPOSIT/"static").glob("*.dat")):
        return len(list((DEPOSIT/"static").glob("*.dat")))
    crc = []
    with zipfile.ZipFile(io.BufferedReader(_HttpRange(DEPOSIT_ZIP), 1 << 20)) as z:
        for info in z.infolist():
            name = info.filename.split("/", 1)[-1]
            if info.filename.startswith("__MACOSX") or info.is_dir():
                continue
            if name.startswith(SINE_PREFIX) and name.endswith(".dat"):
                crc.append(dict(file=Path(name).name, crc32=info.CRC, bytes=info.file_size))
            elif name in wanted:
                wanted[name].write_bytes(z.read(info))
            elif name.startswith(STATIC_PREFIX) and name.endswith(".dat"):
                (DEPOSIT/"static"/Path(name).name).write_bytes(z.read(info))
    pd.DataFrame(crc).sort_values("file").to_csv(DEPOSIT/CRC_NAME, index=False)
    return len(list((DEPOSIT/"static").glob("*.dat")))


def _changed(new: pd.DataFrame, committed: Path, key: str, columns: list) -> list:
    """Keys whose checksum differs from the committed inventory."""
    if not committed.exists():
        return []
    old = pd.read_csv(committed).set_index(key)
    return [f"{k}:{c}" for k, row in new.set_index(key).iterrows() for c in columns
            if k in old.index and old.loc[k, c] != row[c]]


def main() -> None:
    """Fetch, verify against the committed inventories, write them."""
    n_runs, n_static = fetch_runs(), fetch_deposit()
    inv, stat, src = inventory_from_cache(), static_inventory_from_cache(), sources_from_cache()
    changed = (_changed(inv, INV, "run", ["sha256_dat", "sha256_coeffs", "sha256_mat"])
               + _changed(stat, STATIC_INV, "run", ["sha256_dat"]) + _changed(src, SOURCES, "file", ["sha256"]))
    if changed:
        sys.exit(f"[fetch-glasgow] checksum changed for {changed}: the files are not the ones this study "
                 "recorded. Nothing written.")
    if len(inv) != n_runs:
        sys.exit(f"[fetch-glasgow] {n_runs} runs listed upstream, {len(inv)} complete in the cache. Nothing written.")
    inv.to_csv(INV, index=False)
    stat.to_csv(STATIC_INV, index=False)
    src.to_csv(SOURCES, index=False)
    print(f"[fetch-glasgow] {len(inv)} sinusoidal runs ({dict(inv['set'].value_counts())}), "
          f"{n_static} static runs, {len(src)} source documents in {CACHE}")


if __name__ == "__main__":
    main()
