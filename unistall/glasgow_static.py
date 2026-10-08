# -*- coding: utf-8 -*-
"""
unistall / glasgow_static.py
----------------------------
The static inputs of the load model for the University of Glasgow NACA 0012
runs (Mach 0.08 to 0.16, chord Reynolds number 1.0 to 2.0 million):
lift-curve slope, zero-lift incidence, static-stall incidence, the normal
force there (C_N1), the separation point f(alpha), the static moment and the
static pressure drag at zero lift. No held-out run is read.

Author: Akosa Samuel Onyejekwe (independent)

TWO SOURCES, built by the same rule and kept apart.

  "static"        (the one the model uses) the static runs of the same model
                  in the same tunnel, from the original deposit
                  (doi:10.5525/gla.researchdata.464, files 1100efgh.dat; see
                  unistall/fetch_glasgow.py). Each is one sweep up from about
                  -5 deg to about 25 deg and back, 128 points, of surface
                  pressures; the forces are formed here by the panel sum the
                  deposit's own code uses (DSplot.m, aerofoil_taps.m), which
                  is how the coefficient files of the sinusoidal runs were
                  made. The sum is checked against those files on the
                  calibration-type runs (results/glasgow_integration_check.csv).
                  Runs are grouped into three Mach stations by MACH_BANDS.
                  A station is built from the static runs made in the days of
                  the sinusoidal tests (from PERIOD_START on); the runs of the
                  earlier weeks are described too and the difference written
                  out.
  "quasi_steady"  the slowest calibration-type sinusoidal runs (k < QUASI_K,
                  all at Mach 0.116 to 0.119) read as sweeps, as
                  unistall/static_model.py reads the NASA quasi-steady frames:
                  the one of widest incidence range, extended above its peak
                  by the one of highest peak incidence. At k = 0.01 these are
                  not static: stall is delayed. They are here to measure how
                  far a static description taken from the sinusoidal runs
                  alone would differ, not to be used.

RULE, the same for every sweep (it is the rule of static_stations.py and
static_model.py):
  slope, zero-lift incidence   least-squares line through C_l of the up
                               stroke between FIT_RANGE_DEG
  static-stall incidence       where C_l of the up stroke is largest
  C_N1                         the measured C_n at that point
  two post-stall curves        "less" separated: the up stroke;
                               "more" separated (default): the up stroke to
                               static stall, the down stroke beyond it
  f                            Kirchhoff inverse with the sweep's own slope and
                               zero-lift incidence (static_model._kirchhoff_f)
  static moment                C_m of the same strokes
  C_D0                         pressure C_d of the up stroke at the zero-lift
                               incidence
A station is the mean over its sweeps: slope, zero-lift incidence, stall
incidence, C_N1 and C_D0 directly; f and C_m along incidence divided by each
sweep's stall incidence, so the stall break is not smeared by its scatter.

BETWEEN STATIONS everything is interpolated in Mach number exactly as
static_model.StaticModel does between its two. The stations differ in
Reynolds number as much as in Mach number (the tunnel speed was changed);
Mach number is only the label the interpolation uses.

WHAT IS NOT IN THESE NUMBERS. The forces are pressure forces: the drag has no
skin friction and the chord force none either. No wall correction is applied.
The measured zero-lift incidence of this symmetric section is +0.2 to +0.6 deg
and is carried as an offset of the incidence scale, as for the NASA data.

Outputs (results/):
  glasgow_integration_check.csv   the panel sum against the stored coefficients
  glasgow_static_runs.csv         every sweep, one row
  glasgow_static_stations.csv     the stations of both sources and periods
  glasgow_static_table.csv        f and static C_m against incidence
  glasgow_static_difference.csv   quasi-steady source, and earlier period,
                                  against the static source
Usage:  python3 -m unistall.glasgow_static
"""
from functools import lru_cache

import numpy as np
import pandas as pd

from unistall import fetch_glasgow as fg
from unistall import naca4
from unistall.paths import RESULTS
from unistall.static_model import (ALPHA_MAX, CURVES, TABLE_STEP_DEG, StaticModel, Station, _kirchhoff_f,
                                   _monotone_xy)
from unistall.strokes import stroke_split

THICKNESS = 0.12
FIT_RANGE_DEG = (-5.0, 8.0)                     # as unistall/static_stations.py
GRID = np.arange(0.0, ALPHA_MAX + 1e-9, 0.5)
MACH_BANDS = ((0.0, 0.10), (0.10, 0.14), (0.14, 0.20))
PERIOD_START = "1991-03-01"                     # date of the first sinusoidal run in the inventory
MIN_START_DEG = -10.0                           # static sweeps starting below this cover negative incidence only
QUASI_K = 0.0125                                # sinusoidal runs slower than this are read as sweeps
MACH_MARGIN = 0.016                             # as static_model.StaticModel.in_range
SOURCES = ("static", "quasi_steady")


# --------------------------------------------------------------------------- #
#  pressure to force: the panel sum of the deposit's DSplot.m / aerofoil_taps.m
# --------------------------------------------------------------------------- #
@lru_cache(maxsize=1)
def panels() -> dict:
    """Transducer positions and panel geometry of model 11, fractions of
    chord: x, y of the 30 transducers, the panel length lp each one stands
    for (between the intersections of the surface tangents at neighbouring
    transducers, closed at the trailing edge x = 1), and the direction of the
    outward normal [rad from the chord line]."""
    x = np.loadtxt(fg.DEPOSIT/"TransducerLocations.dat", comments="%", usecols=range(fg.N_TRANSDUCERS))[
        fg.MODEL_NUMBER - 1]
    n_upper = int(np.where(np.diff(x) < 0)[0][-1]) + 2            # transducers run upper TE -> LE -> lower TE
    side = np.r_[np.ones(n_upper), -np.ones(len(x) - n_upper)]
    y = side*naca4.thickness(x, THICKNESS)
    slope = side*5*THICKNESS*(naca4.A0/(2*np.sqrt(x)) + naca4.A1 + 2*naca4.A2*x + 3*naca4.A3*x**2
                              + 4*naca4.A4_OPEN*x**3)
    y_te = float(naca4.thickness(1.0, THICKNESS))
    xl, yl = np.empty(len(x) + 1), np.empty(len(x) + 1)
    xl[0], yl[0], xl[-1], yl[-1] = 1.0, y_te, 1.0, -y_te
    xl[1:-1] = ((y[1:] - y[:-1]) + (x[:-1]*slope[:-1] - x[1:]*slope[1:]))/(slope[:-1] - slope[1:])
    yl[1:-1] = y[:-1] + slope[:-1]*(xl[1:-1] - x[:-1])
    normal = np.arctan(slope) + np.pi/2 + np.where(side < 0, np.pi, 0.0)
    return dict(x=x, y=y, lp=np.hypot(np.diff(xl), np.diff(yl)), normal=normal)


def forces(rows: np.ndarray) -> dict:
    """Incidence [deg], C_n, C_t (positive towards the leading edge), C_m
    (quarter chord, nose-up) and pressure C_l, C_d from the sample rows of a
    database file (fetch_glasgow.read_dat)."""
    p = panels()
    cp = rows[:, 1:-1]
    fy = -cp*p["lp"]*np.sin(p["normal"])
    fx = -cp*p["lp"]*np.cos(p["normal"])
    a, cn, ct = rows[:, -1], fy.sum(axis=1), -fx.sum(axis=1)
    cm = (-fy*(p["x"] - fg.PITCH_AXIS_X_C) + fx*p["y"]).sum(axis=1)
    ar = np.radians(a)
    return dict(a=a, cn=cn, ct=ct, cm=cm, cl=cn*np.cos(ar) + ct*np.sin(ar), cd=cn*np.sin(ar) - ct*np.cos(ar))


def integration_check() -> pd.DataFrame:
    """The panel sum against the coefficients stored with each
    calibration-type sinusoidal run: largest and RMS difference per run."""
    inv = pd.read_csv(fg.INV)
    out = []
    for run in inv[inv["set"] == "calibration"].run:
        _, rows = fg.read_dat(fg.RUNS/f"{run}.dat")
        f, co = forces(rows), fg.read_coeffs(fg.RUNS/f"{run}_coeffs.dat")
        d = {q: f[q] - co[:, j] for q, j in (("cn", 2), ("ct", 3), ("cm", 4))}
        out.append(dict(run=run, **{f"max_abs_d{q}": float(np.abs(v).max()) for q, v in d.items()},
                        **{f"rms_d{q}": float(np.sqrt(np.mean(v**2))) for q, v in d.items()}))
    return pd.DataFrame(out)


# --------------------------------------------------------------------------- #
#  one sweep
# --------------------------------------------------------------------------- #
def _strokes(f: dict) -> dict:
    """Up and down stroke of one record, each as (alpha, C_n, C_m, C_l, C_d)
    sorted on incidence with repeated incidences merged."""
    s = stroke_split(f["a"])
    out = {}
    for st in ("up", "down"):
        m = s == st
        cols = [_monotone_xy(f["a"][m], f[q][m]) for q in ("cn", "cm", "cl", "cd")]
        out[st] = (cols[0][0], *[c[1] for c in cols])
    return out


def describe(up: tuple, down: tuple, label: dict) -> dict:
    """The static description of one sweep from its up and down stroke (each
    as _strokes gives it), with `label` copied into the result."""
    a, cn, cm, cl, cd = up
    lin = (a >= FIT_RANGE_DEG[0]) & (a <= FIT_RANGE_DEG[1])
    slope, intercept = np.polyfit(a[lin], cl[lin], 1)
    a0 = float(-intercept/slope)
    i = int(np.argmax(cl))
    cn_less, cm_less = np.interp(GRID, a, cn), np.interp(GRID, a, cm)
    beyond = GRID > a[i]
    cn_more = np.minimum(np.where(beyond, np.interp(GRID, down[0], down[1]), cn_less), cn_less)
    cm_more = np.where(beyond, np.interp(GRID, down[0], down[2]), cm_less)
    return dict(**label, slope_per_deg=float(slope), alpha0_deg=a0, alpha_stall_deg=float(a[i]),
                CL_max=float(cl[i]), CN1=float(cn[i]), CD0=float(np.interp(a0, a, cd)),
                CM_at_zero_lift=float(np.interp(a0, a, cm)), alpha_max_deg=float(a.max()),
                points_in_fit=int(lin.sum()),
                f={"more": _kirchhoff_f(GRID, cn_more, slope, a0), "less": _kirchhoff_f(GRID, cn_less, slope, a0)},
                cm={"more": cm_more, "less": cm_less})


@lru_cache(maxsize=1)
def static_sweeps() -> list:
    """Every static run of the deposit that covers positive incidence,
    described."""
    out = []
    for r in pd.read_csv(fg.STATIC_INV).itertuples():
        if r.start_nominal_deg < MIN_START_DEG:
            continue
        _, rows = fg.read_dat(fg.DEPOSIT/"static"/f"{r.run}.dat")
        st = _strokes(forces(rows))
        out.append(describe(st["up"], st["down"], dict(source="static", run=str(r.run), M=float(r.M), Re=float(r.Re),
                                                        date=r.date, in_period=bool(r.date >= PERIOD_START))))
    return out


def quasi_steady_runs() -> pd.DataFrame:
    """The calibration-type sinusoidal runs slower than QUASI_K, from the
    inventory (conditions only)."""
    inv = pd.read_csv(fg.INV)
    return inv[(inv["set"] == "calibration") & (inv.k < QUASI_K)].reset_index(drop=True)


@lru_cache(maxsize=1)
def quasi_steady_sweep() -> dict:
    """The quasi-steady description: the slow calibration-type run of widest
    incidence range, extended above its peak by the slow calibration-type run
    of highest peak incidence."""
    q = quasi_steady_runs()
    main = q.sort_values(["amp_deg", "k"], ascending=[False, True]).iloc[0]
    high = q.sort_values("peak_alpha_deg", ascending=False).iloc[0]
    s1, s2 = (_strokes({**(fr := fg.load_run(int(r.run))), "a": fr["acl"]}) for r in (main, high))
    hi = s2["up"][0] > s1["up"][0].max()
    up = tuple(np.r_[u1, u2[hi]] for u1, u2 in zip(s1["up"], s2["up"], strict=True))
    a_down = np.r_[s1["down"][0], s2["down"][0]]
    cols = [_monotone_xy(a_down, np.r_[d1, d2]) for d1, d2 in zip(s1["down"][1:], s2["down"][1:], strict=True)]
    down = (cols[0][0], *[c[1] for c in cols])
    return describe(up, down, dict(source="quasi_steady", run=f"{int(main.run)}+{int(high.run)}", M=float(main.M),
                                   Re=float(main.Re), date=main.date, in_period=True))


# --------------------------------------------------------------------------- #
#  stations
# --------------------------------------------------------------------------- #
SCALARS = ("slope_per_deg", "alpha0_deg", "alpha_stall_deg", "CL_max", "CN1", "CD0", "CM_at_zero_lift")


def combine(sweeps: list) -> dict:
    """The mean of several sweep descriptions (see the module text), with the
    standard deviation of each scalar over the sweeps."""
    mean = {k: float(np.mean([s[k] for s in sweeps])) for k in ("M", "Re", *SCALARS)}
    sd = {f"sd_{k}": float(np.std([s[k] for s in sweeps], ddof=1)) if len(sweeps) > 1 else float("nan")
          for k in SCALARS}
    def along(key: str, c: str) -> np.ndarray:
        """Mean of a tabulated curve along incidence scaled to the mean stall incidence."""
        return np.mean([np.interp(GRID*s["alpha_stall_deg"]/mean["alpha_stall_deg"], GRID, s[key][c])
                        for s in sweeps], axis=0)
    return dict(**mean, **sd, n_sweeps=len(sweeps), runs=" ".join(s["run"] for s in sweeps),
                f={c: along("f", c) for c in CURVES}, cm={c: along("cm", c) for c in CURVES})


def _band(M: float) -> int:
    return next(i for i, (lo, hi) in enumerate(MACH_BANDS) if lo <= M < hi)


@lru_cache(maxsize=None)
def station_table(source: str = "static", period: str = "in_period") -> tuple:
    """The stations of one source, in order of Mach number, each as combine
    returns it. `period` applies to the static source: "in_period" (static
    runs from PERIOD_START on, the default), "earlier" or "all"."""
    if source == "quasi_steady":
        return (combine([quasi_steady_sweep()]),)
    keep = {"in_period": lambda s: s["in_period"], "earlier": lambda s: not s["in_period"], "all": lambda s: True}[period]
    sweeps = [s for s in static_sweeps() if keep(s)]
    groups = [[s for s in sweeps if _band(s["M"]) == i] for i in range(len(MACH_BANDS))]
    return tuple(combine(g) for g in groups if g)


def _as_station(c: dict, source: str) -> Station:
    return Station(c["M"], c["slope_per_deg"], c["alpha0_deg"], c["alpha_stall_deg"], c["CN1"], GRID, c["f"], c["cm"],
                   f"Glasgow {source}: runs {c['runs']}")


class GlasgowStatic(StaticModel):
    """Static inputs at Mach number M [-] for the Glasgow runs, with the
    interface unistall.dsmodel reads (slope_per_deg [1/deg], CN_alpha
    [1/rad], alpha0_deg, alpha_stall_deg [deg], CN1, f_at, cm_at, in_range)
    and CD0, the static pressure drag at zero lift. `curve` is the post-stall
    curve ("more" or "less" separated); `source` is "static" or
    "quasi_steady"."""

    def __init__(self, M: float, curve: str = "more", source: str = "static") -> None:
        table = station_table(source)
        sts = [_as_station(c, source) for c in table]
        j = int(np.clip(np.searchsorted([s.M for s in sts], M) - 1, 0, max(len(sts) - 2, 0)))
        a, b = sts[j], sts[min(j + 1, len(sts) - 1)]
        self.M, self.curve, self.source = float(M), curve, source
        self.in_range = bool(sts[0].M - MACH_MARGIN <= M <= sts[-1].M + MACH_MARGIN)
        w = float(np.clip((M - a.M)/(b.M - a.M), 0.0, 1.0)) if b.M > a.M else 0.0
        self._a, self._b, self.w = a, b, w

        def lerp(x: float, y: float) -> float:
            """Linear interpolation in Mach number between the two stations."""
            return (1.0 - w)*x + w*y
        self.slope_per_deg = lerp(a.slope_per_deg, b.slope_per_deg)
        self.CN_alpha = self.slope_per_deg*180.0/np.pi
        self.alpha0_deg = lerp(a.alpha0_deg, b.alpha0_deg)
        self.alpha_stall_deg = lerp(a.alpha_stall_deg, b.alpha_stall_deg)
        self.CN1 = lerp(a.CN1, b.CN1)
        self.CD0 = lerp(table[j]["CD0"], table[min(j + 1, len(sts) - 1)]["CD0"])
        self._grid = np.arange(0.0, ALPHA_MAX + TABLE_STEP_DEG/2, TABLE_STEP_DEG)
        self._ftab = {c: self._between("f", self._grid, c) for c in CURVES}
        self._cmtab = {c: self._between("cm_static", self._grid, c) for c in CURVES}


# --------------------------------------------------------------------------- #
#  tables
# --------------------------------------------------------------------------- #
def _rows_of_stations() -> pd.DataFrame:
    rows = []
    for source, period in (("static", "in_period"), ("static", "earlier"), ("static", "all"),
                           ("quasi_steady", "in_period")):
        for c in station_table(source, period):
            rows.append(dict(source=source, period=period if source == "static" else "sinusoidal tests",
                             used_by_model=bool(source == "static" and period == "in_period"),
                             **{k: round(v, 5) for k, v in c.items() if isinstance(v, float)},
                             n_sweeps=c["n_sweeps"], runs=c["runs"]))
    return pd.DataFrame(rows)


def _difference(a: dict, b: dict, what: str) -> dict:
    """Station b less station a, scalar by scalar and along the tables."""
    out = dict(comparison=what, M_reference=round(a["M"], 4), M_other=round(b["M"], 4),
               **{f"d_{k}": round(b[k] - a[k], 5) for k in SCALARS},
               pct_slope=round(100.0*(b["slope_per_deg"] - a["slope_per_deg"])/a["slope_per_deg"], 2),
               pct_CN1=round(100.0*(b["CN1"] - a["CN1"])/a["CN1"], 2))
    for c in CURVES:
        df, dm = b["f"][c] - a["f"][c], b["cm"][c] - a["cm"][c]
        out.update({f"max_abs_df_{c}": round(float(np.abs(df).max()), 4), f"rms_df_{c}": round(float(np.sqrt(np.mean(df**2))), 4),
                    f"alpha_of_max_df_{c}_deg": float(GRID[int(np.abs(df).argmax())]),
                    f"max_abs_dCm_{c}": round(float(np.abs(dm).max()), 4)})
    return out


def differences() -> pd.DataFrame:
    """How far the quasi-steady source, and the static runs of the earlier
    weeks, differ from the static source the model uses."""
    used = station_table("static", "in_period")
    mid = min(used, key=lambda c: abs(c["M"] - quasi_steady_sweep()["M"]))
    rows = [_difference(mid, station_table("quasi_steady")[0], "quasi-steady sweeps less static runs")]
    for early in station_table("static", "earlier"):
        ref = min(used, key=lambda c, e=early: abs(c["M"] - e["M"]))
        rows.append(_difference(ref, early, "static runs before the period less static runs in it"))
    return pd.DataFrame(rows)


def main() -> None:
    """Write the result files named in the module text."""
    chk = integration_check()
    chk.round(6).to_csv(RESULTS/"glasgow_integration_check.csv", index=False)
    print("[glasgow-static] panel sum against stored coefficients, calibration-type runs: largest difference "
          + ", ".join(f"{q} {chk[f'max_abs_d{q}'].max():.4f}" for q in ("cn", "ct", "cm")))
    sweeps = [*static_sweeps(), quasi_steady_sweep()]
    pd.DataFrame([{k: (round(v, 5) if isinstance(v, float) else v) for k, v in s.items() if k not in ("f", "cm")}
                  for s in sweeps]).to_csv(RESULTS/"glasgow_static_runs.csv", index=False)
    st = _rows_of_stations()
    st.to_csv(RESULTS/"glasgow_static_stations.csv", index=False)
    print(st[["source", "period", "M", "slope_per_deg", "alpha0_deg", "alpha_stall_deg", "CL_max", "CN1", "CD0",
              "n_sweeps"]].to_string(index=False))
    tab = []
    for source in SOURCES:
        for c in station_table(source):
            for i, al in enumerate(GRID):
                tab.append(dict(source=source, M=round(c["M"], 4), alpha_deg=al,
                                **{f"f_{k}": round(float(c["f"][k][i]), 4) for k in CURVES},
                                **{f"Cm_{k}": round(float(c["cm"][k][i]), 4) for k in CURVES}))
    pd.DataFrame(tab).to_csv(RESULTS/"glasgow_static_table.csv", index=False)
    d = differences()
    d.to_csv(RESULTS/"glasgow_static_difference.csv", index=False)
    print(d.T.to_string())


if __name__ == "__main__":
    main()
