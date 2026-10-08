# -*- coding: utf-8 -*-
"""
unistall / glasgow_blind.py
---------------------------
The blind test of the load model on the University of Glasgow NACA 0012 runs:
the forecast on the calibration-type runs, the go rule, the seal, and the
single scoring of the held-out runs.

Author: Akosa Samuel Onyejekwe (independent)

THE ORDER, each step a command of this module:

  targets    write data/glasgow_targets.json: the targets (those of
             data/targets.json, "held_out_validation", copied), the
             procedure and the go rule. Refuses to replace the file.
  forecast   on the calibration-type runs only: cross-validation by the
             procedure of unistall/calibrate.py (same cost, same bounds, same
             starts, same folds, same selection rule), with the static inputs
             of unistall/glasgow_static.py. Writes
             results/glasgow_forecast.csv (every left-out run scored once,
             for each candidate set of fitted constants),
             results/glasgow_forecast_summary.csv and
             results/glasgow_verdict.json (the go rule applied).
             `--static=quasi_steady` repeats it with the other static source
             into files of that suffix; those take no part in the rule.
  constants  fit the chosen set of constants on all calibration-type runs,
             at the search march as every fold of the forecast does, and
             write results/glasgow_constants.json.
  seal       compute the SHA-256 of every file that fixes the test and write
             results/glasgow_seal.json with the date, the forecast numbers
             and the verdict. Refuses if the verdict is no-go, and refuses to
             replace a seal.
  score      refuses unless the seal exists, every checksum still matches,
             the verdict is go and no blind score exists. Then solves every
             held-out run once at its own measured conditions, scores it with
             unistall/metrics.py, writes results/glasgow_blind_scores.csv and
             results/glasgow_blind_summary.csv, and records the date in the
             seal. If the verdict is no-go this command is never run: the
             held-out runs stay unscored and the forecast is what is reported.

THE GUARD. Measured coefficients of a run are read through `load` only. It
raises HeldOutAccess for a held-out run unless the purpose is "score" and the
scoring has been opened by `score` after its checks. Under it,
fetch_glasgow.load_run refuses a held-out run on its own.

THE GROUPS (from the conditions in data/glasgow_inventory.csv, never from a
result):
  primary                    Mach number inside the range of the static
                             stations and measured peak incidence at or
                             below dsmodel.PEAK_ALPHA_MAX_COMPARED (25 deg).
                             The fit, the forecast, the go rule and the
                             targets use this group only.
  beyond_static_lift_range   peak incidence above 25 deg: scored and
                             reported, no target applies.

Usage:  python3 -m unistall.glasgow_blind targets|forecast|constants|seal|score
"""
import datetime
import hashlib
import json
import sys
from contextlib import contextmanager
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from unistall import calibrate as cal
from unistall import dsmodel as dm
from unistall import fetch_glasgow as fg
from unistall import glasgow_static as gs
from unistall import metrics as mt
from unistall.paths import DATA, PACKAGE, RESULTS

AUTHOR = "Akosa Samuel Onyejekwe (independent)"
TARGETS = DATA/"glasgow_targets.json"
NASA_TARGETS = DATA/"targets.json"
FORECAST = RESULTS/"glasgow_forecast.csv"
FORECAST_SUMMARY = RESULTS/"glasgow_forecast_summary.csv"
VERDICT = RESULTS/"glasgow_verdict.json"
CONSTANTS = RESULTS/"glasgow_constants.json"
SEAL = RESULTS/"glasgow_seal.json"
SCORES = RESULTS/"glasgow_blind_scores.csv"
SUMMARY = RESULTS/"glasgow_blind_summary.csv"
ATTACHED_FILE = RESULTS/"attached_moment_factor.json"
SOURCE_FILES = ("dsmodel.py", "attached_flow.py", "static_model.py", "metrics.py", "strokes.py",
                "glasgow_static.py", "fetch_glasgow.py", "glasgow_blind.py")
STATIC_FILES = ("glasgow_static_stations.csv", "glasgow_static_table.csv")
CANDIDATES = ("none", "sensitive only", "all")          # in order of size, as calibrate.select
RULE_MEASURES = ("nRMS_CL", "nRMS_CM", "nRMS_CD")
PROCESSES = 2                                           # worker processes of the forecast and of the fit (two cores)
TARGET_MEASURES = (("mean_nRMS_CL", "mean_nRMS_CL_max", "at most"), ("mean_nRMS_CM", "mean_nRMS_CM_max", "at most"),
                   ("mean_nRMS_CD", "mean_nRMS_CD_max", "at most"),
                   ("mean_abs_dalpha_CLmax_deg", "mean_abs_dalpha_CLmax_deg_max", "at most"),
                   ("mean_abs_dalpha_Mstall_deg", "mean_abs_dalpha_Mstall_deg_max", "at most"),
                   ("mean_Xi_sign_agree", "damping_sign_agreement_min", "at least"))


class NoGo(RuntimeError):
    """The forecast does not meet the go rule: the held-out runs are not scored."""


class SealError(RuntimeError):
    """The seal is missing, broken, or the scoring has already been done."""


# --------------------------------------------------------------------------- #
#  the guard
# --------------------------------------------------------------------------- #
@contextmanager
def _scoring_open() -> object:
    """Opens the held-out runs for the duration of the single scoring."""
    fg._HELD_OUT_OPEN = True
    try:
        yield
    finally:
        fg._HELD_OUT_OPEN = False


def load(run: int, purpose: str) -> dict:
    """The measured loop of run number `run` (fetch_glasgow.load_run).
    `purpose` is "forecast", "constants", "seal" or "score"; a held-out run is
    returned only for "score", and only while `score` has the scoring open."""
    if fg.split_of(run) != "calibration" and (purpose != "score" or not fg._HELD_OUT_OPEN):
        raise fg.HeldOutAccess(f"run {run} is held out and was asked for with purpose '{purpose}'")
    return fg.load_run(run)


# --------------------------------------------------------------------------- #
#  runs and groups
# --------------------------------------------------------------------------- #
def group_of(M: float, peak_alpha_deg: float, mach_range: tuple) -> str:
    """The group of a run from its conditions: Mach number M [-], measured
    peak incidence [deg], and the Mach range (low, high) of the static
    stations."""
    if not mach_range[0] - gs.MACH_MARGIN <= M <= mach_range[1] + gs.MACH_MARGIN:
        return "outside_static_mach_range"
    return "primary" if peak_alpha_deg <= dm.PEAK_ALPHA_MAX_COMPARED else "beyond_static_lift_range"


def runs_table(set_name: str) -> pd.DataFrame:
    """The runs of one set ("calibration" or "held_out") from the inventory,
    with the group of each."""
    inv = pd.read_csv(fg.INV)
    inv = inv[inv["set"] == set_name].reset_index(drop=True)
    st = gs.station_table("static")
    rng = (st[0]["M"], st[-1]["M"])
    inv["group"] = [group_of(r.M, r.peak_alpha_deg, rng) for r in inv.itertuples()]
    return inv


# --------------------------------------------------------------------------- #
#  the model on one run
# --------------------------------------------------------------------------- #
def base_constants() -> dict:
    """Literature values, with the attached-flow moment factor measured on
    the NASA calibration loops (results/attached_moment_factor.json). The
    factor acts in proportion to the compressible weight, 0 at and below
    Mach 0.10; it is held fixed here."""
    base = dict(dm.DEFAULTS)
    if ATTACHED_FILE.exists():
        with open(ATTACHED_FILE) as fh:
            base["cm_unsteady"] = json.load(fh)["cm_unsteady"]
    return base


def solve_run(fr: dict, static: gs.GlasgowStatic, consts: dict, search: bool = False) -> dict:
    """The model on one run at its own measured conditions (mean incidence,
    amplitude, k, Mach number, chord 0.55 m, pitch axis at the quarter
    chord), with the static pressure drag at zero lift of its static inputs.
    `search=True` marches at the coarse resolution of the search stage."""
    n, cyc = dm.march_resolution(fr["k"], cal.SEARCH_DS, cal.SEARCH_SETTLE) if search else dm.march_resolution(fr["k"])
    return dm.solve(fr["a0"], fr["da"], fr["k"], fr["M"], consts={**consts, "CD0": static.CD0}, static=static,
                    n_per_cycle=n, n_cycles=cyc, chord=fg.CHORD_M, a_sound=fr["U"]/fr["M"],
                    x_pitch=fg.PITCH_AXIS_X_C)


def score_rows(x: np.ndarray, names: list, frames: list, statics: list, base: dict, search: bool = False) -> list:
    """unistall.metrics.score_frame of every run in `frames` with the
    constants `names` set to `x` over `base`."""
    c = {**base, **dict(zip(names, x, strict=True))}
    return [mt.score_frame(fr, solve_run(fr, sm, c, search), set_name="calibration")
            for fr, sm in zip(frames, statics, strict=True)]


def residuals(x: np.ndarray, names: list, frames: list, statics: list, base: dict, search: bool = False) -> np.ndarray:
    """The vector of normalised loop errors, unrounded: the cost of
    unistall/calibrate.py is the sum of its squares."""
    return np.array([r[q] for r in score_rows(x, names, frames, statics, base, search) for q in cal.QUANTITIES])


# --------------------------------------------------------------------------- #
#  cross-validated forecast (calibration-type runs only)
# --------------------------------------------------------------------------- #
_JOB = {}


def _best_fit(names: list, frames: list, statics: list, base: dict, rng: np.random.Generator, starts: int) -> np.ndarray:
    """Bounded least squares of the constants `names` from `starts` starts at
    the coarse march; the lowest cost is kept (as calibrate._cv_fold)."""
    if not names:
        return np.array([])
    spec = {f[0]: f for f in cal.FIT}
    lo, hi = np.array([spec[n][1] for n in names]), np.array([spec[n][2] for n in names])
    x0s = [np.array([dm.DEFAULTS[n] for n in names])] + [lo + (hi - lo)*rng.random(len(names)) for _ in range(starts - 1)]
    sols = [least_squares(residuals, np.clip(x0, lo, hi), bounds=(lo, hi), method="trf", x_scale=hi - lo,
                          diff_step=cal.DIFF_STEP, max_nfev=cal.MAX_EVALUATIONS, args=(names, frames, statics, base, True))
            for x0 in x0s]
    return min(sols, key=lambda r: cal.cost(r.fun)).x


def _sensitive(x: np.ndarray, names: list, frames: list, statics: list, base: dict, threshold: float,
               search: bool) -> list:
    """The constants whose +/- SENSITIVITY_PCT change moves the cost by at
    least `threshold` per cent."""
    j0 = cal.cost(residuals(x, names, frames, statics, base, search))
    keep = []
    for i, n in enumerate(names):
        change = []
        for pct in (-cal.SENSITIVITY_PCT, cal.SENSITIVITY_PCT):
            xi = x.copy()
            xi[i] *= 1.0 + pct/100.0
            change.append(abs(100.0*(cal.cost(residuals(xi, names, frames, statics, base, search)) - j0)/j0))
        if max(change) >= threshold:
            keep.append(n)
    return keep


def _cv_fold(fold: int) -> list:
    """One fold: every candidate fitted on the training runs only (the
    sensitive set decided from them alone), each left-out run scored once per
    candidate at the reporting resolution."""
    frames, statics, folds, base, thr = (_JOB[k] for k in ("frames", "statics", "folds", "base", "threshold"))
    test = sorted(folds[fold].tolist())
    train = [i for i in range(len(frames)) if i not in test]
    tr = ([frames[i] for i in train], [statics[i] for i in train])
    te = ([frames[i] for i in test], [statics[i] for i in test])
    rng = np.random.default_rng(cal.SEED + 1 + fold)
    every = [f[0] for f in cal.FIT]
    x_all = _best_fit(every, *tr, base, rng, cal.CV_STARTS)
    sens = _sensitive(x_all, every, *tr, base, thr, True)
    rows = []
    for label, names, x in (("all", every, x_all), ("sensitive only", sens, _best_fit(sens, *tr, base, rng, cal.CV_STARTS)),
                            ("none", [], np.array([]))):
        for r in score_rows(x, names, *te, base):
            rows.append(dict(constants_fitted=label, fold=fold, constants=" ".join(names),
                             **{n: round(float(v), 4) for n, v in zip(names, x, strict=True)}, **r))
    return rows


def forecast_table(static_source: str = "static") -> pd.DataFrame:
    """The left-out scores of every primary calibration-type run for each
    candidate, by CV_FOLDS-fold cross-validation with the seed and the folds
    rule of calibrate.select."""
    runs = runs_table("calibration")
    runs = runs[runs.group == "primary"].reset_index(drop=True)
    frames = [load(int(r), "forecast") for r in runs.run]
    statics = [gs.GlasgowStatic(fr["M"], source=static_source) for fr in frames]
    with open(NASA_TARGETS) as fh:
        thr = json.load(fh)["calibration"]["drop_constant_if_cost_change_below_pct"]
    perm = np.random.default_rng(cal.SEED).permutation(len(frames))
    _JOB.update(frames=frames, statics=statics, base=base_constants(), threshold=thr,
                folds=[perm[i::cal.CV_FOLDS] for i in range(cal.CV_FOLDS)])
    with Pool(min(PROCESSES, cal.CV_FOLDS)) as pool:
        rows = [row for fold in pool.map(_cv_fold, range(cal.CV_FOLDS), chunksize=1) for row in fold]
    cv = pd.DataFrame(rows)
    cv["left_out_cost"] = sum(cv[q]**2 for q in cal.QUANTITIES)
    cv["group"] = "primary"
    return cv


def choose(cv: pd.DataFrame) -> tuple:
    """The candidate kept by the one-standard-error rule of
    data/targets.json ("calibration", "selection_rule"), and the record of
    its steps."""
    by_fold = cv.pivot_table(index="fold", columns="constants_fitted", values="left_out_cost", aggfunc="sum")
    return cal._one_standard_error_choice(by_fold, list(CANDIDATES))


def verdict(chosen_rows: pd.DataFrame, targets: dict) -> dict:
    """The go rule of data/glasgow_targets.json applied to the left-out
    scores of the chosen candidate (`chosen_rows`, primary group)."""
    factor, held = targets["go_rule"]["factor"], targets["held_out_validation"]
    measures = {}
    for q in RULE_MEASURES:
        mean, limit = float(chosen_rows[q].mean()), factor*held[f"mean_{q}_max"]
        measures[f"mean_{q}"] = dict(forecast=round(mean, 4), target=held[f"mean_{q}_max"], limit=round(limit, 4),
                                     within=bool(mean <= limit))
    return dict(go=bool(all(m["within"] for m in measures.values())), factor=factor, n_runs=int(len(chosen_rows)),
                measures=measures)


def forecast(static_source: str = "static") -> dict:
    """Run the forecast and write its files; returns the verdict record."""
    suffix = "" if static_source == "static" else f"_{static_source}"
    cv = forecast_table(static_source)
    chosen, steps = choose(cv)
    cv["chosen"] = cv.constants_fitted == chosen
    cv.to_csv(RESULTS/f"glasgow_forecast{suffix}.csv", index=False)
    rows = cv[cv.chosen].reset_index(drop=True)
    summ = mt.summarise(rows, {"forecast_primary": np.ones(len(rows), bool)})
    summ.to_csv(RESULTS/f"glasgow_forecast_summary{suffix}.csv", index=False)
    with open(TARGETS) as fh:
        v = verdict(rows, json.load(fh))
    record = dict(author=AUTHOR, written_on=datetime.date.today().isoformat(), static_source=static_source,
                  chosen=chosen, selection_steps=steps,
                  left_out_cost={c: round(float(cv[cv.constants_fitted == c].left_out_cost.sum()), 6) for c in CANDIDATES},
                  takes_part_in_go_rule=bool(static_source == "static"), **v)
    with open(RESULTS/f"glasgow_verdict{suffix}.json", "w") as fh:
        json.dump(record, fh, indent=2)
        fh.write("\n")
    print(f"[glasgow-blind] forecast ({static_source}): kept '{chosen}'; "
          + "; ".join(f"{k} {m['forecast']} (limit {m['limit']})" for k, m in v["measures"].items())
          + f" -> {'GO' if v['go'] else 'NO-GO'}")
    return record


# --------------------------------------------------------------------------- #
#  the constants used for the blind scoring (calibration-type runs only)
# --------------------------------------------------------------------------- #
def _search(job: tuple) -> dict:
    i, x0 = job
    lo, hi = _JOB["lo"], _JOB["hi"]
    sol = least_squares(residuals, np.clip(x0, lo, hi), bounds=(lo, hi), method="trf", x_scale=hi - lo,
                        diff_step=cal.DIFF_STEP, max_nfev=cal.MAX_EVALUATIONS,
                        args=(_JOB["names"], _JOB["frames"], _JOB["statics"], _JOB["base"], True))
    return dict(start=i, x=sol.x, cost=cal.cost(sol.fun))


def _fit(names: list) -> tuple:
    """The constants `names` fitted on the runs in _JOB as each fold of the
    forecast fits them, with the starts of calibrate.fit: N_STARTS seeded
    starts, bounded least squares at the coarse march, the lowest cost kept.
    Returns (constants, search costs of every start)."""
    if not names:
        return {}, []
    spec = {f[0]: f for f in cal.FIT}
    lo, hi = np.array([spec[n][1] for n in names]), np.array([spec[n][2] for n in names])
    _JOB.update(names=names, lo=lo, hi=hi)
    rng = np.random.default_rng(cal.SEED)
    starts = [np.array([dm.DEFAULTS[n] for n in names])] + [lo + (hi - lo)*rng.random(len(names))
                                                            for _ in range(cal.N_STARTS - 1)]
    with Pool(PROCESSES) as pool:
        found = pool.map(_search, list(enumerate(starts)), chunksize=1)
    best = min(found, key=lambda r: r["cost"])
    return ({n: round(float(v), 6) for n, v in zip(names, best["x"], strict=True)},
            [round(float(r["cost"]), 6) for r in found])


def constants() -> dict:
    """Fit the set of constants the forecast kept on all primary
    calibration-type runs and write results/glasgow_constants.json. "all":
    every constant of calibrate.FIT; "sensitive only": those the sensitivity
    of that fit keeps, refitted; "none": literature values."""
    with open(VERDICT) as fh:
        chosen = json.load(fh)["chosen"]
    with open(NASA_TARGETS) as fh:
        thr = json.load(fh)["calibration"]["drop_constant_if_cost_change_below_pct"]
    runs = runs_table("calibration")
    runs = runs[runs.group == "primary"].reset_index(drop=True)
    frames = [load(int(r), "constants") for r in runs.run]
    statics = [gs.GlasgowStatic(fr["M"]) for fr in frames]
    base = base_constants()
    _JOB.update(frames=frames, statics=statics, base=base)
    every = [f[0] for f in cal.FIT]
    fitted, costs = ({}, []) if chosen == "none" else _fit(every)
    if chosen == "sensitive only":
        keep = _sensitive(np.array([fitted[n] for n in every]), every, frames, statics, base, thr, True)
        fitted, costs = _fit(keep)
    cost = cal.cost(residuals(np.array(list(fitted.values())), list(fitted), frames, statics, base))
    result = dict(author=AUTHOR, written_on=datetime.date.today().isoformat(), model="tabulated", constants=fitted,
                  fixed={k: v for k, v in base.items() if k not in fitted and k != "CD0"},
                  CD0="per run, the static pressure drag at zero lift of unistall.glasgow_static (not fitted)",
                  chosen=chosen, cost=round(cost, 8), search_costs=costs,
                  cost_definition="sum over runs of nRMS_CL^2 + nRMS_CM^2 + nRMS_CD^2, at the reporting resolution",
                  runs=[int(r) for r in runs.run], n_runs=int(len(runs)), static_source="static",
                  procedure="bounds, starts, seed, algorithm and search march of unistall.calibrate; fitted at the "
                            "search march, as every fold of the forecast is, and the cost evaluated once at the "
                            "reporting resolution")
    with open(CONSTANTS, "w") as fh:
        json.dump(result, fh, indent=2)
        fh.write("\n")
    print(f"[glasgow-blind] constants ({chosen}): {fitted}, cost {cost:.5f} on {len(runs)} runs")
    return result


# --------------------------------------------------------------------------- #
#  the seal
# --------------------------------------------------------------------------- #
def sealed_files() -> dict:
    """Every file that fixes the test, by a short name."""
    files = {"constants": CONSTANTS, "targets": TARGETS, "split": fg.SPLIT, "inventory": fg.INV,
             "static_inventory": fg.STATIC_INV, "forecast": FORECAST, "verdict": VERDICT}
    files.update({f"results/{n}": RESULTS/n for n in STATIC_FILES})
    files.update({f"unistall/{n}": PACKAGE/n for n in SOURCE_FILES})
    return files


def checksums(files: dict) -> dict:
    """SHA-256 of each file of `files` (name -> path). A missing file raises
    SealError."""
    missing = [n for n, p in files.items() if not Path(p).exists()]
    if missing:
        raise SealError(f"cannot seal or verify, missing: {missing}")
    return {n: hashlib.sha256(Path(p).read_bytes()).hexdigest() for n, p in files.items()}


def seal(files: dict | None = None, forecast_file: Path = FORECAST, targets_file: Path = TARGETS,
         seal_file: Path = SEAL) -> dict:
    """Write the seal. The verdict is computed here from the forecast file
    and the targets file; NoGo is raised, and nothing written, if it is
    no-go. An existing seal is not replaced."""
    if Path(seal_file).exists():
        raise SealError(f"{seal_file} exists: a seal is written once")
    files = sealed_files() if files is None else files
    sums = checksums(files)
    cv = pd.read_csv(forecast_file)
    with open(targets_file) as fh:
        v = verdict(cv[cv.chosen & (cv.group == "primary")], json.load(fh))
    if not v["go"]:
        raise NoGo("the forecast does not meet the go rule; no seal is written and the held-out runs are not scored: "
                   + json.dumps(v["measures"]))
    record = dict(author=AUTHOR, sealed_on=datetime.date.today().isoformat(),
                  statement="No held-out Glasgow run had been scored when this seal was written.",
                  verdict=v, sha256=sums, scored_on=None)
    with open(seal_file, "w") as fh:
        json.dump(record, fh, indent=2)
        fh.write("\n")
    return record


def check_seal(files: dict | None = None, seal_file: Path = SEAL, scores_file: Path = SCORES) -> dict:
    """The seal record, if the scoring may proceed. Raises SealError if the
    seal is missing, a sealed file has changed, the verdict is not go, the
    scoring date is set, or the score file exists."""
    if not Path(seal_file).exists():
        raise SealError(f"{seal_file} does not exist: run seal first")
    with open(seal_file) as fh:
        record = json.load(fh)
    now = checksums(sealed_files() if files is None else files)
    changed = sorted(n for n in set(now) | set(record["sha256"]) if now.get(n) != record["sha256"].get(n))
    if changed:
        raise SealError(f"changed since the seal: {changed}")
    if not record["verdict"]["go"]:
        raise SealError("the sealed verdict is not go")
    if record.get("scored_on") or Path(scores_file).exists():
        raise SealError("the held-out runs have been scored already; they are scored once")
    return record


# --------------------------------------------------------------------------- #
#  the single scoring
# --------------------------------------------------------------------------- #
def _cache_matches_inventory(runs: pd.DataFrame) -> None:
    """The cached files of `runs` have the SHA-256 the inventory records."""
    bad = [int(r.run) for r in runs.itertuples()
           if fg.sha256_of(fg.RUNS/f"{r.run}.dat") != r.sha256_dat or fg.sha256_of(fg.RUNS/f"{r.run}_coeffs.dat") != r.sha256_coeffs]
    if bad:
        raise SealError(f"cached files differ from the inventory for runs {bad}")


def targets_table(summary: pd.DataFrame, held: dict) -> pd.DataFrame:
    """Each target against the value of the primary group in `summary`
    (metrics.summarise)."""
    rows = []
    for metric, key, kind in TARGET_MEASURES:
        s = summary[(summary.group == "primary") & (summary.metric == metric)].iloc[0]
        met = s.value <= held[key] if kind == "at most" else s.value >= held[key]
        rows.append(dict(group="primary", metric=metric, value=s.value, ci95_lo=s.ci95_lo, ci95_hi=s.ci95_hi,
                         n_frames=int(s.n_frames), target=held[key], kind=kind, met=bool(met)))
    return pd.DataFrame(rows)


def score() -> pd.DataFrame:
    """Score the held-out runs, once (see the module text)."""
    record = check_seal()
    runs = runs_table("held_out")
    _cache_matches_inventory(runs)
    consts = dm.load_constants(CONSTANTS)
    with open(TARGETS) as fh:
        held = json.load(fh)["held_out_validation"]
    rows = []
    with _scoring_open():
        for r in runs.itertuples():
            fr = load(int(r.run), "score")
            row = mt.score_frame(fr, solve_run(fr, gs.GlasgowStatic(fr["M"]), consts), set_name="held_out")
            rows.append(dict(**row, group=r.group, peak_alpha_deg=r.peak_alpha_deg, Re=r.Re))
    table = pd.DataFrame(rows)
    table.to_csv(SCORES, index=False)
    summ = mt.summarise(table, {g: (table.group == g).values for g in ("primary", "beyond_static_lift_range")
                                if (table.group == g).any()})
    out = summ.merge(targets_table(summ, held)[["group", "metric", "target", "kind", "met"]],
                     on=["group", "metric"], how="left")
    out.to_csv(SUMMARY, index=False)
    record["scored_on"] = datetime.date.today().isoformat()
    record["n_runs_scored"] = int(len(table))
    with open(SEAL, "w") as fh:
        json.dump(record, fh, indent=2)
        fh.write("\n")
    print(out[out.target.notna()].to_string(index=False))
    return out


# --------------------------------------------------------------------------- #
#  the targets file
# --------------------------------------------------------------------------- #
def go_factor() -> dict:
    """The factor of the go rule, from the NASA results as they stand: how the
    held-out error of the NASA loops compared with the error forecast for
    them by cross-validation, and how wide the interval of a mean loop error
    is. See data/glasgow_targets.json, "go_rule", for the reading."""
    q = list(cal.QUANTITIES)
    cv = pd.read_csv(RESULTS/"cross_validation.csv")
    with open(RESULTS/"calibrated_constants.json") as fh:
        nasa = json.load(fh)
    chosen, n_cal = nasa["selection"]["chosen"], int(nasa["n_frames"])
    left_out = float(np.sqrt(cv[cv.constants_fitted == chosen].left_out_cost.sum()/(len(q)*n_cal)))
    tab = pd.read_csv(RESULTS/"validation_table.csv")
    prim = tab[(tab.model == "tabulated") & (tab["set"] == "held_out") & tab.in_mach_range]
    held_out = float(np.sqrt((prim[q]**2).sum(axis=1).mean()/len(q)))
    s = pd.read_csv(RESULTS/"validation_summary.csv")
    s = s[(s.model == "tabulated") & (s.group == "primary") & s.metric.isin([f"mean_{m}" for m in q])]
    half = float(((s.ci95_hi - s.ci95_lo)/(2.0*s.value)).max())
    ratio = held_out/left_out
    return dict(factor=round((1.0 + half)/ratio, 2), held_out_over_left_out=round(ratio, 4),
                nasa_left_out_rms_error=round(left_out, 4), nasa_held_out_rms_error=round(held_out, 4),
                nasa_left_out_candidate=chosen, nasa_calibration_loops=n_cal, nasa_held_out_loops=int(len(prim)),
                widest_relative_half_interval_of_a_mean=round(half, 4),
                files_used={n: hashlib.sha256((RESULTS/n).read_bytes()).hexdigest()
                            for n in ("cross_validation.csv", "calibrated_constants.json", "validation_table.csv",
                                      "validation_summary.csv")})


def _procedure() -> dict:
    return dict(
        split="data/glasgow_split.json: run number divisible by 3 -> calibration, otherwise held_out",
        conditions="every run is solved at its own measured mean incidence and amplitude ((max+min)/2 and (max-min)/2 "
                   "of the incidence column), its stored reduced frequency and Mach number, chord 0.55 m, pitch axis at "
                   "the quarter chord, as a pure sine",
        measured_loads="C_l and pressure C_d formed from the stored C_n and C_t; C_m about the quarter chord as stored; "
                       "128 points per run, the average of 10 cycles",
        static_inputs="unistall/glasgow_static.py, source 'static': the static runs of the same model dated from "
                      f"{gs.PERIOD_START} on (doi:{fg.DEPOSIT_DOI}), three Mach stations, interpolated in Mach number; "
                      "slope, zero-lift incidence, static-stall incidence, C_N1, f(alpha), static C_m and the pressure "
                      "C_D0 by the rule of static_stations.py and static_model.py; post-stall curve 'more' separated, "
                      "the default of the model. None is fitted. No held-out run enters them.",
        constants="literature values of dsmodel.DEFAULTS and the attached-flow moment factor of "
                  "results/attached_moment_factor.json (NASA calibration loops) held fixed; C_D0 from the static runs; "
                  "Tp, Tf0, Tv0, Tvl, eta candidates for fitting within the bounds of calibrate.FIT",
        fitted_on="the primary calibration-type runs only",
        cost="sum over runs of nRMS_CL^2 + nRMS_CM^2 + nRMS_CD^2 (unistall/metrics.py), as unistall/calibrate.py",
        selection=f"{cal.CV_FOLDS}-fold cross-validation inside the calibration-type runs, seed {cal.SEED}, "
                  f"{cal.CV_STARTS} starts a fold, candidates none / sensitive only / all, the one-standard-error "
                  "rule of data/targets.json ('calibration', 'selection_rule')",
        forecast="the left-out scores of the kept candidate: each primary calibration-type run scored once, at the "
                 "reporting resolution, by constants fitted without it",
        final_constants=f"the kept set refitted on all primary calibration-type runs from {cal.N_STARTS} seeded starts at "
                        f"the search march of unistall.calibrate (step {cal.SEARCH_DS} semichords), the lowest cost "
                        "kept. This is how every fold of the forecast fits, so the forecast describes these "
                        "constants. The refinement of every start at the reporting resolution that calibrate.fit adds "
                        "is left out: one evaluation of the cost on these runs at that resolution is 1.9 million "
                        "steps of the march. Every score, in the forecast and in the blind test, is at the reporting "
                        "resolution.",
        measures="those of data/targets.json 'held_out_validation', computed by unistall/metrics.py: mean nRMS of C_L, "
                 "C_M, C_D; mean absolute error of the incidence of C_L,max and of moment stall; cycle-damping sign "
                 "agreement; means over runs with bootstrap 95 % intervals",
        scoring="unistall.glasgow_blind score, once, after the seal; a no-go verdict means it is never run",
        static_source_uncertainty="the forecast is repeated with the 'quasi_steady' static source on the "
                                  "calibration-type runs only; it takes no part in the go rule and no held-out run is "
                                  "scored with it")


def write_targets(model_results_existing: str) -> dict:
    """Write data/glasgow_targets.json. `model_results_existing` states what
    model results existed for Glasgow runs at that moment. Refuses to replace
    the file."""
    if TARGETS.exists():
        raise SealError(f"{TARGETS} exists: the targets are written once")
    with open(NASA_TARGETS) as fh:
        nasa = json.load(fh)
    held = dict(nasa["held_out_validation"])
    held["group"] = ("held-out Glasgow NACA 0012 runs inside the Mach range of the static stations whose measured peak "
                     "incidence is at or below 25 deg, the end of the static lift data (group 'primary')")
    g = go_factor()
    f = g["factor"]
    record = dict(
        author=AUTHOR, frozen_on=datetime.date.today().isoformat(),
        rule="Targets, procedure and go rule are stated here before any held-out Glasgow run is scored and are not "
             "changed after a result is seen.",
        statement="No held-out Glasgow run had been scored, plotted against the model or otherwise evaluated when this "
                  "file was written. " + model_results_existing,
        targets_copied_from="data/targets.json, 'held_out_validation' (same measures, same limits)",
        held_out_validation=held,
        procedure=_procedure(),
        go_rule=dict(
            status="binding",
            statement=f"GO if, in the cross-validated forecast on the primary calibration-type runs (static source "
                      f"'static', the candidate kept by the selection rule), the mean left-out nRMS of C_L, of C_M and "
                      f"of C_D are each at most {f} times the corresponding target: C_L {round(f*held['mean_nRMS_CL_max'], 4)}, "
                      f"C_M {round(f*held['mean_nRMS_CM_max'], 4)}, C_D {round(f*held['mean_nRMS_CD_max'], 4)}. "
                      "Otherwise NO-GO.",
            factor=f, measures=[f"mean_{m}" for m in RULE_MEASURES],
            if_go="the held-out runs are scored once by unistall.glasgow_blind score, after the seal, and the result "
                  "is reported against the targets whatever it is",
            if_no_go="the held-out runs are never scored; the case is reported as a stated limit of the model, with "
                     "the forecast numbers",
            reason_for_the_factor="The forecast is taken as near the targets if the held-out means it implies could "
                                  "still meet them. On the NASA loops the held-out error was "
                                  f"{g['held_out_over_left_out']} of the error forecast by cross-validation (pooled "
                                  "RMS of the three normalised loop errors), and the 95 % interval of a mean loop "
                                  f"error was up to {round(100*g['widest_relative_half_interval_of_a_mean'], 1)} % of "
                                  "the mean either side. factor = (1 + that half-width) / that ratio, rounded to two "
                                  "decimals. The measures of incidence and of damping are reported in the forecast "
                                  "but are not in the rule: the fit does not minimise them.",
            basis=g))
    with open(TARGETS, "w") as fh:
        json.dump(record, fh, indent=2)
        fh.write("\n")
    return record


COMMANDS = {"forecast": lambda: forecast(next((a.split("=")[1] for a in sys.argv if a.startswith("--static=")), "static")),
            "constants": constants, "seal": seal, "score": score}

if __name__ == "__main__":
    command = sys.argv[1] if len(sys.argv) > 1 else ""
    if command == "targets":
        write_targets(" ".join(sys.argv[2:]))
    elif command in COMMANDS:
        COMMANDS[command]()
    else:
        sys.exit(__doc__)
