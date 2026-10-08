# -*- coding: utf-8 -*-
"""
unistall / calibrate.py
-------------------
Fit the dynamic constants of the load model (unistall/dsmodel.py) on the
CALIBRATION frames, by a stated cost and a stated algorithm, and write the
result to results/calibrated_constants.json, which the model reads.

Author: Akosa Samuel Onyejekwe (independent)

WHAT IS FITTED. Only the constants in FIT below. Slope, zero-lift angle, the
onset level C_N1 and the separation function come from static measurements
and are not touched. Every other constant keeps the value in dsmodel.DEFAULTS
and is listed in the output as "literature / fixed".

DATA. The usable dynamic NACA 0012 frames of the calibration set that lie
inside the Mach range covered by static data (data/data_manifest.csv). No
held-out frame is opened by this script.

COST. For each frame, the RMS error in C_L, C_M and C_D over the measured
loop, each divided by the measured range of that quantity (unistall/metrics.py).
The cost is the sum of the squares of those three numbers over all frames:
        J = sum_frames ( nRMS_CL^2 + nRMS_CM^2 + nRMS_CD^2 ).

ALGORITHM. Bounded nonlinear least squares (scipy.optimize.least_squares,
trust-region reflective) on the vector of those normalised errors, unrounded.
SEARCH: from N_STARTS points (the literature values, then points drawn
uniformly inside the bounds by numpy's default generator with seed SEED), at
a coarser march (step SEARCH_DS semichords). POLISH: every search result is
refined at the resolution every reported result uses, and the lowest is the
calibration. All starts are written out, so the
spread between them is visible.

Usage:
  python3 -m unistall.calibrate --attached      measure the unsteady-moment factor on the
                                            attached-flow calibration loops (run first)
  python3 -m unistall.calibrate --fit           fit and write calibrated_constants.json
  python3 -m unistall.calibrate --sensitivity   +/- 10 % on each fitted constant ->
                                            results/sensitivity.csv
  add --model=reference to any of them to calibrate the reference model instead
  python3 -m unistall.calibrate --onset         cost against the stall-onset level (run after --select)
  python3 -m unistall.calibrate --curve         cost with each post-stall static curve
  python3 -m unistall.calibrate --repeat        run the first start again and compare
  python3 -m unistall.calibrate --select        choose which constants to keep fitted by
                                            cross-validation inside the calibration
                                            set (run after the two above)
"""
import sys, json
from pathlib import Path
import numpy as np
import pandas as pd
from multiprocessing import Pool
from scipy.optimize import least_squares, minimize_scalar

HERE = Path(__file__).resolve().parent
from unistall import dsmodel as dm
from unistall import metrics as mt
from unistall.static_model import StaticModel
from unistall.paths import DATA, RESULTS

#        name    lower  upper
FIT = [("Tp",    0.3,   4.0),
       ("Tf0",   0.5,  10.0),
       ("Tv0",   1.0,  15.0),
       ("Tvl",   3.0,  13.0),
       ("eta",   0.6,   1.0)]
SEED = 20261006
N_STARTS = 8
SEARCH_DS, SEARCH_SETTLE = 0.32, 120.0  # march of the search stage: step and settling travel, semichords
POLISH_STARTS = N_STARTS                # every search result is refined at full resolution
DIFF_STEP = 1e-2                        # relative step of the finite-difference Jacobian
MAX_EVALUATIONS = 150                   # per start, per stage
CV_FOLDS, CV_STARTS, PROCESSES = 5, 2, 4
SENSITIVITY_PCT = 10
QUANTITIES = ("nRMS_CL", "nRMS_CM", "nRMS_CD")

# Which model is being calibrated. "tabulated" is unistall/dsmodel.py with the
# measured separation table; "reference" is the same model with the fitted
# exponential separation law (unistall/reference_lb.py). Both go through exactly
# the same cost, algorithm, starts, bounds and frames.
MODEL = "reference" if "--model=reference" in sys.argv else "tabulated"
SUFFIX = "" if MODEL == "tabulated" else "_reference"
CONSTANTS_FILE = RESULTS/f"calibrated_constants{SUFFIX}.json"
ATTACHED_FILE = RESULTS/"attached_moment_factor.json"
ATTACHED_BOUNDS = (0.3, 1.2)


def _static(M):
    if MODEL == "reference":
        from unistall.reference_lb import ExponentialStatic
        return ExponentialStatic(M)
    return StaticModel(M)


def calibration_frames():
    m = mt.manifest_frames(sets=("calibration",))
    m = m[m.within_static_mach_range].reset_index(drop=True)
    frames = [mt.load_frame(mt.FRAME_CACHE/f"{f}.mat") for f in m.frame]
    statics = [_static(fr["M"]) for fr in frames]
    return m, frames, statics


def residuals(x, names, frames, statics, base=None, search=False):
    """The vector of normalised loop errors, unrounded. `search=True` marches
    at the coarser resolution of the search stage; otherwise at the resolution
    every reported result uses (dsmodel.march_resolution)."""
    c = dict(base or dm.DEFAULTS)
    c.update(dict(zip(names, x, strict=False)))
    out = []
    for fr, sm in zip(frames, statics, strict=False):
        n, cyc = dm.march_resolution(fr["k"], SEARCH_DS, SEARCH_SETTLE) if search else dm.march_resolution(fr["k"])
        o = dm.solve(fr["a0"], fr["da"], fr["k"], fr["M"], consts=c, static=sm, n_per_cycle=n, n_cycles=cyc)
        r = mt.score_frame(fr, o)
        out += [r[q] for q in QUANTITIES]
    return np.array(out)


def cost(res):
    return float(np.sum(res**2))


def _unlagged_static(out_full, out_static, sm):
    """The model moment with the static table read at the geometric incidence
    instead of the lagged effective one: the unsteady terms of `out_full`
    (its moment less that of `out_static`, the run with those terms removed)
    added to the static moment at the incidence itself."""
    o = dict(out_full)
    static_here = np.array([sm.cm_at(float(a)) for a in out_full["alpha_deg"]])
    o["CM"] = np.asarray(out_full["CM"]) - np.asarray(out_static["CM"]) + static_here
    return o


ONSET_MARGIN_DEG = 1.0        # up-stroke points closer than this to the model's onset incidence are left out
BOOT_FACTOR = 2000            # bootstrap resamples of the loops for the interval of the factor


def _attached_portion_fit(frames: list, statics: list) -> dict:
    """The unsteady-moment factor from the attached part of every calibration
    loop, with a bootstrap interval over loops.

    For each loop the model is run with literature stall constants. The
    measured moment points used are those on the up-stroke at least
    ONSET_MARGIN_DEG below the incidence at which the model reaches stall
    onset, or every point of the loop if the model never reaches onset on it.
    The moment is linear in the factor (C_M = C_M at factor 0 + factor times
    the unsteady terms), so two runs a loop give the least-squares factor in
    closed form; each loop has equal weight and its errors are divided by its
    measured moment range. Returns the factor, its 95 % interval, and the
    loops and points used."""
    rng = np.random.default_rng(SEED)
    num, den, used = [], [], []
    for fr, sm in zip(frames, statics, strict=True):
        if dm.af.compressible_weight(fr["M"]) < 1.0:
            continue                                        # the factor acts in full only on the compressible branch
        out = [dm.solve(fr["a0"], fr["da"], fr["k"], fr["M"], consts=dict(dm.DEFAULTS, cm_unsteady=g_), static=sm)
               for g_ in (0.0, 1.0)]
        onset = out[1]["onset_alpha_deg"]
        stroke = mt.stroke_split(fr["acm"])
        keep = np.ones(len(fr["acm"]), bool) if not np.isfinite(onset) else \
            (stroke == "up") & (fr["acm"] <= onset - ONSET_MARGIN_DEG)
        if keep.sum() < 3:
            continue
        m0, m1 = ([np.interp(a, *mt._branches(o["alpha_deg"], o["alpha_dot"], o["CM"])[st])
                   for a, st in zip(fr["acm"][keep], stroke[keep], strict=True)] for o in out)
        m0, d = np.array(m0), np.array(m1) - np.array(m0)
        r = float(np.ptp(fr["cm"]))
        num.append(float(np.sum(d*(fr["cm"][keep] - m0))/(keep.sum()*r*r)))
        den.append(float(np.sum(d*d)/(keep.sum()*r*r)))
        used.append(dict(frame=fr["frame"], k=round(fr["k"], 3), points=int(keep.sum()),
                         model_reaches_onset=bool(np.isfinite(onset)),
                         factor_this_loop=round(num[-1]/den[-1], 3) if den[-1] > 0.0 else None))
    num, den = np.array(num), np.array(den)
    pick = rng.integers(0, len(num), size=(BOOT_FACTOR, len(num)))
    boot = num[pick].sum(axis=1)/den[pick].sum(axis=1)
    lo, hi = np.percentile(boot, [2.5, 97.5])
    return dict(cm_unsteady=round(float(num.sum()/den.sum()), 3), interval_95=[round(float(lo), 3), round(float(hi), 3)],
                n_loops=len(used), n_points=int(sum(u["points"] for u in used)),
                onset_margin_deg=ONSET_MARGIN_DEG, loops=used)


def attached_moment():
    """An empirical factor on the unsteady moment, from attached flow.

    Before stall onset the model, with its unsteady moment terms at full
    strength, gives a wider moment loop than is measured. The factor
    cm_unsteady is the single number that minimises the moment error on the
    attached part of every calibration loop (_attached_portion_fit), given
    with a bootstrap interval over loops. It is then held fixed while the
    stall constants are fitted. The loops whose peak incidence is below static
    stall are also examined one by one, as before. It is a correction
    of this model, not a measurement of the flow, and what it corrects is
    measured here loop by loop:

      damping_static_part      damping the model gives with the unsteady
                               moment terms removed; it comes only from
                               reading the static moment table at the lagged
                               effective incidence
      damping_unlagged_static  damping with the terms at full strength and
                               the static table read at the incidence itself
      best_factor_this_loop    the factor that this loop alone would choose
      peak_vortex_normal_force, peak_CN_prime_over_CN1
                               whether the model's leading-edge criterion is
                               reached on a loop that is below static stall
    """
    m, frames, _statics = calibration_frames()
    statics = [StaticModel(fr["M"]) for fr in frames]
    pick = [i for i, (fr, sm) in enumerate(zip(frames, statics, strict=False))
            if abs(fr["a0"]) + abs(fr["da"]) < sm.alpha_stall_deg]

    def run(i, g):
        c = dict(dm.DEFAULTS, cm_unsteady=g)
        return dm.solve(frames[i]["a0"], frames[i]["da"], frames[i]["k"], frames[i]["M"], consts=c, static=statics[i])

    def scores(g):
        return [mt.score_frame(frames[i], run(i, g)) for i in pick]

    def best_for(only):
        r = minimize_scalar(lambda g: sum(mt.score_frame(frames[i], run(i, g))["nRMS_CM"]**2 for i in only),
                            bounds=ATTACHED_BOUNDS, method="bounded", options=dict(xatol=1e-3))
        return round(float(r.x), 3)
    base = _attached_portion_fit(frames, statics)
    g = base["cm_unsteady"]
    before, after, none = scores(1.0), scores(g), scores(0.0)
    full = [run(i, 1.0) for i in pick]
    unlagged = [mt.score_frame(frames[i], _unlagged_static(o, run(i, 0.0), statics[i]))
                for i, o in zip(pick, full, strict=True)]
    own = [best_for([i]) for i in pick]
    closed = [f for f, a in zip(own, after, strict=True) if a["moment_loop_closed"]]

    def num(v):
        return None if not np.isfinite(v) else round(float(v), 4)
    result = dict(author="Akosa Samuel Onyejekwe (independent)", cm_unsteady=g, bounds=list(ATTACHED_BOUNDS),
                  rule="the attached part of every calibration loop: measured moment points that the model, with "
                       "literature stall constants, places before stall onset (see basis)",
                  basis=base, factor_from_the_loops_below_static_stall_alone=best_for(pick),
                  range_over_loops_with_a_closed_moment_loop=[min(closed), max(closed)],
                  frames=[dict(frame=m.frame[i], M=frames[i]["M"], k=frames[i]["k"],
                               moment_loop_closed=bool(a["moment_loop_closed"]),
                               damping_measured=num(a["Xi_exp"]), damping_without_factor=num(b["Xi_model"]),
                               damping_with_factor=num(a["Xi_model"]), damping_static_part=num(z["Xi_model"]),
                               damping_unlagged_static=num(u["Xi_model"]),
                               nRMS_CM_without_factor=num(b["nRMS_CM"]), nRMS_CM_with_factor=num(a["nRMS_CM"]),
                               nRMS_CM_unlagged_static=num(u["nRMS_CM"]), best_factor_this_loop=f,
                               peak_vortex_normal_force=num(np.abs(o["CN_vortex"]).max()),
                               peak_CN_prime_over_CN1=num(np.abs(o["CN_prime"]).max()/statics[i].CN1))
                          for i, b, a, z, u, f, o in zip(pick, before, after, none, unlagged, own, full, strict=True)])
    json.dump(result, open(ATTACHED_FILE, "w"), indent=2)
    open(ATTACHED_FILE, "a").write("\n")
    print(f"[calibrate] unsteady-moment factor from the attached part of {base['n_loops']} loops "
          f"({base['n_points']} points): {g}, 95 % interval {base['interval_95'][0]} to {base['interval_95'][1]}; "
          f"the {len(pick)} loops below static stall alone give {result['factor_from_the_loops_below_static_stall_alone']} "
          f"(closed loops singly {min(closed)} to {max(closed)})")
    for r in result["frames"]:
        print("            " + ", ".join(f"{k} {v}" for k, v in r.items()))
    return result


def base_constants():
    """Literature values, with the attached-flow moment factor if measured."""
    base = dict(dm.DEFAULTS)
    if ATTACHED_FILE.exists():
        base["cm_unsteady"] = json.load(open(ATTACHED_FILE))["cm_unsteady"]
    return base


_FIT = {}


def _least_squares(x0, search):
    names, lo, hi = _FIT["names"], _FIT["lo"], _FIT["hi"]
    return least_squares(residuals, np.clip(x0, lo, hi), bounds=(lo, hi), method="trf", x_scale=hi - lo,
                         diff_step=DIFF_STEP, max_nfev=MAX_EVALUATIONS,
                         args=(names, _FIT["frames"], _FIT["statics"], _FIT["base"], search))


def _search(job):
    i, x0 = job
    sol = _least_squares(x0, True)
    return dict(start=i, x0=x0, x=sol.x, cost=cost(sol.fun), evaluations=int(sol.nfev), status=int(sol.status))


def _polish(run):
    sol = _least_squares(run["x"], False)
    return dict(start=run["start"], x=sol.x, cost=cost(sol.fun), evaluations=int(sol.nfev), status=int(sol.status))


def fit(names=None, base=None, write=True):
    """Two stages. SEARCH: bounded least squares from N_STARTS seeded starts at
    the coarser march, to find the basins. POLISH: every one is refined at
    the resolution every reported result uses; the lowest is the calibration."""
    base = base or base_constants()
    spec = [f for f in FIT if names is None or f[0] in names]
    names = [f[0] for f in spec]
    lo = np.array([f[1] for f in spec])
    hi = np.array([f[2] for f in spec])
    m, frames, statics = calibration_frames()
    _FIT.update(names=names, lo=lo, hi=hi, frames=frames, statics=statics, base=base)
    rng = np.random.default_rng(SEED)
    starts = [np.array([dm.DEFAULTS[n] for n in names])]
    starts += [lo + (hi - lo)*rng.random(len(names)) for _ in range(N_STARTS - 1)]
    with Pool(PROCESSES) as pool:
        found = pool.map(_search, list(enumerate(starts)), chunksize=1)
        best_two = sorted(found, key=lambda r: r["cost"])[:POLISH_STARTS]
        polished = {r["start"]: r for r in pool.map(_polish, best_two, chunksize=1)}
    rows = []
    for r in found:
        pol = polished.get(r["start"])
        rows.append(dict(start=r["start"], **{f"start_{n}": round(float(v), 4) for n, v in zip(names, r["x0"], strict=False)},
                         **{f"search_{n}": float(v) for n, v in zip(names, r["x"], strict=False)},
                         search_cost=r["cost"], search_evaluations=r["evaluations"], search_status=r["status"],
                         polished=pol is not None,
                         **{n: float(v) for n, v in zip(names, (pol or r)["x"], strict=False)},
                         cost=pol["cost"] if pol else np.nan,
                         polish_evaluations=pol["evaluations"] if pol else 0))
        print(f"[calibrate] start {r['start']}: search cost {r['cost']:.5f}"
              + (f", polished {pol['cost']:.5f}" if pol else "") + " at "
              + ", ".join(f"{n} {v:.4f}" for n, v in zip(names, (pol or r)["x"], strict=False)), flush=True)
    runs = pd.DataFrame(rows)
    best = runs.loc[runs.cost.idxmin()]
    x = np.array([best[n] for n in names])
    c0 = cost(residuals(np.array([dm.DEFAULTS[n] for n in names]), names, frames, statics, base))
    result = dict(
        author="Akosa Samuel Onyejekwe (independent)", model=MODEL,
        constants={n: round(float(v), 6) for n, v in zip(names, x, strict=False)},
        fixed={k: v for k, v in base.items() if k not in names},
        cost=round(float(best.cost), 8), cost_at_literature_values=round(c0, 8),
        cost_definition="sum over frames of nRMS_CL^2 + nRMS_CM^2 + nRMS_CD^2, at the reporting resolution",
        algorithm="scipy.optimize.least_squares, method trf, bounded; search at a coarser march, "
                  "then every start refined at the reporting resolution",
        starts=N_STARTS, seed=SEED, bounds={f[0]: [f[1], f[2]] for f in spec},
        march=dict(search_step_semichords=SEARCH_DS, search_settle_semichords=SEARCH_SETTLE,
                   reporting_step_semichords=dm.DS_MAX, reporting_settle_semichords=dm.SETTLE_SEMICHORDS,
                   jacobian_relative_step=DIFF_STEP),
        frames=list(m.frame), n_frames=int(len(m)),
        spread_over_starts={n: [round(float(runs[f"search_{n}"].min()), 4), round(float(runs[f"search_{n}"].max()), 4)]
                            for n in names},
        spread_of_search_cost=[round(float(runs.search_cost.min()), 5), round(float(runs.search_cost.max()), 5)],
        starts_within_1pct_of_best_search=int((runs.search_cost <= 1.01*runs.search_cost.min()).sum()),
        polished_costs=[round(float(v), 6) for v in runs.cost.dropna()],
        polished_constants_differ_by={n: round(float(runs[runs.polished][n].max() - runs[runs.polished][n].min()), 4)
                                      for n in names})
    if write:
        runs.to_csv(RESULTS/f"calibration_runs{SUFFIX}.csv", index=False)
        json.dump(result, open(CONSTANTS_FILE, "w"), indent=2)
        open(CONSTANTS_FILE, "a").write("\n")
        print(f"[calibrate] best cost {best.cost:.5f} (literature values {c0:.5f}); search costs "
              f"{result['spread_of_search_cost']}")
    return result, runs


def sensitivity():
    cal = json.load(open(CONSTANTS_FILE))
    names = list(cal["constants"])
    x = np.array([cal["constants"][n] for n in names])
    base = {**dm.DEFAULTS, **cal["fixed"]}
    m, frames, statics = calibration_frames()
    r0 = residuals(x, names, frames, statics, base)
    j0 = cost(r0)
    def mean(r):
        return {q: float(np.mean(r[i::3])) for i, q in enumerate(QUANTITIES)}
    rows = [dict(constant="(calibrated)", change_pct=0, cost=round(j0, 6), cost_change_pct=0.0,
                 **{f"mean_{q}": round(v, 4) for q, v in mean(r0).items()})]
    for i, name in enumerate(names):
        for pct in (-10, 10):
            xi = x.copy()
            xi[i] *= 1.0 + pct/100.0
            r = residuals(xi, names, frames, statics, base)
            j = cost(r)
            rows.append(dict(constant=name, change_pct=pct, cost=round(j, 6),
                             cost_change_pct=round(100.0*(j - j0)/j0, 3),
                             **{f"mean_{q}": round(v, 4) for q, v in mean(r).items()}))
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS/f"sensitivity{SUFFIX}.csv", index=False)
    worst = df[df.constant != "(calibrated)"].groupby("constant").cost_change_pct.apply(lambda s: s.abs().max())
    print("[calibrate] largest cost change for +/- 10 %:\n" + worst.round(3).to_string())
    return df, worst


def curve_choice():
    """Which post-stall static curve the model follows: the cost on the
    calibration loops with each, at the literature stall constants."""
    m, frames, _ = calibration_frames()
    base = base_constants()
    rows = []
    for curve in ("more", "less"):
        statics = [StaticModel(fr["M"], curve) for fr in frames]
        r = residuals(np.array([]), [], frames, statics, base).reshape(-1, len(QUANTITIES))
        rows.append(dict(post_stall_curve=curve + " separated", cost=round(cost(r), 5),
                         **{f"mean_{q}": round(float(v), 4) for q, v in zip(QUANTITIES, r.mean(axis=0), strict=False)}))
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS/"static_curve_choice.csv", index=False)
    print(df.to_string(index=False))
    return df


ONSET_FACTORS = (0.92, 0.96, 1.0, 1.04, 1.08)


def onset_level():
    """Is the stall-onset level, taken from the static-stall normal force,
    where the dynamic loops want it? The calibration cost with the onset
    level multiplied by each factor in ONSET_FACTORS, every other constant as
    calibrated. A minimum at 1.0 says the static value needs no adjustment."""
    cal = json.load(open(CONSTANTS_FILE))
    m, frames, statics = calibration_frames()
    base = {**dm.DEFAULTS, **cal["fixed"], **cal["constants"]}
    rows = []
    for f in ONSET_FACTORS:
        r = residuals(np.array([f]), ["cn1_factor"], frames, statics, base).reshape(-1, len(QUANTITIES))
        rows.append(dict(onset_level_factor=f, cost=round(cost(r), 5),
                         **{f"mean_{q}": round(float(v), 4) for q, v in zip(QUANTITIES, r.mean(axis=0), strict=False)}))
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS/f"onset_level_check{SUFFIX}.csv", index=False)
    print(df.to_string(index=False))
    return df


def repeatability():
    """Run the first search start again and compare with the recorded run."""
    T = json.load(open(DATA/"targets.json"))["calibration"]["repeatability_max_abs_difference"]
    runs = pd.read_csv(RESULTS/f"calibration_runs{SUFFIX}.csv").set_index("start")
    cal = json.load(open(CONSTANTS_FILE))
    cal = cal.get("first_stage_record", cal)
    names = list(cal["constants"])
    spec = {f[0]: f for f in FIT}
    m, frames, statics = calibration_frames()
    _FIT.update(names=names, lo=np.array([spec[n][1] for n in names]), hi=np.array([spec[n][2] for n in names]),
                frames=frames, statics=statics, base={**dm.DEFAULTS, **cal["fixed"]})
    again = _search((0, np.array([dm.DEFAULTS[n] for n in names])))
    diff = max(abs(float(v) - float(runs.loc[0, f"search_{n}"])) for n, v in zip(names, again["x"], strict=False))
    out = dict(start=0, max_abs_difference_in_constants=diff, cost_recorded=float(runs.loc[0, "search_cost"]),
               cost_repeated=again["cost"], target=T, met=bool(diff <= T))
    json.dump(out, open(RESULTS/f"calibration_repeatability{SUFFIX}.json", "w"), indent=2)
    print(f"[calibrate] repeat of start 0: largest difference in any constant {diff:.2e} (target {T})")
    return out


_CV = {}


def _cv_fold(fold):
    """One fold, with everything that depends on data taken from its training
    loops only: the fit of every constant (CV_STARTS starts, the best on the
    training loops kept), the sensitivity of the training cost at that fit,
    the set of constants that sensitivity keeps, and the fit of that set. The
    left-out loops are used once, to score each candidate."""
    frames, statics, folds, base, thr = (_CV[k] for k in ("frames", "statics", "folds", "base", "threshold"))
    test = sorted(folds[fold].tolist())
    train = [i for i in range(len(frames)) if i not in test]

    def take(idx):
        return [frames[i] for i in idx], [statics[i] for i in idx]
    spec = {f[0]: f for f in FIT}
    rng = np.random.default_rng(SEED + 1 + fold)

    def best_fit(names):
        if not names:
            return np.array([])
        lo = np.array([spec[n][1] for n in names])
        hi = np.array([spec[n][2] for n in names])
        starts = [np.array([dm.DEFAULTS[n] for n in names])] + [lo + (hi - lo)*rng.random(len(names))
                                                                 for _ in range(CV_STARTS - 1)]
        sols = [least_squares(residuals, np.clip(x0, lo, hi), bounds=(lo, hi), method="trf", x_scale=hi - lo,
                              diff_step=DIFF_STEP, max_nfev=MAX_EVALUATIONS, args=(names, *take(train), base, True))
                for x0 in starts]
        return min(sols, key=lambda r: cost(r.fun)).x
    every = [f[0] for f in FIT]
    x_all = best_fit(every)
    j0 = cost(residuals(x_all, every, *take(train), base, True))
    worst = {}
    for i, n in enumerate(every):
        changes = []
        for pct in (-SENSITIVITY_PCT, SENSITIVITY_PCT):
            xi = x_all.copy()
            xi[i] *= 1.0 + pct/100.0
            changes.append(abs(100.0*(cost(residuals(xi, every, *take(train), base, True)) - j0)/j0))
        worst[n] = max(changes)
    sensitive = [n for n in every if worst[n] >= thr]
    rows = []
    for label, names, x in (("all", every, x_all), ("sensitive only", sensitive, best_fit(sensitive)), ("none", [], np.array([]))):
        rows.append(dict(constants_fitted=label, fold=fold, left_out_cost=cost(residuals(x, names, *take(test), base)),
                         constants=" ".join(names), **{n: round(float(v), 4) for n, v in zip(names, x, strict=False)}))
    return rows


def _one_standard_error_choice(by_fold, order):
    """Walk the candidates from the smallest set of fitted constants to the
    largest. A candidate replaces the one kept so far only if its mean
    left-out cost over the folds is lower by more than one standard error of
    the paired fold differences. Returns (chosen, the record of each step)."""
    chosen, steps = order[0], []
    for cand in order[1:]:
        d = (by_fold[chosen] - by_fold[cand]).values          # improvement of the larger set, fold by fold
        mean, se = float(d.mean()), float(d.std(ddof=1)/np.sqrt(len(d)))
        take = bool(mean > se)
        steps.append(dict(smaller=chosen, larger=cand, mean_improvement_per_fold=round(mean, 6),
                          standard_error=round(se, 6), larger_kept=take))
        if take:
            chosen = cand
    return chosen, steps


def select():
    """Decide which constants to fit by cross-validation inside the
    calibration set, and write the final constants.

    Candidates: every constant; only those whose +/- 10 % change moves the
    cost by at least the threshold in targets.json; none (literature values).
    The calibration loops are split into CV_FOLDS folds with SEED. For each
    fold the candidates are fitted on the other folds, the "sensitive" set
    being decided from those training loops alone, and scored on the fold
    left out (_cv_fold). The candidates are taken in order of size and a
    larger one replaces the next smaller only if its left-out cost is lower by
    more than one standard error of the fold-by-fold differences
    (_one_standard_error_choice); if the sensitive set is kept, the set used for the final fit is the
    one the sensitivity of the whole calibration set gives. No held-out loop
    is opened."""
    thr = json.load(open(DATA/"targets.json"))["calibration"]["drop_constant_if_cost_change_below_pct"]
    sens = pd.read_csv(RESULTS/f"sensitivity{SUFFIX}.csv")
    worst = sens[sens.constant != "(calibrated)"].groupby("constant").cost_change_pct.apply(lambda v: v.abs().max())
    first = json.load(open(CONSTANTS_FILE))
    first = first.get("first_stage_record", first)
    every = list(first["constants"])
    sensitive = [n for n in every if worst[n] >= thr]
    candidates = {"all": every, "sensitive only": sensitive, "none": []}
    m, frames, statics = calibration_frames()
    perm = np.random.default_rng(SEED).permutation(len(frames))
    _CV.update(frames=frames, statics=statics, base=base_constants(), threshold=thr,
               folds=[perm[i::CV_FOLDS] for i in range(CV_FOLDS)])
    with Pool(min(PROCESSES, CV_FOLDS)) as pool:
        cv = pd.DataFrame([row for rows in pool.map(_cv_fold, range(CV_FOLDS), chunksize=1) for row in rows])
    cv.to_csv(RESULTS/f"cross_validation{SUFFIX}.csv", index=False)
    by_fold = cv.pivot(index="fold", columns="constants_fitted", values="left_out_cost")
    totals = by_fold.sum()
    chosen, steps = _one_standard_error_choice(by_fold, ["none", "sensitive only", "all"])
    others = [c for c in candidates if c != chosen]
    print("[calibrate] left-out cost by candidate:\n" + totals.round(4).to_string() + f"\n[calibrate] kept: {chosen}; "
          + "; ".join(f"better than '{c}' in {int((by_fold[chosen] < by_fold[c]).sum())} of {CV_FOLDS} folds" for c in others))
    if candidates[chosen] == every:
        result, runs = dict(first), pd.read_csv(RESULTS/f"calibration_runs{SUFFIX}.csv")
    else:
        result, runs = fit(names=candidates[chosen], write=False)
    runs.to_csv(RESULTS/f"calibration_runs_final{SUFFIX}.csv", index=False)
    result["first_stage_record"] = {k: first[k] for k in first if k != "first_stage_record"}
    result["selection"] = dict(
        method=f"{CV_FOLDS}-fold cross-validation inside the calibration set", chosen=chosen,
        candidates={c: dict(constants=n, left_out_cost=round(float(totals[c]), 6)) for c, n in candidates.items()},
        starts_per_fold=CV_STARTS, rule=json.load(open(DATA/"targets.json"))["calibration"]["selection_rule"],
        steps=steps,
        folds_in_which_chosen_is_better={c: int((by_fold[chosen] < by_fold[c]).sum()) for c in others},
        sensitive_set_by_fold={int(r.fold): r.constants for r in cv[cv.constants_fitted == "sensitive only"].itertuples()},
        returned_to_literature_values=[n for n in every if n not in candidates[chosen]],
        largest_cost_change_pct_for_10pct={n: round(float(worst[n]), 3) for n in every})
    json.dump(result, open(CONSTANTS_FILE, "w"), indent=2)
    open(CONSTANTS_FILE, "a").write("\n")
    return result


if __name__ == "__main__":
    if "--attached" in sys.argv:
        attached_moment()
    elif "--select" in sys.argv:
        select()
    elif "--onset" in sys.argv:
        onset_level()
    elif "--curve" in sys.argv:
        curve_choice()
    elif "--repeat" in sys.argv:
        repeatability()
    elif "--fit" in sys.argv:
        fit()
    elif "--sensitivity" in sys.argv:
        sensitivity()
    else:
        sys.exit(__doc__)
