# -*- coding: utf-8 -*-
"""
unistall / validate.py
------------------
Score the calibrated load model and the reference model on every usable
measured frame, each at its own measured conditions and at the
reporting resolution, and compare the held-out results with the frozen targets.

Author: Akosa Samuel Onyejekwe (independent)

Outputs (results/):
  validation_table.csv      one row per frame and model, every metric
  validation_summary.csv    means with bootstrap 95 % intervals by group
  validation_targets.csv    each frozen target against the measured value
  comparison.csv            tabulated model against reference model, frame by
                            frame, with the paired difference and its interval
  comparison_same_constants.csv  the same with one set of dynamic constants
  load_difference_per_frame.csv  difference between the loads the two
                            separation laws predict with identical constants
  equivalence.csv           those differences against the equivalence margin
  validation_intervals_by_condition.csv  the held-out means with intervals
                            from resampling loops and from resampling
                            groups of loops that repeat one condition
  split_fixed_frames.csv    the frames whose set was fixed by name, and the
                            set the frame-number rule would have given each
  split_departures.csv      the held-out measures with those frames moved to
                            the rule's set, and with them left out
  validation_by_stall_depth.csv  the held-out errors grouped by stall depth
  validation_damping.csv    cycle damping: correlation, error and sign counts
  validation_moment_stall.csv    loops with a moment stall measured, predicted, both, neither
  validation_other_airfoil.csv   the AMES-01 frames, scored with that
                            aerofoil's own static sweep and the NACA 0012
                            dynamic constants
  validation_static_mach.csv     the static model against held-out
                            quasi-steady sweeps at Mach numbers between the
                            two stations it was built from

The held-out frames are opened here and nowhere upstream of here.

GROUPS
  primary        held-out, inside the static Mach range, peak
                 incidence inside the static lift range (<= 25 deg)
                 -- the group the frozen targets apply to
  primary_le_20  the same, peaking at or below 20 deg (static moment and drag
                 measured over the whole excursion)
  primary_20_25  the same, peaking between 20 and 25 deg
  outside_mach   held-out frames below Mach 0.20 (reported, no claim made)
  calibration    the frames the constants were fitted on
"""
import json, time
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
from unistall import dsmodel as dm
from unistall import metrics as mt
from unistall import reference_lb as ref
from unistall.static_model import StaticModel, Station, _frame_sweep, _kirchhoff_f, ALPHA_MAX
from unistall.strokes import stroke_split
from unistall.paths import DATA, RESULTS

T = json.load(open(DATA/"targets.json"))
# degrees by which the peak incidence passes static stall: below the first the
# loop is attached or marginal, at or above the second it is in deep stall
DEPTH_LIGHT_DEG, DEPTH_DEEP_DEG = 2.0, 6.0
REF_CONSTANTS = RESULTS/"calibrated_constants_reference.json"


STALL_CONSTANTS = ("Tp", "Tf0", "Tv0", "Tvl", "eta")


def _reference_constants():
    return dm.load_constants(REF_CONSTANTS)


def _literature_constants():
    """The calibrated set with the stall constants put back to their
    literature values (the attached-flow moment factor is kept)."""
    return {**dm.load_constants(), **{n: dm.DEFAULTS[n] for n in STALL_CONSTANTS}}


MACH_CALIBRATED = 0.285      # every calibration loop is at or above this Mach number
NEAR = dict(M=0.012, k=0.006, alpha0_deg=0.35, amp_deg=0.35)   # a held-out loop this close to a calibration loop repeats it


def near_calibration(man):
    """Held-out loops whose conditions repeat those of a calibration loop
    (every difference inside NEAR), as a boolean Series indexed by frame."""
    cal = man[(man["set"] == "calibration") & man.within_static_mach_range]
    def repeats(r):
        return bool(np.any(np.all([np.abs(cal[c].values - r[c]) <= tol for c, tol in NEAR.items()], axis=0)))
    return man.apply(repeats, axis=1) & (man["set"] == "held_out")


def groups(df):
    ho = df["set"] == "held_out"
    prim = ho & df.in_mach_range
    return {"primary": prim, "primary_le_20": prim & ~df.beyond_static_moment_drag_range,
            "primary_20_25": prim & df.beyond_static_moment_drag_range,
            "primary_at_calibrated_mach": prim & (df.M >= MACH_CALIBRATED),
            "primary_below_calibrated_mach": prim & (df.M < MACH_CALIBRATED),
            "primary_not_repeating_calibration": prim & ~df.near_calibration,
            "outside_mach": ho & ~df.in_mach_range, "calibration": (df["set"] == "calibration") & df.in_mach_range}


def score_all():
    man = pd.read_csv(mt.MANIFEST).set_index("frame")
    models = {"tabulated": dm.frame_runner(dm.load_constants())}
    if REF_CONSTANTS.exists():
        models["reference"] = ref.frame_runner(_reference_constants())
    tabs, cpu = [], {}
    for name, run in models.items():
        t0 = time.process_time()
        df = mt.score_manifest(run)
        cpu[name] = (time.process_time() - t0)/len(df)            # per loop, march and scoring together
        df.insert(0, "model", name)
        tabs.append(df)
    tab = pd.concat(tabs, ignore_index=True)
    tab["in_mach_range"] = tab.frame.map(man.within_static_mach_range).astype(bool)
    tab["quality_flags"] = tab.frame.map(man["flags"]).fillna("")
    tab["near_calibration"] = tab.frame.map(near_calibration(man)).astype(bool)
    return tab, cpu


def summaries(tab):
    out = []
    for name, df in tab.groupby("model", sort=False):
        s = mt.summarise(df.reset_index(drop=True), groups(df.reset_index(drop=True)))
        s.insert(0, "model", name)
        out.append(s)
    return pd.concat(out, ignore_index=True)


def targets_table(summ):
    """Each target fixed in advance against the value measured on the primary
    held-out group."""
    W = T["held_out_validation"]

    def val(metric):
        return float(summ[(summ.model == "tabulated") & (summ.group == "primary")
                          & (summ.metric == metric)].value.iloc[0])
    rows = []
    for metric in ("mean_nRMS_CL", "mean_nRMS_CM", "mean_nRMS_CD",
                   "mean_abs_dalpha_CLmax_deg", "mean_abs_dalpha_Mstall_deg"):
        new = val(metric)
        rows.append(dict(measure=metric, target=W[metric + "_max"], kind="at most", measured=round(new, 4),
                         met=bool(new <= W[metric + "_max"])))
    new = val("mean_Xi_sign_agree")
    rows.append(dict(measure="damping_sign_agreement", target=W["damping_sign_agreement_min"],
                     kind="at least", measured=round(new, 4), met=bool(new >= W["damping_sign_agreement_min"])))
    return pd.DataFrame(rows)


def _paired(a, b, cpu=None):
    """Tabulated (a) minus reference (b), frame by frame, on the primary
    held-out group. Both tables are indexed by frame."""
    g = groups(a.reset_index())
    prim = a.index[g["primary"].values]
    rng = np.random.default_rng(mt.BOOT_SEED)
    cols = ["nRMS_CL", "nRMS_CM", "nRMS_CD", "CLmax_err_pct", "dalpha_CLmax_deg", "dalpha_Mstall_deg"]
    per = pd.DataFrame({"frame": prim})
    rows = []
    for c in cols:
        xa = a.loc[prim, c].astype(float).abs().values
        xb = b.loc[prim, c].astype(float).abs().values
        per[f"{c}_tabulated"] = xa
        per[f"{c}_reference"] = xb
        both = np.isfinite(xa) & np.isfinite(xb)             # the loops on which both models have the measure
        xa, xb = xa[both], xb[both]
        d = xa - xb
        m = d[rng.integers(0, len(d), size=(mt.BOOT_N, len(d)))].mean(axis=1)
        lo, hi = np.percentile(m, [2.5, 97.5])
        rows.append(dict(metric="mean |" + c + "|", tabulated=round(float(np.nanmean(xa)), 4),
                         reference=round(float(np.nanmean(xb)), 4),
                         paired_difference=round(float(d.mean()), 4), ci95_lo=round(float(lo), 4),
                         ci95_hi=round(float(hi), 4), n_frames=len(d),
                         verdict=("tabulated better" if hi < 0 else "reference better" if lo > 0
                                  else "no difference at 95 %")))
    sa = a.loc[prim, "Xi_sign_agree"].mean()
    sb = b.loc[prim, "Xi_sign_agree"].mean()
    rows.append(dict(metric="damping sign agreement", tabulated=round(float(sa), 4), reference=round(float(sb), 4),
                     paired_difference=round(float(sa - sb), 4), n_frames=len(prim), verdict=""))
    if cpu is not None:
        rows.append(dict(metric="CPU seconds per loop (march and scoring)", tabulated=round(cpu["tabulated"], 4),
                         reference=round(cpu["reference"], 4),
                         paired_difference=round(cpu["tabulated"] - cpu["reference"], 4),
                         verdict="same machine, same run"))
    return pd.DataFrame(rows), per


def comparison(tab, cpu):
    """Each model with the constants its own calibration gave it."""
    return _paired(tab[tab.model == "tabulated"].set_index("frame"),
                   tab[tab.model == "reference"].set_index("frame"), cpu)


def comparison_same_constants(tab):
    """Both models with one set of dynamic constants, so that the only
    difference left is the separation law. Done three times: with the
    literature values, with the constants fitted for the tabulated model and
    with those fitted for the reference model."""
    flags = tab[tab.model == "tabulated"].set_index("frame")[["in_mach_range", "near_calibration"]]
    out = []
    for label, consts in _constant_sets().items():
        pair = []
        for runner in (dm.frame_runner(consts), ref.frame_runner(consts)):
            df = mt.score_manifest(runner)
            for col in flags:
                df[col] = df.frame.map(flags[col]).astype(bool)
            pair.append(df.set_index("frame"))
        res = _paired(*pair)[0]
        res.insert(0, "constants", label)
        out.append(res)
    return pd.concat(out, ignore_index=True)


EQUIVALENCE = T["equivalence_of_separation_laws"]


def _constant_sets():
    return {"literature values": _literature_constants(), "fitted for the tabulated model": dm.load_constants(),
            "fitted for the reference model": _reference_constants()}


def load_differences(tab):
    """Difference between the loads the two separation laws predict with
    identical constants, loop by loop on the primary held-out group: RMS and
    largest absolute difference of each coefficient over the cycle, raw and
    divided by the measured range of that coefficient (the normalisation of
    the loop errors). Returns one row per set of constants and loop."""
    t = tab[tab.model == "tabulated"].reset_index(drop=True)
    prim = set(t.frame[groups(t)["primary"].values])
    rows = []
    for label, consts in _constant_sets().items():
        run_tab, run_ref = dm.frame_runner(consts), ref.frame_runner(consts)
        for r in mt.manifest_frames().itertuples():
            if r.frame not in prim:
                continue
            fr = mt.load_frame(mt.FRAME_CACHE/f"{r.frame}.mat")
            a, b = run_tab(fr), run_ref(fr)
            if not np.allclose(a["alpha_deg"], b["alpha_deg"], atol=1e-12):
                raise RuntimeError(f"{r.frame}: the two models were not marched on the same incidences")
            row = dict(constants=label, frame=r.frame, M=round(fr["M"], 3), k=round(fr["k"], 3),
                       alpha0_deg=round(fr["a0"], 2), amp_deg=round(fr["da"], 2))
            for q in ("CL", "CM", "CD"):
                d = np.asarray(a[q], float) - np.asarray(b[q], float)
                rng = float(np.ptp(fr[q.lower()]))
                row[f"rms_d{q}"] = float(np.sqrt(np.mean(d**2)))
                row[f"max_abs_d{q}"] = float(np.abs(d).max())
                row[f"mean_d{q}"] = float(d.mean())
                row[f"rms_d{q}_over_range"] = row[f"rms_d{q}"]/rng
                row[f"max_abs_d{q}_over_range"] = row[f"max_abs_d{q}"]/rng
            rows.append(row)
    return pd.DataFrame(rows)


def equivalence_test(diff, same):
    """The equivalence margin of data/targets.json against (a) the mean
    normalised RMS load difference, by the upper end of its bootstrap 95 %
    interval, and (b) the paired difference in loop error, by its whole 95 %
    interval. One row per set of constants and coefficient."""
    rng = np.random.default_rng(mt.BOOT_SEED)
    rows = []
    for label, d in diff.groupby("constants", sort=False):
        for q in ("CL", "CM", "CD"):
            margin = EQUIVALENCE[f"margin_{q}"]
            mean, lo, hi, n = mt._boot(d[f"rms_d{q}_over_range"].values, rng)
            s = same[(same.constants == label) & (same.metric == f"mean |nRMS_{q}|")].iloc[0]
            worst = d.loc[d[f"max_abs_d{q}_over_range"].idxmax()]
            rows.append(dict(
                constants=label, coefficient=q, n_frames=n, margin=margin,
                mean_rms_difference_over_range=round(mean, 5), ci95_lo=round(lo, 5), ci95_hi=round(hi, 5),
                loads_within_margin=bool(hi < margin),
                mean_rms_difference_raw=round(float(d[f"rms_d{q}"].mean()), 5),
                mean_signed_difference_raw=round(float(d[f"mean_d{q}"].mean()), 5),
                worst_loop_rms_over_range=round(float(d[f"rms_d{q}_over_range"].max()), 5),
                worst_loop_max_abs_over_range=round(float(worst[f"max_abs_d{q}_over_range"]), 5),
                worst_loop_max_abs_raw=round(float(d[f"max_abs_d{q}"].max()), 5),
                worst_loop_frame=worst.frame,
                loop_error_difference=s.paired_difference, loop_error_ci95_lo=s.ci95_lo,
                loop_error_ci95_hi=s.ci95_hi,
                loop_error_within_margin=bool(-margin < s.ci95_lo and s.ci95_hi < margin),
                loop_error_bias=s.verdict))
    out = pd.DataFrame(rows)
    out["equivalent"] = out.loads_within_margin & out.loop_error_within_margin
    return out


SPLIT_METRICS = ("mean_nRMS_CL", "mean_nRMS_CM", "mean_nRMS_CD", "mean_abs_dalpha_CLmax_deg",
                 "mean_abs_dalpha_Mstall_deg", "mean_Xi_sign_agree")


def split_departures(tab):
    """The frames whose set was fixed by name in data/split.json, which of
    them the frame-number rule would have assigned differently, and the
    headline held-out measures three ways: with the split as fixed, with
    those frames moved to the set the rule gives (constants unchanged, so a
    frame moved out of calibration has been fitted), and with them left out.
    Returns (frames, measures)."""
    split = json.load(open(DATA/"split.json"))
    rule = split["rule_for_frames_added_later"]["naca0012_calibration_if"]
    if rule != "frame_number % 3 == 0":
        raise RuntimeError("the rule in data/split.json is not the one this check implements")
    man = pd.read_csv(mt.MANIFEST).set_index("frame")
    rows = []
    for name, lst in split["frames"].items():
        for e in lst:
            if e["airfoil"] != "NACA 0012":
                continue
            by_rule = "calibration" if int(e["frame"].split("_")[1]) % 3 == 0 else "held_out"
            rows.append(dict(frame=e["frame"], set_as_fixed=name, set_by_rule=by_rule,
                             departs_from_rule=bool(by_rule != name), M=man.M[e["frame"]], k=man.k[e["frame"]],
                             alpha0_deg=man.alpha0_deg[e["frame"]], amp_deg=man.amp_deg[e["frame"]]))
    frames = pd.DataFrame(rows)
    moved = set(frames.frame[frames.departs_from_rule])
    t = tab[tab.model == "tabulated"].reset_index(drop=True)
    by_rule = t["set"].where(~t.frame.isin(moved), t["set"].map({"calibration": "held_out", "held_out": "calibration"}))
    masks = {"split as fixed": (t["set"] == "held_out") & t.in_mach_range,
             "departing frames moved to the set the rule gives": (by_rule == "held_out") & t.in_mach_range,
             "departing frames left out": (t["set"] == "held_out") & t.in_mach_range & ~t.frame.isin(moved)}
    s = mt.summarise(t, masks)
    measures = s[s.metric.isin(SPLIT_METRICS)].reset_index(drop=True)
    return frames, measures


def clustered_intervals(tab):
    """The headline held-out means with two bootstrap intervals: resampling
    the loops one by one, as everywhere else, and resampling groups of loops
    that repeat one condition (every difference inside NEAR), which does not
    count a repeated condition as independent evidence. One row per measure."""
    t = tab[tab.model == "tabulated"].reset_index(drop=True)
    p = t[groups(t)["primary"].values].reset_index(drop=True)
    cluster = -np.ones(len(p), int)
    for i in range(len(p)):
        if cluster[i] >= 0:
            continue
        same = np.all([np.abs(p[c].values - p[c].values[i]) <= tol for c, tol in NEAR.items()], axis=0)
        cluster[same & (cluster < 0)] = i
    ids = np.unique(cluster)
    rng = np.random.default_rng(mt.BOOT_SEED)
    rows = []
    for col, how in (("nRMS_CL", "mean"), ("nRMS_CM", "mean"), ("nRMS_CD", "mean"), ("dalpha_CLmax_deg", "abs"),
                     ("dalpha_Mstall_deg", "abs"), ("Xi_sign_agree", "mean")):
        v = p[col].astype(float).abs().values if how == "abs" else p[col].astype(float).values
        ok = np.isfinite(v)
        mean, lo, hi, n = mt._boot(v, rng)
        means = []
        for _ in range(mt.BOOT_N):
            take = rng.choice(ids, size=len(ids))
            pick = np.concatenate([np.where((cluster == c) & ok)[0] for c in take])
            means.append(v[pick].mean() if len(pick) else np.nan)
        clo, chi = np.nanpercentile(means, [2.5, 97.5])
        rows.append(dict(measure=("mean_abs_" if how == "abs" else "mean_") + col, n_loops=n, n_conditions=len(ids),
                         value=round(mean, 4), ci95_lo_by_loop=round(lo, 4), ci95_hi_by_loop=round(hi, 4),
                         ci95_lo_by_condition=round(float(clo), 4), ci95_hi_by_condition=round(float(chi), 4)))
    return pd.DataFrame(rows)


def stall_depth_breakdown(tab):
    """Where the error sits: the primary held-out frames grouped by how far
    the peak incidence goes beyond static stall at the frame's Mach number.
    A diagnosis, not a target: the targets apply to the whole primary group."""
    t = tab[tab.model == "tabulated"].reset_index(drop=True)
    t = t[groups(t)["primary"].values].copy()
    over = np.array([r.alpha0_deg + r.amp_deg - StaticModel(r.M).alpha_stall_deg for r in t.itertuples()])
    t["beyond_static_stall_deg"] = over
    lo, hi = DEPTH_LIGHT_DEG, DEPTH_DEEP_DEG
    bins = [(f"attached or marginal (peak less than {lo:.0f} deg beyond static stall)", over < lo),
            (f"light stall ({lo:.0f} to {hi:.0f} deg beyond)", (over >= lo) & (over < hi)),
            (f"deep stall ({hi:.0f} deg or more beyond)", over >= hi)]
    rows = []
    for label, m in bins:
        d = t[m]
        if not len(d):
            continue
        rows.append(dict(group=label, n_frames=len(d), nRMS_CL=round(d.nRMS_CL.mean(), 4),
                         nRMS_CM=round(d.nRMS_CM.mean(), 4), nRMS_CD=round(d.nRMS_CD.mean(), 4),
                         damping_sign_agreement=round(float(d.Xi_sign_agree.mean()), 4),
                         mean_abs_measured_damping=round(float(d.Xi_exp.abs().mean()), 4)))
    return pd.DataFrame(rows)


def moment_stall_counts(tab):
    """How many primary held-out loops show a moment stall on the up-stroke in
    the measurement, in the model, in both or in neither. The incidence error
    is defined only where both show one, so the misses are counted here."""
    rows = []
    for model, t in tab.groupby("model", sort=False):
        t = t.reset_index(drop=True)
        t = t[groups(t)["primary"].values]
        e, m = t.Mstall_exp_found.astype(bool), t.Mstall_model_found.astype(bool)
        both = t[e & m]
        rows.append(dict(model=model, n_frames=len(t), in_both=int((e & m).sum()),
                         measured_only_missed_by_model=int((e & ~m).sum()),
                         model_only_not_measured=int((~e & m).sum()), in_neither=int((~e & ~m).sum()),
                         mean_abs_incidence_error_deg_where_both=round(float(both.dalpha_Mstall_deg.abs().mean()), 4),
                         mean_abs_error_in_minimum_moment=round(float((t.CMmin_model - t.CMmin_exp).abs().mean()), 4),
                         mean_measured_minimum_moment=round(float(t.CMmin_exp.mean()), 4)))
    return pd.DataFrame(rows)


def damping_uncertainty(frame, amp_deg):
    """Uncertainty of a measured cycle damping from the stated uncertainty u
    of the measured moment (uncertainty.U_CM and U_DIGIT_CM combined).
    Returns (sigma, bound): sigma if the error is independent from point to
    point; bound if it has one sign on the up-stroke and the other on the
    down-stroke, which is the worst case for a loop integral. A constant
    offset in C_M does not change the damping."""
    from unistall.uncertainty import U_CM, U_DIGIT_CM
    fr = mt.load_frame(mt.FRAME_CACHE/f"{frame}.mat")
    a = np.radians(np.append(fr["acm"], fr["acm"][0]))
    w = 0.5*(np.roll(a[:-1], -1) - np.roll(a[:-1], 1))            # weight of each point in the loop integral
    scale = np.hypot(U_CM, U_DIGIT_CM)/(np.pi*np.radians(amp_deg)**2)
    return float(scale*np.sqrt(np.sum(w**2))), float(scale*np.sum(np.abs(np.diff(a))))


def damping_table(tab):
    """What the model does and does not get right about cycle damping on the
    primary held-out loops: how well the predicted value follows the measured
    one, how often each is negative, and on how many loops the measured sign
    is itself certain."""
    t = tab[tab.model == "tabulated"].reset_index(drop=True)
    t = t[groups(t)["primary"].values].copy()
    t = t[t.moment_loop_closed.astype(bool)]              # a damping is measured only on a closed loop
    over = np.array([r.alpha0_deg + r.amp_deg - StaticModel(r.M).alpha_stall_deg for r in t.itertuples()])
    unc = np.array([damping_uncertainty(r.frame, r.amp_deg) for r in t.itertuples()])
    t["Xi_sigma"], t["Xi_bound"] = unc[:, 0], unc[:, 1]
    sets = [("all primary held-out loops with a closed moment loop", np.ones(len(t), bool)),
            ("measured damping positive", (t.Xi_exp > 0).values),
            ("measured damping negative", (t.Xi_exp < 0).values),
            (f"deep stall ({DEPTH_DEEP_DEG:.0f} deg or more beyond static stall)", over >= DEPTH_DEEP_DEG),
            ("measured sign certain if the moment error is random (beyond 2 sigma)", (t.Xi_exp.abs() > 2*t.Xi_sigma).values),
            ("measured sign certain even if the moment error follows the stroke", (t.Xi_exp.abs() > t.Xi_bound).values)]
    rows = []
    for label, m in sets:
        d = t[m]
        err = d.Xi_model - d.Xi_exp
        rows.append(dict(group=label, n_frames=len(d), measured_negative=int((d.Xi_exp < 0).sum()),
                         predicted_negative=int((d.Xi_model < 0).sum()),
                         sign_agreement=round(float(d.Xi_sign_agree.mean()), 4) if len(d) else np.nan,
                         sign_agreement_if_always_positive=round(float((d.Xi_exp > 0).mean()), 4) if len(d) else np.nan,
                         correlation=round(float(np.corrcoef(d.Xi_model, d.Xi_exp)[0, 1]), 3) if len(d) > 2 else np.nan,
                         rms_error=round(float(np.sqrt(np.mean(err**2))), 4) if len(d) else np.nan,
                         mean_error=round(float(err.mean()), 4) if len(d) else np.nan,
                         mean_measured=round(float(d.Xi_exp.mean()), 4) if len(d) else np.nan,
                         mean_predicted=round(float(d.Xi_model.mean()), 4) if len(d) else np.nan,
                         mean_uncertainty_random=round(float(d.Xi_sigma.mean()), 4) if len(d) else np.nan,
                         mean_uncertainty_stroke=round(float(d.Xi_bound.mean()), 4) if len(d) else np.nan))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
#  second aerofoil: its own static sweep, the NACA 0012 dynamic constants
# --------------------------------------------------------------------------- #
AMES_SWEEP = "frame_26020"       # quasi-steady sweep, M 0.30, -5 to 25 deg (NASA TM-84245 Table 12)


class SingleStation(StaticModel):
    """StaticModel for an aerofoil with static data at one Mach number."""
    def __init__(self, st, M, mach_tolerance=0.016):
        self.M, self.curve, self.w = float(M), "more", 0.0
        self._a = self._b = st
        self.in_range = bool(abs(M - st.M) <= mach_tolerance)
        self.slope_per_deg, self.alpha0_deg = st.slope_per_deg, st.alpha0_deg
        self.CN_alpha = st.slope_per_deg*180.0/np.pi
        self.alpha_stall_deg, self.CN1 = st.alpha_stall_deg, st.CN1
        self._grid = np.arange(0.0, ALPHA_MAX + 0.01, 0.02)
        self._ftab = {c: st.f(self._grid, c) for c in ("more", "less")}
        self._cmtab = {c: st.cm_static(self._grid, c) for c in ("more", "less")}


def ames_station():
    a, cn, cm, fr = _frame_sweep(AMES_SWEEP, "up")
    s = stroke_split(fr["acl"]) == "up"
    au, cl = fr["acl"][s], fr["cl"][s]
    lin = (au >= -5.1) & (au <= 8.0)
    slope, ic = np.polyfit(au[lin], cl[lin], 1)
    a0 = -ic/slope
    i = int(np.argmax(cn[a <= 20.0]))
    grid = np.arange(0.0, ALPHA_MAX + 1e-9, 0.5)
    cng = np.interp(grid, a, cn)
    f = _kirchhoff_f(grid, cng, slope, a0)
    cmg = np.interp(grid, a, cm)
    return Station(fr["M"], float(slope), float(a0), float(a[i]), float(cn[i]), grid,
                   {"more": f, "less": f}, {"more": cmg, "less": cmg}, f"quasi-steady sweep {AMES_SWEEP}")


def other_airfoil():
    st = ames_station()
    m = mt.manifest_frames(airfoil="AMES-01", sets=("other_airfoil",))
    consts = dm.load_constants()
    rows = []
    for r in m.itertuples():
        fr = mt.load_frame(mt.FRAME_CACHE/f"{r.frame}.mat")
        own = SingleStation(st, fr["M"])
        if not own.in_range:
            continue
        for label, sm in (("own static sweep", own), ("NACA 0012 static data", StaticModel(fr["M"]))):
            o = dm.solve(fr["a0"], fr["da"], fr["k"], fr["M"], consts=consts, static=sm)
            row = mt.score_frame(fr, o, airfoil="AMES-01", set_name="other_airfoil")
            row["static_inputs"] = label
            rows.append(row)
    df = pd.DataFrame(rows)
    return df, st


def static_mach_check():
    """Held-out quasi-steady sweeps at Mach numbers between the two stations:
    does the interpolated static model put static stall in the right place?"""
    m = mt.manifest_frames(sets=("held_out",), dynamic=False)
    rows = []
    for r in m.itertuples():
        if not r.within_static_mach_range:
            continue
        a, cn, cm, fr = _frame_sweep(r.frame, "up")
        sm = StaticModel(fr["M"])
        i = int(np.argmax(cn))
        g = np.arange(4.0, min(a.max(), ALPHA_MAX), 0.05)
        cnm = sm.CN_static(g)
        j = int(np.argmax(cnm))
        pre = a <= a[i]
        rms = float(np.sqrt(np.mean((sm.CN_static(a[pre]) - cn[pre])**2)))
        rows.append(dict(frame=r.frame, M=round(fr["M"], 3), CN_max_measured=round(float(cn[i]), 3),
                         CN_max_model=round(float(cnm[j]), 3),
                         CN_max_error_pct=round(100*(cnm[j] - cn[i])/cn[i], 2),
                         alpha_stall_measured=round(float(a[i]), 2), alpha_stall_model=round(float(g[j]), 2),
                         alpha_stall_error_deg=round(float(g[j] - a[i]), 2),
                         RMS_CN_up_to_stall=round(rms, 4)))
    return pd.DataFrame(rows)


if __name__ == "__main__":
    tab, cpu = score_all()
    tab.round(4).to_csv(RESULTS/"validation_table.csv", index=False)
    summ = summaries(tab)
    summ.to_csv(RESULTS/"validation_summary.csv", index=False)
    tg = targets_table(summ)
    ms = moment_stall_counts(tab).set_index("model").loc["tabulated"]
    tg["note"] = ""
    tg.loc[tg.measure == "mean_abs_dalpha_Mstall_deg", "note"] = (
        f"judged on the {int(ms.in_both)} loops where measurement and model both show a moment stall; the model "
        f"misses a measured one on {int(ms.measured_only_missed_by_model)} loops and predicts an unmeasured one on "
        f"{int(ms.model_only_not_measured)}")
    tg.to_csv(RESULTS/"validation_targets.csv", index=False)
    print(tg.to_string(index=False))
    if "reference" in set(tab.model):
        comp, per = comparison(tab, cpu)
        comp.to_csv(RESULTS/"comparison.csv", index=False)
        per.round(4).to_csv(RESULTS/"comparison_per_frame.csv", index=False)
        print(comp.to_string(index=False))
        same = comparison_same_constants(tab)
        same.to_csv(RESULTS/"comparison_same_constants.csv", index=False)
        print(same.to_string(index=False))
        diff = load_differences(tab)
        diff.round(6).to_csv(RESULTS/"load_difference_per_frame.csv", index=False)
        eq = equivalence_test(diff, same)
        eq.to_csv(RESULTS/"equivalence.csv", index=False)
        print(eq.to_string(index=False))
    sp_frames, sp_measures = split_departures(tab)
    sp_frames.to_csv(RESULTS/"split_fixed_frames.csv", index=False)
    sp_measures.to_csv(RESULTS/"split_departures.csv", index=False)
    print(sp_frames.to_string(index=False))
    print(sp_measures.to_string(index=False))
    clustered = clustered_intervals(tab)
    clustered.to_csv(RESULTS/"validation_intervals_by_condition.csv", index=False)
    print(clustered.to_string(index=False))
    depth = stall_depth_breakdown(tab)
    depth.to_csv(RESULTS/"validation_by_stall_depth.csv", index=False)
    print(depth.to_string(index=False))
    stall = moment_stall_counts(tab)
    stall.to_csv(RESULTS/"validation_moment_stall.csv", index=False)
    print(stall.to_string(index=False))
    damp = damping_table(tab)
    damp.to_csv(RESULTS/"validation_damping.csv", index=False)
    print(damp.to_string(index=False))
    oa, st = other_airfoil()
    oa.to_csv(RESULTS/"validation_other_airfoil.csv", index=False)
    if len(oa):
        print(f"[validate] AMES-01 (static stall {st.alpha_stall_deg:.1f} deg, C_N1 {st.CN1:.3f}, slope "
              f"{st.slope_per_deg:.4f}/deg), {oa.frame.nunique()} frames:")
        print(oa.groupby("static_inputs")[["nRMS_CL", "nRMS_CM", "nRMS_CD"]].mean().round(3).to_string())
        print(oa.groupby("static_inputs").apply(lambda d: d.CLmax_err_pct.abs().mean()).round(2).to_string())
    sm = static_mach_check()
    sm.to_csv(RESULTS/"validation_static_mach.csv", index=False)
    print(sm.to_string(index=False))
