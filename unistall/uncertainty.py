# -*- coding: utf-8 -*-
"""
unistall / uncertainty.py
---------------------
How much of the model-minus-measurement difference on the held-out frames can
be attributed to uncertainty in the inputs, and how much cannot.

Author: Akosa Samuel Onyejekwe (independent)

Four sources are considered, the first three propagated, one at a time, through the calibrated model on
the primary held-out group (inside Mach 0.20-0.30, peaking at or below 25 deg):

  CALIBRATION   every start of the calibration whose final cost is within
                COST_WINDOW of the best is an equally defensible set of
                constants; the model is re-run with each.
  STATIC DATA   the onset level C_N1 is moved by +/- the report's stated
                uncertainty in maximum lift (0.03), and the lift-curve slope
                by +/- its stated uncertainty (0.003 per degree)
                (NASA TM-84245 Vol. 1 Table 8, last line).
  MOMENT FACTOR the empirical factor on the unsteady moment is moved to the
                two ends of its 95 % bootstrap interval over calibration loops
                (results/attached_moment_factor.json); the stall constants
                are not refitted.
  FIT RANGE     the lift-curve slope and the zero-lift incidence are moved by
                half their spread over five choices of the incidence range
                the static line is fitted on (fit_range_spread).
  POST-STALL    the other of the two measured post-stall static states is
  CURVE         used. This is a choice of the model, made on the calibration
                loops, not an uncertainty of an input; it is reported beside
                the band and kept out of it.
  COMBINED      the sources above are varied one at a time; the combined
                half-width is the root sum of squares of the largest
                departure each produces, which assumes they are independent.
  MEASUREMENT   a model that reproduced the true loads exactly would still
                show a normalised RMS error of (uncertainty)/(measured range)
                against data carrying that uncertainty. The report states
                0.03 in C_L near maximum lift and 0.005 in C_M at M = 0.30;
                DIGITISATION adds what the secondary source's tracing of the
                printed loops adds, estimated from the agreement of its
                extreme values with the values printed in Volume 2.

Outputs: results/uncertainty_variants.csv (group means for every variant) and
results/uncertainty_summary.csv (band per metric, and the measurement floor).
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
from unistall import dsmodel as dm
from unistall import metrics as mt
from unistall.static_model import StaticModel
from unistall.paths import RESULTS

COST_WINDOW = 0.02            # starts within 2 % of the best cost
U_CN1, U_SLOPE = 0.03, 0.003  # NASA TM-84245 Vol. 1 Table 8, nominal uncertainties
U_CL, U_CM = 0.03, 0.005      # stated measurement uncertainty at M = 0.30
U_CD = 0.01                   # pressure drag: the report calls it markedly less certain than C_L or C_M; no figure is given, this one is assumed
U_DIGIT_CL, U_DIGIT_CM = 0.005, 0.005   # secondary-source tracing, from extreme-value agreement with Vol. 2
SOURCES = ("calibration", "static data", "fit range", "moment factor")   # varied one at a time and combined
METRICS = ["nRMS_CL", "nRMS_CM", "nRMS_CD", "CLmax_err_pct", "dalpha_CLmax_deg", "dalpha_Mstall_deg"]


def primary_frames():
    m = mt.manifest_frames(sets=("held_out",))
    return m[m.within_static_mach_range].reset_index(drop=True)


class _Shifted(StaticModel):
    """The static inputs with the stall level, the slope or the zero-lift
    incidence moved, or with the other post-stall curve. The separation table
    is kept as it is: it was derived with the unshifted slope and zero-lift
    incidence, so a shift here moves the attached loads and the onset level
    against a fixed table."""
    def __init__(self, M, dCN1=0.0, dslope=0.0, dalpha0=0.0, curve="more"):
        super().__init__(M, curve)
        self.CN1 += dCN1
        self.slope_per_deg += dslope
        self.CN_alpha = self.slope_per_deg*180.0/np.pi
        self.alpha0_deg += dalpha0


FIT_RANGES_DEG = ((-5.0, 8.0), (-5.0, 6.0), (-3.0, 8.0), (-5.0, 10.0), (0.0, 8.0))   # the first is the one used


def fit_range_spread():
    """Half the spread of the lift-curve slope [1/deg] and of the zero-lift
    incidence [deg] over FIT_RANGES_DEG, largest over the two static sweeps:
    how far those two inputs move with the choice of fit range."""
    from unistall import static_stations as st
    from unistall.strokes import stroke_split
    slope, zero = 0.0, 0.0
    for name in st.SWEEPS:
        fr = mt.load_frame(mt.FRAME_CACHE/f"{name}.mat")
        up = stroke_split(fr["acl"]) == "up"
        a, cl = fr["acl"][up], fr["cl"][up]
        fits = [np.polyfit(a[(a >= lo) & (a <= hi)], cl[(a >= lo) & (a <= hi)], 1) for lo, hi in FIT_RANGES_DEG]
        slopes = np.array([f[0] for f in fits])
        zeros = np.array([-f[1]/f[0] for f in fits])
        slope, zero = max(slope, 0.5*float(np.ptp(slopes))), max(zero, 0.5*float(np.ptp(zeros)))
    return slope, zero


def run_variant(frames, consts, dCN1=0.0, dslope=0.0, dalpha0=0.0, curve="more"):
    rows = []
    for fr in frames:
        o = dm.solve(fr["a0"], fr["da"], fr["k"], fr["M"], consts=consts,
                     static=_Shifted(fr["M"], dCN1, dslope, dalpha0, curve))
        rows.append(mt.score_frame(fr, o))
    df = pd.DataFrame(rows)
    out = {m: float(df[m].astype(float).abs().mean()) for m in METRICS}
    out["Xi_sign_agree"] = float(df.Xi_sign_agree.mean())
    return out, df


def main():
    m = primary_frames()
    frames = [mt.load_frame(mt.FRAME_CACHE/f"{f}.mat") for f in m.frame]
    base = dm.load_constants()
    runs = pd.read_csv(RESULTS/"calibration_runs_final.csv")
    names = list(json.load(open(dm.CONSTANTS_FILE))["constants"])
    near = runs[runs.cost <= (1.0 + COST_WINDOW)*runs.cost.min()]
    d_slope, d_zero = fit_range_spread()
    variants = [("calibrated", "baseline", base, {})]
    for r in near.itertuples():
        variants.append((f"calibration start {r.start}", "calibration", {**base, **{n: getattr(r, n) for n in names}}, {}))
    variants += [("C_N1 + 0.03", "static data", base, dict(dCN1=U_CN1)), ("C_N1 - 0.03", "static data", base, dict(dCN1=-U_CN1)),
                 ("slope + 0.003/deg", "static data", base, dict(dslope=U_SLOPE)),
                 ("slope - 0.003/deg", "static data", base, dict(dslope=-U_SLOPE))]
    variants += [(f"slope {s:+.4f}/deg (fit range)", "fit range", base, dict(dslope=s)) for s in (d_slope, -d_slope)]
    variants += [(f"zero-lift incidence {z:+.3f} deg (fit range)", "fit range", base, dict(dalpha0=z))
                 for z in (d_zero, -d_zero)]
    lo, hi = json.load(open(RESULTS/"attached_moment_factor.json"))["basis"]["interval_95"]
    variants += [(f"unsteady-moment factor {g}", "moment factor", {**base, "cm_unsteady": g}, {}) for g in (lo, hi)]
    variants += [("the other post-stall static curve (less separated)", "post-stall curve", base, dict(curve="less"))]
    rows = []
    base_df = None
    for label, source, c, shift in variants:
        out, df = run_variant(frames, c, **shift)
        if label == "calibrated":
            base_df = df
        rows.append(dict(variant=label, source=source, **{k: round(v, 4) for k, v in out.items()}))
    var = pd.DataFrame(rows)
    var.to_csv(RESULTS/"uncertainty_variants.csv", index=False)
    b = var.iloc[0]
    summ = []
    for metric in METRICS + ["Xi_sign_agree"]:
        row = dict(metric=metric, calibrated=b[metric])
        half = []
        for source in SOURCES:
            v = var[var.source.isin([source, "baseline"])][metric]
            row[f"{source}_min"] = round(float(v.min()), 4)
            row[f"{source}_max"] = round(float(v.max()), 4)
            half.append(max(abs(float(v.max()) - b[metric]), abs(float(v.min()) - b[metric])))
        allv = var[var.source != "post-stall curve"][metric]
        row["all_inputs_min"], row["all_inputs_max"] = round(float(allv.min()), 4), round(float(allv.max()), 4)
        row["combined_half_width"] = round(float(np.sqrt(np.sum(np.square(half)))), 4)
        other = var[var.source == "post-stall curve"][metric]
        row["other_post_stall_curve"] = round(float(other.iloc[0]), 4)
        summ.append(row)
    summ = pd.DataFrame(summ)
    # measurement floor: normalised RMS a perfect model would show
    floor = {}
    for q, u in (("CL", np.hypot(U_CL, U_DIGIT_CL)), ("CM", np.hypot(U_CM, U_DIGIT_CM)), ("CD", U_CD)):
        floor[f"nRMS_{q}"] = float((u/base_df[f"range_{q}"]).mean())
    summ["measurement_floor"] = summ.metric.map(floor).round(4)
    summ["calibrated_minus_floor"] = (summ.calibrated - summ.measurement_floor).round(4)
    summ.to_csv(RESULTS/"uncertainty_summary.csv", index=False)
    pd.set_option("display.width", 220)
    print(f"[uncertainty] {len(frames)} held-out frames; {len(near)} calibration starts within "
          f"{100*COST_WINDOW:.0f} % of the best cost")
    print(summ.to_string(index=False))


if __name__ == "__main__":
    main()
