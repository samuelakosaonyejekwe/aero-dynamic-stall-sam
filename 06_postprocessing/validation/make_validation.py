# -*- coding: utf-8 -*-
# Run from the repository root:  PYTHONPATH=. python3 06_postprocessing/validation/make_validation.py
"""
06_postprocessing / validation / make_validation.py
---------------------------------------------------
The comparison of the load model with measurement, set out beside the case
study: every file of this folder is written here from the files of this
solver (results/, data/, the measured frames in the cache), and the loops of
the figure are marched here with unistall.dsmodel at the calibrated
constants. Nothing is typed in.

Author: Akosa Samuel Onyejekwe (independent)

Standing: the held-out figures in this folder are the same numbers as those
of results/ (results/validation_table.csv, results/validation_targets.csv).
They are not a second assessment. A target the files give as not met is
reported as not met.

Writes, beside this file:
  validation_static.csv        the static inputs the model reads, and its
                               slow-sweep drag where results/limit_static_drag.csv
                               gives it, against the digitised static points
                               of data/static_naca0012_M030_*.csv
  fig_validation_static.png    the same as a figure
  validation_loops_summary.csv one row per scored loop of
                               results/validation_table.csv (tabulated model)
  validation_targets.csv       results/validation_targets.csv, passed through
  fig_validation_loops.png     measured against predicted loops of lift,
                               moment and drag for the loops `chosen_loops`
                               gives
  loop_<frame>_CL.csv, _CM.csv, _CD.csv   for those loops, the model read at
                               the measured incidences on the same stroke
  calibration_constants.csv    the constants of results/calibrated_constants.json
  PROVENANCE.md                where the measured frames come from
  README.md                    what each file is
  figure_record.csv            resolution and smallest lettering of the figures

Units: incidences in degrees, coefficients dimensionless, time constants in
semichords of travel.
"""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import project_meta as pm
from unistall import dsmodel, fetch_frames, metrics
from unistall.static_model import StaticModel
from unistall.strokes import stroke_split
from unistall.style import INK, INK_SOFT, PALETTE, apply_style

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
RESULTS, DATA = ROOT/"results", ROOT/"data"
VTABLE = RESULTS/"validation_table.csv"
MODEL = "tabulated"                      # the model of the case study in results/validation_table.csv
STATIC_MACH = 0.30                       # Mach number of the digitised static points
RANK_BY = "nRMS_CL"                      # the loops of the figure are chosen by this loop error
SUMMARY_COLUMNS = ["frame", "set", "M", "k", "alpha0_deg", "amp_deg", "nRMS_CL", "nRMS_CM", "nRMS_CD",
                   "CLmax_err_pct", "dalpha_CLmax_deg", "dalpha_Mstall_deg", "Xi_exp", "Xi_model", "Xi_sign_agree",
                   "in_mach_range"]
LOADS = (("CL", "cl", "acl", "lift coefficient  $C_L$"), ("CM", "cm", "acm", "moment coefficient about c/4  $C_M$"),
         ("CD", "cd", "acd", "drag coefficient  $C_D$"))
NOTE_PT = 9.0                            # lettering of the notes: the loops figure is printed reduced
RECORD = []                              # (file, dpi, smallest lettering) of each figure written


def save(fig: plt.Figure, name: str) -> None:
    """Write a figure, enter it in the figure record and close it."""
    RECORD.append((name, float(plt.rcParams["savefig.dpi"]), pm.smallest_lettering(fig)))
    fig.savefig(HERE/name)
    plt.close(fig)


# --------------------------------------------------------------------------- #
#  static
# --------------------------------------------------------------------------- #
def static_comparison() -> pd.DataFrame:
    """The digitised static points beside what the model reads at the same
    incidence. For lift the model's static input is a NORMAL force, set
    beside the measured lift as it is; for the moment it is the static
    moment; for drag the model has no static input, and the value is its
    slow-sweep drag from results/limit_static_drag.csv."""
    sm = StaticModel(STATIC_MACH)
    drag = pd.read_csv(RESULTS/"limit_static_drag.csv")
    frames = []
    for quantity, col, u_col, model_of, what in (
            ("CL", "Cl", "u_Cl", lambda a: sm.CN_static(a), "static normal force the model reads (CN, not CL)"),
            ("CM", "Cm_c4", "u_Cm", lambda a: sm.cm_static(a), "static moment the model reads"),
            ("CD", "Cd", "u_Cd", None, "drag of the model in a slow sweep (results/limit_static_drag.csv)")):
        d = pd.read_csv(DATA/f"static_naca0012_M030_{quantity}.csv")
        a = d["alpha_deg"].to_numpy(float)
        if model_of is None:
            model = np.interp(a, drag["alpha_deg"], drag["CD_model_slow_sweep"])
        else:
            model = np.asarray(model_of(a), float)
        frames.append(pd.DataFrame(dict(quantity=quantity, alpha_deg=a, measured=d[col].to_numpy(float),
                                        u_measured=d[u_col].to_numpy(float), model=np.round(model, 5),
                                        model_minus_measured=np.round(model - d[col].to_numpy(float), 5),
                                        model_value_is=what, measured_source=d["source"])))
    return pd.concat(frames, ignore_index=True)


def static_figure(table: pd.DataFrame) -> None:
    """The static comparison: measured points with their stated uncertainty,
    the model as a broken line."""
    fig, axs = plt.subplots(1, 3, figsize=(9.2, 4.6), constrained_layout=True)
    labels = {"CL": "$C_L$ measured;  $C_N$ read by the model", "CM": "$C_M$ about c/4", "CD": "$C_D$"}
    for ax, (q, sub) in zip(axs, table.groupby("quantity", sort=False), strict=True):
        ax.errorbar(sub["alpha_deg"], sub["measured"], yerr=sub["u_measured"], fmt="o", color=PALETTE[1], ms=4.5,
                    mfc="none", mew=1.2, elinewidth=0.9, capsize=2, label="measured (digitised)")
        m = sub.dropna(subset=["model"])
        ax.plot(m["alpha_deg"], m["model"], color=PALETTE[0], lw=1.9, ls="--", marker="s", ms=3.5,
                label="model")
        ax.set_xlabel("incidence  α  [deg]")
        ax.set_ylabel(labels[q])
    axs[1].legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2, frameon=True, fontsize=9)
    fig.suptitle("Static data at Mach %.2f against what the model reads" % STATIC_MACH, fontweight="bold",
                 color=INK, fontsize=12.5)
    fig.supxlabel("Left: the model reads a static normal force from quasi-steady sweeps of the experiment, which are "
                  "not corrected for the tunnel walls; the measured lift is from\nthe report's corrected static figure. "
                  "Middle: static moment. Right: the model has no static drag input; the line is its drag in a slow "
                  "sweep, which is\npressure drag where the measured points are total drag from a wake survey. "
                  "Sources of the points: data/static_naca0012_M030_*.csv.", fontsize=NOTE_PT, color=INK_SOFT)
    save(fig, "fig_validation_static.png")


# --------------------------------------------------------------------------- #
#  loops
# --------------------------------------------------------------------------- #
def scored_loops() -> pd.DataFrame:
    """The rows of results/validation_table.csv for the model of the case study."""
    df = pd.read_csv(VTABLE)
    return df[df["model"] == MODEL].reset_index(drop=True)


def chosen_loops(table: pd.DataFrame) -> list:
    """The loops of the figure, by a fixed rule, as (frame, why) pairs: of
    the held-out loops inside the static Mach range, the ones with the
    lowest, the median and the highest loop error in lift (the median being
    the middle one of the ranked list); and the calibration loop the case
    study sets beside Case A."""
    h = table[(table["set"] == "held_out") & table["in_mach_range"]].sort_values([RANK_BY, "frame"])
    picks = [(h.iloc[0]["frame"], "held-out loop with the lowest loop error in lift"),
             (h.iloc[len(h)//2]["frame"], "held-out loop with the median loop error in lift"),
             (h.iloc[-1]["frame"], "held-out loop with the highest loop error in lift")]
    frame = pm.CASES[pm.CASE_A]["measured_frame"]
    return picks + [(frame, "the CALIBRATION loop set beside Case A (not an independent test)")]


def on_measured_points(fr: dict, out: dict, key: str, a_key: str, y_key: str) -> pd.DataFrame:
    """The model read at the measured incidences of one load, each on the
    stroke the measured point belongs to."""
    a, y = fr[a_key], fr[y_key]
    branches = metrics._branches(out["alpha_deg"], out["alpha_dot"], out[key])
    strokes = stroke_split(a)
    model = np.array([np.interp(av, *branches[st]) for st, av in zip(strokes, a, strict=True)])
    return pd.DataFrame(dict(alpha_deg=np.round(a, 4), stroke=strokes, measured=np.round(y, 5),
                             model=np.round(model, 5), model_minus_measured=np.round(model - y, 5)))


def loops(table: pd.DataFrame) -> list:
    """March the chosen loops at their own measured conditions with the
    calibrated constants, write their per-point tables, and return
    (frame, why, measured frame, model solution, row of the table)."""
    run = dsmodel.frame_runner(dsmodel.load_constants())
    done = []
    for frame, why in chosen_loops(table):
        path = metrics.FRAME_CACHE/f"{frame}.mat"
        if not path.exists():
            raise FileNotFoundError(f"{path} missing: run python3 -m unistall.fetch_frames first")
        fr = metrics.load_frame(path)
        out = run(fr)
        for key, y_key, a_key, _ in LOADS:
            on_measured_points(fr, out, key, a_key, y_key).to_csv(HERE/f"loop_{frame}_{key}.csv", index=False)
        done.append((frame, why, fr, out, table[table["frame"] == frame].iloc[0]))
    return done


def loops_figure(done: list) -> None:
    """Measured against predicted loops: one row a loop, one column a load."""
    fig, axs = plt.subplots(len(done), 3, figsize=(9.2, 2.7*len(done) + 0.8), constrained_layout=True)
    for row, (frame, why, fr, out, score) in zip(axs, done, strict=True):
        for ax, (key, y_key, a_key, label) in zip(row, LOADS, strict=True):
            ax.plot(out["alpha_deg"], out[key], color=PALETTE[0], lw=1.9, label="predicted (load model)")
            ax.plot(fr[a_key], fr[y_key], "o", color=PALETTE[1], ms=3.8, mfc="none", mew=1.1, label="measured")
            ax.set_ylabel("%s\n(loop error %.3f)" % (label, score["nRMS_" + key]), fontsize=9.5)
            ax.tick_params(labelsize=9)
        row[1].set_title("%s: %s\nM = %.3f, k = %.3f, α = %.1f° ± %.1f°"
                         % (frame.replace("_", " "), why, fr["M"], fr["k"], fr["a0"], fr["da"]),
                         fontsize=9, color=INK, fontweight="normal")
    for ax in axs[-1]:
        ax.set_xlabel("incidence  α  [deg]")
    axs[-1][1].legend(loc="upper center", bbox_to_anchor=(0.5, -0.3), ncol=2, frameon=True, fontsize=9)
    fig.suptitle("Measured against predicted loops", fontweight="bold", color=INK, fontsize=12.5)
    fig.supxlabel("Loops chosen by rule from results/validation_table.csv: lowest, median and highest loop error in "
                  "lift of the held-out loops inside the static Mach range,\nand the calibration loop of Case A. Each "
                  "is marched at its own measured conditions with the calibrated constants. The loop error is the RMS "
                  "difference\nover the measured range, as in that file.", fontsize=NOTE_PT, color=INK_SOFT)
    save(fig, "fig_validation_loops.png")


# --------------------------------------------------------------------------- #
#  constants and notes
# --------------------------------------------------------------------------- #
def calibration_constants() -> pd.DataFrame:
    """The constants of results/calibrated_constants.json: fitted ones with
    their bounds and literature values, fixed ones with theirs."""
    cal = json.loads((RESULTS/"calibrated_constants.json").read_text(encoding="utf-8"))
    rows = []
    for status, group in (("fitted", cal["constants"]), ("fixed", cal["fixed"])):
        for name, value in group.items():
            lo, hi = cal.get("bounds", {}).get(name, (np.nan, np.nan))
            rows.append(dict(constant=name, value=value, status=status, literature_value=dsmodel.DEFAULTS.get(name),
                             lower_bound=lo, upper_bound=hi))
    return pd.DataFrame(rows)


def provenance() -> str:
    """The text of PROVENANCE.md."""
    n = len(list(metrics.FRAME_CACHE.glob("frame_*.mat")))
    return f"""# Provenance of the measured data

The measured loops are third-party data. They are fetched by `unistall/fetch_frames.py` from the repository
named in that script at the fixed commit `{fetch_frames.SOURCE_COMMIT}`, each file is checked against the
SHA-256 recorded in `data/frames_inventory.csv`, and they are kept in a local cache
(`{metrics.FRAME_CACHE.relative_to(ROOT)}`, {n} frame files at the time of writing). They are not redistributed in
this repository, and none is copied into this folder: the `loop_<frame>_*.csv` tables here hold the measured
points of {len(chosen_loops(scored_loops()))} loops only, beside the model's values.

The measurements themselves are those of the oscillating-aerofoil experiment reported in
{pm.experiment_citation()} (entries `mccroskey1982v1` and `mcalister1982v2` of `docs/references.bib`).

The static points are digitised or transcribed from that report; each row of `data/static_naca0012_M030_*.csv`
names the figure or table and the page it comes from.
"""


def readme(done: list) -> str:
    """The text of README.md."""
    chosen = "\n".join(f"- `{frame}`: {why}" for frame, why, *_ in done)
    return f"""**Standing: the held-out results here are the same numbers as those of `results/`, not a second assessment.**

# 06_postprocessing/validation

Every file is written by `make_validation.py` from the files of this solver. Run it from the repository root:
`PYTHONPATH=. python3 06_postprocessing/validation/make_validation.py`.

| File | What it is |
|---|---|
| `validation_static.csv`, `fig_validation_static.png` | the digitised static points at Mach {STATIC_MACH:.2f} beside what the model reads at the same incidence: its static normal force (not lift) and static moment, and for drag its slow-sweep value from `results/limit_static_drag.csv` |
| `validation_loops_summary.csv` | one row per scored loop of `results/validation_table.csv` for the {MODEL} model: conditions, loop errors, peak-lift error, stall-incidence errors, measured and predicted damping |
| `validation_targets.csv` | `results/validation_targets.csv`, passed through: each target with the measured value and whether it is met, as that file gives it |
| `fig_validation_loops.png` | measured against predicted loops for the loops listed below |
| `loop_<frame>_CL.csv`, `_CM.csv`, `_CD.csv` | for those loops, the model read at the measured incidences, each point on its own stroke |
| `calibration_constants.csv` | the constants of `results/calibrated_constants.json`, fitted and fixed, with literature values and bounds |
| `PROVENANCE.md` | where the measured frames come from |
| `figure_record.csv` | resolution and smallest lettering of the two figures |

The loops of the figure are chosen by a fixed rule (`chosen_loops`), not picked by eye:

{chosen}

A target that `validation_targets.csv` gives as not met is not met; nothing here changes that.
"""


def main() -> None:
    """Write every file of the folder."""
    apply_style()
    static = static_comparison()
    static.to_csv(HERE/"validation_static.csv", index=False)
    static_figure(static)
    table = scored_loops()
    table[SUMMARY_COLUMNS].to_csv(HERE/"validation_loops_summary.csv", index=False)
    pd.read_csv(RESULTS/"validation_targets.csv").to_csv(HERE/"validation_targets.csv", index=False)
    for old in HERE.glob("loop_*.csv"):
        old.unlink()
    done = loops(table)
    loops_figure(done)
    calibration_constants().to_csv(HERE/"calibration_constants.csv", index=False)
    (HERE/"PROVENANCE.md").write_text(provenance(), encoding="utf-8")
    (HERE/"README.md").write_text(readme(done), encoding="utf-8")
    pm.record_figures(HERE/"figure_record.csv", RECORD)
    print("[validation] %d scored loops summarised; loops drawn: %s"
          % (len(table), ", ".join(frame for frame, *_ in done)))


if __name__ == "__main__":
    main()
