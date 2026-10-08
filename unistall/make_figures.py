# -*- coding: utf-8 -*-
"""
unistall / make_figures.py
----------------------
Every figure in results/figures/, drawn from the result files and the model.

Author: Akosa Samuel Onyejekwe (independent)

  fig01_static_inputs       static normal force, separation point and moment
                            at the two Mach stations, table and fitted law
  fig02_attached_flow       attached-flow loads against Theodorsen's solution
  fig03_attached_moment     moment loops below stall: inviscid theory and the
                            measured unsteady-moment factor
  fig04_calibration_loops   three calibration loops (lowest, median and
                            highest cost), each with its RMS errors
  fig05_held_out_loops      three held-out loops chosen by rule (lowest,
                            median and highest lift error), both separation
                            laws, at print size
  fig06_effect_of_constants one held-out deep-stall loop with the calibrated
                            constants, the literature stall constants, and
                            the unsteady moment terms at full strength
  fig07_dynamic_and_static  predicted dynamic loops against the static curves
                            the model reads and the measured static points,
                            deep stall and light stall, with direction arrows
  fig08_design_space        maps of peak lift, minimum moment and peak drag
                            over mean incidence and reduced frequency
  fig09_damping_map         cycle damping over the same plane as contour
                            lines, with the measured loops marked by the sign
                            of their damping and by their set
  fig10_mach_trends         maps of the same peaks over reduced frequency and
                            Mach number, inside the range the measured loops
                            span, with the loops measured at each condition
  fig11_error_by_condition  loop error against stall depth and frequency
  fig12_separation_laws     tabulated against fitted law: paired differences
  fig13_cycle_damping       predicted against measured cycle damping
  fig14_calibration         cost of every start, sensitivity, cross-validation
                            and the cost surface
  fig15_model_states        the model's internal states over one cycle
  fig16_accuracy_summary    held-out errors, their 95 % intervals, the
                            measurement floor and the targets
  fig17_convergence         loads and loop errors against step size
  fig18_state_space         the state-space form against the indicial march

Two of them need model runs that no other script makes; those are written to
results/design_space.csv and results/cost_surface.csv. Each of the two tables
carries, in the column `inputs_sha256`, a signature of the files and the grid
it was computed from, and is computed again when that signature is not the
current one (or when --refresh is given).

Every figure is written by `_save`, which stores in the PNG (text chunk
INPUTS_KEY) the SHA-256 of each result and data file the figure read.
`stale_inputs` compares them with the files as they are now, so a figure
older than a file it was drawn from can be found (tests/test_figures.py).

Loops shown in figs 4 to 6 are chosen by rule from the model's own scores
(cost or lift error), not by eye; the rule is printed on each figure.

Usage:  python3 -m unistall.make_figures            all figures
        python3 -m unistall.make_figures fig07      the ones whose name starts so
        python3 -m unistall.make_figures --refresh  compute the two tables again
"""
import sys, json, os, hashlib
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.patches import Patch
from PIL import Image
from PIL.PngImagePlugin import PngInfo

from unistall import dsmodel as dm
from unistall import metrics as mt
from unistall import reference_lb as ref
from unistall import calibrate as cal
from unistall.static_model import StaticModel, stations, F_MIN
from unistall.style import apply_style, PALETTE, INK, INK_SOFT, CMAP_FIELD, HATCHES
from unistall.paths import DATA, RESULTS, FIGURES, ROOT

BLUE, RED, GREEN, ORANGE, PURPLE = PALETTE[:5]
MEASURED = dict(marker="o", ms=3.2, mfc=RED, mec=RED, ls="none")       # measured points, everywhere
MAP_MACH, MAP_AMP = 0.30, 5.0                       # condition of the design-space maps
MAP_MACH_BAND, MAP_AMP_BAND = 0.02, 0.3             # a measured loop this close in Mach and amplitude [deg] belongs to it
MAP_MEAN_SPAN, MAP_MEAN_STEP = (4.0, 18.0), 0.5     # deg; widened to the measured loops if they reach further
MAP_K_SPAN, MAP_K_STEP = (0.01, 0.20), 0.01
SURFACE_TVL_SPAN, SURFACE_TVL_N = (5.0, 13.0), 9    # semichords; the calibrated value is added to each grid
SURFACE_TP_SPAN, SURFACE_TP_N = (0.5, 3.5), 7
TREND_MATCH_DEG = 0.6                               # a loop this close in mean and amplitude is "at" a trend condition
STATE_CASE = dict(alpha_mean_deg=15.0, alpha_amp_deg=10.0, k=0.10, M=0.30)
STATIC_CASES = (("deep stall", dict(alpha_mean_deg=10.0, alpha_amp_deg=10.0, k=0.10, M=0.30)),
                ("light stall", dict(alpha_mean_deg=10.0, alpha_amp_deg=5.0, k=0.10, M=0.30)))
DEPTH_LIGHT, DEPTH_DEEP = 2.0, 6.0
QUANTITIES = (("cl", "CL", "$C_L$"), ("cm", "CM", "$C_{M,c/4}$"), ("cd", "CD", "$C_D$"))
# one load, one colour, one marker, one line style and one hatch, wherever the three appear together
SERIES = (("CL", "lift", BLUE, "o", "-", HATCHES[0]), ("CM", "moment", RED, "s", "--", HATCHES[1]),
          ("CD", "drag", GREEN, "^", ":", HATCHES[2]))
# the set a measured loop belongs to: filled for calibration, open for held out
SET_STYLE = {"calibration": dict(mfc=INK_SOFT, label="calibration loop"),
             "held_out": dict(mfc="white", label="held-out loop")}
LABEL_CLEAR, LABEL_CLEAR_LINE = 0.07, 0.03          # clear space round a contour label, as a fraction of the axes
INPUTS_KEY = "unistall-inputs"                      # PNG text chunk: {file: sha256} of what the figure read
SIGNATURE = "inputs_sha256"                         # column of the two tables this script owns
WORKERS = os.cpu_count() or 1
REFRESH = "--refresh" in sys.argv                   # compute the two tables this script owns again
# what unistall/static_model.py reads to build the static inputs of any model run
STATIC_FILES = (RESULTS/"static_polar_naca0012_M030.csv", RESULTS/"static_by_mach_naca0012.csv",
                DATA/"static_naca0012_M030_CL.csv", mt.FRAME_CACHE/"frame_13308.mat", mt.FRAME_CACHE/"frame_12300.mat")


def journal_style():
    """Serif type and thin lines, as in a printed journal figure; no black."""
    apply_style()
    plt.rcParams.update({"font.family": "serif", "font.serif": ["STIXGeneral", "DejaVu Serif"],
                         "mathtext.fontset": "stix", "font.size": 10, "axes.titlesize": 10,
                         "axes.titleweight": "normal", "axes.labelsize": 10.5, "legend.fontsize": 9.5,
                         "lines.linewidth": 1.3, "axes.linewidth": 0.8, "grid.alpha": 0.35,
                         "xtick.direction": "in", "ytick.direction": "in"})


# --------------------------------------------------------------------------- #
#  what a figure was drawn from
# --------------------------------------------------------------------------- #
_READ = {}                 # file -> sha256, for every file read since the last figure was saved


def file_hash(path: Path) -> str:
    """SHA-256 of the file at `path`, as hexadecimal text."""
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _name(path):
    """A file's name in the record: relative to the repository where it lies inside it."""
    p = Path(path).resolve()
    return p.relative_to(ROOT).as_posix() if p.is_relative_to(ROOT) else p.as_posix()


def _use(path):
    """Note that the figure being drawn reads `path`, with the file's hash at this moment."""
    _READ[_name(path)] = file_hash(path)
    return Path(path)


def figure_inputs(png: Path) -> dict:
    """{file: sha256} stored in the figure `png` when it was saved; empty if it carries no record."""
    with Image.open(png) as im:
        return json.loads(getattr(im, "text", {}).get(INPUTS_KEY, "{}"))


def stale_inputs(png: Path) -> list:
    """Files the figure `png` was drawn from that have changed since. A file
    that is absent now is not compared."""
    return sorted(f for f, h in figure_inputs(png).items() if (ROOT/f).exists() and file_hash(ROOT/f) != h)


def _save(fig, name):
    """Write one figure, with the record of the files it read; the only way a figure is written."""
    if not _READ:
        raise RuntimeError(f"{name}: no input file was recorded, so the figure cannot be dated against its inputs")
    info = PngInfo()
    info.add_text(INPUTS_KEY, json.dumps(_READ, sort_keys=True))
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURES/f"{name}.png", dpi=220, pil_kwargs={"pnginfo": info})
    plt.close(fig)
    _READ.clear()
    print(f"[figures] {name}.png")


def _csv(name):
    """A result table, by its file name in results/ or by its path."""
    return pd.read_csv(_use(RESULTS/name))


def _json(path):
    with open(_use(path)) as fh:
        return json.load(fh)


def _targets():
    return _json(DATA/"targets.json")


def _frame(name):
    return mt.load_frame(_use(mt.FRAME_CACHE/f"{name}.mat"))


def _static(M, law="tabulated"):
    """The static inputs at Mach number M with the tabulated or the fitted separation law."""
    for f in STATIC_FILES:
        _use(f)
    return StaticModel(M) if law == "tabulated" else ref.ExponentialStatic(M)


def _constants(law="tabulated"):
    """The calibrated constants of the model with the tabulated or the fitted separation law."""
    return dm.load_constants(_use(dm.CONSTANTS_FILE if law == "tabulated" else RESULTS/"calibrated_constants_reference.json"))


def _solve(alpha_mean_deg, alpha_amp_deg, k, M, consts=None, law="tabulated"):
    """One cycle of the model; angles in degrees. `consts=None` is the calibrated set of that law."""
    return dm.solve(alpha_mean_deg, alpha_amp_deg, k, M, consts=_constants(law) if consts is None else consts,
                    static=_static(M, law))


def _solve_frame(fr, consts=None, law="tabulated"):
    """The model at the measured conditions of one loaded frame."""
    return _solve(fr["a0"], fr["da"], fr["k"], fr["M"], consts, law)


def _manifest(sets=("calibration", "held_out")):
    """The measured loops the study uses that lie inside the static Mach range."""
    _use(mt.MANIFEST)
    m = mt.manifest_frames(sets=sets)
    return m[m.within_static_mach_range].reset_index(drop=True)


def _table(model="tabulated"):
    t = _csv("validation_table.csv")
    return t[t.model == model].reset_index(drop=True)


def _primary():
    t = _table()
    t = t[(t["set"] == "held_out") & t.in_mach_range].reset_index(drop=True)
    t["beyond_stall_deg"] = [r.alpha0_deg + r.amp_deg - _static(r.M).alpha_stall_deg for r in t.itertuples()]
    return t


# --------------------------------------------------------------------------- #
#  drawing helpers
# --------------------------------------------------------------------------- #
def _letter(i):
    return f"({chr(97 + i)})"


def _legend_below(fig, axs, ncol=3, extra=()):
    """One legend under the panels, from the labelled artists of every axes in
    `axs` (and the handles in `extra`), each label once. The room kept for it
    follows from the number of its rows, so it never reaches the axis labels."""
    uniq = {}
    for ax in np.atleast_1d(axs).ravel():
        h, lab = ax.get_legend_handles_labels()
        uniq.update({k: v for k, v in zip(lab, h, strict=True) if k not in uniq})
    uniq.update({h.get_label(): h for h in extra})
    size = plt.rcParams["legend.fontsize"]
    rows = int(np.ceil(len(uniq)/ncol))
    pad = (1.5*rows + 1.0)*size/72.0/fig.get_figheight()           # rows of text and a margin, in inches
    fig.get_layout_engine().set(rect=(0, pad, 1, 1 - pad))
    fig.legend(uniq.values(), uniq.keys(), loc="lower center", ncol=ncol, frameon=False)


def _condition(fr):
    return (f"M = {fr['M']:.2f}, k = {fr['k']:.3f}, " + rf"$\alpha$ = {fr['a0']:.0f}° ± {fr['da']:.0f}°")


def _arrows(ax, x, y, color, n=3):
    """Direction of travel round a loop, as n small arrowheads."""
    for i in np.linspace(0, len(x) - 2, n + 2).astype(int)[1:-1]:
        ax.annotate("", xy=(x[i + 1], y[i + 1]), xytext=(x[i], y[i]),
                    arrowprops=dict(arrowstyle="-|>", color=color, lw=0.0, mutation_scale=11))


def _bars(ax, x, height, width, color, hatch, label):
    """Bars told apart by hatching as well as by colour."""
    return ax.bar(x, height, width, facecolor=color, edgecolor="white", linewidth=0.0, hatch=hatch, label=label)


def _field(fig, ax, x, y, z, label, reverse=False):
    """A filled map of z over (x, y) on a scale that survives greyscale, with
    its level lines and a scale bar."""
    levels = ticker.MaxNLocator(nbins=9).tick_values(float(np.nanmin(z)), float(np.nanmax(z)))
    cs = ax.contourf(x, y, z, levels=levels, cmap=CMAP_FIELD + ("_r" if reverse else ""))
    ax.contour(x, y, z, levels=levels, colors=[INK_SOFT], linewidths=0.45, linestyles="-")
    fig.colorbar(cs, ax=ax, label=label)
    return cs


def _log_ticks(axis, lo, hi):
    """Label a logarithmic axis from lo to hi fully: 1, 2 and 5 in every decade
    where it spans few decades (a single decade tick would otherwise be the
    only label), and every decade where it spans many."""
    if np.log10(hi/lo) < 2.5:
        axis.set_major_locator(ticker.LogLocator(base=10.0, subs=(1.0, 2.0, 5.0), numticks=15))
        axis.set_major_formatter(ticker.FuncFormatter(lambda v, _pos: f"{v:g}"))
    else:
        axis.set_major_locator(ticker.LogLocator(base=10.0, numticks=15))
    axis.set_minor_formatter(ticker.NullFormatter())


def _label_points(cs, ax, marks):
    """Where to label each contour level: the point of the level's own line
    that is farthest from the markers at `marks` (rows of x, y), from the
    labels already placed and from the lines of the other levels, and clear
    of the frame. A level with no clear point is left unlabelled."""
    (x0, x1), (y0, y1) = ax.get_xlim(), ax.get_ylim()

    def unit(p: np.ndarray) -> np.ndarray:
        """Data coordinates as fractions of the axes."""
        return (np.asarray(p, float).reshape(-1, 2) - (x0, y0))/(x1 - x0, y1 - y0)

    def nearest(u: np.ndarray, others: np.ndarray) -> np.ndarray:
        """Distance from each row of u to the nearest row of others."""
        if not len(others):
            return np.full(len(u), np.inf)
        return np.sqrt(((u[:, None, :] - others[None, :, :])**2).sum(axis=2)).min(axis=1)
    lines = [unit(np.vstack(s)) if len(s) else np.empty((0, 2)) for s in cs.allsegs]
    taken, out = unit(marks), []
    for i, u in enumerate(lines):
        u = u[((u > LABEL_CLEAR) & (u < 1.0 - LABEL_CLEAR)).all(axis=1)]
        if not len(u):
            continue
        rest = np.vstack([np.empty((0, 2))] + [v for j, v in enumerate(lines) if j != i])
        room = np.minimum(nearest(u, taken)/LABEL_CLEAR, nearest(u, rest)/LABEL_CLEAR_LINE)
        best = int(np.argmax(room))
        if room[best] >= 1.0:
            out.append(tuple(u[best]*(x1 - x0, y1 - y0) + (x0, y0)))
            taken = np.vstack([taken, u[best]])
    return out


# --------------------------------------------------------------------------- #
#  static inputs and attached flow
# --------------------------------------------------------------------------- #
def fig01_static_inputs():
    g = np.arange(0.0, 25.01, 0.1)
    pol = _csv("static_polar_naca0012_M030.csv")
    measured_to = {stations()[-1].M: float(pol.alpha_deg[pol.CM_measured].max())}      # the polar is the upper station
    fig, axs = plt.subplots(1, 3, figsize=(11.5, 4.2), constrained_layout=True)
    for j, (st, col, mk) in enumerate(zip(stations(), (BLUE, RED), ("o", "s"), strict=False)):
        sm, ex = _static(st.M), _static(st.M, "fitted")
        tag = f"M = {st.M:.3f}"

        def sym(n: int, j: int = j, mk: str = mk) -> dict:
            """Sparse symbols of this station, staggered between its curves."""
            return dict(marker=mk, ms=3.6, mfc="white", mew=0.8, markevery=(8 + 14*n + 7*j, 42))
        axs[0].plot(g, sm.CN_static(g), "-", color=col, label=f"{tag}, tabulated", **sym(0))
        axs[0].plot(g, ex.CN_static(g), "--", color=col, lw=1.0, label=f"{tag}, fitted law", **sym(1))
        axs[1].plot(g, sm.f(g), "-", color=col, label=f"{tag}, tabulated", **sym(0))
        axs[1].plot(g, ex.f(g), "--", color=col, lw=1.0, label=f"{tag}, fitted law", **sym(1))
        axs[1].plot(g, sm.f(g, "less"), ":", color=col, lw=1.1, label=f"{tag}, less separated post-stall state", **sym(2))
        end = measured_to.get(st.M, g[-1])
        axs[2].plot(g[g <= end], sm.cm_static(g[g <= end]), "-", color=col, label=f"{tag}, tabulated", **sym(0))
        if end < g[-1]:
            axs[2].plot(g[g >= end], sm.cm_static(g[g >= end]), "-.", color=col, lw=1.0,
                        label=f"{tag}, not measured beyond {end:.0f}°: held at the last value")
    axs[1].axhline(F_MIN, color=INK_SOFT, lw=1.0, ls=(0, (1, 2.5)), label=f"floor on $f$ ({F_MIN:g})")
    for ax, yl, ttl in zip(axs, ("static $C_N$", "separation point $f$", "static $C_{M,c/4}$"),
                           ("(a) Static normal force", "(b) Separation point", "(c) Static moment"), strict=False):
        ax.set_xlabel(r"incidence $\alpha$, deg")
        ax.set_ylabel(yl)
        ax.set_title(ttl, loc="left")
    _legend_below(fig, axs, ncol=3)
    _save(fig, "fig01_static_inputs")


def fig02_attached_flow():
    a = _csv("attached_checks.csv")
    a = a[a.check == "Theodorsen"]
    T = _targets()["attached_flow"]
    cols = (("lift_amp_err_pct", "lift amplitude error, %", T["lift_and_moment_amplitude_error_max_pct"]),
            ("lift_phase_err_deg", "lift phase error, deg", T["lift_and_moment_phase_error_max_deg"]),
            ("moment_amp_err_pct", "moment amplitude error, %", T["lift_and_moment_amplitude_error_max_pct"]),
            ("moment_phase_err_deg", "moment phase error, deg", T["lift_and_moment_phase_error_max_deg"]))
    keep = [v for v in a.variant.unique() if "simplified" not in v]
    fig, axs = plt.subplots(1, 4, figsize=(13, 3.9), constrained_layout=True)
    for i, (ax, (c, yl, lim)) in enumerate(zip(axs, cols, strict=False)):
        ax.axhspan(-lim, lim, facecolor=GREEN, alpha=0.12, edgecolor=GREEN, hatch="...", lw=0.0, label="target band")
        for v, col, mk, ls in zip(keep, (BLUE, ORANGE, RED), ("o", "s", "^"), ("-", "--", ":"), strict=False):
            d = a[a.variant == v]
            ax.plot(d.k, d[c], ls, marker=mk, ms=4.5, color=col, label=v)
        ax.set_xlabel("reduced frequency $k$")
        ax.set_ylabel(yl)
        ax.set_title(f"{_letter(i)} {yl.split(',')[0].capitalize()}", loc="left")
    fig.suptitle(f"Pitching about the quarter chord against Theodorsen's solution, M = {T['theodorsen_M']:g}", fontsize=10)
    _legend_below(fig, axs[0], ncol=2)
    _save(fig, "fig02_attached_flow")


def fig03_attached_moment():
    rec = _json(cal.ATTACHED_FILE)
    g = rec["cm_unsteady"]
    fig, axs = plt.subplots(1, len(rec["frames"]), figsize=(4.3*len(rec["frames"]), 4.4), constrained_layout=True)
    for i, (ax, r) in enumerate(zip(np.atleast_1d(axs), rec["frames"], strict=False)):
        fr = _frame(r["frame"])
        ax.plot(fr["acm"], fr["cm"], label="measured", **MEASURED)
        for factor, ls, col, lab in ((1.0, "--", ORANGE, "unsteady moment terms at full strength"),
                                     (g, "-", BLUE, f"with the empirical factor {g:.2f}")):
            o = _solve_frame(fr, consts=dict(dm.DEFAULTS, cm_unsteady=factor))
            ax.plot(o["alpha_deg"], o["CM"], ls, color=col, label=lab)
        ax.set_xlabel(r"$\alpha$, deg")
        ax.set_ylabel("$C_{M,c/4}$")
        meas = ("measured loop not closed, so it has none" if r["damping_measured"] is None
                else f"measured {r['damping_measured']:.3f}")
        level = r["peak_CN_prime_over_CN1"]
        onset = "reached: the model sheds a vortex" if level >= 1.0 else "not reached"
        ax.set_title(f"{_letter(i)} frame {r['frame'][6:]}, {_condition(fr)}\n"
                     f"cycle damping: full strength {r['damping_without_factor']:.3f}, with factor "
                     f"{r['damping_with_factor']:.3f};\n{meas}\n"
                     f"stall onset {onset} (peak $C'_N$ = {level:.2f} $C_{{N1}}$)", loc="left", fontsize=8.5)
    fig.suptitle("Calibration loops whose peak incidence is below static stall. Both model curves use the literature "
                 "stall constants (dsmodel.DEFAULTS), as the fit of the factor does", fontsize=9.5)
    _legend_below(fig, np.atleast_1d(axs)[0])
    _save(fig, "fig03_attached_moment")


# --------------------------------------------------------------------------- #
#  measured and predicted loops
# --------------------------------------------------------------------------- #
def _loop_rows(axs, picks, runners, size=9):
    """One row per (label, frame name, scores): measured points and one curve
    per runner (law, line style, colour, legend label). Where there are two
    runners the largest difference between their curves is printed."""
    for i, (label, name, score) in enumerate(picks):
        fr = _frame(name)
        outs = [(_solve_frame(fr, law=law), ls, col, lab) for law, ls, col, lab in runners]
        for j, (ax, (q, Q, yl)) in enumerate(zip(axs[i], QUANTITIES, strict=False)):
            ax.plot(fr["a" + q], fr[q], label="measured", **MEASURED)
            for o, ls, col, lab in outs:
                ax.plot(o["alpha_deg"], o[Q], ls, color=col, label=lab)
            ax.set_xlabel(r"$\alpha$, deg")
            ax.set_ylabel(yl)
            note = f"{_letter(3*i + j)} loop error in {yl} = {score['nRMS_' + Q]:.3f}"
            if len(outs) == 2:
                note += f"\ncurves differ by at most {float(np.max(np.abs(outs[0][0][Q] - outs[1][0][Q]))):.3f}"
            head = f"{label}: frame {name[6:]}, {_condition(fr)}" if j == 0 else ""
            ax.set_title(f"{head}\n{note}", loc="left", fontsize=size)


def fig04_calibration_loops():
    c = _table()
    c = c[(c["set"] == "calibration") & c.in_mach_range].copy()
    c["cost"] = c.nRMS_CL**2 + c.nRMS_CM**2 + c.nRMS_CD**2
    c = c.sort_values("cost").reset_index(drop=True)
    picks = [(lab, c.frame[j], c.iloc[j]) for lab, j in
             (("calibration, lowest cost", 0), ("calibration, median cost", len(c)//2), ("calibration, highest cost", len(c) - 1))]
    fig, axs = plt.subplots(3, 3, figsize=(11, 9.6), constrained_layout=True)
    _loop_rows(axs, picks, [("tabulated", "-", BLUE, "model")])
    fig.suptitle(f"Chosen by rule among the {len(c)} calibration loops inside the static Mach range: lowest, median and "
                 "highest cost of the calibrated model", fontsize=10)
    _legend_below(fig, axs[0, 0])
    _save(fig, "fig04_calibration_loops")


def fig05_held_out_loops():
    p = _primary().sort_values("nRMS_CL").reset_index(drop=True)
    n = len(p)
    picks = [(f"held out, {lab} lift error", p.frame[j], p.iloc[j]) for lab, j in
             (("lowest", 0), ("median", n//2), ("highest", n - 1))]
    with plt.rc_context({"font.size": 7.5, "axes.labelsize": 8, "legend.fontsize": 7.5, "lines.linewidth": 1.0}):
        fig, axs = plt.subplots(len(picks), 3, figsize=(7.2, 7.9), constrained_layout=True)
        _loop_rows(axs, picks, [("tabulated", "-", BLUE, "model, tabulated separation"),
                                ("fitted", "--", ORANGE, "model, fitted separation law")], size=6.8)
        fig.suptitle(f"Chosen by rule among the {n} primary held-out loops: lowest, median and highest lift error "
                     "of the model with tabulated separation", fontsize=7.5)
        _legend_below(fig, axs[0, 0])
        _save(fig, "fig05_held_out_loops")


def fig06_effect_of_constants():
    p = _primary()
    deep = p[p.beyond_stall_deg >= DEPTH_DEEP].sort_values("nRMS_CL").reset_index(drop=True)
    name = deep.frame[len(deep)//2]                       # the deep-stall loop of median lift error
    fr = _frame(name)
    cals = _constants()
    outs = [(_solve_frame(fr, consts=c), ls, col, lab) for c, ls, col, lab in (
        (cals, "-", BLUE, "calibrated"),
        ({**cals, **{n: dm.DEFAULTS[n] for n in ("Tp", "Tf0", "Tv0", "Tvl", "eta")}}, "--", INK_SOFT,
         "literature stall constants"),
        (dict(cals, cm_unsteady=1.0), ":", ORANGE, "unsteady moment terms at full strength"))]
    fig, axs = plt.subplots(1, 3, figsize=(11.5, 4.1), constrained_layout=True)
    for i, (ax, (q, Q, yl)) in enumerate(zip(axs, QUANTITIES, strict=False)):
        ax.plot(fr["a" + q], fr[q], label="measured", **MEASURED)
        for o, ls, col, lab in outs:
            ax.plot(o["alpha_deg"], o[Q], ls, color=col, label=lab)
        peak = fr[q].min() if q == "cm" else fr[q].max()
        ax.axhline(peak, color=RED, lw=0.8, ls=(0, (6, 2, 1, 2)), label="measured peak (in $C_{M,c/4}$: the break)")
        same = float(np.max(np.abs(outs[2][0][Q] - outs[0][0][Q]))) < 5e-4
        ax.set_xlabel(r"$\alpha$, deg")
        ax.set_ylabel(yl)
        ax.set_title(f"{_letter(i)} {yl}" + (": full-strength curve lies on the calibrated one" if same else ""),
                     loc="left", fontsize=9)
    fig.suptitle(f"The deep-stall held-out loop of median lift error ({len(deep)} loops peak {DEPTH_DEEP:g}° or more "
                 f"beyond static stall): frame {name[6:]}, {_condition(fr)}", fontsize=10)
    _legend_below(fig, axs[0], ncol=3)
    _save(fig, "fig06_effect_of_constants")


def _static_points(M, lo, hi):
    """The measured static points between lo and hi degrees, where the case's
    Mach number is that of the digitised static polar: per load, a list of
    (incidence, value, marker, legend label). Empty per load otherwise."""
    out = {"CL": [], "CM": [], "CD": []}
    if abs(M - stations()[-1].M) > MAP_MACH_BAND:
        return out
    for Q, col in (("CL", "Cl"), ("CM", "Cm_c4"), ("CD", "Cd")):
        d = pd.read_csv(_use(DATA/f"static_naca0012_M030_{Q}.csv"))
        d = d[(d.alpha_deg >= lo) & (d.alpha_deg <= hi)]
        if Q != "CD":
            out[Q].append((d.alpha_deg, d[col], "o", "static, measured points"))
            continue
        wake = d.kind.str.startswith("total")
        out[Q].append((d.alpha_deg[wake], d[col][wake], "D", "static drag, measured: wake survey (total drag)"))
        out[Q].append((d.alpha_deg[~wake], d[col][~wake], "v", "static drag, measured: surface pressure (form drag)"))
    return out


def fig07_dynamic_and_static():
    g = np.arange(0.0, 20.01, 0.25)
    fig, axs = plt.subplots(len(STATIC_CASES), 3, figsize=(11.5, 3.9*len(STATIC_CASES)), constrained_layout=True)
    for i, (label, case) in enumerate(STATIC_CASES):
        o = _solve(**case)
        sm = _static(case["M"])
        curve = {"CL": sm.CN_static(g)*np.cos(np.radians(g)), "CM": sm.cm_static(g)}
        points = _static_points(case["M"], g[0], g[-1])
        for j, (ax, (_q, Q, yl)) in enumerate(zip(axs[i], QUANTITIES, strict=False)):
            ax.plot(o["alpha_deg"], o[Q], "-", color=BLUE, label="model, dynamic")
            _arrows(ax, o["alpha_deg"], o[Q], BLUE)
            if Q in curve:
                ax.plot(g, curve[Q], "--", color=INK_SOFT, label=r"static input of the model ($C_N\cos\alpha$; $C_{M,c/4}$)")
            for a, y, mk, lab in points[Q]:
                ax.plot(a, y, mk, ms=4.2, mfc="white", mec=INK_SOFT, mew=1.0, ls="none", label=lab)
            ax.set_xlabel(r"angle of attack $\alpha$, deg")
            ax.set_ylabel(yl)
            ax.set_title(_letter(3*i + j) + (" static drag: measured, not an input of the model's drag"
                                             if Q == "CD" else ""), loc="left", fontsize=9)
        axs[i, 0].set_title(rf"{_letter(3*i)} {label}: $\alpha$ = {case['alpha_mean_deg']:.0f}° ± "
                            f"{case['alpha_amp_deg']:.0f}°, k = {case['k']:.2f}, M = {case['M']:.2f}", loc="left")
    _legend_below(fig, axs[0], ncol=3)
    _save(fig, "fig07_dynamic_and_static")


# --------------------------------------------------------------------------- #
#  the two tables this script owns
# --------------------------------------------------------------------------- #
_BUILT = set()             # the owned tables computed in this run, so --refresh computes each once


def _signature(paths, grid):
    """One hash of the files a table is computed from and of the grid it is computed on."""
    record = [{_name(p): file_hash(p) for p in paths if Path(p).exists()}, grid]
    return hashlib.sha256(json.dumps(record, sort_keys=True).encode()).hexdigest()


def _owned(name, build, signature):
    """A table this script owns: read if it was computed from the present
    inputs on the present grid, computed again otherwise."""
    f = RESULTS/name
    if f.exists() and (not REFRESH or name in _BUILT):
        df = pd.read_csv(f)
        if SIGNATURE in df and (df[SIGNATURE] == signature).all():
            _use(f)
            return df
    df = build()
    _BUILT.add(name)
    _use(f)
    return df


def _in_parallel(point, jobs):
    """point(job) for every job, on all processors, in the order of the jobs."""
    with ProcessPoolExecutor(max_workers=WORKERS) as pool:
        return list(pool.map(point, jobs, chunksize=max(1, len(jobs)//(8*WORKERS))))


def map_loops() -> pd.DataFrame:
    """The measured loops at the condition of the design-space maps: rows of
    the manifest within MAP_MACH_BAND of MAP_MACH [-] and MAP_AMP_BAND of the
    amplitude MAP_AMP [deg]."""
    m = _manifest()
    return m[((m.M - MAP_MACH).abs() <= MAP_MACH_BAND) & ((m.amp_deg - MAP_AMP).abs() <= MAP_AMP_BAND)].reset_index(drop=True)


def map_grid() -> tuple:
    """(mean incidences [deg], reduced frequencies [-]) of the design-space
    maps: the spans MAP_MEAN_SPAN and MAP_K_SPAN, widened where a measured
    loop of `map_loops` lies beyond them, so the maps cover every loop that is
    marked on them."""
    m = map_loops()
    a = [*MAP_MEAN_SPAN, *m.alpha0_deg]
    k = [*MAP_K_SPAN, *m.k]
    means = np.arange(np.floor(min(a)), np.ceil(max(a)) + MAP_MEAN_STEP/2, MAP_MEAN_STEP)
    inner = np.arange(np.ceil(min(k)/MAP_K_STEP - 1e-9), np.floor(max(k)/MAP_K_STEP + 1e-9) + 0.5)*MAP_K_STEP
    return means.round(4), np.unique(np.r_[min(k), inner, max(k)].round(4))


def _map_signature():
    means, ks = map_grid()
    return _signature([dm.CONSTANTS_FILE, mt.MANIFEST, *STATIC_FILES], [means.tolist(), ks.tolist(), MAP_MACH, MAP_AMP])


def _map_point(job):
    """One point of the design-space maps: (mean incidence [deg], reduced frequency)."""
    a0, k = job
    o = dm.solve(a0, MAP_AMP, k, MAP_MACH)
    return dict(alpha_mean_deg=a0, alpha_amp_deg=MAP_AMP, M=MAP_MACH, k=k,
                CL_max=round(float(o["CL"].max()), 4), CM_min=round(float(o["CM"].min()), 4),
                CD_max=round(float(o["CD"].max()), 4),
                cycle_damping=round(mt.cycle_damping(o["alpha_deg"], o["CM"], MAP_AMP), 5))


def design_space() -> pd.DataFrame:
    """Peaks of the predicted loads and the cycle damping over mean incidence
    [deg] and reduced frequency [-] at MAP_MACH and MAP_AMP, on `map_grid`;
    written to results/design_space.csv with the signature of its inputs."""
    means, ks = map_grid()
    df = pd.DataFrame(_in_parallel(_map_point, [(float(a0), float(k)) for a0 in means for k in ks]))
    df[SIGNATURE] = _map_signature()
    df.to_csv(RESULTS/"design_space.csv", index=False)
    return df


def _design_space():
    return _owned("design_space.csv", design_space, _map_signature())


def surface_grid() -> tuple:
    """(T_VL values, T_p values) [semichords] of the cost surface: even grids
    over SURFACE_TVL_SPAN and SURFACE_TP_SPAN, widened to hold the calibrated
    values, with the calibrated value of each added, so the calibrated point
    is itself a point of the surface."""
    c = dm.load_constants()

    def axis(span: tuple, n: int, value: float) -> np.ndarray:
        """An even grid over the span and the value, with the value in place of a grid point it nearly coincides with."""
        even = np.linspace(min(span[0], value), max(span[1], value), n)
        even = even[np.abs(even - value) > 0.25*(even[1] - even[0])]
        return np.unique(np.r_[even, value].round(6))
    return axis(SURFACE_TVL_SPAN, SURFACE_TVL_N, c["Tvl"]), axis(SURFACE_TP_SPAN, SURFACE_TP_N, c["Tp"])


def _surface_signature():
    tvl, tp = surface_grid()
    frames = [mt.FRAME_CACHE/f"{f}.mat" for f in mt.manifest_frames(sets=("calibration",)).frame]
    return _signature([dm.CONSTANTS_FILE, mt.MANIFEST, *STATIC_FILES, *frames], [tvl.tolist(), tp.tolist()])


_CALIBRATION = {}          # the calibration loops and constants, loaded once per process


def _surface_point(job):
    """The calibration cost at one (T_VL, T_p) [semichords], other constants as calibrated."""
    if not _CALIBRATION:
        _m, frames, statics = cal.calibration_frames()
        _CALIBRATION.update(frames=frames, statics=statics, base=dm.load_constants())
    c = _CALIBRATION
    return dict(Tvl=job[0], Tp=job[1], cost=cal.cost(cal.residuals(np.array(job), ["Tvl", "Tp"], c["frames"],
                                                                   c["statics"], c["base"])))


def cost_surface() -> pd.DataFrame:
    """Calibration cost over the vortex travel time T_VL and the pressure lag
    T_p [semichords] on `surface_grid`, other constants as calibrated; written
    to results/cost_surface.csv with the signature of its inputs."""
    tvl, tp = surface_grid()
    _CALIBRATION.clear()
    df = pd.DataFrame(_in_parallel(_surface_point, [(float(a), float(b)) for a in tvl for b in tp]))
    df[SIGNATURE] = _surface_signature()
    df.to_csv(RESULTS/"cost_surface.csv", index=False)
    return df


def stale_tables() -> list:
    """Names of the two tables this script owns that were not computed from
    the present inputs on the present grid (or are absent)."""
    out = []
    for name, signature in (("design_space.csv", _map_signature), ("cost_surface.csv", _surface_signature)):
        f = RESULTS/name
        df = pd.read_csv(f) if f.exists() else pd.DataFrame()
        if SIGNATURE not in df or not (df[SIGNATURE] == signature()).all():
            out.append(name)
    return out


# --------------------------------------------------------------------------- #
#  maps
# --------------------------------------------------------------------------- #
def fig08_design_space():
    d = _design_space()
    fig, axs = plt.subplots(1, 3, figsize=(13.2, 4.3), constrained_layout=True)
    for i, (ax, (q, lab, name)) in enumerate(zip(axs, (("CL_max", "$C_{L,max}$", "Peak lift"),
                                                       ("CM_min", "$C_{M,min}$", "Minimum moment"),
                                                       ("CD_max", "$C_{D,max}$", "Peak drag")), strict=False)):
        Z = d.pivot(index="k", columns="alpha_mean_deg", values=q)
        _field(fig, ax, Z.columns.values, Z.index.values, Z.values, lab, reverse=q == "CM_min")
        ax.set_xlabel(r"mean incidence $\alpha_{mean}$, deg")
        ax.set_ylabel("reduced frequency $k$")
        ax.set_title(f"{_letter(i)} {name} {lab}", loc="left")
    fig.suptitle(rf"Predicted peaks at M = {MAP_MACH:.2f}, amplitude {MAP_AMP:.0f}°; the model is run every "
                 rf"{MAP_MEAN_STEP:g}° in mean incidence and every {MAP_K_STEP:g} in $k$ "
                 f"({d.alpha_mean_deg.nunique()} × {d.k.nunique()} points). Lighter is larger in magnitude", fontsize=10)
    _save(fig, "fig08_design_space")


def damping_map_loops() -> pd.DataFrame:
    """The measured loops marked on the damping map: those of `map_loops` with
    a closed moment loop, with their measured damping (Xi_exp) and their set."""
    t = _table()
    t = t[t.frame.isin(set(map_loops().frame)) & t.Xi_exp.notna()]
    return t[["frame", "set", "k", "alpha0_deg", "Xi_exp"]].reset_index(drop=True)


def fig09_damping_map():
    d = _design_space()
    Z = d.pivot(index="alpha_mean_deg", columns="k", values="cycle_damping")
    x, y, z = Z.columns.values, Z.index.values, Z.values
    t = damping_map_loops()
    fig, ax = plt.subplots(figsize=(7.8, 6.0), constrained_layout=True)
    levels = ticker.MaxNLocator(nbins=9).tick_values(float(z.min()), float(z.max()))
    extra = []
    if z.min() < 0:
        ax.contourf(x, y, z, levels=[2*float(z.min()) - 1.0, 0.0], colors=["#e9eef4"], hatches=["////"])
        extra.append(Patch(facecolor="#e9eef4", edgecolor=INK, hatch="////", lw=0.0,
                           label="predicted damping negative (the flow feeds the motion)"))
    cs = ax.contour(x, y, z, levels=[v for v in levels if abs(v) > 1e-12], colors=[BLUE], linewidths=1.0,
                    negative_linestyles="dashed")
    if z.min() < 0 < z.max():
        ax.contour(x, y, z, levels=[0.0], colors=[INK], linewidths=2.0)
        ax.plot([], [], "-", color=INK, lw=2.0, label="predicted zero damping")
    ax.plot([], [], "-", color=BLUE, lw=1.0, label=f"predicted damping, every {levels[1] - levels[0]:g} (dashed: negative)")
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(y.min(), y.max())
    for text in ax.clabel(cs, fmt=lambda v: f"{v:g}".replace("-", "\u2212"), fontsize=8.5, colors=[INK],
                          manual=_label_points(cs, ax, t[["k", "alpha0_deg"]].values)):
        text.set_bbox(dict(facecolor="white", edgecolor="none", pad=1.2))          # readable on the hatching
    for sign, mk, word in ((1, "o", "positive"), (-1, "v", "negative")):
        for name, style in SET_STYLE.items():
            m = (np.sign(t.Xi_exp) == sign) & (t["set"] == name)
            ax.plot(t[m].k, t[m].alpha0_deg, mk, ms=7.5, mfc=style["mfc"], mec=INK, mew=1.3, ls="none", clip_on=False,
                    zorder=4, label=f"{style['label']}, measured damping {word}")
    ax.set_xlabel("reduced frequency $k$")
    ax.set_ylabel(r"mean incidence $\alpha_{mean}$, deg")
    ax.set_title(rf"Predicted cycle damping at M = {MAP_MACH:.2f}, amplitude {MAP_AMP:.0f}°, and the {len(t)} measured "
                 "loops at that condition", loc="left")
    _legend_below(fig, ax, ncol=2, extra=extra)
    _save(fig, "fig09_damping_map")


def compared_ranges() -> tuple:
    """((lowest, highest reduced frequency), (lowest, highest Mach number)) of
    the measured loops the model is compared with: the manifest's usable
    NACA 0012 loops inside the static Mach range."""
    m = _manifest()
    return (float(m.k.min()), float(m.k.max())), (float(m.M.min()), float(m.M.max()))


def trend_field() -> pd.DataFrame:
    """The rows of results/trend_sweep.csv inside `compared_ranges`; nothing
    outside them is drawn."""
    (k0, k1), (M0, M1) = compared_ranges()
    t = _csv("trend_sweep.csv")
    eps = 1e-9
    return t[(t.k >= k0 - eps) & (t.k <= k1 + eps) & (t.M >= M0 - eps) & (t.M <= M1 + eps)].reset_index(drop=True)


def trend_loops(alpha_mean_deg: float, alpha_amp_deg: float) -> pd.DataFrame:
    """The measured loops within TREND_MATCH_DEG [deg] of a mean incidence and
    an amplitude [deg], whatever their Mach number and frequency."""
    m = _manifest()
    return m[((m.alpha0_deg - alpha_mean_deg).abs() <= TREND_MATCH_DEG)
             & ((m.amp_deg - alpha_amp_deg).abs() <= TREND_MATCH_DEG)].reset_index(drop=True)


def fig10_mach_trends():
    t = trend_field()
    (k0, k1), (M0, M1) = compared_ranges()
    conds = list(t.groupby(["alpha_mean_deg", "alpha_amp_deg"]).groups)
    quantities = (("CL_max", "$C_{L,max}$"), ("CM_min", "$C_{M,min}$"), ("CD_max", "$C_{D,max}$"))
    fig, axs = plt.subplots(len(conds), 3, figsize=(12.5, 3.6*len(conds)), constrained_layout=True, squeeze=False)
    for i, (a0, a1) in enumerate(conds):
        d = t[(t.alpha_mean_deg == a0) & (t.alpha_amp_deg == a1)]
        loops = trend_loops(a0, a1)
        for j, (ax, (q, lab)) in enumerate(zip(axs[i], quantities, strict=False)):
            Z = d.pivot(index="M", columns="k", values=q)
            _field(fig, ax, Z.columns.values, Z.index.values, Z.values, lab, reverse=q == "CM_min")
            for (name, style), mk in zip(SET_STYLE.items(), ("D", "o"), strict=True):
                s = loops[loops["set"] == name]
                ax.plot(s.k, s.M, mk, ms=6, mfc=style["mfc"], mec=INK, mew=1.1, ls="none", zorder=4,
                        label=f"{style['label']} measured at this mean incidence and amplitude")
            ax.set_xlim(k0 - 0.04*(k1 - k0), k1 + 0.04*(k1 - k0))      # room for the symbols at the ends of the span
            ax.set_ylim(M0 - 0.04*(M1 - M0), M1 + 0.04*(M1 - M0))
            ax.set_xlabel("reduced frequency $k$")
            ax.set_ylabel("Mach number")
            ax.set_title(rf"{_letter(3*i + j)} {lab}, $\alpha$ = {a0:.0f}° ± {a1:.0f}°: {len(loops)} measured loops",
                         loc="left")
    fig.suptitle(f"The model inside the span of all the measured loops, k = {k0:g} to {k1:g} and M = {M0:g} to {M1:g}. "
                 "Away from the symbols no loop was measured at that mean incidence and amplitude; the map there is "
                 "the model alone", fontsize=10)
    _legend_below(fig, axs[0, 0], ncol=2)
    _save(fig, "fig10_mach_trends")


# --------------------------------------------------------------------------- #
#  scores
# --------------------------------------------------------------------------- #
def fig11_error_by_condition():
    p = _primary()
    depth = _csv("validation_by_stall_depth.csv")
    fig, axs = plt.subplots(1, 3, figsize=(12.5, 4.0), constrained_layout=True)
    for ax, x, xl in ((axs[0], p.beyond_stall_deg, "peak incidence beyond static stall, deg"),
                      (axs[1], p.k, "reduced frequency $k$")):
        for Q, lab, col, mk, _ls, _hatch in SERIES:
            ax.plot(x, p[f"nRMS_{Q}"], mk, ms=4.5, mfc="none", mec=col, mew=0.9, label=lab)
        ax.set_xlabel(xl)
        ax.set_ylabel("normalised loop error")
    x = np.arange(len(depth))
    for j, (Q, lab, col, _mk, _ls, hatch) in enumerate(SERIES):
        _bars(axs[2], x + (j - 1)*0.26, depth[f"nRMS_{Q}"], 0.26, col, hatch, f"{lab} (bars)")
    axs[2].set_xticks(x, [f"{g.split(' (')[0]}\n({n} loops)" for g, n in zip(depth.group, depth.n_frames, strict=False)])
    axs[2].set_ylabel("mean normalised loop error")
    axs[2].set_xlabel("stall-depth group of the held-out loops")
    for ax, ttl in zip(axs, ("(a) Against stall depth", "(b) Against frequency", "(c) By stall-depth group"), strict=False):
        ax.set_title(ttl, loc="left")
    _legend_below(fig, axs, ncol=6)
    _save(fig, "fig11_error_by_condition")


def fig12_separation_laws():
    own = _csv("comparison.csv")
    same = _csv("comparison_same_constants.csv")
    own.insert(0, "constants", "each model's own calibration")
    d = pd.concat([own, same], ignore_index=True)
    metrics = [("mean |nRMS_CL|", "Lift loop error", "mean normalised loop error in $C_L$ [-]"),
               ("mean |nRMS_CM|", "Moment loop error", "mean normalised loop error in $C_{M,c/4}$ [-]"),
               ("mean |nRMS_CD|", "Drag loop error", "mean normalised loop error in $C_D$ [-]"),
               ("mean |CLmax_err_pct|", "Peak lift", "mean error in $C_{L,max}$, % of the measured peak"),
               ("mean |dalpha_CLmax_deg|", "Incidence of maximum lift", "mean error in that incidence, deg"),
               ("mean |dalpha_Mstall_deg|", "Incidence of moment stall", "mean error in that incidence, deg")]
    sets = list(d.constants.unique())
    fig, axs = plt.subplots(2, 3, figsize=(12.5, 6.0), constrained_layout=True)
    for n, (ax, (m, name, quantity)) in enumerate(zip(axs.ravel(), metrics, strict=False)):
        for i, (s, col, mk) in enumerate(zip(sets, (RED, BLUE, GREEN, PURPLE), ("o", "s", "^", "D"), strict=False)):
            r = d[(d.constants == s) & (d.metric == m)].iloc[0]
            ax.errorbar(r.paired_difference, i, xerr=[[r.paired_difference - r.ci95_lo], [r.ci95_hi - r.paired_difference]],
                        fmt=mk, ms=5, color=col, capsize=3.5, label=f"constants: {s}")
        ax.axvline(0.0, color=INK_SOFT, lw=0.9)
        ax.set_yticks([])
        ax.set_ylim(len(sets) - 0.5, -0.5)
        ax.set_xlabel(f"difference, tabulated minus fitted, in the\n{quantity}")
        ax.set_title(f"{_letter(n)} {name}", loc="left")
    fig.suptitle("Paired difference over the held-out loops with its 95 % interval; left of zero the tabulated law is "
                 "the more accurate", fontsize=10)
    _legend_below(fig, axs[0, 0], ncol=2)
    _save(fig, "fig12_separation_laws")


def fig13_cycle_damping():
    p = _primary()
    p = p[p.Xi_exp.notna() & p.Xi_model.notna()]
    groups = (("attached or marginal", p.beyond_stall_deg < DEPTH_LIGHT, "o", BLUE),
              ("light stall", (p.beyond_stall_deg >= DEPTH_LIGHT) & (p.beyond_stall_deg < DEPTH_DEEP), "s", ORANGE),
              ("deep stall", p.beyond_stall_deg >= DEPTH_DEEP, "^", RED))

    def span(v: pd.Series) -> tuple:
        """The range of the data and of zero, with a margin."""
        lo, hi = min(float(v.min()), 0.0), max(float(v.max()), 0.0)
        return lo - 0.08*(hi - lo), hi + 0.08*(hi - lo)
    (x0, x1), (y0, y1) = span(p.Xi_exp), span(p.Xi_model)
    fig, ax = plt.subplots(figsize=(6.2, 5.6), constrained_layout=True)
    wrong = dict(facecolor=RED, alpha=0.10, edgecolor=RED, hatch="xxx", lw=0.0)
    ax.fill_between([x0, 0], 0, y1, **wrong)
    ax.fill_between([0, x1], y0, 0, label="sign predicted wrongly", **wrong)
    both = (max(x0, y0), min(x1, y1))
    ax.plot(both, both, "-", color=INK_SOFT, lw=0.9, label="exact agreement")
    ax.axhline(0, color=INK_SOFT, lw=0.7)
    ax.axvline(0, color=INK_SOFT, lw=0.7)
    for lab, m, mk, col in groups:
        ok = f"{float(p[m].Xi_sign_agree.mean()):.0%}" if m.any() else "none"
        ax.plot(p[m].Xi_exp, p[m].Xi_model, mk, ms=6, mfc="none", mec=col, mew=1.2,
                label=f"{lab}: sign correct on {ok} of {int(m.sum())}")
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_xlabel("measured cycle damping [-]")
    ax.set_ylabel("predicted cycle damping [-]")
    ax.set_title(f"The {len(p)} primary held-out loops with a closed moment loop", loc="left")
    _legend_below(fig, ax, ncol=1)
    _save(fig, "fig13_cycle_damping")


# --------------------------------------------------------------------------- #
#  calibration
# --------------------------------------------------------------------------- #
LAWS = (("", "tabulated", BLUE, "o", HATCHES[0]), ("_reference", "fitted law", ORANGE, "s", HATCHES[1]))


def _starts_panel(ax):
    for suffix, lab, col, mk, _hatch in LAWS:
        r = _csv(f"calibration_runs{suffix}.csv").sort_values("search_cost").reset_index(drop=True)
        ax.plot(range(1, len(r) + 1), r.search_cost, mk + ("-" if suffix == "" else "--"), ms=5, color=col,
                label=f"{lab}, at the end of the search")
        pol = r[r.polished]
        ax.plot(pol.index + 1, pol.cost, mk, ms=11, mfc="none", mec=col, mew=1.6, label=f"{lab}, refined at full resolution")
    ax.set_xlabel("start, in order of the cost at the end of its search")
    ax.set_ylabel("calibration cost")
    ax.set_title("(a) Cost reached from each start", loc="left")
    ax.legend(frameon=False, fontsize=8.5)


def _sensitivity_panel(ax):
    pct = _targets()["calibration"]["sensitivity_perturbation_pct"]
    s = _csv("sensitivity.csv")
    base = float(s.cost[s.constant == "(calibrated)"].iloc[0])
    s = s[s.constant != "(calibrated)"]
    names = list(s.constant.unique())
    for sign, col, hatch in ((-1, BLUE, HATCHES[0]), (1, RED, HATCHES[1])):
        e = s[s.change_pct == sign*pct].set_index("constant").reindex(names)
        _bars(ax, np.arange(len(names)) + 0.2*sign, 100*(e.cost - base)/base, 0.4, col, hatch,
              f"{sign*pct:+g} % on the constant")
    ax.set_xticks(range(len(names)), names)
    ax.set_xlabel("constant changed, the others as calibrated")
    ax.set_ylabel("change in cost, %")
    ax.set_title("(b) Sensitivity of the cost", loc="left")
    ax.legend(frameon=False, loc="upper left")


def _cross_validation_panel(ax):
    for suffix, lab, col, mk, hatch in LAWS:
        cv = _csv(f"cross_validation{suffix}.csv")
        order = [x for x in ("all", "sensitive only", "none") if x in set(cv.constants_fitted)]
        tot = cv.groupby("constants_fitted").left_out_cost.sum()
        off = -0.2 if suffix == "" else 0.2
        ax.bar(np.arange(len(order)) + off, tot[order], 0.4, facecolor=col, alpha=0.55, edgecolor="white", lw=0.0,
               hatch=hatch, label=f"{lab}, total")
        for n, name in enumerate(order):
            folds = cv[cv.constants_fitted == name].left_out_cost.values
            ax.plot(np.full(len(folds), n + off), folds, mk, ms=4, mfc="white", mec=INK, mew=0.9, ls="none",
                    label=f"{lab}, single folds" if n == 0 else None)
        ax.set_xticks(range(len(order)), [f"constants fitted:\n{x}" for x in order])
    ax.set_ylabel("cost on the left-out calibration loops")
    ax.set_title("(c) Cross-validation inside the calibration set", loc="left")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2, fontsize=8.5)


def _surface_panel(fig, ax):
    s = _owned("cost_surface.csv", cost_surface, _surface_signature())
    z = s.pivot(index="Tp", columns="Tvl", values="cost")
    _field(fig, ax, z.columns.values, z.index.values, z.values, "calibration cost", reverse=True)
    ax.plot(s.Tvl, s.Tp, ".", ms=2.2, color=INK_SOFT, ls="none", label="computed points")
    k = _constants()
    ax.plot(k["Tvl"], k["Tp"], "*", ms=15, mfc="white", mec=INK, mew=1.2, ls="none", label="calibrated constants")
    low = s.loc[s.cost.idxmin()]
    if abs(low.Tvl - k["Tvl"]) > 1e-5 or abs(low.Tp - k["Tp"]) > 1e-5:
        ax.plot(low.Tvl, low.Tp, "X", ms=9, mfc="white", mec=INK, mew=1.1, ls="none", label="lowest computed cost")
        where = "the lowest computed cost is not at the calibrated point"
    else:
        where = "the lowest computed cost is at the calibrated point"
    ax.set_xlabel("vortex travel time $T_{vl}$, semichords")
    ax.set_ylabel("pressure lag $T_p$, semichords")
    ax.set_title(f"(d) Cost surface: {where}", loc="left")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=3, fontsize=8.5)


def fig14_calibration():
    fig, axs = plt.subplots(2, 2, figsize=(11.6, 8.6), constrained_layout=True)
    _starts_panel(axs[0, 0])
    _sensitivity_panel(axs[0, 1])
    _cross_validation_panel(axs[1, 0])
    _surface_panel(fig, axs[1, 1])
    _save(fig, "fig14_calibration")


def fig15_model_states():
    o = _solve(**STATE_CASE)
    ph = np.linspace(0.0, 360.0, len(o["alpha_deg"]))
    fig, axs = plt.subplots(5, 1, figsize=(8.8, 10.8), sharex=True, constrained_layout=True)
    axs[0].plot(ph, o["alpha_deg"], color=BLUE, label=r"incidence $\alpha$")
    axs[0].plot(ph, o["alpha_f_deg"], "--", color=RED, label=r"lagged effective incidence $\alpha_f$")
    axs[0].set_ylabel("incidence, deg")
    axs[1].plot(ph, o["CN"], color=BLUE, label="$C_N$")
    axs[1].plot(ph, o["CN_prime"], "--", color=RED, label="lagged $C'_N$")
    axs[1].axhline(o["meta"]["CN1"], color=GREEN, ls=":", label="onset level $C_{N1}$")
    axs[1].set_ylabel("normal force [-]")
    axs[2].plot(ph, o["f_static"], "--", color=RED, label="static $f'$ at $\\alpha_f$")
    axs[2].plot(ph, o["f_sep"], color=BLUE, label="lagged $f''$")
    axs[2].set_ylabel("separation point [-]")
    axs[3].plot(ph, o["CN_vortex"], color=PURPLE, label="vortex lift $C_N^v$")
    axs[3].plot(ph, o["tau_v"]/max(o["tau_v"].max(), 1e-9), ":", color=ORANGE, label="vortex clock, scaled to 1")
    axs[3].set_ylabel("vortex [-]")
    axs[4].plot(ph, o["CL"], color=BLUE, label="$C_L$")
    axs[4].plot(ph, 10*o["CM"], "--", color=RED, label="10 × $C_{M,c/4}$")
    axs[4].plot(ph, o["CD"], ":", color=GREEN, label="$C_D$")
    axs[4].set_ylabel("loads [-]")
    axs[4].set_xlabel("phase of the cycle, deg")
    c = STATE_CASE
    titles = (rf"incidence; $\alpha$ = {c['alpha_mean_deg']:.0f}° ± {c['alpha_amp_deg']:.0f}°, k = {c['k']}, M = {c['M']}",
              "normal force", "separation point",
              f"vortex: the clock starts {o['n_sheddings']} times in the cycle, {o['n_repeated_sheddings']} of them by "
              "repeated shedding", "loads")
    for i, (ax, ttl) in enumerate(zip(axs, titles, strict=True)):
        ax.set_title(f"{_letter(i)} {ttl}", loc="left", fontsize=9.5)
        ax.legend(frameon=False, loc="center left", bbox_to_anchor=(1.01, 0.5))
    _save(fig, "fig15_model_states")


def fig16_accuracy_summary():
    tg = _csv("validation_targets.csv").set_index("measure")
    unc = _csv("uncertainty_summary.csv").set_index("metric")
    summ = _csv("validation_summary.csv")
    ci = summ[(summ.model == "tabulated") & (summ.group == "primary")].set_index("metric")
    items = (("mean_nRMS_CL", "nRMS_CL", "lift"), ("mean_nRMS_CM", "nRMS_CM", "moment"), ("mean_nRMS_CD", "nRMS_CD", "drag"))
    look = {True: ("target met", GREEN, HATCHES[0]), False: ("target not met", ORANGE, HATCHES[1])}
    fig, ax = plt.subplots(figsize=(7.2, 5.0), constrained_layout=True)
    ticks = []
    for i, (key, u, lab) in enumerate(items):
        v = tg.measured[key]
        state, col, hatch = look[bool(tg.met[key])]
        _bars(ax, i, v, 0.5, col, hatch, f"mean error, {state}")
        ax.errorbar(i, v, yerr=[[v - ci.ci95_lo[key]], [ci.ci95_hi[key] - v]], color=INK, capsize=5, lw=1.2,
                    label="95 % interval of the mean over loops")
        ax.plot([i - 0.36, i + 0.36], [tg.target[key]]*2, "-", color=RED, lw=2.0, label="target (at most)")
        ax.plot([i - 0.36, i + 0.36], [unc.measurement_floor[u]]*2, "--", color=BLUE, lw=1.6,
                label="error from measurement uncertainty alone")
        ticks.append(f"{lab}\n{v:.3f} against {tg.target[key]:g}\n{state}")
    ax.set_xticks(range(len(items)), ticks)
    ax.set_ylabel("mean normalised loop error, held-out loops [-]")
    _legend_below(fig, ax, ncol=2)
    _save(fig, "fig16_accuracy_summary")


def fig17_convergence():
    st = _csv("convergence_steps.csv")
    summ = _csv("convergence_summary.csv").set_index("study")
    T = _targets()["numerical"]["loop_error_change_max"]
    st = st[st.step_factor < st.step_factor.max()]              # the finest march is what the others are compared with
    fig, axs = plt.subplots(1, 2, figsize=(10.4, 4.5), constrained_layout=True)
    hidden, total = 0, 0
    for Q, _lab, col, mk, ls, _hatch in SERIES:
        v = st[f"max_abs_d{Q}"]
        hidden, total = hidden + int((v <= 0).sum()), total + len(v)
        axs[0].loglog(st.step_semichords[v > 0], v[v > 0], mk, ms=3.5, mfc="none", mec=col, mew=0.8, label=f"$C_{Q[1]}$")
        g = st.groupby("step_factor")[f"d_nRMS_{Q}"].max()
        g = g[g > 0]
        axs[1].loglog(g.index, g.values, mk + ls, ms=5, color=col, label=f"$C_{Q[1]}$")
    axs[0].axvline(dm.DS_MAX, color=INK_SOFT, ls=":", lw=1.2, label="largest step used")
    axs[0].set_xlabel("step, semichords of travel")
    axs[0].set_ylabel("largest change in the load over the cycle [-]")
    axs[0].set_title(f"(a) Loads against step size, each calibration loop\n{hidden} of {total} changes are zero to the "
                     "six decimals on file and are not drawn", loc="left", fontsize=9.5)
    axs[1].axhline(T, color=RED, lw=1.4, ls="-.", label=f"target ({T:g})")
    axs[1].axvline(1.0, color=INK_SOFT, ls=":", lw=1.2, label="resolution used")
    axs[1].set_xlabel("steps per cycle, as a multiple of those used")
    axs[1].set_ylabel("largest change in loop error, any loop [-]")
    axs[1].set_title("(b) Loop error against step count\n(change from the finest march)", loc="left", fontsize=9.5)
    for ax in axs:
        _log_ticks(ax.xaxis, *ax.get_xlim())
        _log_ticks(ax.yaxis, *ax.get_ylim())
        ax.legend(frameon=False, fontsize=8.5, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=5)
    cyc = summ.loc["cycles marched"]
    fig.suptitle(f"Cycles marched: against {cyc.default_compared_with} the loop error of any of the {int(cyc.n_frames)} "
                 f"calibration loops changes by at most {cyc.worst_change_in_loop_error:g} (target {T:g})", fontsize=10)
    _save(fig, "fig17_convergence")


def _state_space_cycle(axs, a, b):
    """Top row of the state-space figure and the difference over the cycle."""
    for i, (ax, (_q, Q, yl)) in enumerate(zip(axs[0], QUANTITIES, strict=False)):
        ax.plot(a["alpha_deg"], a[Q], "-", color=BLUE, lw=1.6, label="indicial march")
        ax.plot(b["alpha_deg"][::12], b[Q][::12], "o", ms=3.2, mfc="none", mec=RED, mew=0.9, label="state-space form")
        ax.set_xlabel(r"$\alpha$, deg")
        ax.set_ylabel(yl)
        ax.set_title(_letter(i), loc="left")
    c = STATE_CASE
    axs[0, 0].set_title(rf"(a) $\alpha$ = {c['alpha_mean_deg']:.0f}° ± {c['alpha_amp_deg']:.0f}°, "
                        f"k = {c['k']}, M = {c['M']}", loc="left")
    ph = np.linspace(0.0, 360.0, len(a["CL"]))
    for Q, _lab, col, _mk, ls, _hatch in SERIES:
        axs[1, 0].plot(ph, b[Q] - a[Q], ls, color=col, label=f"$C_{Q[1]}$")
    axs[1, 0].set_xlabel("phase of the cycle, deg")
    axs[1, 0].set_ylabel("state-space minus indicial [-]")
    axs[1, 0].set_title("(d) Difference over the same cycle", loc="left")


def fig18_state_space():
    from unistall import statespace
    chk = _csv("statespace_check.csv")
    T = _targets()["numerical"]["state_space_loop_error_difference_max"]
    a = _solve(**STATE_CASE)
    b = statespace.solve(**STATE_CASE, consts=_constants(), static=_static(STATE_CASE["M"]))
    fig, axs = plt.subplots(2, 3, figsize=(12.5, 7.6), constrained_layout=True)
    _state_space_cycle(axs, a, b)
    x = np.arange(len(chk))
    for j, (Q, _lab, col, _mk, _ls, hatch) in enumerate(SERIES):
        _bars(axs[1, 1], x + (j - 1)*0.27, (chk[f"nRMS_{Q}_state_space"] - chk[f"nRMS_{Q}_indicial"]).abs(), 0.27,
              col, hatch, f"$C_{Q[1]}$")
    axs[1, 1].axhline(T, color=INK_SOFT, ls="--", lw=1.2, label=f"target ({T:g})")
    axs[1, 1].set_xticks(x, [f[6:] for f in chk.frame], rotation=90, fontsize=6.5)
    axs[1, 1].set_xlabel("calibration loop (frame number)")
    axs[1, 1].set_ylabel("difference in loop error [-]")
    axs[1, 1].set_title(f"(e) Loop error, each of the {len(chk)} calibration loops", loc="left")
    axs[1, 2].plot(chk.damping_indicial, chk.damping_state_space, "o", ms=5, mfc="none", mec=BLUE, mew=1.1, ls="none",
                   label="one calibration loop")
    both = pd.concat([chk.damping_indicial, chk.damping_state_space])
    lim = [min(float(both.min()), 0.0) - 0.02, float(both.max()) + 0.02]
    axs[1, 2].plot(lim, lim, "-", color=INK_SOFT, lw=0.9, label="exact agreement")
    axs[1, 2].set_xlabel("cycle damping, indicial march [-]")
    axs[1, 2].set_ylabel("cycle damping, state-space form [-]")
    axs[1, 2].set_title(f"(f) Cycle damping, the same {len(chk)} loops", loc="left")
    for ax in (axs[0, 0], axs[1, 0], axs[1, 1], axs[1, 2]):
        ax.legend(frameon=False, fontsize=8.5, loc="upper center", bbox_to_anchor=(0.5, -0.24), ncol=4)
    _save(fig, "fig18_state_space")


FIGURES_ALL = (fig01_static_inputs, fig02_attached_flow, fig03_attached_moment, fig04_calibration_loops,
               fig05_held_out_loops, fig06_effect_of_constants, fig07_dynamic_and_static, fig08_design_space,
               fig09_damping_map, fig10_mach_trends, fig11_error_by_condition, fig12_separation_laws,
               fig13_cycle_damping, fig14_calibration, fig15_model_states, fig16_accuracy_summary,
               fig17_convergence, fig18_state_space)


def main(argv: list) -> None:
    """Draw the figures whose names start with a word of `argv` (all of them if there is none)."""
    journal_style()
    wanted = [a for a in argv if not a.startswith("--")]
    for make in FIGURES_ALL:
        if not wanted or any(make.__name__.startswith(w) for w in wanted):
            _READ.clear()
            make()


if __name__ == "__main__":
    main(sys.argv[1:])
