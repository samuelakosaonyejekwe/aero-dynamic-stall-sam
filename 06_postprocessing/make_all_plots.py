"""
06_postprocessing / make_all_plots.py
-------------------------------------
Generates the 2-D figures for the case study from the solution data, all
written to 06_postprocessing/plots/. The 3-D figures (fig3d_*.png) are NOT
produced here -- they come from make_3d_plots.py, and this stage deliberately
excludes them from its own count.

Outputs, by family. <case> is A_validation or B_application, <phase> is one of
rise|peak|fall|dsv and <deg> the incidence in whole degrees, matching the field
files run_case.py writes:
  hyst_cl_<case>.png / hyst_cd_<case>.png / hyst_cm_<case>.png
                                    hysteresis loops with stroke-direction arrows
  timehist_loads_<case>.png         unsteady loads against time
  states_<case>.png                 indicial state variables
  static_polar_calibration.png      model against the published reference polar
  convergence_residuals.png         cycle-to-cycle convergence (periodicity)
  timestep_refinement.png           time-step refinement (step adequacy)
  cp_distribution_<case>.png        surface Cp at the written phases
  contour_Cp_<case>_<phase>_a<deg>.png     pressure field
  contour_speed_stream_<case>_<phase>_a<deg>.png   speed magnitude + streamlines
  contour_vorticity_<case>_<phase>_a<deg>.png      vorticity
  contour_Mach_<case>_<phase>_a<deg>.png           local Mach number
  contour_Tstatic_<case>_<phase>_a<deg>.png        static temperature
  contour_Trecovery_<case>_<phase>_a<deg>.png      recovery (skin) temperature
  contour_vectors_<case>_<phase>_a<deg>.png        velocity vectors (quiver)
  temperature_profile_<case>.png    surface recovery temperature vs x/c, upper
                                    and lower surface, at the PEAK-incidence
                                    phase only (hence no phase in the name)
LAYOUT RULES (enforced): constrained_layout everywhere, colorbars on their own
axes, titles padded, legends in clear regions -> text never overlaps a figure.
No black is ever used (shared aero_style).
"""
import sys, glob
import numpy as np
import pandas as pd
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))          # for aero_style only: this stage reads
                                       # 05_solution CSVs and never imports the solver
from aero_style import (apply_style, PALETTE, INK, INK_SOFT,
                        CMAP_PRESSURE, CMAP_CP, CMAP_TEMP, CMAP_VORT)
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm as _SymLogNorm
from matplotlib.patches import Polygon as MplPoly
apply_style()
plt.rcParams["figure.constrained_layout.use"] = True

SOL = ROOT/"05_solution"; OUT = HERE/"plots"; OUT.mkdir(exist_ok=True)
SETUP = ROOT/"03_model_setup"; GEO = ROOT/"01_geometry"/"naca0012_coordinates.csv"
AF = pd.read_csv(GEO)

# case metadata is READ from 03_model_setup, never restated here
_flow = pd.read_csv(SETUP/"flow_conditions.csv").set_index("parameter")
_kin  = pd.read_csv(SETUP/"kinematics.csv").set_index("case_id")
def _case(col, label):
    f = _flow[col]; kr = _kin.loc[f["case_id"]]
    return dict(c=float(f["chord_c"]), U=float(f["freestream_velocity_U"]),
                M=float(f["freestream_mach_M"]),
                label=label % (float(f["freestream_mach_M"]), float(kr["reduced_freq_k"]),
                               float(kr["alpha_mean_deg"]), float(kr["alpha_amp_deg"])))
CASES = {"A_validation":  _case("case_A_validation",
                                "Case A — validation rig (NACA0012, M=%.2f, k=%.2f, α=%.0f°±%.0f°)"),
         "B_application": _case("case_B_application",
                                "Case B — rotor retreating blade (r/R=0.75, M=%.2f, k=%.3f, α=%.0f°±%.0f°)")}

# Every per-case figure carried a title that did not say which case it was, so
# the A and B versions of the hysteresis loops, the time histories, the internal
# states, the Cp distributions, the contours and the temperature profiles were
# titled identically -- 56 of the 90 published figures in indistinguishable
# pairs, and each PNG also ships on its own. The label that resolves it was
# already being built from the published flow conditions and kinematics at the
# top of this file, and then used nowhere. It is used here.
def case_title(ax, main, cs, pad=22):
    ax.set_title(main, pad=pad)
    ax.text(0.5, 1.012, CASES[cs]["label"], transform=ax.transAxes, ha="center",
            va="bottom", fontsize=8.5, color=INK_SOFT)


AIRFOIL_FC = "#e3e9f0"

def airfoil_patch(ax, c, fc=AIRFOIL_FC):
    poly = np.column_stack([AF["x_over_c"].values*c, AF["y_over_c"].values*c])
    ax.add_patch(MplPoly(poly, closed=True, facecolor=fc, edgecolor=INK, lw=1.3, zorder=5))

def stroke_arrows(ax, x, y, n=7, color=INK_SOFT):
    idx = np.linspace(2, len(x)-3, n).astype(int)
    for i in idx:
        ax.annotate("", xy=(x[i+1], y[i+1]), xytext=(x[i-1], y[i-1]),
                    arrowprops=dict(arrowstyle="-|>", color=color, lw=1.1), zorder=6)

def save(fig, name):
    fig.savefig(OUT/name); plt.close(fig)

# ============================================================ 1. HYSTERESIS
for cs, meta in CASES.items():
    th = pd.read_csv(SOL/f"time_history_{cs}.csv")
    a = th["alpha_deg"].values
    stat = pd.read_csv(SETUP/"static_polar_reference.csv")
    for var, lab, fname, cmcol in [("CL", "Lift coefficient $C_L$", "cl", "Cl"),
                                   ("CD", "Drag coefficient $C_D$", "cd", "Cd"),
                                   ("CM_c4", "Pitching moment $C_{M,c/4}$", "cm", "Cm_c4")]:
        fig, ax = plt.subplots(figsize=(6.2, 5))
        ax.plot(a, th[var], color=PALETTE[0], lw=2.2, label="UNISTALL (dynamic)")
        ax.plot(stat["alpha_deg"], stat[cmcol], "o--", color=PALETTE[1], ms=4,
                lw=1.3, label="static (published ref.)")
        stroke_arrows(ax, a, th[var].values)
        ax.set_xlabel("angle of attack  α  [deg]"); ax.set_ylabel(lab)
        case_title(ax, "Dynamic-stall hysteresis loop — "
                   + var.replace("_c4", " (c/4)"), cs)
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2, frameon=True)
        save(fig, f"hyst_{fname}_{cs}.png")

# ============================================================ 2. TIME HISTORIES
for cs, meta in CASES.items():
    th = pd.read_csv(SOL/f"time_history_{cs}.csv")
    t = th["time_s"].values*1e3
    # combined 4-panel loads vs time
    fig, axs = plt.subplots(4, 1, figsize=(8.2, 9), sharex=True)
    series = [("alpha_deg", "α [deg]", PALETTE[3]),
              ("CL", "$C_L$", PALETTE[0]),
              ("CD", "$C_D$", PALETTE[2]),
              ("CM_c4", "$C_{M,c/4}$", PALETTE[1])]
    for ax, (col, lab, col_c) in zip(axs, series):
        ax.plot(t, th[col], color=col_c, lw=2)
        ax.set_ylabel(lab)
    axs[-1].set_xlabel("time  [ms]")
    case_title(axs[0], "Unsteady load time histories (one converged cycle)", cs)
    save(fig, f"timehist_loads_{cs}.png")

    # state variables: separation point, vortex normal force, CN'
    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    ax.plot(th["phase_deg"], th["f_separation"], color=PALETTE[4], lw=2,
            label="separation point  f''")
    ax.plot(th["phase_deg"], th["CN_vortex"], color=PALETTE[1], lw=2,
            label="vortex normal force  $C_N^{v}$")
    ax.plot(th["phase_deg"], th["CN_prime"]/th["CN_prime"].max(), color=PALETTE[6],
            lw=1.5, ls="--", label="$C_N'$ (norm.)")
    ax.set_xlabel("cycle phase  ωt  [deg]"); ax.set_ylabel("state value")
    ax.set_ylim(-0.05, 1.18)
    case_title(ax, "UIBS internal states — separation & dynamic-stall vortex", cs)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=3, frameon=True)
    save(fig, f"states_{cs}.png")

# ============================================================ 3. STATIC POLAR
mp = pd.read_csv(SOL/"model_static_polar.csv")
stat = pd.read_csv(SETUP/"static_polar_reference.csv")
fig, axs = plt.subplots(1, 2, figsize=(11, 4.6))
# The model polar and f(alpha) run to 22 deg while the reference they are fitted
# to stops at ALPHA_CAL_MAX. Both curves used to be drawn as one solid line, so
# the last two degrees -- the EXTRAPOLATED branch of f(alpha), which no reported
# case uses -- looked exactly as calibrated as the rest. Drawn dashed beyond the
# calibration limit and stated, the same standard the response surface is held to.
_ACAL = float(stat["alpha_deg"].max())
_cal, _ext = mp["alpha_deg"] <= _ACAL, mp["alpha_deg"] >= _ACAL
axs[0].plot(mp["alpha_deg"][_cal], mp["Cl_model"][_cal], color=PALETTE[0], lw=2.2,
            label="UNISTALL model")
axs[0].plot(mp["alpha_deg"][_ext], mp["Cl_model"][_ext], color=PALETTE[0], lw=2.2,
            ls=(0, (4, 2)), label=f"extrapolated beyond {_ACAL:.0f}°")
axs[0].plot(stat["alpha_deg"], stat["Cl"], "s", color=PALETTE[1], ms=5, label="published ref.")
axs[0].set_xlabel("α [deg]"); axs[0].set_ylabel("$C_L$"); axs[0].legend(fontsize=9)
axs[0].set_title("Static lift polar (calibration)", pad=10)
axs[1].plot(mp["alpha_deg"][_cal], mp["f_sep"][_cal], color=PALETTE[2], lw=2.2)
axs[1].plot(mp["alpha_deg"][_ext], mp["f_sep"][_ext], color=PALETTE[2], lw=2.2,
            ls=(0, (4, 2)))
axs[1].axvline(_ACAL, color=INK_SOFT, lw=0.9, ls=":")
# ha="right", inside the axes. Anchored left of the line it ran off the right
# edge of the frame -- the calibration limit is at 20 deg and the axis stops at
# 22, so a left-anchored two-word line had nowhere to go. The band above the
# curve to the LEFT of the limit is empty (f is below 0.2 there), so it fits.
axs[1].text(_ACAL - 0.4, 0.97, f"calibrated to {_ACAL:.0f}°;\ndashed is extrapolation",
            transform=axs[1].get_xaxis_transform(), va="top", ha="right",
            fontsize=8.5, color=INK_SOFT)
axs[1].set_xlabel("α [deg]"); axs[1].set_ylabel("separation point  f")
axs[1].set_title("Calibrated static separation  f(α)", pad=10)
save(fig, "static_polar_calibration.png")

# ============================================================ 4. CONVERGENCE
# Once the cycle repeats, |ΔCLmax| is exactly 0 and simply disappears from a log
# axis -- which left 4 of 6 points off the chart and made it look broken. Floor
# the zeros and mark them as converged-to-machine-precision instead.
FLOOR = 1e-8
# Once both cases reach the floor their curves and markers coincide EXACTLY, so
# whichever is drawn second hid the other: cycles 5 and 6 showed one case where
# there are two. Distinct dash patterns and marker shapes let the lower curve
# show through the upper one, without moving either off its true value.
_LS = ["-", (0, (6, 3))]
_MK = ["o", "s"]
fig, ax = plt.subplots(figsize=(7.5, 4.6))
ncyc = 0
for i, cs in enumerate(CASES):
    rdf = pd.read_csv(SOL/"convergence"/f"residuals_{cs}.csv")
    x = rdf["cycle"].values; r = rdf["peakCL_residual"].values
    ncyc = max(ncyc, int(x.max()))
    m = ~np.isnan(r)
    conv = m & (r <= 0)
    rp = np.where(conv, FLOOR, r)
    ax.semilogy(x[m], rp[m], ls=_LS[i % 2], color=PALETTE[i], lw=2,
                label=cs.replace("_", " "))
    ax.plot(x[m & ~conv], rp[m & ~conv], _MK[i % 2], color=PALETTE[i], ms=6)
    ax.plot(x[conv], rp[conv], _MK[i % 2], ms=7 - 1.5*i, mfc="white",
            mec=PALETTE[i], mew=1.5)
ax.axhline(FLOOR, color=INK_SOFT, lw=0.8, ls=":")
ax.set_ylim(FLOOR/4, None); ax.set_xlim(0.7, ncyc + 0.3)
ax.set_xticks(range(1, ncyc + 1))
ax.set_xlabel("cycle number"); ax.set_ylabel("peak-$C_L$ residual  |ΔC$_{L,max}$|")
ax.set_title("Cycle-to-cycle convergence: the limit cycle is reached", pad=10)
# place the note in the empty mid-right band: at the bottom it sat on top of
# the converged floor line and its markers
ax.text(0.98, 0.34, "open symbols: ΔC$_{L,max}$ = 0 to double precision\n"
        "(cycle 1 has no predecessor)", transform=ax.transAxes, ha="right",
        va="top", fontsize=8.5, color=INK_SOFT,
        bbox=dict(boxstyle="round,pad=0.35", fc="white", ec=INK_SOFT, alpha=0.9))
ax.legend(loc="upper right")
save(fig, "convergence_residuals.png")

# --------------------------------------------------- 4b. TIME-STEP REFINEMENT
# The figure above shows the march reaches a limit cycle. It does NOT show the
# step resolves that cycle, and until this study was added nothing did.
tdf = pd.read_csv(SOL/"convergence"/"timestep_refinement.csv")
fig, (axL, axR) = plt.subplots(1, 2, figsize=(9.6, 4.3))
_cols = [("pct_from_finest_CL_max", "$C_{L,max}$"),
         ("pct_from_finest_CM_min", "$C_{M,min}$"),
         ("pct_from_finest_CD_max", "$C_{D,max}$")]
_n = tdf["steps_per_cycle"].values
_rep = tdf[tdf["reported_resolution"]].iloc[0]
for i, (col, lab) in enumerate(_cols):
    y = tdf[col].values.copy()
    y[y <= 0] = np.nan                      # the reference row is exactly zero
    axL.loglog(_n, y, ls=_LS[i % 2], marker=_MK[i % 2], color=PALETTE[i],
               lw=2, ms=6, label=lab)
# Deviation is measured against the FINEST run, which carries the same onset
# quantisation as every other row. So the curve cannot fall below that jitter and
# turns up again past ~1440 -- that upturn is the floor of the measurement, not
# the scheme diverging. Draw the floor so the figure cannot be read the wrong way.
_fine = tdf[tdf["steps_per_cycle"] > _rep["steps_per_cycle"]]
_floor = float(max(_fine[c].max() for c, _ in _cols))
_lo = float(min(np.nanmin(tdf[c].replace(0.0, np.nan)) for c, _ in _cols))/2.5
axL.set_ylim(_lo, None)
axL.axhspan(_lo, _floor, color=INK_SOFT, alpha=0.10, lw=0)
axL.axhline(_floor, color=INK_SOFT, lw=0.8, ls="--")
# the note goes under the axes, not in them: inside, it covered the 720-step
# marker, which is the one point the figure exists to show
fig.text(0.5, 0.015, "Shaded band is the onset-quantisation floor: the reference run is quantised too, so a "
         "difference below it says which step the trigger fired on, not how accurate the march is.",
         ha="center", va="bottom", fontsize=8.2, color=INK_SOFT)
axL.annotate("reference", xy=(_n[-1], _lo), xytext=(0, 4), textcoords="offset points",
             ha="center", va="bottom", fontsize=8.0, color=INK_SOFT, rotation=90)
for _ax in (axL, axR):
    _ax.axvline(_rep["steps_per_cycle"], color=INK_SOFT, lw=1.0, ls=":")
    _ax.set_xscale("log"); _ax.set_xticks(_n); _ax.set_xticklabels([str(int(v)) for v in _n])
    _ax.set_xticks([], minor=True)
    _ax.set_xlabel("steps per cycle")
axL.annotate("used here", xy=(_rep["steps_per_cycle"], 1.0), xycoords=("data", "axes fraction"),
             xytext=(-4, -6), textcoords="offset points", ha="right", va="top",
             fontsize=8.5, color=INK_SOFT)
axL.set_ylabel("deviation from finest run  [%]")
axL.set_title("Time-step refinement", pad=10)
axL.legend(loc="upper right", fontsize=9)

axR.plot(_n, tdf["stall_onset_alpha_deg"].values, ls=_LS[0], marker=_MK[0],
         color=PALETTE[0], lw=2, ms=6)
axR.set_ylabel(r"stall onset  $\alpha$  [deg]")
axR.set_title("Onset is detected at a discrete step", pad=10)
axR.text(0.97, 0.95, "the trigger $C_N' \\geq C_{N1}$ is tested once per step,\n"
         "so onset is quantised — this, not a scheme error,\n"
         "sets the residual scatter at the finest steps",
         transform=axR.transAxes, ha="right", va="top", fontsize=8.5,
         color=INK_SOFT,
         bbox=dict(boxstyle="round,pad=0.35", fc="white", ec=INK_SOFT, alpha=0.9))
fig.get_layout_engine().set(rect=(0, 0.055, 1, 0.945))
save(fig, "timestep_refinement.png")

# ============================================================ 5. CP DISTRIBUTION
for cs, meta in CASES.items():
    cp = pd.read_csv(SOL/f"cp_distribution_{cs}.csv")
    tags = cp["phase_tag"].unique()
    fig, ax = plt.subplots(figsize=(7.6, 5))
    for i, tg in enumerate(tags):
        sub = cp[cp["phase_tag"] == tg]
        adeg = sub["alpha_deg"].iloc[0]
        ax.plot(sub["x_c"], sub["Cp"], color=PALETTE[i % len(PALETTE)], lw=1.8,
                label=f"{tg}  (α={adeg:.1f}°)")
    ax.invert_yaxis()
    ax.set_xlabel("x/c"); ax.set_ylabel("$C_p$")
    case_title(ax, "Surface pressure coefficient at cycle phases", cs)
    # The last few control points carry a discretisation oscillation, not a flow
    # feature: the panels there are the shortest on the body and the imposed
    # circulation is not the Kutta one. Switching the dynamic-stall vortex off
    # barely changes it, so it is the panelling. Say so on the figure, with the
    # size read from the published metric rather than restated here.
    _mt = pd.read_csv(SOL/f"metrics_{cs}.csv").set_index("metric")["value"]
    _zone = float(_mt["Cp_TE_panel_oscillation_zone_x_c"])
    ax.axvspan(_zone, 1.0, color=INK_SOFT, alpha=0.12, lw=0)
    fig.get_layout_engine().set(rect=(0, 0.05, 1, 0.95))
    fig.text(0.5, 0.008, "shaded x/c > %.3f: panel-discretisation oscillation where the "
             "reconstruction closes the open trailing edge, |C_p| up to %.1f — not a flow feature"
             % (_zone, float(_mt["Cp_TE_panel_oscillation_max_abs"])),
             ha="center", va="bottom", fontsize=8, color=INK_SOFT)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2, fontsize=9, frameon=True)
    save(fig, f"cp_distribution_{cs}.png")

# ============================================================ field helpers
def load_field(path):
    df = pd.read_csv(path)
    xu = np.unique(df["x_m"].values); yu = np.unique(df["y_m"].values)
    nx, ny = len(xu), len(yu)
    flds = {c: df[c].values.reshape(ny, nx) for c in df.columns if c not in ("x_m", "y_m")}
    return xu, yu, flds

field_files = sorted(glob.glob(str(SOL/"field_*.csv")))

# Common Cp scale for the contour panels, DERIVED from the published fields
# rather than typed in. It was the literal (-5.0, 1.0) with a comment naming the
# deepest value any field reached; deriving the vortex made the fields deeper
# (-7.92) and a hardcoded floor would have silently clipped the very feature the
# figures are for. +1.0 is the stagnation bound Cp cannot exceed, which
# verify_invariants asserts; the floor is rounded down to the next half so the
# colorbar ticks stay round.
_cpmin = min(float(pd.read_csv(f, usecols=["Cp"])["Cp"].min()) for f in field_files)
CP_SCALE = (float(np.floor(_cpmin*2.0)/2.0), 1.0)
print(f"[plots] common Cp scale {CP_SCALE} (deepest published field value {_cpmin:.2f})")

def _symlog_ticks(lin, top):
    """0 and one tick per DECADE out to the peak -- a readable set for a
    symmetric-log bar, instead of one tick per contour level.

    The decades are true powers of ten, not multiples of the linear threshold:
    starting at the threshold put ticks at 0.72 and 7.2, which the formatter
    then printed as "1" and "7". A tick may be rounded in position, never in
    its label."""
    first = int(np.ceil(np.log10(lin)))
    last = int(np.floor(np.log10(top)))
    dec = [10.0**k for k in range(first, last+1)] or [top]
    return [-t for t in reversed(dec)] + [0.0] + dec


def contour_plot(xu, yu, Z, title, cbar_label, cmap, c, fname,
                 lines=False, levels=24, vector=None, stream=None, vlim=None,
                 note=None, case_key=None, norm=None, cticks=None):
    X, Y = np.meshgrid(xu, yu)
    fig, ax = plt.subplots(figsize=(7.8, 5.2))
    # The reconstruction masks the body AND the one ring of cells touching it
    # (the surface sheet is not resolved within a cell of the wall). contourf
    # draws nothing over NaN, so that ring let the white page show through as a
    # ragged halo around the aerofoil in all 56 contour figures. Painting the
    # axes background with the aerofoil fill makes the unresolved ring read as
    # part of the body, which is the honest reading: there is no field there.
    ax.set_facecolor(AIRFOIL_FC)
    # `levels` may arrive as an explicit array (the symmetric-log vorticity
    # scale), in which case it is used as given and vlim does not apply.
    if np.ndim(levels) > 0:
        lv = levels
    else:
        lv = np.linspace(vlim[0], vlim[1], levels) if vlim else levels
    cf = ax.contourf(X, Y, Z, levels=lv, cmap=cmap, extend="both", norm=norm)
    if lines:
        ax.contour(X, Y, Z, levels=12, colors=[INK_SOFT], linewidths=0.4, alpha=0.6)
    if stream is not None:
        u, v = stream
        xe = np.linspace(xu.min(), xu.max(), len(xu))
        ye = np.linspace(yu.min(), yu.max(), len(yu))
        ax.streamplot(xe, ye, np.nan_to_num(u), np.nan_to_num(v), density=1.1,
                      color=INK_SOFT, linewidth=0.6, arrowsize=0.8)
    if vector is not None:
        u, v = vector
        sk = (slice(None, None, 9), slice(None, None, 9))
        ax.quiver(X[sk], Y[sk], u[sk], v[sk], color=INK, scale_units="xy",
                  angles="xy", width=0.003)
    airfoil_patch(ax, c)
    cb = fig.colorbar(cf, ax=ax, pad=0.02, fraction=0.046)
    # Tick the colorbar on ROUND values. A filled-contour colorbar ticks every
    # contour level by default, so the percentile-clipped temperature and Mach
    # scales came out labelled 291.6854, 291.5930, ... -- four decimals of
    # spurious precision on a 0.65 K range, and long enough to crowd the axis.
    from matplotlib.ticker import MaxNLocator
    cb.locator = MaxNLocator(nbins=7, steps=[1, 2, 2.5, 5, 10])
    cb.update_ticks()
    # Drop any tick the locator puts OUTSIDE the colour limits. extend="both"
    # gives the bar a triangle at each end, and a round tick just beyond the
    # scale is drawn inside that triangle, hard against the last real tick: on
    # the percentile-clipped temperature bars, whose whole range can be 0.5 K,
    # the two labels then print on top of each other (291.2 over 291.1 on the
    # case-B recovery-temperature fields). Every one of the 24 percentile-scaled
    # bars was placing at least one tick out of range this way. The colour
    # limits are read back off the mappable so this holds for the fixed scales
    # too, not just the percentile ones.
    if cticks is not None:
        # An explicit set, for the symmetric-log vorticity bar. Left to itself
        # the colorbar puts a tick on every contour level, and the log-spaced
        # levels bunch at both ends: the labels printed on top of each other
        # (75000 over 50000 over 25000, twice).
        cb.set_ticks(cticks)
    else:
        _lo, _hi = cf.get_clim()
        cb.set_ticks([t for t in cb.get_ticks() if _lo - 1e-9 <= t <= _hi + 1e-9])
    cb.set_label(cbar_label)
    ax.set_aspect("equal"); ax.grid(False)
    ax.set_xlim(xu.min(), xu.max()); ax.set_ylim(yu.min(), yu.max())
    ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]")
    if case_key is None:
        ax.set_title(title, pad=10)
    else:
        case_title(ax, title, case_key)
    if note:
        # constrained_layout does not reserve space for figure-level text, so
        # shrink the layout rect first -- otherwise the note lands on top of the
        # x-axis label, which is exactly what it did on first writing. Same fix
        # as the figure legend in validate_nasa_real.py.
        fig.get_layout_engine().set(rect=(0, 0.045, 1, 0.955))
        fig.text(0.5, 0.008, note, ha="center", va="bottom", fontsize=8,
                 color=INK_SOFT)
    save(fig, fname)

# ============================================================ 6. CONTOURS
for ff in field_files:
    base = Path(ff).stem            # field_A_validation_peak_a20
    parts = base.split("_")
    cs = parts[1] + "_" + parts[2]
    tag = parts[3]; adeg = parts[4]
    c = CASES[cs]["c"]
    xu, yu, F = load_field(ff)
    pre = f"{cs}_{tag}_{adeg}"
    U = CASES[cs]["U"]
    # pressure -- one scale for all eight fields (set above), so the panels can
    # be compared with each other; the colorbar extends at both ends and nothing
    # is clipped in the data.
    contour_plot(xu, yu, F["Cp"], f"Pressure coefficient $C_p$ — {tag} ({adeg.replace('a','α=')}°)",
                 "$C_p$", CMAP_CP, c, f"contour_Cp_{pre}.png", lines=True, vlim=CP_SCALE,
                 case_key=cs)
    # velocity magnitude + streamlines
    contour_plot(xu, yu, F["speed_ms"], f"Velocity magnitude + streamlines — {tag}",
                 "|V| [m/s]", CMAP_PRESSURE, c, f"contour_speed_stream_{pre}.png",
                 stream=(F["u_ms"], F["v_ms"]), vlim=(0, 1.7*U), case_key=cs)
    # velocity vectors
    contour_plot(xu, yu, F["speed_ms"], f"Velocity vector field — {tag}",
                 "|V| [m/s]", CMAP_PRESSURE, c, f"contour_vectors_{pre}.png",
                 vector=(F["u_ms"], F["v_ms"]), vlim=(0, 1.7*U), case_key=cs)
    # vorticity (DSV). The scale is a robust percentile, and it has to be: the
    # dynamic-stall vortex's core is far more intense than the bound sheet
    # around the body, so an unclipped scale renders everything except the
    # vortex as one flat colour. The clip is now STATED on the figure, the same
    # standard the 3-D surfaces are held to -- otherwise the vortex reads as a
    # small feature at the colorbar's value when it is orders of magnitude past
    # the end of it.
    _w = F["vorticity_1s"]
    # A LINEAR scale cannot show this field. The bound sheet on the body sits at
    # tens of 1/s while the dynamic-stall vortex core reaches ~1e5, so clipping
    # to the 98th percentile (33 1/s) rendered the vortex -- the feature the
    # figure is named after -- as a four-pixel dot, with a note admitting the
    # core was 2600x off the end of the bar. A symmetric log scale shows both:
    # the sheet keeps its structure and the core is on the same bar rather than
    # past it. The linear threshold is the 90th percentile of |w|, so the noise
    # floor stays linear and only the real structure is logarithmic.
    _wa = np.abs(_w[np.isfinite(_w)])
    _lin = max(float(np.nanpercentile(_wa, 90)), 1e-3)
    _top = float(np.nanmax(_wa))
    _dec = np.log10(_top/_lin)
    _lv = np.concatenate([-_lin*np.logspace(_dec, 0, 9), [0.0],
                          _lin*np.logspace(0, _dec, 9)])
    contour_plot(xu, yu, _w,
                 f"Vorticity (dynamic-stall vortex) — {tag}",
                 "ω_z [1/s]", CMAP_VORT, c, f"contour_vorticity_{pre}.png",
                 levels=_lv, norm=_SymLogNorm(linthresh=_lin, vmin=-_top, vmax=_top),
                 cticks=_symlog_ticks(_lin, _top),
                 note="symmetric-log scale: linear within ±%.2g 1/s (the bound sheet), "
                      "logarithmic beyond, to the vortex core at %.3g 1/s. Nothing is clipped."
                      % (_lin, np.nanmin(_w)), case_key=cs)
    # Local Mach and the two temperature fields have long thin tails at the
    # vortex core. Auto-scaling to the full range put ~99 % of the domain into
    # one or two colour bands, so those three plots came out essentially blank;
    # clip the scale to robust percentiles instead (colorbars still extend).
    def rlim(Z, lo=1.5, hi=98.5):
        a = np.nanpercentile(Z, [lo, hi])
        return (float(a[0]), float(a[1])) if a[1] > a[0] else None
    contour_plot(xu, yu, F["Mach_local"], f"Local Mach number — {tag}",
                 "$M_{local}$", CMAP_PRESSURE, c, f"contour_Mach_{pre}.png",
                 lines=True, vlim=rlim(F["Mach_local"]), case_key=cs)
    contour_plot(xu, yu, F["T_static_K"], f"Static air temperature — {tag}",
                 "T [K]", CMAP_TEMP, c, f"contour_Tstatic_{pre}.png",
                 lines=True, vlim=rlim(F["T_static_K"]), case_key=cs)
    contour_plot(xu, yu, F["T_recovery_K"], f"Recovery (skin) temperature — {tag}",
                 "$T_r$ [K]", CMAP_TEMP, c, f"contour_Trecovery_{pre}.png",
                 lines=True, vlim=rlim(F["T_recovery_K"]), case_key=cs)

# ============================================================ 7. TEMPERATURE PROFILE
# surface recovery temperature vs x/c at the 'peak' phase for each case
for cs in CASES:
    pk = [f for f in field_files if cs in f and "peak" in f]
    if not pk: continue
    xu, yu, F = load_field(pk[0]); c = CASES[cs]["c"]
    # sample recovery T just above & below surface along x
    from scipy.interpolate import RegularGridInterpolator
    # interior of the airfoil is NaN by construction; leave it NaN so the
    # profile breaks rather than plotting a fabricated value there
    itp = RegularGridInterpolator((yu, xu), F["T_recovery_K"],
                                  bounds_error=False, fill_value=np.nan)
    xq = np.linspace(0.02*c, 0.98*c, 120)
    # UNITS. yt comes from y_over_c and is NON-dimensional; xq, the offset and
    # the interpolator grid are all in metres. Adding the two directly put the
    # sample point 0.17 c above the surface instead of 0.03 c -- a factor 5.7 --
    # so this "surface" profile was reading the field well outside the boundary
    # layer while still looking perfectly smooth.
    yt = np.interp(xq/c, AF["x_over_c"][:len(AF)//2][::-1],
                   AF["y_over_c"][:len(AF)//2][::-1])*c        # -> metres
    # STANDOFF, derived from the grid rather than chosen. The reconstruction
    # blanks the body and the one ring of cells touching it, and a linear
    # interpolation returns NaN if EITHER bracketing row is blank -- so a sample
    # must clear the surface by three rows, not by a round number of chords. At
    # the flat 0.03c it used to use (2.1 rows here) 12 of the 240 samples fell in
    # that ring and the published profile broke in two places with nothing on the
    # figure to say why. The 0.03c intent is kept as a floor; on this grid the
    # three-row bound is what binds, and the title prints whichever applied.
    _dy = float(yu[1] - yu[0])
    OFF = max(0.03*c, 3.0*_dy)                                 # standoff [m]
    fig, ax = plt.subplots(figsize=(7.8, 4.6))
    for sgn, lab, col in [(+1.0, "upper surface", PALETTE[1]),
                          (-1.0, "lower surface", PALETTE[0])]:
        Tq = itp(np.column_stack([sgn*(yt + OFF), xq]))
        ax.plot(xq/c, Tq, color=col, lw=2, label=lab)
    ax.set_xlabel("x/c"); ax.set_ylabel("recovery temperature  $T_r$ [K]")
    case_title(ax, "Recovery (skin) temperature %.1f%%c off the surface — peak incidence"
               % (100*OFF/c), cs)
    ax.legend(loc="best")
    save(fig, f"temperature_profile_{cs}.png")

# count what THIS script wrote, not everything sitting in plots/ (which also
# holds the fig3d_* figures written by make_3d_plots.py)
print("[plots] 2D figures done:",
      len([f for f in OUT.glob("*.png") if not f.name.startswith("fig3d_")]))
