# -*- coding: utf-8 -*-
# Run from the repository root:  PYTHONPATH=. python3 06_postprocessing/make_all_plots.py
"""
06_postprocessing / make_all_plots.py
-------------------------------------
The two-dimensional figures of the case study, drawn from the files in
05_solution/ and 03_model_setup/ and written to 06_postprocessing/plots/.

Author: Akosa Samuel Onyejekwe (independent)

Standing: the load loops, time histories and states show the prediction of
the load model. The surface pressures and every field map show the
potential-flow reconstruction of unistall.flowfield, an illustration: not a
flow solution, incompressible with a linearised (Karman-Tsien) correction of
the pressure, quasi-steady, compared with no measurement. Local Mach number
and temperatures follow from the corrected pressure by the isentropic
relations; where the corrected pressure is below the critical value they are
not physical, and the figures hatch that region and print its size.

<case> is a key of project_meta.CASES; <phase> is rise, peak, dsv or fall
and <deg> the incidence in whole degrees, as in the field files.
  hyst_cl_<case>.png, hyst_cm_<case>.png, hyst_cd_<case>.png   load loops
  timehist_loads_<case>.png          loads against time over one cycle
  states_<case>.png                  separation point, vortex force and clock
  static_inputs.png                  static curves the solver reads
  convergence_residuals.png          change of the cycle with cycles marched
  timestep_refinement.png            loads against steps per cycle
  cp_distribution_<case>.png         reconstructed surface pressure
  temperature_profile_<case>.png     static and recovery temperature along
                                     the surface
  contour_<quantity>_<case>_<phase>_a<deg>.png   reconstructed field, with
      <quantity> one of Cp, Mach, Tstatic, Trecovery, speed_stream,
      vorticity, vectors

Also written, beside this script:
  overlay_record.csv   for each case, whether the measured loop was drawn
  masked_region_record.csv   for each figure of a corrected quantity, the
                       part of the surface and of the field beyond the
                       critical pressure that its note prints
  figure_record.csv    resolution and smallest lettering of every figure

Every colour scale spans the values it colours and says so. The field maps
are drawn in axes fixed to the section. They are written with a 256-colour
palette to keep the files small.
"""
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import AsinhNorm, Normalize, PowerNorm
from matplotlib.patches import Polygon
from PIL import Image
from scipy.ndimage import distance_transform_edt

import project_meta as pm
from unistall import flowfield, metrics, naca4
from unistall.style import CMAP_CP, CMAP_PRESSURE, CMAP_TEMP, CMAP_VORT, INK, INK_SOFT, PALETTE, apply_style

ROOT = Path(__file__).resolve().parent.parent
HERE = Path(__file__).resolve().parent
SOL, SETUP = ROOT/"05_solution", ROOT/"03_model_setup"
OUT = HERE/"plots"
SECTION = pd.read_csv(ROOT/"01_geometry"/"naca0012_coordinates.csv")
OTHER_SCRIPT_PREFIXES = ("fig3d_", "cp_phase_map_")     # figures of make_3d_plots.py, left alone here

BODY_FILL = "#e3e9f0"
MEASURED_LABEL = "measured, NASA TM-84245 frame %s, a calibration loop,\nat its own conditions: M = %.3f, k = %.3f, α = %.1f° ± %.1f°"
MODEL_LABEL = "%s (load model) at the nominal condition of the case:\nM = %.2f, k = %.3f, α = %.0f° ± %.0f°"
PHASE_TEXT = {"mean_up": "mean incidence, upstroke", "rise": "rise (upstroke)", "peak": "greatest lift",
              "dsv": "greatest vortex lift", "fall": "fall (downstroke)", "mean_down": "mean incidence, downstroke"}
FIELD_PHASES = ("rise", "peak", "dsv", "fall")
WINDOW = (-0.75, 1.75, -0.85, 0.85)      # part of the field drawn, in chords
VECTOR_STRIDE = 6                        # every n-th grid node carries an arrow
RESIDUAL_FLOOR = 1e-11                   # where a change of exactly zero is drawn on the log axis
MAP_DPI = 150                            # field maps: 7.4 in wide
MAP_COLOURS = 256                        # size of the palette the field maps are stored with
NOTE_PT = pm.NOTE_PT                     # lettering of the notes under a figure
NOTE_WIDTH = 118                         # characters to a line of a note under a field map
CRITICAL_COLOUR = PALETTE[5]             # the line of the critical pressure coefficient
MARK_COLOUR = PALETTE[3]                 # the separation point of the load model


SPEED_MAP, TEMP_MAP = pm.trimmed_colormap(CMAP_PRESSURE), pm.trimmed_colormap(CMAP_TEMP)

# quantity -> (file tag, column, title, colour-bar label, colour map, kind)
QUANTITIES = [
    ("Cp", "Cp", "Pressure coefficient", "$C_p$ (Kármán-Tsien corrected)", CMAP_CP, "fill"),
    ("Mach", "Mach", "Local Mach number", "local Mach number  [-]", SPEED_MAP, "fill"),
    ("Tstatic", "T_static_K", "Static temperature", "static temperature  [K]", TEMP_MAP, "fill"),
    ("Trecovery", "T_recovery_K", "Recovery temperature", "recovery temperature  [K]", TEMP_MAP, "fill"),
    ("speed_stream", "speed_over_U", "Speed and streamlines", "|V| / U  (incompressible)", SPEED_MAP, "stream"),
    ("vorticity", "vorticity_c_over_U", "Vorticity", r"$\omega c / U$  (negative = clockwise)", CMAP_VORT,
     "vorticity"),
    ("vectors", "speed_over_U", "Velocity vectors", "|V| / U  (incompressible)", SPEED_MAP, "vectors"),
]
CORRECTED = ("Cp", "Mach", "T_static_K", "T_recovery_K")     # columns that follow from the corrected pressure
RECORD = []                              # (file, dpi, smallest lettering) of every figure written
MASKED = []                              # rows of masked_region_record.csv


# --------------------------------------------------------------------------- #
#  inputs and shared pieces
# --------------------------------------------------------------------------- #
def load_cases() -> dict:
    """The cases of 03_model_setup with the per-instant table of the
    reconstruction added."""
    cases, air, _ = pm.read_setup()
    return {name: dict(c, air=air, recon=pd.read_csv(SOL/f"reconstruction_{name}.csv").set_index("phase"))
            for name, c in cases.items()}


def heading(fig: plt.Figure, ax: plt.Axes, title: str, case: dict | None = None) -> None:
    """Figure title, with the case line under it."""
    fig.suptitle(title, fontweight="bold", color=INK, fontsize=12.5)
    if case is not None:
        ax.set_title(case["line"], fontsize=NOTE_PT, color=INK_SOFT, fontweight="normal", pad=8)


def footnote(fig: plt.Figure, text: str) -> None:
    """Small note under the figure."""
    fig.supxlabel(text, fontsize=NOTE_PT, color=INK_SOFT)


def wrapped(*paragraphs: str, width: int = NOTE_WIDTH) -> str:
    """Paragraphs wrapped to `width` characters, one under the other."""
    return "\n".join(textwrap.fill(p, width) for p in paragraphs if p)


def legend_below(ax: plt.Axes, ncol: int = 2, drop: float = -0.16) -> None:
    """Legend under the axes, clear of the curves."""
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, drop), ncol=ncol, frameon=True, fontsize=9)


def save(fig: plt.Figure, name: str, palette: bool = False) -> None:
    """Write the figure to the plots folder, enter it in the figure record
    and close it. With `palette` it is written at MAP_DPI and stored with a
    palette of MAP_COLOURS colours."""
    dpi = MAP_DPI if palette else float(plt.rcParams["savefig.dpi"])
    RECORD.append((name, dpi, pm.smallest_lettering(fig)))
    fig.savefig(OUT/name, dpi=dpi)
    if palette:
        with Image.open(OUT/name) as im:
            small = im.convert("RGB").quantize(MAP_COLOURS, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.NONE)
        small.save(OUT/name, optimize=True)
    plt.close(fig)


def stroke_arrows(ax: plt.Axes, x: np.ndarray, y: np.ndarray, n: int = 8) -> None:
    """Arrow heads along a loop, showing the direction it is traced."""
    for i in np.linspace(len(x)*0.04, len(x)*0.96, n).astype(int):
        ax.annotate("", xy=(x[i + 2], y[i + 2]), xytext=(x[i - 2], y[i - 2]),
                    arrowprops=dict(arrowstyle="-|>", color=INK_SOFT, lw=1.1), zorder=6)


# --------------------------------------------------------------------------- #
#  loads
# --------------------------------------------------------------------------- #
def measured_loop(name: str, case: dict) -> tuple:
    """(measured loop or None, row of the overlay record). A case with a
    measured frame whose file is not in the cache is drawn without it; the
    record says so, a warning is printed, and the documents read the record."""
    frame = case["measured_frame"]
    if not frame:
        return None, dict(case=name, measured_frame="", drawn=False, note="the case has no measured loop")
    path = metrics.FRAME_CACHE/(frame + ".mat")
    if not path.exists():
        print(f"[plots] WARNING: {path} is absent (run python3 -m unistall.fetch_frames); the loops of {name} are "
              "drawn WITHOUT the measured points")
        return None, dict(case=name, measured_frame=frame, drawn=False,
                          note="the frame file was not in the cache when the figures were drawn")
    meas = metrics.load_frame(path)
    return meas, dict(case=name, measured_frame=frame, drawn=True, note="drawn", M=round(meas["M"], 4),
                      k=round(meas["k"], 4), alpha_mean_deg=round(meas["a0"], 3), alpha_amp_deg=round(meas["da"], 3))


def hysteresis(name: str, case: dict, th: pd.DataFrame) -> dict:
    """Lift, moment and drag against incidence over the cycle. The legend
    gives the conditions of the curve and of the measured points, which are
    not the same. Returns the row of the overlay record."""
    meas, record = measured_loop(name, case)
    a = th["alpha_deg"].to_numpy()
    for tag, col, label in (("cl", "CL", "lift coefficient  $C_L$"),
                            ("cm", "CM_c4", "moment coefficient about c/4  $C_M$"),
                            ("cd", "CD", "drag coefficient  $C_D$")):
        fig, ax = plt.subplots(figsize=(5.6, 5.4))
        ax.plot(a, th[col], color=PALETTE[0], lw=2.2,
                label=MODEL_LABEL % (case["solver"], case["M"], case["k"], case["a_mean"], case["a_amp"]))
        stroke_arrows(ax, a, th[col].to_numpy())
        if meas is not None:
            ax.plot(meas["a" + tag], meas[tag], "o", color=PALETTE[1], ms=4.5, mfc="none", mew=1.3,
                    label=MEASURED_LABEL % (case["measured_frame"].split("_")[1], meas["M"], meas["k"], meas["a0"],
                                            meas["da"]))
        ax.set_xlabel("incidence  α  [deg]")
        ax.set_ylabel(label)
        heading(fig, ax, "Load loop: " + label.split("  ")[0], case)
        legend_below(ax, ncol=1, drop=-0.17)
        save(fig, f"hyst_{tag}_{name}.png")
    return record


def time_histories(name: str, case: dict, th: pd.DataFrame) -> None:
    """Incidence and the three loads against time over the cycle."""
    t = th["time_s"].to_numpy()*1e3
    fig, axs = plt.subplots(4, 1, figsize=(6.4, 7.6), sharex=True)
    for ax, (col, label, colour) in zip(axs, (("alpha_deg", "α [deg]", PALETTE[3]), ("CL", "$C_L$", PALETTE[0]),
                                              ("CD", "$C_D$", PALETTE[2]), ("CM_c4", "$C_M$ (c/4)", PALETTE[1])),
                                        strict=True):
        ax.plot(t, th[col], color=colour, lw=2.0)
        ax.set_ylabel(label)
    axs[-1].set_xlabel("time  [ms]")
    heading(fig, axs[0], "Loads over one cycle", case)
    save(fig, f"timehist_loads_{name}.png")


def states(name: str, case: dict, th: pd.DataFrame) -> None:
    """The states of the load model against cycle phase."""
    ph = th["phase_deg"]
    cn1 = case["static_inputs"]["CN1_onset"]
    fig, axs = plt.subplots(2, 1, figsize=(6.6, 7.4), sharex=True)
    axs[0].plot(ph, th["f_static"], color=PALETTE[6], lw=1.6, ls="--", label="static separation point  f'")
    axs[0].plot(ph, th["f_separation"], color=PALETTE[4], lw=2.0, label="lagged separation point  f''")
    axs[0].plot(ph, th["CN_vortex"], color=PALETTE[1], lw=2.0, ls="-.", label="vortex normal force  $C_N^v$")
    axs[0].set_ylabel("separation point,  $C_N^v$")
    axs[1].plot(ph, th["CN_prime"]/cn1, color=PALETTE[0], lw=2.0, label="lagged normal force  $C'_N / C_{N1}$")
    axs[1].plot(ph, th["tau_v_semichords"]/case["Tvl"], color=PALETTE[3], lw=2.0, ls="--",
                label="vortex clock  $τ_v / T_{VL}$")
    axs[1].axhline(1.0, color=INK_SOFT, lw=1.0, ls=":", label="onset level and end of the chord (1)")
    axs[1].set_ylabel("ratio")
    axs[1].set_xlabel("cycle phase  ωt  [deg]")
    axs[0].legend(loc="lower center", bbox_to_anchor=(0.5, 1.12), ncol=2, frameon=True, fontsize=9)
    legend_below(axs[1], ncol=2, drop=-0.2)
    fig.suptitle("States of the load model over one cycle", fontweight="bold", color=INK, fontsize=12.5)
    axs[0].set_title(case["line"], fontsize=NOTE_PT, color=INK_SOFT, fontweight="normal", pad=6)
    save(fig, f"states_{name}.png")


def static_inputs(cases: dict) -> None:
    """The static normal force, moment and separation point the solver reads."""
    tab = pd.read_csv(SETUP/"static_inputs.csv")
    fig, axs = plt.subplots(1, 3, figsize=(9.2, 3.9))
    for i, name in enumerate(cases):
        sub = tab[tab["case"] == name]
        label = "M = %.2f  (%s)" % (cases[name]["M"], pm.short_label(name))
        for ax, col in zip(axs, ("CN_static", "CM_static_c4", "f_separation"), strict=True):
            ax.plot(sub["alpha_deg"], sub[col], color=PALETTE[i], lw=2.0, ls=("-", "--")[i % 2], label=label)
    for ax, label in zip(axs, ("static normal force  $C_N$", "static moment about c/4  $C_M$",
                               "separation point  f"), strict=True):
        ax.set_xlabel("incidence  α  [deg]")
        ax.set_ylabel(label)
    legend_below(axs[1], drop=-0.24)
    fig.suptitle("Static inputs of the load model at the Mach numbers of the two cases", fontweight="bold",
                 color=INK, fontsize=12.5)
    save(fig, "static_inputs.png")


# --------------------------------------------------------------------------- #
#  convergence
# --------------------------------------------------------------------------- #
def convergence_residuals(cases: dict) -> None:
    """Largest change of lift and moment over the cycle when one more cycle is marched."""
    fig, ax = plt.subplots(figsize=(6.6, 5.2))
    reported = []
    for i, name in enumerate(cases):
        df = pd.read_csv(SOL/"convergence"/f"residuals_{name}.csv").dropna()
        reported.append("%s: %d" % (pm.short_label(name), cases[name]["cycles"]))
        for col, sym, ls, tone in (("max_abs_change_CL", "o", "-", 0), ("max_abs_change_CM", "s", "--", 1)):
            y = df[col].to_numpy()
            zero = y <= 0.0
            ax.semilogy(df["cycles_marched"], np.where(zero, RESIDUAL_FLOOR, y), ls=ls, marker=sym, ms=6,
                        lw=1.8, color=PALETTE[2*i + tone], mfc="white",
                        label="%s, %s" % (pm.short_label(name), col.split("_")[-1]
                                          .replace("CL", "$C_L$").replace("CM", "$C_M$")))
            ax.semilogy(df["cycles_marched"][~zero], y[~zero], sym, ms=6, color=PALETTE[2*i + tone])
    ax.axhline(RESIDUAL_FLOOR, color=INK_SOFT, lw=0.8, ls=":")
    ax.set_xlabel("cycles marched  n")
    ax.set_ylabel("largest change over the cycle from n − 1 to n")
    fig.suptitle("Cycle-to-cycle change of the loads", fontweight="bold", color=INK, fontsize=12.5)
    legend_below(ax, ncol=2)
    footnote(fig, wrapped("Open symbols on the dotted line: no change at all. Cycles marched for the reported "
                          "results: " + ", ".join(reported) + ". A change that reappears after the loop has settled "
                          "comes from the sign of the pitch rate at a turning point, which rounding can place on "
                          "either side of zero.", width=100))
    save(fig, "convergence_residuals.png")


def timestep_refinement(cases: dict) -> None:
    """Difference of the peak loads from the finest march, against steps per cycle."""
    tab = pd.read_csv(SOL/"convergence"/"timestep_refinement.csv")
    order = pd.read_csv(SOL/"convergence"/"timestep_order.csv")
    fig, axs = plt.subplots(1, len(cases), figsize=(9.2, 4.6), sharey=True)
    for ax, name in zip(axs, cases, strict=True):
        sub = tab[tab["case"] == name].iloc[:-1]            # the finest march is the reference
        for i, (q, label) in enumerate((("CL_max", "$C_{L,max}$"), ("CM_min", "$C_{M,min}$"),
                                        ("CD_max", "$C_{D,max}$"))):
            ax.loglog(sub["steps_per_cycle"], sub["pct_from_finest_" + q], marker="os^"[i], ms=6, lw=1.8,
                      ls=("-", "--", "-.")[i], color=PALETTE[i], label=label)
        used = int(tab[(tab["case"] == name) & tab["reported"]]["steps_per_cycle"].iloc[0])
        ax.axvline(used, color=INK_SOFT, lw=1.1, ls=":", label="steps used (%d)" % used)
        ax.set_xticks(sub["steps_per_cycle"])
        ax.set_xticklabels([str(int(v)) for v in sub["steps_per_cycle"]], fontsize=9)
        ax.set_xticks([], minor=True)
        ax.set_xlabel("steps per cycle")
        ax.set_title(cases[name]["label"], fontsize=10, color=INK, fontweight="normal")
        legend_below(ax, ncol=2, drop=-0.2)
    finest = tab.groupby("case")["steps_per_cycle"].max()
    seen = int(order["order_observed"].sum())
    axs[0].set_ylabel("difference from the finest march  [%]")
    fig.suptitle("Refinement of the time step", fontweight="bold", color=INK, fontsize=12.5)
    footnote(fig, wrapped("Reference: " + ", ".join("%d steps per cycle (%s)" % (finest[n], pm.short_label(n))
                                                    for n in cases)
                          + ". The curves are not straight lines of one slope: an order of convergence is observed "
                          "in %d of the %d load quantities examined (05_solution/convergence/timestep_order.csv)."
                          % (seen, len(order)), width=125))
    save(fig, "timestep_refinement.png")


# --------------------------------------------------------------------------- #
#  reconstructed surface pressure
# --------------------------------------------------------------------------- #
def footprint_sentence(recon: pd.DataFrame) -> str:
    """What the vortex does to the surface pressure beneath it, at the
    instants at which it is over the chord."""
    over = recon[recon["vortex_over_chord"] == 1]
    if over.empty:
        return "The vortex marker is over the chord at none of these instants."
    parts = ["%s %+.2f (x/c = %.2f)" % (PHASE_TEXT[tag], r["dCp_under_vortex"], r["vortex_x_over_c"])
             for tag, r in over.iterrows()]
    raised = int((over["dCp_under_vortex"] > 0.0).sum())
    return ("Change of surface $C_p$ beneath the vortex marker when it is put in: " + "; ".join(parts)
            + ". It RAISES the pressure beneath it at %d of these %d instants: it does not reproduce the suction "
            "footprint of a dynamic-stall vortex." % (raised, len(over)))


def masked_sentence(recon: pd.DataFrame, tags: object) -> str:
    """The part of the surface beyond the critical pressure at each of the
    instants `tags`, in words."""
    parts = ["%s %.1f %%" % (PHASE_TEXT[t], 100.0*recon.loc[t, "beyond_critical_surface_fraction"]) for t in tags]
    return "Part of the surface beyond the critical pressure: " + "; ".join(parts) + "."


def record_masked(fname: str, case_name: str, recon: pd.DataFrame, tags: object) -> None:
    """Enter in the masked-region record what the note of a figure prints."""
    for t in tags:
        r = recon.loc[t]
        MASKED.append(dict(file=fname, case=case_name, phase=t,
                           surface_fraction=r["beyond_critical_surface_fraction"],
                           field_area_over_c2=r.get("beyond_critical_field_area_over_c2", np.nan),
                           field_fraction=r.get("beyond_critical_field_fraction", np.nan)))


def cp_distribution(name: str, case: dict, cp: pd.DataFrame) -> None:
    """Reconstructed surface pressure, corrected for compressibility, at the
    instants written by the solution stage. Both panels have suction
    upwards; the band beyond the critical pressure is hatched."""
    r = case["recon"]
    cp_star = float(r["Cp_critical"].iloc[0])
    tags = list(cp["phase"].unique())
    fig = plt.figure(figsize=(9.2, 8.0))
    gs = fig.add_gridspec(2, 2, height_ratios=(1.0, 0.16))
    axs, low = [fig.add_subplot(gs[0, i]) for i in range(2)], fig.add_subplot(gs[1, :])
    low.axis("off")
    for i, tag in enumerate(tags):
        sub = cp[cp["phase"] == tag]
        label = "%s, α = %.1f°" % (PHASE_TEXT[tag], sub["alpha_deg"].iloc[0])
        for ax in axs:
            for surface, ls in (("upper", "-"), ("lower", "--")):
                q = sub[sub["surface"] == surface].sort_values("x_over_c")
                ax.plot(q["x_over_c"], q["Cp"], color=PALETTE[i % len(PALETTE)], lw=1.7, ls=ls,
                        label=label if surface == "upper" else None)
    floor = float(cp["Cp"].min())
    for ax, title in zip(axs, ("whole range", "detail above the critical value"), strict=True):
        ax.axhspan(cp_star, floor - 1.0, facecolor="none", edgecolor=CRITICAL_COLOUR, hatch="///", lw=0.0, alpha=0.5)
        ax.axhline(cp_star, color=CRITICAL_COLOUR, lw=1.3, ls=(0, (1, 1)),
                   label="critical $C_p^*$ = %.2f; hatched beyond it: not physical" % cp_star)
        ax.set_xlabel("x/c")
        ax.set_ylabel("$C_p$ (Kármán-Tsien corrected; suction upwards)")
        ax.set_title(title, fontsize=10, color=INK, fontweight="normal")
    axs[0].set_ylim(1.2, floor - 0.6)
    axs[1].set_ylim(1.1, cp_star - 0.4)
    low.legend(*axs[0].get_legend_handles_labels(), loc="center", ncol=3, frameon=True, fontsize=9)
    fig.suptitle("Reconstructed surface pressure\n" + case["line"].split("\n")[0], fontweight="bold", color=INK,
                 fontsize=12)
    footnote(fig, wrapped(
        "Solid: upper surface; dashed: lower surface. " + flowfield.STANDING_NOTE,
        "Beyond the dotted line the corrected $C_p$ is below the critical value %.2f (local Mach number above one): "
        "the correction does not hold there and the values are not physical; the corrected value is not allowed "
        "below %.1f (%.0f %% of the vacuum value). " % (cp_star, floor, 100.0*flowfield.VACUUM_CAP)
        + masked_sentence(r, tags) + " Lift returned by this pressure at greatest lift: %+.3f %% from the lift given, "
        "with the sheet carrying %.3f of the circulation of the lift (Prandtl-Glauert factor %.3f); the "
        "incompressible pressure with the whole circulation returns it within %+.2f %%."
        % (r.loc["peak", "closure_error_pct"], r.loc["peak", "circulation_factor"],
           r.loc["peak", "prandtl_glauert_factor"], r.loc["peak", "closure_error_pct_incompressible"]),
        footprint_sentence(r), width=150))
    record_masked(f"cp_distribution_{name}.png", name, r, tags)
    save(fig, f"cp_distribution_{name}.png")


def temperature_profile(name: str, case: dict, cp: pd.DataFrame) -> None:
    """Static and recovery temperature along the surface at the four field
    instants, from the corrected pressure by the isentropic relations. The
    points beyond the critical pressure are left out."""
    r, air = case["recon"], case["air"]
    T0 = air.T_inf*(1.0 + 0.5*(air.gamma - 1.0)*case["M"]**2)
    fig = plt.figure(figsize=(9.2, 7.2))
    gs = fig.add_gridspec(2, 2, height_ratios=(1.0, 0.16))
    axs, low = [fig.add_subplot(gs[0, i]) for i in range(2)], fig.add_subplot(gs[1, :])
    low.axis("off")
    for i, tag in enumerate(FIELD_PHASES):
        sub = cp[cp["phase"] == tag]
        for ax, col in zip(axs, ("T_static_K", "T_recovery_K"), strict=True):
            for surface, ls in (("upper", "-"), ("lower", "--")):
                q = sub[sub["surface"] == surface].sort_values("x_over_c")
                ax.plot(q["x_over_c"], q[col].where(q["beyond_critical"] == 0), color=PALETTE[i], lw=1.7, ls=ls,
                        label="%s, α = %.1f°" % (PHASE_TEXT[tag], sub["alpha_deg"].iloc[0])
                        if surface == "upper" else None)
    for ax, label in zip(axs, ("static temperature  [K]", "recovery temperature  [K]"), strict=True):
        ax.axhline(T0, color=INK_SOFT, lw=1.0, ls=":", label="total temperature %.2f K" % T0)
        ax.axhline(air.T_inf, color=INK_SOFT, lw=1.0, ls="-.", label="stream temperature %.2f K" % air.T_inf)
        ax.set_xlabel("x/c")
        ax.set_ylabel(label)
    low.legend(*axs[0].get_legend_handles_labels(), loc="center", ncol=3, frameon=True, fontsize=9)
    fig.suptitle("Reconstructed surface temperatures\n" + case["line"].split("\n")[0], fontweight="bold", color=INK,
                 fontsize=12)
    footnote(fig, wrapped(
        "Solid: upper surface; dashed: lower surface. From the corrected pressure by the isentropic relations with "
        "the stagnation state of the stream; recovery factor %.3f (%s). The recovery temperature is that of an "
        "adiabatic wall under a turbulent boundary layer, which the reconstruction does not have. The curves stop "
        "where the corrected pressure is below the critical value: the values there are not physical. "
        % (air.recovery, flowfield.RECOVERY_SOURCE) + masked_sentence(r, FIELD_PHASES),
        flowfield.STANDING_NOTE, width=150))
    record_masked(f"temperature_profile_{name}.png", name, r, FIELD_PHASES)
    save(fig, f"temperature_profile_{name}.png")


# --------------------------------------------------------------------------- #
#  reconstructed fields
# --------------------------------------------------------------------------- #
def drawn(F: dict, col: str = "") -> np.ndarray:
    """Nodes outside the section and inside the drawn part of the field; for
    a column that follows from the corrected pressure, those of them where
    that pressure is not below the critical value."""
    keep = F["valid"] if col in CORRECTED else F["outside"]
    return keep & (F["x"] >= WINDOW[0]) & (F["x"] <= WINDOW[1]) & (F["y"] >= WINDOW[2]) & (F["y"] <= WINDOW[3])


def colour_scales(fields: dict, case: dict) -> dict:
    """The colour scale of each quantity, common to the instants of one case,
    as (norm, lowest, highest). Each spans the values it colours: the nodes
    drawn outside the section, over all the instants, and for a quantity
    that follows from the corrected pressure those of them outside the
    hatched region beyond the critical pressure. The scale of Cp is linear
    near zero and logarithmic beyond; that of the speed is a square root;
    those of the Mach number and the temperatures are linear; that of the
    vorticity is linear and symmetric, set by the largest
    core vorticity of the vortex at these instants as its circulation and
    core radius give it, not by a sampled node."""
    scales = {}
    for _, col, _, _, _, kind in QUANTITIES:
        vals = np.concatenate([F[col][drawn(F, col)] for F in fields.values()])
        lo, hi = float(vals.min()), float(vals.max())
        if kind == "vorticity":
            top = float(case["recon"].loc[list(fields), "vortex_peak_vorticity_c_over_U"].max()) or 1.0
            scales[col] = (Normalize(-top, top), -top, top)
        elif col == "Cp":
            scales[col] = (AsinhNorm(1.0, vmin=lo, vmax=hi, clip=True), lo, hi)
        elif col in CORRECTED:
            scales[col] = (Normalize(lo, hi, clip=True), lo, hi)
        else:
            scales[col] = (PowerNorm(0.5, vmin=lo, vmax=hi), lo, hi)
    return scales


def carried_inside(Z: np.ndarray, outside: np.ndarray) -> np.ndarray:
    """Z with each node inside the section given the value of the nearest
    node outside it, so that the colours can be drawn up to the outline."""
    nearest = distance_transform_edt(~outside, return_distances=False, return_indices=True)
    return Z[tuple(nearest)]


def draw_quantity(ax: plt.Axes, F: dict, col: str, cmap: object, kind: str, norm: Normalize) -> object:
    """Draw one quantity on the axes; returns the mappable for the colour
    bar. The colours run up to the outline, which is drawn over them: within
    one grid spacing of it they are carried over from the nearest node
    outside the section (the values on the surface itself are those of the
    surface-pressure figure)."""
    X, Y, Z = F["x"], F["y"], carried_inside(F[col], F["outside"])
    beyond = np.where(F["outside"], F["beyond_critical"], 0.0)
    if col in CORRECTED and (beyond > 0.5).any():
        ax.contourf(X, Y, beyond, levels=[0.5, 1.5], colors=["white"], alpha=0.75, hatches=["////"], zorder=4)
        ax.contour(X, Y, beyond, levels=[0.5], colors=[CRITICAL_COLOUR], linewidths=1.2, linestyles=[(0, (1, 1))],
                   zorder=4)
    if kind == "vectors":
        sk = (slice(None, None, VECTOR_STRIDE), slice(None, None, VECTOR_STRIDE))
        keep = F["outside"][sk]
        speed = np.hypot(F["u_ms"], F["v_ms"])[sk][keep]
        length = 0.8*VECTOR_STRIDE*float(X[0, 1] - X[0, 0])          # arrow of a node moving at U
        return ax.quiver(X[sk][keep], Y[sk][keep], F["u_ms"][sk][keep]/speed, F["v_ms"][sk][keep]/speed, Z[sk][keep],
                         cmap=cmap, norm=norm, angles="xy", scale_units="xy", scale=1.0/length, width=0.0032)
    shown = ax.pcolormesh(X, Y, Z, cmap=cmap, norm=norm, shading="gouraud", rasterized=True)
    if kind == "stream":
        gx = np.linspace(X[0, 0], X[0, -1], X.shape[1])              # exactly uniform, as streamplot needs
        gy = np.linspace(Y[0, 0], Y[-1, 0], Y.shape[0])
        ax.streamplot(gx, gy, np.ma.masked_where(~F["outside"], F["u_ms"]),
                      np.ma.masked_where(~F["outside"], F["v_ms"]), density=1.3, color=INK, linewidth=0.55,
                      arrowsize=0.7)
    return shown


def mark_model(ax: plt.Axes, r: pd.Series) -> None:
    """Mark on the section the separation point of the load model and, where
    there is one in view, the centre of the vortex marker."""
    side = 1.0 if r["transpiration_side"] >= 0 else -1.0
    f = float(np.clip(r["f_separation"], 0.0, 1.0))
    y = side*float(naca4.closed_thickness(np.array([f]), naca4.thickness_ratio(pm.SECTION))[0])
    ax.plot([f], [y], marker="v" if side > 0 else "^", ms=11, mfc=MARK_COLOUR, mec="white", mew=1.2, ls="none",
            zorder=7, label="separation point of the load model, x/c = %.2f" % f)
    if r["vortex_peak_vorticity_c_over_U"] > 0.0 and WINDOW[0] < r["vortex_x_over_c"] < WINDOW[1]:
        ax.plot([r["vortex_x_over_c"]], [r["vortex_y_over_c"]], marker="+", ms=12, mew=2.0, color=INK, ls="none",
                zorder=7, label="centre of the vortex marker (clock at %.2f)" % r["vortex_clock_tau_over_Tvl"])


def field_note(kind: str, col: str, lo: float, hi: float, case: dict, r: pd.Series, n_inst: int) -> str:
    """The note under a field map: what the scale spans, the region beyond
    the critical pressure and its size, the separated region, the axes and
    the standing of the field."""
    cp_star = case["critical"]["Cp_critical"]
    if kind == "vorticity":
        scale = ("Colour scale: linear, symmetric, set by the largest core vorticity of the vortex marker over the %d "
                 "instants of this case (|ωc/U| = %.1f, from its circulation and core radius). The vorticity is zero "
                 "outside that core: the bound sheet lies on the surface itself and the transpiration sheds none; a "
                 "map of the singularities, with no boundary-layer or wake vorticity in it." % (n_inst, hi))
    else:
        scale = ("Colour scale spans the values at the nodes drawn outside the section%s over the %d instants of this "
                 "case: %.2f to %.2f. Within one grid spacing of the outline the colour is carried over from the "
                 "nearest node outside." % (" and outside the hatched region" if col in CORRECTED else "", n_inst,
                                            lo, hi))
    if col in CORRECTED:
        scale += (" " + flowfield.CRITICAL_NOTE + " Critical $C_p$ %.2f. At this instant the region is %.1f %% of the "
                  "surface and %.4f c² of the stored field (%.3f %% of it)."
                  % (cp_star, 100.0*r["beyond_critical_surface_fraction"], r["beyond_critical_field_area_over_c2"],
                     100.0*r["beyond_critical_field_fraction"]))
    if col in ("T_static_K", "T_recovery_K"):
        scale += (" Temperatures from the corrected pressure by the isentropic relations; recovery factor %.3f (%s)."
                  % (case["air"].recovery, flowfield.RECOVERY_SOURCE))
    if col == "T_recovery_K":
        scale += (" A recovery temperature belongs to an adiabatic wall: away from the wall the map shows what such "
                  "a wall would take at the local state.")
    if kind == "vectors":
        scale = "Arrows show direction only; colour shows speed. " + scale
    return wrapped(scale, "The separated region is not drawn: the triangle marks where the load model puts the "
                   "separation point (f = %.2f), and the flow shown aft of it is not separated flow. "
                   % r["f_separation"] + flowfield.AXES_NOTE, flowfield.STANDING_NOTE)


def field_figure(F: dict, quantity: tuple, case: dict, phase: str, alpha: float, scale: tuple, fname: str,
                 n_inst: int) -> None:
    """One map of one quantity at one instant."""
    tag, col, title, label, cmap, kind = quantity
    norm, lo, hi = scale
    r = case["recon"].loc[phase]
    fig, ax = plt.subplots(figsize=(7.4, 7.5))
    ax.set_facecolor("white")
    ax.grid(False)
    mappable = draw_quantity(ax, F, col, cmap, kind, norm)
    ax.add_patch(Polygon(np.column_stack([SECTION["x_over_c"], SECTION["y_over_c"]]), closed=True,
                         facecolor=BODY_FILL, edgecolor=INK, lw=1.2, zorder=5))
    mark_model(ax, r)
    ax.set_xlim(WINDOW[0], WINDOW[1])
    ax.set_ylim(WINDOW[2], WINDOW[3])
    ax.set_aspect("equal")
    ax.set_xlabel("x/c")
    ax.set_ylabel("y/c")
    ax.legend(loc="lower left", fontsize=NOTE_PT, frameon=True, framealpha=0.92)
    bar = fig.colorbar(mappable, ax=ax, pad=0.02, fraction=0.05, shrink=0.9)
    bar.set_label(label)
    if col == "Cp":
        pm.cp_colour_bar(bar, lo, hi)
        bar.set_label(label + ", suction upwards")
    bar.ax.tick_params(labelsize=9)
    heading(fig, ax, "%s at %s, α = %.1f°" % (title, PHASE_TEXT[phase], alpha), case)
    footnote(fig, field_note(kind, col, lo, hi, case, r, n_inst))
    if col in CORRECTED:
        record_masked(fname, "", case["recon"], [phase])
    save(fig, fname, palette=True)


def field_figures(name: str, case: dict, cp: pd.DataFrame) -> int:
    """Every field map of one case; returns how many were written."""
    alpha = cp.groupby("phase")["alpha_deg"].first()
    paths = {ph: next(SOL.glob(f"field_{name}_{ph}_a*.csv.gz")) for ph in FIELD_PHASES}
    fields = {ph: flowfield.read_field(p, case["chord"], case["U"]) for ph, p in paths.items()}
    scales = colour_scales(fields, case)
    for ph, F in fields.items():
        suffix = paths[ph].name[len("field_"):-len(".csv.gz")]
        for q in QUANTITIES:
            field_figure(F, q, case, ph, float(alpha[ph]), scales[q[1]], f"contour_{q[0]}_{suffix}.png", len(fields))
    return len(fields)*len(QUANTITIES)


def main() -> None:
    """Draw every two-dimensional figure."""
    apply_style()
    plt.rcParams["figure.constrained_layout.use"] = True
    plt.rcParams["hatch.color"] = INK_SOFT
    plt.rcParams["hatch.linewidth"] = 0.6
    OUT.mkdir(exist_ok=True)
    for old in OUT.glob("*.png"):
        if not old.name.startswith(OTHER_SCRIPT_PREFIXES):
            old.unlink()
    cases = load_cases()
    n_maps, overlay = 0, []
    for name, case in cases.items():
        th = pd.read_csv(SOL/f"time_history_{name}.csv")
        cp = pd.read_csv(SOL/f"cp_distribution_{name}.csv")
        overlay.append(hysteresis(name, case, th))
        time_histories(name, case, th)
        states(name, case, th)
        cp_distribution(name, case, cp)
        temperature_profile(name, case, cp)
        n_maps += field_figures(name, case, cp)
    static_inputs(cases)
    convergence_residuals(cases)
    timestep_refinement(cases)
    pd.DataFrame(overlay).to_csv(HERE/"overlay_record.csv", index=False)
    pd.DataFrame(MASKED).drop(columns="case").to_csv(HERE/"masked_region_record.csv", index=False)
    pm.record_figures(HERE/"figure_record.csv", RECORD)
    print("[plots] %d field maps and %d other figures written to %s"
          % (n_maps, len(RECORD) - n_maps, OUT.relative_to(ROOT)))


if __name__ == "__main__":
    main()
