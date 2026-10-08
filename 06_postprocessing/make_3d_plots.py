# -*- coding: utf-8 -*-
# Run from the repository root:  PYTHONPATH=. python3 06_postprocessing/make_3d_plots.py
"""
06_postprocessing / make_3d_plots.py
------------------------------------
The figures of the case study that show a quantity over two variables,
written to 06_postprocessing/plots/. Run 04_solver/run_case.py first.

Author: Akosa Samuel Onyejekwe (independent)

  fig3d_response_surface.png      peak lift and cycle damping of the load
                                  model over mean incidence and reduced
                                  frequency (05_solution/response_surface.csv),
                                  as surfaces. The damping is model output
                                  only; the figure prints how often the model
                                  gets the sign of the damping right on the
                                  held-out loops (results/validation_damping.csv).
  cp_phase_map_<case>.png         reconstructed upper-surface Cp over the
                                  whole chord and the whole cycle, as a map
  fig3d_cp_phase_surface.png      the same for Case A as a surface
  fig3d_field_surface_Cp.png      reconstructed Cp field as a surface, Case A
                                  at greatest lift
  fig3d_field_surface_speed.png   the same for the speed
  fig3d_section_vectors.png       the section drawn with span, with the
                                  reconstructed velocity on its near end plane

Standing: all but the first show the potential-flow reconstruction of
unistall.flowfield, an illustration: not a flow solution, incompressible with
a linearised (Karman-Tsien) correction of the pressure, quasi-steady,
compared with no measurement. In every figure of the case study the pressure
coefficient is shown with suction upwards (on an axis or on a colour bar),
and a surface of the corrected pressure is cut off at the critical value,
beyond which it is not physical.
"""
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import AsinhNorm
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

import project_meta as pm
from unistall import dsmodel, flowfield
from unistall.style import CMAP_CP, CMAP_PRESSURE, INK, INK_SOFT, PALETTE, apply_style

ROOT = Path(__file__).resolve().parent.parent
HERE = Path(__file__).resolve().parent
SOL = ROOT/"05_solution"
GEO = ROOT/"01_geometry"/"naca0012_coordinates.csv"
DAMPING = ROOT/"results"/"validation_damping.csv"
OUT = HERE/"plots"
PREFIXES = ("fig3d_", "cp_phase_map_")   # the figures this script owns

CASE = pm.CASE_A                         # the case of the surface and vector figures
N_PHASES = 121                           # instants of the cycle in the Cp map
SURFACE_PHASES = 41                      # instants of the cycle in the Cp surface
SPAN_OVER_C = 0.75                       # span drawn in the pictorial view, in chords
VECTOR_STRIDE = 9                        # every n-th grid node carries an arrow
VECTOR_WINDOW = (-0.35, 1.55, 0.55)      # x0, x1 and half-height of the arrows' window, in chords
ARROW_OVER_C = 0.16                      # length of the longest arrow, in chords
NOTE_WIDTH = 118                         # characters to a line of a note
NOTE_PT = pm.NOTE_PT
ALL_LOOPS = "all primary held-out loops with a closed moment loop"
RECORD = []                              # (file, dpi, smallest lettering) of every figure written


def save(fig: plt.Figure, name: str, **kw: object) -> None:
    """Write a figure, enter it in the figure record and close it."""
    RECORD.append((name, float(plt.rcParams["savefig.dpi"]), pm.smallest_lettering(fig)))
    fig.savefig(OUT/name, **kw)
    plt.close(fig)


def axes3d(fig: plt.Figure, rect: list) -> plt.Axes:
    """A 3-D axes with light panes and room for its labels."""
    ax = fig.add_axes(rect, projection="3d")
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.set_edgecolor(INK_SOFT)
        axis.pane.set_alpha(0.04)
        axis.line.set_color(INK_SOFT)
        axis.labelpad = 9
    ax.tick_params(pad=2, labelsize=9)
    return ax


# --------------------------------------------------------------------------- #
def damping_statement(df: pd.DataFrame) -> str:
    """How the damping of the surface stands against the held-out loops."""
    d = pd.read_csv(DAMPING).set_index("group").loc[ALL_LOOPS]
    return ("Cycle damping is model output only. It is negative at %d of these %d conditions. On the %d held-out "
            "measured loops with a closed moment loop the model gets the SIGN of the damping right on a fraction "
            "%.2f (a prediction of positive damping everywhere would score %.2f); %d of those loops measure negative "
            "damping. The damping surface is no evidence of what a measurement would show."
            % (int((df["cycle_damping_Xi"] < 0.0).sum()), len(df), int(d["n_frames"]), d["sign_agreement"],
               d["sign_agreement_if_always_positive"], int(d["measured_negative"])))


def response_surface() -> None:
    """Peak lift and cycle damping over mean incidence and reduced frequency."""
    df = pd.read_csv(SOL/"response_surface.csv")
    grid = {q: df.pivot(index="alpha_mean_deg", columns="reduced_freq_k", values=q)
            for q in ("CL_max", "cycle_damping_Xi")}
    A, K = np.meshgrid(grid["CL_max"].index, grid["CL_max"].columns, indexing="ij")
    fig = plt.figure(figsize=(9.2, 5.6))
    for i, (q, label) in enumerate((("CL_max", "peak lift  $C_{L,max}$"),
                                    ("cycle_damping_Xi", "cycle damping  Ξ  (model output only)"))):
        ax = axes3d(fig, [0.0 + 0.5*i, 0.24, 0.44, 0.62])
        ax.plot_surface(A, K, grid[q].to_numpy(), cmap=pm.trimmed_colormap(CMAP_PRESSURE), edgecolor=INK_SOFT,
                        linewidth=0.3, alpha=0.95)
        ax.set_xlabel("mean incidence [deg]")
        ax.set_ylabel("reduced frequency k")
        ax.set_zlabel(label)
        ax.view_init(elev=24, azim=-60)
    fig.text(0.5, 0.975, "Response of the load model to mean incidence and reduced frequency", ha="center", va="top",
             fontsize=12.5, fontweight="bold", color=INK)
    fig.text(0.5, 0.925, "NACA %s, M = %.2f, pitch amplitude %.0f° (the Mach number and amplitude of Case B); a "
             "steady stream and prescribed pitch, not a rotor calculation"
             % (pm.SECTION, df["mach_M"].iloc[0], df["alpha_amp_deg"].iloc[0]), ha="center", va="top",
             fontsize=NOTE_PT, color=INK_SOFT)
    fig.text(0.5, 0.02, "%d conditions. The highest peak incidence is %.0f°; the measured loops the model was "
             "compared with reach %.0f°.\n" % (len(df), df["peak_alpha_deg"].max(), dsmodel.PEAK_ALPHA_MAX_COMPARED)
             + textwrap.fill(damping_statement(df), 132), ha="center", va="bottom", fontsize=NOTE_PT, color=INK_SOFT)
    save(fig, "fig3d_response_surface.png", bbox_inches=fig.bbox_inches)


def surface_through_cycle(name: str, case: dict, air: flowfield.Air, n_panel: int, n_phases: int) -> tuple:
    """(x/c of the upper-surface control points, phase [deg], corrected
    upper-surface Cp as an array of shape (phases, stations)), from the
    stored time history."""
    th = pd.read_csv(SOL/f"time_history_{name}.csv")
    rows = th.iloc[np.unique(np.linspace(0, len(th) - 1, n_phases).astype(int))]
    Z, x_up = [], None
    for r in rows.itertuples():
        rec = flowfield.solve(GEO, case["chord"], case["U"], r.alpha_deg, r.CL, r.CN_vortex,
                              r.tau_v_semichords/case["Tvl"], r.f_separation, n_panel=n_panel, mach=case["M"],
                              air=air)
        st = flowfield.surface_state(rec, case["M"], air)
        order = np.argsort(st["x_over_c"][st["upper"]])
        x_up = st["x_over_c"][st["upper"]][order]
        Z.append(st["Cp"][st["upper"]][order])
    return x_up, rows["phase_deg"].to_numpy(), np.array(Z)


def notes(*paragraphs: str) -> str:
    """Paragraphs wrapped for the note under a figure."""
    return "\n".join(textwrap.fill(p, NOTE_WIDTH) for p in paragraphs if p)


def finish(fig: plt.Figure, title: str, subtitle: str, note: str, name: str) -> None:
    """Title, case line and note of a 3-D figure, then write it."""
    fig.text(0.5, 0.975, title, ha="center", va="top", fontsize=12.5, fontweight="bold", color=INK)
    fig.text(0.5, 0.925, subtitle, ha="center", va="top", fontsize=NOTE_PT, color=INK_SOFT)
    fig.text(0.5, 0.015, note, ha="center", va="bottom", fontsize=NOTE_PT, color=INK_SOFT)
    save(fig, name, bbox_inches=fig.bbox_inches)


def cp_phase_map(name: str, case: dict, air: flowfield.Air, n_panel: int) -> None:
    """Reconstructed upper-surface Cp against chordwise position and cycle
    phase, over the whole chord, as a map. The region beyond the critical
    pressure is hatched."""
    x_up, phase, Z = surface_through_cycle(name, case, air, n_panel, N_PHASES)
    cp_star = case["critical"]["Cp_critical"]
    ok = Z >= cp_star
    lo, hi = float(Z[ok].min()), float(Z[ok].max())
    fig, ax = plt.subplots(figsize=(7.0, 6.4), constrained_layout=True)
    ax.grid(False)
    shown = ax.pcolormesh(x_up, phase, Z, cmap=CMAP_CP, shading="gouraud", rasterized=True,
                          norm=AsinhNorm(1.0, vmin=lo, vmax=hi, clip=True))
    if not ok.all():
        ax.contourf(x_up, phase, (~ok).astype(float), levels=[0.5, 1.5], colors=["white"], alpha=0.75,
                    hatches=["////"])
        ax.contour(x_up, phase, (~ok).astype(float), levels=[0.5], colors=[INK], linewidths=1.2,
                   linestyles=[(0, (1, 1))])
    ax.set_xscale("function", functions=(np.sqrt, np.square))
    ax.set_xlim(0.0, 1.0)
    ax.set_xticks([0.0, 0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_xticklabels(["0", "0.01", "0.05", "0.1", "0.2", "0.4", "0.6", "0.8", "1"], fontsize=9)
    ax.set_xlabel("x/c  (square-root scale, to open up the leading edge)")
    ax.set_ylabel("cycle phase  ωt  [deg]")
    bar = fig.colorbar(shown, ax=ax, pad=0.02)
    pm.cp_colour_bar(bar, lo, hi)
    bar.set_label("upper-surface $C_p$ (Kármán-Tsien corrected; suction upwards)")
    bar.ax.tick_params(labelsize=9)
    fig.suptitle("Reconstructed upper-surface pressure through the cycle", fontweight="bold", color=INK,
                 fontsize=12.5)
    ax.set_title(case["line"], fontsize=NOTE_PT, color=INK_SOFT, fontweight="normal", pad=8)
    fig.supxlabel(notes(
        "%d chordwise stations at %d instants; the colour scale spans the values outside the hatched region, "
        "%.1f to %.2f. %s Critical $C_p$ %.2f; the region is %.1f %% of the points shown."
        % (Z.shape[1], Z.shape[0], lo, hi, flowfield.CRITICAL_NOTE, cp_star, 100.0*float((~ok).mean())),
        flowfield.STANDING_NOTE), fontsize=NOTE_PT, color=INK_SOFT)
    save(fig, f"cp_phase_map_{name}.png")


def cp_phase_surface(case: dict, air: flowfield.Air, n_panel: int) -> None:
    """Reconstructed upper-surface Cp against chordwise position and cycle
    phase as a surface, suction upwards, over the whole chord, cut off at
    the critical value."""
    x_up, phase, Z = surface_through_cycle(CASE, case, air, n_panel, SURFACE_PHASES)
    cp_star = case["critical"]["Cp_critical"]
    cut = float((Z < cp_star).mean())
    XC, PH = np.meshgrid(np.sqrt(x_up), phase)
    fig = plt.figure(figsize=(7.4, 7.4))
    ax = axes3d(fig, [0.02, 0.29, 0.9, 0.62])
    ax.plot_surface(XC, PH, np.maximum(Z, cp_star), cmap=CMAP_CP, edgecolor=INK_SOFT, linewidth=0.15, alpha=0.97,
                    rstride=1, cstride=4)
    ticks = np.array([0.0, 0.05, 0.2, 0.5, 1.0])
    ax.set_xticks(np.sqrt(ticks))
    ax.set_xticklabels(["%g" % t for t in ticks])
    ax.set_xlabel("x/c  (square-root scale)")
    ax.set_ylabel("cycle phase ωt [deg]")
    ax.set_zlabel("upper-surface $C_p$ (suction upwards)")
    ax.set_zlim(1.0, cp_star)
    ax.view_init(elev=26, azim=-52)
    finish(fig, "Reconstructed upper-surface pressure through the cycle", case["line"],
           notes("%d chordwise stations over the whole chord at %d instants; Kármán-Tsien corrected $C_p$. The "
                 "surface is cut off at the critical value %.2f, beyond which the correction does not hold and "
                 "the values are not physical (%.1f %% of the points shown)."
                 % (Z.shape[1], Z.shape[0], cp_star, 100.0*cut), flowfield.STANDING_NOTE),
           "fig3d_cp_phase_surface.png")


def field_surfaces(case: dict, F: dict, alpha: float) -> None:
    """The corrected Cp field and the incompressible speed field at greatest
    lift, each as a surface over the plane. The Cp axis has suction upwards
    and the surface is cut off at the critical value; the speed is shown
    over its whole range."""
    cp_star = case["critical"]["Cp_critical"]
    beyond = int((F["outside"] & (F["beyond_critical"] > 0.5)).sum())
    for col, cmap, label, tag, flip, text in (
            ("Cp", CMAP_CP, "$C_p$ (suction upwards)", "Cp", True,
             "Kármán-Tsien corrected $C_p$, cut off at the critical value %.2f, beyond which the correction does not "
             "hold and the values are not physical (%d of %d nodes)." % (cp_star, beyond, int(F["outside"].sum()))),
            ("speed_over_U", pm.trimmed_colormap(CMAP_PRESSURE), "|V| / U", "speed", False,
             "Speed of the incompressible solution, over its whole range.")):
        Z = np.where(F["outside"], F[col], np.nan)
        if flip:
            Z = np.maximum(Z, cp_star)
        lo, hi = float(np.nanmin(Z)), float(np.nanmax(Z))
        fig = plt.figure(figsize=(7.4, 7.4))
        ax = axes3d(fig, [0.02, 0.29, 0.9, 0.62])
        ax.plot_surface(F["x"], F["y"], Z, cmap=cmap, vmin=lo, vmax=hi, linewidth=0, antialiased=True,
                        rstride=1, cstride=1)
        ax.set_zlim((hi, lo) if flip else (lo, hi))
        ax.set_xlabel("x/c")
        ax.set_ylabel("y/c")
        ax.set_zlabel(label)
        ax.view_init(elev=40, azim=-58)
        finish(fig, "Reconstructed field as a surface: %s at greatest lift, α = %.1f°" % (label.split(" (")[0], alpha),
               case["line"],
               notes(text + " Range shown: %.3g to %.3g. The gap is the section." % (lo, hi) + " " + flowfield.AXES_NOTE,
                     flowfield.STANDING_NOTE), f"fig3d_field_surface_{tag}.png")


def section_vectors(case: dict, F: dict, alpha: float) -> None:
    """The section drawn with span and the reconstructed velocity on its
    near end plane. The arrows lie on the plane nearest the observer, so
    none passes behind the section, and they are drawn last, on top."""
    sec = pd.read_csv(GEO)
    sx, sy = sec["x_over_c"].to_numpy(), sec["y_over_c"].to_numpy()
    fig = plt.figure(figsize=(7.4, 7.4))
    ax = axes3d(fig, [0.02, 0.29, 0.9, 0.62])
    ax.computed_zorder = False
    for z, order in ((SPAN_OVER_C, 1), (0.0, 3)):
        ax.add_collection3d(Poly3DCollection([list(zip(sx, np.full_like(sx, z), sy, strict=True))],
                                             facecolor="#c6d7ea", edgecolor=INK, linewidths=1.4, zorder=order))
    for i in np.linspace(0, len(sx) - 1, 12).astype(int):
        ax.plot([sx[i], sx[i]], [0.0, SPAN_OVER_C], [sy[i], sy[i]], color=INK_SOFT, lw=0.5, alpha=0.6, zorder=2)
    sk = (slice(None, None, VECTOR_STRIDE), slice(None, None, VECTOR_STRIDE))
    x, y, u, v = F["x"][sk], F["y"][sk], F["u_ms"][sk]/case["U"], F["v_ms"][sk]/case["U"]
    keep = F["outside"][sk] & (x > VECTOR_WINDOW[0]) & (x < VECTOR_WINDOW[1]) & (np.abs(y) < VECTOR_WINDOW[2])
    speed = np.hypot(u[keep], v[keep])
    cap = float(np.percentile(speed, 95))
    shrink = np.minimum(1.0, cap/np.maximum(speed, 1e-9))*ARROW_OVER_C/cap
    ax.quiver(x[keep], np.zeros(keep.sum()), y[keep], u[keep]*shrink, np.zeros(keep.sum()), v[keep]*shrink,
              color=PALETTE[0], linewidth=0.9, arrow_length_ratio=0.32, zorder=4)
    ax.set_xlim(VECTOR_WINDOW[0], VECTOR_WINDOW[1])
    ax.set_ylim(0.0, SPAN_OVER_C)
    ax.set_zlim(-VECTOR_WINDOW[2], VECTOR_WINDOW[2])
    ax.set_xlabel("x/c")
    ax.set_ylabel("span / c")
    ax.set_zlabel("y/c")
    ax.view_init(elev=22, azim=-62)
    finish(fig, "Section with the reconstructed velocity at greatest lift, α = %.1f°" % alpha, case["line"],
           notes("Arrows: incompressible velocity on the end plane nearest the observer, in front of the section; "
                 "the longest stands for %.2f U and faster ones are drawn at that length. The span is drawn for the "
                 "picture only: the calculation is two-dimensional. " % cap + flowfield.AXES_NOTE,
                 flowfield.STANDING_NOTE), "fig3d_section_vectors.png")


def main() -> None:
    """Draw the response surface, the Cp map of each case, and the surface
    and vector figures of Case A."""
    apply_style()
    plt.rcParams["hatch.color"] = INK_SOFT
    plt.rcParams["hatch.linewidth"] = 0.6
    OUT.mkdir(exist_ok=True)
    for old in OUT.glob("*.png"):
        if old.name.startswith(PREFIXES):
            old.unlink()
    cases, air, cfg = pm.read_setup()
    n_panel = cfg["field_reconstruction"]["read_by_the_driver"]["n_panels"]
    response_surface()
    for name, case in cases.items():
        cp_phase_map(name, case, air, n_panel)
    case = cases[CASE]
    cp_phase_surface(case, air, n_panel)
    cp = pd.read_csv(SOL/f"cp_distribution_{CASE}.csv")
    alpha = float(cp[cp["phase"] == "peak"]["alpha_deg"].iloc[0])
    F = flowfield.read_field(next(SOL.glob(f"field_{CASE}_peak_a*.csv.gz")), case["chord"], case["U"])
    field_surfaces(case, F, alpha)
    section_vectors(case, F, alpha)
    pm.record_figures(HERE/"figure_record.csv", RECORD)
    print("[3d] %d figures written to %s" % (len(RECORD), OUT.relative_to(ROOT)))


if __name__ == "__main__":
    main()
