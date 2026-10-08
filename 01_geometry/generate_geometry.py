# -*- coding: utf-8 -*-
# Run from the repository root:  PYTHONPATH=. python3 01_geometry/generate_geometry.py
"""
01_geometry / generate_geometry.py
----------------------------------
The NACA 0012 section of the case study: surface coordinates, a table of its
geometric properties and two figures.

Author: Akosa Samuel Onyejekwe (independent)

The outline is the closed-edge form of the four-digit thickness distribution
(unistall.naca4.outline): the thickness is zero at the trailing edge. The
mesh, the panels of the flow-field reconstruction and the drawings use the
same outline.

Writes, beside this file:
  naca0012_coordinates.csv      x/c, y/c from the trailing edge along the upper
                                surface to the leading edge and back along the
                                lower surface to the same trailing-edge point
  section_geometry_summary.csv  form of the trailing edge, thickness, its
                                position, nose radius, trailing-edge gap,
                                enclosed area, number of points
  fig_geometry_profile.png      the profile
  fig_geometry_thickness.png    the thickness distribution
  figure_record.csv             resolution and smallest lettering of the figures
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import brentq

import project_meta as pm
from unistall import naca4
from unistall.style import INK, INK_SOFT, PALETTE, apply_style

HERE = Path(__file__).resolve().parent
CODE = pm.SECTION
THICKNESS_RATIO = naca4.thickness_ratio(CODE)
FILL = "#cfdcec"
COORDINATE_DECIMALS = 9              # so that the file lies on the outline to well within 1e-6 of the chord
RECORD = []                          # (file, dpi, smallest lettering) of each figure written


def summary_table(X: np.ndarray, Y: np.ndarray, yt: np.ndarray) -> pd.DataFrame:
    """Geometric properties of the section; lengths as fractions of chord unless a unit is given."""
    x_tmax = brentq(naca4.closed_thickness_slope, 0.05, 0.9, args=(THICKNESS_RATIO,))
    t_max = 2.0*float(naca4.closed_thickness(x_tmax, THICKNESS_RATIO))
    area = 0.5*abs(float(np.dot(X, np.roll(Y, -1)) - np.dot(Y, np.roll(X, -1))))
    rows = [("section", "NACA " + CODE, "-"),
            ("trailing_edge", "closed: x^4 coefficient %s of the four-digit thickness form" % naca4.A4_CLOSED, "-"),
            ("max_thickness_t_c", round(t_max, 5), "-"),
            ("x_at_max_thickness_x_c", round(x_tmax, 4), "-"),
            ("LE_radius_r_c", round(naca4.NOSE_RADIUS_FACTOR*THICKNESS_RATIO**2, 5), "-"),
            ("TE_gap_percent_chord", round(200.0*float(yt[-1]), 3), "%"),
            ("enclosed_area_A_c2", round(area, 5), "-"),
            ("n_surface_points", len(X), "-")]
    return pd.DataFrame(rows, columns=["property", "value", "units"])


def profile_figure(X: np.ndarray, Y: np.ndarray, path: Path) -> None:
    """The profile, to scale."""
    fig, ax = plt.subplots(figsize=(9, 2.9))
    ax.fill(X, Y, color=FILL)
    ax.plot(X, Y, color=INK, lw=1.8)
    ax.axhline(0.0, color=INK_SOFT, lw=0.7, ls="--")
    ax.set_xlabel("x/c")
    ax.set_ylabel("y/c")
    ax.set_title("NACA %s section (closed trailing edge)" % CODE, pad=10)
    ax.set_aspect("equal")
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-THICKNESS_RATIO, THICKNESS_RATIO)
    RECORD.append((path.name, float(plt.rcParams["savefig.dpi"]), pm.smallest_lettering(fig)))
    fig.savefig(path)
    plt.close(fig)


def thickness_figure(x: np.ndarray, yt: np.ndarray, summary: pd.DataFrame, path: Path) -> None:
    """Thickness against chordwise position, with the maximum marked."""
    val = summary.set_index("property")["value"]
    x_tmax, t_max = float(val["x_at_max_thickness_x_c"]), float(val["max_thickness_t_c"])
    fig, ax = plt.subplots(figsize=(8, 4.2))
    ax.plot(x, 2.0*yt, color=PALETTE[0], lw=2.0, label="thickness  2 $y_t$/c")
    ax.plot([x_tmax], [t_max], "o", color=PALETTE[1], ms=7,
            label="maximum %.2f %% of chord at x/c = %.3f" % (100.0*t_max, x_tmax))
    ax.set_xlabel("x/c")
    ax.set_ylabel("thickness / chord")
    ax.set_title("NACA %s thickness distribution" % CODE, pad=10)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2, frameon=True)
    RECORD.append((path.name, float(plt.rcParams["savefig.dpi"]), pm.smallest_lettering(fig)))
    fig.savefig(path)
    plt.close(fig)


def main() -> None:
    """Write the coordinates, the summary and the two figures."""
    apply_style()
    X, Y, x, yt = naca4.outline(CODE, pm.N_SURFACE_POINTS)
    pd.DataFrame({"node_id": np.arange(1, len(X) + 1), "x_over_c": np.round(X, COORDINATE_DECIMALS),
                  "y_over_c": np.round(Y, COORDINATE_DECIMALS) + 0.0}).to_csv(HERE/"naca0012_coordinates.csv",
                                                                             index=False)
    summary = summary_table(X, Y, yt)
    summary.to_csv(HERE/"section_geometry_summary.csv", index=False)
    profile_figure(X, Y, HERE/"fig_geometry_profile.png")
    thickness_figure(x, yt, summary, HERE/"fig_geometry_thickness.png")
    pm.record_figures(HERE/"figure_record.csv", RECORD)
    val = summary.set_index("property")["value"]
    print("[geometry] %d surface points; t/c = %s at x/c = %s; area/c^2 = %s"
          % (len(X), val["max_thickness_t_c"], val["x_at_max_thickness_x_c"], val["enclosed_area_A_c2"]))


if __name__ == "__main__":
    main()
