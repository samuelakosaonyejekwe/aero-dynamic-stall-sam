# -*- coding: utf-8 -*-
"""
02_mesh / generate_mesh.py
--------------------------
A structured O-grid around the NACA 0012 section, in chord units.

Standing: an illustration. No equation of the flow is solved on this grid and
no load depends on it.

Purpose. The load model (package unistall) has no spatial mesh and does not read
this grid. The grid exists so that the reconstructed flow field can be
evaluated at its nodes, and for illustration. It is not a CFD mesh that has
been run: no flow solution has been computed on it, and its wall spacing is
not sized for a boundary layer.

Construction. The wall line runs trailing edge -> upper surface -> leading
edge -> lower surface -> trailing edge, at x = (1 + cos(phi))/2 with
phi = theta + EDGE_BIAS sin(theta) and equal steps of theta: nodes cluster at
both edges, and the bias moves nodes from the rear half of the section to the
front half. The section is the closed-edge NACA 0012 of
unistall.naca4.closed_thickness, the outline of 01_geometry, of the panels of
the flow-field reconstruction and of the drawings. Each wall node is
marched outward along a direction that turns from the smoothed wall normal to
the direction of its far-field node, with geometric growth of the step. The
outer boundary is a circle about mid-chord with equally spaced nodes. There is
no wake cut. Index i = 0 and i = N_WRAP - 1 are the same points (the seam).

Quality. The trailing edge is a sharp wedge, and an O-grid has to turn round
it at one node: the two cells that touch that node have a corner of about
172 degrees, which sets the maximum skewness and the minimum orthogonality
angle of the grid. The quality table gives the measures both for all cells
and with those two cells left out. The script stops without writing anything
if any cell has a non-positive signed area.

Outputs (all in this folder)
  mesh_nodes.csv             i, j, x_over_c, y_over_c, wall_distance_over_c
  mesh_quality_metrics.csv   metric, value, unit
  mesh_radial_spacing.csv    the wall-normal spacing law
  fig_mesh_full.png          whole domain and near field
  fig_mesh_le_zoom.png       leading-edge region
  fig_mesh_te_zoom.png       trailing-edge region
  fig_mesh_wall_spacing.png  layer height and distance from the wall

Run from the repository root:  PYTHONPATH=. python3 02_mesh/generate_mesh.py

Author: Akosa Samuel Onyejekwe (independent)
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
RECORD = []                 # (file, dpi, smallest lettering) of each figure written

THICKNESS_RATIO = naca4.thickness_ratio(pm.SECTION)
N_WRAP = 193                # nodes round the section, seam node counted twice
N_NORMAL = 81               # nodes from the wall to the far field
FARFIELD_RADIUS = 20.0      # chords, about mid-chord
FIRST_HEIGHT = 1.0e-3       # chords, height of the first layer
NORMAL_SMOOTHING_PASSES = 5
BLEND_EXPONENT = 0.75       # direction weight = (distance / radius) ** exponent
EDGE_BIAS = 0.3             # 0: cosine spacing in x; > 0: more nodes on the front half
CENTRE = (0.5, 0.0)


def wall_line():
    """Wall nodes, counter-clockwise from the trailing edge at (1, 0)."""
    theta = np.linspace(0.0, 2.0*np.pi, N_WRAP)
    x = 0.5*(1.0 + np.cos(theta + EDGE_BIAS*np.sin(theta)))
    y = np.sign(np.sin(theta))*naca4.closed_thickness(x, THICKNESS_RATIO)
    x[0] = x[-1] = 1.0
    y[0] = y[-1] = 0.0
    y[(N_WRAP - 1)//2] = 0.0
    return theta, x, y


def wall_normals(x, y):
    """Outward unit normals of the closed wall loop, smoothed along the wall."""
    px, py = x[:-1], y[:-1]                       # the unique nodes of the loop
    tx = np.roll(px, -1) - np.roll(px, 1)
    ty = np.roll(py, -1) - np.roll(py, 1)
    nx, ny = ty, -tx                              # outward for a counter-clockwise loop
    norm = np.hypot(nx, ny)
    nx, ny = nx/norm, ny/norm
    for _ in range(NORMAL_SMOOTHING_PASSES):
        nx = 0.25*np.roll(nx, 1) + 0.5*nx + 0.25*np.roll(nx, -1)
        ny = 0.25*np.roll(ny, 1) + 0.5*ny + 0.25*np.roll(ny, -1)
        norm = np.hypot(nx, ny)
        nx, ny = nx/norm, ny/norm
    return np.append(nx, nx[0]), np.append(ny, ny[0])


def growth_ratio(n_layers, first, total):
    """Ratio r of a geometric series of n_layers steps, first step given, summing to total."""
    return brentq(lambda r: first*(r**n_layers - 1.0)/(r - 1.0) - total, 1.0 + 1e-6, 2.0)


def normal_coordinate():
    """Distance from the wall of each layer, and the growth ratio."""
    ratio = growth_ratio(N_NORMAL - 1, FIRST_HEIGHT, FARFIELD_RADIUS)
    steps = FIRST_HEIGHT*ratio**np.arange(N_NORMAL - 1)
    return np.concatenate([[0.0], np.cumsum(steps)]), ratio


def build_grid():
    """Node coordinates X[i, j], Y[i, j], the wall line, the layer distances and the growth ratio."""
    theta, xw, yw = wall_line()
    nx, ny = wall_normals(xw, yw)
    dist, ratio = normal_coordinate()
    xf = CENTRE[0] + FARFIELD_RADIUS*np.cos(theta)
    yf = CENTRE[1] + FARFIELD_RADIUS*np.sin(theta)
    reach = np.hypot(xf - xw, yf - yw)            # wall node to its far-field node
    fx, fy = (xf - xw)/reach, (yf - yw)/reach
    frac = dist/dist[-1]
    X = np.zeros((N_WRAP, N_NORMAL))
    Y = np.zeros((N_WRAP, N_NORMAL))
    for j in range(N_NORMAL):
        w = frac[j]**BLEND_EXPONENT
        dx, dy = (1.0 - w)*nx + w*fx, (1.0 - w)*ny + w*fy
        norm = np.hypot(dx, dy)
        length = dist[j] + (reach - dist[-1])*frac[j]**2   # ends exactly on the circle
        X[:, j] = xw + length*dx/norm
        Y[:, j] = yw + length*dy/norm
    return dict(X=X, Y=Y, xw=xw, yw=yw, dist=dist, ratio=ratio)


def cell_measures(X, Y):
    """Signed area, aspect ratio and corner angles (degrees) of every quadrilateral cell."""
    px = np.stack([X[:-1, :-1], X[:-1, 1:], X[1:, 1:], X[1:, :-1]])   # counter-clockwise corners
    py = np.stack([Y[:-1, :-1], Y[:-1, 1:], Y[1:, 1:], Y[1:, :-1]])
    area = 0.5*np.sum(px*np.roll(py, -1, axis=0) - py*np.roll(px, -1, axis=0), axis=0)
    ex, ey = np.roll(px, -1, axis=0) - px, np.roll(py, -1, axis=0) - py
    edge = np.hypot(ex, ey)
    ux, uy = ex/edge, ey/edge
    cosine = -(ux*np.roll(ux, 1, axis=0) + uy*np.roll(uy, 1, axis=0))
    angle = np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))
    return area, edge.max(axis=0)/edge.min(axis=0), angle


def quality_table(grid):
    """The quality measures of the grid as a table of metric, value, unit."""
    X, Y, xw, yw, dist = grid["X"], grid["Y"], grid["xw"], grid["yw"], grid["dist"]
    area, aspect, angle = cell_measures(X, Y)
    deviation = np.abs(angle - 90.0).max(axis=0)
    skew = deviation/90.0
    wall_step = np.hypot(np.diff(xw), np.diff(yw))
    outer = np.hypot(X[:, -1] - CENTRE[0], Y[:, -1] - CENTRE[1])
    worst = np.unravel_index(np.argmax(skew), skew.shape)
    rest = np.ones(skew.shape, bool)
    rest[[0, -1], 0] = False                      # the two cells at the trailing-edge node
    rows = [
        ("topology", "O-grid, no wake cut, seam at the trailing edge", "-"),
        ("status", "not a CFD mesh that has been run; not read by the load model", "-"),
        ("nodes_wrap_i", N_WRAP, "-"),
        ("nodes_normal_j", N_NORMAL, "-"),
        ("total_nodes", N_WRAP*N_NORMAL, "-"),
        ("unique_nodes", (N_WRAP - 1)*N_NORMAL, "-"),
        ("total_cells", (N_WRAP - 1)*(N_NORMAL - 1), "-"),
        ("farfield_radius", round(float(outer.mean()), 4), "chord"),
        ("farfield_radius_deviation_max", float("%.1e" % np.abs(outer - FARFIELD_RADIUS).max()), "chord"),
        ("first_layer_height", float("%.4g" % dist[1]), "chord"),
        ("wall_normal_growth_ratio", round(float(grid["ratio"]), 4), "-"),
        ("wall_spacing_min", float("%.3g" % wall_step.min()), "chord"),
        ("wall_spacing_max", float("%.3g" % wall_step.max()), "chord"),
        ("wall_spacing_at_leading_edge", float("%.3g" % wall_step[(N_WRAP - 1)//2]), "chord"),
        ("wall_spacing_at_trailing_edge", float("%.3g" % wall_step[0]), "chord"),
        ("max_skewness_equiangle", round(float(skew.max()), 3), "-"),
        ("max_skewness_cell_i_j", "%d %d" % worst, "-"),
        ("max_skewness_without_trailing_edge_cells", round(float(skew[rest].max()), 3), "-"),
        ("mean_skewness_equiangle", round(float(skew.mean()), 3), "-"),
        ("cells_with_skewness_above_0.5", round(100.0*float(np.mean(skew > 0.5)), 2), "percent"),
        ("min_orthogonality_angle", round(float(90.0 - deviation.max()), 1), "deg"),
        ("min_orthogonality_angle_without_trailing_edge_cells",
         round(float(90.0 - deviation[rest].max()), 1), "deg"),
        ("max_aspect_ratio", round(float(aspect.max()), 1), "-"),
        ("mean_aspect_ratio", round(float(aspect.mean()), 1), "-"),
        ("min_cell_area_signed", "%.3e" % area.min(), "chord^2"),
        ("inverted_cells", int(np.sum(area <= 0.0)), "-"),
    ]
    return pd.DataFrame(rows, columns=["metric", "value", "unit"])


def write_tables(grid, metrics):
    """Write the node list, the quality table and the wall-normal spacing law."""
    X, Y, dist = grid["X"], grid["Y"], grid["dist"]
    ii, jj = np.meshgrid(np.arange(N_WRAP), np.arange(N_NORMAL), indexing="ij")
    wall = np.hypot(X - X[:, :1], Y - Y[:, :1])
    pd.DataFrame({"i": ii.ravel(), "j": jj.ravel(), "x_over_c": X.ravel(), "y_over_c": Y.ravel(),
                  "wall_distance_over_c": wall.ravel()}
                 ).to_csv(HERE/"mesh_nodes.csv", index=False, float_format="%.9g")
    metrics.to_csv(HERE/"mesh_quality_metrics.csv", index=False)
    pd.DataFrame({"layer_j": np.arange(N_NORMAL), "normal_coordinate_over_c": dist,
                  "layer_height_over_c": np.concatenate([[np.nan], np.diff(dist)])}
                 ).to_csv(HERE/"mesh_radial_spacing.csv", index=False, float_format="%.9g")


def draw_grid(ax, grid, every_i, every_j, lw):
    """Draw every every_i-th normal line and every every_j-th wrap line, and say which."""
    X, Y = grid["X"], grid["Y"]
    for j in range(0, N_NORMAL, every_j):
        ax.plot(X[:, j], Y[:, j], color=INK_SOFT, lw=lw)
    if (N_NORMAL - 1) % every_j:
        ax.plot(X[:, -1], Y[:, -1], color=INK_SOFT, lw=lw)
    for i in range(0, N_WRAP - 1, every_i):
        ax.plot(X[i, :], Y[i, :], color=INK_SOFT, lw=lw)
    ax.fill(grid["xw"], grid["yw"], color=PALETTE[0], zorder=3)
    ax.set_aspect("equal")
    ax.grid(False)
    ax.set_xlabel("x / c")
    ax.set_ylabel("y / c")
    shown = "all grid lines drawn" if every_i == every_j == 1 else (
        "drawn: one line in %d round the section, one in %d outward" % (every_i, every_j))
    ax.text(0.5, -0.16, shown, transform=ax.transAxes, ha="center", va="top", fontsize=9.5, color=INK_SOFT)


def write_figure(fig, name):
    """Write a figure and enter it in the figure record."""
    RECORD.append((name, float(plt.rcParams["savefig.dpi"]), pm.smallest_lettering(fig)))
    fig.savefig(HERE/name)
    plt.close(fig)


def figure_full(grid):
    """The whole domain beside the near field."""
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 5.2))
    draw_grid(axes[0], grid, 4, 4, 0.35)
    r = 1.05*FARFIELD_RADIUS
    axes[0].set_xlim(CENTRE[0] - r, CENTRE[0] + r)
    axes[0].set_ylim(-r, r)
    axes[0].set_title("Whole domain, radius %.0f chords" % FARFIELD_RADIUS)
    draw_grid(axes[1], grid, 2, 2, 0.4)
    axes[1].set_xlim(-1.0, 2.0)
    axes[1].set_ylim(-1.5, 1.5)
    axes[1].set_title("Near field")
    fig.suptitle("O-grid round the NACA %s section: %d x %d nodes (illustration; no flow solved on it)"
                 % (pm.SECTION, N_WRAP, N_NORMAL), color=INK, fontweight="bold", fontsize=11.5)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    write_figure(fig, "fig_mesh_full.png")


def figure_zoom(grid, name, title, xlim, ylim):
    """A close view with every grid line drawn."""
    fig, ax = plt.subplots(figsize=(4.4, 4.5))
    draw_grid(ax, grid, 1, 1, 0.5)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_title(title, pad=12, fontsize=11)
    write_figure(fig, name)


def figure_spacing(grid):
    """Layer height and distance from the wall against the layer index."""
    dist = grid["dist"]
    layer = np.arange(1, N_NORMAL)
    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    h1, = ax.semilogy(layer, np.diff(dist), color=PALETTE[1], lw=2, marker="o", ms=3, label="layer height")
    h2, = ax.semilogy(layer, dist[1:], color=PALETTE[2], lw=1.8, ls="--", label="distance from the wall")
    ax.set_xlabel("layer index j")
    ax.set_ylabel("length / c")
    ax.set_title("Wall-normal spacing: first layer %.1e c, growth ratio %.4f" % (dist[1], grid["ratio"]))
    ax.legend(handles=[h1, h2], loc="upper left")
    write_figure(fig, "fig_mesh_wall_spacing.png")


def main():
    apply_style()
    grid = build_grid()
    metrics = quality_table(grid)
    inverted = int(metrics.set_index("metric").loc["inverted_cells", "value"])
    if inverted:
        raise SystemExit("[mesh] %d inverted cells: grid not valid, nothing written" % inverted)
    write_tables(grid, metrics)
    figure_full(grid)
    figure_zoom(grid, "fig_mesh_le_zoom.png", "Leading-edge region", (-0.06, 0.14), (-0.10, 0.10))
    figure_zoom(grid, "fig_mesh_te_zoom.png", "Trailing-edge region (no wake cut)", (0.90, 1.10), (-0.10, 0.10))
    figure_spacing(grid)
    pm.record_figures(HERE/"figure_record.csv", RECORD)
    print(metrics.to_string(index=False))


if __name__ == "__main__":
    main()
