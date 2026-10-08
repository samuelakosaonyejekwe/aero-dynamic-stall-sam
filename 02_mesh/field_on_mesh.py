# -*- coding: utf-8 -*-
"""
02_mesh / field_on_mesh.py
--------------------------
Evaluate the reconstructed flow field at the nodes of the O-grid, for both
cases at the instant of peak lift. Run from the repository root:

    PYTHONPATH=. python3 02_mesh/field_on_mesh.py

Author: Akosa Samuel Onyejekwe (independent)

Standing: an illustration. This is the one use the grid has. The load model
has no spatial mesh; the field is the potential-flow reconstruction of
unistall/flowfield.py, drawn round the lift the model predicts. No equation
of the flow is solved on the grid.

The script marches nothing. The instant is the row of greatest lift of
05_solution/time_history_<case>.csv, and the conditions are those of
03_model_setup. Each node takes its value from one of two places:
  surface     the nodes of the wall line (j = 0): the values on the surface
              itself, unistall.flowfield.surface_state, interpolated along
              each surface to the node;
  field       every node further from the wall than
              unistall.flowfield.smoothing_radius: the field evaluated there.
The nodes between, where the field shows the junctions of the panels, carry
no value (source = none).

Outputs: 05_solution/field_<case>_peak_on_mesh.csv.gz (one row per node) and
         02_mesh/fig_mesh_field_<case>.png (pressure coefficient on the grid).
"""
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import AsinhNorm

import project_meta as pm
from unistall import flowfield
from unistall.style import CMAP_CP, INK_SOFT, apply_style

ROOT = Path(__file__).resolve().parents[1]
GEO = ROOT/"01_geometry"/"naca0012_coordinates.csv"
SOL = ROOT/"05_solution"
HERE = ROOT/"02_mesh"
VIEW_CHORDS = (-0.6, 1.8, -1.0, 1.0)           # window of the figure, in chords
RECORD = []                                    # (file, dpi, smallest lettering) of each figure written
GZIP = {"method": "gzip", "compresslevel": 9, "mtime": 0}


def peak_instant(name: str, case: dict, air: flowfield.Air, cfg: dict) -> tuple:
    """(reconstruction, stored row) at the instant of greatest lift of the
    stored time history of a case."""
    th = pd.read_csv(SOL/f"time_history_{name}.csv")
    row = th.loc[th["CL"].idxmax()]
    rec = flowfield.solve(GEO, case["chord"], case["U"], float(row["alpha_deg"]), float(row["CL"]),
                          float(row["CN_vortex"]), float(row["tau_v_semichords"])/case["Tvl"],
                          float(row["f_separation"]),
                          n_panel=cfg["field_reconstruction"]["read_by_the_driver"]["n_panels"], mach=case["M"],
                          air=air)
    return rec, row


def wall_values(rec: flowfield.Reconstruction, case: dict, air: flowfield.Air, x_over_c: np.ndarray,
                y_over_c: np.ndarray) -> tuple:
    """(Cp, speed [m/s], flag of the critical pressure) at wall nodes: the
    surface values interpolated in x/c along the surface each node lies on. A node on the chord line (the
    leading or the trailing edge) takes the mean of the two surfaces."""
    st = flowfield.surface_state(rec, case["M"], air)
    side = {}
    for upper in (True, False):
        pick = st["upper"] == upper
        order = np.argsort(st["x_over_c"][pick])
        xs = st["x_over_c"][pick][order]
        side[upper] = tuple(np.interp(x_over_c, xs, np.asarray(st[q], float)[pick][order])
                            for q in ("Cp", "speed", "beyond_critical"))
    on_upper, on_chord = y_over_c > 0.0, y_over_c == 0.0
    out = []
    for q in (0, 1, 2):
        both = 0.5*(side[True][q] + side[False][q])
        out.append(np.where(on_chord, both, np.where(on_upper, side[True][q], side[False][q])))
    return out[0], out[1], out[2] >= 0.5


def field_on_nodes(name: str, case: dict, air: flowfield.Air, cfg: dict, nodes: pd.DataFrame) -> tuple:
    """The reconstructed field of one case at peak lift, at every mesh node.
    Returns (table with one row per node, incidence at that instant [deg],
    near-wall distance over chord)."""
    rec, row = peak_instant(name, case, air, cfg)
    near = flowfield.smoothing_radius(rec)/case["chord"]
    table = nodes[["i", "j", "x_over_c", "y_over_c", "wall_distance_over_c"]].copy()
    wall = (table["j"] == 0).to_numpy()
    far = (table["wall_distance_over_c"] > near).to_numpy() & ~wall
    Cp, speed, beyond = np.full(len(table), np.nan), np.full(len(table), np.nan), np.zeros(len(table), bool)
    f = flowfield.evaluate(rec, table["x_over_c"].to_numpy()[far]*case["chord"],
                           table["y_over_c"].to_numpy()[far]*case["chord"], case["M"], air, blank_inside=False)
    Cp[far], speed[far], beyond[far] = f["Cp"], f["speed"], f["beyond_critical"]
    Cp[wall], speed[wall], beyond[wall] = wall_values(rec, case, air, table["x_over_c"].to_numpy()[wall],
                                        table["y_over_c"].to_numpy()[wall])
    table["Cp"], table["speed_ms"] = np.round(Cp, 4), np.round(speed, 2)
    table["beyond_critical"] = beyond.astype(int)
    table["source"] = np.where(wall, "surface", np.where(far, "field", "none"))
    return table, float(row["alpha_deg"]), near


def figure(name: str, label: str, table: pd.DataFrame, alpha_deg: float, near: float) -> None:
    """Pressure coefficient on the grid, with the grid lines drawn."""
    ni, nj = int(table.i.max()) + 1, int(table.j.max()) + 1
    X = table.x_over_c.values.reshape(ni, nj)
    Y = table.y_over_c.values.reshape(ni, nj)
    Cp = np.where(table.beyond_critical.values == 1, np.nan, table.Cp.values).reshape(ni, nj)
    lo, hi = float(np.nanmin(Cp)), float(np.nanmax(Cp))
    n_beyond = int(table.beyond_critical.sum())
    fig, ax = plt.subplots(figsize=(7.4, 6.4), constrained_layout=True)
    ax.grid(False)
    pc = ax.pcolormesh(X, Y, np.ma.masked_invalid(Cp), cmap=CMAP_CP, norm=AsinhNorm(1.0, vmin=lo, vmax=hi),
                       shading="nearest", rasterized=True)
    ax.plot(X[:, ::4], Y[:, ::4], color=INK_SOFT, lw=0.25, alpha=0.6)
    ax.plot(X[::4].T, Y[::4].T, color=INK_SOFT, lw=0.25, alpha=0.6)
    bar = fig.colorbar(pc, ax=ax, label="$C_p$, suction upwards (scale spans the data: %.2f to %.2f)" % (lo, hi))
    pm.cp_colour_bar(bar, lo, hi)
    bar.ax.tick_params(labelsize=9)
    ax.set_xlim(VIEW_CHORDS[0], VIEW_CHORDS[1])
    ax.set_ylim(VIEW_CHORDS[2], VIEW_CHORDS[3])
    ax.set_aspect("equal")
    ax.set_xlabel("x / c")
    ax.set_ylabel("y / c")
    ax.set_title(f"{label}: reconstructed field at the grid nodes, peak lift, "
                 rf"$\alpha$ = {alpha_deg:.1f}°", loc="left", fontsize=10.5)
    fig.supxlabel(textwrap.fill(
        "Reconstruction, not a flow solution; no equation of the flow is solved on this grid. Kármán-Tsien "
        f"corrected $C_p$. Wall line: values on the surface itself. Nodes nearer the wall than {near:.4f} c carry "
        f"no value (white), nor do the {n_beyond} nodes where the corrected pressure is below the critical value "
        "(not physical). Every fourth grid line drawn. " + flowfield.AXES_NOTE, 112),
                  fontsize=8.5, color=INK_SOFT)
    RECORD.append((f"fig_mesh_field_{name}.png", 200.0, pm.smallest_lettering(fig)))
    fig.savefig(HERE/f"fig_mesh_field_{name}.png", dpi=200)
    plt.close(fig)


def main() -> None:
    """Write the field at the mesh nodes and its picture for both cases."""
    apply_style()
    nodes = pd.read_csv(HERE/"mesh_nodes.csv")
    cases, air, cfg = pm.read_setup()
    for name, case in cases.items():
        table, alpha, near = field_on_nodes(name, case, air, cfg, nodes)
        table.to_csv(SOL/f"field_{name}_peak_on_mesh.csv.gz", index=False, compression=GZIP)
        figure(name, case["label"], table, alpha, near)
        counts = table["source"].value_counts()
        pm.record_figures(HERE/"figure_record.csv", RECORD)
        print(f"[mesh field] {name}: {len(table)} nodes: {counts.get('surface', 0)} on the wall line from the "
              f"surface values, {counts.get('field', 0)} from the field, {counts.get('none', 0)} nearer the wall "
              f"than {near:.4f} c without a value; alpha {alpha:.2f} deg")


if __name__ == "__main__":
    main()
