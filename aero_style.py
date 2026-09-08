"""
aero_style.py  — shared plotting style for the whole UNISTALL case study.
HARD RULE: never use black. All text/axes/ticks/lines use deep navy (#1f3350)
or coloured tones; all colormaps are black-free.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from cycler import cycler

INK      = "#1f3350"   # deep navy — replaces black for all text/axes
INK_SOFT = "#4a5d78"   # muted navy-grey for secondary elements
GRID     = "#9fb3c8"   # light blue-grey grid

# qualitative palette (no black)
PALETTE = ["#2e6fb7", "#c0392b", "#27895f", "#e08a1e",
           "#7b4fb5", "#d14f8c", "#1aa3a3", "#8a6d3b"]

# black-free colormaps for fields
CMAP_PRESSURE = "turbo"      # velocity / pressure magnitude
CMAP_CP       = "Spectral_r" # surface pressure coefficient
CMAP_TEMP     = "plasma"     # temperature (low end deep blue, NOT black)
CMAP_VORT     = "coolwarm"   # vorticity / signed fields

def apply_style():
    plt.rcParams.update({
        "font.size": 11,
        "figure.dpi": 150,
        "savefig.dpi": 150,
        "savefig.bbox": "tight",
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "text.color":       INK,
        "axes.labelcolor":  INK,
        "axes.titlecolor":  INK,
        "axes.edgecolor":   INK_SOFT,
        "axes.titleweight": "bold",
        "xtick.color":      INK,
        "ytick.color":      INK,
        "axes.grid": True,
        "grid.color": GRID,
        "grid.alpha": 0.45,
        "grid.linewidth": 0.7,
        "axes.prop_cycle": cycler(color=PALETTE),
        "legend.framealpha": 0.9,
        "legend.edgecolor": INK_SOFT,
    })
    return PALETTE
