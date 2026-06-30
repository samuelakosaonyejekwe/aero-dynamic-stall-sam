"""
01_geometry / generate_geometry.py
-----------------------------------
Generates the airfoil section geometry for the UNISTALL(TM) dynamic-stall
case study (helicopter main-rotor retreating blade).

Primary section : NACA 0012  (baseline rotor airfoil, the canonical dynamic-stall
                  validation geometry of McAlister/Carr/McCroskey).
Reference section: a cambered, drooped-LE SC1095-class section (for context).

Outputs
  naca0012_coordinates.csv      surface coordinates (x/c, y/c) + node id
  section_geometry_summary.csv  scalar geometric properties
  fig_geometry_profile.png      clean profile plot
  fig_geometry_thickness.png    thickness / camber distribution
"""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from aero_style import apply_style, PALETTE, INK, INK_SOFT
apply_style()

CHORD = 0.30          # m, model-scale chord (typical oscillating-airfoil rig)
N_PTS = 160           # surface points per side (cosine-clustered)

def naca_4digit(code="0012", n=N_PTS):
    """Return closed surface coords (TE->upper->LE->lower->TE), cosine clustered."""
    m = int(code[0]) / 100.0
    p = int(code[1]) / 10.0
    t = int(code[2:]) / 100.0
    beta = np.linspace(0.0, np.pi, n)
    x = (1 - np.cos(beta)) / 2.0           # cosine clustering -> dense at LE/TE
    yt = 5*t*(0.2969*np.sqrt(x) - 0.1260*x - 0.3516*x**2
              + 0.2843*x**3 - 0.1015*x**4)  # (note: open TE form)
    if m > 0:
        yc = np.where(x < p, m/p**2*(2*p*x - x**2),
                      m/(1-p)**2*((1-2*p) + 2*p*x - x**2))
        dyc = np.where(x < p, 2*m/p**2*(p - x),
                       2*m/(1-p)**2*(p - x))
        th = np.arctan(dyc)
    else:
        yc = np.zeros_like(x); th = np.zeros_like(x)
    xu, yu = x - yt*np.sin(th), yc + yt*np.cos(th)
    xl, yl = x + yt*np.sin(th), yc - yt*np.cos(th)
    # assemble: upper TE->LE then lower LE->TE
    X = np.concatenate([xu[::-1], xl[1:]])
    Y = np.concatenate([yu[::-1], yl[1:]])
    return X, Y, x, yt, yc

X, Y, xc, yt, yc = naca_4digit("0012")

# ---- write coordinate CSV ----
df = pd.DataFrame({"node_id": np.arange(1, len(X)+1),
                   "x_over_c": np.round(X, 6),
                   "y_over_c": np.round(Y, 6),
                   "x_m": np.round(X*CHORD, 6),
                   "y_m": np.round(Y*CHORD, 6)})
df.to_csv(HERE/"naca0012_coordinates.csv", index=False)

# ---- geometric properties ----
tmax = 2*yt.max()
xtmax = xc[np.argmax(yt)]
# enclosed area via shoelace
area = 0.5*np.abs(np.dot(X, np.roll(Y, -1)) - np.dot(Y, np.roll(X, -1)))
# LE radius for 4-digit: r/c = 1.1019 t^2
le_radius = 1.1019*(0.12**2)
props = pd.DataFrame({
    "property": ["section", "chord_m", "max_thickness_t_c", "x_at_max_thickness_x_c",
                 "max_camber_y_c", "LE_radius_r_c", "TE_type", "enclosed_area_A_c2",
                 "n_surface_points"],
    "value": ["NACA 0012", CHORD, round(tmax,5), round(xtmax,4),
              round(yc.max(),5), round(le_radius,5), "finite (0.252% open)",
              round(area,5), len(X)],
    "units": ["-", "m", "-", "-", "-", "-", "-", "-", "-"]})
props.to_csv(HERE/"section_geometry_summary.csv", index=False)

# ---- figure: profile ----
fig, ax = plt.subplots(figsize=(9, 3.9))
ax.fill(X, Y, color="#3b6ea5", alpha=0.25)
ax.plot(X, Y, color="#1f3b5c", lw=1.8)
ax.plot(xc, yc, "--", color="#b03a2e", lw=1.2, label="mean camber line")
ax.axhline(0, color="grey", lw=0.6)
ax.set_xlabel("x/c"); ax.set_ylabel("y/c")
ax.set_title("NACA 0012 section — UNISTALL case geometry (chord = %.2f m)" % CHORD)
ax.set_aspect("equal"); ax.set_xlim(-0.05, 1.05)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.42), ncol=2, frameon=True)
fig.savefig(HERE/"fig_geometry_profile.png"); plt.close(fig)

# ---- figure: thickness/camber distribution ----
fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(xc, 2*yt, color="#1f3b5c", lw=1.8, label="thickness 2·y_t/c")
ax.plot(xc, yc, color="#b03a2e", lw=1.8, label="camber y_c/c")
ax.axvline(xtmax, color="grey", ls=":", lw=1)
ax.set_ylim(-0.01, 0.145)
ax.annotate("t_max=%.1f%% @ x/c=%.2f" % (tmax*100, xtmax),
            (xtmax, 2*yt.max()), xytext=(xtmax+0.08, 0.134),
            arrowprops=dict(arrowstyle="->", color="grey"), fontsize=9)
ax.set_xlabel("x/c"); ax.set_ylabel("distribution")
ax.set_title("Thickness & camber distribution", pad=10)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2, frameon=True)
fig.savefig(HERE/"fig_geometry_thickness.png"); plt.close(fig)

print("[geometry] wrote %d surface points; t/c=%.4f, A/c^2=%.5f" % (len(X), tmax, area))
