"""
06_postprocessing / make_3d_plots.py
------------------------------------
3D engineering visualisations:
  * fig3d_response_surface.png   CL_max response surface over (mean-α, k)  [+CSV]
  * fig3d_cp_phase_surface.png   Cp(x/c, cycle-phase) carpet surface
  * fig3d_field_surface_*.png    pressure / speed field as a 3D surface
  * fig3d_section_vectors_*.png   pictorial blade section, pressure-coloured,
                                  with 3D velocity vectors
No black; clean layouts; labels kept clear of the surfaces.
"""
import sys, glob
import numpy as np
import pandas as pd
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT/"04_solver"))
from aero_style import (apply_style, PALETTE, INK, INK_SOFT,
                        CMAP_PRESSURE, CMAP_CP, CMAP_TEMP)
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa
import unistall_solver as us
apply_style()

SOL = ROOT/"05_solution"; OUT = HERE/"plots"; OUT.mkdir(exist_ok=True)
SETUP = ROOT/"03_model_setup"; GEO = ROOT/"01_geometry"/"naca0012_coordinates.csv"
AF = pd.read_csv(GEO)

def tidy3d(ax):
    for pane in (ax.xaxis, ax.yaxis, ax.zaxis):
        pane.pane.set_edgecolor(INK_SOFT); pane.pane.set_alpha(0.04)
    ax.grid(True)
    ax.xaxis.labelpad = 14; ax.yaxis.labelpad = 14; ax.zaxis.labelpad = 10
    ax.tick_params(pad=4)

# calibrate solver (same as run_case)
import json
cfg = json.load(open(SETUP/"solver_config.json"))
CNALPHA = cfg["lift_curve_slope_CNalpha_per_rad"]
consts = dict(**cfg["indicial_circulatory"], **cfg["time_constants_semichords"])
consts.update({k: v for k, v in cfg["calibrated_constants"].items() if k != "comment"})
stat = pd.read_csv(SETUP/"static_polar_reference.csv")
f_static = us.calibrate_separation(stat["alpha_deg"], stat["Cl"], stat["Cd"], CNALPHA)

# ====================================================== 1. RESPONSE SURFACE
means = np.linspace(6, 16, 6)
ks    = np.linspace(0.04, 0.16, 6)
CLmax = np.zeros((len(means), len(ks)))
CMmin = np.zeros_like(CLmax)
for i, am in enumerate(means):
    for j, kk in enumerate(ks):
        o = us.solve_dynamic_stall(am, 8.0, kk, 0.30, 0.30, 102.0, f_static,
                                   CNalpha=CNALPHA, consts=consts,
                                   n_per_cycle=360, n_cycles=3)
        CLmax[i, j] = o["CL"].max(); CMmin[i, j] = o["CM"].min()
M, K = np.meshgrid(means, ks, indexing="ij")
pd.DataFrame({"alpha_mean_deg": M.ravel(), "reduced_freq_k": K.ravel(),
              "CL_max": CLmax.ravel().round(4), "CM_min": CMmin.ravel().round(4)}
             ).to_csv(SOL/"response_surface.csv", index=False)

fig = plt.figure(figsize=(8.5, 6.2)); ax = fig.add_subplot(111, projection="3d")
surf = ax.plot_surface(M, K, CLmax, cmap=CMAP_PRESSURE, edgecolor=INK_SOFT,
                       linewidth=0.3, antialiased=True, alpha=0.95)
ax.set_xlabel("\nmean α [deg]"); ax.set_ylabel("\nreduced freq k")
ax.set_zlabel("dynamic $C_{L,max}$")
ax.set_title("Dynamic-stall lift response surface  $C_{L,max}(α_{mean}, k)$", pad=18)
cb = fig.colorbar(surf, ax=ax, pad=0.10, shrink=0.6); cb.set_label("$C_{L,max}$")
ax.view_init(elev=24, azim=-60); tidy3d(ax)
fig.savefig(OUT/"fig3d_response_surface.png", bbox_inches="tight"); plt.close(fig)

# ====================================================== 2. Cp(x/c, phase) SURFACE
nph = 16
o = us.solve_dynamic_stall(10.0, 10.0, 0.10, 0.30, 0.30, 102.0, f_static,
                           CNalpha=CNALPHA, consts=consts, n_per_cycle=720, n_cycles=6)
idxs = np.linspace(0, len(o["alpha_deg"])-1, nph).astype(int)
xc_ref = None; Zsurf = []; phases = []
for ii in idxs:
    xoc, cp, upper = us.surface_cp(GEO, 0.30, 102.0, 0.30, o["alpha_deg"][ii],
                                   o["CL"][ii], o["CNv"][ii], o["tau_v"][ii]/consts["Tvl"],
                                   nx_grid=180, ny_grid=140)
    # upper-surface Cp resampled on a common x/c grid
    xu = np.linspace(0.02, 0.98, 60)
    cu = np.interp(xu, np.sort(xoc[upper]), cp[upper][np.argsort(xoc[upper])])
    Zsurf.append(cu); phases.append(o["phase_deg"][ii])
Zsurf = np.array(Zsurf)
XC, PH = np.meshgrid(xu, np.array(phases), indexing="xy")
fig = plt.figure(figsize=(8.5, 6.2)); ax = fig.add_subplot(111, projection="3d")
s = ax.plot_surface(XC, PH, Zsurf, cmap=CMAP_CP, edgecolor=INK_SOFT, linewidth=0.2, alpha=0.95)
ax.set_xlabel("\nx/c"); ax.set_ylabel("\ncycle phase ωt [deg]"); ax.set_zlabel("upper $C_p$")
ax.invert_zaxis()
ax.set_title("Upper-surface $C_p$ evolution through the cycle (Case A)", pad=18)
cb = fig.colorbar(s, ax=ax, pad=0.10, shrink=0.6); cb.set_label("$C_p$")
ax.view_init(elev=26, azim=-52); tidy3d(ax)
fig.savefig(OUT/"fig3d_cp_phase_surface.png", bbox_inches="tight"); plt.close(fig)

# ====================================================== 3. FIELD AS 3D SURFACE
def load_field(path):
    df = pd.read_csv(path)
    xu = np.unique(df["x_m"].values); yu = np.unique(df["y_m"].values)
    nx, ny = len(xu), len(yu)
    flds = {c: df[c].values.reshape(ny, nx) for c in df.columns if c not in ("x_m", "y_m")}
    return xu, yu, flds

peakA = [f for f in glob.glob(str(SOL/"field_A_validation_peak*.csv"))][0]
xu, yu, F = load_field(peakA)
X, Y = np.meshgrid(xu, yu)
for key, cmap, lab, fn in [("Cp", CMAP_CP, "$C_p$", "Cp"),
                           ("speed_ms", CMAP_PRESSURE, "|V| [m/s]", "speed")]:
    Z = np.array(F[key], dtype=float)
    fig = plt.figure(figsize=(8.5, 6.0)); ax = fig.add_subplot(111, projection="3d")
    s = ax.plot_surface(X, Y, np.nan_to_num(Z, nan=np.nanmin(Z)), cmap=cmap,
                        linewidth=0, antialiased=True, alpha=0.96, rstride=2, cstride=2)
    ax.set_xlabel("\nx [m]"); ax.set_ylabel("\ny [m]"); ax.set_zlabel(lab)
    ax.set_title(f"3D field surface — {lab}  (Case A, peak incidence)", pad=18)
    cb = fig.colorbar(s, ax=ax, pad=0.10, shrink=0.6); cb.set_label(lab)
    ax.view_init(elev=40, azim=-58); tidy3d(ax)
    fig.savefig(OUT/f"fig3d_field_surface_{fn}.png", bbox_inches="tight"); plt.close(fig)

# ====================================================== 4. PICTORIAL SECTION + VECTORS
import matplotlib.cm as cm
from matplotlib.colors import Normalize
xu, yu, F = load_field(peakA)
X, Y = np.meshgrid(xu, yu)
span = 0.45   # m, pseudo-span for the pictorial extrusion
fig = plt.figure(figsize=(9, 6.2)); ax = fig.add_subplot(111, projection="3d")
# extrude airfoil section at two span stations, colour by chordwise Cp proxy
afx = AF["x_over_c"].values*0.30; afy = AF["y_over_c"].values*0.30
for z in (0.0, span):
    ax.plot(afx, np.full_like(afx, z), afy, color=INK, lw=1.6)
# connect LE & TE to suggest the blade surface
for frac in np.linspace(0, 1, 12):
    k = int(frac*(len(afx)-1))
    ax.plot([afx[k], afx[k]], [0, span], [afy[k], afy[k]], color=INK_SOFT, lw=0.5, alpha=0.6)
# 3D velocity vectors on a mid-span plane
sk = (slice(None, None, 10), slice(None, None, 10))
xs = X[sk]; ys = Y[sk]; us_ = F["u_ms"][sk]; vs = F["v_ms"][sk]
mask = ~np.isnan(us_)
zc = span/2
ax.quiver(xs[mask], np.full(mask.sum(), zc), ys[mask],
          us_[mask], np.zeros(mask.sum()), vs[mask],
          length=0.0016, normalize=False, color=PALETTE[0], linewidth=0.9)
ax.set_xlim(X.min(), X.max()); ax.set_zlim(Y.min(), Y.max()); ax.set_ylim(0, span)
ax.set_xlabel("\nx [m]"); ax.set_ylabel("\nspan z [m]"); ax.set_zlabel("y [m]")
ax.set_title("Pictorial blade section with reconstructed velocity field (peak incidence)",
             pad=20)
ax.view_init(elev=18, azim=-68); tidy3d(ax)
ax.set_box_aspect((1.4, 0.8, 0.7))
fig.savefig(OUT/"fig3d_section_vectors.png", bbox_inches="tight"); plt.close(fig)

print("[3d] response surface + Cp-phase surface + field surfaces + pictorial section done")
