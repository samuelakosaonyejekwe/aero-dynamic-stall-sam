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
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import unistall_solver as us
apply_style()

SOL = ROOT/"05_solution"; OUT = HERE/"plots"; OUT.mkdir(exist_ok=True)
SETUP = ROOT/"03_model_setup"; GEO = ROOT/"01_geometry"/"naca0012_coordinates.csv"
AF = pd.read_csv(GEO)

def tidy3d(ax):
    for pane in (ax.xaxis, ax.yaxis, ax.zaxis):
        pane.pane.set_edgecolor(INK_SOFT); pane.pane.set_alpha(0.04)
    ax.grid(True)
    ax.xaxis.labelpad = 12; ax.yaxis.labelpad = 12; ax.zaxis.labelpad = 8
    ax.tick_params(pad=3)

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
pd.DataFrame({"alpha_mean_deg": M.ravel().round(3), "reduced_freq_k": K.ravel().round(4),
              "CL_max": CLmax.ravel().round(4), "CM_min": CMmin.ravel().round(4)}
             ).to_csv(SOL/"response_surface.csv", index=False)

fig = plt.figure(figsize=(8.5, 6.2)); ax = fig.add_subplot(111, projection="3d")
surf = ax.plot_surface(M, K, CLmax, cmap=CMAP_PRESSURE, edgecolor=INK_SOFT,
                       linewidth=0.3, antialiased=True, alpha=0.95)
ax.set_xlabel("mean α [deg]"); ax.set_ylabel("reduced freq k")
ax.set_zlabel("dynamic $C_{L,max}$")
ax.set_title("Dynamic-stall lift response surface  $C_{L,max}(α_{mean}, k)$", pad=18)
cb = fig.colorbar(surf, ax=ax, pad=0.10, shrink=0.6); cb.set_label("$C_{L,max}$")
ax.view_init(elev=24, azim=-60); tidy3d(ax)
fig.savefig(OUT/"fig3d_response_surface.png", bbox_inches="tight", pad_inches=0.35); plt.close(fig)

# ====================================================== 2. Cp(x/c, phase) SURFACE
nph = 16
o = us.solve_dynamic_stall(10.0, 10.0, 0.10, 0.30, 0.30, 102.0, f_static,
                           CNalpha=CNALPHA, consts=consts, n_per_cycle=720, n_cycles=6)
idxs = np.linspace(0, len(o["alpha_deg"])-1, nph).astype(int)
xcp = np.linspace(0.02, 0.98, 60)      # common upper-surface x/c grid
Zsurf = []; phases = []
for ii in idxs:
    xoc, cp, upper = us.surface_cp(GEO, 0.30, 102.0, 0.30, o["alpha_deg"][ii],
                                   o["CL"][ii], o["CNv"][ii], o["tau_v"][ii]/consts["Tvl"],
                                   nx_grid=180, ny_grid=140)
    # upper-surface Cp resampled on the common x/c grid
    cu = np.interp(xcp, np.sort(xoc[upper]), cp[upper][np.argsort(xoc[upper])])
    Zsurf.append(cu); phases.append(o["phase_deg"][ii])
Zsurf = np.array(Zsurf)
XC, PH = np.meshgrid(xcp, np.array(phases), indexing="xy")
fig = plt.figure(figsize=(8.5, 6.2)); ax = fig.add_subplot(111, projection="3d")
s = ax.plot_surface(XC, PH, Zsurf, cmap=CMAP_CP, edgecolor=INK_SOFT, linewidth=0.2, alpha=0.95)
ax.set_xlabel("x/c"); ax.set_ylabel("cycle phase ωt [deg]"); ax.set_zlabel("upper $C_p$")
ax.invert_zaxis()
ax.set_title("Upper-surface $C_p$ evolution through the cycle (Case A)", pad=18)
cb = fig.colorbar(s, ax=ax, pad=0.10, shrink=0.6); cb.set_label("$C_p$")
ax.view_init(elev=26, azim=-52); tidy3d(ax)
fig.savefig(OUT/"fig3d_cp_phase_surface.png", bbox_inches="tight", pad_inches=0.35); plt.close(fig)

# ====================================================== 3. FIELD AS 3D SURFACE
def load_field(path):
    df = pd.read_csv(path)
    xu = np.unique(df["x_m"].values); yu = np.unique(df["y_m"].values)
    nx, ny = len(xu), len(yu)
    flds = {c: df[c].values.reshape(ny, nx) for c in df.columns if c not in ("x_m", "y_m")}
    return xu, yu, flds

_pk = sorted(glob.glob(str(SOL/"field_A_validation_peak*.csv")))
if not _pk:
    raise SystemExit("[3d] no field_A_validation_peak*.csv — run 04_solver/run_case.py first")
peakA = _pk[0]
xu, yu, F = load_field(peakA)
X, Y = np.meshgrid(xu, yu)
for key, cmap, lab, fn in [("Cp", CMAP_CP, "$C_p$", "Cp"),
                           ("speed_ms", CMAP_PRESSURE, "|V| [m/s]", "speed")]:
    Z = np.array(F[key], dtype=float)
    # Leave the masked airfoil interior as NaN so it renders as a hole. Filling
    # it with the global minimum (the Cp clip floor) punched a canyon several
    # units deep that set the z-scale and flattened the actual field.
    lo, hi = np.nanpercentile(Z, [0.5, 99.5])
    Zc = np.clip(Z, lo, hi)
    fig = plt.figure(figsize=(8.5, 6.0)); ax = fig.add_subplot(111, projection="3d")
    s = ax.plot_surface(X, Y, Zc, cmap=cmap, vmin=lo, vmax=hi,
                        linewidth=0, antialiased=True, alpha=0.96, rstride=2, cstride=2)
    ax.set_zlim(lo, hi)
    ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]"); ax.set_zlabel(lab)
    ax.set_title(f"3D field surface — {lab}  (Case A, peak incidence)", pad=18)
    cb = fig.colorbar(s, ax=ax, pad=0.10, shrink=0.6); cb.set_label(lab)
    ax.view_init(elev=40, azim=-58); tidy3d(ax)
    fig.savefig(OUT/f"fig3d_field_surface_{fn}.png", bbox_inches="tight", pad_inches=0.35); plt.close(fig)

# ====================================================== 4. PICTORIAL SECTION + VECTORS
xu, yu, F = load_field(peakA)
X, Y = np.meshgrid(xu, yu)
span = 0.45   # m, pseudo-span for the pictorial extrusion
fig = plt.figure(figsize=(9, 6.2)); ax = fig.add_subplot(111, projection="3d")
# extrude airfoil section at two span stations, filled so it reads as a solid
afx = AF["x_over_c"].values*0.30; afy = AF["y_over_c"].values*0.30
# fill BOTH end sections so the extrusion reads as a solid blade rather than
# two thin sticks lost among the arrows
for z in (0.0, span):
    ax.add_collection3d(Poly3DCollection(
        [list(zip(afx, np.full_like(afx, z), afy))], facecolor="#c6d7ea",
        edgecolor=INK, linewidths=1.6, alpha=1.0, zorder=10))
for frac in np.linspace(0, 1, 12):
    k = int(frac*(len(afx)-1))
    ax.plot([afx[k], afx[k]], [0, span], [afy[k], afy[k]], color=INK_SOFT, lw=0.5, alpha=0.6)
# 3-D velocity vectors on a mid-span plane. Sampling the whole domain every
# 10th node produced a solid mat of arrows across the entire box that buried the
# section; restrict them to a window around the aerofoil and thin them out.
sk = (slice(None, None, 12), slice(None, None, 12))
xs = X[sk]; ys = Y[sk]; us_ = F["u_ms"][sk]; vs = F["v_ms"][sk]
win = (~np.isnan(us_)) & (xs > -0.18) & (xs < 0.52) & (np.abs(ys) < 0.20)
zc = span/2
# Cap the plotted magnitude. Arrow length is proportional to |V|, and a handful
# of near-surface cells carry several times the freestream speed, so unclipped
# they drew metre-long arrows off the top of the axes.
uu, vv = us_[win].copy(), vs[win].copy()
mag = np.hypot(uu, vv)
cap = np.nanpercentile(mag, 96)
scale = np.where(mag > cap, cap/np.maximum(mag, 1e-9), 1.0)
uu *= scale; vv *= scale
ax.quiver(xs[win], np.full(win.sum(), zc), ys[win],
          uu, np.zeros(win.sum()), vv,
          length=0.0026, normalize=False, color=PALETTE[0], linewidth=1.1,
          arrow_length_ratio=0.30, zorder=2)
ax.set_xlim(X.min(), X.max()); ax.set_zlim(Y.min(), Y.max()); ax.set_ylim(0, span)
ax.set_xlabel("x [m]"); ax.set_ylabel("span z [m]"); ax.set_zlabel("y [m]")
ax.set_title("Pictorial blade section with reconstructed velocity field (peak incidence)",
             pad=20)
ax.set_yticks(np.linspace(0, span, 4))
ax.view_init(elev=18, azim=-68); tidy3d(ax)
ax.set_box_aspect((1.4, 0.8, 0.7))
fig.savefig(OUT/"fig3d_section_vectors.png", bbox_inches="tight", pad_inches=0.35); plt.close(fig)

print("[3d] response surface + Cp-phase surface + field surfaces + pictorial section done")
