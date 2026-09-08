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
                        CMAP_PRESSURE, CMAP_CP)
import matplotlib.pyplot as plt
# NOTE: mpl_toolkits.mplot3d.Axes3D used to be imported here purely for its
# side effect of registering the '3d' projection. That has been built in
# since matplotlib 3.2 and requirements.txt pins >=3.6, so the import was
# dead. Removed and the stage re-run to confirm the 3-D axes still build.
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

# ---- Case A conditions and numerics: READ from 03_model_setup, never restated.
#      This script used to carry its own copies (c=0.30, U=102.0, M=0.30,
#      alpha=10+-10, k=0.10, 360x3 steps) while the README claimed every
#      plotting script reads the setup stage. They are now one source. ----
_flow = pd.read_csv(SETUP/"flow_conditions.csv").set_index("parameter")["case_A_validation"]
_kin  = pd.read_csv(SETUP/"kinematics.csv").set_index("case_id").loc["A_validation_rig"]
C_A = float(_flow["chord_c"]); U_A = float(_flow["freestream_velocity_U"])
M_A = float(_flow["freestream_mach_M"])
AM_A = float(_kin["alpha_mean_deg"]); AA_A = float(_kin["alpha_amp_deg"])
K_A = float(_kin["reduced_freq_k"])
NPC = cfg["numerics"]["steps_per_cycle"]; NCYC = cfg["numerics"]["n_cycles"]
# the sweep runs many points, so it uses a coarser march than the reported cases;
# both numbers are declared here rather than buried as literals
SWEEP_NPC, SWEEP_NCYC = NPC//2, 3

# ====================================================== 1. RESPONSE SURFACE
# sweep grid. Amplitude is held at the case-A value so the surface is a sweep in
# (mean incidence, reduced frequency) about the reported case, not about a third
# unrelated condition; the CSV records the amplitude with the data.
# The upper bound on mean incidence is set by the calibration, not by taste.
# f(alpha) is fitted by inverse Kirchhoff to static_polar_reference.csv, which is
# tabulated to ALPHA_CAL_MAX; past that, calibrate_separation() decays f towards
# full separation on an assumption rather than on data. With the amplitude held
# at the case-A value the peak incidence is mean + amp, so the largest mean that
# keeps every point on calibrated data is ALPHA_CAL_MAX - amp. This grid used to
# run to a mean of 16 deg, i.e. a 26 deg peak, which put half of the 36 published
# points on the extrapolated branch. They are no longer computed: an uncalibrated
# result does not become safe by being labelled.
ALPHA_CAL_MAX = float(stat["alpha_deg"].max())
MEAN_MAX      = ALPHA_CAL_MAX - AA_A
means = np.linspace(6, MEAN_MAX, 6)
ks    = np.linspace(0.04, 0.16, 6)
CLmax = np.zeros((len(means), len(ks)))
CMmin = np.zeros_like(CLmax)
for i, am in enumerate(means):
    for j, kk in enumerate(ks):
        o = us.solve_dynamic_stall(am, AA_A, kk, M_A, C_A, U_A, f_static,
                                   CNalpha=CNALPHA, consts=consts,
                                   n_per_cycle=SWEEP_NPC, n_cycles=SWEEP_NCYC)
        CLmax[i, j] = o["CL"].max(); CMmin[i, j] = o["CM"].min()
M, K = np.meshgrid(means, ks, indexing="ij")
# within_calibration is kept as a published invariant rather than a warning: the
# grid is bounded so every point is True, and the assert below fails the build if
# a future change to the sweep or to the polar ever breaks that.
PEAK = M + AA_A
IN_CAL = PEAK <= ALPHA_CAL_MAX
pd.DataFrame({"alpha_mean_deg": M.ravel().round(3), "reduced_freq_k": K.ravel().round(4),
              "alpha_amp_deg": AA_A, "mach_M": M_A, "chord_m": C_A,
              "steps_per_cycle": SWEEP_NPC, "n_cycles": SWEEP_NCYC,
              "peak_alpha_deg": PEAK.ravel().round(2),
              "within_calibration": IN_CAL.ravel(),
              "CL_max": CLmax.ravel().round(4), "CM_min": CMmin.ravel().round(4)}
             ).to_csv(SOL/"response_surface.csv", index=False)

# reserve a band at the bottom for the calibration note, so it cannot land on
# the x tick labels the way an axes-coordinate annotation did
fig = plt.figure(figsize=(8.5, 6.6))
ax = fig.add_axes([0.0, 0.13, 1.0, 0.82], projection="3d")
surf = ax.plot_surface(M, K, CLmax, cmap=CMAP_PRESSURE, edgecolor=INK_SOFT,
                       linewidth=0.3, antialiased=True, alpha=0.95)
ax.set_xlabel("mean α [deg]"); ax.set_ylabel("reduced freq k")
ax.set_zlabel("dynamic $C_{L,max}$")
ax.set_title("Dynamic-stall lift response surface  $C_{L,max}(α_{mean}, k)$\n"
             "(NACA 0012, M=%.2f, α amplitude %.0f°)" % (M_A, AA_A), pad=18)
cb = fig.colorbar(surf, ax=ax, pad=0.10, shrink=0.6); cb.set_label("$C_{L,max}$")
# state the bound the grid is drawn inside, so the limit travels with the figure
assert IN_CAL.all(), "response surface must not contain extrapolated points"
_note = ("Every point lies on calibrated data: the mean incidence stops at %.0f° so the peak\n"
         "(mean + %.0f° amplitude) never exceeds the %.0f° the static polar is tabulated to.\n"
         "Both reported cases peak at exactly %.0f°."
         % (MEAN_MAX, AA_A, ALPHA_CAL_MAX, ALPHA_CAL_MAX))
fig.text(0.5, 0.015, _note, ha="center", va="bottom", fontsize=8.5, color=INK_SOFT)
ax.view_init(elev=24, azim=-60); tidy3d(ax)
fig.savefig(OUT/"fig3d_response_surface.png", bbox_inches="tight", pad_inches=0.35); plt.close(fig)

# ====================================================== 2. Cp(x/c, phase) SURFACE
nph = 16
o = us.solve_dynamic_stall(AM_A, AA_A, K_A, M_A, C_A, U_A, f_static,
                           CNalpha=CNALPHA, consts=consts,
                           n_per_cycle=NPC, n_cycles=NCYC)
idxs = np.linspace(0, len(o["alpha_deg"])-1, nph).astype(int)
xcp = np.linspace(0.02, 0.98, 60)      # common upper-surface x/c grid
Zsurf = []; phases = []
for ii in idxs:
    xoc, cp, upper = us.surface_cp(GEO, C_A, U_A, M_A, o["alpha_deg"][ii],
                                   o["CL"][ii], o["CNv"][ii], o["tau_v"][ii]/consts["Tvl"])
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
    # The cells immediately outside the body carry the regularised surface sheet
    # and reach ~3x the free-stream speed at peak incidence -- real for an
    # unseparated potential field, but they would set the whole z-scale. They are
    # clipped to robust percentiles, and the clip is STATED on the figure: the
    # flat-topped walls around the aerofoil hole are clipped values, not data.
    lo, hi = np.nanpercentile(Z, [0.5, 99.5])
    Zc = np.clip(Z, lo, hi)
    n_clipped = int(np.sum((Z < lo) | (Z > hi)))
    z_true = (float(np.nanmin(Z)), float(np.nanmax(Z)))
    # reserve a band at the bottom for the clip note, so it cannot land on the
    # x/y tick labels the way an axes-coordinate annotation did
    fig = plt.figure(figsize=(8.5, 6.4))
    ax = fig.add_axes([0.0, 0.12, 1.0, 0.84], projection="3d")
    # stride 1, NOT 2. plot_surface builds a quad from points (i, i+2) when
    # rstride=2, so a quad whose skipped middle point is masked is still drawn:
    # it bridges straight across the aerofoil and renders as a thin vertical fin
    # sticking out of the hole. At stride 1 the hole follows the mask exactly.
    s = ax.plot_surface(X, Y, Zc, cmap=cmap, vmin=lo, vmax=hi,
                        linewidth=0, antialiased=True, alpha=0.96,
                        rstride=1, cstride=1)
    ax.set_zlim(lo, hi)
    ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]"); ax.set_zlabel(lab)
    ax.set_title(f"3D field surface — {lab}  (Case A, peak incidence)", pad=18)
    # No colorbar here: height and colour encode the SAME variable, so a colorbar
    # put a second identical scale, with the same label, immediately beside the
    # z-axis. The z-axis is the scale.
    ax.view_init(elev=40, azim=-58); tidy3d(ax)
    fig.text(0.5, 0.015,
             "z clipped to the 0.5–99.5 percentile band [%.3g, %.3g] — %d of %d "
             "cells clipped, true range %.3g to %.3g.\nThe flat walls around the "
             "aerofoil hole are clipped near-surface cells, not a plateau in the field."
             % (lo, hi, n_clipped, Z.size - int(np.isnan(Z).sum()), z_true[0], z_true[1]),
             ha="center", va="bottom", fontsize=8.5, color=INK_SOFT)
    fig.savefig(OUT/f"fig3d_field_surface_{fn}.png", bbox_inches="tight", pad_inches=0.35); plt.close(fig)

# ====================================================== 4. PICTORIAL SECTION + VECTORS
xu, yu, F = load_field(peakA)
X, Y = np.meshgrid(xu, yu)
span = 0.45   # m, pseudo-span for the pictorial extrusion
fig = plt.figure(figsize=(9, 6.2)); ax = fig.add_subplot(111, projection="3d")
# extrude airfoil section at two span stations, filled so it reads as a solid
afx = AF["x_over_c"].values*C_A; afy = AF["y_over_c"].values*C_A
# fill BOTH end sections so the extrusion reads as a solid blade rather than
# two thin sticks lost among the arrows
for z in (0.0, span):
    ax.add_collection3d(Poly3DCollection(
        [list(zip(afx, np.full_like(afx, z), afy))], facecolor="#c6d7ea",
        edgecolor=INK, linewidths=1.6, alpha=1.0, zorder=10))
for frac in np.linspace(0, 1, 12):
    k = int(frac*(len(afx)-1))
    ax.plot([afx[k], afx[k]], [0, span], [afy[k], afy[k]], color=INK_SOFT, lw=0.5, alpha=0.6)
# 3-D velocity vectors on a mid-span plane.
# This figure was unreadable: the "window" that was supposed to restrict the
# arrows spanned -0.18..0.52 m in x and +-0.20 m in y, i.e. almost the whole
# 0.9 m x 0.72 m domain, and length=0.0026 against |V| ~ 150 m/s drew 0.39 m
# arrows -- 43 % of the axis -- so the arrows overlapped into a solid mat that
# buried the section. Both are now set FROM the geometry and the data: the
# window is a fixed number of chords around the section, and the arrow length is
# chosen so the longest arrow drawn is a set fraction of a chord.
WIN_X = (-0.35*C_A, 1.55*C_A)          # chords around the section
WIN_Y = 0.55*C_A
ARROW_MAX = 0.18*C_A                   # longest arrow drawn, in metres
sk = (slice(None, None, 9), slice(None, None, 9))
xs = X[sk]; ys = Y[sk]; us_ = F["u_ms"][sk]; vs = F["v_ms"][sk]
win = ((~np.isnan(us_)) & (xs > WIN_X[0]) & (xs < WIN_X[1]) & (np.abs(ys) < WIN_Y))
zc = span/2
uu, vv = us_[win].copy(), vs[win].copy()
mag = np.hypot(uu, vv)
# Cap the plotted magnitude: a handful of near-surface cells carry several times
# the freestream speed and would otherwise set the arrow scale for everything.
cap = float(np.nanpercentile(mag, 95))
scale = np.where(mag > cap, cap/np.maximum(mag, 1e-9), 1.0)
uu *= scale; vv *= scale
LEN = ARROW_MAX/cap                    # metres of arrow per m/s
ax.quiver(xs[win], np.full(win.sum(), zc), ys[win],
          uu, np.zeros(win.sum()), vv,
          length=LEN, normalize=False, color=PALETTE[0], linewidth=0.9,
          arrow_length_ratio=0.32, zorder=2, alpha=0.85)
# a scale bar, so the arrows mean something
ax.text2D(0.02, 0.03, "arrows: velocity on the mid-span plane, "
          "longest = %.0f m/s (%.2f c)" % (cap, ARROW_MAX/C_A),
          transform=ax.transAxes, fontsize=8.5, color=INK_SOFT)
ax.set_xlim(WIN_X[0], WIN_X[1]); ax.set_zlim(-WIN_Y, WIN_Y); ax.set_ylim(0, span)
ax.set_xlabel("x [m]"); ax.set_ylabel("span z [m]"); ax.set_zlabel("y [m]")
ax.set_title("Pictorial blade section with reconstructed velocity field (peak incidence)",
             pad=20)
ax.set_yticks(np.linspace(0, span, 4))
ax.view_init(elev=18, azim=-68); tidy3d(ax)
ax.set_box_aspect((1.4, 0.8, 0.7))
fig.savefig(OUT/"fig3d_section_vectors.png", bbox_inches="tight", pad_inches=0.35); plt.close(fig)

print("[3d] response surface + Cp-phase surface + field surfaces + pictorial section done")
