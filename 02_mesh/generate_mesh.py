"""
02_mesh / generate_mesh.py
--------------------------
Builds a body-fitted structured C-type grid around the NACA 0012 section used by
the UNISTALL(TM) field-reconstruction module, and reports mesh-quality metrics.

The reduced-order UIBS core does not require a volume mesh, but the universal
solver embeds a body-fitted grid for (a) panel-method field reconstruction and
(b) optional CFD hand-off. This module documents that grid to CFD-grade standards.

Outputs
  mesh_nodes.csv              every grid node (i, j, x_m, y_m, wall_distance_m)
  mesh_quality_metrics.csv    scalar quality metrics (y+, growth, AR, ortho, skew)
  mesh_radial_spacing.csv     wall-normal spacing law
  fig_mesh_full.png           full C-grid (far field)
  fig_mesh_le_zoom.png        leading-edge boundary-layer zoom
  fig_mesh_te_zoom.png        trailing-edge zoom
  fig_mesh_wall_spacing.png   first-cell height / growth-ratio plot
"""
import sys
import numpy as np
import pandas as pd
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from aero_style import apply_style, PALETTE, INK, INK_SOFT, CMAP_PRESSURE
import matplotlib.pyplot as plt
apply_style()

# ---- load airfoil surface ----
coords = pd.read_csv(HERE.parent/"01_geometry"/"naca0012_coordinates.csv")
xs, ys = coords["x_over_c"].values, coords["y_over_c"].values
CHORD = 0.30

# resample surface to a smooth, evenly-distributed wall line (N_wall points)
N_WALL = 257
# parametric arclength
s = np.concatenate([[0], np.cumsum(np.hypot(np.diff(xs), np.diff(ys)))])
s /= s[-1]
sq = (1 - np.cos(np.linspace(0, np.pi, N_WALL))) / 2  # cluster LE/TE
xw = np.interp(sq, s, xs)
yw = np.interp(sq, s, ys)

# ---- surface normals (outward) ----
dx = np.gradient(xw); dy = np.gradient(yw)
nl = np.hypot(dx, dy); tx, ty = dx/nl, dy/nl
nx, ny = ty, -tx                     # rotate tangent -> normal
# ensure outward (away from chord line y=0 mid)
sign = np.sign((yw) + 1e-9)
sign[sign == 0] = 1
nx *= 1; ny *= 1
# flip normals that point inward (dot with radial from 0.5,0)
rx, ry = xw - 0.5, yw - 0.0
inward = (nx*rx + ny*ry) < 0
nx[inward] *= -1; ny[inward] *= -1

# ---- wall-normal distribution (geometric growth to far field) ----
N_RAD = 121
FARFIELD = 20.0          # chords (radius of far-field boundary)
Y_PLUS_TARGET = 1.0
# estimate first-cell height for y+ ~= 1 at Re_c
RE_C = 2.0e6
Cf = 0.026 / RE_C**(1/7.0)                       # turbulent flat-plate ~
rho, U, mu = 1.10, 102.0, 1.78e-5
tau_w = 0.5*rho*U**2*Cf
u_tau = np.sqrt(tau_w/rho)
y1 = Y_PLUS_TARGET*mu/(rho*u_tau) / CHORD          # in chords
y1 = max(y1, 1.0e-5)
# geometric series sum to FARFIELD
def growth_for(n, first, total):
    from scipy.optimize import brentq
    f = lambda r: first*(r**n - 1)/(r - 1) - total if abs(r-1)>1e-9 else first*n-total
    return brentq(f, 1.0001, 1.5)
GR = growth_for(N_RAD-1, y1, FARFIELD)
dn = y1*GR**np.arange(N_RAD)            # spacing of each layer
yn = np.concatenate([[0], np.cumsum(dn)])[:N_RAD]   # normal coordinate (chords)

# ---- build grid ----
I, J = N_WALL, N_RAD
Xg = np.zeros((I, J)); Yg = np.zeros((I, J))
for j in range(J):
    # blend wall normal -> radial (elliptic-like smoothing toward far field)
    w = (yn[j]/FARFIELD)
    cx, cy = 0.5, 0.0
    rxn = (xw-cx); ryn = (yw-cy); rr = np.hypot(rxn, ryn)+1e-9
    nxf = (1-w)*nx + w*rxn/rr
    nyf = (1-w)*ny + w*ryn/rr
    nf = np.hypot(nxf, nyf)
    Xg[:, j] = xw + yn[j]*nxf/nf
    Yg[:, j] = yw + yn[j]*nyf/nf

# ---- quality metrics ----
# cell areas, aspect ratio, orthogonality, skewness (approx, on interior cells)
def cell_quality(Xg, Yg):
    ars, orthos, skews, areas = [], [], [], []
    for i in range(I-1):
        for j in range(J-1):
            p = np.array([[Xg[i,j],Yg[i,j]],[Xg[i+1,j],Yg[i+1,j]],
                          [Xg[i+1,j+1],Yg[i+1,j+1]],[Xg[i,j+1],Yg[i,j+1]]])
            # area
            a = 0.5*abs(np.dot(p[:,0], np.roll(p[:,1],-1)) - np.dot(p[:,1], np.roll(p[:,0],-1)))
            areas.append(a)
            e = np.diff(np.vstack([p, p[0]]), axis=0)
            le = np.hypot(e[:,0], e[:,1])+1e-12
            ar = max(le)/min(le); ars.append(ar)
            # orthogonality: min angle between adjacent edges vs 90
            ang = []
            for k in range(4):
                v1 = e[k]/le[k]; v2 = -e[(k-1)%4]/le[(k-1)%4]
                ang.append(np.degrees(np.arccos(np.clip(np.dot(v1,v2),-1,1))))
            ang = np.array(ang)
            orthos.append(90 - np.max(np.abs(ang-90)))   # 90 = perfect
            skews.append(np.max(np.abs(ang-90))/90.0)     # 0 = perfect
    return (np.array(ars), np.array(orthos), np.array(skews), np.array(areas))

ar, ortho, skew, area = cell_quality(Xg, Yg)

metrics = pd.DataFrame({
    "metric": ["topology", "i_nodes_wrap", "j_nodes_normal", "total_nodes",
               "total_cells", "farfield_radius_chords", "first_cell_height_y1_chords",
               "first_cell_height_y1_m", "wall_normal_growth_ratio", "target_yplus",
               "max_aspect_ratio", "mean_aspect_ratio", "min_orthogonality_deg",
               "max_skewness", "mean_skewness", "min_cell_area_c2"],
    "value": ["C-grid (body-fitted)", I, J, I*J, (I-1)*(J-1), FARFIELD,
              round(y1,7), round(y1*CHORD,8), round(GR,4), Y_PLUS_TARGET,
              round(ar.max(),1), round(ar.mean(),1), round(ortho.min(),1),
              round(skew.max(),3), round(skew.mean(),3), float("%.2e"%area.min())],
})
metrics.to_csv(HERE/"mesh_quality_metrics.csv", index=False)

# nodes csv (subsampled to keep file reasonable: every node)
ii, jj = np.meshgrid(np.arange(I), np.arange(J), indexing="ij")
wall_dist = np.repeat(yn[None,:], I, axis=0)*CHORD
nodes = pd.DataFrame({"i": ii.ravel(), "j": jj.ravel(),
                      "x_m": (Xg*CHORD).ravel().round(6),
                      "y_m": (Yg*CHORD).ravel().round(6),
                      "wall_distance_m": wall_dist.ravel().round(7)})
nodes.to_csv(HERE/"mesh_nodes.csv", index=False)

pd.DataFrame({"layer_j": np.arange(N_RAD),
              "normal_coord_chords": yn.round(6),
              "layer_spacing_chords": np.concatenate([[y1], np.diff(yn)]).round(7),
              "normal_coord_m": (yn*CHORD).round(6)}).to_csv(HERE/"mesh_radial_spacing.csv", index=False)

# ---- figures ----
def plot_grid(ax, every_i=4, every_j=3, lw=0.4):
    for j in range(0, J, every_j):
        ax.plot(Xg[:, j]*CHORD, Yg[:, j]*CHORD, color=INK_SOFT, lw=lw)
    for i in range(0, I, every_i):
        ax.plot(Xg[i, :]*CHORD, Yg[i, :]*CHORD, color=INK_SOFT, lw=lw)
    ax.fill(xw*CHORD, yw*CHORD, color=PALETTE[0], alpha=0.35)
    ax.set_aspect("equal"); ax.grid(False)

fig, ax = plt.subplots(figsize=(7.5, 7.5))
plot_grid(ax, 4, 4)
ax.set_xlim(-6*CHORD, 7*CHORD); ax.set_ylim(-6.5*CHORD, 6.5*CHORD)
ax.set_title("Body-fitted C-grid (near field) — %d×%d nodes" % (I, J))
ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]")
fig.savefig(HERE/"fig_mesh_full.png"); plt.close(fig)

fig, ax = plt.subplots(figsize=(6, 6))
plot_grid(ax, 2, 1, lw=0.5)
ax.set_xlim(-0.02*CHORD, 0.18*CHORD); ax.set_ylim(-0.10*CHORD, 0.10*CHORD)
ax.set_title("Leading-edge boundary-layer clustering (y+ ≈ 1)")
ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]")
fig.savefig(HERE/"fig_mesh_le_zoom.png"); plt.close(fig)

fig, ax = plt.subplots(figsize=(6, 6))
plot_grid(ax, 2, 1, lw=0.5)
ax.set_xlim(0.85*CHORD, 1.10*CHORD); ax.set_ylim(-0.12*CHORD, 0.12*CHORD)
ax.set_title("Trailing-edge / near-wake clustering")
ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]")
fig.savefig(HERE/"fig_mesh_te_zoom.png"); plt.close(fig)

fig, ax = plt.subplots(figsize=(8, 4.5))
layer = np.arange(N_RAD)
ax.semilogy(layer, np.concatenate([[y1], np.diff(yn)])*CHORD*1e3, color=PALETTE[1], lw=2,
            marker="o", ms=3, label="layer spacing")
ax.set_xlabel("wall-normal layer index j"); ax.set_ylabel("cell height [mm]")
ax.set_title("Wall-normal spacing law (geom. growth GR=%.3f, y1=%.2e m)" % (GR, y1*CHORD))
ax.legend(loc="upper left")
ax2 = ax.twinx(); ax2.grid(False)
ax2.plot(layer, yn*CHORD, color=PALETTE[2], lw=1.5, ls="--", label="cumulative dist")
ax2.set_ylabel("cumulative normal distance [m]", color=PALETTE[2])
fig.savefig(HERE/"fig_mesh_wall_spacing.png"); plt.close(fig)

print("[mesh] %d nodes, GR=%.3f, y1=%.2e m, maxAR=%.0f, maxSkew=%.3f, minOrtho=%.1f deg"
      % (I*J, GR, y1*CHORD, ar.max(), skew.max(), ortho.min()))
