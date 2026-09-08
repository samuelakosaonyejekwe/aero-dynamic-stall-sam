"""
02_mesh / generate_mesh.py
--------------------------
Builds a body-fitted structured O-type grid around the NACA 0012 section used by
the UNISTALL(TM) field-reconstruction module, and reports mesh-quality metrics.

The reduced-order UIBS core does not require a volume mesh, but the universal
solver embeds a body-fitted grid for (a) panel-method field reconstruction and
(b) optional CFD hand-off. This module generates that grid and reports the usual
CFD quality metrics for it.

Quality caveat -- every number below is written to mesh_quality_metrics.csv by
this script, so it cannot drift away from the grid it describes. Read them
before reusing this grid for CFD. The worst cells sit near the trailing edge,
where the wrap turns through the largest angle. The loads reported by this study
do not depend on this grid -- the UIBS core is meshless and the field
reconstruction is panel-based -- so the grid is a documentation and hand-off
artefact, not a production CFD mesh.

Two defects that an unsigned quality metric cannot see were found and fixed here:
  * INVERTED CELLS. The previous grid contained 22 folded (negative-Jacobian)
    cells, in two single columns either side of the trailing-edge seam, spanning
    wall-normal layers 95-105. min_cell_area was computed with abs(), so a
    folded cell reported a healthy positive area and the fault was invisible in
    the published metrics. The cause was noise in the wall-normal directions:
    they were finite-differenced from a piecewise-LINEAR resample of the
    coordinate CSV at a trailing-edge spacing finer than the source data, which
    put a ~4e-4 staircase into the normal direction -- enough for neighbouring
    marching lines to cross several chords out. Fixed by taking the normals
    before the trailing-edge merge (the merge itself put a spurious near-normal
    first segment at i=1), setting the seam normal to the trailing-edge
    bisector, and applying one 3-point smoothing pass to the normal field.
    inverted_cells is now a published metric and is asserted to be zero.
  * BACKWARDS WALL CLUSTERING. The wall points were cosine-clustered in
    ARCLENGTH, whose two ends are both at the trailing edge, so the leading
    edge -- the one place a boundary-layer grid most needs resolution, and the
    place fig_mesh_le_zoom.png claims to show clustered -- carried the coarsest
    spacing on the body: 9.6x the first trailing-edge segment, and 53x the
    finest trailing-edge segment. The distribution is now cosine on each
    surface separately, clustering at the leading edge as well; the published
    wall_spacing_at_LE/TE metrics now stand at 0.12x.

Measured effect of the two fixes, old -> new (all rows of mesh_quality_metrics.csv):
    inverted cells               22    -> 0
    max aspect ratio             3546  -> 952
    cells with AR > 1000         0.111 % -> 0.000 %
    max skewness                 0.765 -> 0.642
    cells with skewness > 0.5    1.621 % -> 1.204 %
    min orthogonality            21.1 deg -> 32.3 deg

Outputs
  mesh_nodes.csv              every grid node (i, j, x_m, y_m, wall_distance_m)
  mesh_quality_metrics.csv    scalar quality metrics (y+, growth, AR, ortho, skew)
  mesh_radial_spacing.csv     wall-normal spacing law
  fig_mesh_full.png           full O-grid (far field)
  fig_mesh_le_zoom.png        leading-edge boundary-layer zoom
  fig_mesh_te_zoom.png        trailing-edge zoom
  fig_mesh_wall_spacing.png   first-cell height / growth-ratio plot

Topology note: the wall line wraps the complete surface and every node is offset
along a normal that is blended into a radial direction at the far field, so the
outer boundary is a circle and there is NO wake cut. That is an O-grid, not a
C-grid; the metrics CSV records it as such.
"""
import sys
import numpy as np
import pandas as pd
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from aero_style import apply_style, PALETTE, INK_SOFT
import matplotlib.pyplot as plt
apply_style()

# ---- case conditions: READ from 03_model_setup, never restated here. They used
#      to be a second hardcoded copy (rho=1.10, U=102.0, mu=1.78e-5, c=0.30) kept
#      "in step by hand"; the setup stage now derives rho from p/(R*T) and U from
#      the Mach test point, so the copies would have disagreed. This is why the
#      pipeline runs 03_model_setup BEFORE 02_mesh -- see run_all.py. ----
_flow = pd.read_csv(HERE.parent/"03_model_setup"/"flow_conditions.csv"
                    ).set_index("parameter")["case_A_validation"]
CHORD = float(_flow["chord_c"])

# ---- load airfoil surface ----
coords = pd.read_csv(HERE.parent/"01_geometry"/"naca0012_coordinates.csv")
xs, ys = coords["x_over_c"].values, coords["y_over_c"].values

# ---- wall line: resample the section to N_WALL points, clustered at BOTH the
#      leading and the trailing edge. A single cosine in arclength clusters at
#      s=0 and s=1, which are BOTH the trailing edge, and starves the leading
#      edge; see the header. Cosine on each surface separately puts refinement
#      where the boundary layer needs it. ----
N_WALL = 257
s = np.concatenate([[0], np.cumsum(np.hypot(np.diff(xs), np.diff(ys)))])
s /= s[-1]
_half = (1 - np.cos(np.linspace(0, np.pi, (N_WALL + 1)//2))) / 2   # dense at both ends
sq = np.concatenate([0.5*_half, 0.5 + 0.5*_half[1:]])              # TE -> LE -> TE
xw = np.interp(sq, s, xs)
yw = np.interp(sq, s, ys)

# ---- surface normals (outward). Taken on the OPEN wall, BEFORE the trailing
#      edge is merged: merging first moves the seam point by half the TE gap,
#      which makes the segment i=0 -> i=1 almost perpendicular to the surface and
#      hands i=1 a normal ~50 deg away from its neighbours'. ----
dx = np.gradient(xw); dy = np.gradient(yw)
nl = np.hypot(dx, dy); tx, ty = dx/nl, dy/nl
nx, ny = ty, -tx                     # rotate tangent -> normal
# flip normals that point inward (dot with radial from mid-chord)
inward = (nx*(xw - 0.5) + ny*yw) < 0
nx = np.where(inward, -nx, nx); ny = np.where(inward, -ny, ny)

# CLOSE the wall curve. The 4-digit section has an open trailing edge (0.252 %c),
# so the wrap starts and ends at two distinct points. Left open, the two branches
# march away from each other downstream and leave an unmeshed wedge several
# chords wide in the mid-field. Merging them makes the wall a genuinely closed
# curve, which is what an O-grid needs; the seam normal is the bisector of the
# two trailing-edge branch normals.
xm, ym = 0.5*(xw[0] + xw[-1]), 0.5*(yw[0] + yw[-1])
xw[0] = xw[-1] = xm
yw[0] = yw[-1] = ym
_bx, _by = nx[0] + nx[-1], ny[0] + ny[-1]; _bn = np.hypot(_bx, _by)
nx[0] = nx[-1] = _bx/_bn; ny[0] = ny[-1] = _by/_bn

# One 3-point smoothing pass over the closed normal field. np.interp resamples
# the coordinate CSV piecewise-LINEARLY, so where the requested wall spacing is
# comparable to the source spacing the normal direction picks up a staircase of
# order 4e-4. That is small, but it is enough to make two neighbouring marching
# lines cross several chords out: it is what produced the 22 folded cells in the
# previous grid. One pass removes it; the assertion below proves the result.
N_SMOOTH = 1
_a, _b = nx[:-1].copy(), ny[:-1].copy()            # unique points of the loop
for _ in range(N_SMOOTH):
    _a = 0.25*np.roll(_a, 1) + 0.5*_a + 0.25*np.roll(_a, -1)
    _b = 0.25*np.roll(_b, 1) + 0.5*_b + 0.25*np.roll(_b, -1)
    _n = np.hypot(_a, _b); _a /= _n; _b /= _n
nx = np.append(_a, _a[0]); ny = np.append(_b, _b[0])

# ---- wall-normal distribution (geometric growth to far field) ----
N_RAD = 121
FARFIELD = 20.0          # chords (radius of far-field boundary)
Y_PLUS_TARGET = 1.0
# estimate first-cell height for y+ ~= 1 at Re_c (case A conditions, read above)
rho = float(_flow["air_density_rho"]); U = float(_flow["freestream_velocity_U"])
mu  = float(_flow["dynamic_viscosity_mu"])
RE_C = rho*U*CHORD/mu                            # derived, never quoted
Cf = 0.026 / RE_C**(1/7.0)                       # turbulent flat-plate ~
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
# Each wall point is given its OWN far-field target, placed on the outer circle
# by arclength so the mapping wall -> circle is one-to-one and continuous. The
# marching direction blends from the wall normal (orthogonality at the wall) to
# the direction of that target. Blending towards a common radial direction from
# mid-chord, as before, sent the two trailing-edge branches apart and opened the
# wedge; blending towards a closed loop of targets closes the wrap at every j.
I, J = N_WALL, N_RAD
CX, CY = 0.5, 0.0
sw = np.concatenate([[0], np.cumsum(np.hypot(np.diff(xw), np.diff(yw)))]); sw /= sw[-1]
th = 2.0*np.pi*sw                                  # TE -> 0, LE -> pi, TE -> 2pi
xf = CX + FARFIELD*np.cos(th); yf = CY + FARFIELD*np.sin(th)
ufx, ufy = xf - xw, yf - yw
uf = np.hypot(ufx, ufy); ufx /= uf; ufy /= uf

Xg = np.zeros((I, J)); Yg = np.zeros((I, J))
for j in range(J):
    w = yn[j]/FARFIELD
    nxf = (1-w)*nx + w*ufx
    nyf = (1-w)*ny + w*ufy
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
            # SIGNED area. Taking abs() here is what hid 22 folded cells in the
            # previous grid: an inverted cell has the opposite orientation to
            # its neighbours but a perfectly healthy |area|.
            a = 0.5*(np.dot(p[:,0], np.roll(p[:,1],-1)) - np.dot(p[:,1], np.roll(p[:,0],-1)))
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

# ---- validity: no cell may be inverted relative to the dominant orientation ----
_orient = np.sign(np.median(area))
n_inverted = int(np.sum(area*_orient <= 0.0))
area_abs = np.abs(area)
# wall spacing, so the leading-edge clustering claim has a number behind it
_dwall = np.hypot(np.diff(xw), np.diff(yw))
_iLE = int(np.argmin(xw))
pct_skew_gt_05 = 100.0*float(np.mean(skew > 0.5))
pct_ar_gt_1000 = 100.0*float(np.mean(ar > 1000.0))

metrics = pd.DataFrame({
    "metric": ["topology", "i_nodes_wrap", "j_nodes_normal", "total_nodes",
               "total_cells", "farfield_radius_chords", "first_cell_height_y1_chords",
               "first_cell_height_y1_m", "wall_normal_growth_ratio", "target_yplus",
               "max_aspect_ratio", "mean_aspect_ratio", "pct_cells_aspect_ratio_gt_1000",
               "min_orthogonality_deg", "max_skewness", "mean_skewness",
               "pct_cells_skewness_gt_0.5", "min_cell_area_c2", "inverted_cells",
               "wall_spacing_at_LE_chords", "wall_spacing_at_TE_chords"],
    "value": ["O-grid (body-fitted, no wake cut)", I, J, I*J, (I-1)*(J-1), FARFIELD,
              round(y1,7), round(y1*CHORD,8), round(GR,4), Y_PLUS_TARGET,
              round(ar.max(),1), round(ar.mean(),1), round(pct_ar_gt_1000,3),
              round(ortho.min(),1), round(skew.max(),3), round(skew.mean(),3),
              round(pct_skew_gt_05,3), float("%.2e"%area_abs.min()), n_inverted,
              float("%.2e"%_dwall[_iLE]), float("%.2e"%_dwall[0])],
})
metrics.to_csv(HERE/"mesh_quality_metrics.csv", index=False)
if n_inverted:
    raise SystemExit("[mesh] FAILED: %d inverted cells — grid is not valid" % n_inverted)

# nodes csv (subsampled to keep file reasonable: every node)
ii, jj = np.meshgrid(np.arange(I), np.arange(J), indexing="ij")
wall_dist = np.repeat(yn[None,:], I, axis=0)*CHORD
nodes = pd.DataFrame({"i": ii.ravel(), "j": jj.ravel(),
                      "x_m": (Xg*CHORD).ravel().round(6),
                      "y_m": (Yg*CHORD).ravel().round(6),
                      "wall_distance_m": wall_dist.ravel().round(7)})
nodes.to_csv(HERE/"mesh_nodes.csv", index=False)

# spacing of layer j is yn[j]-yn[j-1]; layer 0 lies ON the wall, so it has none
layer_spacing = np.concatenate([[np.nan], np.diff(yn)])
pd.DataFrame({"layer_j": np.arange(N_RAD),
              "normal_coord_chords": yn.round(6),
              "layer_spacing_chords": layer_spacing.round(7),
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
ax.set_title("Body-fitted O-grid (near field) — %d×%d nodes" % (I, J))
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
l1, = ax.semilogy(layer, layer_spacing*CHORD*1e3, color=PALETTE[1], lw=2,
                  marker="o", ms=3, label="layer spacing")
ax.set_xlabel("wall-normal layer index j"); ax.set_ylabel("cell height [mm]")
ax.set_title("Wall-normal spacing law (geom. growth GR=%.3f, y1=%.2e m)" % (GR, y1*CHORD))
ax2 = ax.twinx(); ax2.grid(False)
l2, = ax2.plot(layer, yn*CHORD, color=PALETTE[2], lw=1.5, ls="--",
               label="cumulative normal distance")
ax2.set_ylabel("cumulative normal distance [m]", color=PALETTE[2])
# single combined legend: both curves live on different axes
ax.legend(handles=[l1, l2], loc="upper left")
fig.savefig(HERE/"fig_mesh_wall_spacing.png"); plt.close(fig)

print("[mesh] %d nodes, GR=%.3f, y1=%.2e m, maxAR=%.0f, maxSkew=%.3f, minOrtho=%.1f deg, "
      "inverted=%d, wall spacing LE/TE=%.2f"
      % (I*J, GR, y1*CHORD, ar.max(), skew.max(), ortho.min(), n_inverted,
         _dwall[_iLE]/_dwall[0]))
