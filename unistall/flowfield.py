# -*- coding: utf-8 -*-
"""
unistall / flowfield.py
-------------------
A potential-flow RECONSTRUCTION of the flow around the section, drawn around
the lift the load model predicts. It is not a flow solution.

Author: Akosa Samuel Onyejekwe (independent)

What it is
  * the closed section as straight panels, cosine-clustered at both edges,
    each carrying a source of constant strength and a vortex sheet of one
    uniform strength; the sheet carries the circulation Gamma = 0.5 C_L U c of
    the lift it is given. The velocity of every panel is taken from the
    closed-form integrals of the constant-strength panel, the same ones in
    the wall condition, on the surface and in the field;
  * a transpiration through the suction surface standing for the separated
    region. The predicted lift is not, in general, the lift at which
    potential flow leaves the trailing edge smoothly; in stall it is far below
    it. The wall condition on the suction-surface panels aft of the load
    model's separation point x/c = f is therefore v_n = v_t w, with w rising
    smoothly from 0 at the separation point to 1 at the trailing edge, and
    the one number v_t is solved so that the pressure is the same at the two
    control points next to the trailing edge (the Kutta condition) while the
    circulation stays Gamma. v_t takes either sign: positive is an outflow,
    negative an inflow (suction into the wall), which is what the condition
    asks for when the predicted lift is above the Kutta value. Where the
    model's flow is attached (f near 1) the transpiration is confined to the
    last ATTACHED_OUTFLOW_FRACTION of the suction surface. It is an
    EQUIVALENT transpiration chosen to satisfy the Kutta condition. It
    displaces the outer flow; it is not a model of the flow inside a
    separated region, and it draws no reversed flow, no shear layer and no
    wake;
  * one Lamb-Oseen vortex, a MARKER of the load model's vortex clock: its
    circulation is 0.5 C_N^v U c, its centre is at
    x/c = VORTEX_X0 + VORTEX_DX tau, y/c = VORTEX_Y0 + VORTEX_DY tau with
    tau = tau_v / T_VL, so that it is over the trailing edge when the clock
    says the vortex leaves the chord (tau = 1) and moves on downstream after
    it, and its core radius is VORTEX_CORE_OVER_C chords. It is included in
    the wall condition. It does not reproduce the suction footprint of a
    dynamic-stall vortex: a clockwise vortex above the surface slows the
    surface flow beneath it, so the pressure under it is raised unless the
    vortex is strong and close enough to reverse that flow at more than the
    local speed (`vortex_footprint` gives the value at each instant);
  * velocity, the incompressible pressure coefficient (steady Bernoulli,
    plus the radial-equilibrium correction inside the vortex core) and
    vorticity, at the nodes of a grid or at any points handed in;
  * a compressibility correction of the pressure. For a stream Mach number
    M > 0 the incompressible pressure coefficient is corrected by the
    Karman-Tsien rule (`karman_tsien`), and the sheet is given the fraction
    of the circulation of the lift at which the corrected surface pressure
    returns the lift given (`solve`; the fraction starts from the
    Prandtl-Glauert factor sqrt(1 - M^2), and `surface_load_closure` reports
    the closure). From the corrected pressure
    the local Mach number, the static temperature and the recovery
    temperature follow by the isentropic relations with the stagnation state
    of the stream (`compressible_state`).

What it is not
  * it is quasi-steady: each instant is a steady flow at that instant's
    incidence and lift. There is no shed wake, no pitch-rate term in the wall
    condition, no unsteady term in the pressure, and the vortex is not
    force-free (it is held where the clock puts it and does not move with
    the flow);
  * it is not a compressible solution. The Karman-Tsien rule is a
    linearised correction applied to an incompressible reconstruction and
    holds only where the local flow is subcritical. Where the corrected
    pressure coefficient is below the critical value C_p* of the stream Mach
    number (`sonic_values`) the local Mach number it gives is above one, no
    shock or supersonic region is represented, and the pressure, Mach number
    and temperatures there are not physical: `beyond_critical` marks those
    points, `critical_coverage` measures them, and the figures hatch them.
    The corrected pressure is not allowed below VACUUM_CAP of the vacuum
    value -2/(gamma M^2). The recovery temperature is that of an adiabatic
    wall under a turbulent boundary layer (recovery factor Pr^(1/3)); the
    reconstruction has no boundary layer, and away from the wall the value
    is only what such a wall would take at the local state. `limit_to_sonic`
    and `thermal_state` are older relations kept for a caller who wants
    them; the reconstruction does not use them unless asked (`limit=True`);
  * it has no boundary layer and no separated shear layer, and its vorticity
    is zero everywhere except in the core of the vortex (the bound sheet
    lies on the surface itself);
  * it is given the lift only. With the transpiration and the vortex the
    surface pressure does not integrate exactly to that lift, and it does not
    return the model's moment or drag; `surface_load_closure` reports all
    three;
  * the constants of the vortex path and core and the shape of the
    transpiration are chosen for illustration, not derived or measured;
  * none of it is compared with a measurement, and no load depends on it: the
    loads come from unistall.dsmodel and are only read here.

Units: lengths in metres, speeds in m/s, angles in degrees at the interface
(radians inside), temperatures in kelvin, vorticity in 1/s, circulation in
m^2/s. Coefficients are dimensionless. Positive circulation is clockwise
(positive lift for a stream in +x), for the sheet and the vortex alike.
"""
import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.path import Path as _Polygon
from scipy.interpolate import CubicSpline

N_PANELS = 320                       # panels around the closed section
NEAR_WALL_PANEL_LENGTHS = 0.6        # distance from the wall, in mean panel lengths, within which the field
                                     # shows the junctions of the panels
DOMAIN_CHORDS = (-1.0, 2.0, -1.2, 1.2)   # field window (x0, x1, y0, y1), in chords
GRID_NX, GRID_NY = 180, 140          # nodes of the field grid
NEAR_WALL_RINGS_BLANKED = 0          # rings of grid nodes blanked around the section
EVALUATION_BLOCK = 1000              # field points evaluated at a time
FAR_PANEL_LENGTHS = 4.0              # beyond this many of its lengths a panel is taken by its two-term expansion

# Transpiration standing for the separated region.
ATTACHED_OUTFLOW_FRACTION = 0.05     # chord fraction ahead of the trailing edge that always carries it

# The vortex, a marker of the load model's vortex clock. Path of its centre:
#   x_v / c = VORTEX_X0 + VORTEX_DX * tau ,  y_v / c = VORTEX_Y0 + VORTEX_DY * tau ,  tau = tau_v / T_VL,
# so that x_v / c = 1 at tau = 1. Core radius VORTEX_CORE_OVER_C chords. All chosen for illustration.
VORTEX_X0, VORTEX_DX = 0.25, 0.75
VORTEX_Y0, VORTEX_DY = 0.10, 0.06
VORTEX_CORE_OVER_C = 0.05
CORE_INTEGRAL_POINTS = 800               # radial stations of the core-pressure integral

# Peak-swirl coefficient of the Lamb-Oseen profile, max over s of (1 - exp(-s^2))/s.
_s = np.linspace(1e-6, 5.0, 200001)
LAMB_OSEEN_PEAK = float(np.max((1.0 - np.exp(-_s*_s))/_s))
del _s

# Sutherland's law for air: reference viscosity [Pa s], reference temperature [K], constant [K]
SUTHERLAND_MU0, SUTHERLAND_T0, SUTHERLAND_S = 1.716e-5, 273.15, 110.4

# What the reconstruction is, in the words every figure note, caption and document uses.
QUASI_STEADY_OMISSIONS = ("no shed wake", "no pitch-rate boundary condition", "no unsteady pressure term",
                          "vortex not force-free")
STANDING_NOTE = ("Potential-flow reconstruction drawn round the predicted lift, not a flow solution: incompressible, "
                 "with a linearised (Karman-Tsien) correction of the pressure that holds only where the local flow "
                 "is subcritical, and quasi-steady (" + ", ".join(QUASI_STEADY_OMISSIONS) + "); no boundary layer, "
                 "no separated shear layer. A transpiration through the suction surface, of either sign, stands for the "
                 "separated region, which is not drawn: no reversed flow and no wake appear. The vortex is a "
                 "marker of the model's vortex clock and does not reproduce a suction footprint. Compared with no "
                 "measurement; no load depends on it.")
AXES_NOTE = "Axes fixed to the section; the stream comes from below at the incidence shown."
CRITICAL_NOTE = ("Hatched: the corrected pressure coefficient is below the critical value (local Mach number above "
                 "one); the correction does not hold there and the values are not physical.")
VACUUM_CAP = 0.95                    # the corrected Cp is not allowed below this fraction of the vacuum value
RECOVERY_SOURCE = "Pr^(1/3), the usual value for a turbulent boundary layer on a flat plate, with Pr = 0.72 for air"


@dataclass(frozen=True)
class Air:
    """Sea-level standard air. gamma [-], R_gas [J/kg/K], prandtl [-],
    T_inf [K], p_inf [Pa]. The speed of sound follows from the first, second
    and fourth."""
    gamma: float = 1.4
    R_gas: float = 287.05
    prandtl: float = 0.72
    T_inf: float = 288.15
    p_inf: float = 101325.0

    @property
    def a_sound(self) -> float:
        """Speed of sound [m/s], sqrt(gamma R T)."""
        return math.sqrt(self.gamma*self.R_gas*self.T_inf)

    @property
    def cp(self) -> float:
        """Specific heat at constant pressure [J/kg/K]."""
        return self.gamma*self.R_gas/(self.gamma - 1.0)

    @property
    def recovery(self) -> float:
        """Recovery factor of a turbulent boundary layer, Pr^(1/3) [-]."""
        return self.prandtl**(1.0/3.0)

    @property
    def rho(self) -> float:
        """Density [kg/m^3] from the gas law."""
        return self.p_inf/(self.R_gas*self.T_inf)

    @property
    def mu(self) -> float:
        """Dynamic viscosity [Pa s] from Sutherland's law."""
        return (SUTHERLAND_MU0*(self.T_inf/SUTHERLAND_T0)**1.5
                * (SUTHERLAND_T0 + SUTHERLAND_S)/(self.T_inf + SUTHERLAND_S))


AIR = Air()


@dataclass(frozen=True)
class _Panels:
    """The panels and their influence at the control points, just outside the
    surface. W is the complex influence matrix of `_influence` there, and
    `solve` the inverse of the wall-condition matrix (normal velocity per
    unit source strength), formed once."""
    xp: np.ndarray
    yp: np.ndarray
    xc: np.ndarray
    yc: np.ndarray
    L: np.ndarray
    tx: np.ndarray
    ty: np.ndarray
    nx: np.ndarray
    ny: np.ndarray
    W: np.ndarray | None
    solve: np.ndarray | None


@dataclass(frozen=True)
class Reconstruction:
    """One instant of the reconstruction. Panel end points (xp, yp) [m],
    control points (xc, yc) [m], lengths L [m], source strengths sigma [m/s],
    sheet strength gam [m/s]; stream speed U [m/s], incidence alpha [rad],
    chord [m], the lift coefficient CL it carries; vortex centre (xv, yv) [m],
    circulation Gv [m^2/s], core radius rc [m] and the speed Ve [m/s] of the
    flow at the vortex centre without the vortex; transpiration speed vt
    [m/s] at the trailing edge (negative: inflow), its shape w [-] on each
    panel (zero where there is none) and the chord fraction x_sep at which
    it starts; `panels`, the influence matrices the instant was solved with;
    `factor`, the fraction of the circulation of the lift CL that the sheet
    carries (1, or the value solved for a stream Mach number, which is then
    kept in `mach`; see `solve`)."""
    xp: np.ndarray
    yp: np.ndarray
    xc: np.ndarray
    yc: np.ndarray
    L: np.ndarray
    sigma: np.ndarray
    gam: float
    U: float
    alpha: float
    chord: float
    CL: float
    xv: float
    yv: float
    Gv: float
    rc: float
    Ve: float
    vt: float
    w: np.ndarray
    x_sep: float
    panels: _Panels | None = field(default=None, repr=False, compare=False)
    factor: float = 1.0
    mach: float = 0.0


# --------------------------------------------------------------------------- #
#  geometry
# --------------------------------------------------------------------------- #
def panel_geometry(coords_csv: Path, chord: float, n_panel: int = N_PANELS) -> tuple:
    """Panel end points (xp, yp) [m] around the CLOSED section, for chord [m].

    `coords_csv` has columns x_over_c, y_over_c running trailing edge, upper
    surface, leading edge, lower surface, trailing edge, the first and last
    point being the one trailing-edge point. The end points lie on a cubic
    spline through those coordinates, parametrised by the distance along the
    outline, and are cosine-spaced in that distance on each surface
    separately, so they cluster at both edges and their number can exceed
    the number of coordinates."""
    df = pd.read_csv(coords_csv)
    x, y = df["x_over_c"].to_numpy(float), df["y_over_c"].to_numpy(float)
    if abs(x[0] - x[-1]) > 1e-9 or abs(y[0] - y[-1]) > 1e-9:
        raise ValueError(f"{coords_csv}: the outline is not closed at the trailing edge")
    s = np.concatenate([[0.0], np.cumsum(np.hypot(np.diff(x), np.diff(y)))])
    s_le = s[int(np.argmin(x))]
    h = (1.0 - np.cos(np.linspace(0.0, np.pi, n_panel//2 + 1)))/2.0
    sq = np.concatenate([s_le*h, s_le + (s[-1] - s_le)*h[1:]])
    xp, yp = CubicSpline(s, x)(sq)*chord, CubicSpline(s, y)(sq)*chord
    xp[-1], yp[-1] = xp[0], yp[0]
    return xp, yp


def _influence(x, y, P):
    """Complex influence of the panels at points (x, y) [m], shape (points,
    panels), from the closed-form integral of a constant-strength panel:
    with z = x + i y and a panel of length L from z1 to z2 at angle theta,

        W = exp(-i theta) ln((z - z1)/(z - z2)) / (2 pi),

    and a panel carrying a source of strength sigma and a clockwise vortex
    sheet of strength gamma induces the velocity u - i v = (sigma + i gamma) W.
    Further from a panel's mid-point zc than FAR_PANEL_LENGTHS of its lengths
    the logarithm is replaced by the first two terms of its expansion in
    L/(z - zc), which there differ from it by less than one part in 10^4."""
    z = (x + 1j*y)[:, None]
    d = z - (P.xc + 1j*P.yc)[None, :]
    d[d == 0.0] = 1.0
    inv = np.reciprocal(d, out=d)
    W = inv*(P.L*(P.tx + 1j*P.ty))[None, :]                   # (z2 - z1)/(z - zc)
    near = np.nonzero(W.real**2 + W.imag**2 > 1.0/FAR_PANEL_LENGTHS**2)
    W *= W
    W *= 1.0/12.0
    W += 1.0
    W *= inv
    W *= (P.L/(2.0*np.pi))[None, :]
    if len(near[0]):
        j = near[1]
        d1 = z[near[0], 0] - (P.xp[j] + 1j*P.yp[j])
        d2 = z[near[0], 0] - (P.xp[j + 1] + 1j*P.yp[j + 1])
        tiny = 1e-12*P.L.min()
        d1, d2 = np.where(d1 == 0.0, tiny, d1), np.where(d2 == 0.0, tiny, d2)
        W[near] = np.log(d1/d2)*(P.tx[j] - 1j*P.ty[j])/(2.0*np.pi)
    return W


def _build_panels(xp, yp):
    """The panels with end points (xp, yp) and their influence at the control
    points. A panel's own velocity there is the limit from outside: half the
    source strength along the outward normal (and with it half the sheet
    strength against the direction a counter-clockwise outline runs in)."""
    xc, yc = 0.5*(xp[:-1] + xp[1:]), 0.5*(yp[:-1] + yp[1:])
    L = np.hypot(np.diff(xp), np.diff(yp))
    tx, ty = np.diff(xp)/L, np.diff(yp)/L
    sense = 1.0 if 0.5*np.sum(xp[:-1]*yp[1:] - xp[1:]*yp[:-1]) > 0.0 else -1.0
    nx, ny = sense*ty, -sense*tx
    W = _influence(xc, yc, _Panels(xp, yp, xc, yc, L, tx, ty, nx, ny, None, None))
    i = np.arange(len(xc))
    W[i, i] = 0.5*(nx - 1j*ny)
    A = np.real(W*(nx + 1j*ny)[:, None])
    return _Panels(xp, yp, xc, yc, L, tx, ty, nx, ny, W, np.linalg.inv(A))


_PANEL_CACHE: dict = {}
_PANEL_CACHE_SIZE = 16


def _panels_of(coords_csv: Path, chord: float, n_panel: int) -> _Panels:
    """The panels of one coordinate file, chord and panel count, built once.
    The key holds the file's modification time, so an edited file is read
    again."""
    path = Path(coords_csv).resolve()
    key = (str(path), path.stat().st_mtime_ns, float(chord), int(n_panel))
    if key not in _PANEL_CACHE:
        if len(_PANEL_CACHE) >= _PANEL_CACHE_SIZE:
            _PANEL_CACHE.clear()
        _PANEL_CACHE[key] = _build_panels(*panel_geometry(path, float(chord), int(n_panel)))
    return _PANEL_CACHE[key]


# --------------------------------------------------------------------------- #
#  transpiration standing for the separated region
# --------------------------------------------------------------------------- #
def outflow_shape(xc_over_c: np.ndarray, upper: np.ndarray, f_sep: float, CL: float) -> tuple:
    """(w, x_sep): the shape of the transpiration on each panel and the chord
    fraction at which it starts.

    xc_over_c are the control points as fractions of chord, `upper` marks the
    upper-surface panels, f_sep is the load model's separation point (fraction
    of chord from the leading edge, 1 = attached) and CL the lift coefficient.
    The transpiration is on the suction surface: the upper one for positive
    lift, the lower one for negative. It starts at x_sep = min(f_sep,
    1 - ATTACHED_OUTFLOW_FRACTION) and w rises from 0 there to 1 at the
    trailing edge as a smooth step."""
    x_sep = float(min(max(f_sep, 0.0), 1.0 - ATTACHED_OUTFLOW_FRACTION))
    s = np.clip((np.asarray(xc_over_c, float) - x_sep)/(1.0 - x_sep), 0.0, 1.0)
    suction = upper if CL >= 0.0 else ~upper
    return np.where(suction, s*s*(3.0 - 2.0*s), 0.0), x_sep


def _edge_velocity(u, v):
    """The surface velocity (u, v) at the two control points next to the
    trailing edge, as a 2 x 2 array: first row the first panel, second row
    the last."""
    return np.array([[u[0], v[0]], [u[-1], v[-1]]])


def _kutta_strength(P, base, unit):
    """The transpiration strength vt at which the pressure is the same at
    the two control points next to the trailing edge. `base` and `unit` are
    the edge velocities (see `_edge_velocity`) without transpiration and per
    unit vt. Equal pressure is equal speed, a quadratic in vt; of its roots
    the one taken is the nearer to the strength at which the two tangential
    speeds are equal and opposite, the flow leaving both surfaces towards
    the edge. Without a real root that strength itself is returned."""
    t_first, t_last = np.array([P.tx[0], P.ty[0]]), np.array([P.tx[-1], P.ty[-1]])
    q0, q1 = base[0] @ t_first + base[1] @ t_last, unit[0] @ t_first + unit[1] @ t_last
    linear = -q0/q1
    a = unit[0] @ unit[0] - unit[1] @ unit[1]
    b = 2.0*(base[0] @ unit[0] - base[1] @ unit[1])
    c = base[0] @ base[0] - base[1] @ base[1]
    if abs(a) < 1e-14*max(abs(b), 1e-300):
        return -c/b
    disc = b*b - 4.0*a*c
    if disc < 0.0:
        return linear
    roots = ((-b + math.sqrt(disc))/(2.0*a), (-b - math.sqrt(disc))/(2.0*a))
    return min(roots, key=lambda r: abs(r - linear))


def _surface_velocity(P, U, alpha, gam, sigma, extra=None):
    """Velocity (u, v) at the control points, just outside the surface."""
    eu, ev = extra if extra is not None else (0.0, 0.0)
    w = P.W @ (sigma + 1j*gam)
    return U*np.cos(alpha) + eu + w.real, U*np.sin(alpha) + ev - w.imag


def _solve_sources(P, U, alpha, gam, w=None, extra=None):
    """Source strengths for the wall condition v_n = vt w, in the presence of
    the sheet of strength gam and of `extra`, a (u, v) pair induced at the
    control points by anything else. With w given, vt is the value that
    satisfies the Kutta condition; without it the wall is solid. Returns
    (sigma, vt)."""
    u0, v0 = _surface_velocity(P, U, alpha, gam, np.zeros_like(P.xc), extra)
    sigma = P.solve @ -(u0*P.nx + v0*P.ny)
    if w is None or not np.any(w):
        return sigma, 0.0
    unit = P.solve @ w                                        # sources of a unit transpiration
    ws, wu = P.W @ sigma, P.W @ unit
    vt = _kutta_strength(P, _edge_velocity(u0 + ws.real, v0 - ws.imag), _edge_velocity(wu.real, -wu.imag))
    return sigma + vt*unit, float(vt)


# --------------------------------------------------------------------------- #
#  induced velocities
# --------------------------------------------------------------------------- #
def _panel_velocity(x, y, P, sigma, gam):
    """Velocity (u, v) induced at points (x, y) [m], arrays of one shape, by
    the sources and the sheet of the panels P."""
    xf, yf = np.ravel(x), np.ravel(y)
    w = np.empty(len(xf), complex)
    strength = sigma + 1j*gam
    for i in range(0, len(xf), EVALUATION_BLOCK):
        part = slice(i, i + EVALUATION_BLOCK)
        w[part] = _influence(xf[part], yf[part], P) @ strength
    return w.real.reshape(np.shape(x)), -w.imag.reshape(np.shape(x))


def _vortex_velocity(x, y, xv, yv, Gv, rc):
    """Lamb-Oseen vortex: (u, v, r^2) at points (x, y)."""
    rx, ry = x - xv, y - yv
    r2 = rx*rx + ry*ry
    safe = np.where(r2 == 0.0, 1.0, r2)
    core = (1.0 - np.exp(-r2/rc**2))*Gv/(2.0*np.pi*safe)
    return core*ry, -core*rx, r2


def vortex_centre(tau_over_Tvl: float) -> tuple:
    """(x/c, y/c) of the vortex centre at vortex clock tau_v / T_VL: on the
    straight path of the constants above, over the trailing edge at 1."""
    t = max(float(tau_over_Tvl), 0.0)
    return VORTEX_X0 + VORTEX_DX*t, VORTEX_Y0 + VORTEX_DY*t


CLOSURE_TOL, CLOSURE_ITERATIONS = 1e-6, 12        # the circulation that returns the lift: tolerance in CL, steps


def _instant(P, chord, U, alpha, CL, CN_vortex, tau_over_Tvl, f_sep, factor, mach):
    """One panel solution with the sheet carrying `factor` times the
    circulation of CL."""
    gam = float(factor*0.5*CL*U*chord/P.L.sum())
    upper = np.arange(len(P.xc)) < len(P.xc)//2
    w, x_sep = outflow_shape(P.xc/chord, upper, f_sep, CL) if f_sep is not None else (np.zeros_like(P.xc), 1.0)
    xv, yv = (chord*q for q in vortex_centre(tau_over_Tvl))
    rc = VORTEX_CORE_OVER_C*chord
    Gv = 0.5*CN_vortex*U*chord if (CN_vortex > 0.0 and tau_over_Tvl > 0.0 and U > 0.0) else 0.0
    extra = _vortex_velocity(P.xc, P.yc, xv, yv, Gv, rc)[:2] if Gv > 0.0 else None
    sigma, vt = _solve_sources(P, U, alpha, gam, w, extra)
    du, dv = _panel_velocity(np.array([xv]), np.array([yv]), P, sigma, gam)
    Ve = float(np.hypot(U*np.cos(alpha) + du[0], U*np.sin(alpha) + dv[0]))
    return Reconstruction(P.xp, P.yp, P.xc, P.yc, P.L, sigma, gam, float(U), alpha, float(chord), float(CL),
                          float(xv), float(yv), float(Gv), float(rc), Ve, float(vt), w, float(x_sep), P,
                          float(factor), float(mach))


def solve(coords_csv: Path, chord: float, U: float, alpha_deg: float, CL: float,
          CN_vortex: float = 0.0, tau_over_Tvl: float = 0.0, f_sep: float | None = None,
          n_panel: int = N_PANELS, mach: float = 0.0, factor: float | None = None,
          air: Air = AIR) -> Reconstruction:
    """Panel solution for one instant.

    chord [m], U [m/s], alpha_deg [deg]; CL is the lift coefficient the
    reconstruction is given; CN_vortex is the vortex part of the normal force
    (vortex circulation 0.5 CN_vortex U chord) and tau_over_Tvl the vortex
    clock as a fraction of T_VL; there is no vortex unless both are
    positive. f_sep is the load model's separation point (fraction of chord):
    with it, the suction surface aft of it carries the transpiration that
    satisfies the Kutta condition (see `outflow_shape`); with None the wall is
    solid everywhere and the trailing edge is left as the given circulation
    makes it.

    The sheet carries `factor` times the circulation 0.5 CL U chord. With
    mach = 0 the factor is 1: the incompressible pressure returns CL but for
    the closure error. With a stream Mach number `mach` > 0 and no factor
    given, the factor is solved so that the surface pressure corrected for
    compressibility at that Mach number (Karman-Tsien, `surface_state`)
    returns CL: it starts from the Prandtl-Glauert value sqrt(1 - mach^2) and
    ends below it where the suction is strong, because the correction grows
    faster than that value there. A factor given is used as it is.

    The panels of a coordinate file and the inverse of their wall-condition
    matrix are built once and kept."""
    alpha = float(np.radians(alpha_deg))
    P = _panels_of(coords_csv, chord, n_panel)
    args = (P, chord, U, alpha, CL, CN_vortex, tau_over_Tvl, f_sep)
    if factor is not None or mach <= 0.0 or abs(CL) < 1e-9:
        return _instant(*args, 1.0 if factor is None else factor, mach)
    s0 = math.sqrt(1.0 - mach*mach)
    rec = _instant(*args, s0, mach)
    e0 = surface_load_closure(rec, mach, air)["dCL"]
    s1 = s0*CL/(CL + e0)
    for _ in range(CLOSURE_ITERATIONS):
        rec = _instant(*args, s1, mach)
        e1 = surface_load_closure(rec, mach, air)["dCL"]
        if abs(e1) <= CLOSURE_TOL or e1 == e0:
            break
        s0, s1, e0 = s1, s1 - e1*(s1 - s0)/(e1 - e0), e1
    return rec


def from_solution(coords_csv: Path, chord: float, U: float, out: dict, i: int, Tvl: float,
                  vortex: bool = True, outflow: bool = True, n_panel: int = N_PANELS,
                  mach: float = 0.0, factor: float | None = None, air: Air = AIR) -> Reconstruction:
    """Panel solution at step i of `out`, the dict unistall.dsmodel.solve
    returns. chord [m], U [m/s], Tvl in semichords (the calibrated T_VL).
    `vortex=False` leaves the vortex out and keeps the same lift;
    `outflow=False` leaves the transpiration out (solid wall everywhere);
    `mach`, `factor` and `air` as in `solve`."""
    cnv = float(out["CN_vortex"][i]) if vortex else 0.0
    return solve(coords_csv, chord, U, float(out["alpha_deg"][i]), float(out["CL"][i]),
                 cnv, float(out["tau_v"][i])/Tvl, float(out["f_sep"][i]) if outflow else None, n_panel, mach,
                 factor, air)


def _panels(rec: Reconstruction) -> _Panels:
    """The panels of a reconstruction (built from its end points if it was
    made without them)."""
    return rec.panels if rec.panels is not None else _build_panels(rec.xp, rec.yp)


def without_vortex(rec: Reconstruction) -> Reconstruction:
    """The same instant solved again with the vortex left out: the same lift,
    the same transpiration shape, its strength solved afresh."""
    P = _panels(rec)
    sigma, vt = _solve_sources(P, rec.U, rec.alpha, rec.gam, rec.w)
    return Reconstruction(rec.xp, rec.yp, rec.xc, rec.yc, rec.L, sigma, rec.gam, rec.U, rec.alpha, rec.chord,
                          rec.CL, rec.xv, rec.yv, 0.0, rec.rc, rec.Ve, float(vt), rec.w, rec.x_sep, P, rec.factor,
                          rec.mach)


def outflow_state(rec: Reconstruction) -> dict:
    """The transpiration as numbers, all dimensionless: vt_over_U (its speed
    at the trailing edge over the stream speed; positive outflow, negative
    inflow), flux_over_Uc (volume flux through the surface per unit span over
    U c), x_sep (chord fraction at which it starts) and side (+1 upper
    surface, -1 lower, 0 none)."""
    upper = np.arange(len(rec.xc)) < len(rec.xc)//2
    side = 0 if not np.any(rec.w) else (1 if np.any(rec.w[upper]) else -1)
    return dict(vt_over_U=rec.vt/rec.U, flux_over_Uc=float(rec.vt*np.sum(rec.w*rec.L)/(rec.U*rec.chord)),
                x_sep=rec.x_sep, side=side)


def circulation_state(rec: Reconstruction) -> dict:
    """The circulations over U c, all dimensionless: `bound` (the sheet, which
    is 0.5 CL times the factor of the reconstruction), `vortex` (0.5 C_N^v when the vortex is present) and `total`,
    their sum, which is what a contour round the section and the vortex
    encloses. The lift the reconstruction is given corresponds to `bound`
    alone."""
    Uc = rec.U*rec.chord
    bound = float(rec.gam*rec.L.sum()/Uc)
    return dict(bound=bound, vortex=rec.Gv/Uc, total=bound + rec.Gv/Uc)


# --------------------------------------------------------------------------- #
#  field quantities
# --------------------------------------------------------------------------- #
def _core_pressure_deficit(r, Gv, rc, U):
    """Cp to add to Bernoulli's value so that the pressure inside the vortex
    core satisfies radial equilibrium, dp/dr = rho v_theta^2 / r. It vanishes
    outside the core, where the two agree."""
    r = np.asarray(r, float)
    if Gv <= 0.0 or rc <= 0.0 or U <= 0.0:
        return np.zeros_like(r)
    rmax = max(float(np.nanmax(r)), 12.0*rc)
    r0 = rc*1e-4
    rr = r0*np.exp(np.linspace(0.0, np.log(rmax/r0), CORE_INTEGRAL_POINTS))   # geometric: resolves any rc
    vt = (Gv/(2.0*np.pi*rr))*(1.0 - np.exp(-rr**2/rc**2))
    g = vt*vt/rr
    inner = np.concatenate([[0.0], np.cumsum(0.5*(g[1:] + g[:-1])*np.diff(rr))])
    tail = (Gv/(2.0*np.pi))**2/(2.0*rmax**2)                  # the integral beyond rmax, in closed form
    cp_equilibrium = -(2.0/U**2)*((inner[-1] + tail) - inner)
    return np.interp(r, rr, cp_equilibrium + (vt/U)**2)


def sonic_values(M: float, air: Air = AIR) -> tuple:
    """(sonic speed [m/s], critical pressure coefficient C_p* [-]) for stream
    Mach number M: the speed at which the isentropic local Mach number is 1
    in a flow of uniform total temperature, and the isentropic pressure
    coefficient at that point."""
    g = air.gamma
    ratio = (2.0/(g + 1.0))*(1.0 + 0.5*(g - 1.0)*M*M)         # sonic over stream temperature
    return float(np.sqrt(g*air.R_gas*air.T_inf*ratio)), float((2.0/(g*M*M))*(ratio**(g/(g - 1.0)) - 1.0))


def limit_to_sonic(u: np.ndarray, v: np.ndarray, Cp: np.ndarray, M: float, air: Air = AIR) -> tuple:
    """A sonic limiter: (u, v, Cp, limited) with the local speed held at the
    sonic speed and Cp held at the critical value C_p* of stream Mach number M.

    It is a limiter, not a compressible solution, and it is not consistent:
    Bernoulli's incompressible Cp reaches C_p* at a speed below the
    isentropic sonic speed, so `limited` [-] is 1 where only Cp is held, 2
    where the speed is held as well and 0 elsewhere. u, v [m/s] keep their
    direction. With M <= 0 nothing is limited. The reconstruction applies it
    only where a caller asks with `limit=True`; the case study does not."""
    u, v, Cp = np.asarray(u, float), np.asarray(v, float), np.asarray(Cp, float)
    if M <= 0.0:
        return u, v, Cp, np.zeros_like(Cp)
    v_sonic, cp_star = sonic_values(M, air)
    speed = np.hypot(u, v)
    fast = speed > v_sonic
    scale = np.where(fast, v_sonic/np.where(fast, speed, 1.0), 1.0)
    held = Cp < cp_star
    return u*scale, v*scale, np.where(held, cp_star, Cp), np.where(fast, 2.0, np.where(held, 1.0, 0.0))


def thermal_state(speed: np.ndarray, M: float, air: Air = AIR) -> tuple:
    """(static temperature [K], adiabatic recovery temperature [K], local Mach
    number [-]) at flow speed `speed` [m/s] for stream Mach number M, from the
    energy equation with uniform total temperature. These are relations of a
    compressible flow; fed with the speeds of the incompressible
    reconstruction they do not give its Mach number or temperature, and the
    case study does not use them."""
    speed = np.asarray(speed, float)
    T0 = air.T_inf*(1.0 + 0.5*(air.gamma - 1.0)*M*M)
    T_static = T0 - speed**2/(2.0*air.cp)
    T_recovery = T0 - (1.0 - air.recovery)*speed**2/(2.0*air.cp)
    return T_static, T_recovery, speed/np.sqrt(air.gamma*air.R_gas*np.maximum(T_static, 1.0))


def vacuum_cp(M: float, air: Air = AIR) -> float:
    """The pressure coefficient of a vacuum at stream Mach number M,
    -2/(gamma M^2): no pressure can be lower."""
    return -2.0/(air.gamma*M*M)


def karman_tsien(Cp: np.ndarray, M: float, air: Air = AIR) -> tuple:
    """(corrected Cp, capped): the Karman-Tsien compressibility correction of
    an incompressible pressure coefficient at stream Mach number M,

        Cp_c = Cp / (beta + M^2 Cp / (2 (1 + beta))),   beta = sqrt(1 - M^2).

    It is a linearised correction and holds only where the local flow is
    subcritical. The rule has a pole at strong suction; the corrected value
    is not allowed below VACUUM_CAP of the vacuum value (`vacuum_cp`), and
    `capped` is True where that bound acts. With M <= 0 the pressure is
    returned as it is."""
    Cp = np.asarray(Cp, float)
    if M <= 0.0:
        return Cp.copy(), np.zeros(Cp.shape, bool)
    beta = math.sqrt(1.0 - M*M)
    floor = VACUUM_CAP*vacuum_cp(M, air)
    den = beta + M*M*Cp/(2.0*(1.0 + beta))
    capped = (den <= 0.0) | (Cp < floor*np.where(den > 0.0, den, 1.0))
    return np.where(capped, floor, Cp/np.where(den > 0.0, den, 1.0)), capped


def compressible_state(Cp: np.ndarray, M: float, air: Air = AIR) -> dict:
    """The local state that goes with a pressure coefficient Cp (corrected
    for compressibility) at stream Mach number M > 0, by the isentropic
    relations with the stagnation state of the stream:

        p/p_inf = 1 + gamma M^2 Cp / 2,
        M_l^2   = 2/(gamma - 1) ((p0/p)^((gamma - 1)/gamma) - 1),
        T       = T0 / (1 + (gamma - 1) M_l^2 / 2),
        T_r     = T (1 + r (gamma - 1) M_l^2 / 2),

    with p0, T0 the stagnation pressure and temperature of the stream and r
    the recovery factor of `air` (Pr^(1/3), turbulent boundary layer).
    Returns Mach [-], T_static and T_recovery [K] and beyond_critical (True
    where Cp is below the critical value, that is where Mach > 1). A
    pressure above the stagnation value gives Mach 0."""
    Cp = np.asarray(Cp, float)
    g = air.gamma
    stag = 1.0 + 0.5*(g - 1.0)*M*M
    p = np.maximum(1.0 + 0.5*g*M*M*Cp, 1e-12)                   # p / p_inf
    m2 = np.maximum((2.0/(g - 1.0))*((stag**(g/(g - 1.0))/p)**((g - 1.0)/g) - 1.0), 0.0)
    T = air.T_inf*stag/(1.0 + 0.5*(g - 1.0)*m2)
    return dict(Mach=np.sqrt(m2), T_static=T, T_recovery=T*(1.0 + air.recovery*0.5*(g - 1.0)*m2),
                beyond_critical=Cp < sonic_values(M, air)[1])


def smoothing_radius(rec: Reconstruction) -> float:
    """Distance from the wall [m] within which the field of the
    constant-strength panels shows their junctions:
    NEAR_WALL_PANEL_LENGTHS mean panel lengths. Nearer the wall than this
    the values on the surface itself (`surface_state`) are the ones to
    read."""
    return float(NEAR_WALL_PANEL_LENGTHS*rec.L.mean())


def _state(rec, u, v, r2, M, air, limit):
    """Speed, pressure and local state from the velocity (u, v) and the
    squared distance r2 to the vortex centre. Cp_incompressible is
    Bernoulli's value; Cp is that value corrected by the Karman-Tsien rule
    for stream Mach number M > 0 (the same value for M <= 0), and the local
    Mach number and temperatures follow from it. With `limit` the older
    sonic limiter is applied to the incompressible values instead and no
    correction is made."""
    cp_i = 1.0 - (u*u + v*v)/rec.U**2 + _core_pressure_deficit(np.sqrt(r2), rec.Gv, rec.rc, rec.U)
    limited = np.zeros_like(cp_i)
    if limit:
        u, v, cp_i, limited = limit_to_sonic(u, v, cp_i, M, air)
    out = dict(u=u, v=v, speed=np.hypot(u, v), Cp_incompressible=cp_i, limited=limited)
    if M > 0.0 and not limit:
        out["Cp"], out["capped"] = karman_tsien(cp_i, M, air)
        out.update(compressible_state(out["Cp"], M, air))
    else:
        nan = np.full(np.shape(cp_i), np.nan)
        out.update(Cp=cp_i.copy(), capped=np.zeros(np.shape(cp_i), bool), Mach=nan, T_static=nan.copy(),
                   T_recovery=nan.copy(), beyond_critical=np.zeros(np.shape(cp_i), bool))
    return out


def inside_section(rec: Reconstruction, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """True at the points (x, y) [m] that lie inside the section."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    return _Polygon(np.column_stack([rec.xp, rec.yp])).contains_points(
        np.column_stack([x.ravel(), y.ravel()])).reshape(x.shape)


def evaluate(rec: Reconstruction, x: np.ndarray, y: np.ndarray, M: float, air: Air = AIR,
             blank_inside: bool = True, limit: bool = False) -> dict:
    """The reconstructed field at arbitrary points (x, y) [m], arrays of one
    shape (grid nodes, mesh nodes, a probe line).

    M is the stream Mach number. Returns arrays of that shape: u, v, speed
    [m/s] and Cp_incompressible [-] (the incompressible solution), vorticity
    [1/s], and for M > 0 the pressure corrected for compressibility and the
    local state that goes with it: Cp [-] (Karman-Tsien), Mach [-], T_static
    and T_recovery [K], `beyond_critical` (True where the corrected Cp is
    below C_p*, the local Mach number above one: the correction does not
    hold there and the values are not physical) and `capped` (True where the
    corrected Cp is held at the vacuum-side bound). With M <= 0, Cp is the
    incompressible value and the Mach number and temperatures are NaN.
    `limited` is zero unless `limit` (see `_state`). Points inside the
    section are NaN when `blank_inside` (False in the flags). Within
    `smoothing_radius` of the wall the field shows the panel junctions;
    `surface_state` gives the values on the surface itself."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    pu, pv = _panel_velocity(x, y, _panels(rec), rec.sigma, rec.gam)
    vu, vv, r2 = _vortex_velocity(x, y, rec.xv, rec.yv, rec.Gv, rec.rc)
    out = _state(rec, rec.U*np.cos(rec.alpha) + pu + vu, rec.U*np.sin(rec.alpha) + pv + vv, r2, M, air, limit)
    out["vorticity"] = -(rec.Gv/(np.pi*rec.rc**2))*np.exp(-r2/rec.rc**2)
    if blank_inside:
        inside = inside_section(rec, x, y)
        for key, a in out.items():
            a[inside] = False if a.dtype == bool else (0.0 if key == "limited" else np.nan)
    return out


def _dilate(mask: np.ndarray, rings: int) -> np.ndarray:
    """Boolean grid mask grown by `rings` nodes in all eight directions. The
    grid does not wrap: a mask touching one edge does not reach the opposite
    one."""
    grown = mask.copy()
    ny, nx = mask.shape
    for _ in range(rings):
        padded = np.pad(grown, 1, constant_values=False)
        grown = np.zeros_like(mask)
        for sy in (0, 1, 2):
            for sx in (0, 1, 2):
                grown |= padded[sy:sy + ny, sx:sx + nx]
    return grown


def reconstruct_field(rec: Reconstruction, M: float, air: Air = AIR, domain: tuple = DOMAIN_CHORDS,
                      nx: int = GRID_NX, ny: int = GRID_NY, limit: bool = False,
                      rings: int = NEAR_WALL_RINGS_BLANKED) -> dict:
    """The field on a uniform grid of nx by ny nodes over `domain`
    (x0, x1, y0, y1) in chords. Returns the dict of `evaluate` with the node
    coordinates X, Y [m] and `inside` added: True at the nodes inside the
    section and at the `rings` rings of nodes next to it. The values at those
    nodes are kept: inside the section they are the field of the
    singularities, with no physical meaning, and are there so that a contour
    map can be drawn up to the outline without a stepped gap."""
    c = rec.chord
    gx, gy = np.linspace(domain[0]*c, domain[1]*c, nx), np.linspace(domain[2]*c, domain[3]*c, ny)
    X, Y = np.meshgrid(gx, gy)
    out = evaluate(rec, X, Y, M, air, blank_inside=False, limit=limit)
    return dict(X=X, Y=Y, inside=_dilate(inside_section(rec, X, Y), rings), **out)


def read_field(path: Path, chord: float, U: float) -> dict:
    """One stored field file as grids: its columns as arrays of shape
    (ny, nx), with x, y in chords, speed_over_U and vorticity_c_over_U
    added, `outside` True at the nodes outside the section and `valid` True
    at those of them where the corrected pressure is not below the critical
    value."""
    df = pd.read_csv(path)
    nx = df["x_m"].nunique()
    F = {c: df[c].to_numpy(float).reshape(-1, nx) for c in df.columns}
    F["x"], F["y"] = F["x_m"]/chord, F["y_m"]/chord
    F["speed_over_U"] = np.hypot(F["u_ms"], F["v_ms"])/U
    F["vorticity_c_over_U"] = F["vorticity_1s"]*chord/U
    F["outside"] = F["inside"] < 0.5
    F["valid"] = F["outside"] & (F["beyond_critical"] < 0.5)
    return F


# --------------------------------------------------------------------------- #
#  surface quantities
# --------------------------------------------------------------------------- #
def surface_state(rec: Reconstruction, M: float, air: Air = AIR, limit: bool = False) -> dict:
    """Values ON the surface, at the panel control points, in panel order
    (trailing edge, upper surface, leading edge, lower surface, trailing edge).

    Returns x_over_c [-], s_over_c (distance along the outline from the
    trailing edge) [-], speed [m/s], Cp_incompressible [-], and for stream
    Mach number M > 0 the corrected Cp [-], Mach [-], T_static and
    T_recovery [K], `beyond_critical` and `capped` as in `evaluate`;
    `limited` (zero unless `limit`), `upper` (True on the upper surface) and
    vn_over_U, the transpiration through the wall [-]."""
    P = _panels(rec)
    du, dv, rv2 = _vortex_velocity(rec.xc, rec.yc, rec.xv, rec.yv, rec.Gv, rec.rc)
    u, v = _surface_velocity(P, rec.U, rec.alpha, rec.gam, rec.sigma, (du, dv))
    out = _state(rec, u, v, rv2, M, air, limit)
    out.pop("u")
    out.pop("v")
    upper = np.arange(len(rec.xc)) < len(rec.xc)//2
    return dict(x_over_c=rec.xc/rec.chord, s_over_c=(np.cumsum(rec.L) - 0.5*rec.L)/rec.chord, upper=upper,
                vn_over_U=rec.vt*rec.w/rec.U, **out)


def surface_cp(rec: Reconstruction, M: float = 0.0, air: Air = AIR) -> tuple:
    """(x/c, Cp, upper) at the panel control points; see `surface_state`.
    With stream Mach number M > 0 the Cp is the one corrected for
    compressibility; with the default 0 it is the incompressible one."""
    st = surface_state(rec, M, air)
    return st["x_over_c"], st["Cp"], st["upper"]


def surface_load_closure(rec: Reconstruction, M: float = 0.0, air: Air = AIR) -> dict:
    """The loads returned by integrating the surface Cp, against the lift
    that pressure should return.

    With M = 0 the incompressible pressure is integrated and compared with
    the lift of the circulation the sheet carries (rec.factor rec.CL). With
    stream Mach number M > 0 the pressure corrected for compressibility is
    integrated and compared with the lift the reconstruction was given,
    rec.CL: this is the closure of the pressures the case study stores, and
    it is meaningful when the instant was solved for the same Mach number.

    Returns CL_given (the lift compared with), CL_from_Cp, dCL (their
    difference) and error_pct (100 dCL / CL_given, NaN at zero lift), all
    dimensionless, and te_jump, the difference in Cp between the two control
    points next to the trailing edge. With a solid wall, no vortex and M = 0
    the difference is discretisation error. Otherwise it also holds what the
    pressure integral does not see (the momentum of the transpiration, the
    force of the free vortex on the section); for M > 0 and an instant whose
    circulation was solved for that Mach number it is the tolerance of that
    solution, and with the circulation fixed at the Prandtl-Glauert factor it
    shows how far the Karman-Tsien correction departs from that factor.

    Also returns CM_from_Cp, the pitching moment coefficient of the surface
    pressure about the quarter chord (nose up positive), and CD_from_Cp, its
    pressure drag coefficient. The reconstruction is given the lift only:
    neither is constrained, and the caller compares them with the model."""
    _, cp, _ = surface_cp(rec, M, air)
    given = rec.CL if M > 0.0 else rec.factor*rec.CL
    xp, yp = rec.xp, rec.yp
    sense = 1.0 if 0.5*np.sum(xp[:-1]*yp[1:] - xp[1:]*yp[:-1]) > 0.0 else -1.0
    fn, fa = sense*cp*np.diff(xp), -sense*cp*np.diff(yp)       # normal and axial force of each panel / (q)
    CN, CA = np.sum(fn)/rec.chord, np.sum(fa)/rec.chord
    cl = float(CN*np.cos(rec.alpha) - CA*np.sin(rec.alpha))
    x_ref = xp.min() + 0.25*rec.chord
    cm = float(-np.sum((rec.xc - x_ref)*fn - rec.yc*fa)/rec.chord**2)
    d = cl - given
    return dict(CL_given=given, CL_from_Cp=cl, dCL=d,
                error_pct=100.0*d/given if given != 0.0 else float("nan"),
                te_jump=float(abs(cp[0] - cp[-1])),
                CM_from_Cp=cm, CD_from_Cp=float(CN*np.sin(rec.alpha) + CA*np.cos(rec.alpha)))


def critical_coverage(rec: Reconstruction, M: float, air: Air = AIR, field: dict | None = None) -> dict:
    """Where the corrected pressure is below the critical value of stream
    Mach number M, that is where the compressibility correction does not
    hold. surface_fraction and surface_fraction_capped: the fractions of the
    surface arc length beyond the critical pressure and held at the
    vacuum-side bound (the second is part of the first). With `field` (from
    `reconstruct_field`): field_area_over_c2, the area of the field grid
    outside the section beyond the critical pressure over chord^2, and
    field_fraction, the same as a fraction of that grid outside the section.
    All dimensionless."""
    st = surface_state(rec, M, air)
    out = dict(surface_fraction=float(np.sum(rec.L[st["beyond_critical"]])/rec.L.sum()),
               surface_fraction_capped=float(np.sum(rec.L[st["capped"]])/rec.L.sum()))
    if field is not None:
        cell = (field["X"][0, 1] - field["X"][0, 0])*(field["Y"][1, 0] - field["Y"][0, 0])/rec.chord**2
        outside = ~field["inside"]
        beyond = field["beyond_critical"] & outside
        out["field_area_over_c2"] = float(np.sum(beyond)*cell)
        out["field_fraction"] = float(np.sum(beyond)/np.sum(outside))
    return out


def limiter_coverage(rec: Reconstruction, M: float, air: Air = AIR, field: dict | None = None) -> dict:
    """Where the incompressible solution, without the compressibility
    correction, passes the critical state of stream Mach number M (the
    measure that goes with `limit_to_sonic`; the case study uses
    `critical_coverage`).

    surface_fraction_held and surface_fraction_sonic: the fractions of the
    surface arc length on which Cp is below C_p* and on which the speed
    exceeds the isentropic sonic speed (the second is part of the first).
    With `field` (from `reconstruct_field`), field_area_held_over_c2 and
    field_area_sonic_over_c2: the same as areas of the field grid outside
    the section over chord^2. All dimensionless."""
    lim = surface_state(rec, M, air, limit=True)["limited"]
    out = dict(surface_fraction_held=float(np.sum(rec.L[lim >= 1.0])/rec.L.sum()),
               surface_fraction_sonic=float(np.sum(rec.L[lim >= 2.0])/rec.L.sum()))
    if field is not None:
        cell = (field["X"][0, 1] - field["X"][0, 0])*(field["Y"][1, 0] - field["Y"][0, 0])/rec.chord**2
        v_sonic, cp_star = sonic_values(M, air)
        outside = ~field["inside"]
        out["field_area_held_over_c2"] = float(np.sum((field["Cp_incompressible"] < cp_star) & outside)*cell)
        out["field_area_sonic_over_c2"] = float(np.sum((field["speed"] > v_sonic) & outside)*cell)
    return out


def vortex_state(rec: Reconstruction) -> dict:
    """The reconstructed vortex as numbers: x_over_c, y_over_c of its centre,
    Gamma_over_Uc (circulation / (U c)), rc_over_c (core radius),
    peak_swirl_over_U (its largest swirl speed), peak_vorticity_c_over_U (the
    size of the vorticity at its centre, times c / U), edge_speed_over_U
    (speed of the panel flow at its centre) and Cp_centre, the pressure
    coefficient at its centre. All dimensionless."""
    c, U = rec.chord, rec.U
    here = evaluate(rec, np.array([rec.xv]), np.array([rec.yv]), 0.0, blank_inside=False)
    swirl = LAMB_OSEEN_PEAK*rec.Gv/(2.0*np.pi*rec.rc)/U if rec.Gv > 0.0 else 0.0
    return dict(x_over_c=rec.xv/c, y_over_c=rec.yv/c, Gamma_over_Uc=rec.Gv/(U*c), rc_over_c=rec.rc/c,
                peak_swirl_over_U=float(swirl), peak_vorticity_c_over_U=float(rec.Gv/(np.pi*rec.rc**2)*c/U),
                edge_speed_over_U=rec.Ve/U, Cp_centre=float(here["Cp"][0]))


def vortex_footprint(rec: Reconstruction) -> dict:
    """What the vortex does to the pressure on the surface beneath it.

    Returns over_chord (True if the vortex is present and its centre is over
    the chord), x_over_c of the upper-surface control point nearest below its
    centre, Cp_under (the surface Cp there), Cp_under_without_vortex (the
    same with the vortex left out, see `without_vortex`) and dCp_under, their
    difference: positive means the vortex raises the pressure beneath it,
    negative that it lowers it (a suction footprint). The three Cp values
    are NaN when over_chord is False."""
    x_over_c = rec.xv/rec.chord
    over = bool(rec.Gv > 0.0 and 0.0 <= x_over_c <= 1.0)
    nan = float("nan")
    if not over:
        return dict(over_chord=False, x_over_c=x_over_c, Cp_under=nan, Cp_under_without_vortex=nan, dCp_under=nan)
    st, bare = surface_state(rec, 0.0), surface_state(without_vortex(rec), 0.0)
    i = int(np.argmin(np.where(st["upper"], np.abs(rec.xc - rec.xv), np.inf)))
    return dict(over_chord=True, x_over_c=float(st["x_over_c"][i]), Cp_under=float(st["Cp"][i]),
                Cp_under_without_vortex=float(bare["Cp"][i]), dCp_under=float(st["Cp"][i] - bare["Cp"][i]))
