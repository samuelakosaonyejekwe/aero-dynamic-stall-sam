# -*- coding: utf-8 -*-
"""
unistall / static_model.py
----------------------
The static inputs of the load model as functions of Mach number: lift-curve
slope, zero-lift angle, the normal force at static stall, and the Kirchhoff
separation function f(alpha). Everything comes from measurements in the
CALIBRATION set or from the NASA report's static figures; nothing is fitted
here and no held-out frame is read.

Author: Akosa Samuel Onyejekwe (independent)

TWO MACH STATIONS, because static stall in this experiment depends strongly on
Mach number (maximum static lift 1.64 at M = 0.215, 1.36 at M = 0.302):

  M = 0.302   separation function from the static polar assembled in
              static_polar_naca0012_M030.csv (NASA TM-84245 Vol. 1 Fig. 16(a),
              Fig. 9, Table 7; Vol. 2 Fig. 4), 0 to 25 deg.
              load scale (slope, zero-lift angle) from quasi-steady sweep
              frame_12102, which is the same measurement chain as the dynamic
              frames (not wall-corrected).
  M = 0.215   separation function, slope and zero-lift angle from quasi-steady
              sweeps frame_13308 (-3 to 17 deg) and frame_12300 (10 to 30 deg,
              M = 0.204), both in the calibration set.

WHICH POST-STALL CURVE. Past static stall the data are not single-valued: two
branches at M = 0.30, an up-sweep and a down-sweep at M = 0.21. Each station
therefore carries two curves, named for the flow state and not for how they
were measured:
  "more"   attached up to static stall, then the more separated post-stall
           state: the lower lift branch at M = 0.30, the down-sweep at
           M = 0.21                                                <- DEFAULT
  "less"   the less separated state: the upper branch, the up-sweep
The default was chosen on the calibration loops, where it gives the lower
cost of the two (results/static_curve_choice.csv, written by
calibrate.py --curve).

THE SEPARATION FUNCTION is the Kirchhoff inverse,
    f = (2 sqrt(C_N / (C_Nalpha (alpha - alpha_0))) - 1)^2 ,   F_MIN <= f <= 1,
evaluated with the slope and zero-lift angle OF THE DATA SET THE C_N CAME FROM,
so that f is a property of the flow state and not of a tunnel correction.

BETWEEN STATIONS the slope, zero-lift angle and stall level are linear in Mach
number, and f is interpolated along incidence normalised by the static-stall
angle, so that the stall break moves with Mach number instead of being smeared
into two steps:
    f(alpha; M) = (1-w) f_a(alpha * as_a/as(M)) + w f_b(alpha * as_b/as(M)).
Outside 0.215 <= M <= 0.302 the nearer station is used unchanged and the
result is flagged as outside the static Mach range.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator

HERE = Path(__file__).resolve().parent
from unistall import metrics as mt
from unistall.strokes import stroke_split
from unistall.paths import DATA, RESULTS

F_MIN = 0.02                  # floor on the separation point
ATTACHED_BELOW_DEG = 2.0      # below this |alpha - alpha_0| the ratio is ill-conditioned: f = 1
ALPHA_MAX = 25.0
MERGE_DEG = 0.25              # points closer than this in incidence are one station
TABLE_STEP_DEG = 0.02         # resolution of the per-Mach lookup tables


def _monotone_xy(a, y):
    """Sort on incidence and average points that share a nominal incidence."""
    o = np.argsort(a, kind="stable")
    a, y = np.asarray(a, float)[o], np.asarray(y, float)[o]
    ga, gy, i = [], [], 0
    while i < len(a):
        j = i
        while j + 1 < len(a) and a[j+1] - a[i] <= MERGE_DEG:
            j += 1
        ga.append(a[i:j+1].mean())
        gy.append(y[i:j+1].mean())
        i = j + 1
    return np.array(ga), np.array(gy)


def _kirchhoff_f(alpha_deg, CN, slope_per_deg, alpha0_deg):
    x = np.asarray(alpha_deg, float) - alpha0_deg
    f = np.ones_like(x)
    m = np.abs(x) >= ATTACHED_BELOW_DEG
    ratio = np.clip(np.asarray(CN, float)[m]/(slope_per_deg*x[m]), 0.0, None)
    f[m] = np.clip(2.0*np.sqrt(ratio) - 1.0, np.sqrt(F_MIN), 1.0)**2
    return f


def _frame_sweep(name, stroke):
    """(alpha, C_N, C_M) of one stroke of a quasi-steady sweep frame, with C_N
    formed from the frame's own C_L and pressure C_D of that stroke."""
    fr = mt.load_frame(mt.FRAME_CACHE/f"{name}.mat")
    s = stroke_split(fr["acl"]) == stroke
    a, cl = _monotone_xy(fr["acl"][s], fr["cl"][s])
    sd = stroke_split(fr["acd"]) == stroke                  # drag of the same stroke as the lift
    ad, cd = _monotone_xy(fr["acd"][sd], fr["cd"][sd])
    am, cm = _monotone_xy(fr["acm"][stroke_split(fr["acm"]) == stroke],
                          fr["cm"][stroke_split(fr["acm"]) == stroke])
    cdi = np.interp(a, ad, cd)
    cn = cl*np.cos(np.radians(a)) + cdi*np.sin(np.radians(a))
    return a, cn, np.interp(a, am, cm), fr


CURVES = ("more", "less")


class Station:
    """Static description at one Mach number: slope [1/deg], zero-lift angle
    and static-stall angle [deg], normal force at static stall, and the
    separation point and static moment on each curve of CURVES, tabulated at
    the incidences `alpha` [deg]."""
    def __init__(self, M: float, slope_per_deg: float, alpha0_deg: float, alpha_stall_deg: float,
                 CN1: float, alpha: np.ndarray, f: dict, cm: dict, source: str) -> None:
        self.M, self.slope_per_deg, self.alpha0_deg = M, slope_per_deg, alpha0_deg
        self.alpha_stall_deg, self.CN1, self.source = alpha_stall_deg, CN1, source
        self.alpha = alpha
        self._f = {c: PchipInterpolator(alpha, f[c], extrapolate=False) for c in CURVES}
        self._cm = {c: PchipInterpolator(alpha, cm[c], extrapolate=False) for c in CURVES}

    def f(self, alpha_deg: np.ndarray, curve: str = "more") -> np.ndarray:
        """Separation point at incidence alpha_deg [deg] on the named curve."""
        q = np.clip(np.abs(np.asarray(alpha_deg, float)), self.alpha[0], self.alpha[-1])
        return np.clip(self._f[curve](q), F_MIN, 1.0)

    def cm_static(self, alpha_deg: np.ndarray, curve: str = "more") -> np.ndarray:
        """Static C_M at |alpha| [deg] on the named curve."""
        q = np.clip(np.abs(np.asarray(alpha_deg, float)), self.alpha[0], self.alpha[-1])
        return self._cm[curve](q)


def _station_030():
    pol = pd.read_csv(RESULTS/"static_polar_naca0012_M030.csv")
    cl = pd.read_csv(DATA/"static_naca0012_M030_CL.csv")
    lin = cl[(cl.alpha_deg > -6) & (cl.alpha_deg < 9)]
    s_fig, ic = np.polyfit(lin.alpha_deg, lin.Cl, 1)       # the static figure's own slope
    a0_fig = -ic/s_fig
    by = pd.read_csv(RESULTS/"static_by_mach_naca0012.csv")
    r = by[by.M > 0.29].iloc[0]
    a = pol.alpha_deg.values
    f = {"more": _kirchhoff_f(a, pol.Cn_lower.values, s_fig, a0_fig),
         "less": _kirchhoff_f(a, pol.Cn_upper.values, s_fig, a0_fig)}
    cm = {"more": pol.Cm_lower.values, "less": pol.Cm_upper.values}
    return Station(float(r.M), float(r.CL_alpha_per_deg), float(r.alpha_zero_lift_deg),
                   float(r.alpha_at_CL_max_deg), float(r.CN_at_CL_max), a, f, cm,
                   "f: NASA TM-84245 Vol. 1 Fig. 16(a) static polar (slope "
                   f"{s_fig:.4f}/deg, alpha_0 {a0_fig:+.2f} deg); scale: frame_12102")


def _station_021():
    by = pd.read_csv(RESULTS/"static_by_mach_naca0012.csv")
    r = by[by.M < 0.25].iloc[0]
    slope, a0 = float(r.CL_alpha_per_deg), float(r.alpha_zero_lift_deg)
    a1, cn1, cm1, _ = _frame_sweep("frame_13308", "up")       # -3 to 17 deg
    a2, cn2, cm2, _ = _frame_sweep("frame_12300", "up")       # 10 to 30 deg, M 0.204
    hi = a2 > a1.max() + MERGE_DEG
    au, cnu, cmu = np.r_[a1, a2[hi]], np.r_[cn1, cn2[hi]], np.r_[cm1, cm2[hi]]
    ad1, cnd1, cmd1, _ = _frame_sweep("frame_13308", "down")
    ad2, cnd2, cmd2, _ = _frame_sweep("frame_12300", "down")
    add_, cnd = _monotone_xy(np.r_[ad1, ad2], np.r_[cnd1, cnd2])
    _, cmd = _monotone_xy(np.r_[ad1, ad2], np.r_[cmd1, cmd2])
    grid = np.arange(0.0, ALPHA_MAX + 1e-9, 0.5)
    cn_less = np.interp(grid, au, cnu)
    cm_less = np.interp(grid, au, cmu)
    # the more separated state: attached up to the static-stall incidence of
    # the up-sweep (stall with increasing incidence happens there), and beyond
    # it the down-sweep, which is the lower of the two post-stall levels
    down = grid > float(r.alpha_at_CL_max_deg)
    cn_more = np.minimum(np.where(down, np.interp(grid, add_, cnd), cn_less), cn_less)
    cm_more = np.where(down, np.interp(grid, add_, cmd), cm_less)
    f = {"more": _kirchhoff_f(grid, cn_more, slope, a0), "less": _kirchhoff_f(grid, cn_less, slope, a0)}
    return Station(float(r.M), slope, a0, float(r.alpha_at_CL_max_deg), float(r.CN_at_CL_max),
                   grid, f, {"more": cm_more, "less": cm_less},
                   "f and scale: quasi-steady sweeps frame_13308 and frame_12300 (calibration set)")


_STATIONS = None


def stations() -> list:
    """The Mach stations, in order of Mach number (built once)."""
    global _STATIONS
    if _STATIONS is None:
        _STATIONS = sorted([_station_021(), _station_030()], key=lambda s: s.M)
    return _STATIONS


class StaticModel:
    """Static inputs at one Mach number, interpolated between the stations."""
    def __init__(self, M: float, curve: str = "more") -> None:
        a, b = stations()
        self.M, self.curve = float(M), curve
        self.in_range = bool(a.M - 0.016 <= M <= b.M + 0.016)
        w = float(np.clip((M - a.M)/(b.M - a.M), 0.0, 1.0))
        self._a, self._b, self.w = a, b, w

        def lerp(x: float, y: float) -> float:
            """Linear interpolation in Mach number between the two stations."""
            return (1.0 - w)*x + w*y
        self.slope_per_deg = lerp(a.slope_per_deg, b.slope_per_deg)
        self.CN_alpha = self.slope_per_deg*180.0/np.pi          # per radian
        self.alpha0_deg = lerp(a.alpha0_deg, b.alpha0_deg)
        self.alpha_stall_deg = lerp(a.alpha_stall_deg, b.alpha_stall_deg)
        self.CN1 = lerp(a.CN1, b.CN1)
        # The march looks f and the static moment up once per time step. Both
        # are tabulated here on a fine grid so each lookup is one linear
        # interpolation; TABLE_STEP_DEG is fine enough that the tabulation
        # changes f by less than 1e-4 anywhere.
        self._grid = np.arange(0.0, ALPHA_MAX + TABLE_STEP_DEG/2, TABLE_STEP_DEG)
        self._ftab = {c: self._between("f", self._grid, c) for c in CURVES}
        self._cmtab = {c: self._between("cm_static", self._grid, c) for c in CURVES}

    def _between(self, what: str, q: np.ndarray, c: str) -> np.ndarray:
        """Station values of f or the static moment, interpolated along
        incidence scaled by the static-stall angle."""
        va = getattr(self._a, what)(q*self._a.alpha_stall_deg/self.alpha_stall_deg, c)
        vb = getattr(self._b, what)(q*self._b.alpha_stall_deg/self.alpha_stall_deg, c)
        return (1.0 - self.w)*va + self.w*vb

    def f(self, alpha_deg: np.ndarray, curve: str | None = None) -> np.ndarray:
        """Separation point at incidence alpha_deg [deg]. The section is
        symmetric and its measured zero-lift incidence is an offset of the
        incidence scale, so f is even about that incidence."""
        return np.interp(self._mirror(alpha_deg), self._grid, self._ftab[curve or self.curve])

    def cm_static(self, alpha_deg: np.ndarray, curve: str | None = None) -> np.ndarray:
        """Static C_M (c/4) at alpha_deg [deg]. The table is for incidence
        above the zero-lift incidence; below it the deviation from the value
        at zero lift is mirrored (symmetric section)."""
        tab = self._cmtab[curve or self.curve]
        zero = np.interp(max(self.alpha0_deg, 0.0), self._grid, tab)
        return zero + np.sign(np.asarray(alpha_deg, float) - self.alpha0_deg)*(
            np.interp(self._mirror(alpha_deg), self._grid, tab) - zero)

    def _lookup(self, tab: np.ndarray, x: float) -> float:
        """Linear interpolation of a table on the uniform incidence grid at
        x [deg], held at its end values outside the grid."""
        u = x/TABLE_STEP_DEG
        if u <= 0.0:
            return float(tab[0])
        i = int(u)
        if i >= len(tab) - 1:
            return float(tab[-1])
        return float(tab[i] + (u - i)*(tab[i + 1] - tab[i]))

    def f_at(self, alpha_deg: float) -> float:
        """f at one incidence [deg] on the default curve; the same value as
        f(), for use inside the time march."""
        return self._lookup(self._ftab[self.curve], self.alpha0_deg + abs(alpha_deg - self.alpha0_deg))

    def cm_at(self, alpha_deg: float) -> float:
        """Static C_M (c/4) at one incidence [deg] on the default curve; the
        same value as cm_static(), for use inside the time march."""
        tab = self._cmtab[self.curve]
        zero = self._lookup(tab, max(self.alpha0_deg, 0.0))
        dev = self._lookup(tab, self.alpha0_deg + abs(alpha_deg - self.alpha0_deg)) - zero
        return zero + (dev if alpha_deg >= self.alpha0_deg else -dev)

    def _mirror(self, alpha_deg: np.ndarray) -> np.ndarray:
        """Incidence reflected about the zero-lift incidence onto the side the
        tables cover."""
        return self.alpha0_deg + np.abs(np.asarray(alpha_deg, float) - self.alpha0_deg)

    def CN_static(self, alpha_deg: np.ndarray, curve: str | None = None) -> np.ndarray:
        """Kirchhoff static normal force rebuilt from f, at alpha_deg [deg]."""
        x = np.radians(np.asarray(alpha_deg, float) - self.alpha0_deg)
        return self.CN_alpha*((1.0 + np.sqrt(self.f(alpha_deg, curve)))/2.0)**2*x


if __name__ == "__main__":
    rows = []
    for st in stations():
        sm = StaticModel(st.M)
        g = np.arange(4.0, ALPHA_MAX + 1e-9, 0.5)
        print(f"[static-model] M {st.M}: slope {st.slope_per_deg}/deg, alpha_0 {st.alpha0_deg:+.2f} deg, "
              f"static stall {st.alpha_stall_deg} deg, C_N1 {st.CN1}; "
              f"f at 10/14/16/18/20/25 deg: "
              + " ".join(f"{float(sm.f(a)):.2f}" for a in (10, 14, 16, 18, 20, 25))
              + f" | model static C_N max {float(sm.CN_static(g).max()):.3f}")
        for a in g:
            rows.append(dict(M=st.M, alpha_deg=a, f_more=round(float(sm.f(a, "more")), 4),
                             f_less=round(float(sm.f(a, "less")), 4),
                             CN_static_model=round(float(sm.CN_static(a)), 4),
                             Cm_static=round(float(sm.cm_static(a)), 4)))
    for M in (0.25, 0.28):
        sm = StaticModel(M)
        g = np.arange(4.0, ALPHA_MAX + 1e-9, 0.25)
        cn = sm.CN_static(g)
        print(f"[static-model] M {M} (interpolated, w = {sm.w:.2f}): slope {sm.slope_per_deg:.4f}/deg, "
              f"static stall {sm.alpha_stall_deg:.1f} deg, C_N1 {sm.CN1:.3f}; "
              f"model static C_N max {cn.max():.3f} at {g[cn.argmax()]:.2f} deg")
    pd.DataFrame(rows).to_csv(RESULTS/"static_model_table.csv", index=False)
