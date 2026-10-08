# -*- coding: utf-8 -*-
"""
unistall / reference_lb.py
----------------------
The reference model: the SAME Leishman-Beddoes load model as unistall/dsmodel.py
with one thing changed -- the static separation point is the standard
three-parameter exponential fit instead of a table of the measured values.

Author: Akosa Samuel Onyejekwe (independent)

    f(alpha) = 1 - 0.3 exp((|a| - alpha_1)/S_1)            |a| <= alpha_1
             = 0.04 + 0.66 exp((alpha_1 - |a|)/S_2)        |a| >  alpha_1
    with a = alpha - alpha_0             (Leishman; [DH19] Eq. 1.33)

alpha_1, S_1 and S_2 are fitted by least squares to the same static
separation-point values the tabulated model uses, at each of the two Mach
stations. Between the stations the two fitted laws are blended exactly as the
tabulated model blends its two tables: each is evaluated at the incidence
scaled by its static-stall angle, and the values are weighted linearly in Mach
number. The law is then held on the same incidence grid as the table. Slope, zero-lift angle, onset level C_N1, the static moment table
and every dynamic equation are identical in the two models, so a difference
between their scores is a difference between the two separation laws and
nothing else.

Usage:
  python3 -m unistall.reference_lb     fit and print the exponential constants and
                                   their residual against the table
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
from unistall import dsmodel as dm
from unistall.static_model import StaticModel, stations, F_MIN, ALPHA_MAX
from unistall.paths import RESULTS

FIT_FROM_DEG, FIT_STEP_DEG = 2.0, 0.25


def exponential_f(x_deg: np.ndarray, alpha1: float, S1: float, S2: float) -> np.ndarray:
    """[DH19] Eq. (1.33); x_deg is incidence measured from zero lift [deg];
    alpha1, S1, S2 in degrees."""
    x = np.abs(np.asarray(x_deg, float))
    return np.where(x <= alpha1,
                    1.0 - 0.3*np.exp((x - alpha1)/S1),
                    0.04 + 0.66*np.exp((alpha1 - x)/S2))


def fit_station(st: object) -> dict:
    """alpha_1, S_1, S_2 for one Mach station, fitted to its tabulated f."""
    a = np.arange(FIT_FROM_DEG, ALPHA_MAX + 1e-9, FIT_STEP_DEG)
    f_tab = st.f(a)
    x = a - st.alpha0_deg
    def res(p: np.ndarray) -> np.ndarray:
        """Fitted minus tabulated separation point for constants p [deg]."""
        return exponential_f(x, *p) - f_tab
    sol = least_squares(res, [st.alpha_stall_deg, 3.0, 2.0],
                        bounds=(FIT_LOWER, FIT_UPPER))
    rms = float(np.sqrt(np.mean(sol.fun**2)))
    return dict(M=st.M, alpha1_deg=float(sol.x[0]), S1_deg=float(sol.x[1]), S2_deg=float(sol.x[2]),
                rms_f=rms, max_abs_f=float(np.max(np.abs(sol.fun))))


_FITS = None
FIT_LOWER, FIT_UPPER = [5.0, 0.02, 0.02], [25.0, 20.0, 20.0]     # alpha_1, S_1, S_2 [deg]


def station_fits() -> list:
    """The exponential constants at each Mach station (fitted once)."""
    global _FITS
    if _FITS is None:
        _FITS = [fit_station(st) for st in stations()]
    return _FITS


class ExponentialStatic(StaticModel):
    """StaticModel with the fitted exponential separation law."""
    def __init__(self, M: float, curve: str = "more") -> None:
        super().__init__(M, curve)
        fa, fb = station_fits()
        a, b = stations()
        self._exp = []
        for st, ft in ((a, fa), (b, fb)):
            self._exp.append((st.alpha_stall_deg, st.alpha0_deg, ft["alpha1_deg"], ft["S1_deg"], ft["S2_deg"]))
        # the law replaces the table on the same incidence grid, so the march
        # and every other look-up read it exactly as they read the table
        self._ftab = {curve: self._law(self._grid)}

    def _law(self, alpha_deg: np.ndarray) -> np.ndarray:
        """The fitted law at alpha_deg [deg]: the two stations' laws, each at
        the incidence scaled by its static-stall angle, blended in Mach."""
        q = self._mirror(alpha_deg)
        out = 0.0
        for wgt, (a_ss, a0, a1, S1, S2) in zip((1.0 - self.w, self.w), self._exp, strict=False):
            out = out + wgt*exponential_f(q*a_ss/self.alpha_stall_deg - a0, a1, S1, S2)
        return np.clip(out, F_MIN, 1.0)


def frame_runner(consts: dict | None = None, **kw: object) -> object:
    """run(fr) for unistall/metrics.py, with the fitted separation law."""
    def run(fr: dict) -> dict:
        """Solve one loaded frame (angles in degrees as stored)."""
        return dm.solve(fr["a0"], fr["da"], fr["k"], fr["M"], consts=consts,
                        static=ExponentialStatic(fr["M"]), **kw)
    return run


if __name__ == "__main__":
    df = pd.DataFrame(station_fits()).round(4)
    df.to_csv(RESULTS/"reference_exponential_fit.csv", index=False)
    print(df.to_string(index=False))
    for st in stations():
        e = ExponentialStatic(st.M)
        t = StaticModel(st.M)
        g = np.arange(4.0, ALPHA_MAX + 1e-9, 0.25)
        print(f"[reference] M {st.M}: static C_N maximum {float(e.CN_static(g).max()):.3f} (exponential) "
              f"against {float(t.CN_static(g).max()):.3f} (table); f at 12/14/16/18/22 deg: "
              + " ".join(f"{float(e.f(a)):.2f}/{float(t.f(a)):.2f}" for a in (12, 14, 16, 18, 22)))
