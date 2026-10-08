# -*- coding: utf-8 -*-
"""
unistall / naca4.py
---------------
The NACA four-digit thickness distribution, in one place, for the geometry,
the mesh, the flow-field reconstruction and the drawings.

Author: Akosa Samuel Onyejekwe (independent)
"""
import numpy as np

#  y_t / c = 5 t (a0 sqrt(x) + a1 x + a2 x^2 + a3 x^3 + a4 x^4)
#  with the OPEN trailing-edge coefficient a4 = -0.1015 (Abbott & von Doenhoff);
#  the closed-edge variant uses -0.1036, which brings the thickness to zero at
#  x = 1. The case study uses the CLOSED form everywhere (`outline`).
A0, A1, A2, A3, A4_OPEN = 0.2969, -0.1260, -0.3516, 0.2843, -0.1015
A4_CLOSED = -0.1036
NOSE_RADIUS_FACTOR = 1.1019          # r/c = 1.1019 t^2 for the four-digit sections


def thickness(x: np.ndarray, t: float) -> np.ndarray:
    """Half-thickness y_t/c at chordwise positions x (fractions of chord) for
    thickness ratio t (0.12 for a NACA 0012)."""
    x = np.asarray(x, float)
    return 5*t*(A0*np.sqrt(x) + A1*x + A2*x**2 + A3*x**3 + A4_OPEN*x**4)


def thickness_ratio(code: str) -> float:
    """Thickness ratio t of a symmetric four-digit section from its
    designation: 0.12 for "0012"."""
    if len(code) != 4 or not code.isdigit() or code[:2] != "00":
        raise ValueError(f"not a symmetric NACA four-digit designation: {code!r}")
    return int(code[2:])/100.0


def closed_thickness(x: np.ndarray, t: float) -> np.ndarray:
    """Half-thickness y_t/c of the closed-edge form (zero at x = 1) at
    chordwise positions x (fractions of chord) for thickness ratio t."""
    x = np.asarray(x, float)
    return 5*t*(A0*np.sqrt(x) + A1*x + A2*x**2 + A3*x**3 + A4_CLOSED*x**4)


def closed_thickness_slope(x: float, t: float) -> float:
    """d(y_t/c)/d(x/c) of the closed-edge form at x (fraction of chord, > 0)."""
    return 5.0*t*(A0/(2.0*np.sqrt(x)) + A1 + 2.0*A2*x + 3.0*A3*x**2 + 4.0*A4_CLOSED*x**3)


def outline(code: str, n_per_side: int) -> tuple:
    """(X, Y, x, yt): the closed outline of the symmetric section `code`
    (fractions of chord), running trailing edge -> upper surface -> leading
    edge -> lower surface -> trailing edge with n_per_side cosine-spaced
    stations on each surface, and those stations x with their half-thickness
    yt. The trailing edge is the single point (1, 0)."""
    x = (1.0 - np.cos(np.linspace(0.0, np.pi, n_per_side)))/2.0
    yt = closed_thickness(x, thickness_ratio(code))
    yt[0] = yt[-1] = 0.0
    return np.concatenate([x[::-1], x[1:]]), np.concatenate([yt[::-1], -yt[1:]]), x, yt
