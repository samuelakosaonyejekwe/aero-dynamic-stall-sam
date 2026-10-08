# -*- coding: utf-8 -*-
"""
unistall / strokes.py
-----------------
Which points of a traced hysteresis loop are on the up-stroke and which on the
down-stroke.

Author: Akosa Samuel Onyejekwe (independent)
"""
import numpy as np

__all__ = ["stroke_split"]


def stroke_split(alpha_deg: np.ndarray) -> np.ndarray:
    """Label each point of a traced loop "up" or "down".

    alpha_deg: incidence of the points in the order they were traced [deg].
    The record is a closed loop that may start anywhere, so it is split on both
    turning points. A monotone sweep comes back as a single stroke.
    """
    a = np.asarray(alpha_deg, float)
    imax, imin = int(np.argmax(a)), int(np.argmin(a))
    s = np.empty(len(a), dtype="<U4")
    if imin <= imax:
        s[:] = "down"
        s[imin:imax+1] = "up"
    else:
        s[:] = "up"
        s[imax+1:imin+1] = "down"
    return s
