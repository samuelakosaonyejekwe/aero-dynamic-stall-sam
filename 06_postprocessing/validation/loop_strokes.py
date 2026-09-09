# -*- coding: utf-8 -*-
"""
06_postprocessing / validation / loop_strokes.py
------------------------------------------------
One answer to one question: given the incidences of a traced hysteresis loop,
which points are on the up-stroke and which on the down-stroke.

There were two answers before this module. validate_nasa_real.py split the
record at its two turning points, which is right for a closed loop that may
start anywhere in the cycle. validate_digitized.py had no answer at all: its
`stroke` column is documented as optional and when it was omitted every point
was labelled "up", so down-stroke points were compared against the up-stroke
model a whole loop-width away (RMS_CL 0.3158 on data that was the solver's own
output scaled by 1.02) and the matched model loop collapsed to exactly zero
area, published as a -100.0 % error. Replacing that with a third method --
the sign of the increment between successive points -- would have left two
implementations that disagree on jittery digitised data, where a single
non-monotone point flips a label. So there is one, and it is the robust one.

Importable with no side effects, unlike the two harnesses that use it.
"""
import numpy as np

__all__ = ["stroke_split"]


def stroke_split(a):
    """Label each point of a traced loop "up" or "down".

    The record is a closed loop that may start anywhere, so split on BOTH
    turning points rather than assuming it opens on the upstroke. A monotone
    sweep has its extremes at the two ends and comes back as a single stroke.
    """
    a = np.asarray(a, float)
    imax, imin = int(np.argmax(a)), int(np.argmin(a))
    s = np.empty(len(a), dtype="<U4")
    if imin <= imax:                       # ... min ... max ...  -> rising between them
        s[:] = "down"; s[imin:imax+1] = "up"
    else:                                  # ... max ... min ...  -> falling between them
        s[:] = "up";   s[imax+1:imin+1] = "down"
    return s
