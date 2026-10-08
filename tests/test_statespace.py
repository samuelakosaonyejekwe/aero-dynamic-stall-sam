"""The state-space form must be the same model as the indicial march.
Author: Akosa Samuel Onyejekwe (independent)"""
import json
import math

import numpy as np
import pytest

from unistall import dsmodel as dm
from unistall import metrics as mt
from unistall import statespace as ss
from unistall.paths import DATA, RESULTS

pytestmark = pytest.mark.skipif(not (mt.FRAME_CACHE/"frame_12102.mat").exists(),
                                reason="measured frames not fetched: run python3 -m unistall.fetch_frames")


def test_a_settled_state_does_not_move():
    """At constant incidence every lag has caught up, so dx/ds is zero."""
    model = ss.StateSpace(0.30, 102.0, 0.61)
    x = model.initial_state(math.radians(6.0), 0.0, 0.0)
    mode = dict(Tf=model.p["Tf0"], Tv=model.p["Tv0"], feed=False)
    assert len(x) == ss.N_STATES
    assert np.max(np.abs(model.rhs(x, math.radians(6.0), 0.0, 0.0, mode))) < 1e-12


def test_state_space_reproduces_the_indicial_loads():
    """Attached, light-stall and deep-stall cycles: the two forms give the same loads."""
    for a0, a1, k, M, tol in ((5.0, 5.0, 0.10, 0.30, 1e-3), (10.0, 5.0, 0.10, 0.30, 5e-3), (10.0, 10.0, 0.10, 0.30, 5e-3)):
        a, b = dm.solve(a0, a1, k, M), ss.solve(a0, a1, k, M)
        for q in ("CL", "CM", "CD"):
            assert float(np.max(np.abs(a[q] - b[q]))) < tol, (a0, a1, q)


def test_state_space_covers_low_mach_number():
    """Below M = 0.20 the incompressible loads are blended in; the two forms still agree."""
    for M in (0.05, 0.12, 0.15):
        a, b = dm.solve(8.0, 6.0, 0.10, M), ss.solve(8.0, 6.0, 0.10, M)
        for q in ("CL", "CM", "CD"):
            assert float(np.max(np.abs(a[q] - b[q]))) < 5e-3, (M, q)


def test_recorded_agreement_meets_its_target():
    import pandas as pd
    target = json.load(open(DATA/"targets.json"))["numerical"]["state_space_loop_error_difference_max"]
    assert float(pd.read_csv(RESULTS/"statespace_summary.csv").worst_difference_in_loop_error.iloc[0]) <= target
