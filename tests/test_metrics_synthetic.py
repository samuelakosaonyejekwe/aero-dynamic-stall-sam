"""The scoring functions on loops whose answers are known in closed form.
Author: Akosa Samuel Onyejekwe (independent)"""
import numpy as np

from unistall import metrics as mt
from unistall.strokes import stroke_split


def _ellipse(amp_deg, cm_amp, lag_rad, n=721, start=0.0):
    """alpha = a1 sin(ph), C_M = -cm_amp sin(ph - lag): the closed integral of
    C_M d(alpha) is pi a1 cm_amp sin(lag), so Xi = -cm_amp sin(lag) / a1 [rad]."""
    ph = start + np.linspace(0.0, 2.0*np.pi, n)
    return 10.0 + amp_deg*np.sin(ph), -cm_amp*np.sin(ph - lag_rad)


def test_cycle_damping_of_an_ellipse():
    for lag, cm_amp, amp in ((0.3, 0.05, 10.0), (-0.2, 0.08, 5.0), (0.0, 0.05, 8.0)):
        a, cm = _ellipse(amp, cm_amp, lag)
        exact = -cm_amp*np.sin(lag)/np.radians(amp)
        assert abs(mt.cycle_damping(a, cm, amp) - exact) < 2e-5
    a, cm = _ellipse(10.0, 0.05, 0.3)
    assert abs(mt.cycle_damping(a[::-1], cm[::-1], 10.0) + mt.cycle_damping(a, cm, 10.0)) < 1e-12     # traversed backwards


def test_cycle_damping_does_not_depend_on_where_the_record_starts():
    a0, c0 = _ellipse(10.0, 0.05, 0.3)
    a1, c1 = _ellipse(10.0, 0.05, 0.3, start=1.234)
    assert abs(mt.cycle_damping(a0, c0, 10.0) - mt.cycle_damping(a1, c1, 10.0)) < 1e-6


def test_stroke_split_on_a_loop_started_anywhere():
    for start in (0.0, 1.0, 2.5, 4.0, 5.5):
        ph = start + np.linspace(0.0, 2.0*np.pi, 200, endpoint=False)
        a = 10.0 + 8.0*np.sin(ph)
        s = stroke_split(a)
        rising = np.cos(ph) > 0.0
        interior = np.abs(np.cos(ph)) > 0.05                # away from the two turning points
        assert np.all((s[interior] == "up") == rising[interior])


def test_stroke_split_of_a_monotone_sweep_is_one_stroke():
    assert set(stroke_split(np.linspace(-5.0, 25.0, 40))) == {"up"}


def test_moment_stall_incidence_of_a_ramp_with_a_break():
    a = np.linspace(0.0, 20.0, 201)
    cm = np.where(a < 14.0, 0.01, 0.01 - 0.05*(a - 14.0))      # falls 0.05 per degree from 14 deg
    found = mt.moment_stall_alpha(a, cm)
    assert abs(found - (14.0 + mt.MSTALL_DROP/0.05)) < 1e-9
    assert np.isnan(mt.moment_stall_alpha(a, np.full_like(a, 0.01)))


def test_loop_error_of_a_known_offset():
    ph = np.linspace(0.0, 2.0*np.pi, 721)
    alpha, rate = 10.0 + 10.0*np.sin(ph), np.cos(ph)
    model = 1.0 + 0.1*np.sin(ph) + 0.05*np.cos(ph)              # different on the two strokes
    pick = np.arange(5, 716, 30)
    assert mt.loop_rms(alpha[pick], model[pick], alpha, rate, model) < 2e-3
    assert abs(mt.loop_rms(alpha[pick], model[pick] + 0.02, alpha, rate, model) - 0.02) < 2e-3


def test_loop_is_closed_needs_both_strokes():
    ph = np.linspace(0.0, 2.0*np.pi, 60, endpoint=False)
    assert mt.loop_is_closed(10.0 + 5.0*np.sin(ph))
    assert not mt.loop_is_closed(np.linspace(5.0, 15.0, 30))
