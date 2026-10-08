"""Static inputs and the dynamic model. These read the calibration-set frames
from the fetch cache (run python3 -m unistall.fetch_frames first).
Author: Akosa Samuel Onyejekwe (independent)"""
import json
from pathlib import Path
import numpy as np
import pytest
from unistall import metrics as mt

PKG = Path(__file__).resolve().parents[1]/"unistall"
from unistall.paths import DATA
T = json.load(open(DATA/"targets.json"))
pytestmark = pytest.mark.skipif(not (mt.FRAME_CACHE/"frame_12102.mat").exists(),
                                reason="measured frames not fetched: run python3 -m unistall.fetch_frames")


def test_static_model_rebuilds_its_stations():
    from unistall.static_model import StaticModel, stations
    g = np.arange(4.0, 25.01, 0.25)
    for st in stations():
        cn = StaticModel(st.M).CN_static(g)
        assert abs(cn.max()/st.CN1 - 1) < 0.02
        f = StaticModel(st.M).f(g)
        assert f.min() >= 0.02 - 1e-12 and f.max() <= 1.0 + 1e-12


def test_no_negative_drag_beyond_the_limit():
    import unistall.dsmodel as dm
    for a0, a1, k in ((9.0, 5.0, 0.20), (10.0, 10.0, 0.10), (5.0, 5.0, 0.10)):
        o = dm.solve(a0, a1, k, 0.30)
        assert float(np.mean(o["CD"][:-1] < 0.0)) <= T["separated_flow"]["CD_negative_fraction_of_cycle_max"]
        assert o["CD"].min() >= T["separated_flow"]["CD_min"]


def test_static_stall_moves_smoothly_with_mach():
    from unistall.static_model import StaticModel
    g = np.arange(4.0, 25.01, 0.05)
    stall = [g[StaticModel(M).CN_static(g).argmax()] for M in (0.215, 0.24, 0.27, 0.302)]
    assert all(b <= a for a, b in zip(stall, stall[1:], strict=False))


def test_exponential_fit_follows_the_table():
    from unistall import reference_lb as ref
    for ft in ref.station_fits():
        assert ft["rms_f"] < 0.03


def test_onset_does_not_depend_on_the_time_step():
    from unistall import dsmodel as dm
    a = dm.solve(10.0, 10.0, 0.10, 0.30)
    b = dm.solve(10.0, 10.0, 0.10, 0.30, n_per_cycle=4*a["meta"]["steps_per_cycle"])
    assert abs(a["onset_alpha_deg"] - b["onset_alpha_deg"]) < T["separated_flow"]["onset_shift_with_4x_steps_max_deg"]


def test_attached_model_is_recovered_without_separation():
    from unistall import dsmodel as dm
    o = dm.solve(3.0, 2.0, 0.05, 0.30, separation=False, vortex=False)
    assert np.all(o["f_sep"] == 1.0) and np.all(o["CN_vortex"] == 0.0)


def test_calibration_file_lists_only_declared_constants():
    from unistall import dsmodel as dm, calibrate as cal
    c = json.load(open(dm.CONSTANTS_FILE))
    assert set(c["constants"]) <= {f[0] for f in cal.FIT}
    for name, lo, hi in cal.FIT:
        if name in c["constants"]:
            assert lo <= c["constants"][name] <= hi


def test_conditions_outside_the_measured_range_are_refused_when_asked():
    """Above Mach 0.30 there are no measured loops to compare with: the model
    says so, and refuses in strict mode."""
    import pytest
    from unistall import dsmodel as dm
    inside = dm.solve(10.0, 5.0, 0.10, 0.30, n_per_cycle=90, n_cycles=2, strict=True)
    assert inside["meta"]["outside_compared_range"] == []
    for kw in (dict(M=0.45), dict(M=0.10), dict(k=0.35), dict(alpha_mean_deg=20.0, alpha_amp_deg=10.0)):
        case = dict(alpha_mean_deg=10.0, alpha_amp_deg=5.0, k=0.10, M=0.30)
        case.update(kw)
        assert dm.solve(**case, n_per_cycle=90, n_cycles=2)["meta"]["outside_compared_range"]
        with pytest.raises(ValueError):
            dm.solve(**case, n_per_cycle=90, n_cycles=2, strict=True)


def test_tables_are_symmetric_about_the_zero_lift_incidence():
    """A symmetric section: f is even and the static moment odd about the
    measured zero-lift incidence."""
    from unistall.static_model import StaticModel
    sm = StaticModel(0.30)
    a0 = sm.alpha0_deg
    for d in (0.5, 3.0, 9.0, 14.0, 18.0):
        assert abs(sm.f_at(a0 + d) - sm.f_at(a0 - d)) < 1e-12
        zero = sm.cm_at(a0)
        assert abs((sm.cm_at(a0 + d) - zero) + (sm.cm_at(a0 - d) - zero)) < 1e-12
        assert abs(sm.f_at(a0 + d) - float(sm.f(a0 + d))) < 1e-12
        assert abs(sm.cm_at(a0 - d) - float(sm.cm_static(a0 - d))) < 1e-12


def test_first_step_returns_the_attached_loads():
    """The march starts from the loads of the attached state, not from zero."""
    from unistall import dsmodel as dm
    o = dm.solve(8.0, 4.0, 0.10, 0.30, n_per_cycle=360, n_cycles=1)
    assert o["CN"][0] > 0.5 and abs(o["CN"][1] - o["CN"][0]) < 0.1    # the impulsive terms enter at the second step


def test_step_is_set_in_semichords():
    """No march uses a step longer than DS_MAX semichords, at any frequency."""
    import math
    from unistall import dsmodel as dm
    for k in (0.005, 0.01, 0.025, 0.05, 0.1, 0.2, 0.3):
        steps, cycles = dm.march_resolution(k)
        assert 2.0*math.pi/(k*steps) <= dm.DS_MAX + 1e-12
        assert (cycles - 1)*2.0*math.pi/k >= min(dm.SETTLE_SEMICHORDS, (dm.CYCLES_MAX - 1)*2.0*math.pi/k) - 1e-9


def test_loads_do_not_depend_on_chord_or_sound_speed():
    """The model is written in semichords of travel, so stalled loops must be
    the same for any chord and speed of sound. A pitch rate taken from finite
    differences fails this at the turning points of the cycle."""
    import numpy as np
    from unistall import dsmodel as dm
    for case in ((10.0, 10.0, 0.10, 0.30), (15.0, 10.0, 0.10, 0.29), (15.0, 10.0, 0.15, 0.22)):
        a = dm.solve(*case, n_per_cycle=720, n_cycles=4)
        b = dm.solve(*case, n_per_cycle=720, n_cycles=4, chord=1.7, a_sound=295.0)
        for q in ("CL", "CM", "CD"):
            assert np.abs(a[q] - b[q]).max() < 1e-10, (case, q)


def test_pitch_rate_is_exactly_zero_at_the_turning_points():
    import numpy as np
    from unistall import dsmodel as dm
    alpha, rate = dm.prescribed_pitch(10.0, 10.0, 3.0, 720, 1440)
    assert rate[180] == 0.0 and rate[540] == 0.0 and rate[900] == 0.0
    assert rate[179] > 0.0 > rate[181]
    assert np.argmax(alpha[:720]) == 180


def test_normal_force_and_moment_are_odd_about_zero_lift():
    """A loop and its mirror image about the zero-lift incidence: the normal
    force changes sign exactly and the moment does about its static value at
    zero lift. Lift and drag are resolved on the experiment's incidence scale
    and are not expected to mirror."""
    import numpy as np
    from unistall import dsmodel as dm
    from unistall.static_model import StaticModel
    sm = StaticModel(0.30)
    a = dm.solve(10.0, 10.0, 0.10, 0.30, static=sm, n_per_cycle=720, n_cycles=4)
    b = dm.solve(2.0*sm.alpha0_deg - 10.0, -10.0, 0.10, 0.30, static=sm, n_per_cycle=720, n_cycles=4)
    assert np.abs(a["CN"] + b["CN"]).max() < 1e-10
    assert np.ptp(a["CM"] + b["CM"]) < 1e-10
    assert abs((a["CM"] + b["CM"]).mean() - 2.0*sm.cm_at(sm.alpha0_deg)) < 1e-10


def test_every_shedding_in_a_cycle_is_counted():
    """A slow deep-stall loop keeps the onset condition long enough for the
    vortex clock to restart; each restart is a shedding."""
    from unistall import dsmodel as dm
    o = dm.solve(15.0, 10.0, 0.05, 0.30)
    assert o["n_repeated_sheddings"] >= 1 and o["n_sheddings"] == 1 + o["n_repeated_sheddings"]


def test_time_constant_multipliers_follow_the_stated_order():
    """Every combination of the flags, against the order of tests written in
    docs/formulation.md (41) and (42), coded here a second time."""
    import itertools
    from unistall import dsmodel as dm
    from unistall.static_model import StaticModel
    m = dm.DynamicStall(0.30, 100.0, 0.61, 1e-4, None, StaticModel(0.30), 0.25, True, True)
    Tvl = m.p["Tvl"]
    for tesf, lesf, on_chord, vrtx, direction, fpp in itertools.product(
            (False, True), (False, True), (False, True), (False, True), ("away", "returning", "turning"), (0.5, 0.9)):
        away, returning = direction == "away", direction == "returning"
        m.fpp = fpp
        if tesf:
            want = 2.0 if returning else 1.0 if not lesf else 2.0 if fpp <= 0.7 else 1.75
        else:
            want = 1.0
            if not lesf:
                want = 0.5
            if vrtx and on_chord:
                want = 0.25
            if away:
                want = 0.75
        assert m._sigma1(tesf, lesf, on_chord, vrtx, away, returning) == want
        for tv in (0.5*Tvl, 1.5*Tvl, 2.5*Tvl):
            want3 = 1.0
            if Tvl <= tv <= 2.0*Tvl:
                want3 = 3.0 if tesf else 4.0
            if vrtx and on_chord:
                want3 = 2.0 if returning else 1.0
            elif returning:
                want3 = 4.0
            assert m._sigma3(tv, tesf, on_chord, vrtx, returning) == want3


def test_formulation_lists_every_constant_and_the_metric_constants():
    """docs/formulation.md section 5 names every entry of dsmodel.DEFAULTS, and
    the constants of the measures it quotes are those of the code."""
    from unistall import dsmodel as dm, metrics as mt
    from unistall.paths import DOCS
    text = (DOCS/"formulation.md").read_text(encoding="utf-8")
    table = text[text.index("## 5. Constants"):text.index("## 5a.")]
    for name in dm.DEFAULTS:
        assert f"`{name}`" in table, name
    assert f"falls {mt.MSTALL_DROP:g} below" in text
    assert f"{mt.CLOSED_FRACTION:g} of its incidence range" in text
    assert f"{mt.BOOT_N} bootstrap" in text
