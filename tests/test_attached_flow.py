"""Attached-flow checks that need no measured data.
Author: Akosa Samuel Onyejekwe (independent)"""
import json
from pathlib import Path
import numpy as np
from unistall import attached_flow as af

T = json.load(open(Path(__file__).resolve().parents[1]/"data"/"targets.json"))["attached_flow"]
CHORD = 0.61
from unistall.dsmodel import A_SOUND


def _response(solver, k, U, npc=1440, ncyc=40):
    w = 2*k*U/CHORD
    t = np.linspace(0, ncyc*2*np.pi/w, npc*ncyc + 1)
    al = np.radians(1.0)*np.sin(w*t)
    o = solver(al, t, U)
    s = slice(-npc - 1, None)
    A = af.first_harmonic(al[s], t[s], w)
    return af.first_harmonic(o["CN"][s], t[s], w)/A, af.first_harmonic(o["CM"][s], t[s], w)/A


def test_marching_scheme_reproduces_theodorsen():
    """The deficiency-function march, run with incompressible constants, must
    agree with Theodorsen's closed form inside the frozen target."""
    for x_pitch in (0.25, 0.40):
        for k in T["theodorsen_k"]:
            CN, CM = _response(lambda al, t, U, x_pitch=x_pitch: af.solve_attached_incompressible(
                al, t, U, CHORD, x_pitch=x_pitch), k, 17.0)
            CLt, CMt, _ = af.theodorsen_pitch(k, x_pitch)
            assert abs(abs(CN)/abs(CLt) - 1)*100 <= T["lift_and_moment_amplitude_error_max_pct"]
            assert abs(np.degrees(np.angle(CN/CLt))) <= T["lift_and_moment_phase_error_max_deg"]
            assert abs(abs(CM)/abs(CMt) - 1)*100 <= T["lift_and_moment_amplitude_error_max_pct"]
            assert abs(np.degrees(np.angle(CM/CMt))) <= T["lift_and_moment_phase_error_max_deg"]


def test_model_meets_the_incompressible_limit():
    """At M = 0.05 the model uses the incompressible form and must agree with
    Theodorsen inside the target, in lift and in moment."""
    M = T["theodorsen_M"]
    for k in T["theodorsen_k"]:
        CN, CM = _response(lambda al, t, U: af.solve_attached(
            al, t, U, CHORD, M, 2*np.pi/np.sqrt(1 - M*M)), k, M*A_SOUND)
        CLt, CMt, _ = af.theodorsen_pitch(k)
        assert abs(abs(CN)/abs(CLt) - 1)*100 <= T["lift_and_moment_amplitude_error_max_pct"]
        assert abs(np.degrees(np.angle(CN/CLt))) <= T["lift_and_moment_phase_error_max_deg"]
        assert abs(abs(CM)/abs(CMt) - 1)*100 <= T["lift_and_moment_amplitude_error_max_pct"]
        assert abs(np.degrees(np.angle(CM/CMt))) <= T["lift_and_moment_phase_error_max_deg"]


def test_compressible_constants_alone_miss_the_incompressible_limit():
    """Why the low-Mach blend exists: without it the compressible indicial
    constants are several per cent from Theodorsen at M = 0.05."""
    M = T["theodorsen_M"]
    worst = 0.0
    for k in T["theodorsen_k"]:
        CN, _ = _response(lambda al, t, U: af.solve_attached(
            al, t, U, CHORD, M, 2*np.pi/np.sqrt(1 - M*M), low_mach_blend=False), k, M*A_SOUND)
        CLt, _, _ = af.theodorsen_pitch(k)
        worst = max(worst, abs(abs(CN)/abs(CLt) - 1)*100)
    assert 5.5 <= worst <= 7.5


def test_blend_leaves_the_compressible_range_untouched():
    t = np.linspace(0, 0.3, 1201)
    al = np.radians(5 + 5*np.sin(60*t))
    for M in (0.20, 0.30):
        a = af.solve_attached(al, t, M*A_SOUND, CHORD, M, 6.245)
        b = af.solve_attached(al, t, M*A_SOUND, CHORD, M, 6.245, low_mach_blend=False)
        assert np.max(np.abs(a["CN"] - b["CN"])) == 0.0 and np.max(np.abs(a["CM"] - b["CM"])) == 0.0


def test_march_depends_on_k_and_mach_only():
    out = []
    for c, a_s in ((0.30, A_SOUND), (1.70, 295.0)):
        M, k = 0.30, 0.10
        U = M*a_s
        w = 2*k*U/c
        t = np.linspace(0, 6*2*np.pi/w, 6*720 + 1)
        out.append(af.solve_attached(np.radians(5.0 + 5.0*np.sin(w*t)), t, U, c, M, 6.245))
    d = max(np.max(np.abs(out[0]["CN"] - out[1]["CN"])), np.max(np.abs(out[0]["CM"] - out[1]["CM"])))
    assert d <= T["chord_speed_independence_max_abs_dCL"]


def test_impulsive_time_constant_is_c_over_a():
    M, U, c = 0.30, 102.0, 0.61
    t = np.linspace(0, 0.1, 201)
    o = af.solve_attached(np.radians(5 + 5*np.sin(60*t)), t, U, c, M, 6.245)
    assert abs(o["meta"]["T_I"] - c/(U/M)) < 1e-15


def test_compressible_starting_values_are_those_of_linear_theory():
    """At the Mach numbers of the study the four starting values after a step
    are the exact ones to 0.2 %, and the stored file is what the check gives."""
    import pandas as pd
    from unistall import check_indicial
    from unistall.paths import RESULTS
    for M in check_indicial.MACH:
        for row in check_indicial.starting_values(M):
            assert abs(row["error_pct"]) < 0.2, row
    stored = pd.read_csv(RESULTS/"indicial_checks.csv")
    fresh = pd.DataFrame([r for M in check_indicial.MACH
                          for r in check_indicial.starting_values(M) + [check_indicial.short_time_lift(M)]])
    assert (stored.error_pct - fresh.error_pct).abs().max() < 1e-3
