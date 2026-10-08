# -*- coding: utf-8 -*-
"""
unistall / check_attached.py
------------------------
The attached-flow acceptance checks of work package 2, written to
results/attached_checks.csv. Nothing here is fitted.

Author: Akosa Samuel Onyejekwe (independent)

  1. THEODORSEN, harmonic pitch about the quarter chord, k = 0.05, 0.10, 0.20,
     at M = 0.05:
       (a) the model, which uses the incompressible form at this Mach number;
       (b) the incompressible form called directly;
       (c) the compressible indicial constants alone, to show how far a set
           fitted for compressible flow sits from the incompressible limit
           and why the blend exists;
       (d) the simplified impulsive-moment form of [DH19] Eq. (1.27), to show
           why it is not used.
  2. CHORD AND SPEED INDEPENDENCE of the Leishman-Beddoes march at fixed (k, M).
  3. MEASURED NO-STALL FRAMES of the calibration set (5 deg +/- 5 deg, M = 0.30):
     normalised RMS error in C_L and C_M, with the static moment taken (i) from
     the report's Table 8 (C_M0 and aerodynamic centre) and (ii) from the
     measured static moment curve of Fig. 9 at the effective incidence.
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
from unistall import attached_flow as af
from unistall import metrics as mt
from unistall.dsmodel import A_SOUND as dm_A_SOUND
from unistall.paths import DATA, RESULTS

T = json.load(open(DATA/"targets.json"))["attached_flow"]
CHORD, A_SOUND = 0.61, dm_A_SOUND       # the tunnel model's chord; speed of sound of dsmodel
NO_STALL_FRAMES = ("frame_10221", "frame_10218")


def _harmonic_response(solver, k, M, npc=2880, ncyc=40, **kw):
    U = M*A_SOUND
    w = 2*k*U/CHORD
    t = np.linspace(0, ncyc*2*np.pi/w, npc*ncyc + 1)
    al = np.radians(1.0)*np.sin(w*t)
    o = solver(al, t, U, **kw)
    s = slice(-npc - 1, None)
    A = af.first_harmonic(al[s], t[s], w)
    return af.first_harmonic(o["CN"][s], t[s], w)/A, af.first_harmonic(o["CM"][s], t[s], w)/A


def theodorsen_rows():
    rows = []
    M = T["theodorsen_M"]
    cna = 2*np.pi/np.sqrt(1 - M*M)
    variants = {
        "the model (low-Mach blend)":
            lambda al, t, U: af.solve_attached(al, t, U, CHORD, M, cna),
        "incompressible form alone":
            lambda al, t, U: af.solve_attached_incompressible(al, t, U, CHORD),
        "compressible indicial constants alone, M = 0.05":
            lambda al, t, U: af.solve_attached(al, t, U, CHORD, M, cna, low_mach_blend=False),
        "compressible constants with the simplified impulsive moment [DH19 Eq. 1.27]":
            lambda al, t, U: af.solve_attached(al, t, U, CHORD, M, cna, moment_alpha_form="quarter",
                                               low_mach_blend=False),
    }
    for name, solver in variants.items():
        for k in T["theodorsen_k"]:
            CN, CM = _harmonic_response(solver, k, T["theodorsen_M"])
            CLt, CMt, _ = af.theodorsen_pitch(k)
            rows.append(dict(check="Theodorsen", variant=name, k=k,
                             lift_amp_err_pct=round(100*(abs(CN)/abs(CLt) - 1), 2),
                             lift_phase_err_deg=round(float(np.degrees(np.angle(CN/CLt))), 2),
                             moment_amp_err_pct=round(100*(abs(CM)/abs(CMt) - 1), 2),
                             moment_phase_err_deg=round(float(np.degrees(np.angle(CM/CMt))), 2)))
    x2 = json.load(open(DATA/"targets.json"))["structural"]["second_pitch_axis_x_c"]
    for k in T["theodorsen_k"]:
        CN, CM = _harmonic_response(lambda al, t, U: af.solve_attached_incompressible(al, t, U, CHORD, x_pitch=x2), k, M)
        CLt, CMt, _ = af.theodorsen_pitch(k, x2)
        rows.append(dict(check="Theodorsen second axis", variant=f"incompressible form, pitch axis at {x2} chord", k=k,
                         lift_amp_err_pct=round(100*(abs(CN)/abs(CLt) - 1), 2),
                         lift_phase_err_deg=round(float(np.degrees(np.angle(CN/CLt))), 2),
                         moment_amp_err_pct=round(100*(abs(CM)/abs(CMt) - 1), 2),
                         moment_phase_err_deg=round(float(np.degrees(np.angle(CM/CMt))), 2)))
    return rows


def independence():
    out = []
    for c, a_s in ((0.30, dm_A_SOUND), (1.70, 295.0)):
        M, k = 0.30, 0.10
        U = M*a_s
        w = 2*k*U/c
        t = np.linspace(0, 6*2*np.pi/w, 6*720 + 1)
        al = np.radians(5.0 + 5.0*np.sin(w*t))
        out.append(af.solve_attached(al, t, U, c, M, 6.245))
    return float(max(np.max(np.abs(out[0]["CN"] - out[1]["CN"])),
                     np.max(np.abs(out[0]["CM"] - out[1]["CM"]))))


def no_stall_rows():
    summ = pd.read_csv(DATA/"static_naca0012_M030_summary.csv").set_index("quantity")["value"]
    cm_static = pd.read_csv(DATA/"static_naca0012_M030_CM.csv")
    pre = cm_static[cm_static.alpha_deg <= 14.2]
    CNa = float(summ["CL_alpha_per_deg"])*180/np.pi
    a0 = np.radians(float(summ["alpha_zero_lift_deg"]))
    rows = []
    for name in NO_STALL_FRAMES:
        fr = mt.load_frame(mt.FRAME_CACHE/f"{name}.mat")
        U = fr["M"]*A_SOUND
        w = 2*fr["k"]*U/CHORD
        npc, ncyc = 1440, 8
        t = np.linspace(0, ncyc*2*np.pi/w, npc*ncyc + 1)
        al = np.radians(fr["a0"] + fr["da"]*np.sin(w*t))
        o = af.solve_attached(al, t, U, CHORD, fr["M"], CNa, alpha0=a0)   # x_ac = 1/4, C_M0 = 0
        s = slice(-npc - 1, None)
        a = al[s]
        CN = o["CN"][s]
        CC = CN*np.tan(a)                       # attached flow: full leading-edge suction
        base = dict(alpha_deg=np.degrees(a), alpha_dot=np.gradient(a, t[s]),
                    CL=CN*np.cos(a) + CC*np.sin(a),
                    CD=CN*np.sin(a) - CC*np.cos(a) + float(summ["CD_min"]))
        for label, cm_s in (
                ("Table 8: C_M0 and aerodynamic centre",
                 float(summ["CM0"]) + o["CN_circ"][s]*(0.25 - float(summ["x_ac_over_c"]))),
                ("Fig. 9: measured static C_M at the effective incidence",
                 np.interp(np.degrees(o["alpha_e"][s] + a0), pre.alpha_deg, pre.Cm_c4))):
            r = mt.score_frame(fr, dict(base, CM=o["CM"][s] + cm_s), set_name="calibration")
            rows.append(dict(check="no-stall frame", variant=label, frame=name, M=r["M"], k=r["k"],
                             nRMS_CL=r["nRMS_CL"], nRMS_CM=r["nRMS_CM"], RMS_CL=r["RMS_CL"],
                             RMS_CM=r["RMS_CM"], range_CM=r["range_CM"],
                             Xi_model=r["Xi_model"], Xi_exp=r["Xi_exp"]))
    return rows


def model_rows():
    """The calibrated model itself on the calibration loops that peak below
    static stall, with the error a perfect model would show from the stated
    measurement uncertainty alone (uncertainty / measured range)."""
    from unistall import dsmodel as dm
    from unistall.uncertainty import U_CL, U_CM, U_DIGIT_CL, U_DIGIT_CM
    frames = [r["frame"] for r in json.load(open(RESULTS/"attached_moment_factor.json"))["frames"]]
    consts = dm.load_constants()
    rows = []
    for name in frames:
        fr = mt.load_frame(mt.FRAME_CACHE/f"{name}.mat")
        r = mt.score_frame(fr, dm.solve(fr["a0"], fr["da"], fr["k"], fr["M"], consts=consts), set_name="calibration")
        rows.append(dict(check="no-stall frame, the model", variant="the calibrated model", frame=name, M=r["M"], k=r["k"],
                         nRMS_CL=r["nRMS_CL"], nRMS_CM=r["nRMS_CM"], RMS_CL=r["RMS_CL"], RMS_CM=r["RMS_CM"],
                         range_CM=r["range_CM"], Xi_model=r["Xi_model"], Xi_exp=r["Xi_exp"],
                         measurement_floor_nRMS_CL=float(np.hypot(U_CL, U_DIGIT_CL)/r["range_CL"]),
                         measurement_floor_nRMS_CM=float(np.hypot(U_CM, U_DIGIT_CM)/r["range_CM"])))
    return rows


if __name__ == "__main__":
    rows = theodorsen_rows()
    d = independence()
    rows.append(dict(check="chord and speed independence", variant="Leishman-Beddoes march",
                     max_abs_difference=d))
    rows += no_stall_rows()
    rows += model_rows()
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS/"attached_checks.csv", index=False)
    th = df[df.check == "Theodorsen"]
    for v, g in th.groupby("variant", sort=False):
        print(f"[attached] {v}: lift {g.lift_amp_err_pct.abs().max():.2f} % / "
              f"{g.lift_phase_err_deg.abs().max():.2f} deg, moment "
              f"{g.moment_amp_err_pct.abs().max():.2f} % / {g.moment_phase_err_deg.abs().max():.2f} deg "
              f"(target {T['lift_and_moment_amplitude_error_max_pct']} % / "
              f"{T['lift_and_moment_phase_error_max_deg']} deg)")
    print(f"[attached] chord/speed independence: {d:.2e} (target {T['chord_speed_independence_max_abs_dCL']:.0e})")
    ns = df[df.check == "no-stall frame"]
    print(ns[["frame", "variant", "nRMS_CL", "nRMS_CM", "RMS_CM", "Xi_model", "Xi_exp"]].to_string(index=False))
