# -*- coding: utf-8 -*-
"""
unistall / structural.py
--------------------
The two things the load model needs before anything structural can be said:
loads in physical units, and a section that is free to respond.

Author: Akosa Samuel Onyejekwe (independent)

  dimensional_loads   lift and drag per unit span [N/m] and pitching moment
                      per unit span [N m/m] from the coefficients
  cycle_work          work done on the section by the aerodynamic moment over
                      one cycle of prescribed pitch [J/m]
  TorsionalSection    a one-degree-of-freedom section in pitch -- inertia,
                      torsional stiffness, structural damping -- driven by the
                      model's own pitching moment, marched together with the
                      load model
  check               the acceptance measures of work package 13, written to
                      results/structural_checks.csv

WHY THE FREE RESPONSE. A cycle-damping coefficient from PRESCRIBED motion says
whether the air does work on the section at that amplitude and frequency. It
does not show that a section left to itself grows or decays. The torsional
section does: it is released from rest at the top of the stroke and its
amplitude is followed over many cycles.

The section is deliberately heavy (its inertia is INERTIA_RATIO times
q c^2 / w^2), so the aerodynamic moment changes its frequency and amplitude
slowly and the comparison with the prescribed-motion damping at the same
frequency is meaningful. It is a check on the load model, not a model of a
blade: a blade section also flaps and bends, and sees a varying stream.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

from unistall import dsmodel as dm
from unistall import metrics as mt
from unistall import statespace as ss
from unistall.static_model import StaticModel

HERE = Path(__file__).resolve().parent
from unistall.paths import DATA, RESULTS
INERTIA_RATIO = 60.0       # I / (q c^2 / w^2): heavy enough that amplitude changes slowly
FREE_CYCLES = 30           # cycles of free response followed
STEPS_PER_CYCLE = 720
RHO = 1.0                  # kg/m^3, nominal; every check below is a ratio and does not depend on it
#  mean incidence [deg], amplitude [deg], reduced frequency, Mach number
CONDITIONS = [(5.0, 5.0, 0.10, 0.30), (9.0, 5.0, 0.10, 0.30), (11.0, 5.0, 0.10, 0.30),
              (12.0, 5.0, 0.05, 0.30), (15.0, 5.0, 0.10, 0.29), (10.0, 8.0, 0.15, 0.25)]


# a second chord, speed of sound and density for the check that the
# coefficients do not depend on them
OTHER_CHORD_M, OTHER_A_SOUND, OTHER_RHO = 1.7, 295.0, 0.9


def dimensional_loads(CL: np.ndarray, CD: np.ndarray, CM: np.ndarray, rho: float, U: float,
                      chord: float) -> dict:
    """Loads per unit span from coefficients.

    rho [kg/m^3], U [m/s], chord [m]. Returns lift and drag [N/m] and pitching
    moment about the quarter chord [N m/m], with the dynamic pressure [Pa]."""
    q = 0.5*rho*U*U
    return dict(q_Pa=q, lift_N_per_m=q*chord*np.asarray(CL), drag_N_per_m=q*chord*np.asarray(CD),
                moment_Nm_per_m=q*chord*chord*np.asarray(CM))


def cycle_work(alpha_rad: np.ndarray, moment_Nm_per_m: np.ndarray) -> float:
    """Work done ON the section by the aerodynamic moment over the record,
    the integral of M d(alpha) along the time history [J/m]. alpha [rad]."""
    return float(np.sum(0.5*(moment_Nm_per_m[1:] + moment_Nm_per_m[:-1])*np.diff(alpha_rad)))


class TorsionalSection:
    """One degree of freedom in pitch about the quarter chord:
        I a_ddot + 2 zeta I w_n a_dot + I w_n^2 (a - a_m) = M_aero - M_trim .
    M_trim is the mean aerodynamic moment of the prescribed cycle, so the
    section's equilibrium is the mean incidence a_m.

    inertia [kg m^2/m], omega_n [rad/s], zeta [-], alpha_mean [rad],
    trim moment [N m/m], q c^2 [N m/m per unit C_M]."""

    def __init__(self, inertia: float, omega_n: float, zeta: float, alpha_mean: float,
                 moment_trim: float, qc2: float) -> None:
        self.I, self.wn, self.zeta = inertia, omega_n, zeta
        self.alpha_mean, self.M_trim, self.qc2 = alpha_mean, moment_trim, qc2

    def acceleration(self, alpha: float, rate: float, CM: float) -> float:
        """Pitch acceleration [rad/s^2] at incidence alpha [rad], rate [rad/s]
        and aerodynamic moment coefficient CM."""
        return ((self.qc2*CM - self.M_trim)/self.I - 2.0*self.zeta*self.wn*rate
                - self.wn**2*(alpha - self.alpha_mean))


def free_response(alpha_mean_deg: float, alpha_amp_deg: float, k: float, M: float,
                  moment_trim_coeff: float, zeta: float = 0.0) -> dict:
    """Release the torsional section from rest at mean + amplitude and march
    it with the load model. Angles in degrees; k = w_n c / (2U).
    moment_trim_coeff is the mean C_M the spring is trimmed against.
    Returns the peak amplitude of each cycle [deg] and the growth per cycle."""
    c, U = dm.CHORD_M, M*dm.A_SOUND
    wn = 2.0*k*U/c
    dt = 2.0*np.pi/wn/STEPS_PER_CYCLE
    qc2 = 0.5*RHO*U*U*c*c
    sec = TorsionalSection(INERTIA_RATIO*qc2/wn**2, wn, zeta, np.radians(alpha_mean_deg),
                           qc2*moment_trim_coeff, qc2)
    a = np.radians(alpha_mean_deg + alpha_amp_deg)
    rate = 0.0
    model = dm.DynamicStall(M, U, c, dt, static=StaticModel(M))
    rec = model.start(a, a, 0.0)
    peaks, means, hi, lo = [], [], -np.inf, np.inf
    for n in range(1, FREE_CYCLES*STEPS_PER_CYCLE + 1):
        rate += dt*sec.acceleration(a, rate, rec["CM"])     # semi-implicit Euler: rate first,
        a += dt*rate                                         # then incidence with the new rate
        rec = model.step(a, rate)
        hi, lo = max(hi, a), min(lo, a)
        if n % STEPS_PER_CYCLE == 0:
            peaks.append(np.degrees(0.5*(hi - lo)))
            means.append(np.degrees(0.5*(hi + lo)))
            hi, lo = -np.inf, np.inf
    peaks = np.array(peaks)
    late = peaks[5:]                                         # let the starting transient pass
    growth = float(np.polyfit(np.arange(len(late)), np.log(late), 1)[0])
    return dict(amplitude_deg=peaks, mean_deg=np.array(means), log_growth_per_cycle=growth)


def free_response_state_space(alpha_mean_deg: float, alpha_amp_deg: float, k: float, M: float,
                              moment_trim_coeff: float, zeta: float = 0.0) -> dict:
    """The same release as free_response, with the section and the load model
    integrated as ONE system: the state vector is the pitch angle, its rate and
    the twelve aerodynamic states of statespace.py, advanced together by the
    classical fourth-order Runge-Kutta method. Angles in degrees;
    k = w_n c / (2U). Returns the peak amplitude of each cycle [deg] and the
    growth per cycle."""
    c, U = dm.CHORD_M, M*dm.A_SOUND
    wn = 2.0*k*U/c
    qc2 = 0.5*RHO*U*U*c*c
    sec = TorsionalSection(INERTIA_RATIO*qc2/wn**2, wn, zeta, np.radians(alpha_mean_deg),
                           qc2*moment_trim_coeff, qc2)
    aero = ss.StateSpace(M, U, c, static=StaticModel(M))
    h = 2.0*np.pi/(k*STEPS_PER_CYCLE)                     # semichords: w_n t = k s
    march = ss.Marcher(aero, h)
    theta, rate = np.radians(alpha_mean_deg + alpha_amp_deg), 0.0      # rate is d(theta)/ds
    march.start(theta, 0.0, 0.0)

    def derivative(y: np.ndarray, mode: dict) -> np.ndarray:
        """d/ds of (theta, theta', aerodynamic states)."""
        th, v, x = y[0], y[1], y[2:]
        CM = aero.loads(x, th, v, 0.0, march.tau)["CM"]
        acc = sec.acceleration(th, v/aero.per_s, CM)*aero.per_s**2
        return np.concatenate([[v, acc], aero.rhs(x, th, v, acc, mode)])
    peaks, hi, lo = [], -np.inf, np.inf
    for n in range(1, FREE_CYCLES*STEPS_PER_CYCLE + 1):
        y0 = np.concatenate([[theta, rate], march.x])
        mode0 = march.mode(theta, rate, march.lesf)
        k1 = derivative(y0, mode0)
        lesf_end = abs(y0[2 + ss.CNP] + h*k1[2 + ss.CNP]) > aero.CN1
        x0 = march.x.copy()
        march.x = x0 + h*k1[2:]
        march.clock(lesf_end, x0[ss.CNP], theta, theta + h*rate, mode0["returning"])
        march.x = x0
        mode = march.mode(theta + h*rate, rate + h*k1[1], lesf_end)
        k1 = derivative(y0, mode)
        k2 = derivative(y0 + 0.5*h*k1, mode)
        k3 = derivative(y0 + 0.5*h*k2, mode)
        k4 = derivative(y0 + h*k3, mode)
        y = y0 + (h/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4)
        theta, rate, march.x = float(y[0]), float(y[1]), y[2:]
        march.finish(x0, lesf_end)
        out = aero.loads(march.x, theta, rate, 0.0, march.tau)
        if out["CN_vortex"]*out["CNfs"] < 0.0:
            march.x[ss.CNV] = 0.0
        hi, lo = max(hi, theta), min(lo, theta)
        if n % STEPS_PER_CYCLE == 0:
            peaks.append(np.degrees(0.5*(hi - lo)))
            hi, lo = -np.inf, np.inf
    peaks = np.array(peaks)
    late = peaks[5:]
    return dict(amplitude_deg=peaks, log_growth_per_cycle=float(np.polyfit(np.arange(len(late)), np.log(late), 1)[0]))


def check() -> pd.DataFrame:
    """The structural measures for every condition in CONDITIONS.

    work_identity_error_pct compares the cycle work with the same integral from
    an independent march at twice the steps (work_from_damping_J_per_m holds
    that second value). dimensional_identity_max_abs compares the coefficients
    recovered from dimensional loads marched at another chord, speed of sound
    and density. free_over_predicted_growth sets the growth rate of the freely
    pitching section beside the rate that the cycle damping of the prescribed
    motion predicts for small changes of amplitude; in stall the damping
    depends on amplitude, so the two differ as the amplitude changes, and
    free_mean_drift_deg shows whether the mean incidence stayed where the
    spring was trimmed."""
    T = json.load(open(DATA/"targets.json"))["structural"]
    rows = []
    for a0, a1, k, M in CONDITIONS:
        o = dm.solve(a0, a1, k, M)
        U, c = M*dm.A_SOUND, dm.CHORD_M
        d = dimensional_loads(o["CL"], o["CD"], o["CM"], RHO, U, c)
        q = d["q_Pa"]
        # the same condition marched for another chord, speed of sound and density:
        # its dimensional loads, made non-dimensional again, must be these coefficients
        c2, a2, rho2 = OTHER_CHORD_M, OTHER_A_SOUND, OTHER_RHO
        o2 = dm.solve(a0, a1, k, M, chord=c2, a_sound=a2)
        d2 = dimensional_loads(o2["CL"], o2["CD"], o2["CM"], rho2, M*a2, c2)
        ident = max(float(np.max(np.abs(d2["lift_N_per_m"]/(d2["q_Pa"]*c2) - d["lift_N_per_m"]/(q*c)))),
                    float(np.max(np.abs(d2["drag_N_per_m"]/(d2["q_Pa"]*c2) - d["drag_N_per_m"]/(q*c)))),
                    float(np.max(np.abs(d2["moment_Nm_per_m"]/(d2["q_Pa"]*c2*c2) - d["moment_Nm_per_m"]/(q*c*c)))))
        xi = mt.cycle_damping(o["alpha_deg"], o["CM"], a1)
        W = cycle_work(np.radians(o["alpha_deg"]), d["moment_Nm_per_m"])
        # the same work from an independent march with twice as many steps
        fine = dm.solve(a0, a1, k, M, n_per_cycle=2*o["meta"]["steps_per_cycle"])
        W_xi = cycle_work(np.radians(fine["alpha_deg"]), q*c*c*fine["CM"])
        fr = free_response(a0, a1, k, M, float(np.mean(o["CM"][:-1])))
        g = fr["log_growth_per_cycle"]
        g_ss = free_response_state_space(a0, a1, k, M, float(np.mean(o["CM"][:-1])))["log_growth_per_cycle"]
        predicted = -np.pi*xi*np.radians(a1)**2*q*c*c/(INERTIA_RATIO*q*c*c*np.radians(a1)**2)
        rows.append(dict(alpha_mean_deg=a0, alpha_amp_deg=a1, k=k, M=M,
                         cycle_damping=round(xi, 4), cycle_work_J_per_m=round(W, 4),
                         work_from_damping_J_per_m=round(W_xi, 4),
                         work_identity_error_pct=round(100.0*abs(W - W_xi)/max(abs(W_xi), 1e-12), 4),
                         free_log_growth_per_cycle=round(g, 5),
                         predicted_log_growth_per_cycle=round(predicted, 5),
                         free_over_predicted_growth=round(g/predicted, 3) if predicted != 0.0 else float("nan"),
                         free_amplitude_first_deg=round(float(fr["amplitude_deg"][0]), 3),
                         free_amplitude_last_deg=round(float(fr["amplitude_deg"][-1]), 3),
                         free_mean_drift_deg=round(float(fr["mean_deg"][-1] - a0), 3),
                         free_log_growth_per_cycle_state_space=round(g_ss, 5),
                         state_space_agrees_in_sign=bool((g_ss > 0) == (g > 0)),
                         free_response="grows" if g > 0 else "decays",
                         prescribed_damping="negative" if xi < 0 else "positive",
                         agree=bool((g > 0) == (xi < 0)),
                         dimensional_identity_max_abs=ident,
                         peak_lift_N_per_m=round(float(d["lift_N_per_m"].max()), 1),
                         min_moment_Nm_per_m=round(float(d["moment_Nm_per_m"].min()), 2)))
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS/"structural_checks.csv", index=False)
    print(df[["alpha_mean_deg", "alpha_amp_deg", "k", "M", "cycle_damping", "work_identity_error_pct",
              "free_log_growth_per_cycle", "free_log_growth_per_cycle_state_space", "predicted_log_growth_per_cycle",
              "free_response", "agree"]].to_string(index=False))
    print(f"[structural] work identity worst {df.work_identity_error_pct.max():.4f} % "
          f"(target {T['cycle_work_identity_max_error_pct']} %); free response agrees with the sign of the "
          f"prescribed damping in {int(df.agree.sum())} of {len(df)} conditions "
          f"(target >= {T['free_response_sign_agreement_min_conditions']}); dimensional identity "
          f"{df.dimensional_identity_max_abs.max():.1e} (target {T['dimensional_load_identity_max_abs']})")
    return df


if __name__ == "__main__":
    check()
