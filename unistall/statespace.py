# -*- coding: utf-8 -*-
"""
unistall / statespace.py
----------------------
The load model in state-space form: the same equations as dsmodel.py, written
as first-order differential equations in a state vector and integrated by the
classical fourth-order Runge-Kutta method.

Author: Akosa Samuel Onyejekwe (independent)

WHY IT EXISTS. dsmodel.py advances each lag with an exponential recurrence
(the indicial form). Every one of those recurrences is the exact solution,
over one step with a linearly varying input, of a first-order equation
    dz/ds = (u - z)/T ,
so the model can equally be written as dx/ds = F(x, alpha, alpha', alpha'').
This module does that. It serves two purposes: a second integration of the
same equations, switching rules and static tables as the march in dsmodel.py,
which checks the integration and nothing else (results/statespace_check.csv;
an error in a rule or a table would be in both), and a form that can be
coupled to a structural equation and integrated as one system
(structural.py).

TIME is s, semichords of travel: s = 2 U t / c. A prime is d/ds.

STATE VECTOR (N_STATES entries; indices in the constants below)
  W1, W2   the three-quarter-chord incidence lagged at rates b1 beta^2, b2 beta^2
           (circulatory deficiency: X_i = A_i (alpha_34 - W_i))
  W5       the pitch rate q lagged at rate b5 beta^2 (circulatory pitch-rate moment)
  ZA, ZQ   d(alpha)/dt and d(q)/dt lagged with T_alpha, T_q (impulsive normal force)
  Z3, Z4   d(alpha)/dt lagged with T_3, T_4 (impulsive moment)
  ZM       d(q)/dt lagged with k_mq^2 T_I (impulsive pitch-rate moment)
  CNP      the potential normal force lagged with T_p (C'_N)
  FPP      the separation point lagged with T_f (f'')
  AM       the effective incidence lagged with c_m T_f0, for the static moment [deg]
  CNV      the vortex lift
  J1, J2   the three-quarter-chord incidence lagged at the two rates of Jones'
           approximation to the Wagner function (incompressible circulatory
           lift; they act only below M = 0.20, where the loads are blended)

WHAT IS NOT CONTINUOUS. The time constants T_f and T_V take one of several
values according to the flow state (separating or reattaching, vortex on the
chord, incidence rising or falling), and the vortex clock starts, restarts
and re-arms at events. Those rules are the same as in dsmodel.py. They are
evaluated at the start of each step and held over it, and the events are
applied at its end, so across a switch the two forms agree to first order in
the step and not better. Between switches the Runge-Kutta steps are fourth
order.

LOW MACH NUMBER. Below attached_flow.M_COMPRESSIBLE the attached-flow loads
are the blend of dsmodel.py: weight w on the compressible loads above and
1 - w on the incompressible ones, which are the two Jones states for the
circulatory part and the apparent-mass terms of thin-aerofoil theory, written
directly in alpha' and alpha''. So the state-space form covers the same Mach
range as the indicial march.
"""
import math

import numpy as np

from unistall import attached_flow as af
from unistall import dsmodel as dm
from unistall.static_model import StaticModel, F_MIN

SUBSTEP_FRACTION = 1.0     # a Runge-Kutta sub-step is at most this multiple of the fastest lag in play
W1, W2, W5, ZA, ZQ, Z3, Z4, ZM, CNP, FPP, AM, CNV, J1, J2 = range(14)
N_STATES = 14
STATE_NAMES = ("W1", "W2", "W5", "ZA", "ZQ", "Z3", "Z4", "ZM", "CNP", "FPP", "AM", "CNV", "J1", "J2")


class StateSpace:
    """The model as dx/ds = F(x, alpha, alpha', alpha''; mode).

    M [-]; U [m/s]; chord [m]; `consts` a dict of constants (None = the
    calibrated set); `static` a StaticModel (built for M if None); x_pitch as a
    fraction of chord. Incidence in radians; a prime is d/ds with s in
    semichords."""

    def __init__(self, M: float, U: float, chord: float, consts: dict | None = None,
                 static: StaticModel | None = None, x_pitch: float = 0.25) -> None:
        self.p = dm.load_constants() if consts is None else {**dm.DEFAULTS, **consts}
        self.sm = static or StaticModel(M)
        self.M, self.U, self.c, self.xp = M, U, chord, x_pitch
        self.a0 = math.radians(self.sm.alpha0_deg)
        self.CNa = self.sm.CN_alpha
        self.CN1 = self.p["cn1_factor"]*self.sm.CN1
        ref = af.AttachedFlow(U, chord, M, self.CNa, 1.0, alpha0=self.a0, x_pitch=x_pitch).meta
        self.beta2 = 1.0 - M*M
        self.per_s = chord/(2.0*U)                       # seconds per semichord
        self.T_a, self.T_q, self.T_I = ref["T_alpha"], ref["T_q"], ref["T_I"]
        self.T_mq = ref["k_mq"]**2*self.T_I
        self.T3, self.T4 = af.B3*ref["K_alphaM"]*self.T_I, af.B4*ref["K_alphaM"]*self.T_I
        self.b5 = af.B5_DEFAULT
        self.Tcm = max(self.p["cm_lag"]*self.p["Tf0"], dm.CM_LAG_FLOOR)
        self.w = af.compressible_weight(M)                # weight of the compressible loads
        self.g = 1.0 - (1.0 - self.p["cm_unsteady"])*self.w
        self.a_axis = 2.0*x_pitch - 1.0                   # pitch axis aft of mid-chord, semichords

    def fastest_lag(self) -> float:
        """The shortest time constant among the states that carry weight
        [semichords]. The impulsive lags count only where the compressible
        loads do (w > 0)."""
        lags = [self.p["Tp"], self.Tcm, 1.0/(af.B2*self.beta2), 1.0/af.JONES_B2]
        if self.w > 0.0:
            lags += [T/self.per_s for T in (self.T_a, self.T_q, self.T3, self.T4, self.T_mq)]
        return min(lags)

    # ---- inputs -----------------------------------------------------------
    def kinematics(self, alpha: float, d1: float, d2: float) -> tuple:
        """From alpha [rad], alpha' and alpha'' (per semichord, per semichord
        squared): (alpha_34, alpha_34', q, K_alpha [rad/s], K_q [1/s])."""
        q = 2.0*d1                                        # alpha_dot c / U
        a34 = alpha + q*(0.75 - self.xp)
        a34_d = d1 + 2.0*d2*(0.75 - self.xp)
        return a34, a34_d, q, d1/self.per_s, 2.0*d2/self.per_s

    def initial_state(self, alpha: float, d1: float, d2: float) -> np.ndarray:
        """The state in which every lag has caught up with its input."""
        a34, _a34d, q, Ka, Kq = self.kinematics(alpha, d1, d2)
        x = np.zeros(N_STATES)
        x[[W1, W2, J1, J2]] = a34
        x[W5] = q
        x[[ZA, Z3, Z4]] = Ka
        x[[ZQ, ZM]] = Kq
        x[CNP] = self.attached(x, a34, q, d2)["CNpot"]
        af_deg = math.degrees(x[CNP]/self.CNa + self.a0)
        x[FPP] = self.sm.f_at(af_deg)
        x[AM] = af_deg
        return x

    # ---- outputs ----------------------------------------------------------
    def attached(self, x: np.ndarray, a34: float, q: float, d2: float) -> dict:
        """Attached-flow loads from the state: circulatory and impulsive parts,
        blended with the incompressible loads below M_COMPRESSIBLE. d2 is
        alpha'' [1/semichord^2], which the apparent-mass terms need."""
        X1, X2 = af.A1*(a34 - x[W1]), af.A2*(a34 - x[W2])
        out = dict(alpha_e=a34 - self.a0 - X1 - X2)
        out["CNc"] = self.CNa*out["alpha_e"]
        out["CNnc"] = (4.0*self.T_a/self.M)*x[ZA] + (self.T_q/self.M)*x[ZQ]
        out["cm_q"] = -(self.CNa/16.0)*(q - af.A5*(q - x[W5]))
        out["cm_nc"] = (-(af.A3*self.T3/self.M)*x[Z3] - (af.A4*self.T4/self.M)*x[Z4]
                        - (7.0*self.T_mq/(12.0*self.M))*x[ZM])
        if self.w < 1.0:
            a = self.a_axis
            inc = dict(alpha_e=a34 - self.a0 - af.JONES_A1*(a34 - x[J1]) - af.JONES_A2*(a34 - x[J2]), cm_q=0.0)
            inc["CNc"] = self.CNa*inc["alpha_e"]
            inc["CNnc"] = math.pi*(q/2.0 - a*d2)
            inc["cm_nc"] = ((math.pi/2.0)*(-(0.5 - a)*q/2.0 - (0.125 + a*a)*d2) + (a + 0.5)*inc["CNc"]/2.0
                            - (inc["CNc"] + inc["CNnc"])*(self.xp - 0.25))
            out = {key: self.w*out[key] + (1.0 - self.w)*inc[key] for key in out}
        out["CNpot"] = out["CNc"] + out["CNnc"]
        return out

    def loads(self, x: np.ndarray, alpha: float, d1: float, d2: float, tau: float) -> dict:
        """C_N, C_C, C_M (c/4) and the quantities behind them, from the state,
        the incidence [rad] and the vortex clock tau [semichords]."""
        p = self.p
        a34, _a34d, q, _Ka, _Kq = self.kinematics(alpha, d1, d2)
        a = self.attached(x, a34, q, d2)
        fpp = min(max(x[FPP], F_MIN), 1.0)
        Kf = ((1.0 + math.sqrt(fpp))/2.0)**2
        af_deg = math.degrees(x[CNP]/self.CNa + self.a0)
        fprime = self.sm.f_at(af_deg)
        cc_angle = alpha if p["cc_geometric"] else a["alpha_e"] + self.a0
        suction = a["CNc"]*math.tan(cc_angle)
        if p["cc_floor"] and suction < 0.0:
            suction = 0.0
        cc = p["eta"]*suction*(math.sqrt(fpp) - p["cc_offset"])
        if abs(x[CNP]) > self.CN1:
            cc *= fpp**max(2.0*(abs(x[CNP]) - self.CN1) + fpp - fprime, 0.0)
        cpv = p["xcp_v"]*(1.0 - math.cos(math.pi*min(tau/p["Tvl"], 1.0))) if tau > 0.0 else 0.0
        cm_q = a["cm_q"]*(Kf if p["cmq_kirchhoff"] else 1.0)
        CM = self.sm.cm_at(x[AM]) + self.g*(cm_q + a["cm_nc"]) - cpv*x[CNV]
        return dict(CN=a["CNc"]*Kf + a["CNnc"] + x[CNV], CC=cc, CM=CM, CN_prime=x[CNP], f_sep=fpp,
                    f_static=fprime, CN_vortex=x[CNV], alpha_f_deg=af_deg, CNfs=a["CNc"]*Kf + a["CNnc"])

    # ---- the equations ----------------------------------------------------
    def rhs(self, x: np.ndarray, alpha: float, d1: float, d2: float, mode: dict) -> np.ndarray:
        """dx/ds for incidence alpha [rad], alpha' and alpha''; `mode` holds the
        time constants in force (Tf, Tv) and whether the vortex is fed."""
        a34, a34_d, q, Ka, Kq = self.kinematics(alpha, d1, d2)
        a = self.attached(x, a34, q, d2)
        dx = np.empty(N_STATES)
        dx[J1], dx[J2] = af.JONES_B1*(a34 - x[J1]), af.JONES_B2*(a34 - x[J2])
        r1, r2 = af.B1*self.beta2, af.B2*self.beta2
        dx[W1], dx[W2] = r1*(a34 - x[W1]), r2*(a34 - x[W2])
        dx[W5] = self.b5*self.beta2*(q - x[W5])
        if self.w > 0.0:
            dx[ZA], dx[ZQ] = self.per_s*(Ka - x[ZA])/self.T_a, self.per_s*(Kq - x[ZQ])/self.T_q
            dx[Z3], dx[Z4] = self.per_s*(Ka - x[Z3])/self.T3, self.per_s*(Ka - x[Z4])/self.T4
            dx[ZM] = self.per_s*(Kq - x[ZM])/self.T_mq
        else:                                             # the impulsive lags carry no weight: left at rest
            dx[[ZA, ZQ, Z3, Z4, ZM]] = 0.0
        dx[CNP] = (a["CNpot"] - x[CNP])/self.p["Tp"]
        af_deg = math.degrees(x[CNP]/self.CNa + self.a0)
        dx[FPP] = (self.sm.f_at(af_deg) - x[FPP])/mode["Tf"]
        dx[AM] = (af_deg - x[AM])/self.Tcm
        dx[CNV] = -x[CNV]/mode["Tv"]
        if mode["feed"]:
            fpp = min(max(x[FPP], F_MIN), 1.0)
            root = math.sqrt(fpp)
            Kf = ((1.0 + root)/2.0)**2
            dKf = ((1.0 + root)/(4.0*root))*dx[FPP] if F_MIN < x[FPP] < 1.0 else 0.0
            dCNc = self.w*self.CNa*(a34_d - af.A1*(a34_d - dx[W1]) - af.A2*(a34_d - dx[W2]))
            dCNc += (1.0 - self.w)*self.CNa*(a34_d - af.JONES_A1*(a34_d - dx[J1]) - af.JONES_A2*(a34_d - dx[J2]))
            dx[CNV] += dCNc*(1.0 - Kf) - a["CNc"]*dKf
        return dx


class Marcher:
    """StateSpace advanced step by step with the switching rules of dsmodel.py.

    `ds` is the step [semichords]. Call `start` once, then `step` with the
    motion, a function of s returning (alpha [rad], alpha', alpha''), and the
    value of s at the start of the step. The step is divided into equal
    Runge-Kutta sub-steps where the fastest lag in play needs it (the
    impulsive lags shorten with Mach number), the switching rules being
    applied once per step as in dsmodel.py."""

    def __init__(self, model: StateSpace, ds: float, vortex: bool = True) -> None:
        self.m, self.ds, self.vortex = model, ds, vortex
        self.substeps = max(1, math.ceil(ds/(SUBSTEP_FRACTION*model.fastest_lag())))
        self.fed = dm.DynamicStall(model.M, model.U, model.c, 1.0, model.p, model.sm, model.xp)   # for its switching rules only

    def start(self, alpha: float, d1: float, d2: float) -> dict:
        """Set the state at the first instant; returns the loads there."""
        self.x = self.m.initial_state(alpha, d1, d2)
        self.tau, self.lesf, self.dfpp = 0.0, False, 0.0
        return self._record(alpha, d1, d2, float("nan"))

    def mode(self, alpha: float, d1: float, lesf: bool) -> dict:
        """Time constants in force over the coming step, by the rules of dsmodel.py."""
        p, x = self.m.p, self.x
        d0 = alpha - self.m.a0
        away, returning = d1*d0 > 0.0, d1*d0 < 0.0
        tesf = self.dfpp < 0.0
        on_chord, vrtx = 0.0 < self.tau <= p["Tvl"], 0.0 < self.tau <= 2.0*p["Tvl"]
        self.fed.fpp = min(max(x[FPP], F_MIN), 1.0)
        Tf = p["Tf0"]/self.fed._sigma1(tesf, lesf, on_chord, vrtx, away, returning)
        Tv = p["Tv0"]/self.fed._sigma3(self.tau, tesf, on_chord, vrtx, returning)
        return dict(Tf=Tf, Tv=Tv, feed=bool(self.vortex and on_chord), returning=returning)

    def clock(self, lesf: bool, cnp_before: float, alpha0: float, alpha1: float, returning: bool) -> float:
        """Advance the vortex clock over the step just taken; returns the onset
        incidence [deg] if the step contains one, else nan."""
        p, onset = self.m.p, float("nan")
        if not self.vortex:
            self.tau = 0.0
        elif self.tau > 0.0:
            self.tau += self.ds
        elif lesf and not self.lesf:
            x0, x1 = abs(cnp_before) - self.m.CN1, abs(self.x[CNP]) - self.m.CN1
            frac = min(max(x1/(x1 - x0) if x1 != x0 else 1.0, 0.0), 1.0)
            self.tau = frac*self.ds
            onset = math.degrees(alpha1 - (alpha1 - alpha0)*frac)
        fpp = min(max(self.x[FPP], F_MIN), 1.0)
        if self.tau >= p["Tvl"] + 2.0*(1.0 - fpp)/p["St_sh"]:
            self.tau = self.ds if lesf else 0.0
        if self.tau > p["Tvl"] and (not lesf) and returning:
            self.tau = 0.0
        return onset

    def step(self, motion: object, s0: float) -> dict:
        """One step from s0 to s0 + ds."""
        m, h, x0 = self.m, self.ds, self.x
        begin, end = motion(s0), motion(s0 + h)
        # the vortex clock is advanced first, as dsmodel.py does, from the
        # onset state predicted at the end of the step by a first-order look
        mode0 = self.mode(begin[0], begin[1], self.lesf)
        k1 = m.rhs(x0, *begin, mode0)
        lesf_end = abs(x0[CNP] + h*k1[CNP]) > m.CN1
        self.x = x0 + h*k1                                   # provisional, for the clock's crossing
        onset = self.clock(lesf_end, x0[CNP], begin[0], end[0], mode0["returning"])
        mode = self.mode(end[0], end[1], lesf_end)
        x, hs = x0, h/self.substeps
        for i in range(self.substeps):
            sa = s0 + i*hs
            a, mid, b = motion(sa), motion(sa + 0.5*hs), motion(sa + hs)
            k1 = m.rhs(x, *a, mode)
            k2 = m.rhs(x + 0.5*hs*k1, *mid, mode)
            k3 = m.rhs(x + 0.5*hs*k2, *mid, mode)
            k4 = m.rhs(x + hs*k3, *b, mode)
            x = x + (hs/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4)
        self.x = x
        self.finish(x0, lesf_end)
        out = self._record(*end, onset)
        if out["CN_vortex"]*out["CNfs"] < 0.0:              # vortex lift never opposes the attached force
            self.x[CNV] = 0.0
            out = self._record(*end, onset)
        return out

    def finish(self, x0: np.ndarray, lesf_end: bool) -> None:
        """Book-keeping after a step that started from state x0 [-]: the
        direction in which the separation point moved, and the onset flag."""
        self.dfpp = self.x[FPP] - x0[FPP]
        self.lesf = lesf_end

    def _record(self, alpha: float, d1: float, d2: float, onset: float) -> dict:
        out = self.m.loads(self.x, alpha, d1, d2, self.tau)
        out.update(alpha=alpha, tau_v=self.tau, onset_alpha_deg=onset)
        return out


def solve(alpha_mean_deg: float, alpha_amp_deg: float, k: float, M: float, consts: dict | None = None,
          static: StaticModel | None = None, n_per_cycle: int | None = None, n_cycles: int | None = None,
          chord: float = dm.CHORD_M, a_sound: float = dm.A_SOUND, x_pitch: float = 0.25,
          vortex: bool = True) -> dict:
    """March alpha = mean + amp sin(k s) in state-space form and return the
    last cycle, with the same keys as dsmodel.solve.

    Angles in degrees; k [-]; M [-]; chord [m]; a_sound [m/s]. The step and the number of cycles default to
    dsmodel.march_resolution(k)."""
    auto = dm.march_resolution(k)
    n_per_cycle, n_cycles = n_per_cycle or auto[0], n_cycles or auto[1]
    U = M*a_sound
    model = StateSpace(M, U, chord, consts, static, x_pitch)
    N = n_per_cycle*n_cycles
    h = 2.0*math.pi/(k*n_per_cycle)
    a_m, a_1 = math.radians(alpha_mean_deg), math.radians(alpha_amp_deg)

    def motion(s: float) -> tuple:
        """(alpha, alpha', alpha'') at s."""
        return a_m + a_1*math.sin(k*s), a_1*k*math.cos(k*s), -a_1*k*k*math.sin(k*s)
    march = Marcher(model, h, vortex)
    recs = [march.start(*motion(0.0))]
    recs += [march.step(motion, n*h) for n in range(N)]
    sl = slice(N - n_per_cycle, N + 1)

    def col(key: str) -> np.ndarray:
        """One recorded quantity over the last cycle."""
        return np.array([r[key] for r in recs[sl]])
    alpha = col("alpha")
    s = h*np.arange(N + 1)[sl]
    CN, CC = col("CN"), col("CC")
    on = col("onset_alpha_deg")
    on = on[np.isfinite(on)]
    tau = col("tau_v")
    restarts = int(np.sum((tau[1:] < tau[:-1]) & (tau[1:] > 0.0)))      # the vortex clock restarted under sustained onset
    return dict(t=(s - s[0])*model.per_s, alpha_deg=np.degrees(alpha), alpha_dot=a_1*k*np.cos(k*s)/model.per_s,
                CL=CN*np.cos(alpha) + CC*np.sin(alpha), CD=CN*np.sin(alpha) - CC*np.cos(alpha) + model.p["CD0"],
                CM=col("CM"), CN=CN, CC=CC, CN_prime=col("CN_prime"), f_sep=col("f_sep"), f_static=col("f_static"),
                CN_vortex=col("CN_vortex"), tau_v=col("tau_v"), alpha_f_deg=col("alpha_f_deg"),
                onset_alpha_deg=float(on[0]) if len(on) else float("nan"),
                n_sheddings=int(len(on)) + restarts, n_repeated_sheddings=restarts,
                meta=dict(M=M, k=k, steps_per_cycle=n_per_cycle, cycles=n_cycles, ds=h, states=STATE_NAMES,
                          substeps=march.substeps,
                          integrator="classical fourth-order Runge-Kutta, fixed step", constants=model.p))
