# -*- coding: utf-8 -*-
"""
unistall / dsmodel.py
-----------------
Dynamic-stall load model for a pitching aerofoil: the Leishman-Beddoes method
with the separation point and the static moment read from measured static data
(a table) instead of fitted curves.

Author: Akosa Samuel Onyejekwe (independent)

Formulation: Damiani & Hayman (2019), NREL/TP-5000-66347, doi:10.2172/1576488
([DH19]); equation numbers below are theirs. The attached-flow part is
unistall/attached_flow.py; the static inputs are unistall/static_model.py.

Per time step n:

  attached flow       C_N^pot = C_N^c + C_N^nc ;  alpha_e                 (1.14)-(1.20)
  pressure lag        C'_N = C_N^pot - D_p ,  time constant T_p           (1.35)
  effective angle     alpha_f = C'_N / C_Nalpha + alpha_0                  (1.34)
  static separation   f' = f(alpha_f ; M)      looked up, not fitted      (text after 1.36)
  boundary-layer lag  f'' = f' - D_f ,  T_f = T_f0 / sigma_1               (1.36)-(1.37)
  separated force     C_N^fs = C_N^c ((1+sqrt f'')/2)^2 + C_N^nc           (1.38)
  onset               leading-edge separation when |C'_N| > C_N1(M)        (1.46)
  vortex feed         C_V = C_N^c (1 - ((1+sqrt f'')/2)^2)                 (1.49)
  vortex force        C_N^v, time constant T_V = T_V0 / sigma_3            (1.47), (1.52)
  vortex clock        tau_v in semichords from onset; the vortex is over
                      the chord while 0 < tau_v <= T_VL                    (1.51)
  repeated shedding   tau_v restarts at T_VL + T_sh,
                      T_sh = 2 (1 - f'') / St_sh                           (1.54)
  normal force        C_N = C_N^fs + C_N^v                                 (1.53)
  chord force         C_C = eta C_N^c tan(alpha) sqrt f'' , alpha the geometric
                      incidence (alpha_e + alpha_0 with cc_geometric = 0),
                      times f''^k2 with k2 = 2(C'_N - C_N1) + f'' - f'
                      once |C'_N| > C_N1                                   (1.21), (1.32), (1.56)
  moment about c/4    C_M = C_M,static(alpha'_f) + C_M^c,q + C_M^nc + C_M^v (1.59)
                      alpha'_f = alpha_f lagged with 0.1 T_f0              (1.43)
                      C_M^v = -x_cp,v (1 - cos(pi tau_v / T_VL)) C_N^v     (1.57)
  lift and drag       C_L = C_N cos a + C_C sin a
                      C_D = C_N sin a - C_C cos a + C_D0                   (1.2)

The multipliers sigma_1 and sigma_3 follow the rules of [DH19] section
2.6.2.3: they speed the separation point up while separation is in progress,
slow it down during reattachment, hold it back while a vortex is over the
chord, and speed the decay of vortex lift once the vortex has left.

DEPARTURES FROM THE SOURCE, each stated (docs/formulation.md gives the reason):
  * eta multiplies the chord force on both sides of the onset threshold. As
    printed, (1.56a) drops it above the threshold, which puts a step of
    (1 - eta) in the chord force at onset.
  * the onset crossing is located inside the time step by linear
    interpolation of C'_N, so the vortex clock starts at the crossing and not
    at the end of the step in which it happened.
  * the leading-edge suction is resolved along the geometric incidence
    (cc_geometric), not the lagged effective one.
  * the circulatory pitch-rate moment is reduced by the Kirchhoff factor
    (cmq_kirchhoff).
  * the unsteady moment terms carry an empirical factor taken from
    loops below static stall (cm_unsteady).
  * below M = 0.20 the attached-flow loads are blended with the
    incompressible ones (attached_flow.py), b_5 = 0.5 and the two-term
    impulsive moment are used there, the exponent k2 of (1.56) is not allowed
    below zero, and vortex lift that would oppose the attached force is set
    to zero.

WHAT IS CALIBRATED AND WHAT IS NOT. Slope, zero-lift angle, C_N1 and f come
from static measurements (static_model.py) and are never fitted. The dynamic
constants in DEFAULTS are starting values; unistall/calibrate.py fits a subset of
them on the calibration frames and writes the result to a file this module
reads through `load_constants`.
"""
import json, math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
from unistall import attached_flow as af
from unistall.static_model import StaticModel, F_MIN
from unistall.paths import RESULTS

CHORD_M = 0.61            # model chord of the NASA TM-84245 experiment
A_SOUND = math.sqrt(1.4*287.05*288.15)   # sea-level standard air, m/s; the coefficients depend on (k, M),
                                         # not on chord or sound speed (tested)
CM_LAG_FLOOR = 1e-6       # keeps the moment lag defined when cm_lag is set to zero

# The conditions the measured loops cover (data/data_manifest.csv, NACA 0012
# frames inside the static Mach range). Outside them the model still runs, but
# nothing it returns has been compared with a measurement.
K_MAX_COMPARED = 0.20            # highest reduced frequency among those frames
PEAK_ALPHA_MAX_COMPARED = 25.0   # deg; the static lift data end here

# Resolution of the march. The lags are measured in semichords of travel, so
# the step is set in semichords, not as a fraction of the cycle: at low reduced
# frequency one cycle is hundreds of semichords long. results/convergence_*.csv
# shows what these two numbers buy.
DS_MAX = 0.02                    # largest step, semichords of travel: the step at which the loop errors,
                                 # extrapolated to zero step, are within the target of data/targets.json
SETTLE_SEMICHORDS = 200.0        # travel marched before the reported cycle
STEPS_MIN, STEPS_MULTIPLE = 720, 360
CYCLES_MIN, CYCLES_MAX = 3, 8


def march_resolution(k: float, ds_max: float = DS_MAX, settle: float = SETTLE_SEMICHORDS) -> tuple:
    """(steps per cycle, cycles marched) for reduced frequency k [-]: the step
    is at most `ds_max` semichords and at least `settle` semichords of travel
    precede the reported (last) cycle."""
    per_cycle = 2.0*math.pi/k                             # semichords in one cycle: s = 2 U t / c
    steps = max(STEPS_MIN, STEPS_MULTIPLE*math.ceil(per_cycle/ds_max/STEPS_MULTIPLE))
    cycles = int(min(max(math.ceil(settle/per_cycle) + 1, CYCLES_MIN), CYCLES_MAX))
    return steps, cycles

# Multipliers on 1/T_f and 1/T_V, and the separation point below which
# separation is accelerated. All from [DH19] section 2.6.2.3; none is fitted.
F_FAST_SEPARATION = 0.7
SIGMA1 = dict(separating_returning=2.0, separating_attached_le=1.0, separating_low_f=2.0,
              separating_high_f=1.75, default=1.0, reattaching=0.5, vortex_on_chord=0.25,
              reattaching_away=0.75)
SIGMA3 = dict(default=1.0, shed_separating=3.0, shed_reattaching=4.0, on_chord_away=1.0,
              on_chord_returning=2.0, returning=4.0)

DEFAULTS = dict(
    Tp=1.7,               # pressure lag, semichords
    Tf0=3.0,              # boundary-layer lag
    Tv0=6.0,              # vortex-lift decay
    Tvl=7.0,              # vortex travel time over the chord; [DH19] gives the range 6-13
    eta=0.95,             # chord-force recovery factor; [DH19] gives 0.85-0.95
    xcp_v=0.20,           # vortex centre-of-pressure travel, [DH19] Eq. (1.57)
    St_sh=0.19,           # Strouhal number of repeated shedding
    cm_lag=0.10,          # fraction of T_f0 that lags alpha_f for the moment, [DH19] Eq. (1.43)
    CD0=0.0072,           # minimum static drag, NASA TM-84245 Vol. 1 Table 8
    cn1_factor=1.0,       # onset level as a multiple of the static-stall normal force
    cmq_kirchhoff=1,      # 1 = circulatory pitch-rate moment scaled by the Kirchhoff factor
    cm_unsteady=1.0,      # factor on the compressible unsteady moment terms (1 = inviscid theory);
                          # fitted on loops below static stall by calibrate.py --attached
    sigma_rules=1,        # 1 = time-constant multipliers of [DH19] 2.6.2.3; 0 = none (T_V halved once the vortex has left)
    cc_geometric=1,       # 1 = chord force uses the geometric incidence, 0 = the lagged effective one
    cc_offset=0.0,        # 0 = Leishman-Beddoes chord force; 0.2 = Gonzalez's form, [DH19] Eq. (1.40)
    cc_floor=1,           # 1 = the leading-edge suction C_N^c tan(alpha) is not allowed below zero (it cannot
                          # push); 0 = as the formula gives it
)
CONSTANTS_FILE = RESULTS/"calibrated_constants.json"


def load_constants(file: Path = CONSTANTS_FILE) -> dict:
    """DEFAULTS overridden by a calibration record, if it exists: first the
    values it held fixed (the attached-flow moment factor among them), then
    the constants it fitted. `file` is the record's path."""
    p = dict(DEFAULTS)
    if Path(file).exists():
        with open(file) as fh:
            cal = json.load(fh)
        p.update({k: v for k, v in cal.get("fixed", {}).items() if k in DEFAULTS})
        p.update(cal["constants"])
    return p


class DynamicStall:
    """The load model advanced one time step at a time, for a motion that is
    not known in advance (a prescribed motion or a structural response).

    M [-], k-independent; U [m/s]; chord [m]; dt [s]; `consts` as in DEFAULTS
    (None = the calibrated set); `static` a StaticModel (built for M if None);
    x_pitch as a fraction of chord. Call `start(alpha0, alpha1, alpha_dot0)`
    once and then `step(alpha, alpha_dot)`; incidences in radians, alpha_dot
    [rad/s] is used for its sign only. Each call returns the loads and states
    at that step as a dict of floats.
    """

    def __init__(self, M: float, U: float, chord: float, dt: float, consts: dict | None = None,
                 static: StaticModel | None = None, x_pitch: float = 0.25,
                 vortex: bool = True, separation: bool = True) -> None:
        self.p = load_constants() if consts is None else {**DEFAULTS, **consts}
        self.sm = static or StaticModel(M)
        self.vortex, self.separation = vortex, separation
        self.a0 = np.radians(self.sm.alpha0_deg)
        # critical normal force for leading-edge separation: the static-stall
        # value times cn1_factor. The factor allows the onset level to sit above
        # the quasi-steady stall value, as it does in the original model.
        self.CN1 = self.p["cn1_factor"]*self.sm.CN1
        self.att = af.AttachedFlow(U, chord, M, self.sm.CN_alpha, dt, alpha0=self.a0, x_pitch=x_pitch)
        self.ds = self.att.meta["ds"]
        # the strength of the factor follows the blend weight w: it multiplies
        # the blended unsteady moment by 1 - (1 - g_M) w, which is g_M at and
        # above M_COMPRESSIBLE and 1 at and below M_INCOMPRESSIBLE, so the
        # incompressible limit (Theodorsen) is untouched
        self.g_unsteady = 1.0 - (1.0 - self.p["cm_unsteady"])*af.compressible_weight(M)
        p = self.p
        self.Etp, self.Etph = math.exp(-self.ds/p["Tp"]), math.exp(-self.ds/(2*p["Tp"]))
        Tcm = max(p["cm_lag"]*p["Tf0"], CM_LAG_FLOOR)
        self.Ecm, self.Ecmh = math.exp(-self.ds/Tcm), math.exp(-self.ds/(2*Tcm))

    def _f_static(self, alpha_f_deg: float) -> float:
        """Static separation point at the effective incidence [deg]."""
        return self.sm.f_at(alpha_f_deg) if self.separation else 1.0

    def start(self, alpha_first: float, alpha_second: float, alpha_dot: float) -> dict:
        """State and loads at the first instant."""
        a = self.att.start(alpha_first, alpha_second)
        self.alpha = alpha_first
        self.CNpot = self.CNp = a["CN"]
        self.af_deg = float(math.degrees(self.CNp/self.sm.CN_alpha + self.a0))
        self.fp = self.fpp = self._f_static(self.af_deg)
        self.fpp_before = None
        self.Dp = self.Df = self.Daf = self.Cv = self.CNv = self.tau = 0.0
        self.lesf = False
        Kf = ((1.0 + math.sqrt(self.fpp))/2.0)**2
        cc = self.p["eta"]*a["CN"]*math.tan(alpha_first)*(math.sqrt(self.fpp) - self.p["cc_offset"])
        return dict(alpha=alpha_first, alpha_dot=alpha_dot, CN=a["CN"]*Kf, CC=cc,
                    CM=self.sm.cm_at(self.af_deg), CN_prime=self.CNp,
                    f_sep=self.fpp, f_static=self.fp, CN_vortex=0.0, tau_v=0.0, alpha_f_deg=self.af_deg,
                    onset_alpha_deg=float("nan"), CN_vortex_removed=0.0, reshed=False)

    def _vortex_clock(self, lesf: bool, CNp: float, alpha: float, returning: bool) -> tuple:
        """Advance tau_v [semichords]; returns (tau_v, onset incidence [deg] or nan).
        The onset is located inside the step by linear interpolation of |C'_N|."""
        p, tv, onset = self.p, self.tau, float("nan")
        if not self.vortex:
            tv = 0.0
        elif tv > 0.0:
            tv += self.ds
        elif lesf and not self.lesf:
            x0, x1 = abs(self.CNp) - self.CN1, abs(CNp) - self.CN1
            frac = min(max(x1/(x1 - x0) if x1 != x0 else 1.0, 0.0), 1.0)   # part of the step past the crossing
            tv = frac*self.ds
            onset = float(math.degrees(alpha - (alpha - self.alpha)*frac))
        Tsh = 2.0*(1.0 - self.fpp)/p["St_sh"]
        self.reshed = False
        if tv >= p["Tvl"] + Tsh:
            self.reshed = bool(lesf)
            tv = self.ds if lesf else 0.0                    # shed again, or re-arm
        if tv > p["Tvl"] and (not lesf) and returning:
            tv = 0.0                                         # vortex gone, incidence falling: re-arm
        return tv, onset

    def _sigma1(self, tesf: bool, lesf: bool, on_chord: bool, vrtx: bool, away: bool, returning: bool) -> float:
        """Multiplier on 1/T_f, [DH19] section 2.6.2.3.1."""
        if not self.p["sigma_rules"]:
            return SIGMA1["default"]
        if tesf:
            if returning:
                return SIGMA1["separating_returning"]
            if not lesf:
                return SIGMA1["separating_attached_le"]
            return SIGMA1["separating_low_f"] if self.fpp <= F_FAST_SEPARATION else SIGMA1["separating_high_f"]
        s1 = SIGMA1["default"]
        if not lesf:
            s1 = SIGMA1["reattaching"]
        if vrtx and on_chord:
            s1 = SIGMA1["vortex_on_chord"]
        if away:
            s1 = SIGMA1["reattaching_away"]
        return s1

    def _sigma3(self, tv: float, tesf: bool, on_chord: bool, vrtx: bool, returning: bool) -> float:
        """Multiplier on 1/T_V, [DH19] section 2.6.2.3.2."""
        if not self.p["sigma_rules"]:
            return SIGMA3["default"] if (vrtx and on_chord) else SIGMA3["on_chord_returning"]
        s3 = SIGMA3["default"]
        if self.p["Tvl"] <= tv <= 2.0*self.p["Tvl"]:
            s3 = SIGMA3["shed_separating"] if tesf else SIGMA3["shed_reattaching"]
        if vrtx and on_chord:
            s3 = SIGMA3["on_chord_returning"] if returning else SIGMA3["on_chord_away"]
        elif returning:
            s3 = SIGMA3["returning"]
        return s3

    def step(self, alpha: float, alpha_dot: float) -> dict:
        """Advance to the next instant: incidence alpha [rad], rate alpha_dot [rad/s]."""
        p, sm = self.p, self.sm
        a = self.att.step(alpha)
        CNc, CNnc = a["CN_circ"], a["CN_imp_alpha"] + a["CN_imp_q"]
        CNpot = CNc + CNnc
        d0 = alpha - self.a0                                 # incidence from zero lift
        away, returning = alpha_dot*d0 > 0.0, alpha_dot*d0 < 0.0
        # ---- pressure lag and effective incidence
        self.Dp = self.Dp*self.Etp + (CNpot - self.CNpot)*self.Etph
        CNp = CNpot - self.Dp
        af_deg = float(math.degrees(CNp/sm.CN_alpha + self.a0))
        fprime = self._f_static(af_deg)
        lesf = abs(CNp) > self.CN1
        tesf = (self.fpp < self.fpp_before) if self.fpp_before is not None else False
        tv, onset = self._vortex_clock(lesf, CNp, alpha, returning)
        on_chord, vrtx = 0.0 < tv <= p["Tvl"], 0.0 < tv <= 2.0*p["Tvl"]
        # ---- boundary-layer lag and separated normal force
        Tf = p["Tf0"]/self._sigma1(tesf, lesf, on_chord, vrtx, away, returning)
        self.Df = self.Df*math.exp(-self.ds/Tf) + (fprime - self.fp)*math.exp(-self.ds/(2*Tf))
        fpp = min(max(fprime - self.Df, F_MIN), 1.0) if self.separation else 1.0
        Kf = ((1.0 + math.sqrt(fpp))/2.0)**2
        CNfs = CNc*Kf + CNnc
        # ---- vortex lift
        Tv = p["Tv0"]/self._sigma3(tv, tesf, on_chord, vrtx, returning)
        Cv = CNc*(1.0 - Kf)
        CNv = self.CNv*math.exp(-self.ds/Tv)
        if self.vortex and on_chord:
            CNv += (Cv - self.Cv)*math.exp(-self.ds/(2*Tv))
        removed = 0.0
        if CNv*CNfs < 0.0:
            removed, CNv = CNv, 0.0                          # never opposes the attached force
        # ---- chord force
        # direction of the leading-edge suction: the lagged effective incidence
        # (potential theory, which gives unsteady thrust on the return stroke)
        # or the geometric incidence (no thrust beyond the quasi-steady value)
        cc_angle = alpha if p["cc_geometric"] else a["alpha_e"] + self.a0
        suction = CNc*math.tan(cc_angle)
        if p["cc_floor"] and suction < 0.0:
            suction = 0.0
        cc = p["eta"]*suction*(math.sqrt(fpp) - p["cc_offset"])
        if lesf:
            cc *= fpp**max(2.0*(abs(CNp) - self.CN1) + fpp - fprime, 0.0)
        # ---- moment: static table at the lagged effective incidence + unsteady terms
        self.Daf = self.Daf*self.Ecm + (af_deg - self.af_deg)*self.Ecmh
        cpv = p["xcp_v"]*(1.0 - math.cos(math.pi*min(tv/p["Tvl"], 1.0))) if tv > 0.0 else 0.0
        # unsteady moment: the circulatory pitch-rate term is carried by the
        # circulation, so it is reduced by separation in the same proportion
        # as the circulatory normal force (Kirchhoff factor) when
        # cmq_kirchhoff is set
        cm_circ_q = a["CM_circ_q"]*(Kf if p["cmq_kirchhoff"] else 1.0)
        CM = (sm.cm_at(af_deg - self.Daf)
              + self.g_unsteady*(cm_circ_q + a["CM_imp_alpha"] + a["CM_imp_q"]) - cpv*CNv)
        self.fpp_before, self.fpp, self.fp = self.fpp, fpp, fprime
        self.alpha, self.CNpot, self.CNp, self.af_deg = alpha, CNpot, CNp, af_deg
        self.Cv, self.CNv, self.tau, self.lesf = Cv, CNv, tv, lesf
        return dict(alpha=alpha, alpha_dot=alpha_dot, CN=CNfs + CNv, CC=cc, CM=CM, CN_prime=CNp,
                    f_sep=fpp, f_static=fprime, CN_vortex=CNv, tau_v=tv, alpha_f_deg=af_deg,
                    onset_alpha_deg=onset, CN_vortex_removed=removed, reshed=self.reshed)


def outside_compared_range(alpha_mean_deg: float, alpha_amp_deg: float, k: float,
                           static: StaticModel) -> list:
    """Reasons, as text, why a condition lies outside what the measured loops
    cover; empty when it lies inside. Angles in degrees; k [-]."""
    why = []
    if not static.in_range:
        why.append(f"Mach {static.M:.3f} is outside the range of the static data")
    if k > K_MAX_COMPARED + 0.005:
        why.append(f"reduced frequency {k:.3f} is above {K_MAX_COMPARED}")
    peak = abs(alpha_mean_deg) + abs(alpha_amp_deg)
    if peak > PEAK_ALPHA_MAX_COMPARED + 1e-9:
        why.append(f"peak incidence {peak:.1f} deg is above {PEAK_ALPHA_MAX_COMPARED} deg")
    return why


def prescribed_pitch(alpha_mean_deg: float, alpha_amp_deg: float, omega: float, n_per_cycle: int, N: int) -> tuple:
    """Incidence [rad] and its rate [rad/s] at steps 0..N of the motion
    mean + amp sin(phase), phase = 2 pi n / n_per_cycle, omega in rad/s.

    The phase is formed from the step number, so the incidence does not depend
    on chord or speed, and the rate is the analytic derivative. At a turning
    point that falls on a step the rate is exactly zero, not rounding noise;
    there the incidence is neither moving away from zero lift nor returning,
    and the switching rules of the model take their default branch."""
    n = np.arange(N + 1)
    phase = 2.0*np.pi*n/n_per_cycle
    cos = np.cos(phase)
    cos[((4*n) % n_per_cycle == 0) & (((4*n)//n_per_cycle) % 2 == 1)] = 0.0
    return (np.radians(alpha_mean_deg + alpha_amp_deg*np.sin(phase)),
            np.radians(alpha_amp_deg)*omega*cos)


def solve(alpha_mean_deg: float, alpha_amp_deg: float, k: float, M: float, consts: dict | None = None,
          static: StaticModel | None = None, n_per_cycle: int | None = None, n_cycles: int | None = None,
          chord: float = CHORD_M, a_sound: float = A_SOUND, x_pitch: float = 0.25,
          vortex: bool = True, separation: bool = True, strict: bool = False) -> dict:
    """March alpha(t) = mean + amp sin(wt) and return the last cycle.

    Angles in degrees; k = w c / (2U) [-]; M [-]; chord [m]; a_sound [m/s].
    `static` is a StaticModel (built for M if not given). `separation=False`
    gives the attached-flow model; `vortex=False` the model without the
    leading-edge vortex. `strict=True` raises ValueError for a condition
    outside the range the measured loops cover; otherwise the reasons are
    returned in meta["outside_compared_range"]. Returns a dict of arrays over
    the last cycle. `n_per_cycle` and `n_cycles` default to march_resolution(k).
    n_sheddings counts every vortex shed in the cycle: each onset and each
    restart of the vortex clock under sustained onset (n_repeated_sheddings).
    CN_vortex_removed is the vortex normal force set to zero at a step because
    it opposed the separated force (zero where the rule did not act)."""
    auto = march_resolution(k)
    n_per_cycle, n_cycles = n_per_cycle or auto[0], n_cycles or auto[1]
    static = static if static is not None else StaticModel(M)
    outside = outside_compared_range(alpha_mean_deg, alpha_amp_deg, k, static)
    if strict and outside:
        raise ValueError("outside the compared range: " + "; ".join(outside))
    U = M*a_sound
    omega = 2.0*k*U/chord
    N = n_per_cycle*n_cycles
    t = np.linspace(0.0, n_cycles*2.0*np.pi/omega, N + 1)
    alpha, rate = prescribed_pitch(alpha_mean_deg, alpha_amp_deg, omega, n_per_cycle, N)
    model = DynamicStall(M, U, chord, float(t[1] - t[0]), consts, static, x_pitch, vortex, separation)
    recs = [model.start(alpha[0], alpha[1], rate[0])] + [model.step(alpha[n], rate[n]) for n in range(1, N + 1)]
    def col(key: str) -> np.ndarray:
        """One recorded quantity over every step, as an array."""
        return np.array([r[key] for r in recs])
    CN, CC = col("CN"), col("CC")
    CL = CN*np.cos(alpha) + CC*np.sin(alpha)
    CD = CN*np.sin(alpha) - CC*np.cos(alpha) + model.p["CD0"]
    s = slice(N - n_per_cycle, N + 1)
    on = col("onset_alpha_deg")[s]
    on = on[np.isfinite(on)]
    sm = model.sm
    return dict(t=t[s] - t[s][0], alpha_deg=np.degrees(alpha[s]), alpha_dot=rate[s],
                CL=CL[s], CD=CD[s], CM=col("CM")[s], CN=CN[s], CC=CC[s], CN_prime=col("CN_prime")[s],
                f_sep=col("f_sep")[s], f_static=col("f_static")[s], CN_vortex=col("CN_vortex")[s],
                tau_v=col("tau_v")[s], alpha_f_deg=col("alpha_f_deg")[s],
                onset_alpha_deg=float(on[0]) if len(on) else float("nan"),
                n_sheddings=int(len(on)) + int(col("reshed")[s][1:].sum()),
                n_repeated_sheddings=int(col("reshed")[s][1:].sum()),
                CN_vortex_removed=col("CN_vortex_removed")[s],
                meta=dict(M=M, k=k, CN1=sm.CN1, CN_alpha=sm.CN_alpha, alpha0_deg=sm.alpha0_deg,
                          in_static_mach_range=sm.in_range, outside_compared_range=outside,
                          steps_per_cycle=n_per_cycle, cycles=n_cycles,
                          ds=model.ds, constants=model.p))


def frame_runner(consts: dict | None = None, curve: str = "more", **kw: object) -> object:
    """run(fr) for unistall/metrics.py: one measured frame at its own measured
    conditions (M, k, mean angle and amplitude [deg] as stored in the frame)."""
    def run(fr: dict) -> dict:
        """Solve one loaded frame (angles in degrees as stored)."""
        return solve(fr["a0"], fr["da"], fr["k"], fr["M"], consts=consts,
                     static=StaticModel(fr["M"], curve), **kw)
    return run
