# -*- coding: utf-8 -*-
"""
unistall / attached_flow.py
-----------------------
Unsteady ATTACHED-flow loads on a pitching aerofoil by the Leishman-Beddoes
indicial method, written from a citable statement of it and with nothing
tuned:

  [DH19]  Damiani, R. & Hayman, G. (2019), "The Unsteady Aerodynamics Module
          for FAST 8", NREL/TP-5000-66347, doi:10.2172/1576488, section 1.1,
          Eqs. (1.5)-(1.29).

Author: Akosa Samuel Onyejekwe (independent)

What is computed, per time step n (equation numbers are those of [DH19]):

  reduced time        ds = 2 U dt / c                                    (1.5b)
  pitch rate          q  = alpha_dot c / U                               (1.7)
  3/4-chord angle     a34 = alpha + q (3/4 - x_p)       x_p = pitch axis / c
  circulatory         X1, X2 deficiency functions on d(a34)              (1.15)
                      alpha_e = a34 - alpha_0 - X1 - X2                  (1.14)
                      C_N^c = C_Nalpha(M) alpha_e                        (1.13)
  impulsive, alpha    T_I = c / a_s                                      (1.11c)
                      k_alpha = 1 / [(1-M) + pi beta M^2 (A1 b1 + A2 b2)]   see (d) below
                      T_alpha = 0.75 k_alpha T_I                         (1.10)
                      C_N^nc,alpha = (4 T_alpha / M)(K_alpha - K'_alpha) (1.18)
  impulsive, q        T_q = 0.75 k_q T_I                                 (1.10)
                      C_N^nc,q = (T_q / M)(K_q - K'_q)                   (1.19)-(1.20)
  moment about c/4    C_M = C_M0 + C_N^c (1/4 - x_ac)                    (1.22)
                          - (C_Nalpha / 16)(q - K'''_q)                  (1.25)-(1.26)
                          + C_M^nc,alpha + C_M^nc,q                      (1.27)-(1.29)

T_I = c/a_s is the classical time constant.

C_Nalpha(M) here is the measured lift-curve slope at the Mach number of the
case. It already contains compressibility, so it is NOT divided by beta again.

THREE PLACES WHERE THE SOURCE IS NOT SELF-CONSISTENT, and what is done:

  (a) Eq. (1.19) carries a minus sign that Eq. (1.20) does not. The plus sign
      of (1.20) is used: it is the one for which a nose-up pitch acceleration
      adds lift, as apparent mass must.
  (b) Eq. (1.25) multiplies the circulatory pitch-rate moment by c/U although
      q is already non-dimensional. The factor is omitted; with it the term
      would have the dimension of time.
  (c) The text sets b5 = 5; Leishman's value is 0.5. `b5` is an argument; the
      default 0.5 is the published Leishman-Beddoes constant.

  (d) Eqs. (1.11a) and (1.11b) print the same expression for k_alpha and
      k_q, with C_Nalpha M^2 beta in the denominator. They cannot both be
      right, and the alpha form can be derived. The exact short-time response
      to a step in incidence (Lomax et al. 1952, valid for 0 <= s <= 2M/(1+M))
      is C_N = (4/M)[1 - (1-M) s / (2M)], so its slope at s = 0 is
      -2(1-M)/M^2. The model's response is
          (4/M) exp(-s/T') + (2 pi/beta)[1 - A1 exp(-b1 beta^2 s) - A2 exp(-b2 beta^2 s)],
      whose slope at s = 0 is -4/(M T') + 2 pi beta (A1 b1 + A2 b2). Equating,
          T' = 2M / [(1-M) + pi beta M^2 (A1 b1 + A2 b2)],
      and since T' = 2M k_alpha in semichords,
          k_alpha = 1 / [(1-M) + pi beta M^2 (A1 b1 + A2 b2)].
      That is the form of Leishman & Beddoes (1989) -- beta to the first
      power -- and it is what is implemented. The pitch-rate constant carries
      2 pi beta M^2 in the same place. The factor 0.75 applied to both in
      (1.10) is Leishman's allowance for effects outside piston theory.

TWO FORMS OF THE IMPULSIVE ALPHA MOMENT are available, because [DH19] gives a
simplification of the original:
    "two_term"   -(1/M)[A3 exp(-t/(b3 K_aM T_I)) + A4 exp(-t/(b4 K_aM T_I))],
                 A3 = 1.5, A4 = -0.5, b3 = 0.25, b4 = 0.1,
                 K_aM = (A3 b4 + A4 b3) / (b3 b4 (1 - M))   -- Leishman & Beddoes
    "quarter"    -C_N^nc,alpha / 4                           -- [DH19] Eq. (1.27)
Which one is used is decided by the attached-flow checks, not by preference.

No low-pass filter is applied to alpha, q or their differences ([DH19] Eq. 1.8):
the motions solved here are prescribed and smooth, so there is no noise to
remove, and a filter would add a phase lag the checks would then measure.
"""
import numpy as np

# published Leishman-Beddoes indicial constants
A1, A2, B1, B2 = 0.30, 0.70, 0.14, 0.53
A3, A4, B3, B4 = 1.5, -0.5, 0.25, 0.1
A5, B5_DEFAULT = 1.0, 0.5

# R. T. Jones' two-term approximation to the Wagner function (incompressible)
JONES_A1, JONES_A2, JONES_B1, JONES_B2 = 0.165, 0.335, 0.0455, 0.3

# LOW-MACH BLEND. The compressible indicial constants above were fitted for
# M >= 0.3 and do not tend to the incompressible limit: at M = 0.05 they are
# 6 % and 6.5 deg from Theodorsen's solution in lift and 13 % and 11 deg in
# moment. Incompressible thin-aerofoil theory (Wagner's function in Jones'
# approximation, with the apparent-mass loads) is exact in that limit. The
# model therefore uses the incompressible form at and below M_INCOMPRESSIBLE,
# the compressible form at and above M_COMPRESSIBLE, and a linear blend of the
# two sets of loads in between. Above M_COMPRESSIBLE nothing is changed.
M_INCOMPRESSIBLE, M_COMPRESSIBLE = 0.10, 0.20


def compressible_weight(M: float) -> float:
    """Weight of the compressible loads at Mach number M [-]: 0 at and below
    M_INCOMPRESSIBLE, 1 at and above M_COMPRESSIBLE, linear in between."""
    return float(np.clip((M - M_INCOMPRESSIBLE)/(M_COMPRESSIBLE - M_INCOMPRESSIBLE), 0.0, 1.0))


class _Incompressible:
    """Incompressible attached-flow loads, one step at a time: Jones' Wagner
    approximation for the circulatory part and the apparent-mass terms of
    thin-aerofoil theory for the rest. Same interface as AttachedFlow."""

    def __init__(self, U: float, c: float, CN_alpha: float, dt: float, alpha0: float, x_pitch: float) -> None:
        self.U, self.c, self.CN_alpha, self.dt, self.alpha0, self.x_pitch = U, c, CN_alpha, dt, alpha0, x_pitch
        ds = 2.0*U*dt/c
        self.E1, self.E1h = np.exp(-JONES_B1*ds), np.exp(-JONES_B1*ds/2)
        self.E2, self.E2h = np.exp(-JONES_B2*ds), np.exp(-JONES_B2*ds/2)
        self.X1 = self.X2 = 0.0
        self.a = 2.0*x_pitch - 1.0                       # pitch axis aft of mid-chord, semichords

    def _loads(self, alpha_e: float, q: float, Q2: float) -> dict:
        a = self.a
        CNc = self.CN_alpha*alpha_e
        CNnc = np.pi*(q/2.0 - a*Q2)
        CN = CNc + CNnc
        CMa = (np.pi/2.0)*(-(0.5 - a)*q/2.0 - (0.125 + a*a)*Q2) + (a + 0.5)*CNc/2.0
        return dict(CN=CN, CM=CMa - CN*(self.x_pitch - 0.25), CN_circ=CNc, CN_imp_alpha=CNnc, CN_imp_q=0.0,
                    CM_circ_q=0.0, CM_imp_alpha=CMa - CN*(self.x_pitch - 0.25), CM_imp_q=0.0,
                    alpha_e=alpha_e, q=q)

    def start(self, alpha_first: float, alpha_second: float) -> dict:
        """Loads at the first instant, from the first two incidences [rad]."""
        Ka = (alpha_second - alpha_first)/self.dt
        self.alpha, self.Ka = alpha_first, Ka
        q = Ka*self.c/self.U
        self.a34 = alpha_first + q*(0.75 - self.x_pitch)
        return self._loads(self.a34 - self.alpha0, q, 0.0)

    def step(self, alpha: float) -> dict:
        """Advance to the next instant, at incidence alpha [rad]."""
        Ka = (alpha - self.alpha)/self.dt
        q = Ka*self.c/self.U
        Q2 = (self.c/(2.0*self.U))**2*(Ka - self.Ka)/self.dt
        a34 = alpha + q*(0.75 - self.x_pitch)
        da = a34 - self.a34
        self.X1 = self.X1*self.E1 + JONES_A1*da*self.E1h
        self.X2 = self.X2*self.E2 + JONES_A2*da*self.E2h
        self.alpha, self.Ka, self.a34 = alpha, Ka, a34
        return self._loads(a34 - self.alpha0 - self.X1 - self.X2, q, Q2)


class AttachedFlow:
    """The attached-flow model advanced one time step at a time.

    U [m/s], c [m], M [-], CN_alpha [1/rad] at this Mach number, dt [s];
    x_ac and x_pitch as fractions of chord; CM0 [-]; alpha0 [rad].
    Call `start` once with the first two incidences, then `step` with each
    later incidence. Every call returns the loads at that step as a dict of
    floats: CN, CM, CN_circ, CN_imp_alpha, CN_imp_q, CM_circ_q, CM_imp_alpha,
    CM_imp_q, alpha_e [rad], q [-].
    """

    def __init__(self, U: float, c: float, M: float, CN_alpha: float, dt: float,
                 x_ac: float = 0.25, CM0: float = 0.0, alpha0: float = 0.0,
                 x_pitch: float = 0.25, b5: float = B5_DEFAULT,
                 moment_alpha_form: str = "two_term", low_mach_blend: bool = True) -> None:
        if moment_alpha_form not in ("two_term", "quarter"):
            raise ValueError(moment_alpha_form)
        self.w = compressible_weight(M) if low_mach_blend else 1.0
        self.inc = _Incompressible(U, c, CN_alpha, dt, alpha0, x_pitch) if self.w < 1.0 else None
        self.U, self.c, self.M, self.CN_alpha, self.dt = U, c, M, CN_alpha, dt
        self.x_ac, self.CM0, self.alpha0, self.x_pitch = x_ac, CM0, alpha0, x_pitch
        self.form = moment_alpha_form
        ds = 2.0*U*dt/c
        beta2 = 1.0 - M*M
        beta = np.sqrt(beta2)
        TI = c/(U/M)
        k_a = 1.0/((1.0 - M) + np.pi*beta*M*M*(A1*B1 + A2*B2))
        k_q = 1.0/((1.0 - M) + 2.0*np.pi*beta*M*M*(A1*B1 + A2*B2))
        self.T_a, self.T_q = 0.75*k_a*TI, 0.75*k_q*TI
        k_mq = 7.0/(15.0*(1.0 - M) + 1.5*CN_alpha*A5*b5*beta*M*M)
        self.T_mq = k_mq*k_mq*TI
        K_aM = (A3*B4 + A4*B3)/(B3*B4*(1.0 - M))
        self.T3, self.T4 = B3*K_aM*TI, B4*K_aM*TI
        def half(T: float) -> tuple:
            """Decay factors over a step and half a step for time constant T [s]."""
            return (np.exp(-dt/T), np.exp(-dt/(2*T)))
        self.E1, self.E1h = np.exp(-B1*beta2*ds), np.exp(-B1*beta2*ds/2)
        self.E2, self.E2h = np.exp(-B2*beta2*ds), np.exp(-B2*beta2*ds/2)
        self.E5, self.E5h = np.exp(-b5*beta2*ds), np.exp(-b5*beta2*ds/2)
        (self.Ea, self.Eah), (self.Eq, self.Eqh) = half(self.T_a), half(self.T_q)
        (self.Em, self.Emh) = half(self.T_mq)
        (self.E3, self.E3h), (self.E4, self.E4h) = half(self.T3), half(self.T4)
        self.X1 = self.X2 = self.X5 = self.Dka = self.Dkq = self.Dmq = self.D3 = self.D4 = 0.0
        self.meta = dict(T_I=TI, k_alpha=k_a, k_q=k_q, T_alpha=self.T_a, T_q=self.T_q,
                         k_mq=k_mq, K_alphaM=K_aM, ds=ds, dt=dt)

    def _blend(self, lb: dict, inc: dict) -> dict:
        """Compressible and incompressible loads combined with weight w."""
        if self.inc is None:
            return lb
        out = {k: self.w*lb[k] + (1.0 - self.w)*inc[k] for k in lb}
        # the quarter-chord static offset is added once, to the blended load
        return out

    def start(self, alpha_first: float, alpha_second: float) -> dict:
        """Set the state at the first instant. The first difference is taken
        from the first two incidences [rad], so no start-up spike enters."""
        return self._blend(self._start(alpha_first, alpha_second),
                           self.inc.start(alpha_first, alpha_second) if self.inc else {})

    def step(self, alpha: float) -> dict:
        """Advance to the next instant, at incidence alpha [rad]."""
        return self._blend(self._step(alpha), self.inc.step(alpha) if self.inc else {})

    def _start(self, alpha_first: float, alpha_second: float) -> dict:
        Ka = (alpha_second - alpha_first)/self.dt
        q = Ka*self.c/self.U
        self.alpha, self.Ka, self.q, self.Kq = alpha_first, Ka, q, 0.0
        self.a34 = alpha_first + q*(0.75 - self.x_pitch)
        alpha_e = self.a34 - self.alpha0
        CNc = self.CN_alpha*alpha_e
        return dict(CN=CNc, CM=self.CM0 + CNc*(0.25 - self.x_ac), CN_circ=CNc, CN_imp_alpha=0.0,
                    CN_imp_q=0.0, CM_circ_q=0.0, CM_imp_alpha=0.0, CM_imp_q=0.0, alpha_e=alpha_e, q=q)

    def _step(self, alpha: float) -> dict:
        Ka = (alpha - self.alpha)/self.dt
        q = Ka*self.c/self.U
        Kq = (q - self.q)/self.dt
        a34 = alpha + q*(0.75 - self.x_pitch)
        da = a34 - self.a34
        self.X1 = self.X1*self.E1 + A1*da*self.E1h
        self.X2 = self.X2*self.E2 + A2*da*self.E2h
        alpha_e = a34 - self.alpha0 - self.X1 - self.X2
        CNc = self.CN_alpha*alpha_e
        dKa, dKq = Ka - self.Ka, Kq - self.Kq
        self.Dka = self.Dka*self.Ea + dKa*self.Eah
        CNnc_a = (4.0*self.T_a/self.M)*(Ka - self.Dka)
        self.Dkq = self.Dkq*self.Eq + dKq*self.Eqh
        CNnc_q = (self.T_q/self.M)*(Kq - self.Dkq)
        self.X5 = self.X5*self.E5 + A5*(q - self.q)*self.E5h
        CMc_q = -(self.CN_alpha/16.0)*(q - self.X5)
        if self.form == "two_term":
            self.D3 = self.D3*self.E3 + dKa*self.E3h
            self.D4 = self.D4*self.E4 + dKa*self.E4h
            CMnc_a = -(A3*self.T3/self.M)*(Ka - self.D3) - (A4*self.T4/self.M)*(Ka - self.D4)
        else:
            CMnc_a = -CNnc_a/4.0
        self.Dmq = self.Dmq*self.Em + dKq*self.Emh
        CMnc_q = -(7.0*self.T_mq/(12.0*self.M))*(Kq - self.Dmq)
        CM = self.CM0 + CNc*(0.25 - self.x_ac) + CMc_q + CMnc_a + CMnc_q
        self.alpha, self.Ka, self.q, self.Kq, self.a34 = alpha, Ka, q, Kq, a34
        return dict(CN=CNc + CNnc_a + CNnc_q, CM=CM, CN_circ=CNc, CN_imp_alpha=CNnc_a,
                    CN_imp_q=CNnc_q, CM_circ_q=CMc_q, CM_imp_alpha=CMnc_a, CM_imp_q=CMnc_q,
                    alpha_e=alpha_e, q=q)


def solve_attached(alpha: np.ndarray, t: np.ndarray, U: float, c: float, M: float,
                   CN_alpha: float, x_ac: float = 0.25, CM0: float = 0.0, alpha0: float = 0.0,
                   x_pitch: float = 0.25, b5: float = B5_DEFAULT,
                   moment_alpha_form: str = "two_term", low_mach_blend: bool = True) -> dict:
    """Attached-flow C_N and C_M(c/4) for a prescribed incidence history.

    alpha [rad] at uniformly spaced times t [s]; U [m/s]; c [m]; M [-];
    CN_alpha [1/rad] at this Mach number; x_ac, x_pitch as fractions of chord;
    CM0 [-]; alpha0 [rad]. Returns a dict of arrays the length of `alpha`.
    """
    alpha = np.asarray(alpha, float)
    t = np.asarray(t, float)
    model = AttachedFlow(U, c, M, CN_alpha, float(t[1] - t[0]), x_ac, CM0, alpha0, x_pitch, b5,
                         moment_alpha_form, low_mach_blend)
    recs = [model.start(alpha[0], alpha[1])] + [model.step(a) for a in alpha[1:]]
    out = {k: np.array([r[k] for r in recs]) for k in recs[0]}
    out["meta"] = model.meta
    return out


# --------------------------------------------------------------------------- #
#  The incompressible limit, by the same marching scheme
# --------------------------------------------------------------------------- #
#  R. T. Jones' two-term approximation to the Wagner function,
#      phi(s) = 1 - 0.165 exp(-0.0455 s) - 0.335 exp(-0.3 s),
#  with the apparent-mass loads of incompressible thin-aerofoil theory. This is
#  Theodorsen's theory in the time domain. It exists here for one purpose: it
#  runs through the SAME deficiency-function recurrences as solve_attached, so
#  its agreement with Theodorsen's closed form verifies the marching scheme
#  independently of any compressible constant.

def solve_attached_incompressible(alpha: np.ndarray, t: np.ndarray, U: float, c: float,
                                  CN_alpha: float = 2.0*np.pi, x_pitch: float = 0.25) -> dict:
    """Incompressible attached-flow C_N and C_M(c/4) for prescribed pitch.

    alpha [rad] at uniformly spaced times t [s]; U [m/s]; c [m]; CN_alpha
    [1/rad]; x_pitch as a fraction of chord."""
    alpha = np.asarray(alpha, float)
    t = np.asarray(t, float)
    N = len(alpha)
    dt = float(t[1] - t[0])
    ds = 2.0*U*dt/c
    ad = np.gradient(alpha, dt)
    add = np.gradient(ad, dt)
    q = ad*c/U
    Q2 = (c/(2.0*U))**2*add
    a = 2.0*x_pitch - 1.0
    a34 = alpha + q*(0.75 - x_pitch)
    E1, E1h = np.exp(-JONES_B1*ds), np.exp(-JONES_B1*ds/2)
    E2, E2h = np.exp(-JONES_B2*ds), np.exp(-JONES_B2*ds/2)
    X1 = X2 = 0.0
    CNc = np.zeros(N)
    CNc[0] = CN_alpha*a34[0]
    for n in range(1, N):
        da = a34[n] - a34[n-1]
        X1 = X1*E1 + JONES_A1*da*E1h
        X2 = X2*E2 + JONES_A2*da*E2h
        CNc[n] = CN_alpha*(a34[n] - X1 - X2)
    CNnc = np.pi*(q/2.0 + a*(-Q2))                        # pi b (U a_dot - b a a_ddot)/U^2
    CN = CNc + CNnc
    CMa = (np.pi/2.0)*(-(0.5 - a)*q/2.0 - (0.125 + a*a)*Q2) + (a + 0.5)*CNc/2.0
    CM = CMa - CN*(x_pitch - 0.25)
    return dict(CN=CN, CM=CM, CN_circ=CNc, CN_imp=CNnc, q=q)


# --------------------------------------------------------------------------- #
#  Theodorsen's solution for harmonic pitch, the reference the model is held to
# --------------------------------------------------------------------------- #
def theodorsen_pitch(k: float, x_pitch: float = 0.25) -> tuple:
    """Complex C_L/alpha and C_M(c/4)/alpha for alpha = alpha_bar exp(i w t),
    incompressible, pitch axis at x_pitch (fraction of chord). k = w c / (2U).

    With b = c/2 and a = 2 x_pitch - 1 (axis aft of mid-chord in semichords):
      L   = pi rho b^2 (U a_dot - b a a_ddot)
            + 2 pi rho U b C(k) [U alpha + b (1/2 - a) a_dot]
      M_a = pi rho b^2 [-U b (1/2 - a) a_dot - b^2 (1/8 + a^2) a_ddot]
            + 2 pi rho U b^2 (a + 1/2) C(k) [U alpha + b (1/2 - a) a_dot]
    and the moment is transferred from the axis to the quarter chord.
    """
    from scipy.special import hankel2
    Ck = hankel2(1, k)/(hankel2(1, k) + 1j*hankel2(0, k))
    a = 2.0*x_pitch - 1.0
    ik = 1j*k
    w34 = 1.0 + (0.5 - a)*ik                              # [U alpha + b(1/2-a) a_dot]/(U alpha)
    CL = np.pi*(ik + a*k*k) + 2.0*np.pi*Ck*w34            # L / (rho U^2 b)
    CMa = (np.pi/2.0)*(-(0.5 - a)*ik + (0.125 + a*a)*k*k) + np.pi*(a + 0.5)*Ck*w34   # M_a/(2 rho U^2 b^2)
    # transfer to c/4: M_c/4 = M_a - L * (x_pitch - 1/4) c, nose-up positive
    CM = CMa - CL*(x_pitch - 0.25)
    return CL, CM, Ck


def first_harmonic(y: np.ndarray, t: np.ndarray, omega: float) -> complex:
    """Complex amplitude Y of the exp(i w t) component of y over whole cycles,
    with y = Re{Y exp(i w t)}."""
    trapz = getattr(np, "trapezoid", None) or np.trapz
    return 2.0*trapz(y*np.exp(-1j*omega*t), t)/(t[-1] - t[0])
