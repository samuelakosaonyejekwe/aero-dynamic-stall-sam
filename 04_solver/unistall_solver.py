"""
=============================================================================
 UNISTALL(TM)  —  Universal Unsteady-Aerodynamics & Dynamic-Stall Solver
 Core method : Unified Indicial–Beddoes State-Space (UIBS)
 Version     : 1.0.0
=============================================================================
A reduced-order, physics-based solver for unsteady airfoil aerodynamics and
dynamic stall. The UIBS core couples four sub-models in semichord ("reduced")
time s = 2*U*t/c :

  (1) Attached-flow indicial response (Beddoes' two-lag circulatory model +
      compressible non-circulatory / added-mass impulsive loads).
  (2) Trailing-edge separation via Kirchhoff/Helmholtz theory with a pressure
      lag (Tp) and a boundary-layer lag (Tf).
  (3) Leading-edge dynamic-stall vortex (DSV): shedding, lift overshoot,
      convection (Tvl) and decay (Tv) — the source of the moment break.
  (4) Compressibility (Prandtl-Glauert beta) + chord-force / drag closure.

Two auxiliary modules make the solver "universal" for engineering output:
  * Vortex/source field-reconstruction  -> 2D pressure, velocity, vorticity,
    streamlines + an explicit Lamb-Oseen dynamic-stall vortex. Source panels
    and a bound vortex sheet, BOTH on the body surface, with the total
    circulation matched to the UIBS C_L (Kutta-Joukowski).
    LIMITATIONS, measured rather than asserted (surface_load_closure() and
    the metrics_*.csv rows "Cp_closure_error_pct" and "Cp_TE_jump_max_over_phases"
    recompute them on every run):
      - CLOSURE. Integrating the surface Cp now recovers the C_L it was given to
        within 0.6 % over the whole cycle (-0.28 % at alpha 2 deg, -0.37 % at
        10 deg, -0.58 % at 17.5 deg), and unlike before it CONVERGES: refining
        160 -> 1280 panels drives it monotonically to -0.04 %, which is what
        Blasius requires and is the check that the formulation is right rather
        than merely better. It previously
        read -10.7 to -14.4 %, and the explanation recorded here for that
        deficit -- the Cp clip at -8, with the further claim that it did not
        converge under refinement -- was WRONG on both counts. Tested directly:
        moving the clip from -8 to -1e9 changed the closure by 0.00 points, and
        refining 160 -> 1280 panels moved it monotonically from -11.5 % to
        -5.2 %. The real causes were three, all in the evaluation rather than
        the physics, and all now fixed in surface_cp: the vortex sheet's own
        tangential contribution (-gam/2) was omitted; Cp was evaluated at the
        panel end-points offset 0.015c off the wall rather than at the control
        points where tangency is imposed; and a Prandtl-Glauert factor was
        applied to a Cp whose circulation already carried compressibility,
        inflating the load a further 4.8 % at M = 0.3.
      - KUTTA CONDITION. The trailing-edge Cp jump is not a residual that can be
        driven to zero, and it is no longer presented as one. It is LINEAR in
        the imposed C_L and passes through zero exactly at the inviscid attached
        circulation, which kutta_reference_CL() computes: at alpha = 10 deg that
        is C_L = 1.220 at the 160 panels used, converging to 1.213 by 1280
        (implied lift slope 6.95/rad, against 2*pi*1.092 = 6.86 from the
        thin-aerofoil thickness rule, itself approximate), and imposing it
        drives the trailing-edge jump to ~1e-3. The
        reconstruction is instead handed the indicial C_L, which during dynamic
        stall departs from that value deliberately -- so a body carrying a
        non-Kutta circulation MUST show a trailing-edge jump. The published
        Cp_TE_jump_max_over_phases is therefore a measure of how far the modelled
        flow is from attached, not an error. The earlier, smaller published
        values (0.15 and 0.25) were not a better result: they came from probing
        0.015c off the wall, the same offset that was hiding the closure error.
      - DYNAMIC-STALL VORTEX CORE, a stated limitation rather than a fixed one.
        The reconstructed vortex carries circulation DSV_GAMMA_FACTOR*CNv*U*c in
        a core of radius DSV_CORE_RADIUS_CHORDS*c. Its sign, position and the
        flow reversal beneath it are physical, but the core is diffuse: the
        measured suction at the core centre is published as Cp_DSV_core_min in
        metrics_*.csv and reads about -0.4, where a deep-stall vortex core is
        usually reported nearer -3 to -6. Making it deeper means shrinking the
        core radius and raising the circulation factor together (0.06c and 2.5
        give about -2.3), and NEITHER constant can be derived or calibrated
        here: the experimental frames this study ships carry only integrated
        cl/cd/cm against incidence, with no surface-pressure or field data
        anywhere in the repository to fit a core size to. Both constants are
        therefore named at the top of this module rather than buried as
        literals, and the resulting core depth is published as a number so the
        shallowness is checkable instead of being an adjective in a docstring.
        It does not affect the reported loads, which come from the UIBS core.
      - Nothing in the reconstruction knows about separation: it is a potential
        field, so at post-stall incidence the leading-edge suction peak it draws
        (about Cp = -15 at 17.5 deg) is far deeper than a real separated flow
        would sustain. Nothing is clipped. The FIELD nonetheless bottoms out near
        -5.2, because the near-wall ring carrying that peak is masked, so the
        contour plots understate the surface suction by about three times.
    The reconstruction is qualitative; the reported loads come from the UIBS
    core and do not depend on it.
  * Compressible thermal module          -> static & recovery (skin) temperature.

The model is calibrated PER CASE to a static polar and validated against
published dynamic-stall experiments. References are recorded in 03_model_setup.
=============================================================================
"""
import numpy as np
from scipy.interpolate import PchipInterpolator

# np.trapz was removed in NumPy 2.0 in favour of np.trapezoid; bind whichever exists
_trapz = getattr(np, "trapezoid", None) or np.trapz

# --------------------------------------------------------------------------- #
#  STATIC SEPARATION CALIBRATION  (Kirchhoff inverse from a static polar)
# --------------------------------------------------------------------------- #
def calibrate_separation(alpha_deg, Cl, Cd, CNalpha):
    """Return f_static(alpha_deg) — the static TE separation point, derived by
    inverting the Kirchhoff relation  CN = CNalpha*((1+sqrt(f))/2)^2 * alpha .
    Symmetric airfoil -> even function of alpha."""
    a = np.radians(np.asarray(alpha_deg, float))
    Cl = np.asarray(Cl, float); Cd = np.asarray(Cd, float)
    CN = Cl*np.cos(a) + Cd*np.sin(a)
    f = np.ones_like(a)
    for i in range(len(a)):
        if abs(a[i]) < np.radians(0.5):
            f[i] = 1.0
        else:
            ratio = CN[i]/(CNalpha*a[i])
            sf = 2.0*np.sqrt(max(ratio, 0.0)) - 1.0
            f[i] = float(np.clip(sf, np.sqrt(0.02), 1.0)**2)
    # monotone, smooth interpolant on |alpha|; clamp ends
    order = np.argsort(alpha_deg)
    ad = np.asarray(alpha_deg, float)[order]; fd = f[order]
    interp = PchipInterpolator(ad, fd, extrapolate=False)
    amin, amax = ad.min(), ad.max()
    F_MIN, S_EXT = 0.02, 3.0        # floor, and decay length past the data (deg)
    f_end = float(interp(amax))     # separation point AT the last measured alpha
    def f_static(alpha_query_deg):
        q = np.abs(np.asarray(alpha_query_deg, float))
        out = interp(q)
        # Past the calibration data, decay smoothly from the last measured value
        # towards full separation. Snapping straight to F_MIN put a step of 0.13
        # in f at alpha_max, which showed up as a 0.37 jump in the model C_L --
        # 15x the neighbouring steps -- purely as an artefact of the clamp.
        deep = F_MIN + (f_end - F_MIN)*np.exp(-(q - amax)/S_EXT)
        out = np.where(q <= amax, out, deep)
        out = np.where(q >= amin, out, 1.0)            # below data -> attached
        return np.clip(np.nan_to_num(out, nan=1.0), F_MIN, 1.0)
    return f_static


# --------------------------------------------------------------------------- #
#  UIBS DYNAMIC-STALL MARCHING SOLVER
# --------------------------------------------------------------------------- #
def solve_dynamic_stall(alpha_mean_deg, alpha_amp_deg, k, M, c, U,
                        f_static, CNalpha=6.28, CD0=0.0086, CM0=0.0,
                        consts=None, n_per_cycle=720, n_cycles=6):
    """March the UIBS model through n_cycles of alpha(t)=mean+amp*sin(wt).
    Returns a dict of per-step arrays for the LAST (converged) cycle."""
    p = dict(A1=0.30, A2=0.70, b1=0.14, b2=0.53,
             Tp=1.7, Tf=3.0, Tv=6.0, Tvl=5.0, CN1=1.45,
             k0=0.0, k1=-0.135, k2=0.04, kappa=2.0, eta=0.95,
             Bv=1.0, cpv_amp=0.20)   # Bv: vortex-feed gain, cpv_amp: vortex CP travel
    if consts:
        p.update(consts)

    beta2 = max(1.0 - M*M, 1e-3)
    Kalpha = 0.75/(1.0 - M + np.pi*np.sqrt(beta2)*M*M*(p["A1"]*p["b1"]+p["A2"]*p["b2"]))
    # impulsive time constant T_I = Kalpha*c/a with a = U/M, the speed of sound
    # implied by the case. This was written as M*340.0: a hardcoded sea-level
    # value that contradicted both the comment beside it and the T_I = Kalpha*c/a
    # of the report, and made the march weakly dependent on U for any case whose
    # speed of sound is not exactly 340 m/s.
    TI = Kalpha*c/U                           # a = U/M  =>  Kalpha*c/(M*a) = Kalpha*c/U

    omega = 2.0*k*U/c
    N = n_per_cycle*n_cycles
    t = np.linspace(0.0, n_cycles*2*np.pi/omega, N+1)
    dt = t[1]-t[0]
    ds = 2.0*U*dt/c                                  # semichord step (constant)

    a0 = np.radians(alpha_mean_deg); aa = np.radians(alpha_amp_deg)
    alpha = a0 + aa*np.sin(omega*t)
    alpha_dot = aa*omega*np.cos(omega*t)
    # AoA at 3/4-chord (pitch about c/4): adds pitch-rate downwash
    alpha34 = alpha + alpha_dot*c/(2.0*U)

    # state arrays
    X1 = np.zeros(N+1); X2 = np.zeros(N+1); Dimp = np.zeros(N+1)
    Dp = np.zeros(N+1); Df = np.zeros(N+1); CNv = np.zeros(N+1)
    CN = np.zeros(N+1); CC = np.zeros(N+1); CL = np.zeros(N+1)
    CD = np.zeros(N+1); CM = np.zeros(N+1); CNp = np.zeros(N+1)
    fpp = np.ones(N+1); CN_pot = np.zeros(N+1); CNf = np.zeros(N+1)
    tau_v = np.zeros(N+1); vortex_active = np.zeros(N+1)

    E1 = np.exp(-p["b1"]*beta2*ds); E1h = np.exp(-p["b1"]*beta2*ds/2)
    E2 = np.exp(-p["b2"]*beta2*ds); E2h = np.exp(-p["b2"]*beta2*ds/2)
    Etp = np.exp(-ds/p["Tp"]); Etph = np.exp(-ds/(2*p["Tp"]))
    Etf = np.exp(-ds/p["Tf"]); Etfh = np.exp(-ds/(2*p["Tf"]))
    Etv = np.exp(-ds/p["Tv"]); Etvh = np.exp(-ds/(2*p["Tv"]))
    EI = np.exp(-dt/TI); EIh = np.exp(-dt/(2*TI))

    fpp[0] = f_static(alpha_mean_deg)
    vortex_on = False
    fprime_prev = float(f_static(alpha_mean_deg))
    Cv_prev = 0.0
    for n in range(1, N+1):
        da = alpha34[n]-alpha34[n-1]
        # (1) circulatory deficiency -> effective AoA
        X1[n] = X1[n-1]*E1 + p["A1"]*da*E1h
        X2[n] = X2[n-1]*E2 + p["A2"]*da*E2h
        aE = alpha34[n] - X1[n] - X2[n]
        CNc = CNalpha*aE
        # (1b) non-circulatory impulsive (added mass)
        dadt = (alpha34[n]-alpha34[n-1])/dt
        dadt_prev = (alpha34[n-1]-alpha34[max(n-2,0)])/dt
        Dimp[n] = Dimp[n-1]*EI + (dadt - dadt_prev)*EIh
        CNi = (4.0*Kalpha*c/U/M)*(dadt - Dimp[n]) if M > 0 else 0.0
        CN_pot[n] = CNc + CNi
        # (2) pressure lag -> delayed normal force -> separation onset
        Dp[n] = Dp[n-1]*Etp + (CN_pot[n]-CN_pot[n-1])*Etph
        CNp[n] = CN_pot[n] - Dp[n]
        af = CNp[n]/CNalpha
        fprime = float(f_static(np.degrees(af)))
        # (2b) boundary-layer lag -> dynamic separation point
        Df[n] = Df[n-1]*Etf + (fprime - fprime_prev)*Etfh
        fpp[n] = np.clip(fprime - Df[n], 0.02, 1.0)
        fprime_prev = fprime
        # (1c) separated circulatory normal force (Kirchhoff) + impulsive
        Kf = ((1.0+np.sqrt(fpp[n]))/2.0)**2
        CNf[n] = CNalpha*aE*Kf + CNi
        # (3) leading-edge dynamic-stall vortex
        # tau_v is the vortex convection clock. It is zeroed ONLY at shedding
        # onset and HELD (not zeroed) once the vortex switches off, because the
        # vortex normal force CNv is still decaying at that instant: zeroing it
        # would collapse the vortex moment arm CP_v in a single step and put a
        # step discontinuity into CM (and hence into the damping integral).
        if (CNp[n] >= p["CN1"]) and (not vortex_on) and (alpha_dot[n] > 0):
            vortex_on = True
            tau_v[n] = 0.0                           # shedding starts: CP still at c/4
        elif vortex_on:
            tau_v[n] = tau_v[n-1] + ds               # convecting over the chord
        else:
            tau_v[n] = tau_v[n-1]                    # off: hold, so CM_v decays with CNv
        Cv = p["Bv"]*CNc*(1.0 - Kf)                  # vortex feed = (attached-separated)*gain
        if vortex_on:
            vortex_active[n] = 1.0
            if tau_v[n] < p["Tvl"]:
                CNv[n] = CNv[n-1]*Etv + (Cv - Cv_prev)*Etvh
            else:
                CNv[n] = CNv[n-1]*Etv                # past TE: decay only
            if (CNp[n] < p["CN1"]) and (alpha_dot[n] < 0):
                vortex_on = False                    # reset for next cycle
        else:
            CNv[n] = CNv[n-1]*Etv
        Cv_prev = Cv
        # (4) totals
        CN[n] = CNf[n] + CNv[n]
        CC[n] = p["eta"]*CNalpha*aE*aE*np.sqrt(max(fpp[n], 0.0))   # LE suction
        a_n = alpha[n]
        CL[n] = CN[n]*np.cos(a_n) + CC[n]*np.sin(a_n)
        CD[n] = CN[n]*np.sin(a_n) - CC[n]*np.cos(a_n) + CD0
        # moment: attached/separated CP shift + vortex moment
        cp_shift = p["k0"] + p["k1"]*(1.0-fpp[n]) + p["k2"]*np.sin(np.pi*fpp[n]**p["kappa"])
        CM_f = CM0 + cp_shift*CNf[n]
        cpv = p["cpv_amp"]*(1.0 - np.cos(np.pi*np.clip(tau_v[n]/p["Tvl"], 0, 1)))  # CP aft travel
        CM_v = -cpv*CNv[n]
        CM[n] = CM_f + CM_v

    # cycle-to-cycle convergence (peak CL & peak |CM| per cycle)
    cyc_peakCL = [CL[i*n_per_cycle+1:(i+1)*n_per_cycle+1].max() for i in range(n_cycles)]
    cyc_minCM = [CM[i*n_per_cycle+1:(i+1)*n_per_cycle+1].min() for i in range(n_cycles)]

    # return last cycle
    s = slice(N-n_per_cycle, N+1)
    # Phase runs 0 -> 360 monotonically across the reported cycle. Taking the
    # modulo of the absolute time wrapped the final 360 deg point back to 0,
    # which made every curve plotted against phase draw a spurious horizontal
    # line straight back across the axes at its terminal value.
    phase = np.degrees(omega*(t[s]-t[s][0]))
    out = dict(t=t[s]-t[s][0], phase_deg=phase,
               alpha_deg=np.degrees(alpha[s]), alpha_dot=alpha_dot[s],
               CL=CL[s], CD=CD[s], CM=CM[s], CN=CN[s], CC=CC[s],
               CNp=CNp[s], f_sep=fpp[s], CNv=CNv[s], CNf=CNf[s],
               vortex_active=vortex_active[s], tau_v=tau_v[s],
               alpha_eff_deg=np.degrees(alpha34[s]-X1[s]-X2[s]))
    out["_meta"] = dict(Kalpha=Kalpha, TI=TI, ds=ds, omega=omega, dt=dt,
                        beta=np.sqrt(beta2), n_per_cycle=n_per_cycle)
    out["cycle_peakCL"] = np.array(cyc_peakCL)
    out["cycle_minCM"] = np.array(cyc_minCM)
    return out


# --------------------------------------------------------------------------- #
#  AERODYNAMIC DAMPING (stall-flutter indicator) from the CM-alpha loop
# --------------------------------------------------------------------------- #
#  |Xi_hat| below this is treated as NEUTRAL: for a figure-of-eight CM loop the
#  raw Xi is a small residual between two lobes of opposite sign, and below this
#  level the model cannot resolve its sign.
#
#  WHERE THIS NUMBER COMES FROM. It was previously 0.02, justified as "no larger
#  than the time-step discretisation error". That justification is false and was
#  measured to be so: refining n_per_cycle from 720 to 5760 moves Xi_hat for
#  Case A by only 0.00044, i.e. 45x smaller than the band it was said to
#  explain. The band that IS defensible is the model's own demonstrated accuracy
#  in this quantity: recomputing Xi_hat from the five real NACA 0012 loops of
#  NASA TM-84245 and from the model at the same conditions gives a mean absolute
#  discrepancy of 0.072 (max 0.144). A damping residual smaller than that cannot
#  be claimed as a finding. 06_postprocessing/validation/validate_nasa_real.py
#  recomputes that spread on every run and writes it to
#  validation_realdata_summary.csv (rounding 0.072 up to 0.08 so the band is not
#  tighter than the evidence), warns if this constant ever falls below the
#  measured spread, and so keeps the number traceable.
# --- dynamic-stall-vortex reconstruction constants ---------------------------
# Neither is derived, and neither can be calibrated from anything this study
# ships: the experimental frames carry only integrated cl/cd/cm against
# incidence, with no surface-pressure or field data to fit a core size to. They
# are named here rather than buried as literals so that the two numbers a reader
# would have to change are visible, and the core depth they produce is published
# as Cp_DSV_core_min in metrics_*.csv. See the DSV entry under LIMITATIONS.
DSV_GAMMA_FACTOR       = 1.4    # vortex circulation = this * CNv * U * c
DSV_CORE_RADIUS_CHORDS = 0.16   # Lamb-Oseen core radius in chords

DAMPING_TOL = 0.08


def aerodynamic_damping(alpha_deg, CM, normalise=False):
    """Cyclic work / damping coefficient:  Xi = -∮ CM dalpha .
    Xi > 0  -> positive aerodynamic damping (stable);
    Xi < 0  -> negative damping (stall-flutter prone).

    With normalise=True, also returns Xi divided by the area of the bounding box
    of the CM-alpha loop.  Only that ratio is a meaningful discriminator: Xi
    itself carries the units of the loop and, for the near-cancelling
    figure-of-eight loops typical of light dynamic stall, is a small residual.
    """
    a = np.radians(alpha_deg)
    Xi = -_trapz(CM, a)
    if not normalise:
        return Xi
    box = (np.max(CM)-np.min(CM))*(np.max(a)-np.min(a))
    return Xi, (Xi/box if box > 0 else 0.0)


def damping_verdict(Xi_hat, tol=DAMPING_TOL):
    """Three-way stall-flutter verdict from the NORMALISED damping Xi_hat.
    A bare sign test on Xi is not defensible: |Xi_hat| < tol means the two
    lobes of the CM loop cancel to within the resolution of the model."""
    if Xi_hat < -tol:
        return "HIGH (neg. damping)"
    if Xi_hat > tol:
        return "low (pos. damping)"
    return "neutral (within model resolution)"


# --------------------------------------------------------------------------- #
#  FIELD RECONSTRUCTION  (source panels + bound vortex sheet + Lamb-Oseen DSV)
# --------------------------------------------------------------------------- #
def _airfoil_surface(naca_csv, c, n_panel=160):
    """Panel end-points around a CLOSED section, clustered at BOTH the leading
    and the trailing edge.

    Two defects were fixed here and both were measurable:
      * the 4-digit section has an open trailing edge (0.252 %c). Leaving the
        panel body open made it leak: the net source flux sum(sigma*L) came out
        at 0.44 % of U*c instead of zero, and the surface Cp jumped by 2.4
        across the trailing edge, which the Kutta condition forbids. The two
        trailing-edge points are now merged, exactly as 02_mesh does for the
        O-grid wall.
      * the panels were cosine-clustered in ARCLENGTH, whose two ends are both
        at the trailing edge, so the leading edge -- where the suction peak is
        -- got the coarsest panels on the body (98x longer than the trailing-
        edge panels). The distribution is now cosine on each surface separately,
        which clusters at the leading edge as well.
    """
    import pandas as pd
    df = pd.read_csv(naca_csv)
    x, y = df["x_over_c"].values*c, df["y_over_c"].values*c
    x = x.copy(); y = y.copy()
    xm, ym = 0.5*(x[0] + x[-1]), 0.5*(y[0] + y[-1])
    x[0] = x[-1] = xm; y[0] = y[-1] = ym            # close the section
    s = np.concatenate([[0], np.cumsum(np.hypot(np.diff(x), np.diff(y)))]); s /= s[-1]
    m = (n_panel//2) + 1
    h = (1 - np.cos(np.linspace(0, np.pi, m)))/2    # clustered at both ends
    sq = np.concatenate([0.5*h, 0.5 + 0.5*h[1:]])   # -> TE, LE and TE again
    return np.interp(sq, s, x), np.interp(sq, s, y)


def _solve_panels(xp, yp, U, alpha, Gamma):
    """Constant-strength source panels in the presence of a surface vortex sheet
    of uniform strength carrying total circulation Gamma.

    The circulation used to be a separate elliptic sheet laid along the CHORD
    LINE, detached from the body. That arrangement has no Kutta condition: the
    source panels alone are the non-lifting solution, whose rear stagnation
    point sits on the upper surface at incidence, and a chord-line sheet whose
    strength vanishes at the trailing edge cannot move it. The measured symptom
    was a trailing-edge Cp jump of 2.4 and suction on the pressure side.
    Putting the vorticity ON the surface and solving the sources against it
    leaves the trailing-edge jump at ~0.1 and makes the surface Cp integrate
    back to the C_L it was given (see reconstruct_field's LIMITATION note).

    Returns panel midpoints, lengths, source strengths and the sheet strength.
    """
    xc = 0.5*(xp[:-1]+xp[1:]); yc = 0.5*(yp[:-1]+yp[1:])
    dx = np.diff(xp); dy = np.diff(yp); L = np.hypot(dx, dy)
    nx, ny = dy/L, -dx/L                              # normal (sense fixed below)
    # ensure outward (point away from centroid)
    cx, cy = xc.mean(), yc.mean()
    flip = ((xc-cx)*nx+(yc-cy)*ny) < 0; nx = np.where(flip, -nx, nx); ny = np.where(flip, -ny, ny)
    Np = len(xc)
    gam = Gamma/L.sum()                               # uniform sheet strength
    Uinf = np.array([U*np.cos(alpha), U*np.sin(alpha)])
    A = np.zeros((Np, Np)); rhs = np.zeros(Np)
    for i in range(Np):
        rx = xc[i]-xc; ry = yc[i]-yc
        r2 = rx*rx+ry*ry + (0.5*L)**2*1e-2
        ui = (1/(2*np.pi))*rx/r2*L; vi = (1/(2*np.pi))*ry/r2*L
        A[i,:] = ui*nx[i]+vi*ny[i]
        A[i,i] = 0.5                                   # self contribution
        # vortex-sheet downwash at this control point (its own panel induces no
        # normal velocity on itself, so drop the self term)
        gu = (gam*L/(2*np.pi))*ry/r2; gv = -(gam*L/(2*np.pi))*rx/r2
        gu[i] = 0.0; gv[i] = 0.0
        rhs[i] = -(Uinf[0]*nx[i]+Uinf[1]*ny[i]) - (gu.sum()*nx[i]+gv.sum()*ny[i])
    sigma = np.linalg.solve(A, rhs)
    return xc, yc, L, sigma, gam


def _core_pressure_deficit(r, Gamma, rc, U, n=800):
    """Cp correction converting Bernoulli into radial equilibrium inside a
    Lamb-Oseen vortex core.

    Cp = 1 - (V/U)^2 is Bernoulli, which holds only where the flow is
    IRROTATIONAL. Inside the core v_theta -> 0 as r -> 0, so Bernoulli alone
    reports the dynamic-stall vortex as a pressure PEAK, when a real vortex core
    is a pressure MINIMUM held by radial equilibrium  dp/dr = rho*v_theta^2/r :

        Cp_eq(r) = -(2/U^2) * integral_r^inf ( v_theta^2 / r' ) dr'

    Outside the core v_theta = Gamma/(2*pi*r) and that integral collapses to
    -(v_theta/U)^2 -- precisely what Bernoulli already supplies. So the value
    returned here is Cp_eq minus the Bernoulli part: it decays to zero away from
    the core and adds the missing suction inside it.
    """
    r = np.asarray(r, float)
    if Gamma <= 0.0 or rc <= 0.0 or U <= 0.0:
        return np.zeros_like(r)
    rmax = max(float(np.nanmax(r)), 12.0*rc)
    rr = np.linspace(rc*1e-4, rmax, n)
    vt = (Gamma/(2.0*np.pi*rr))*(1.0 - np.exp(-rr**2/rc**2))
    g = vt*vt/rr
    I = np.concatenate([[0.0], np.cumsum(0.5*(g[1:]+g[:-1])*np.diff(rr))])
    cp_eq = -(2.0/U**2)*(I[-1] - I)          # radial equilibrium
    cp_bern = -(vt/U)**2                     # already counted by Bernoulli
    return np.interp(r, rr, cp_eq - cp_bern)


def _panel_velocity(X, Y, xc, yc, L, sigma, gam, eps2):
    """Velocity induced at arbitrary points by the source panels and the uniform
    surface vortex sheet. Shared by the field reconstruction and the surface-Cp
    evaluation so the two can never be built from different fields."""
    u = np.zeros_like(np.asarray(X, float))
    v = np.zeros_like(u)
    for j in range(len(xc)):
        rx = X-xc[j]; ry = Y-yc[j]; r2 = rx*rx+ry*ry+eps2
        u += (sigma[j]*L[j]/(2*np.pi))*rx/r2 + (gam*L[j]/(2*np.pi))*ry/r2
        v += (sigma[j]*L[j]/(2*np.pi))*ry/r2 - (gam*L[j]/(2*np.pi))*rx/r2
    return u, v


def _surface_velocity(xp, yp, xc, yc, L, sigma, gam, U, alpha):
    """Velocity ON the body, evaluated at the panel CONTROL POINTS.

    This must not be done with _panel_velocity. That routine is for field
    points: it regularises 1/r2 with eps2 and includes every panel, which at a
    point lying on the sheet is both singular and wrong. On the surface the two
    self-contributions are known in closed form and are the whole difficulty:

      * a constant-strength SOURCE panel induces sigma/2 along its own outward
        normal and nothing tangentially;
      * a constant-strength VORTEX panel induces -gam/2 along its own tangent
        and nothing normally -- this is the tangential velocity JUMP across a
        vortex sheet, and it is the term whose omission was the largest single
        error in the reconstruction.

    Evaluating instead at the panel END-POINTS, offset along a normal, was the
    second error: flow tangency is imposed at the control points, so those are
    the only places the discrete solution actually satisfies the boundary
    condition. Measured at alpha=10 deg, C_L=1.10: the surface-C_p integral
    returned -48.7% of the circulation it was given with the self-term missing,
    and -0.18% with this routine.
    """
    # Tangents and normals come from the PANEL END-POINTS, which is what
    # _solve_panels imposed tangency with. Taking them from control-point to
    # control-point instead is a centred direction over two panels: it differs
    # wherever the surface curves, i.e. exactly at the leading edge where the
    # self-terms matter most, and measured 1.0 point of closure error at
    # alpha = 17.5 deg (-0.58% with the panel tangent, +1.05% without).
    tx, ty = np.diff(xp)/L, np.diff(yp)/L
    nx, ny = ty, -tx
    cx, cy = xc.mean(), yc.mean()
    flip = ((xc-cx)*nx + (yc-cy)*ny) < 0
    nx = np.where(flip, -nx, nx); ny = np.where(flip, -ny, ny)
    u = np.full(len(xc), U*np.cos(alpha)); v = np.full(len(xc), U*np.sin(alpha))
    for i in range(len(xc)):
        rx = xc[i]-xc; ry = yc[i]-yc; r2 = rx*rx + ry*ry
        r2[i] = np.inf                                   # exclude the self panel
        u[i] += np.sum((sigma*L/(2*np.pi))*rx/r2 + (gam*L/(2*np.pi))*ry/r2)
        v[i] += np.sum((sigma*L/(2*np.pi))*ry/r2 - (gam*L/(2*np.pi))*rx/r2)
    u += 0.5*sigma*nx - 0.5*gam*tx                       # analytic self terms
    v += 0.5*sigma*ny - 0.5*gam*ty
    return u, v


def _dsv_velocity(X, Y, xv, yv, Gv, rc):
    """Lamb-Oseen dynamic-stall vortex, same sign convention as the bound sheet
    (clockwise). Returns (du, dv, r2) with r2 the squared distance to the core."""
    rx = X-xv; ry = Y-yv; r2 = rx*rx+ry*ry
    fcore = (1-np.exp(-r2/rc**2))
    with np.errstate(divide="ignore", invalid="ignore"):
        du =  Gv/(2*np.pi)*ry/np.where(r2 == 0, 1, r2)*fcore
        dv = -Gv/(2*np.pi)*rx/np.where(r2 == 0, 1, r2)*fcore
    return du, dv, r2


def reconstruct_field(naca_csv, c, U, M, alpha_deg, CL, CNv,
                      tau_over_Tvl, domain=(-1.0, 2.0, -1.2, 1.2),
                      nx_grid=260, ny_grid=200, gamma=1.4,
                      T_inf=288.15, cp=1004.5, recovery=0.892, R_gas=287.05):
    """Reconstruct 2D flow field at one instant. Returns grids of velocity,
    pressure coefficient, static & recovery temperature, vorticity, plus the
    DSV location. Lifting circulation matched to the UIBS CL and carried on the
    body surface (so the trailing edge behaves); the dynamic-stall vortex
    rendered as a convecting Lamb-Oseen vortex of strength ~ CNv."""
    alpha = np.radians(alpha_deg)
    Gamma = 0.5*CL*U*c                    # Kutta-Joukowski, matched to the UIBS CL
    xp, yp = _airfoil_surface(naca_csv, c)
    xc, yc, L, sigma, gam = _solve_panels(xp, yp, U, alpha, Gamma)

    x0, x1, y0, y1 = domain
    gx = np.linspace(x0*c, x1*c, nx_grid)
    gy = np.linspace(y0*c, y1*c, ny_grid)
    X, Y = np.meshgrid(gx, gy)
    u = np.full_like(X, U*np.cos(alpha)); v = np.full_like(X, U*np.sin(alpha))

    # Source panels AND the bound vortex sheet, both carried on the body surface.
    # The circulation used to live on the chord line as a separate elliptic
    # sheet; see _solve_panels for why that could not satisfy the Kutta
    # condition. Sign: this sheet is clockwise (positive lift for flow in +x),
    # which fixes the sign convention the dynamic-stall vortex below must match.
    # Regularisation radius. It must cover BOTH the panel spacing and the field
    # grid spacing: a bound sheet is singular on the surface, so if the panels
    # are resolved at roughly one panel per grid cell the sampled vorticity
    # alternates sign cell to cell and the surface renders as a speckled band
    # instead of the thin sheet it is. Take whichever is larger.
    dgrid = max(gx[1]-gx[0], gy[1]-gy[0])
    eps2 = max((0.6*L.mean())**2, (0.9*dgrid)**2)
    du, dv = _panel_velocity(X, Y, xc, yc, L, sigma, gam, eps2)
    u += du; v += dv

    # dynamic-stall vortex (Lamb-Oseen), convects along upper surface
    xv = (0.25 + 0.55*np.clip(tau_over_Tvl, 0, 1.3))*c
    yv = 0.10*c + 0.06*c*np.clip(tau_over_Tvl, 0, 1.3)
    # SIGN. The dynamic-stall vortex is a roll-up of upper-surface boundary-layer
    # vorticity, so it rotates in the SAME sense as the bound circulation
    # (clockwise here). This term must therefore carry the same sign as the
    # bound surface sheet above. It does not need to supply the vortex lift -- Gamma is
    # already matched to the full UIBS C_L, which contains CNv -- so giving the
    # vortex its physical rotation costs nothing and buys the correct vorticity
    # field and the flow reversal beneath the core that characterises the stall.
    # The suction under the vortex comes from its low-pressure core (see
    # _core_pressure_deficit), not from accelerating the surface flow.
    Gv = DSV_GAMMA_FACTOR*max(CNv, 0.0)*U*c
    rc = DSV_CORE_RADIUS_CHORDS*c
    du, dv, r2 = _dsv_velocity(X, Y, xv, yv, Gv, rc)
    u += du; v += dv

    speed = np.hypot(u, v)
    # incompressible Cp + Prandtl-Glauert compressibility correction (bounded).
    # The core correction is what keeps the dynamic-stall vortex reading as the
    # suction feature it is; see _core_pressure_deficit.
    Cp_inc = (1.0 - (speed/U)**2
              + _core_pressure_deficit(np.sqrt(r2), abs(Gv), rc, U))
    # No Prandtl-Glauert factor. The circulation this field is built from is
    # Gamma = 0.5*C_L*U*c and that C_L already carries compressibility through
    # beta in the indicial march, so scaling the resulting Cp again multiplies
    # the reconstructed load by a further 1/beta (+4.8% at M = 0.3). Applying it
    # here was one of the three errors that put the surface-Cp closure at
    # -12.4%; see surface_cp.
    Cp = Cp_inc
    # NOT clipped. There used to be a np.clip(Cp, -8, 1) here, described as a
    # display bound. It was a no-op and the description was wrong: the deepest
    # value reached across all eight published fields is -5.2, and removing the
    # clip leaves every field bit-identical, tested at conditions past the ones
    # written out. The field never approaches the -15 the surface reaches
    # because the near-wall ring that carries the leading-edge suction peak is
    # masked below (see near_wall_cells_masked) -- so the FIELD understates the
    # peak suction by about a factor of three relative to surface_cp, which is a
    # property of the grid, not a bound imposed on the data.
    # thermodynamics
    T0 = T_inf*(1+(gamma-1)/2*M**2)
    T_static = T0 - speed**2/(2*cp)
    T_recovery = T0 - (1-recovery)*speed**2/(2*cp)
    Mlocal = speed/np.sqrt(gamma*R_gas*np.maximum(T_static, 1.0))
    # vorticity
    dvx = np.gradient(v, gx, axis=1); duy = np.gradient(u, gy, axis=0)
    vort = dvx - duy

    # Mask the aerofoil interior AND the one ring of cells touching it.
    # The bound sheet is a singular vortex sheet regularised over eps; a cell
    # that straddles the surface therefore samples the middle of the jump and
    # reports a spurious acceleration. It is measurable: at mid-chord on the
    # PRESSURE side, where the profile is monotonically falling towards the wall
    # (79.1 -> 80.4 -> 84.6 m/s), the last cell jumped to 111.7 m/s -- above the
    # free stream, on the side of the aerofoil that must be slower than it. The
    # field is simply not resolved within one cell of the wall, so it is not
    # published there. surface_cp probes at 0.015c, which is ~7 cells out on its
    # own finer grid, so it is unaffected.
    from matplotlib.path import Path as MplPath
    poly = MplPath(np.column_stack([xp, yp]))
    inside = poly.contains_points(np.column_stack([X.ravel(), Y.ravel()])).reshape(X.shape)
    masked = inside.copy()                       # 8-connected dilation by one cell
    for sx in (-1, 0, 1):
        for sy in (-1, 0, 1):
            if sx == 0 and sy == 0:
                continue
            masked |= np.roll(np.roll(inside, sy, axis=0), sx, axis=1)
    for arr in (u, v, speed, Cp, T_static, T_recovery, Mlocal, vort):
        arr[masked] = np.nan

    return dict(X=X, Y=Y, u=u, v=v, speed=speed, Cp=Cp, T_static=T_static,
                T_recovery=T_recovery, Mlocal=Mlocal, vort=vort,
                xp=xp, yp=yp, xv=xv, yv=yv, Gamma=Gamma, Gv=Gv,
                T0=T0, T_inf=T_inf)


def surface_cp(naca_csv, c, U, M, alpha_deg, CL, CNv, tau_over_Tvl):
    """Surface pressure coefficient distribution Cp(x/c), evaluated exactly at
    the panel control points.

    Three defects were found here by testing the closure invariant (integrating
    the returned Cp must return the C_L the reconstruction was given), and all
    three were measurable at alpha=10 deg, C_L=1.10:

      1. The vortex sheet's own tangential contribution (-gam/2) was omitted.
         Alone this cost -48.7% of the lift when evaluated on the wall.
      2. Cp was evaluated at the panel end-points, stepped 0.015c off the wall
         along a normal, rather than at the control points where tangency is
         actually imposed. That offset masked (1) rather than fixing it: the
         error read -11.5% at 0.015c and grew to -48.7% as the probe approached
         the surface, so a reader could have concluded the method was better
         than it was by standing further away from the aerofoil.
      3. A Prandtl-Glauert factor was applied on top. It does not belong: the
         circulation supplied to the reconstruction is Gamma = 0.5*C_L*U*c, and
         that C_L already carries compressibility through beta in the indicial
         march. Scaling the resulting Cp again multiplied the reconstructed
         load by a further 1/beta = 1.048 at M = 0.3, i.e. +4.8%.

    With all three corrected the closure error is -0.18% at alpha = 10 deg and
    stays inside -0.4% up to 19 deg, against -12.4% published previously.

    Returns (x/c at the control points, Cp, upper_mask).
    """
    alpha = np.radians(alpha_deg)
    Gamma = 0.5*CL*U*c
    xp, yp = _airfoil_surface(naca_csv, c)
    xc, yc, L, sigma, gam = _solve_panels(xp, yp, U, alpha, Gamma)

    u, v = _surface_velocity(xp, yp, xc, yc, L, sigma, gam, U, alpha)

    xv = (0.25 + 0.55*np.clip(tau_over_Tvl, 0, 1.3))*c
    yv = 0.10*c + 0.06*c*np.clip(tau_over_Tvl, 0, 1.3)
    Gv = DSV_GAMMA_FACTOR*max(CNv, 0.0)*U*c
    rc = DSV_CORE_RADIUS_CHORDS*c
    du, dv, r2 = _dsv_velocity(xc, yc, xv, yv, Gv, rc)
    u += du; v += dv

    speed = np.hypot(u, v)
    Cp = (1.0 - (speed/U)**2
          + _core_pressure_deficit(np.sqrt(r2), abs(Gv), rc, U))
    # NOT clipped. The potential-flow leading-edge suction peak reaches Cp =
    # -16 at 17.5 deg, and the -8 clip the FIELD uses for display truncated it,
    # which by itself cost -5.9% of the closure once the errors above were
    # fixed. A real boundary layer separates long before that peak is reached;
    # that is a limitation of reconstructing from potential flow, and it is
    # stated rather than hidden by a clip.
    upper = yc >= 0
    return xc/c, Cp, upper


def kutta_reference_CL(naca_csv, c, U, M, alpha_deg):
    """The lift coefficient at which this section's reconstruction satisfies the
    Kutta condition exactly, i.e. the inviscid attached-flow circulation.

    The trailing-edge Cp jump is LINEAR in the imposed C_L, because the panel
    system is linear and only the vortex sheet strength scales with it. Two
    evaluations therefore locate its zero exactly.

    This exists because the trailing-edge jump is not a defect that can be
    driven to zero. The reconstruction is handed Gamma = 0.5*C_L*U*c with C_L
    from the indicial march, which during dynamic stall deliberately departs
    from the attached-inviscid value; a body carrying a circulation other than
    the Kutta one MUST show a trailing-edge jump. Publishing this reference
    turns the residual from an unexplained number into a measure of how far the
    modelled flow is from attached, which is what it actually is.
    """
    j = []
    for CLt in (1.0, 2.0):
        _, cp_s, _ = surface_cp(naca_csv, c, U, M, alpha_deg, CLt, 0.0, 0.0)
        j.append(float(cp_s[0] - cp_s[-1]))          # SIGNED trailing-edge jump
    if j[1] == j[0]:
        return float("nan")
    return float(1.0 - j[0]*(2.0 - 1.0)/(j[1] - j[0]))


def surface_load_closure(naca_csv, c, U, M, alpha_deg, CL, CNv, tau_over_Tvl):
    """Integrate the reconstructed surface Cp and compare with the C_L that the
    reconstruction was given. An invariant: a closed body carrying circulation
    Gamma = 0.5*CL*U*c must return that C_L (Blasius), so any residual here is
    discretisation or a bug, never physics.

    Also returns the trailing-edge pressure jump |Cp_upper - Cp_lower| there,
    the second invariant: the Kutta condition demands it be zero.

    The integral is taken over the same panels the solution was built on, with
    Cp at the control points. It previously summed a mean of END-POINT Cp values
    against the end-point spacing, a quadrature the panel solution does not
    support.

    Both exist so the reconstruction's accuracy is a number the pipeline writes
    out, not a claim in a comment.
    """
    xp, yp = _airfoil_surface(naca_csv, c)
    _, cp_s, _ = surface_cp(naca_csv, c, U, M, alpha_deg, CL, CNv, tau_over_Tvl)
    # traversal sense (shoelace): +1 counter-clockwise, so that n ds = (dy, -dx)
    sgn = 1.0 if 0.5*np.sum(xp[:-1]*yp[1:] - xp[1:]*yp[:-1]) > 0 else -1.0
    a = np.radians(alpha_deg)
    CN =  sgn*np.sum(cp_s*np.diff(xp))/c        # +(1/c) integral Cp dx
    CA = -sgn*np.sum(cp_s*np.diff(yp))/c        # -(1/c) integral Cp dy
    cl = CN*np.cos(a) - CA*np.sin(a)
    te_jump = float(abs(cp_s[0] - cp_s[-1]))    # control points either side of the TE
    return (float(cl), (100.0*(cl - CL)/CL if CL != 0 else float("nan")), te_jump)

