"""
=============================================================================
 UNISTALL(TM)  —  Universal Unsteady-Aerodynamics & Dynamic-Stall Solver
 Core method : Unified Indicial-Beddoes State-Space (UIBS)
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
    streamlines + an explicit Lamb-Oseen dynamic-stall vortex.
    LIMITATION: the reconstructed vortex carries circulation 1.4*CNv*U*c in a
    core of radius 0.16c, giving a peak swirl of roughly 0.2*U. Its sign,
    position and the flow reversal beneath it are physical, but the core is
    weaker and more diffuse than a measured dynamic-stall vortex, so the core
    suction it produces is shallow. The reconstruction is qualitative; the
    reported loads come from the UIBS core and do not depend on it.
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
    TI = Kalpha*c/(M*340.0 if M > 0 else U)   # impulsive time const (s), a=U/M

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
#  |Xi_hat| below this is treated as neutral: for a figure-of-eight CM loop the
#  raw Xi is a small residual between two lobes of opposite sign, and at this
#  level it is no larger than the time-step discretisation error.
DAMPING_TOL = 0.02


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
    import pandas as pd
    df = pd.read_csv(naca_csv)
    x, y = df["x_over_c"].values*c, df["y_over_c"].values*c
    # resample to n_panel evenly along arclength
    s = np.concatenate([[0], np.cumsum(np.hypot(np.diff(x), np.diff(y)))]); s/=s[-1]
    sq = (1-np.cos(np.linspace(0, np.pi, n_panel+1)))/2
    return np.interp(sq, s, x), np.interp(sq, s, y)

def _solve_sources(xp, yp, U, alpha):
    """Constant-strength source panels: enforce flow tangency (point-source
    approximation, regularized self term)."""
    xc = 0.5*(xp[:-1]+xp[1:]); yc = 0.5*(yp[:-1]+yp[1:])
    dx = np.diff(xp); dy = np.diff(yp); L = np.hypot(dx, dy)
    nx, ny = dy/L, -dx/L                              # outward normal (CW airfoil)
    # ensure outward (point away from centroid)
    cx, cy = xc.mean(), yc.mean()
    flip = ((xc-cx)*nx+(yc-cy)*ny) < 0; nx[flip]*=-1; ny[flip]*=-1
    Np = len(xc)
    Uinf = np.array([U*np.cos(alpha), U*np.sin(alpha)])
    A = np.zeros((Np, Np)); rhs = np.zeros(Np)
    for i in range(Np):
        rx = xc[i]-xc; ry = yc[i]-yc
        r2 = rx*rx+ry*ry + (0.5*L)**2*1e-2
        ui = (1/(2*np.pi))*rx/r2*L; vi = (1/(2*np.pi))*ry/r2*L
        A[i,:] = ui*nx[i]+vi*ny[i]
        A[i,i] = 0.5                                   # self contribution
        rhs[i] = -(Uinf[0]*nx[i]+Uinf[1]*ny[i])
    sigma = np.linalg.solve(A, rhs)
    return xc, yc, L, sigma

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


def reconstruct_field(naca_csv, c, U, M, alpha_deg, CL, CNv,
                      tau_over_Tvl, domain=(-1.0, 2.0, -1.2, 1.2),
                      nx_grid=260, ny_grid=200, gamma=1.4,
                      T_inf=288.15, cp=1004.5, recovery=0.892):
    """Reconstruct 2D flow field at one instant. Returns grids of velocity,
    pressure coefficient, static & recovery temperature, vorticity, plus the
    DSV location. Lifting circulation matched to the UIBS CL; the dynamic-stall
    vortex rendered as a convecting Lamb-Oseen vortex of strength ~ CNv."""
    alpha = np.radians(alpha_deg)
    xp, yp = _airfoil_surface(naca_csv, c)
    xc, yc, L, sigma = _solve_sources(xp, yp, U, alpha)

    x0, x1, y0, y1 = domain
    gx = np.linspace(x0*c, x1*c, nx_grid)
    gy = np.linspace(y0*c, y1*c, ny_grid)
    X, Y = np.meshgrid(gx, gy)
    u = np.full_like(X, U*np.cos(alpha)); v = np.full_like(X, U*np.sin(alpha))

    # source-panel induced velocity (point-source per panel, smoothed)
    eps2 = (0.6*L.mean())**2
    for j in range(len(xc)):
        rx = X-xc[j]; ry = Y-yc[j]; r2 = rx*rx+ry*ry+eps2
        u += (sigma[j]*L[j]/(2*np.pi))*rx/r2
        v += (sigma[j]*L[j]/(2*np.pi))*ry/r2

    # bound circulation: smooth (elliptic) vortex sheet along chord (y=0), 0..c.
    # Elliptic loading is finite at both ends -> no leading-edge singularity,
    # giving a clean field; total circulation is matched to the UIBS lift.
    Gamma = 0.5*CL*U*c
    nb = 80
    xb = np.linspace(0.01*c, 0.99*c, nb)
    w = np.sqrt(np.clip((xb/c)*(1-xb/c), 0, None))       # elliptic loading
    w /= _trapz(w, xb); dGam = Gamma*w*np.gradient(xb)
    epsb2 = (0.06*c)**2
    for j in range(nb):
        rx = X-xb[j]; ry = Y-0.0; r2 = rx*rx+ry*ry+epsb2
        u +=  dGam[j]/(2*np.pi)*ry/r2
        v += -dGam[j]/(2*np.pi)*rx/r2

    # dynamic-stall vortex (Lamb-Oseen), convects along upper surface
    xv = (0.25 + 0.55*np.clip(tau_over_Tvl, 0, 1.3))*c
    yv = 0.10*c + 0.06*c*np.clip(tau_over_Tvl, 0, 1.3)
    # SIGN. The dynamic-stall vortex is a roll-up of upper-surface boundary-layer
    # vorticity, so it rotates in the SAME sense as the bound circulation
    # (clockwise here). This term must therefore carry the same sign as the
    # bound sheet above. It does not need to supply the vortex lift -- Gamma is
    # already matched to the full UIBS C_L, which contains CNv -- so giving the
    # vortex its physical rotation costs nothing and buys the correct vorticity
    # field and the flow reversal beneath the core that characterises the stall.
    # The suction under the vortex comes from its low-pressure core (see
    # _core_pressure_deficit), not from accelerating the surface flow.
    Gv = 1.4*max(CNv, 0.0)*U*c
    rc = 0.16*c
    rx = X-xv; ry = Y-yv; r2 = rx*rx+ry*ry
    fcore = (1-np.exp(-r2/rc**2))
    with np.errstate(divide="ignore", invalid="ignore"):
        u +=  Gv/(2*np.pi)*ry/np.where(r2 == 0, 1, r2)*fcore
        v += -Gv/(2*np.pi)*rx/np.where(r2 == 0, 1, r2)*fcore

    speed = np.hypot(u, v)
    # incompressible Cp + Prandtl-Glauert compressibility correction (bounded).
    # The core correction is what keeps the dynamic-stall vortex reading as the
    # suction feature it is; see _core_pressure_deficit.
    Cp_inc = (1.0 - (speed/U)**2
              + _core_pressure_deficit(np.sqrt(r2), abs(Gv), rc, U))
    Cp = Cp_inc/np.sqrt(1-M**2) if M > 0 else Cp_inc
    Cp = np.clip(Cp, -8.0, 1.0)               # keep field physical & clean
    # thermodynamics
    T0 = T_inf*(1+(gamma-1)/2*M**2)
    T_static = T0 - speed**2/(2*cp)
    T_recovery = T0 - (1-recovery)*speed**2/(2*cp)
    Mlocal = speed/np.sqrt(gamma*287.05*np.maximum(T_static, 1.0))
    # vorticity
    dvx = np.gradient(v, gx, axis=1); duy = np.gradient(u, gy, axis=0)
    vort = dvx - duy

    # mask airfoil interior
    from matplotlib.path import Path as MplPath
    poly = MplPath(np.column_stack([xp, yp]))
    inside = poly.contains_points(np.column_stack([X.ravel(), Y.ravel()])).reshape(X.shape)
    for arr in (u, v, speed, Cp, T_static, T_recovery, Mlocal, vort):
        arr[inside] = np.nan

    return dict(X=X, Y=Y, u=u, v=v, speed=speed, Cp=Cp, T_static=T_static,
                T_recovery=T_recovery, Mlocal=Mlocal, vort=vort,
                xp=xp, yp=yp, xv=xv, yv=yv, Gamma=Gamma, Gv=Gv,
                T0=T0, T_inf=T_inf)


def surface_cp(naca_csv, c, U, M, alpha_deg, CL, CNv, tau_over_Tvl,
               nx_grid=400, ny_grid=300):
    """Surface pressure coefficient distribution Cp(x/c) upper & lower."""
    fld = reconstruct_field(naca_csv, c, U, M, alpha_deg, CL, CNv, tau_over_Tvl,
                            nx_grid=nx_grid, ny_grid=ny_grid)
    xp, yp = fld["xp"], fld["yp"]
    # sample just outside the surface along outward normals
    dx = np.gradient(xp); dy = np.gradient(yp); Ln = np.hypot(dx, dy)
    nx, ny = dy/Ln, -dx/Ln
    cx, cy = xp.mean(), yp.mean()
    flip = ((xp-cx)*nx+(yp-cy)*ny) < 0; nx[flip]*=-1; ny[flip]*=-1
    off = 0.015*c
    xs, ys = xp+nx*off, yp+ny*off
    from scipy.interpolate import RegularGridInterpolator
    gx = fld["X"][0,:]; gy = fld["Y"][:,0]
    Cpf = np.nan_to_num(fld["Cp"], nan=0.0)
    itp = RegularGridInterpolator((gy, gx), Cpf, bounds_error=False, fill_value=0.0)
    cp = itp(np.column_stack([ys, xs]))
    upper = yp >= 0
    return xp/c, cp, upper
