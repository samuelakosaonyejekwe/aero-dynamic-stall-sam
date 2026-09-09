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
    the metrics_*.csv rows "Cp_closure_error_pct",
    "Cp_closure_worst_dCL_cycle" and "Cp_TE_jump_max_over_phases"
    recompute them on every run):
      - CLOSURE. Integrating the surface Cp now recovers the C_L it was given to
        -0.28 % at alpha 2 deg (C_L 0.22), -0.37 % at 10 deg (C_L 1.10) and
        -0.58 % at 17.5 deg (C_L 1.91) with
        the dynamic-stall vortex switched off -- which is the published
        Cp_closure_error_pct, and the DSV-OFF measurement is the only one that
        is a closure error at all.
        WITH the vortex present the surface integral does not return the imposed
        C_L exactly, and it should not: a free vortex near the body exerts a
        real force on it, so that part of the residual is the vortex's induced
        lift, not an error. It is published as its own row,
        DSV_induced_lift_dCL, and at the derived circulation it is small
        (-0.0025 for Case A). The two were conflated while the figure was quoted
        as -0.8 % (Case A) / -0.6 % (Case B) and read as a discretisation bound.
        Refinement settles which is which: the DSV-off closure converges to
        zero, the with-DSV one converges to the interaction force. Those are single instants and are NOT a bound on
        the cycle; an earlier revision claimed "within 0.6 % over the whole
        cycle", which was never measured. What IS measured over the cycle, on
        every run, is the worst ABSOLUTE residual with the vortex off:
        Cp_closure_worst_dCL_cycle, i.e. its percentage of each case's own C_L,max
        (Cp_closure_worst_dCL_pct_of_CLmax). It is reported in C_L counts rather
        than as the worst instantaneous percentage because the cycle passes
        through C_L = 0.09, where a residual of 0.0014 reads as +1.6 % purely
        from the small denominator. Unlike before the closure also CONVERGES:
        refining 160 -> 1280 panels drives it monotonically to -0.04 %, which is
        what Blasius requires and is the check that the formulation is right
        rather than merely better. It previously read -10.7 to -14.4 %, and
        the explanation recorded here for that
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
        drives the trailing-edge jump to ~1e-3. The reconstruction is instead
        handed the indicial C_L, which during dynamic
        stall departs from that value deliberately -- so a body carrying a
        non-Kutta circulation MUST show a trailing-edge jump. The published
        Cp_TE_jump_max_over_phases is therefore a measure of how far the modelled
        flow is from attached, not an error. The earlier, smaller published
        values (0.15 and 0.25) were not a better result: they came from probing
        0.015c off the wall, the same offset that was hiding the closure error.
      - DYNAMIC-STALL VORTEX. Both of its properties are DERIVED; neither is
        chosen. Its circulation is Kutta-Joukowski on its own normal-force
        contribution, Gamma_v = 0.5*CNv*U*c -- the same relation the bound sheet
        uses -- so it carries exactly the CNv/CN share of the circulation the
        lift implies and the budget closes. Its core is the radius at which that
        circulation swirls at the edge speed of the shear layer that rolls it
        up, rc = LAMB_OSEEN_PEAK*Gamma_v/(2*pi*Ve), with Ve taken off the panel
        solution at the vortex centre with the vortex excluded. See the
        constants block for the derivation, and for the shear-layer-flux
        alternative that was tried and rejected for breaking the budget.
        Published in metrics_*.csv, and all of it recomputed by
        verify_invariants: Gamma_v = 0.121*U*c (DSV_circulation_over_Uc),
        rc = 0.0112c (DSV_core_radius_chords), peak swirl 1.10*U
        (DSV_peak_swirl_over_U) and a core suction of Cp = -4.22
        (Cp_DSV_core_min), evaluated AT the vortex centre by dsv_core_cp rather
        than sampled off a grid. Measurements put a dynamic-stall core at -3 to
        -6, so the depth lands where it should without having been fitted there.
        The pair this replaces -- 1.4*CNv*U*c in a 0.16c core -- gave -0.36.
        The core is SMALL, and the published field grid does not resolve it:
        rc is 0.79 of a cell, published as DSV_core_radius_cells. The depth is
        still exact because it is evaluated off the grid, and the vortex is
        still represented correctly in the field because its vorticity is added
        in closed form rather than differenced (integrating the published
        vorticity around the core returns its circulation to 0.4 %).
        It does NOT reverse the flow at the wall beneath it: it induces 0.121*U
        upstream there against a local 1.098*U. Three places in this study
        asserted a reversal the fields never showed; the numbers that settle it
        are now published rather than the adjective.
        It does not affect the reported loads, which come from the UIBS core.
      - Nothing in the reconstruction knows about separation: it is a potential
        field, so at post-stall incidence the leading-edge suction peak it draws
        (about Cp = -15.4 at 17.5 deg) is far deeper than a real separated flow
        would sustain. Nothing is clipped. The FIELD nonetheless bottoms out near
        -5.2, because the near-wall ring carrying that peak is masked, so the
        contour plots understate the surface suction by about three times. That
        -5.2 is a property of the PUBLISHED GRID and not a converged value, and
        it cannot become one: the masked ring is one cell wide, so refining the
        grid moves it closer to a singularity the potential solution genuinely
        has. Measured at the Case-A dsv instant, the field minimum goes -3.3
        (110x85), -5.2 (220x170), -6.5 (440x340), -8.4 (880x680). Any quantity
        that must not depend on the grid is therefore evaluated off it -- the
        surface Cp at the panel control points (surface_cp) and the vortex-core
        depth at the vortex centre (dsv_core_cp).
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
#  MODEL/DISCRETISATION CONSTANTS -- the single source of truth.
#
#  03_model_setup/generate_setup.py IMPORTS these when it writes
#  solver_config.json, instead of restating them. Every one of them used to be a
#  literal in this file AND a literal in that config: f_min, n_panels,
#  domain_chords, the default grid and near_wall_cells_masked. The two copies
#  happened to agree, but nothing made them agree -- the config would have gone
#  on describing a reconstruction the solver had stopped performing.
# --------------------------------------------------------------------------- #
F_MIN = 0.02                       # floor on the separation point f
N_PANELS_DEFAULT = 160             # panels around the closed section
DOMAIN_CHORDS = (-1.0, 2.0, -1.2, 1.2)      # field domain, in chords
GRID_NX_DEFAULT, GRID_NY_DEFAULT = 260, 200  # field grid when none is given
NEAR_WALL_CELLS_MASKED = 1         # rings of cells blanked around the body

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
            f[i] = float(np.clip(sf, np.sqrt(F_MIN), 1.0)**2)
    # monotone, smooth interpolant on |alpha|; clamp ends
    order = np.argsort(alpha_deg)
    ad = np.asarray(alpha_deg, float)[order]; fd = f[order]
    interp = PchipInterpolator(ad, fd, extrapolate=False)
    amin, amax = ad.min(), ad.max()
    S_EXT = 3.0                     # decay length past the data (deg); floor is F_MIN
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
        # Below the data -> attached. This branch does NOT fire for the polar
        # this study ships, which is tabulated from alpha = 0 and is queried on
        # |alpha|, so q >= amin always. It is a live guard rather than a dead
        # line: fed a polar starting at 4 deg it returns f = 1.0 at 0 deg
        # instead of the NaN PchipInterpolator(extrapolate=False) would give,
        # which is the physically right answer below the fitted range. Tested
        # both ways rather than assumed.
        out = np.where(q >= amin, out, 1.0)
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
    # IMPULSIVE TIME CONSTANT. T_I = Kalpha*c/U, i.e. Kalpha*c/(M*a) with
    # a = U/M the speed of sound implied by the case. Expressed in the semichord
    # time the rest of the march runs in, dt/T_I = ds/(2*Kalpha), so the
    # impulsive lag is a fixed 2*Kalpha semichords and the march stays chord-
    # and speed-independent (test_frame_independence in the validation stage
    # asserts that).
    #
    # This is NOT the classical Leishman-Beddoes T_I = c/a: it is larger by 1/M
    # (3.3x at M = 0.30), and the amplitude 4*Kalpha*c/(U*M) below carries the
    # same extra 1/M against the classical 4*Kalpha*c/U. The choice is stated
    # rather than silent because it is a real departure, and it is the
    # convention the calibrated constants in solver_config.json were fitted
    # with: swapping both terms to the classical scaling and re-running the five
    # real NACA 0012 frames with these same constants moves the held-out mean
    # peak-lift error from 1.7 % to 2.9 %. Re-deriving it would require
    # re-calibrating against data this study does not ship, so the implemented
    # form is kept and documented (report section 4.3 states this T_I, not c/a).
    #
    # The denominator was once the literal M*340.0 -- a hardcoded sea-level
    # speed of sound that made the march weakly dependent on U for any case
    # whose speed of sound is not 340 m/s. U = M*a removes that.
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
        fpp[n] = np.clip(fprime - Df[n], F_MIN, 1.0)
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
DAMPING_TOL = 0.08

# --- dynamic-stall vortex: DERIVED, not chosen -------------------------------
# The reconstructed vortex used to carry DSV_GAMMA_FACTOR*CNv*U*c in a core of
# DSV_CORE_RADIUS_CHORDS*c, with both numbers picked rather than derived, and
# the comment here said they could not be derived from anything this study
# ships. That was wrong on both counts.
#
#   CIRCULATION -- Kutta-Joukowski, the SAME relation the bound sheet uses.
#   The sheet carries Gamma = 0.5*C_L*U*c; the vortex's own contribution to the
#   normal force is CNv, so it carries Gamma_v = 0.5*CNv*U*c. The two then sum
#   to the total circulation the lift implies, i.e. the vortex holds exactly the
#   CNv/CN share of it (12.7 % at the Case-A vortex maximum) and the circulation
#   budget closes. It is not free to be larger: sizing the vortex from the
#   shear-layer vorticity flux instead -- Ve^2/2 sustained over the Tvl feeding
#   clock, which gives 1.79*U*c -- was tried and rejected during the audit that
#   derived this, because 1.79*U*c is nearly twice the whole circulation of a
#   section carrying C_L = 1.91 and inflated the body's integrated lift by 60 %.
#
#   CORE RADIUS -- the radius at which that circulation swirls at the speed of
#   the shear layer that rolled it up. A separated shear layer's velocity scale
#   is the edge speed Ve, read off the panel solution AT the vortex centre with
#   the vortex's own field excluded, so
#
#       rc = LAMB_OSEEN_PEAK*Gamma_v/(2*pi*Ve).
#
# Neither expression contains a fitted number. For Case A at the vortex maximum
# they give Gamma_v = 0.121*U*c and rc = 0.0112c, and hence a core suction of
# Cp = -4.22 -- against the -3 to -6 that dynamic-stall measurements report, and
# without having been fitted to it. The chosen pair gave -0.36, an order of
# magnitude too shallow.
#
# WHAT THE VORTEX THEN IS, and is not. It is small and intense rather than
# broad and weak: rc = 0.0112c is 0.79 of a cell of the published field grid, so
# the FIELD does not resolve the core even though its depth is right. That is
# published as DSV_core_radius_cells and is the same kind of statement the
# leading-edge suction peak already carries -- the surface reaches -15 where the
# field bottoms out far shallower, a property of the grid. The vortex's own
# vorticity is added to the field in CLOSED FORM for exactly this reason, so it
# is still sampled correctly: integrating the published vorticity around the
# core returns its circulation to 0.4 % on that grid.
# It does not reverse the flow at the wall beneath it, and that is a measured
# comparison rather than an assertion: it induces 0.121*U upstream there against
# a local 1.098*U (DSV_induced_at_wall_over_U against DSV_edge_speed_over_U).
# Three places in this study once asserted a reversal the fields never showed.
# It does not affect the reported loads, which come from the UIBS core.
#
# LAMB_OSEEN_PEAK is the peak-swirl coefficient of the Lamb-Oseen profile,
# max_u (1 - exp(-u^2))/u, computed rather than quoted. It is NOT 0.7152: that
# number, which stood here briefly, is 1 - exp(-u*^2), the fraction of the
# circulation enclosed at the peak-swirl radius, not the swirl coefficient.
_u = np.linspace(1e-6, 5.0, 200001)
LAMB_OSEEN_PEAK = float(np.max((1.0 - np.exp(-_u*_u))/_u))     # 0.638173
del _u
# Literature default for the vortex-convection clock, in semichords. It is the
# ONE source for that default: 03_model_setup imports it into solver_config.json
# rather than restating it, exactly as it does F_MIN and the panel count. The
# CALIBRATED value overrides it and is what run_case passes in.
TVL_DEFAULT = 5.0


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
def _airfoil_surface(naca_csv, c, n_panel=N_PANELS_DEFAULT):
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


def _solve_panels(xp, yp, U, alpha, Gamma, extra=None):
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

    `extra` is an optional (u, v) pair of velocities induced at the control
    points by anything else in the field -- in practice the dynamic-stall
    vortex. It belongs in the boundary condition: a free vortex added to the
    solution AFTERWARDS drives flow straight through the body, because nothing
    ever asked the sources to cancel its normal component. That was survivable
    while the vortex was a hundredth of the bound circulation; with the derived
    vortex, at 1.79*U*c, it put reversed flow on the PRESSURE side, under the
    aerofoil, where the vortex's own field simply passed through the section.

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
        if extra is not None:                          # e.g. the dynamic-stall
            rhs[i] -= extra[0][i]*nx[i]+extra[1][i]*ny[i]   # vortex, see below
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
    # GEOMETRIC, not uniform. The core is now derived and can be small against
    # rmax; a uniform grid from 0 to rmax then steps clean over it and returns
    # no core suction at all. Log spacing resolves any rc. Checked against the
    # old uniform grid at the old core size: same answer to 1e-5.
    r0 = rc*1e-4
    rr = r0*np.exp(np.linspace(0.0, np.log(rmax/r0), n))
    vt = (Gamma/(2.0*np.pi*rr))*(1.0 - np.exp(-rr**2/rc**2))
    g = vt*vt/rr
    I = np.concatenate([[0.0], np.cumsum(0.5*(g[1:]+g[:-1])*np.diff(rr))])
    # ANALYTIC TAIL beyond rmax. Outside the core v_theta = Gamma/(2 pi r), so
    # the remaining integral is exactly (Gamma/2pi)^2/(2 rmax^2). Without it the
    # correction did not cancel Bernoulli in the far field but left a small
    # POSITIVE residual, which put one cell of one published field at Cp =
    # 1.0005 -- over the stagnation bound the study asserts. With the tail the
    # two agree at rmax by construction and the residual is zero to rounding.
    tail = (Gamma/(2.0*np.pi))**2/(2.0*rmax**2)
    cp_eq = -(2.0/U**2)*((I[-1] + tail) - I)   # radial equilibrium
    cp_bern = -(vt/U)**2                       # already counted by Bernoulli
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
    and -0.37% with this routine. (-0.18% stood here for a while: that was the
    value BEFORE the panel-tangent fix described just below, which is worth a
    point of closure at 17.5 deg and moved this condition too.)
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


def _solve_with_dsv(naca_csv, c, U, alpha, CL, CNv, tau_over_Tvl, Tvl):
    """Panel solution WITH the dynamic-stall vortex in the boundary condition.

    Two steps, because the two depend on each other: the vortex's circulation is
    derived from the edge speed at its centre, which comes from the panel
    solution, and the panel solution must cancel the vortex's normal velocity on
    the body. Solve without the vortex, derive it, then solve again with it --
    one Picard step, and the matrix is unchanged so only the right-hand side is
    rebuilt. The second solve is what keeps the section a streamline; without it
    the vortex's field passes straight through the aerofoil.

    Returns (xp, yp, xc, yc, L, sigma, gam, xv, yv, Gv, rc, Ve).
    """
    Gamma = 0.5*CL*U*c
    xp, yp = _airfoil_surface(naca_csv, c)
    xc, yc, L, sigma, gam = _solve_panels(xp, yp, U, alpha, Gamma)
    xv, yv, Gv, rc, Ve = _dsv_state(xc, yc, L, sigma, gam, U, alpha, c,
                                    CNv, tau_over_Tvl, Tvl)
    if Gv > 0.0:
        du, dv, _ = _dsv_velocity(xc, yc, xv, yv, Gv, rc)
        xc, yc, L, sigma, gam = _solve_panels(xp, yp, U, alpha, Gamma, extra=(du, dv))
    return xp, yp, xc, yc, L, sigma, gam, xv, yv, Gv, rc, Ve


def _dsv_state(xc, yc, L, sigma, gam, U, alpha, c, CNv, tau_over_Tvl, Tvl=None):
    """Position, circulation and core radius of the reconstructed dynamic-stall
    vortex, all derived (see the constants block). Returns (xv, yv, Gv, rc, Ve).

    Tvl is accepted so that every entry point in this module takes the same
    arguments; the derivation does not use it.

    One definition, shared by reconstruct_field, surface_cp and dsv_core_cp.
    The four lines it replaces were written out three times, so a change to the
    vortex had to be made in three places or the field, the surface Cp and the
    published core depth would have described three different vortices.
    """
    t = float(np.clip(tau_over_Tvl, 0.0, 1.3))
    xv = (0.25 + 0.55*t)*c
    yv = (0.10 + 0.06*t)*c
    # edge speed at the vortex centre, with the vortex's OWN field excluded --
    # this is the speed of the shear layer that rolls it up
    du, dv = _panel_velocity(np.array([xv]), np.array([yv]), xc, yc, L, sigma, gam,
                             (0.6*L.mean())**2)
    Ve = float(np.hypot(U*np.cos(alpha) + du[0], U*np.sin(alpha) + dv[0]))
    if CNv <= 0.0 or U <= 0.0 or Ve <= 0.0:
        return xv, yv, 0.0, 1e-6*c, Ve
    Gv = 0.5*CNv*U*c                              # Kutta-Joukowski, as the sheet
    rc = max(LAMB_OSEEN_PEAK*Gv/(2.0*np.pi*Ve), 1e-6*c)
    return xv, yv, Gv, rc, Ve


def _dsv_velocity(X, Y, xv, yv, Gv, rc):
    """Lamb-Oseen dynamic-stall vortex, same sign convention as the bound sheet
    (clockwise). Returns (du, dv, r2) with r2 the squared distance to the core."""
    rx = X-xv; ry = Y-yv; r2 = rx*rx+ry*ry
    fcore = (1-np.exp(-r2/rc**2))
    with np.errstate(divide="ignore", invalid="ignore"):
        du =  Gv/(2*np.pi)*ry/np.where(r2 == 0, 1, r2)*fcore
        dv = -Gv/(2*np.pi)*rx/np.where(r2 == 0, 1, r2)*fcore
    return du, dv, r2


# Fallback air properties, used only when a caller does not pass the values
# 03_model_setup/material_thermo_properties.csv publishes (run_case.py always
# does). They are DERIVED here for the same reason they are derived there:
# the literals that used to sit in this signature, cp = 1004.5 and
# recovery = 0.892, were the very pair that stage removed for contradicting
# their own definitions (gamma*R/(gamma-1) = 1004.68, Pr^(1/3) = 0.8963).
_GAMMA_DEF, _RGAS_DEF, _PR_DEF = 1.4, 287.05, 0.72
_CP_DEF  = _GAMMA_DEF*_RGAS_DEF/(_GAMMA_DEF - 1.0)
_REC_DEF = _PR_DEF**(1.0/3.0)


def reconstruct_field(naca_csv, c, U, M, alpha_deg, CL, CNv,
                      tau_over_Tvl, Tvl=TVL_DEFAULT, domain=DOMAIN_CHORDS,
                      nx_grid=GRID_NX_DEFAULT, ny_grid=GRID_NY_DEFAULT,
                      gamma=_GAMMA_DEF,
                      T_inf=288.15, cp=_CP_DEF, recovery=_REC_DEF,
                      R_gas=_RGAS_DEF):
    """Reconstruct 2D flow field at one instant. Returns grids of velocity,
    pressure coefficient, static & recovery temperature, vorticity, plus the
    DSV location. Lifting circulation matched to the UIBS CL and carried on the
    body surface (so the trailing edge behaves); the dynamic-stall vortex
    rendered as a convecting Lamb-Oseen vortex of strength ~ CNv."""
    alpha = np.radians(alpha_deg)
    # Gamma = 0.5*CL*U*c (Kutta-Joukowski, matched to the UIBS CL) is imposed
    # inside _solve_with_dsv, which also puts the vortex in the wall boundary
    # condition -- see there for why that matters.
    (xp, yp, xc, yc, L, sigma, gam,
     xv, yv, Gv, rc, _Ve) = _solve_with_dsv(naca_csv, c, U, alpha, CL, CNv,
                                            tau_over_Tvl, Tvl)
    Gamma = 0.5*CL*U*c

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
    # SIGN. The dynamic-stall vortex is a roll-up of upper-surface boundary-layer
    # vorticity, so it rotates in the SAME sense as the bound circulation
    # (clockwise here). This term must therefore carry the same sign as the
    # bound surface sheet above. It does not need to supply the vortex lift -- Gamma is
    # already matched to the full UIBS C_L, which contains CNv -- so giving the
    # vortex its physical rotation costs nothing and buys the correct vorticity
    # field. It does NOT turn the flow over at the wall beneath the core, and
    # this comment used to say it did: the vortex induces 0.121*U upstream
    # there against a local 1.098*U (both published in metrics_*.csv), so the
    # comparison is a number now rather than an assertion.
    # The suction under the vortex comes from its low-pressure core (see
    # _core_pressure_deficit), not from accelerating the surface flow.
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
    # VORTICITY. The panel field is differenced numerically; the dynamic-stall
    # vortex's own vorticity is added in CLOSED FORM. Differencing it too was
    # wrong once the core became derived: at rc = 0.011c the core is under a
    # cell across, and a centred difference of a near-singular field returns a
    # four-lobed rosette of alternating sign -- a stencil artefact that reads as
    # a feature. Lamb-Oseen vorticity is omega = -Gamma_v/(pi rc^2) exp(-r^2/rc^2)
    # in this module's sign convention (the sheet and the vortex are both
    # clockwise), so it is sampled exactly at the nodes instead.
    dvx = np.gradient(v - dv, gx, axis=1); duy = np.gradient(u - du, gy, axis=0)
    vort = dvx - duy
    if Gv > 0.0:
        vort = vort - (Gv/(np.pi*rc*rc))*np.exp(-r2/(rc*rc))

    # Mask the aerofoil interior AND the one ring of cells touching it.
    # The bound sheet is a singular vortex sheet regularised over eps; a cell
    # that straddles the surface therefore samples the middle of the jump and
    # reports a spurious acceleration. It is measurable: at mid-chord on the
    # PRESSURE side, where the profile is monotonically falling towards the wall
    # (79.1 -> 80.4 -> 84.6 m/s), the last cell jumped to 111.7 m/s -- above the
    # free stream, on the side of the aerofoil that must be slower than it. The
    # field is simply not resolved within one cell of the wall, so it is not
    # published there. surface_cp is unaffected because it does not read this grid
    # at all: it evaluates the panel solution exactly at the control points, with
    # the source and vortex self-terms in closed form. (It used to probe 0.015c
    # off the wall, which is where this comment's "~7 cells out" came from; that
    # offset is gone.)
    from matplotlib.path import Path as MplPath
    poly = MplPath(np.column_stack([xp, yp]))
    inside = poly.contains_points(np.column_stack([X.ravel(), Y.ravel()])).reshape(X.shape)
    masked = inside.copy()                       # 8-connected dilation
    for _ring in range(NEAR_WALL_CELLS_MASKED):  # the count the config publishes
        _seed = masked.copy()
        for sx in (-1, 0, 1):
            for sy in (-1, 0, 1):
                if sx == 0 and sy == 0:
                    continue
                masked |= np.roll(np.roll(_seed, sy, axis=0), sx, axis=1)
    for arr in (u, v, speed, Cp, T_static, T_recovery, Mlocal, vort):
        arr[masked] = np.nan

    return dict(X=X, Y=Y, u=u, v=v, speed=speed, Cp=Cp, T_static=T_static,
                T_recovery=T_recovery, Mlocal=Mlocal, vort=vort,
                xp=xp, yp=yp, xv=xv, yv=yv, Gamma=Gamma, Gv=Gv,
                T0=T0, T_inf=T_inf)


def surface_cp(naca_csv, c, U, M, alpha_deg, CL, CNv, tau_over_Tvl,
               Tvl=TVL_DEFAULT):
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

    With all three corrected the closure error is -0.37% at alpha = 10 deg and
    falls monotonically from -0.28% at 2 deg to -0.62% at 19 deg, against -12.4%
    published previously. Both figures here read -0.18% and "inside -0.4% up to
    19 deg" until this audit: they were the pre-panel-tangent values, and they
    contradicted the -0.37% and -0.58% this module's own header already
    published for the same conditions. The header is the measured set.

    Returns (x/c at the control points, Cp, upper_mask).
    """
    alpha = np.radians(alpha_deg)
    (xp, yp, xc, yc, L, sigma, gam,
     xv, yv, Gv, rc, _Ve) = _solve_with_dsv(naca_csv, c, U, alpha, CL, CNv,
                                            tau_over_Tvl, Tvl)
    u, v = _surface_velocity(xp, yp, xc, yc, L, sigma, gam, U, alpha)
    du, dv, r2 = _dsv_velocity(xc, yc, xv, yv, Gv, rc)
    u += du; v += dv

    speed = np.hypot(u, v)
    Cp = (1.0 - (speed/U)**2
          + _core_pressure_deficit(np.sqrt(r2), abs(Gv), rc, U))
    # NOT clipped. The potential-flow leading-edge suction peak reaches Cp =
    # -15.4 at 17.5 deg. Neither this nor the field is clipped any more, and
    # reinstating the old -8 display clip HERE would be a real error rather than
    # a cosmetic one: measured, it takes the closure at 17.5 deg from -0.58 % to
    # -5.86 %. (At 10 deg the peak never reaches -8, so the same clip costs
    # nothing there -- which is why it looked harmless.) A real boundary layer
    # separates long before that peak is reached; that is a limitation of
    # reconstructing from potential flow, and it is stated rather than hidden.
    upper = yc >= 0
    return xc/c, Cp, upper


def dsv_vortex_state(naca_csv, c, U, M, alpha_deg, CL, CNv, tau_over_Tvl,
                     Tvl=TVL_DEFAULT):
    """The reconstructed dynamic-stall vortex, as published numbers.

    Returns a dict with, all derived and none of them read off a grid:
      Gamma_over_Uc  circulation / (U*c)
      rc_chords      Lamb-Oseen core radius, in chords
      peak_swirl_over_U   the vortex's own maximum swirl / U -- by construction
                          the edge speed of the shear layer that rolls it up
      induced_at_wall_over_U   the UPSTREAM velocity the vortex induces on the
                          wall directly beneath its centre, / U. This is the
                          quantity that decides whether the reconstruction shows
                          the flow reversal that characterises dynamic stall:
                          the vortex reverses the flow there when this exceeds
                          the local speed the reversal has to overcome, and it
                          falls off as 1/(2*pi*yv), so a strong vortex standing
                          well off the surface can still fail to reverse it.
      edge_speed_over_U    Ve/U at the vortex centre, vortex excluded.

    peak_swirl_over_U used to be computed with a coefficient of 0.7152. That is
    the fraction of the circulation ENCLOSED at the peak-swirl radius, not the
    peak-swirl coefficient, which is LAMB_OSEEN_PEAK = 0.638173; the published
    figure was 12 % high. Both are computed here rather than quoted.
    """
    alpha = np.radians(alpha_deg)
    (_xp, _yp, _xc, _yc, _L, _sig, _gam,
     xv, yv, Gv, rc, Ve) = _solve_with_dsv(naca_csv, c, U, alpha, CL, CNv,
                                           tau_over_Tvl, Tvl)
    if U <= 0.0:
        return dict(Gamma_over_Uc=float("nan"), rc_chords=float("nan"),
                    peak_swirl_over_U=float("nan"),
                    induced_at_wall_over_U=float("nan"),
                    edge_speed_over_U=float("nan"))
    return dict(Gamma_over_Uc=float(Gv/(U*c)),
                rc_chords=float(rc/c),
                peak_swirl_over_U=float(LAMB_OSEEN_PEAK*Gv/(2.0*np.pi*rc)/U)
                                  if Gv > 0 else 0.0,
                induced_at_wall_over_U=float(Gv/(2.0*np.pi*yv)/U) if yv > 0 else 0.0,
                edge_speed_over_U=float(Ve/U))


def dsv_core_cp(naca_csv, c, U, M, alpha_deg, CL, CNv, tau_over_Tvl,
                Tvl=TVL_DEFAULT):
    """Cp at the CENTRE of the reconstructed dynamic-stall vortex, evaluated
    exactly at that point.

    This is what metrics_*.csv publishes as Cp_DSV_core_min. It used to be read
    off the reconstructed field as the value at whichever grid node happened to
    lie nearest the vortex centre, which made a grid-independent quantity look
    as though it were not: the same instant of Case A gave -0.425 on a 110x85
    grid, -0.381 at 220x170, -0.358 at 440x340 and -0.368 at 880x680 -- drift
    and non-monotonicity that were entirely the sampling point moving, not the
    reconstruction changing. Evaluated here there is one number.

    The panel regularisation used is the panel-spacing one only; the field's
    eps2 also carries a grid term, which has no meaning away from a grid, and
    the core sits at least 0.10c off the surface where neither matters.

    Returns (xv, yv, Cp_at_the_centre).
    """
    alpha = np.radians(alpha_deg)
    (xp, yp, xc, yc, L, sigma, gam,
     xv, yv, Gv, rc, _Ve) = _solve_with_dsv(naca_csv, c, U, alpha, CL, CNv,
                                            tau_over_Tvl, Tvl)
    X = np.array([xv]); Y = np.array([yv])
    du_p, dv_p = _panel_velocity(X, Y, xc, yc, L, sigma, gam, (0.6*L.mean())**2)
    du_v, dv_v, r2 = _dsv_velocity(X, Y, xv, yv, Gv, rc)
    u = U*np.cos(alpha) + du_p + du_v
    v = U*np.sin(alpha) + dv_p + dv_v
    speed = np.hypot(u, v)
    cp = 1.0 - (speed/U)**2 + _core_pressure_deficit(np.sqrt(r2), abs(Gv), rc, U)
    return float(xv), float(yv), float(cp[0])


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


def surface_load_closure(naca_csv, c, U, M, alpha_deg, CL, CNv, tau_over_Tvl,
                         Tvl=TVL_DEFAULT):
    """Integrate the reconstructed surface Cp and compare with the C_L that the
    reconstruction was given.

    An invariant ONLY WITH THE VORTEX OFF (CNv = 0): a closed body carrying
    circulation Gamma = 0.5*CL*U*c and nothing else must return that C_L
    (Blasius), so the residual is then discretisation or a bug, never physics,
    and it converges to zero under panel refinement.

    With the dynamic-stall vortex present that is no longer true and must not be
    read as an error. A free vortex of circulation 1.79*U*c standing 0.16c off
    the body exerts a real force on it, so the residual is the vortex's induced
    lift; it converges under refinement to that force, not to zero. run_case
    therefore measures Cp_closure_error_pct with the vortex OFF and publishes
    the difference the vortex makes separately, as DSV_induced_lift_dCL.

    Also returns the trailing-edge pressure jump |Cp_upper - Cp_lower| there.
    That one is NOT an invariant and must not be read as one: it is zero only if
    the imposed circulation is the inviscid attached value, which is what
    kutta_reference_CL() computes. The reconstruction is deliberately handed the
    indicial C_L instead, so the jump is a measure of how far the modelled flow
    is from attached. It is linear in the imposed C_L, with a measured
    jump/|CL - CL_kutta| of 1.91-2.02 over alpha = 2-19 deg (falling
    monotonically with incidence), at the 160 panels the study uses.

    The integral is taken over the same panels the solution was built on, with
    Cp at the control points. It previously summed a mean of END-POINT Cp values
    against the end-point spacing, a quadrature the panel solution does not
    support.

    Both exist so the reconstruction's accuracy is a number the pipeline writes
    out, not a claim in a comment.
    """
    xp, yp = _airfoil_surface(naca_csv, c)
    _, cp_s, _ = surface_cp(naca_csv, c, U, M, alpha_deg, CL, CNv, tau_over_Tvl, Tvl)
    # traversal sense (shoelace): +1 counter-clockwise, so that n ds = (dy, -dx)
    sgn = 1.0 if 0.5*np.sum(xp[:-1]*yp[1:] - xp[1:]*yp[:-1]) > 0 else -1.0
    a = np.radians(alpha_deg)
    CN =  sgn*np.sum(cp_s*np.diff(xp))/c        # +(1/c) integral Cp dx
    CA = -sgn*np.sum(cp_s*np.diff(yp))/c        # -(1/c) integral Cp dy
    cl = CN*np.cos(a) - CA*np.sin(a)
    te_jump = float(abs(cp_s[0] - cp_s[-1]))    # control points either side of the TE
    return (float(cl), (100.0*(cl - CL)/CL if CL != 0 else float("nan")), te_jump)

