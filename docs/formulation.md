# Formulation of the load model

Author: Akosa Samuel Onyejekwe (independent)

Every equation and every constant that `unistall/attached_flow.py`,
`unistall/static_model.py` and `unistall/dsmodel.py` evaluate, in the order they are
evaluated. Source of the Leishman-Beddoes equations: Damiani & Hayman (2019),
NREL/TP-5000-66347, doi:10.2172/1576488 ("[DH19]"); their equation numbers are
given in brackets. Where the implementation departs from that source it says so.

## 1. Inputs

| Symbol | Meaning | Source |
|---|---|---|
| alpha(t) | incidence, alpha = alpha_m + alpha_1 sin(omega t) | case |
| U, c, M | speed, chord, Mach number; a_s = U/M | case |
| k | reduced frequency, omega c / (2U) | case |
| x_p | pitch axis, fraction of chord (0.25) | case |
| C_Nalpha(M) | lift-curve slope | static data, section 4 |
| alpha_0(M) | zero-lift incidence | static data, section 4 |
| C_N1(M) | normal force at static stall | static data, section 4 |
| f(alpha; M) | static separation point | static data, section 4 |
| C_M,s(alpha; M) | static moment about c/4 | static data, section 4 |

## 2. Attached flow

    (1)  ds = 2 U dt / c                                        [1.5b]
    (2)  beta^2 = 1 - M^2
    (3)  K_alpha,n = (alpha_n - alpha_{n-1}) / dt               [1.7]
    (4)  q_n = K_alpha,n c / U ;  K_q,n = (q_n - q_{n-1}) / dt  [1.7, 1.8]
    (5)  alpha_34,n = alpha_n + q_n (3/4 - x_p)

Circulatory part:

    (6)  X_1,n = X_1,n-1 exp(-b_1 beta^2 ds) + A_1 exp(-b_1 beta^2 ds/2) d(alpha_34)   [1.15]
    (7)  X_2,n = X_2,n-1 exp(-b_2 beta^2 ds) + A_2 exp(-b_2 beta^2 ds/2) d(alpha_34)   [1.15]
    (8)  alpha_e,n = alpha_34,n - alpha_0 - X_1,n - X_2,n       [1.14]
    (9)  C_N^c = C_Nalpha alpha_e                               [1.13]

Impulsive part:

    (10) T_I = c / a_s                                          [1.11c]
    (11) k_alpha = 1 / [(1 - M) + pi beta M^2 (A_1 b_1 + A_2 b_2)]
    (12) k_q     = 1 / [(1 - M) + 2 pi beta M^2 (A_1 b_1 + A_2 b_2)]
    (13) T_alpha = 0.75 k_alpha T_I ;  T_q = 0.75 k_q T_I       [1.10]
    (14) K'_alpha,n = K'_alpha,n-1 exp(-dt/T_alpha) + (K_alpha,n - K_alpha,n-1) exp(-dt/(2 T_alpha))   [1.18]
    (15) C_N^nc,alpha = (4 T_alpha / M)(K_alpha - K'_alpha)     [1.18]
    (16) K'_q,n = K'_q,n-1 exp(-dt/T_q) + (K_q,n - K_q,n-1) exp(-dt/(2 T_q))    [1.19]
    (17) C_N^nc,q = (T_q / M)(K_q - K'_q)                       [1.19, 1.20]
    (18) C_N^pot = C_N^c + C_N^nc,alpha + C_N^nc,q              [1.20]

Equation (11) is derived, not copied: the exact short-time response to a step
in incidence is C_N = (4/M)[1 - (1-M) s / (2M)] (Lomax et al.), and equating
its initial slope to that of (4/M) exp(-s/T') plus the circulatory function
gives T' = 2M k_alpha with k_alpha as in (11). [DH19] prints the same
expression for k_alpha and k_q.

Low Mach number. The constants above are for compressible flow and do not
tend to the incompressible limit: used alone at M = 0.05 they differ from
Theodorsen's solution by up to 6.4 % in lift amplitude, 12.9 % in moment
amplitude and 11 degrees in moment phase (`attached_checks.csv`). Below M = 0.20 the attached-flow loads are therefore blended with the
incompressible loads marched by the same recurrences (Jones' approximation to
the Wagner function, phi(s) = 1 - 0.165 exp(-0.0455 s) - 0.335 exp(-0.3 s),
with the apparent-mass terms of thin-aerofoil theory):

    (18a) w = (M - 0.10)/(0.20 - 0.10), limited to [0, 1]
    (18b) X = w X_compressible + (1 - w) X_incompressible   for every attached-flow load X

At and above M = 0.20, w = 1 and equations (1) to (25) are unchanged. The
values measured against Theodorsen's solution are in `attached_checks.csv`.

At the Mach numbers of the study the compressible branch is checked against
exact linear theory (`python3 -m unistall.check_indicial`,
`results/indicial_checks.csv`): the four starting values 4/M, -1/M, 1/M and
-7/(12M) are returned to better than 0.1 %. The lift after a step in
incidence follows the exact short-time result (4/M)[1 - (1-M)s/(2M)] only
near s = 0: the factor 0.75 in (13) and the single exponential leave it above
the exact value by the end of the interval 0 <= s <= 2M/(1+M), by the amount
in the file (about a tenth at Mach 0.30 and more at Mach 0.22). No
oscillatory compressible solution is used; at these Mach numbers the
frequency response rests on this check and on the measured loops below stall.

Unsteady moment terms about the quarter chord:

    (19) K'''_q,n = K'''_q,n-1 exp(-b_5 beta^2 ds) + A_5 exp(-b_5 beta^2 ds/2) (q_n - q_{n-1})   [1.26]
    (20) C_M^c,q = -(C_Nalpha / 16)(q - K'''_q)                 [1.25, without its factor c/U]
    (21) K_alphaM = (A_3 b_4 + A_4 b_3) / [b_3 b_4 (1 - M)]
    (22) C_M^nc,alpha = -(A_3 T_3 / M)(K_alpha - D_3) - (A_4 T_4 / M)(K_alpha - D_4),
         T_j = b_j K_alphaM T_I, D_j lagged as in (14) with T_j
    (23) k_mq = 7 / [15 (1 - M) + 1.5 C_Nalpha A_5 b_5 beta M^2] [1.29]
    (24) K''_q,n = K''_q,n-1 exp(-dt/(k_mq^2 T_I)) + (K_q,n - K_q,n-1) exp(-dt/(2 k_mq^2 T_I))    [1.29]
    (25) C_M^nc,q = -(7 k_mq^2 T_I / (12 M))(K_q - K''_q)       [1.29]

Equation (22) is the two-term form of Leishman and Beddoes. [DH19] gives the
simplification C_M^nc,alpha = -C_N^nc,alpha / 4 [1.27]; against Theodorsen's
solution at M = 0.05 that form over-predicts the moment amplitude by 49 %, so
it is not used.

## 3. Separated flow and dynamic stall

    (26) D_p,n = D_p,n-1 exp(-ds/T_p) + (C_N^pot_n - C_N^pot_{n-1}) exp(-ds/(2 T_p))   [1.35]
    (27) C'_N = C_N^pot - D_p                                   [1.35]
    (28) alpha_f = C'_N / C_Nalpha + alpha_0                    [1.34]
    (29) f' = f(alpha_f ; M)          table look-up, section 4
    (30) T_f = T_f0 / sigma_1                                   [1.37]
    (31) D_f,n = D_f,n-1 exp(-ds/T_f) + (f'_n - f'_{n-1}) exp(-ds/(2 T_f))             [1.36]
    (32) f'' = f' - D_f ,  F_MIN <= f'' <= 1                    [1.36]
    (33) C_N^fs = C_N^c ((1 + sqrt f'')/2)^2 + C_N^nc           [1.38]

Leading-edge separation and the vortex:

    (34) onset when |C'_N| first exceeds C_N1(M)                [1.46]
    (35) tau_v: 0 before onset; from the onset, tau_v advances by ds each step   [1.51]
    (36) T_sh = 2 (1 - f'') / St_sh ; tau_v restarts when tau_v >= T_VL + T_sh    [1.54]
    (37) C_V = C_N^c (1 - ((1 + sqrt f'')/2)^2)                 [1.49]
    (38) T_V = T_V0 / sigma_3                                   [1.48]
    (39) 0 < tau_v <= T_VL:  C_N^v,n = C_N^v,n-1 exp(-ds/T_V) + (C_V,n - C_V,n-1) exp(-ds/(2 T_V))   [1.47]
         otherwise:          C_N^v,n = C_N^v,n-1 exp(-ds/T_V)   [1.52]
    (40) C_N = C_N^fs + C_N^v                                   [1.53]

The onset crossing in (34) is located inside the time step by linear
interpolation of |C'_N|, and tau_v starts from the crossing. The vortex clock
is returned to zero when the vortex has left the chord (tau_v > T_VL), the
onset level is no longer exceeded and incidence is moving towards zero lift.
C_N^v is set to zero if it would oppose C_N^fs.

Multipliers ([DH19] section 2.6.2.3). With "separating" meaning f'' fell over
the previous step, "LE" meaning |C'_N| > C_N1, "on chord" meaning
0 < tau_v <= T_VL, and "away" meaning incidence is moving away from zero lift:

The tests are applied in this order, and the first that holds decides
("returning" means incidence is moving towards zero lift; at a turning point
the incidence is neither away nor returning):

    (41) sigma_1:
         if separating:
             2      returning
             1      otherwise, not LE
             2      otherwise, LE and f'' <= 0.7
             1.75   otherwise, LE and f'' > 0.7
         if not separating, the last of these that holds:
             1      (default)
             0.5    not LE
             0.25   vortex on chord
             0.75   away
    (42) sigma_3, the last of these that holds:
             1      (default)
             3      T_VL <= tau_v <= 2 T_VL and separating ; 4 if reattaching there
             2      vortex on chord and returning ; 1 vortex on chord and not returning
             4      no vortex on chord and returning

`tests/test_static_and_model.py` enumerates every combination of the flags
against this list.

Chord force, lift and drag:

    (43) C_C = eta C_N^c tan(alpha) sqrt f''                    [1.21, 1.32, with alpha for alpha_e + alpha_0]
    (44) when |C'_N| > C_N1:  C_C is multiplied by f''^k2,
         k2 = max(2 (|C'_N| - C_N1) + f'' - f', 0)              [1.56]
    (45) C_L = C_N cos alpha + C_C sin alpha                    [1.2]
    (46) C_D = C_N sin alpha - C_C cos alpha + C_D0             [1.2]

In (43) the leading-edge suction is resolved along the geometric incidence
alpha. [DH19] uses the lagged effective incidence alpha_e + alpha_0, which is
what potential theory gives; on the return stroke of a stalled cycle that form
keeps the suction pointing forward of where the measured force points, and the
drag falls below zero for part of the cycle, which the measurements do not
show (their lowest value is -0.009). With alpha the drag meets the limit set
in `targets.json` on every calibration loop (`separated_checks.csv`) and the
drag error is slightly lower; on the loops that cross zero lift the chord
force and the drag are listed in `results/limit_zero_lift.csv`, including one
held-out loop on which the drag is below zero for just over 0.05 of the cycle.

Normal force and moment are odd about the zero-lift incidence alpha_0 (the
moment about its static value there), to rounding. Lift and drag are resolved
with the incidence of the experiment's own scale in (43), (45) and (46),
because that is how the measured lift and drag were resolved; a loop and its
mirror image about alpha_0 therefore differ in drag, by the amounts in
`results/limit_mirror_pair.csv` (about 0.012 at most, for alpha_0 = 0.23
degrees).

The slow-sweep limit of (46) is not the measured static drag. The static data
are total drag from a wake survey; the oscillating loops are pressure drag,
and eta is fitted to them. `results/limit_static_drag.csv` sets the two side
by side: below stall the model's drag is higher than the static total drag,
and through stall the two differ by up to about 0.1. C_D0 is the minimum of
the static total drag and is kept as a constant; the file gives the
comparison without it as well.
The lagged form is kept as the option `cc_geometric = 0`.

In (44) eta is kept on both sides of the threshold; as printed, [1.56a] drops
it above the threshold, which would put a step in the chord force at onset.

Moment about the quarter chord:

    (47) alpha'_f = alpha_f lagged as in (31) with time constant c_m T_f0   [1.43]
    (48) x_v = x_cp,v (1 - cos(pi min(tau_v / T_VL, 1)))        [1.57]
    (49) C_M = C_M,s(alpha'_f ; M) + g (K_f C_M^c,q + C_M^nc,alpha + C_M^nc,q) - x_v C_N^v   [1.59]
         K_f = ((1 + sqrt f'')/2)^2
         g   = 1 - (1 - g_M) w ,  w from (18a)

The circulatory pitch-rate moment is carried by the circulation, so in (49)
it is reduced by separation in the same proportion as the circulatory normal
force in (33). K_f = 1 in attached flow.

The factor g_M is empirical. Before stall onset the model gives a wider moment
loop than is measured; g_M is the one number that minimises the moment error
over the attached part of every calibration loop, that is, the measured moment
points that the model with literature stall constants places before onset
(`results/attached_moment_factor.json` gives the value, its 95 % bootstrap
interval over loops, and the loops and points used). It is fixed before the
stall constants are fitted. Its strength follows the blend weight: the blended
unsteady moment is multiplied by 1 - (1 - g_M) w, which is g_M at and above
Mach 0.20 and 1 at and below Mach 0.10, so the incompressible limit checked
against Theodorsen's solution is unchanged. Between the two limits it is not
g_M on the compressible part alone; no measured loop tests that range.

What g_M corrects. The first term of (49), the static table read at the lagged
incidence, opens a loop of its own, and on the loops that peak below static
stall it supplies part of the excess damping (the result file gives that part
loop by loop). That lag is largely physical: the attached part of the static
moment is the normal force acting ahead of the quarter chord, and moving it
from the lagged incidence to the circulatory incidence changes the damping
little. What remains is an excess of the inviscid unsteady moment over
measurement in attached flow at Mach 0.30, which g_M removes. It rests on one
Mach number and is a correction of this model, not a property of the flow;
`results/uncertainty_summary.csv` carries its interval into the held-out
measures.
With g_M = 1, (49) is the form of [DH19].

## 4. Static inputs

Two Mach stations, both from the calibration set or the NASA report's static
figures (`unistall/static_model.py`):

| | M = 0.215 | M = 0.302 |
|---|---|---|
| slope, per degree | 0.1180 | 0.1170 |
| alpha_0, degrees | -0.01 | +0.24 |
| static-stall incidence alpha_ss, degrees | 16.2 | 13.6 |
| C_N1 | 1.598 | 1.336 |
| source of slope, alpha_0, C_N1 | quasi-steady sweep frame 13308 | quasi-steady sweep frame 12102 |
| source of f and C_M,s | sweeps 13308 and, above 17 degrees, 12300 (Mach 0.204) | NASA TM-84245 Vol. 1 Fig. 16(a), Fig. 9, Table 7; Vol. 2 Fig. 4 |

    (50) f = (2 sqrt(C_N,s / (C_Nalpha,s (alpha - alpha_0,s))) - 1)^2 ,  F_MIN <= f <= 1,
         f = 1 for |alpha - alpha_0,s| < 2 degrees

with the slope and zero-lift incidence of the data set the static C_N,s came
from. Past static stall the measurements are not single-valued: two branches
at M = 0.302, an up-sweep and a down-sweep at M = 0.215. The model follows the
attached curve up to static stall and the more separated state beyond it (the
lower lift branch at M = 0.302, the down-sweep at M = 0.215) at both stations.
This choice was made on the calibration loops, where it gives a lower cost
than the less separated state or the mean of the two.
Beyond 20 degrees at M = 0.302, where static drag was not measured,
C_N,s = C_L,s / cos(alpha).

Between the stations, with w = (M - 0.215)/(0.302 - 0.215) limited to [0, 1]:

    (51) C_Nalpha, alpha_0, alpha_ss, C_N1 :  (1 - w) X_a + w X_b
    (52) f(alpha; M) = (1 - w) f_a(alpha alpha_ss,a / alpha_ss) + w f_b(alpha alpha_ss,b / alpha_ss)
    (53) C_M,s(alpha; M) likewise

## 4a. The fitted separation law of the reference model

The reference model (`unistall/reference_lb.py`) differs from the model above
in (29) only. In place of the table it uses the three-constant exponential law
of Leishman and Beddoes, [DH19] Eq. (1.33), with x = |alpha - alpha_0|:

    (53a) f = 1 - 0.3 exp((x - alpha_1)/S_1)        x <= alpha_1
          f = 0.04 + 0.66 exp((alpha_1 - x)/S_2)    x >  alpha_1

alpha_1, S_1 and S_2 are fitted by bounded least squares to the tabulated f of
each Mach station on a grid of 0.25 degrees from 2 degrees up, and
interpolated in Mach number as in (51) and (52)
(`results/reference_exponential_fit.csv` gives the constants and the residual
of each fit). Everything else is shared: slope, alpha_0, C_N1, the static
moment table, the constants and the march. The reference model is therefore
the present model with one input replaced, not an independent implementation
of the published model, and it is fitted to the table, not to the static data
directly.

## 5. Constants

| Constant | Value | Status | Source |
|---|---|---|---|
| A_1, A_2 | 0.30, 0.70 | fixed | Leishman & Beddoes, as given in [DH19] |
| b_1, b_2 | 0.14, 0.53 | fixed | as above |
| A_3, A_4 | 1.5, -0.5 | fixed | Leishman & Beddoes impulsive moment |
| b_3, b_4 | 0.25, 0.1 | fixed | as above |
| A_5, b_5 | 1.0, 0.5 | fixed | Leishman & Beddoes ([DH19] prints b_5 = 5) |
| factor 0.75 in (13) | 0.75 | fixed | [DH19] Eq. (1.10) |
| x_cp,v (`xcp_v`) | 0.20 | fixed | [DH19] Eq. (1.57) |
| St_sh (`St_sh`) | 0.19 | fixed | Strouhal number of repeated shedding, [DH19] Eq. (1.54) |
| c_m (`cm_lag`) | 0.10 | fixed | [DH19] Eq. (1.43) |
| C_D0 (`CD0`) | 0.0072 | fixed | NASA TM-84245 Vol. 1 Table 8 |
| F_MIN | 0.02 | fixed | floor on the separation point |
| chord-force offset (`cc_offset`) | 0 | fixed | 0 gives (43); 0.2 gives Gonzalez's form sqrt f'' - 0.2, [DH19] Eq. (1.40), kept as an option and not used |
| blend limits in (18a) | 0.10, 0.20 | fixed | below the lowest Mach number of the static data |
| g_M (`cm_unsteady`) | see `results/attached_moment_factor.json` | empirical, from loops below static stall | `unistall/calibrate.py --attached` |
| onset-level factor (`cn1_factor`) | 1 | fixed | multiplies C_N1 in (34); 1 is the static value (`results/onset_level_check.csv`) |
| `cmq_kirchhoff` | 1 | fixed | 1 applies K_f to the circulatory pitch-rate moment in (49); 0 does not |
| `sigma_rules` | 1 | fixed | 1 uses (41) and (42); 0 sets both multipliers to 1 (T_V halved once the vortex has left the chord) |
| `cc_floor` | 1 | fixed | 1 holds the leading-edge suction C_N^c tan(alpha) in (43) at zero where the formula would make it negative (the lagged circulatory force and the incidence of opposite sign near zero lift); 0 leaves it as the formula gives it |
| `cc_geometric` | 1 | fixed | 1 uses alpha in (43); 0 uses alpha_e + alpha_0 as [DH19] |
| alpha_1, S_1, S_2 of (53a) | see `results/reference_exponential_fit.csv` | fitted to the table, reference model only | `unistall/reference_lb.py` |
| drop in C_M that marks moment stall, (53g) | 0.05 | fixed | `unistall/metrics.py` |
| T_p, T_f0, T_V0, T_VL, eta (`Tp`, `Tf0`, `Tv0`, `Tvl`, `eta`) | see `calibrated_constants.json` | fitted or literature | `unistall/calibrate.py`; the file lists which are fitted and which, if any, were returned to literature values (1.7, 3.0, 6.0, 7.0, 0.95) |

## 5a. Calibration cost

For each calibration loop j and each coefficient X in {C_L, C_M, C_D}:

    (53b) e_X,j = sqrt( mean_i (X_model(alpha_i) - X_i)^2 ) / (max_i X_i - min_i X_i)

where (alpha_i, X_i) are the measured points of the loop and the model is read
at the same incidence on the same stroke (up-stroke points against the model's
up-stroke, down-stroke against down-stroke; the measured record is split at
its two turning points). The cost is

    (53c) J = sum_j ( e_CL,j^2 + e_CM,j^2 + e_CD,j^2 )

minimised over the fitted constants by bounded least squares on the vector of
the e (`unistall/calibrate.py`), from 8 seeded starts at a coarse march, each
then refined at the resolution of (54) and (55). Which constants stay fitted
is decided by five-fold cross-validation inside the calibration loops under
the rule recorded in `data/targets.json`.

## 5b. Measures reported

    (53d) loop error of X                 e_X of (53b), averaged over loops
    (53e) peak-lift error                 100 (max C_L,model - max C_L,measured) / max C_L,measured
    (53f) incidence of maximum lift       alpha at max C_L, model minus measured
    (53g) incidence of moment stall       the up-stroke incidence at which C_M first falls 0.05 below its
                                          attached level (the median of C_M over the first third of the
                                          up-stroke range), by linear interpolation; model minus measured;
                                          defined only where both show one
    (53h) cycle damping                   Xi = -(1/(pi alpha_1^2)) closed integral of C_M d(alpha), alpha in radians;
                                          taken from a measured record only if each stroke spans at least
                                          0.8 of its incidence range
    (53i) damping sign agreement          fraction of loops on which Xi of model and measurement have the same sign

Means over loops carry a 95 % interval from 10000 bootstrap resamples of the
loops (`unistall/metrics.py`). The measured points are a tracing of the
printed loops, about 38 points a loop, so (53f) and (53g) are resolved no
finer than the spacing of those points.

## 6. Numerics

There is no spatial mesh: the model is a set of first-order lags marched in
time. The linear lags use exponential recurrences with the mid-step factor
exp(-ds/2T) on the increment of the input (exact for the decay, second order
in the step for a linearly varying input), so the march has no stability
limit, but the switching of the time-constant multipliers, the
vortex leaving the chord and repeated shedding are resolved only to one step,
and convergence with the step is first order. The lags are measured in
semichords of travel, so the step is set in semichords and not as a fraction
of the cycle:

    (54) steps per cycle = 2 pi / (k ds_max), rounded up to a multiple of 360 and at least 720,  ds_max = 0.02
    (55) cycles marched  = the number needed for 200 semichords of travel before the reported cycle, between 3 and 8

`results/convergence_summary.csv` gives the change in every loop error when
the step count is multiplied by four and when four cycles are added, and
`results/convergence_extrapolated.csv` the error of the default step against
the value at zero step, by Richardson extrapolation from three resolutions
(first order assumed where a switching event makes the three values
non-monotone). The step is the one at which that error is inside the target. The
search stage of the calibration marches at ds = 0.32; the calibrated constants
are then refined at the resolution of (54) and (55).

The first step returns the loads of the attached state at the first incidence.
The tables of f and the static moment are even and odd about the measured
zero-lift incidence, which for this symmetric section is an offset of the
incidence scale.

## 6a. State-space form

The same model is also written as first-order differential equations
(`unistall/statespace.py`). With s the travel in semichords, a prime d/ds, and
each lag written on the quantity it delays, the state vector and its
equations are

    (56) W_i' = b_i beta^2 (alpha_34 - W_i),  i = 1, 2            X_i = A_i (alpha_34 - W_i)
    (57) W_5' = b_5 beta^2 (q - W_5)                              C_M^c,q = -(C_Nalpha/16)(q - A_5 (q - W_5))
    (58) Z_a' = (c/2U)(K_alpha - Z_a)/T_alpha ,  Z_q' = (c/2U)(K_q - Z_q)/T_q
         C_N^nc = (4 T_alpha/M) Z_a + (T_q/M) Z_q
    (59) Z_3' = (c/2U)(K_alpha - Z_3)/T_3 ,  Z_4' likewise with T_4 ,  Z_m' = (c/2U)(K_q - Z_m)/(k_mq^2 T_I)
         C_M^nc = -(A_3 T_3/M) Z_3 - (A_4 T_4/M) Z_4 - (7 k_mq^2 T_I/(12 M)) Z_m
    (60) C'_N ' = (C_N^pot - C'_N)/T_p
    (61) f'' ' = (f(alpha_f ; M) - f'')/T_f
    (62) alpha_m' = (alpha_f - alpha_m)/(c_m T_f0)                static moment read at alpha_m
    (63) C_N^v ' = -C_N^v/T_V + [vortex on chord] d/ds( C_N^c (1 - K_f) )

Each recurrence of sections 2 and 3 approximates the solution of the matching
equation over one step with a linearly varying input to second order in the
step, so the two forms are
the same model. Equations (56) to (63) are integrated by the classical
fourth-order Runge-Kutta method with the step of (54). T_f and T_V take the
values of (41) and (42), evaluated at the start of each step and held over
it, and the vortex clock of (35) and (36) is advanced at the end of the step;
across a switch the two forms therefore agree to first order in the step.
`results/statespace_summary.csv` gives the measured difference on the
calibration loops. Below M = 0.20 the blend of (18a) is carried in the same
way: two more states, the three-quarter-chord incidence lagged at the two
rates of Jones' approximation, give the incompressible circulatory lift, and
the apparent-mass terms are written directly in alpha' and alpha''. Where a
lag is shorter than the step (the impulsive lags shorten with Mach number)
the step is divided into equal Runge-Kutta sub-steps.

Coupled to the torsional section of `unistall/structural.py`, the pitch
angle and its rate join the state vector and the equations are integrated as
one system.

## 6b. Switching rules: what they do on the measured loops

`python3 -m unistall.check_limits` counts two rules on every usable loop
(`results/limit_vortex_rules.csv`). The vortex clock restarts under sustained
onset, (1.54), on most stalled loops, and many times in a cycle on the slowest
ones; `n_sheddings` counts every restart. The vortex normal force is set to
zero in a step where it would oppose the separated force; the file gives the
number of such steps and the largest force removed, which is small beside the
loop errors.

## 7. Range

Static data are at Mach 0.215 and 0.302. Every calibration loop and all but
four of the held-out loops are at Mach 0.285 or above, so the comparison with
measurement is in effect at Mach 0.29 to 0.30; the four held-out loops between
Mach 0.22 and 0.285 probe the interpolation in Mach number and are predicted
worse (`results/validation_summary.csv`). Reduced frequency runs up to 0.20
and peak incidence up to 25 degrees on the NACA 0012 section, in sinusoidal
pitch only; other motions are not compared with measurement. `dsmodel.solve`
reports any condition outside the range of the static data, that reduced
frequency or that incidence in `meta["outside_compared_range"]` and raises an
error for it when called with `strict=True`.
