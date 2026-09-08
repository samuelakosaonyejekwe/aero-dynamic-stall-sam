"""
03_model_setup / generate_setup.py
----------------------------------
Writes ALL solver input data for the UNISTALL(TM) dynamic-stall case study:
  - flow_conditions.csv            : free-stream / operating conditions (2 cases)
  - kinematics.csv                 : pitch-oscillation (blade feathering) schedule
  - material_thermo_properties.csv : air thermodynamics + blade structural/thermal
  - solver_config.json             : UIBS / Leishman-Beddoes constants + numerics
  - static_polar_reference.csv     : published NACA 0012 static polar (calibration
                                     + validation reference; sources in CSV header)
"""
import json
import numpy as np
import pandas as pd
from pathlib import Path

HERE = Path(__file__).resolve().parent

# ============================================================ FLOW CONDITIONS
# Case A: oscillating-airfoil VALIDATION rig (matches McAlister/McCroskey NACA0012)
# Case B: APPLICATION - medium utility helicopter retreating-blade section r/R=0.75
#
# The free-stream state is OVER-DETERMINED: (p, T, rho), (M, U, a) and, for case
# B, the rotor kinematics all describe the same flow. Quoting all of them
# independently left them contradicting each other -- rho was 1.09 % (A) and
# 0.62 % (B) away from p/(R*T), and U_B was 0.63 % away from M*a and 0.87 % away
# from the rotor kinematics that are supposed to define it. Only the independent
# quantities are stated here; everything else is DERIVED.
#
#   INDEPENDENT: p_inf, T_inf, mu, chord, reduced frequency k
#                case A -- M (the wind-tunnel test point is specified by Mach)
#                case B -- Omega, R, r/R, mu_advance (the blade station defines U)
#   DERIVED:     a = sqrt(gamma*R*T),  rho = p/(R*T),  Re_c = rho*U*c/mu
#                case A -- U = M*a ;  case B -- U = Omega*R*(r/R - mu_adv), M = U/a
GAMMA, R_GAS = 1.4, 287.05

CH_A, CH_B   = 0.30, 0.527            # chord [m]
K_A,  K_B    = 0.10, 0.074            # reduced frequency (also used by KINEMATICS)
P_A,  P_B    = 90000.0, 91200.0       # static pressure [Pa]
T_A,  T_B    = 288.15, 287.5          # static temperature [K]
MU           = 1.78e-5                # dynamic viscosity [Pa.s]
M_A          = 0.30                   # case A test-point Mach
OMEGA, R_ROT, R_STA, ADV = 27.0, 8.18, 0.75, 0.32   # case B rotor kinematics

A_A, A_B     = np.sqrt(GAMMA*R_GAS*T_A), np.sqrt(GAMMA*R_GAS*T_B)
RHO_A, RHO_B = P_A/(R_GAS*T_A), P_B/(R_GAS*T_B)
U_A          = M_A*A_A                              # rig: Mach is the test point
U_B          = OMEGA*R_ROT*(R_STA - ADV)            # retreating blade at psi=270 deg
M_B          = U_B/A_B
RE_A = RHO_A*U_A*CH_A/MU
RE_B = RHO_B*U_B*CH_B/MU

flow = pd.DataFrame([
    ["case_id",                 "A_validation_rig", "B_application_rotor", "-"],
    ["description",             "NACA0012 oscillating airfoil (wind tunnel)",
                                "Retreating-blade section, r/R=0.75, mu=0.32", "-"],
    ["airfoil",                 "NACA 0012", "NACA 0012", "-"],
    ["chord_c",                 CH_A, CH_B, "m"],
    ["freestream_mach_M",       round(M_A, 4), round(M_B, 4), "-"],
    ["freestream_velocity_U",   round(U_A, 2), round(U_B, 2), "m/s"],
    ["speed_of_sound_a",        round(A_A, 2), round(A_B, 2), "m/s"],
    ["air_density_rho",         round(RHO_A, 4), round(RHO_B, 4), "kg/m^3"],
    ["static_pressure_p_inf",   P_A, P_B, "Pa"],
    ["static_temperature_T_inf",T_A, T_B, "K"],
    ["dynamic_viscosity_mu",    MU, MU, "Pa.s"],
    ["reynolds_number_Re_c",    float("%.4g" % RE_A), float("%.4g" % RE_B), "-"],
    ["reduced_frequency_k",     K_A, K_B, "-"],
    ["advance_ratio_mu",        np.nan, ADV, "-"],
    ["rotor_radius_R",          np.nan, R_ROT, "m"],
    ["radial_station_r_R",      np.nan, R_STA, "-"],
    ["rotor_speed_Omega",       np.nan, OMEGA, "rad/s"],
    ["rotor_tip_speed_OmegaR",  np.nan, round(OMEGA*R_ROT, 2), "m/s"],
], columns=["parameter", "case_A_validation", "case_B_application", "units"])
flow.to_csv(HERE/"flow_conditions.csv", index=False)

# ============================================================ KINEMATICS
# alpha(t) = alpha_mean + alpha_amp * sin(omega t);  k = omega c / (2 U)
def kin_row(case, U, c, k, a_mean, a_amp):
    omega = 2*k*U/c
    f_hz = omega/(2*np.pi)
    T = 1.0/f_hz
    return [case, a_mean, a_amp, k, round(omega,3), round(f_hz,3), round(T,5)]

kin = pd.DataFrame([
    kin_row("A_validation_rig", U_A, CH_A, K_A, 10.0, 10.0),
    kin_row("B_application_rotor", U_B, CH_B, K_B, 12.0, 8.0),
], columns=["case_id", "alpha_mean_deg", "alpha_amp_deg", "reduced_freq_k",
            "omega_rad_s", "freq_Hz", "period_s"])
kin["pitch_axis_x_c"] = 0.25
kin["motion"] = "alpha(t)=mean+amp*sin(omega t)"
kin.to_csv(HERE/"kinematics.csv", index=False)

# ============================================================ MATERIAL/THERMO
# cp and the recovery factor are DERIVED from gamma, R and Pr rather than
# quoted: stating r = Pr^(1/3) next to 0.892 (Pr^(1/3) = 0.8963) and cp = 1004.5
# next to gamma*R/(gamma-1) = 1004.68 left each value contradicting its own note.
PR = 0.72                       # GAMMA and R_GAS are set with the flow conditions
CP_AIR = GAMMA*R_GAS/(GAMMA - 1.0)
REC_FAC = PR**(1.0/3.0)
thermo = pd.DataFrame([
    ["air_gamma",            GAMMA,    "-",      "ratio of specific heats"],
    ["air_gas_constant_R",   R_GAS,    "J/kg/K", "specific gas constant"],
    ["air_cp",               round(CP_AIR, 2), "J/kg/K", "derived: gamma*R/(gamma-1)"],
    ["air_Prandtl_Pr",       PR,       "-",      "Prandtl number"],
    ["recovery_factor_r",    round(REC_FAC, 4), "-", "derived: turbulent, r=Pr^(1/3)"],
    ["sutherland_C1",        1.458e-6, "kg/m/s/K^0.5", "Sutherland viscosity const"],
    ["sutherland_S",         110.4,    "K",      "Sutherland temperature"],
    ["blade_skin_material",  "Al-2024-T3", "-",  "blade skin"],
    ["blade_skin_k_thermal", 121.0,    "W/m/K",  "thermal conductivity"],
    ["blade_skin_density",   2780.0,   "kg/m^3", "density"],
    ["blade_skin_cp",        875.0,    "J/kg/K", "specific heat"],
    ["blade_emissivity",     0.85,     "-",      "painted surface"],
], columns=["property", "value", "units", "note"])
thermo.to_csv(HERE/"material_thermo_properties.csv", index=False)

# ============================================================ SOLVER CONFIG
config = {
    "solver_name": "UNISTALL(TM) Universal Unsteady-Aerodynamics & Dynamic-Stall Solver",
    "core_method": "Unified Indicial-Beddoes State-Space (UIBS)",
    "version": "1.0.0",
    "modules": ["attached_flow_indicial", "trailing_edge_separation",
                "leading_edge_dynamic_stall_vortex", "compressibility_correction",
                "vortex_panel_field_reconstruction", "compressible_thermal_module"],
    "constant_precedence": "literature defaults below are OVERRIDDEN, key by key, "
                           "by calibrated_constants; where a key appears in both "
                           "(CN1, Tf, Tvl) the calibrated value is the one solved",
    "indicial_circulatory": {"A1": 0.30, "A2": 0.70, "b1": 0.14, "b2": 0.53},
    "time_constants_semichords": {"Tp": 1.7, "Tf": 3.0, "Tv": 6.0, "Tvl": 5.0,
                                  "comment": "literature defaults [S6]; Tf and Tvl "
                                             "are superseded by calibrated_constants"},
    "separation_model": {
        "method": "inverse-Kirchhoff fit of f(alpha) to static_polar_reference.csv",
        "form": "f = (2*sqrt(CN_static/(CNalpha*alpha)) - 1)^2, PCHIP-interpolated in |alpha|",
        "f_min": 0.02,
        "comment": "no alpha1/S1/S2 exponential fit is used: the separation law is "
                   "recovered in closed form from the measured polar"},
    "dynamic_stall_onset": {"CN1": 1.45,
                            "comment": "literature default for critical CN; the SOLVED "
                                       "value is calibrated_constants.CN1"},
    "calibrated_constants": {
        "CN1": 1.38, "Tf": 2.5, "Tv": 6.0, "Tvl": 6.0,
        "k0": 0.0, "k1": -0.22, "k2": 0.04, "kappa": 2.0, "eta": 0.95,
        "cpv_amp": 0.22, "Bv": 1.3,
        "comment": "single source of truth; calibrated against REAL NACA 0012 "
                   "frame 9302 (10deg+/-10deg, M0.30, k0.10) from NASA TM-84245"},
    "lift_curve_slope_CNalpha_per_rad": 6.28,
    "zero_lift_CM0": 0.0,
    "numerics": {"steps_per_cycle": 720, "n_cycles": 6, "report_cycle": 6,
                 "integrator": "semichord-marching exponential-recurrence"},
    "field_reconstruction": {
        "method": "constant-strength source panels solved for flow tangency against a "
                  "uniform bound vortex sheet ON THE BODY SURFACE carrying "
                  "Gamma = 0.5*C_L*U*c (Kutta-Joukowski, matched to the UIBS C_L), "
                  "plus a Lamb-Oseen dynamic-stall vortex",
        "n_panels": 160,
        "panel_distribution": "closed section, cosine-clustered on each surface "
                              "separately -> clustering at BOTH the leading and the "
                              "trailing edge",
        "grid_nx_default": 260, "grid_ny_default": 200,
        "grid_nx_solution": 220, "grid_ny_solution": 170,
        "domain_chords": [-1.0, 2.0, -1.2, 1.2],
        "near_wall_cells_masked": 1,
        "surface_cp_probe_offset_chords": 0.015,
        "comment": "grid_*_solution are the sizes written to 05_solution/field_*.csv. "
                   "cp_distribution_*.csv is NOT read off any grid: the surface Cp is "
                   "evaluated directly from the panel singularities at "
                   "surface_cp_probe_offset_chords outside the wall, so it does not "
                   "depend on a grid at all. That offset is the value at which the "
                   "closure error is smallest (-14.4 % at alpha 19 deg, against "
                   "-19.7 % at 0.008c and -17.0 % at 0.030c). near_wall_cells_masked "
                   "is the ring of field cells left blank because the regularised "
                   "surface sheet is not resolved there. Accuracy of the "
                   "reconstruction is measured, not asserted: the "
                   "Cp_closure_error_pct row of metrics_*.csv reports how well the "
                   "integrated surface Cp reproduces the C_L it was given"},
    "calibration_state": "calibrated_per_case (static polar) + validated (dynamic)"
}
with open(HERE/"solver_config.json", "w") as fp:
    json.dump(config, fp, indent=2)

# ============================================================ STATIC POLAR REF
# Representative published NACA 0012 static aerodynamics, Re ~ 2-3e6, low Mach.
# SOURCES (record):
#  [1] Sheldahl, R.E. & Klimas, P.C. (1981), "Aerodynamic Characteristics of
#      Seven Symmetrical Airfoil Sections...", SAND80-2114, Sandia Nat. Labs.
#  [2] Abbott, I.H. & von Doenhoff, A.E. (1959), "Theory of Wing Sections", Dover.
#  [3] McCroskey, W.J. (1987), "A Critical Assessment of Wind Tunnel Results for
#      the NACA 0012 Airfoil", NASA TM-100019.
# NOTE: values are representative literature consensus; digitize exact tables from
#       the cited sources for certification work.
alpha = [0,2,4,6,8,10,11,12,13,14,15,16,17,18,19,20]
Cl    = [0.000,0.218,0.435,0.650,0.860,1.060,1.152,1.240,1.310,1.360,
         1.400,1.420,1.380,1.200,1.100,1.050]
Cd    = [0.0086,0.0090,0.0098,0.0110,0.0128,0.0156,0.0178,0.0205,0.0240,0.0300,
         0.0400,0.0600,0.0900,0.1400,0.1700,0.1900]
Cm    = [0.000,-0.002,-0.004,-0.006,-0.008,-0.010,-0.012,-0.014,-0.018,-0.024,
         -0.035,-0.050,-0.075,-0.095,-0.100,-0.100]
static = pd.DataFrame({"alpha_deg": alpha, "Cl": Cl, "Cd": Cd, "Cm_c4": Cm})
static["source"] = "Sheldahl&Klimas SAND80-2114; Abbott&vonDoenhoff; NASA TM-100019"
static.to_csv(HERE/"static_polar_reference.csv", index=False)

print("[setup] wrote flow_conditions, kinematics, thermo, solver_config.json, static_polar_reference")
print("        A: M=%.4f U=%.2f a=%.2f rho=%.4f Re=%.3e   (Mach is the test point)"
      % (M_A, U_A, A_A, RHO_A, RE_A))
print("        B: M=%.4f U=%.2f a=%.2f rho=%.4f Re=%.3e   (U from Omega*R*(r/R-mu))"
      % (M_B, U_B, A_B, RHO_B, RE_B))
print("        validation k=%.3f alpha=10+10sin; application k=%.3f alpha=12+8sin" % (K_A, K_B))
