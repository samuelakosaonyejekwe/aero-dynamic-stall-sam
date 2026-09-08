"""
04_solver / run_case.py
-----------------------
Drives the UNISTALL(TM) solver for both case-study configurations and writes
all solution data to 05_solution/.  Outputs:
  time_history_<case>.csv      per-step unsteady loads & states (last cycle)
  cp_distribution_<case>.csv   surface Cp(x/c) at several phase angles
  field_<case>_ph<deg>.csv     reconstructed 2D fields at key phases
  model_static_polar.csv       quasi-steady model polar (for validation)
  metrics_<case>.csv           engineering scalar metrics (deterministic)
  runtime_environment.csv      the machine, and the CPU time the march took on it
  convergence/residuals_<case>.csv   cycle-to-cycle convergence
"""
import sys, json, time, platform
import numpy as np
import pandas as pd
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
import unistall_solver as us

SETUP = ROOT/"03_model_setup"
SOL   = ROOT/"05_solution"
GEO   = ROOT/"01_geometry"/"naca0012_coordinates.csv"
(SOL/"convergence").mkdir(parents=True, exist_ok=True)

cfg = json.load(open(SETUP/"solver_config.json"))
CNALPHA = cfg["lift_curve_slope_CNalpha_per_rad"]
NPC  = cfg["numerics"]["steps_per_cycle"]
NCYC = cfg["numerics"]["n_cycles"]
consts = dict(**cfg["indicial_circulatory"], **cfg["time_constants_semichords"])
# ---- calibrated constants (single source of truth; calibrated to REAL NACA0012
#      frame 9302 from NASA TM-84245 — see 06_postprocessing/validation) ----
consts.update({k: v for k, v in cfg["calibrated_constants"].items() if k != "comment"})

# ---- calibrate static separation to published polar ----
stat = pd.read_csv(SETUP/"static_polar_reference.csv")
f_static = us.calibrate_separation(stat["alpha_deg"], stat["Cl"], stat["Cd"], CNALPHA)

# ---- closure constants: taken from the config / reference polar, never inlined,
#      so the quasi-steady polar below and the dynamic march use the same values
ETA = consts["eta"]                                        # LE-suction efficiency
CD0 = float(stat.loc[stat["alpha_deg"].abs().idxmin(), "Cd"])   # zero-lift drag

# ---- model quasi-steady polar (validation) ----
aq = np.linspace(0, 22, 89)
fq = f_static(aq)
CN_qs = CNALPHA*((1+np.sqrt(fq))/2)**2*np.radians(aq)
CC_qs = ETA*CNALPHA*np.radians(aq)**2*np.sqrt(fq)
CL_qs = CN_qs*np.cos(np.radians(aq)) + CC_qs*np.sin(np.radians(aq))
CD_qs = CN_qs*np.sin(np.radians(aq)) - CC_qs*np.cos(np.radians(aq)) + CD0
pd.DataFrame({"alpha_deg": aq.round(3), "Cl_model": CL_qs.round(4),
              "Cd_model": CD_qs.round(4), "f_sep": fq.round(4)}
             ).to_csv(SOL/"model_static_polar.csv", index=False)

# ---- air thermodynamics: READ from 03_model_setup. These were previously only
#      rendered into the report while the solver used its own hardcoded defaults,
#      so editing the published "solver input" changed nothing. ----
_th = pd.read_csv(SETUP/"material_thermo_properties.csv").set_index("property")["value"]
THERMO = dict(gamma=float(_th["air_gamma"]), cp=float(_th["air_cp"]),
              recovery=float(_th["recovery_factor_r"]),
              R_gas=float(_th["air_gas_constant_R"]))

# ---- case definitions: READ from 03_model_setup, never restated here, so the
#      documented inputs and the inputs actually solved can never drift apart ----
flow = pd.read_csv(SETUP/"flow_conditions.csv").set_index("parameter")
kin  = pd.read_csv(SETUP/"kinematics.csv").set_index("case_id")

# Sampling incidences are capped just inside the top of the static polar the
# separation law was calibrated on. Past that the model is extrapolating (see
# calibrate_separation's decay branch), so a "peak incidence" field sampled there
# would be showing an extrapolated state. Two different magic caps, 19.5 and
# 19.0, used to sit inline with no reason given; both are now this one number.
ALPHA_SAMPLE_MAX = float(stat["alpha_deg"].max()) - 1.0     # 20 deg polar -> 19 deg

def _case(col):
    """Assemble one case from flow_conditions.csv + kinematics.csv."""
    f = flow[col]
    k_row = kin.loc[f["case_id"]]              # flow_conditions carries the kinematics key
    return dict(a_mean=float(k_row["alpha_mean_deg"]), a_amp=float(k_row["alpha_amp_deg"]),
                k=float(k_row["reduced_freq_k"]), M=float(f["freestream_mach_M"]),
                c=float(f["chord_c"]), U=float(f["freestream_velocity_U"]),
                T_inf=float(f["static_temperature_T_inf"]))

CASES = {"A_validation":  _case("case_A_validation"),
         "B_application": _case("case_B_application")}

def phase_index(out, target_alpha_deg, upstroke=True):
    """index of nearest matching alpha on up/down stroke."""
    a = out["alpha_deg"]; ad = out["alpha_dot"]
    mask = (ad > 0) if upstroke else (ad < 0)
    idx = np.where(mask)[0]
    j = idx[np.argmin(np.abs(a[idx]-target_alpha_deg))]
    return j

# One short throwaway march before any timing. The first call into the solver
# pays interpreter and SciPy/NumPy warm-up, which made Case A read 1.23 s CPU
# against Case B's 0.60 s for identical work -- a 2x difference that is entirely
# cold start and would have been published as if it were a property of the case.
us.solve_dynamic_stall(10.0, 10.0, 0.10, 0.30, 0.30, 102.0, f_static,
                       CNalpha=CNALPHA, consts=consts, n_per_cycle=60, n_cycles=1)

summary_rows = []
solve_cpu = {}
for name, C in CASES.items():
    # ---- CPU time of the march. NOT wall time: on this machine the same code
    #      measured 0.57 s idle and 3.96 s under load, and it was that
    #      load-contaminated wall time that made the report's "< 1 s per
    #      operating point" look false when it is true. CPU time is the
    #      load-independent measure. It is still machine-dependent, so it is
    #      written to runtime_environment.csv beside the machine that produced
    #      it, NOT to metrics_*.csv, which must stay bit-reproducible. ----
    _t0 = time.process_time()
    out = us.solve_dynamic_stall(C["a_mean"], C["a_amp"], C["k"], C["M"], C["c"], C["U"],
                                 f_static, CNalpha=CNALPHA, CD0=CD0,
                                 CM0=cfg["zero_lift_CM0"], consts=consts,
                                 n_per_cycle=NPC, n_cycles=NCYC)
    solve_s = time.process_time() - _t0
    solve_cpu[name] = solve_s
    # ---- time history CSV ----
    th = pd.DataFrame({
        "time_s": out["t"], "phase_deg": out["phase_deg"], "alpha_deg": out["alpha_deg"],
        "alpha_dot_rad_s": out["alpha_dot"], "alpha_eff_deg": out["alpha_eff_deg"],
        "CL": out["CL"], "CD": out["CD"], "CM_c4": out["CM"], "CN": out["CN"],
        "CC": out["CC"], "CN_prime": out["CNp"], "CN_vortex": out["CNv"],
        "CN_attached": out["CNf"], "f_separation": out["f_sep"],
        "vortex_active": out["vortex_active"], "tau_v_semichords": out["tau_v"],
    }).round(6)
    # time_s needs more than 6 decimals. The step is 1.28e-4 s (Case A), so
    # rounding the time base to 1 us leaves a 0.78 % spread on a step that is
    # uniform by construction -- enough that differentiating the published
    # history gives 0.2 % noise in alpha_dot that is not in the solution. The
    # loads are coefficients of order 1 and 6 decimals is ample for them.
    th["time_s"] = out["t"].round(12)
    th.to_csv(SOL/f"time_history_{name}.csv", index=False)

    # ---- engineering metrics ----
    a = out["alpha_deg"]; up = out["alpha_dot"] > 0
    CLmax = out["CL"].max(); iCL = out["CL"].argmax()
    CMmin = out["CM"].min(); iCM = out["CM"].argmin()
    CDmax = out["CD"].max()
    # ---- CD_min and the cycle mean. The instantaneous drag goes NEGATIVE for
    #      about a third of the cycle: on the downstroke the effective incidence
    #      lags 1-2 deg above the geometric one, so the leading-edge suction term
    #      CC*cos(alpha) outweighs CN*sin(alpha). That is the unsteady-thrust
    #      behaviour a chord-force model is meant to reproduce, and the CYCLE
    #      MEAN stays positive (no net propulsion) -- but reporting only CD_max
    #      left a reader with a plainly negative C_D loop and nothing to check
    #      it against. Both bounds and the mean are published. ----
    CDmin = out["CD"].min()
    CDmean = float(us._trapz(out["CD"], out["t"])/(out["t"][-1]-out["t"][0]))
    xi, xi_hat = us.aerodynamic_damping(a, out["CM"], normalise=True)
    verdict = us.damping_verdict(xi_hat)                # three-way, tolerance-banded
    # hysteresis loop areas
    loopCL = np.abs(us._trapz(out["CL"], np.radians(a)))
    # ---- dynamic-stall onset. The model sheds the leading-edge vortex when the
    #      LAGGED normal force CN' reaches CN1 on the upstroke, so that is the
    #      onset. This was computed from CN = CNf + CNv instead, a different
    #      quantity, which reported 11.22 deg for Case A when the vortex actually
    #      starts at 12.50 deg. Read it straight off the solver's own vortex
    #      switch so the metric and the criterion can never disagree again. ----
    if np.any(out["vortex_active"] > 0):
        onset = float(a[int(np.argmax(out["vortex_active"] > 0))])
    else:
        onset = np.nan
    CL_static_max = CL_qs.max()
    # ---- reconstruction closure: does the reconstructed surface Cp integrate
    #      back to the C_L it was handed? Published so the field-reconstruction
    #      accuracy is a measured number rather than a claim in a docstring. ----
    _j = int(np.argmax(out["CL"]))
    _clcp, cp_closure_pct, _ = us.surface_load_closure(
        GEO, C["c"], C["U"], C["M"], out["alpha_deg"][_j], out["CL"][_j],
        out["CNv"][_j], out["tau_v"][_j]/consts["Tvl"])
    # ---- Kutta residual. Reported as the WORST of the phases written to
    #      cp_distribution_*.csv, not just the one at peak lift.
    #
    #      This is NOT an error that can be driven to zero, and the published
    #      CL_kutta_inviscid row is what makes that checkable. The trailing-edge
    #      Cp jump is linear in the imposed C_L and vanishes exactly at the
    #      inviscid attached circulation; the reconstruction is instead handed
    #      the indicial C_L, which during dynamic stall departs from that value
    #      on purpose. The jump is therefore a measure of how far the modelled
    #      flow is from attached. A reader can verify the claim by imposing
    #      CL_kutta_inviscid, which drives the jump to ~1e-3. ----
    cp_te_jump = 0.0
    for _tgt, _ups in [(C["a_mean"], True),
                       (C["a_mean"]+C["a_amp"]*0.7, True),
                       (min(C["a_mean"]+C["a_amp"], ALPHA_SAMPLE_MAX), True),
                       (C["a_mean"], False)]:
        _k = phase_index(out, _tgt, upstroke=_ups)
        _, _, _tj = us.surface_load_closure(
            GEO, C["c"], C["U"], C["M"], out["alpha_deg"][_k], out["CL"][_k],
            out["CNv"][_k], out["tau_v"][_k]/consts["Tvl"])
        cp_te_jump = max(cp_te_jump, _tj)
    cl_kutta = us.kutta_reference_CL(GEO, C["c"], C["U"], C["M"],
                                     float(out["alpha_deg"][_j]))
    # ---- depth of the reconstructed dynamic-stall-vortex core. Published so
    #      that "the core is diffuse" is a number a reader can check rather than
    #      an adjective: it is the Cp at the vortex centre at the instant of
    #      peak vortex loading. Neither DSV constant can be calibrated from the
    #      data this study ships (integrated cl/cd/cm only), so the shallowness
    #      is reported, not tuned away. ----
    _v = int(np.argmax(out["CNv"]))
    _F = us.reconstruct_field(GEO, C["c"], C["U"], C["M"], out["alpha_deg"][_v],
                              out["CL"][_v], out["CNv"][_v],
                              out["tau_v"][_v]/consts["Tvl"], T_inf=C["T_inf"])
    _d = np.hypot(_F["X"]-_F["xv"], _F["Y"]-_F["yv"])
    cp_dsv_core = float(_F["Cp"][np.unravel_index(np.nanargmin(_d), _d.shape)])
    met = pd.DataFrame({
        "metric": ["CL_max_dynamic", "alpha_at_CLmax_deg", "CL_max_static",
                   "dynamic_overshoot_ratio", "CM_min(c/4)", "alpha_at_CMmin_deg",
                   "CD_max", "CD_min", "CD_cycle_mean",
                   "stall_onset_alpha_deg", "aero_damping_Xi",
                   "aero_damping_Xi_normalised", "damping_neutral_band",
                   "stall_flutter_risk",
                   "CL_hysteresis_loop_area", "Cp_closure_error_pct",
                   "Cp_TE_jump_max_over_phases", "CL_kutta_inviscid",
                   "Cp_DSV_core_min",
                   "reduced_frequency_k", "mach_M", "mean_alpha_deg", "amp_alpha_deg"],
        "value": [round(CLmax,3), round(a[iCL],2), round(CL_static_max,3),
                  round(CLmax/CL_static_max,3), round(CMmin,3), round(a[iCM],2),
                  round(CDmax,3), round(CDmin,3), round(CDmean,4),
                  round(float(onset),2), round(xi,5),
                  round(xi_hat,4), us.DAMPING_TOL, verdict,
                  round(loopCL,4), round(cp_closure_pct,1), round(cp_te_jump,3),
                  round(cl_kutta,3), round(cp_dsv_core,3),
                  C["k"], C["M"], C["a_mean"], C["a_amp"]],
    })
    met.to_csv(SOL/f"metrics_{name}.csv", index=False)
    summary_rows.append([name, round(CLmax,3), round(a[iCL],2), round(CMmin,3),
                         round(CDmax,3), round(xi,5), round(xi_hat,4), verdict])

    # ---- convergence ----
    pd.DataFrame({"cycle": np.arange(1, len(out["cycle_peakCL"])+1),
                  "peak_CL": out["cycle_peakCL"].round(5),
                  "min_CM": out["cycle_minCM"].round(5),
                  "peakCL_residual": np.abs(np.concatenate(
                      [[np.nan], np.diff(out["cycle_peakCL"])])).round(6)
                  }).to_csv(SOL/"convergence"/f"residuals_{name}.csv", index=False)

    # ---- surface Cp distributions at phases ----
    cp_rows = []
    phase_specs = [(C["a_mean"], True, "mean_up"),
                   (C["a_mean"]+C["a_amp"]*0.7, True, "pre_stall_up"),
                   (min(C["a_mean"]+C["a_amp"], ALPHA_SAMPLE_MAX), True, "near_peak"),
                   (C["a_mean"], False, "mean_down")]
    for tgt, ups, tag in phase_specs:
        j = phase_index(out, tgt, upstroke=ups)
        xoc, cp, upper = us.surface_cp(GEO, C["c"], C["U"], C["M"],
                                       out["alpha_deg"][j], out["CL"][j], out["CNv"][j],
                                       out["tau_v"][j]/consts["Tvl"])
        for k_ in range(len(xoc)):
            cp_rows.append([tag, round(out["alpha_deg"][j],2),
                            "upper" if upper[k_] else "lower",
                            round(xoc[k_],4), round(float(cp[k_]),4)])
    pd.DataFrame(cp_rows, columns=["phase_tag","alpha_deg","surface","x_c","Cp"]
                 ).to_csv(SOL/f"cp_distribution_{name}.csv", index=False)

    # ---- 2D reconstructed fields at key phases ----
    # rise/peak/fall are selected on ANGLE OF ATTACK. On their own they miss the
    # instant the dynamic-stall vortex is actually strongest -- the phenomenon
    # this study is about -- so the vortex maximum is written as a fourth field.
    field_specs = [(C["a_mean"]+C["a_amp"]*0.5, True, "rise"),
                   (min(C["a_mean"]+C["a_amp"], ALPHA_SAMPLE_MAX), True, "peak"),
                   (C["a_mean"]+C["a_amp"]*0.5, False, "fall"),
                   (None, None, "dsv")]
    fphases = []
    for tgt, ups, tag in field_specs:
        j = int(np.argmax(out["CNv"])) if tag == "dsv" else phase_index(out, tgt, upstroke=ups)
        fld = us.reconstruct_field(GEO, C["c"], C["U"], C["M"], out["alpha_deg"][j],
                                   out["CL"][j], out["CNv"][j],
                                   out["tau_v"][j]/consts["Tvl"], T_inf=C["T_inf"],
                                   nx_grid=220, ny_grid=170, **THERMO)
        dff = pd.DataFrame({
            "x_m": fld["X"].ravel().round(5), "y_m": fld["Y"].ravel().round(5),
            "u_ms": fld["u"].ravel().round(3), "v_ms": fld["v"].ravel().round(3),
            "speed_ms": fld["speed"].ravel().round(3), "Cp": fld["Cp"].ravel().round(4),
            "T_static_K": fld["T_static"].ravel().round(3),
            "T_recovery_K": fld["T_recovery"].ravel().round(3),
            "Mach_local": fld["Mlocal"].ravel().round(4),
            "vorticity_1s": fld["vort"].ravel().round(2)})
        adeg = round(out["alpha_deg"][j],1)
        dff.to_csv(SOL/f"field_{name}_{tag}_a{adeg:.0f}.csv", index=False)
        fphases.append((tag, adeg, fld["xv"], fld["yv"]))
    print(f"[run] {name}: CLmax={CLmax:.2f}@{a[iCL]:.1f}deg CMmin={CMmin:.3f} "
          f"onset={onset:.2f}deg Cp-closure={cp_closure_pct:+.1f}% "
          f"Cp-TEjump(max)={cp_te_jump:.3f} cpu={solve_s:.2f}s "
          f"CNvmax={out["CNv"].max():.3f}@a{out["alpha_deg"][int(np.argmax(out["CNv"]))]:.1f} "
          f"CDmax={CDmax:.3f} Xi={xi:.5f} (norm {xi_hat:+.4f} -> {verdict}) "
          f"fields={[f[0] for f in fphases]}")

pd.DataFrame(summary_rows, columns=["case","CL_max","alpha_CLmax_deg","CM_min",
             "CD_max","aero_damping_Xi","aero_damping_Xi_norm","flutter_risk"]
             ).to_csv(SOL/"summary_all_cases.csv", index=False)
# the timings in metrics_*.csv mean nothing without the machine they were taken
# on, so it is recorded alongside them rather than left implicit
_cpu_A = solve_cpu["A_validation"]
pd.DataFrame({"property": ["python", "platform", "processor", "steps_per_cycle",
                           "n_cycles", "solve_cpu_time_s_case_A",
                           "solve_cpu_time_s_case_B", "solve_cpu_ms_per_cycle_case_A",
                           "solve_cpu_us_per_step", "timing_note"],
              "value": [platform.python_version(), platform.platform(),
                        platform.processor() or "unknown", NPC, NCYC,
                        round(_cpu_A, 2), round(solve_cpu["B_application"], 2),
                        round(1000.0*_cpu_A/NCYC, 0),
                        round(1e6*_cpu_A/(NPC*NCYC), 0),
                        "CPU time, not wall time: wall time on this machine varied "
                        "0.57-3.96 s for identical code purely from background load. "
                        "Machine-dependent, so it lives here and not in metrics_*.csv, "
                        "which stays bit-reproducible."]}
             ).to_csv(SOL/"runtime_environment.csv", index=False)
print("[run] done. solution written to 05_solution/")
