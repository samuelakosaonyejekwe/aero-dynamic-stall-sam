"""
04_solver / run_case.py
-----------------------
Drives the UNISTALL(TM) solver for both case-study configurations and writes
all solution data to 05_solution/.  Outputs:
  time_history_<case>.csv      per-step unsteady loads & states (last cycle)
  cp_distribution_<case>.csv   surface Cp(x/c) at several phase angles
  field_<case>_ph<deg>.csv     reconstructed 2D fields at key phases
  model_static_polar.csv       quasi-steady model polar (for validation)
  metrics_<case>.csv           engineering scalar metrics
  convergence/residuals_<case>.csv   cycle-to-cycle convergence
"""
import sys, json
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
consts = dict(**cfg["indicial_circulatory"], **cfg["time_constants_semichords"])
# ---- calibrated constants (single source of truth; calibrated to REAL NACA0012
#      frame 9302 from NASA TM-84245 — see 06_postprocessing/validation) ----
consts.update({k: v for k, v in cfg["calibrated_constants"].items() if k != "comment"})

# ---- calibrate static separation to published polar ----
stat = pd.read_csv(SETUP/"static_polar_reference.csv")
f_static = us.calibrate_separation(stat["alpha_deg"], stat["Cl"], stat["Cd"], CNALPHA)

# ---- model quasi-steady polar (validation) ----
aq = np.linspace(0, 22, 89)
fq = f_static(aq)
CN_qs = CNALPHA*((1+np.sqrt(fq))/2)**2*np.radians(aq)
CC_qs = 0.95*CNALPHA*np.radians(aq)**2*np.sqrt(fq)
CL_qs = CN_qs*np.cos(np.radians(aq)) + CC_qs*np.sin(np.radians(aq))
CD_qs = CN_qs*np.sin(np.radians(aq)) - CC_qs*np.cos(np.radians(aq)) + 0.0086
pd.DataFrame({"alpha_deg": aq.round(3), "Cl_model": CL_qs.round(4),
              "Cd_model": CD_qs.round(4), "f_sep": fq.round(4)}
             ).to_csv(SOL/"model_static_polar.csv", index=False)

# ---- case definitions ----
flow = pd.read_csv(SETUP/"flow_conditions.csv").set_index("parameter")
kin  = pd.read_csv(SETUP/"kinematics.csv").set_index("case_id")

CASES = {
 "A_validation": dict(col="case_A_validation",
                      a_mean=10.0, a_amp=10.0, k=0.10, M=0.30, c=0.30, U=102.0,
                      T_inf=288.15, phases=[8, 14, 18, 24]),   # deg of cycle... (phase angle wt)
 "B_application": dict(col="case_B_application",
                      a_mean=12.0, a_amp=8.0, k=0.074, M=0.28, c=0.527, U=95.8,
                      T_inf=287.5, phases=[10, 16, 19, 5]),
}

def phase_index(out, target_alpha_deg, upstroke=True):
    """index of nearest matching alpha on up/down stroke."""
    a = out["alpha_deg"]; ad = out["alpha_dot"]
    mask = (ad > 0) if upstroke else (ad < 0)
    idx = np.where(mask)[0]
    j = idx[np.argmin(np.abs(a[idx]-target_alpha_deg))]
    return j

summary_rows = []
for name, C in CASES.items():
    out = us.solve_dynamic_stall(C["a_mean"], C["a_amp"], C["k"], C["M"], C["c"], C["U"],
                                 f_static, CNalpha=CNALPHA, consts=consts,
                                 n_per_cycle=720, n_cycles=6)
    # ---- time history CSV ----
    th = pd.DataFrame({
        "time_s": out["t"], "phase_deg": out["phase_deg"], "alpha_deg": out["alpha_deg"],
        "alpha_dot_rad_s": out["alpha_dot"], "alpha_eff_deg": out["alpha_eff_deg"],
        "CL": out["CL"], "CD": out["CD"], "CM_c4": out["CM"], "CN": out["CN"],
        "CC": out["CC"], "CN_prime": out["CNp"], "CN_vortex": out["CNv"],
        "CN_attached": out["CNf"], "f_separation": out["f_sep"],
        "vortex_active": out["vortex_active"], "tau_v_semichords": out["tau_v"],
    }).round(6)
    th.to_csv(SOL/f"time_history_{name}.csv", index=False)

    # ---- engineering metrics ----
    a = out["alpha_deg"]; up = out["alpha_dot"] > 0
    CLmax = out["CL"].max(); iCL = out["CL"].argmax()
    CMmin = out["CM"].min(); iCM = out["CM"].argmin()
    CDmax = out["CD"].max()
    xi = us.aerodynamic_damping(a, out["CM"])           # net cyclic damping
    # hysteresis loop areas
    loopCL = np.abs(np.trapz(out["CL"], np.radians(a)))
    onset = a[up][np.argmax(out["CN"][up] >= consts["CN1"])] if np.any(out["CN"][up] >= consts["CN1"]) else np.nan
    CL_static_max = CL_qs.max()
    met = pd.DataFrame({
        "metric": ["CL_max_dynamic", "alpha_at_CLmax_deg", "CL_max_static",
                   "dynamic_overshoot_ratio", "CM_min(c/4)", "alpha_at_CMmin_deg",
                   "CD_max", "stall_onset_alpha_deg", "aero_damping_Xi",
                   "stall_flutter_risk", "CL_hysteresis_loop_area",
                   "reduced_frequency_k", "mach_M", "mean_alpha_deg", "amp_alpha_deg"],
        "value": [round(CLmax,3), round(a[iCL],2), round(CL_static_max,3),
                  round(CLmax/CL_static_max,3), round(CMmin,3), round(a[iCM],2),
                  round(CDmax,3), round(float(onset),2), round(xi,4),
                  "HIGH (neg. damping)" if xi < 0 else "low (pos. damping)",
                  round(loopCL,4), C["k"], C["M"], C["a_mean"], C["a_amp"]],
    })
    met.to_csv(SOL/f"metrics_{name}.csv", index=False)
    summary_rows.append([name, round(CLmax,3), round(a[iCL],2), round(CMmin,3),
                         round(CDmax,3), round(xi,4),
                         "HIGH" if xi < 0 else "low"])

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
                   (min(C["a_mean"]+C["a_amp"], 19.5), True, "near_peak"),
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
    field_specs = [(C["a_mean"]+C["a_amp"]*0.5, True, "rise"),
                   (min(C["a_mean"]+C["a_amp"], 19.0), True, "peak"),
                   (C["a_mean"]+C["a_amp"]*0.5, False, "fall")]
    fphases = []
    for tgt, ups, tag in field_specs:
        j = phase_index(out, tgt, upstroke=ups)
        fld = us.reconstruct_field(GEO, C["c"], C["U"], C["M"], out["alpha_deg"][j],
                                   out["CL"][j], out["CNv"][j],
                                   out["tau_v"][j]/consts["Tvl"], T_inf=C["T_inf"],
                                   nx_grid=220, ny_grid=170)
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
          f"CDmax={CDmax:.3f} Xi={xi:.4f} fields={[f[0] for f in fphases]}")

pd.DataFrame(summary_rows, columns=["case","CL_max","alpha_CLmax_deg","CM_min",
             "CD_max","aero_damping_Xi","flutter_risk"]
             ).to_csv(SOL/"summary_all_cases.csv", index=False)
print("[run] done. solution written to 05_solution/")
