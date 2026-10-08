# Tables

Author: Akosa Samuel Onyejekwe (independent)

Written by `python3 -m unistall.make_tables` from the result files; no number is typed by hand.

## Table 1. Measured loops used

NACA 0012, NASA TM-84245 (`results/validation_table.csv`, `data/data_manifest.csv`).

| Set | Inside the static Mach range | Loops | Mach | Reduced frequency | Mean incidence, deg | Amplitude, deg |
|---|---|---|---|---|---|---|
| calibration | no | 4 | 0.04 to 0.18 | 0.099 to 0.248 | 10 to 15 | 10 to 10 |
| calibration | yes | 20 | 0.29 to 0.30 | 0.010 to 0.201 | 4 to 15 | 5 to 10 |
| held out | no | 15 | 0.04 to 0.18 | 0.010 to 0.284 | 8 to 15 | 5 to 10 |
| held out | yes | 45 | 0.22 to 0.30 | 0.010 to 0.202 | 3 to 15 | 5 to 10 |

## Table 2. Static inputs at the two Mach stations

(`results/static_by_mach_naca0012.csv`)

| Mach | Lift-curve slope, per deg | Zero-lift incidence, deg | Maximum static C L | Static-stall incidence, deg | Stall level CN1 | Source sweep |
|---|---|---|---|---|---|---|
| 0.215 | 0.1180 | -0.01 | 1.644 | 16.2 | 1.598 | 13308 |
| 0.302 | 0.1170 | 0.24 | 1.364 | 13.6 | 1.336 | 12102 |

## Table 3. Constants

(`results/calibrated_constants.json`, `results/calibrated_constants_reference.json`)

| Constant | Literature | Tabulated model | Fitted-law model | Status | Range over starts |
|---|---|---|---|---|---|
| pressure lag T_p | 1.700 | 1.652 | 1.602 | fitted; the two refined starts differ by 0.835 | 1.25 to 2.03 |
| boundary-layer lag T_f | 3.000 | 6.278 | 6.314 | fitted; the two refined starts differ by 2.533 | 4.61 to 7.36 |
| vortex decay T_v | 6.000 | 6.318 | 7.269 | fitted; the two refined starts differ by 5.902 | 2.65 to 9.30 |
| vortex travel time T_vl | 7.000 | 8.392 | 8.274 | fitted; the two refined starts differ by 5.928 | 6.84 to 12.99 |
| chord-force recovery eta | 0.950 | 0.907 | 0.909 | fitted; the two refined starts differ by 0.002 | 0.91 to 0.91 |
| unsteady-moment factor g_M | 1.000 | 0.728 | 0.728 | empirical, from loops below static stall |  |

## Table 4. Unsteady moment below static stall

Calibration loops whose peak incidence is below static stall; the empirical factor 0.728 is the value that minimises their moment loop error, and the loops with a closed moment loop would alone choose 0.598 and 0.707. The static part is the damping the model gives with the unsteady moment terms removed, which comes from reading the static table at the lagged incidence (`results/attached_moment_factor.json`).

| Frame | Reduced frequency k | Damping, measured | Damping, full strength | Damping, with factor | Damping, static part | Damping, full strength, static table unlagged | Moment error, full strength | Moment error, with factor | Moment error, full strength, static table unlagged | Factor this loop would choose | Peak vortex normal force | Peak lagged normal force / onset level |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| frame_7110 | 0.100 | 0.160 | 0.212 | 0.165 | 0.042 | 0.170 | 0.139 | 0.118 | 0.088 | 0.707 | 0.016 | 1.054 |
| frame_10218 | 0.010 |  | 0.022 | 0.017 | 0.004 | 0.018 | 0.098 | 0.100 | 0.097 | 1.199 | 0.000 | 0.851 |
| frame_10221 | 0.099 | 0.140 | 0.200 | 0.153 | 0.030 | 0.170 | 0.171 | 0.117 | 0.128 | 0.598 | 0.000 | 0.801 |

## Table 5. Accuracy on the held-out loops

(`results/validation_targets.csv`, `results/validation_summary.csv`)

| Measure | Loops | Measured | 95 % interval | Target | Met |
|---|---|---|---|---|---|
| Loop error, lift | 45 | 0.107 | 0.091 to 0.124 | at most 0.10 | no |
| Loop error, moment | 45 | 0.196 | 0.173 to 0.220 | at most 0.15 | no |
| Loop error, drag | 45 | 0.153 | 0.131 to 0.178 | at most 0.20 | yes |
| Incidence of maximum lift, deg | 45 | 0.572 | 0.413 to 0.774 | at most 0.72 | yes |
| Incidence of moment stall, deg | 29 | 0.622 | 0.392 to 0.884 | at most 1.00 | yes |
| Fraction of loops with the damping sign correct | 45 | 0.600 | 0.466 to 0.733 | at least 0.80 | no |

## Table 6. Accuracy by stall depth, by Mach number and without repeated conditions

(`results/validation_by_stall_depth.csv`, `results/validation_summary.csv`)

| Group | Loops | Lift loop error | Moment loop error | Drag loop error | Damping sign agreement | Mean abs measured damping |
|---|---|---|---|---|---|---|
| attached or marginal (peak less than 2 deg beyond static stall) | 19 | 0.093 | 0.241 | 0.209 | 0.474 | 0.098 |
| light stall (2 to 6 deg beyond) | 7 | 0.131 | 0.155 | 0.152 | 0.429 | 0.056 |
| deep stall (6 deg or more beyond) | 19 | 0.112 | 0.166 | 0.096 | 0.789 | 0.125 |

| Group | Loops | Lift | Moment | Drag | Incidence of maximum lift, deg | Incidence of moment stall, deg | Damping sign correct |
|---|---|---|---|---|---|---|---|
| all primary held-out loops | 45 | 0.107 | 0.196 | 0.153 | 0.572 | 0.622 | 0.600 |
| Mach 0.285 and above (where the calibration loops are) | 41 | 0.107 | 0.198 | 0.157 | 0.451 | 0.429 | 0.561 |
| below Mach 0.285 | 4 | 0.102 | 0.175 | 0.112 | 1.819 | 1.826 | 1.000 |
| not repeating a calibration condition | 39 | 0.106 | 0.198 | 0.158 | 0.562 | 0.632 | 0.590 |
| held out, below the static Mach range (no claim) | 15 | 0.146 | 0.437 | 0.411 | 1.286 | 1.569 | 0.533 |
| calibration loops | 20 | 0.082 | 0.187 | 0.183 | 0.413 | 0.281 | 0.722 |

## Table 6a. Moment stall: measured, predicted, missed

The incidence error of moment stall in Table 5 is defined only for loops where both the measurement and the model show one on the up-stroke (`results/validation_moment_stall.csv`).

| Model | Loops | Moment stall in measurement and model | Measured, missed by the model | Predicted, not measured | In neither | Mean incidence error where both, deg | Mean abs error in minimum moment | Mean measured minimum moment |
|---|---|---|---|---|---|---|---|---|
| tabulated | 45 | 29 | 3 | 0 | 13 | 0.622 | 0.040 | -0.162 |
| reference | 45 | 29 | 3 | 0 | 13 | 0.648 | 0.039 | -0.162 |

## Table 6b. Intervals when repeated conditions are not counted as independent

The held-out means with a 95 % interval from resampling the loops one by one and from resampling groups of loops that repeat one condition (`results/validation_intervals_by_condition.csv`).

| Measure | N loops | N conditions | Value | Ci95 lo by loop | Ci95 hi by loop | Ci95 lo by condition | Ci95 hi by condition |
|---|---|---|---|---|---|---|---|
| mean_nRMS_CL | 45 | 38 | 0.107 | 0.091 | 0.124 | 0.090 | 0.125 |
| mean_nRMS_CM | 45 | 38 | 0.196 | 0.172 | 0.220 | 0.168 | 0.223 |
| mean_nRMS_CD | 45 | 38 | 0.153 | 0.130 | 0.177 | 0.128 | 0.178 |
| mean_abs_dalpha_CLmax_deg | 45 | 38 | 0.572 | 0.413 | 0.774 | 0.414 | 0.788 |
| mean_abs_dalpha_Mstall_deg | 29 | 38 | 0.622 | 0.400 | 0.873 | 0.384 | 0.883 |
| mean_Xi_sign_agree | 45 | 38 | 0.600 | 0.444 | 0.733 | 0.442 | 0.761 |

## Table 7. Tabulated against fitted separation law

Mean absolute error on the held-out loops; difference is tabulated minus fitted (`results/comparison.csv`, `results/comparison_same_constants.csv`).

| Constants | Measure | Tabulated | Fitted law | Difference | 95 % interval | Verdict |
|---|---|---|---|---|---|---|
| each model's own calibration | \|nRMS_CL\| | 0.1067 | 0.1056 | +0.0011 | -0.0007 to +0.0027 | no difference at 95 % |
| each model's own calibration | \|nRMS_CM\| | 0.1959 | 0.1963 | -0.0004 | -0.0017 to +0.0009 | no difference at 95 % |
| each model's own calibration | \|nRMS_CD\| | 0.1526 | 0.1512 | +0.0014 | +0.0002 to +0.0024 | reference better |
| each model's own calibration | \|CLmax_err_pct\| | 4.5311 | 4.3475 | +0.1836 | +0.0910 to +0.2789 | reference better |
| each model's own calibration | \|dalpha_CLmax_deg\| | 0.5723 | 0.5392 | +0.0332 | -0.0154 to +0.0781 | no difference at 95 % |
| each model's own calibration | \|dalpha_Mstall_deg\| | 0.6216 | 0.6476 | -0.0259 | -0.0327 to -0.0191 | tabulated better |
| literature values | \|nRMS_CL\| | 0.1332 | 0.1315 | +0.0017 | +0.0000 to +0.0033 | reference better |
| literature values | \|nRMS_CM\| | 0.2161 | 0.2145 | +0.0017 | +0.0002 to +0.0034 | reference better |
| literature values | \|nRMS_CD\| | 0.1641 | 0.1633 | +0.0008 | -0.0004 to +0.0021 | no difference at 95 % |
| literature values | \|CLmax_err_pct\| | 4.7574 | 4.8025 | -0.0451 | -0.1834 to +0.1085 | no difference at 95 % |
| literature values | \|dalpha_CLmax_deg\| | 0.8938 | 0.9116 | -0.0178 | -0.0697 to +0.0183 | no difference at 95 % |
| literature values | \|dalpha_Mstall_deg\| | 0.7817 | 0.7871 | -0.0054 | -0.0107 to -0.0009 | tabulated better |
| fitted for the tabulated model | \|nRMS_CL\| | 0.1067 | 0.1054 | +0.0014 | -0.0006 to +0.0030 | no difference at 95 % |
| fitted for the tabulated model | \|nRMS_CM\| | 0.1959 | 0.1945 | +0.0014 | +0.0002 to +0.0027 | reference better |
| fitted for the tabulated model | \|nRMS_CD\| | 0.1526 | 0.1520 | +0.0005 | -0.0007 to +0.0016 | no difference at 95 % |
| fitted for the tabulated model | \|CLmax_err_pct\| | 4.5311 | 4.4242 | +0.1069 | +0.0051 to +0.2133 | reference better |
| fitted for the tabulated model | \|dalpha_CLmax_deg\| | 0.5723 | 0.5628 | +0.0095 | -0.0298 to +0.0433 | no difference at 95 % |
| fitted for the tabulated model | \|dalpha_Mstall_deg\| | 0.6216 | 0.6228 | -0.0011 | -0.0028 to +0.0004 | no difference at 95 % |
| fitted for the reference model | \|nRMS_CL\| | 0.1058 | 0.1056 | +0.0002 | -0.0032 to +0.0027 | no difference at 95 % |
| fitted for the reference model | \|nRMS_CM\| | 0.1969 | 0.1963 | +0.0006 | -0.0015 to +0.0024 | no difference at 95 % |
| fitted for the reference model | \|nRMS_CD\| | 0.1510 | 0.1512 | -0.0002 | -0.0028 to +0.0016 | no difference at 95 % |
| fitted for the reference model | \|CLmax_err_pct\| | 4.4817 | 4.3475 | +0.1342 | +0.0278 to +0.2430 | reference better |
| fitted for the reference model | \|dalpha_CLmax_deg\| | 0.5295 | 0.5392 | -0.0097 | -0.0427 to +0.0103 | no difference at 95 % |
| fitted for the reference model | \|dalpha_Mstall_deg\| | 0.6442 | 0.6476 | -0.0034 | -0.0078 to -0.0003 | tabulated better |

## Table 7a. Loads predicted by the two separation laws with identical constants

RMS difference between the two predictions over the cycle, divided by the measured range of the coefficient and averaged over the held-out loops, against the margin recorded beforehand in `data/targets.json` (`results/equivalence.csv`).

| Constants | Coefficient | Margin | Mean RMS difference / range | 95 % interval | Worst loop | Largest difference at any instant (coefficient) | Loads within margin | Loop errors within margin |
|---|---|---|---|---|---|---|---|---|
| literature values | CL | 0.010 | 0.0116 | 0.0094 to 0.0140 | 0.0401 | 0.0775 | no | yes |
| literature values | CM | 0.015 | 0.0061 | 0.0048 to 0.0076 | 0.0225 | 0.0210 | yes | yes |
| literature values | CD | 0.020 | 0.0096 | 0.0079 to 0.0115 | 0.0266 | 0.0372 | yes | yes |
| fitted for the tabulated model | CL | 0.010 | 0.0101 | 0.0079 to 0.0125 | 0.0391 | 0.0731 | no | yes |
| fitted for the tabulated model | CM | 0.015 | 0.0052 | 0.0042 to 0.0063 | 0.0143 | 0.0244 | yes | yes |
| fitted for the tabulated model | CD | 0.020 | 0.0086 | 0.0069 to 0.0104 | 0.0260 | 0.0322 | yes | yes |
| fitted for the reference model | CL | 0.010 | 0.0108 | 0.0082 to 0.0138 | 0.0470 | 0.2430 | no | yes |
| fitted for the reference model | CM | 0.015 | 0.0074 | 0.0049 to 0.0114 | 0.0805 | 0.1049 | yes | yes |
| fitted for the reference model | CD | 0.020 | 0.0093 | 0.0073 to 0.0116 | 0.0382 | 0.0832 | yes | yes |

## Table 8. Cycle damping

(`results/validation_damping.csv`)

| Group | Loops | Measured negative | Predicted negative | Sign agreement | Sign agreement of a prediction of positive damping everywhere | Correlation | Rms error | Mean error | Mean measured | Mean predicted | Mean uncertainty random | Mean uncertainty stroke |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| all primary held-out loops with a closed moment loop | 45 | 15 | 3 | 0.600 | 0.667 | 0.330 | 0.120 | 0.029 | 0.064 | 0.093 | 0.015 | 0.085 |
| measured damping positive | 30 | 0 | 3 | 0.900 | 1.000 | 0.522 | 0.089 | -0.030 | 0.125 | 0.095 | 0.013 | 0.078 |
| measured damping negative | 15 | 15 | 0 | 0.000 | 0.000 | 0.241 | 0.165 | 0.147 | -0.058 | 0.089 | 0.018 | 0.098 |
| deep stall (6 deg or more beyond static stall) | 19 | 3 | 1 | 0.789 | 0.842 | 0.554 | 0.090 | -0.048 | 0.115 | 0.067 | 0.010 | 0.066 |
| measured sign certain if the moment error is random (beyond 2 sigma) | 35 | 9 | 3 | 0.657 | 0.743 | 0.435 | 0.114 | 0.007 | 0.083 | 0.089 | 0.014 | 0.081 |
| measured sign certain even if the moment error follows the stroke | 21 | 3 | 0 | 0.857 | 0.857 | 0.684 | 0.099 | -0.026 | 0.132 | 0.106 | 0.012 | 0.069 |

## Table 9. Attached flow against Theodorsen's solution

Largest error over the reduced frequencies tested (`results/attached_checks.csv`).

| Variant | Lift amplitude, % | Lift phase, deg | Moment amplitude, % | Moment phase, deg |
|---|---|---|---|---|
| the model (low-Mach blend) | 1.86 | 0.62 | 0.01 | 0.06 |
| incompressible form alone | 1.69 | 0.63 | 0.00 | 0.00 |
| compressible indicial constants alone, M = 0.05 | 6.43 | 6.54 | 12.90 | 11.18 |
| compressible constants with the simplified impulsive moment [DH19 Eq. 1.27] | 6.43 | 6.54 | 49.76 | 9.02 |

## Table 9a. Compressible attached flow against exact linear theory

Starting values after a step and the exact short-time lift, at the Mach numbers of the study, where Theodorsen's solution does not apply (`results/indicial_checks.csv`).

| Check | Mach number | Quantity | Model | Theory | Error pct |
|---|---|---|---|---|---|
| starting value | 0.220 | C_N per unit step in incidence | 18.1772 | 18.1818 | -0.03 |
| starting value | 0.220 | C_M about c/4 per unit step in incidence | -4.5446 | -4.5454 | -0.02 |
| starting value | 0.220 | C_N per unit step in pitch rate | 4.5442 | 4.5454 | -0.03 |
| starting value | 0.220 | C_M about c/4 per unit step in pitch rate | -2.6498 | -2.6515 | -0.06 |
| short-time lift after a step in incidence | 0.220 | largest difference over 0 <= s <= 0.361 semichords (at s = 0.361) | 8.0916 | 6.5592 | +23.36 |
| starting value | 0.300 | C_N per unit step in incidence | 13.3309 | 13.3333 | -0.02 |
| starting value | 0.300 | C_M about c/4 per unit step in incidence | -3.3329 | -3.3333 | -0.01 |
| starting value | 0.300 | C_N per unit step in pitch rate | 3.3327 | 3.3333 | -0.02 |
| starting value | 0.300 | C_M about c/4 per unit step in pitch rate | -1.9437 | -1.9444 | -0.04 |
| short-time lift after a step in incidence | 0.300 | largest difference over 0 <= s <= 0.462 semichords (at s = 0.462) | 6.8346 | 6.1544 | +11.05 |

## Table 9b. The model beside the published Leishman-Beddoes model

Two NACA 0012 cases at Mach 0.3 and reduced frequency 0.1 for which Leishman and Crouse (1989) plot their model against the measurements. Loop errors of the published curves and of this model against the same traced measured points, and the RMS difference between the two models over the cycle as a fraction of the measured range (`results/published_reference.csv`; the curves are tracings of printed figures).

| Case | Coefficient | Measured range | Loop error published model | Loop error this model | Rms difference between models over range | Extreme published model | Extreme this model | Extreme measured |
|---|---|---|---|---|---|---|---|---|
| 9.9 +/- 9.9 deg, k 0.1, M 0.3 | CN | 1.782 | 0.0381 | 0.0775 | 0.0598 | 1.810 | 1.820 | 1.830 |
| 9.9 +/- 9.9 deg, k 0.1, M 0.3 | CM | 0.294 | 0.0746 | 0.1128 | 0.0969 | -0.263 | -0.286 | -0.229 |
| 15 +/- 10 deg, k 0.1, M 0.3 | CN | 1.756 | 0.0718 | 0.0810 | 0.0661 | 1.948 | 1.933 | 2.089 |
| 15 +/- 10 deg, k 0.1, M 0.3 | CM | 0.395 | 0.1097 | 0.1266 | 0.0815 | -0.348 | -0.320 | -0.314 |

## Table 10. Second aerofoil (Ames A-01)

Dynamic constants unchanged (`results/validation_other_airfoil.csv`).

| Static inputs | Loops | Loop error, lift | Loop error, moment | Loop error, drag | Peak-lift error, % |
|---|---|---|---|---|---|
| NACA 0012 static data | 38 | 0.223 | 0.372 | 0.261 | 14.070 |
| own static sweep | 38 | 0.216 | 0.390 | 0.308 | 9.487 |

## Table 11. Every target against its measured value

(`data/targets.json`, `results/targets_scoreboard.csv`)

| Group | Kind of target | Target | Kind | Value | Measured | Error implied by the quoted uncertainty | Met | Source |
|---|---|---|---|---|---|---|---|---|
| attached_flow | physical | lift_and_moment_amplitude_error_max_pct | at most | 3 | 1.86 |  | yes | attached_checks.csv |
| attached_flow | physical | lift_and_moment_phase_error_max_deg | at most | 3 | 0.62 |  | yes | attached_checks.csv |
| attached_flow | physical | chord_speed_independence_max_abs_dCL | at most | 1e-10 | 5.04e-14 |  | yes | attached_checks.csv |
| attached_flow | physical | no_stall_frame_nRMS_CL_max | at most | 0.05 | 0.0425 |  | yes | attached_checks.csv |
| attached_flow | physical | no_stall_frame_nRMS_CM_max | at most | 0.1 | 0.1235 |  | no | attached_checks.csv |
| separated_flow | physical | CD_negative_fraction_of_cycle_max | at most | 0.05 | 0 |  | yes | separated_checks.csv |
| separated_flow | physical | CD_min | at least | -0.01 | 0.0021 |  | yes | separated_checks.csv |
| separated_flow | physical | onset_shift_with_4x_steps_max_deg | at most | 0.02 | 0.0011 |  | yes | separated_checks.csv |
| calibration | numerical | repeatability_max_abs_difference | at most | 1e-06 | 1.11e-16 |  | yes | calibration_repeatability.json |
| held_out_validation | physical | mean_nRMS_CL_max | at most | 0.1 | 0.1067 | 0.028 | no | validation_targets.csv |
| held_out_validation | physical | mean_nRMS_CM_max | at most | 0.15 | 0.1959 | 0.059 | no | validation_targets.csv |
| held_out_validation | physical | mean_nRMS_CD_max | at most | 0.2 | 0.1526 | 0.059 | yes | validation_targets.csv |
| held_out_validation | physical | mean_abs_dalpha_CLmax_deg_max | at most | 0.7238 | 0.5723 |  | yes | validation_targets.csv |
| held_out_validation | physical | mean_abs_dalpha_Mstall_deg_max | at most | 1 | 0.6216 |  | yes | validation_targets.csv: judged on the 29 loops where measurement and model both show a moment stall; the model misses a measured one on 3 loops and predicts an unmeasured one on 0 |
| held_out_validation | physical | damping_sign_agreement_min | at least | 0.8 | 0.6 |  | no | validation_targets.csv |
| structural | numerical | cycle_work_identity_max_error_pct | at most | 1 | 0.3814 |  | yes | structural_checks.csv |
| structural | numerical | free_response_sign_agreement_min_conditions | at least | 3 | 6 |  | yes | structural_checks.csv |
| structural | numerical | dimensional_load_identity_max_abs | at most | 1e-09 | 1.138e-15 |  | yes | structural_checks.csv |
| structural | numerical | second_axis_theodorsen_amplitude_max_pct | at most | 3 | 1.69 |  | yes | attached_checks.csv |
| structural | numerical | second_axis_theodorsen_phase_max_deg | at most | 3 | 0.63 |  | yes | attached_checks.csv |
| code_quality | software | regression_max_abs_difference | at most | 1e-12 | 0 |  | yes | code_quality.csv |
| code_quality | software | import_path_edits_load_model_path_max | at most | 0 | 0 |  | yes | code_quality.csv |
| code_quality | software | module_level_executable_lines_per_load_model_file_max | at most | 30 | 1 |  | yes | code_quality.csv |
| code_quality | software | longest_package_function_lines_max | at most | 80 | 56 |  | yes | code_quality.csv |
| code_quality | software | unused_parameters_max | at most | 0 | 0 |  | yes | code_quality.csv |
| code_quality | software | public_functions_typed_and_unit_documented_pct_min | at least | 100 | 100 |  | yes | code_quality.csv |
| code_quality | software | lint_findings_max | at most | 0 | 0 |  | yes | code_quality.csv |
| code_quality | software | silently_skipped_checks_max | at most | 0 | 0 |  | yes | code_quality.csv |
| code_quality | software | revision_history_phrases_in_package_max | at most | 0 | 0 |  | yes | code_quality.csv |
| numerical | numerical | loop_error_change_max | at most | 0.002 | 0.000602 |  | yes | convergence_summary.csv |
| numerical | numerical | state_space_loop_error_difference_max | at most | 0.002 | 0.000398 |  | yes | statespace_summary.csv |
| equivalence_of_separation_laws | physical | margin_CL | at most | 0.01 | 0.01402 |  | no | equivalence.csv |
| equivalence_of_separation_laws | physical | margin_CM | at most | 0.015 | 0.01143 |  | yes | equivalence.csv |
| equivalence_of_separation_laws | physical | margin_CD | at most | 0.02 | 0.01159 |  | yes | equivalence.csv |

## Table 12. Step size and cycles marched

The default resolution against a march with four times as many steps and one with four more cycles, on the calibration loops (`results/convergence_summary.csv`).

| Study | Default resolution compared with | Loops | Largest change in a loop error | Largest change in lift | Largest change in moment | Largest change in drag | Largest error against zero step (extrapolated) | Target | Met |
|---|---|---|---|---|---|---|---|---|---|
| step size | 4 times as many steps | 20 | 0.0008 | 0.0062 | 0.0023 | 0.0017 | 0.001 | 0.002 | yes |
| cycles marched | 4 more cycles | 20 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |  | 0.002 | yes |

## Table 12a. State-space form against the indicial march

The same model integrated as differential equations, on the calibration loops (`results/statespace_summary.csv`).

| Quantity | Value |
|---|---|
| n_frames | 24 |
| n_frames_below_mach_0_20 | 4 |
| lowest_mach | 0.035 |
| states | 14 |
| integrator | classical fourth-order Runge-Kutta, fixed step |
| worst_difference_in_loop_error | 0.000398 |
| mean_difference_in_loop_error | 7.2e-05 |
| worst_difference_in_CL | 0.005225 |
| worst_difference_in_CM | 0.004843 |
| worst_difference_in_CD | 0.001692 |
| worst_difference_in_CL_pct_of_range | 0.5275 |
| worst_difference_in_CM_pct_of_range | 1.9679 |
| worst_difference_in_CD_pct_of_range | 3.6068 |
| worst_difference_in_damping | 0.00027 |
| cpu_ratio_state_space_to_indicial | 4.94 |
| target_difference_in_loop_error | 0.002 |
| met | True |

## Table 12b. The model at its limits

Written by `python3 -m unistall.check_limits`. A loop and its mirror image about the zero-lift incidence: sums and differences that are zero for a model symmetric about it (`results/limit_mirror_pair.csv`).

| Mean incidence, deg | Amplitude, deg | Reduced frequency k | Mach number | Mean incidence, deg | Largest sum of normal force, loop and mirror image | Largest sum of moment, loop and mirror image | Largest difference of chord force, loop and mirror image | Largest sum of lift, loop and mirror image | Largest difference of drag, loop and mirror image | Range of drag over the cycle |
|---|---|---|---|---|---|---|---|---|---|---|
| 10.000 | 10.000 | 0.100 | 0.300 | 0.234 | 2.0e-14 | 0.0103 | 0.0112 | 0.0040 | 0.0126 | 0.552 |
| 15.000 | 10.000 | 0.100 | 0.290 | 0.205 | 2.1e-14 | 0.0097 | 0.0095 | 0.0047 | 0.0125 | 0.690 |
| 9.000 | 5.000 | 0.100 | 0.300 | 0.234 | 1.5e-14 | 0.0103 | 0.0104 | 0.0020 | 0.0023 | 0.068 |
| 12.000 | 8.000 | 0.050 | 0.220 | 0.004 | 2.8e-14 | 0.0058 | 0.0002 | 0.0000 | 0.0002 | 0.332 |

Drag of a slow sweep against the static drag data at Mach 0.30 (`results/limit_static_drag.csv`).

| Incidence, deg | Static drag, measured | Kind | Drag, model slow sweep | Drag, model slow sweep, without CD0 | Error | Error without CD0 |
|---|---|---|---|---|---|---|
| -5.000 | 0.008 | total drag, wake survey | 0.012 | 0.005 | 0.004 | -0.004 |
| -2.000 | 0.007 | total drag, wake survey | 0.008 | 0.001 | 0.001 | -0.006 |
| 0.000 | 0.007 | total drag, wake survey | 0.007 | -0.000 | 0.000 | -0.007 |
| 2.000 | 0.007 | total drag, wake survey | 0.008 | 0.001 | 0.001 | -0.006 |
| 5.000 | 0.009 | total drag, wake survey | 0.012 | 0.005 | 0.003 | -0.004 |
| 8.000 | 0.010 | total drag, wake survey | 0.019 | 0.012 | 0.009 | 0.002 |
| 10.000 | 0.012 | total drag, wake survey | 0.026 | 0.018 | 0.014 | 0.006 |
| 12.000 | 0.017 | total drag, wake survey | 0.035 | 0.028 | 0.018 | 0.011 |
| 14.000 | 0.029 | total drag, wake survey | 0.093 | 0.086 | 0.064 | 0.057 |
| 9.980 | 0.018 | pressure (form) drag, excludes skin friction; not wall-corrected | 0.025 | 0.018 | 0.007 | 0.000 |
| 11.100 | 0.023 | pressure (form) drag, excludes skin friction; not wall-corrected | 0.029 | 0.022 | 0.006 | -0.001 |
| 12.080 | 0.028 | pressure (form) drag, excludes skin friction; not wall-corrected | 0.036 | 0.029 | 0.008 | 0.001 |
| 13.060 | 0.034 | pressure (form) drag, excludes skin friction; not wall-corrected | 0.058 | 0.051 | 0.024 | 0.017 |
| 13.530 | 0.044 | pressure (form) drag, excludes skin friction; not wall-corrected | 0.070 | 0.063 | 0.026 | 0.019 |
| 14.030 | 0.173 | pressure (form) drag, excludes skin friction; not wall-corrected | 0.095 | 0.088 | -0.078 | -0.085 |
| 14.540 | 0.194 | pressure (form) drag, excludes skin friction; not wall-corrected | 0.155 | 0.148 | -0.039 | -0.046 |
| 15.010 | 0.199 | pressure (form) drag, excludes skin friction; not wall-corrected | 0.186 | 0.179 | -0.013 | -0.020 |
| 15.480 | 0.104 | pressure (form) drag, excludes skin friction; not wall-corrected | 0.217 | 0.209 | 0.113 | 0.105 |
| 15.990 | 0.116 | pressure (form) drag, excludes skin friction; not wall-corrected | 0.229 | 0.222 | 0.113 | 0.106 |
| 17.870 | 0.199 | pressure (form) drag, excludes skin friction; not wall-corrected | 0.269 | 0.262 | 0.070 | 0.063 |
| 20.000 | 0.378 | pressure (form) drag, excludes skin friction; not wall-corrected | 0.316 | 0.309 | -0.062 | -0.069 |

Loops that cross zero lift (`results/limit_zero_lift.csv`).

| Frame | Set | Mach number | Reduced frequency k | Mean incidence, deg | Amplitude, deg | Lowest chord force | Lowest drag | Fraction of cycle with chord force negative | Fraction of cycle with drag negative |
|---|---|---|---|---|---|---|---|---|---|
| frame_9221 | held_out | 0.302 | 0.010 | 9.900 | 9.900 | 0.000 | 0.007 | 0.000 | 0.000 |
| frame_9222 | calibration | 0.302 | 0.024 | 9.900 | 9.900 | 0.000 | 0.007 | 0.000 | 0.000 |
| frame_9223 | held_out | 0.302 | 0.048 | 10.000 | 9.900 | 0.000 | 0.007 | 0.000 | 0.000 |
| frame_9302 | calibration | 0.302 | 0.096 | 9.800 | 9.900 | 0.000 | 0.004 | 0.000 | 0.000 |
| frame_9307 | held_out | 0.302 | 0.145 | 9.900 | 9.900 | 0.000 | -0.001 | 0.000 | 0.049 |
| frame_10218 | calibration | 0.300 | 0.010 | 5.000 | 5.000 | 0.000 | 0.007 | 0.000 | 0.000 |
| frame_10221 | calibration | 0.301 | 0.099 | 5.000 | 5.000 | 0.000 | 0.006 | 0.000 | 0.000 |
| frame_10222 | held_out | 0.301 | 0.198 | 5.000 | 5.000 | 0.000 | 0.002 | 0.000 | 0.000 |
| frame_10303 | held_out | 0.301 | 0.099 | 5.000 | 10.000 | 0.000 | 0.002 | 0.000 | 0.000 |
| frame_10305 | calibration | 0.301 | 0.099 | 3.700 | 10.000 | 0.000 | 0.002 | 0.000 | 0.000 |
| frame_10309 | held_out | 0.301 | 0.099 | 3.000 | 10.000 | 0.000 | 0.002 | 0.000 | 0.000 |

Switching rules on every usable loop (`results/limit_vortex_rules.csv`).

| Set | Loops | Loops with a vortex shed | Loops with repeated shedding | Most sheddings in a cycle | Loops on which the sign rule removes vortex lift | Largest vortex normal force removed |
|---|---|---|---|---|---|---|
| calibration | 20 | 18 | 15 | 31 | 17 | 0.0016 |
| held-out | 45 | 43 | 42 | 20 | 43 | 0.0026 |
| all | 65 | 61 | 57 | 31 | 60 | 0.0026 |

## Table 13. Post-stall static curve

Calibration cost with each of the two measured post-stall states, at the literature stall constants (`results/static_curve_choice.csv`).

| Post-stall static curve | Cost | Mean lift loop error | Mean moment loop error | Mean drag loop error |
|---|---|---|---|---|
| more separated | 2.143 | 0.095 | 0.197 | 0.194 |
| less separated | 3.330 | 0.147 | 0.266 | 0.213 |
