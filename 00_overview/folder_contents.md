# Contents of the numbered folders

Written by `00_overview/write_folder_contents.py`. Each row is a kind of file: the pattern of its name, what it is, and what writes it. n is the number of files of that kind present when this list was written.

## 00_overview

| Files | n | What they are | Written by |
|---|---|---|---|
| `case_definition.md` | 1 | the definition of the study in words, every number from 03_model_setup | `03_model_setup/generate_setup.py` |
| `folder_contents.md` | 1 | this list | `00_overview/write_folder_contents.py` |
| `write_folder_contents.py` | 1 | the script that writes this list | written by hand |

## 01_geometry

| Files | n | What they are | Written by |
|---|---|---|---|
| `generate_geometry.py` | 1 | script of the section geometry | written by hand |
| `naca*_coordinates.csv` | 1 | coordinates of the closed section outline, fractions of chord | `01_geometry/generate_geometry.py` |
| `section_geometry_summary.csv` | 1 | thickness, nose radius, form of the trailing edge, area | `01_geometry/generate_geometry.py` |
| `fig_geometry_*.png` | 2 | the profile and its thickness distribution | `01_geometry/generate_geometry.py` |
| `figure_record.csv` | 1 | resolution and smallest lettering of the folder's figures, from which the report works out the size of the lettering as printed | the figure script(s) of the folder |

## 02_mesh

| Files | n | What they are | Written by |
|---|---|---|---|
| `figure_record.csv` | 1 | resolution and smallest lettering of the folder's figures, from which the report works out the size of the lettering as printed | the figure script(s) of the folder |
| `README.md` | 1 | standing of the grid: an illustration, no flow equation solved on it | written by hand |
| `generate_mesh.py` | 1 | script of the O-grid | written by hand |
| `field_on_mesh.py` | 1 | script that evaluates the reconstructed field at the grid nodes | written by hand |
| `mesh_*.csv` | 3 | grid nodes, quality measures and the wall-normal spacing law | `02_mesh/generate_mesh.py` |
| `fig_mesh_field_*.png` | 2 | corrected pressure coefficient of the reconstruction at the grid nodes, at peak lift | `02_mesh/field_on_mesh.py` |
| `fig_mesh_*.png` | 4 | pictures of the grid | `02_mesh/generate_mesh.py` |

## 03_model_setup

| Files | n | What they are | Written by |
|---|---|---|---|
| `generate_setup.py` | 1 | script of the conditions of the two cases | written by hand |
| `flow_conditions.csv` | 1 | chord, Mach number, stream speed, air state, Reynolds number and reduced frequency of each case, each row with its source | `03_model_setup/generate_setup.py` |
| `kinematics.csv` | 1 | the prescribed pitching motion of each case | `03_model_setup/generate_setup.py` |
| `station_condition.csv` | 1 | how the Mach number and reduced frequency of Case B follow from the rotor | `03_model_setup/generate_setup.py` |
| `air_properties.csv` | 1 | the gas properties the cases are run with | `03_model_setup/generate_setup.py` |
| `static_inputs.csv` | 1 | the static normal force, moment and separation point the model reads at each case's Mach number; they come from measurement | `03_model_setup/generate_setup.py` |
| `solver_config.json` | 1 | solver name, calibrated constants, march of each case, settings of the reconstruction | `03_model_setup/generate_setup.py` |

## 04_solver

| Files | n | What they are | Written by |
|---|---|---|---|
| `README.md` | 1 | what the folder holds: drivers only; the model is the package unistall/ | written by hand |
| `run_case.py` | 1 | driver of the two cases: writes 05_solution | written by hand |
| `unistall_solver.py` | 1 | command-line entry to the packaged solver for one condition | written by hand |

## 05_solution

| Files | n | What they are | Written by |
|---|---|---|---|
| `README.md` | 1 | standing of the files: loads are the model's prediction, fields a reconstruction | written by hand |
| `time_history_*.csv` | 2 | loads and model states at every step of the reported cycle | `04_solver/run_case.py` |
| `metrics_*.csv` | 2 | scalar results of a case | `04_solver/run_case.py` |
| `summary_all_cases.csv` | 1 | one row of main results per case | `04_solver/run_case.py` |
| `held_out_error_bands.csv` | 1 | error of the model in each headline quantity on the held-out loops | `04_solver/run_case.py` |
| `static_station_spread.csv` | 1 | Case B with the static inputs of each Mach station | `04_solver/run_case.py` |
| `model_static_polar_*.csv` | 2 | lift, drag, moment and separation point of the model in a slow sweep: its quasi-steady limit, not a measurement | `04_solver/run_case.py` |
| `response_surface.csv` | 1 | peak lift, minimum moment and damping over mean incidence and reduced frequency | `04_solver/run_case.py` |
| `runtime_environment.csv` | 1 | machine, library versions, CPU time of the march | `04_solver/run_case.py` |
| `cp_distribution_*.csv` | 2 | reconstructed surface pressure, local Mach number and temperatures at six instants | `04_solver/run_case.py` |
| `reconstruction_*.csv` | 2 | the reconstruction instant by instant as numbers | `04_solver/run_case.py` |
| `cycle_scan_*.csv` | 2 | the reconstruction at instants spread over the cycle | `04_solver/run_case.py` |
| `field_*_on_mesh.csv.gz` | 2 | the reconstruction at the nodes of the O-grid, gzip-compressed CSV | `02_mesh/field_on_mesh.py` |
| `field_*.csv.gz` | 8 | the reconstructed field on a rectangular grid at one instant, gzip-compressed CSV | `04_solver/run_case.py` |
| `convergence/residuals_*.csv` | 2 | change of the cycle with the number of cycles marched | `04_solver/run_case.py` |
| `convergence/timestep_*.csv` | 2 | loads against the time step; observed order and discretisation uncertainty | `04_solver/run_case.py` |
| `convergence/panel_convergence_*.csv` | 3 | the reconstruction against the number of panels | `04_solver/run_case.py` |
| `convergence/reconstruction_timing.csv` | 1 | CPU time of the reconstruction | `04_solver/run_case.py` |

## 06_postprocessing

| Files | n | What they are | Written by |
|---|---|---|---|
| `figure_record.csv` | 2 | resolution and smallest lettering of the folder's figures, from which the report works out the size of the lettering as printed | the figure script(s) of the folder |
| `validation/figure_record.csv` | 0 | the same for the two validation figures | `06_postprocessing/validation/make_validation.py` |
| `README.md` | 1 | standing of the figures | written by hand |
| `make_all_plots.py` | 1 | script of the two-dimensional figures | written by hand |
| `make_3d_plots.py` | 1 | script of the surfaces, maps over two variables and the vector view | written by hand |
| `overlay_record.csv` | 1 | whether the measured loop was drawn with each case | `06_postprocessing/make_all_plots.py` |
| `masked_region_record.csv` | 1 | the region beyond the critical pressure printed on each figure of a corrected quantity | `06_postprocessing/make_all_plots.py` |
| `plots/fig3d_*.png` | 5 | response surface; surfaces of the reconstructed pressure and speed; section with velocity vectors | `06_postprocessing/make_3d_plots.py` |
| `plots/cp_phase_map_*.png` | 2 | reconstructed upper-surface pressure over chord and cycle | `06_postprocessing/make_3d_plots.py` |
| `plots/contour_*.png` | 56 | maps of the reconstructed field: pressure, local Mach number, static and recovery temperature, speed, vorticity, vectors | `06_postprocessing/make_all_plots.py` |
| `plots/*.png` | 17 | load loops, time histories, model states, static inputs, convergence, surface pressure and surface temperatures | `06_postprocessing/make_all_plots.py` |
| `validation/make_validation.py` | 1 | script of the validation folder | written by hand |
| `validation/README.md` | 1 | what each file of the validation folder is | `06_postprocessing/validation/make_validation.py` |
| `validation/PROVENANCE.md` | 1 | where the measured frames come from | `06_postprocessing/validation/make_validation.py` |
| `validation/validation_*.csv` | 3 | static comparison, one row per scored loop, and the targets with their verdicts; the same numbers as results/ | `06_postprocessing/validation/make_validation.py` |
| `validation/loop_*.csv` | 12 | the model read at the measured points of the loops drawn | `06_postprocessing/validation/make_validation.py` |
| `validation/calibration_constants.csv` | 1 | the calibrated and fixed constants | `06_postprocessing/validation/make_validation.py` |
| `validation/fig_validation_*.png` | 2 | static comparison; measured against predicted loops | `06_postprocessing/validation/make_validation.py` |

## 07_report

| Files | n | What they are | Written by |
|---|---|---|---|
| `build_report.py` | 1 | the one builder of the Word file and the three PDFs | written by hand |
| `rpt_*.py` | 10 | the parts of the builder: data, text, layout, album, dossier | written by hand |
| `case.docx` | 1 | the report as a Word document | `07_report/build_report.py` |
| `UNISTALL_*.pdf` | 3 | the report, the plots album and the data dossier | `07_report/build_report.py` |
| `_equations/eq_*.png` | 24 | the equations of the report | `07_report/build_report.py` |
| `report_numbers.json` | 1 | every number quoted in the report with the file it was read from | `07_report/build_report.py` |
| `figure_lettering.csv` | 1 | smallest lettering of each figure of the report as printed | `07_report/build_report.py` |
| `page_fill.csv` | 1 | how full each page of the three PDFs is | `07_report/build_report.py` |

## 08_engineering_drawings

| Files | n | What they are | Written by |
|---|---|---|---|
| `figure_record.csv` | 1 | resolution and smallest lettering of the folder's figures, from which the report works out the size of the lettering as printed | the figure script(s) of the folder |
| `draw_engineering.py` | 1 | script of the drawing sheets | written by hand |
| `sheet?_*.png` | 4 | four illustrative sheets, not for manufacture | `08_engineering_drawings/draw_engineering.py` |
| `dimensions.csv` | 1 | every dimension printed on a sheet, with its source | `08_engineering_drawings/draw_engineering.py` |
| `sheet_record.csv` | 1 | what each sheet states: scale, viewing direction of section A-A, body coordinates | `08_engineering_drawings/draw_engineering.py` |

## Where to find what

For a reader who expects the layout of a flow-solver case:

| Looking for | Here |
|---|---|
| fluid or material properties | 03_model_setup/air_properties.csv holds the gas properties |
| aerofoil polar or tabulated section data | 03_model_setup/static_inputs.csv holds the static inputs the model reads; they come from measurement. 05_solution/model_static_polar_<case>.csv is the model's own slow-sweep polar |
| solver executable or solver script | the model is the package unistall/; 04_solver/unistall_solver.py runs it for one condition and 04_solver/run_case.py for the two cases |
| mesh and solution on the mesh | 02_mesh is an illustration only; no flow equation is solved. 05_solution/field_*_on_mesh.csv.gz is the reconstruction evaluated at its nodes |
| field or contour data | 05_solution/field_*.csv.gz, compressed (gzip CSV), read with pandas.read_csv |
| validation data and comparison | 06_postprocessing/validation, the same numbers as results/ |
| separate report, album and dossier builders | one builder, 07_report/build_report.py, writes the Word file and the three PDFs |

## The two case names

File names carry the keys `A_validation` and `B_application`. They are labels of the two conditions and nothing more. Case A is the tunnel condition of a measured loop that was used in the calibration, so the comparison drawn for it is a fit and not an independent test; Case B is a section in a steady stream at a blade station's Mach number and reduced frequency, not a rotor calculation. The comparison of the model with loops held out of the calibration is in `results/` and in `06_postprocessing/validation/`.
