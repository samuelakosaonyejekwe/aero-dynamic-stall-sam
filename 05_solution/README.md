**Standing: the loads are the prediction of the load model; the fields and surface pressures are a reconstruction, not a solution of the flow equations.**

# 05_solution

Written by `04_solver/run_case.py` (and `field_*_peak_on_mesh.csv.gz` by `02_mesh/field_on_mesh.py`).

- `time_history_<case>.csv`, `metrics_<case>.csv`, `summary_all_cases.csv`, `response_surface.csv`: the load
  model's prediction. `held_out_error_bands.csv` gives the error of the model in each headline quantity on the
  held-out measured loops, and `static_station_spread.csv` what the blend of static inputs between the two
  Mach stations is worth in Case B.
- `cp_distribution_<case>.csv`, `field_<case>_*.csv.gz`, `reconstruction_<case>.csv`, `cycle_scan_<case>.csv`:
  the potential-flow reconstruction of `unistall/flowfield.py`, drawn round the predicted lift. It is
  incompressible, with a linearised (Kármán-Tsien) correction of the pressure from which the local Mach number
  and the static and recovery temperatures follow by the isentropic relations; where the corrected pressure is
  below the critical value (column `beyond_critical`) those values are not physical. It is quasi-steady, it has
  no boundary layer, no separated shear layer and no wake, its surface pressure does not return the model's
  moment or drag, and it is compared with no measurement. No load depends on it.
- `convergence/`: cycle-to-cycle change, refinement of the time step with the observed order and the
  discretisation uncertainty, the reconstruction against the number of panels, and its CPU time.

Case B is an aerofoil at a rotor-blade-station condition in a steady stream with prescribed pitch. It is not
a rotor calculation.
