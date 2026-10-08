**Standing: the load loops show the model's prediction; every field map and surface-pressure figure shows a reconstruction, not a flow solution.**

# 06_postprocessing

`make_all_plots.py` and `make_3d_plots.py` draw the figures of `plots/` from `05_solution` and
`03_model_setup`. Run them from the repository root with `PYTHONPATH=.`.

- Load loops, time histories and states: the load model. The Case A loops carry the measured points of a
  calibration loop, taken at its own conditions, which the legend prints beside the nominal conditions of the
  curve; `overlay_record.csv` says whether the points were drawn.
- `cp_distribution_*`, `temperature_profile_*`, `cp_phase_map_*`, `contour_*` and the `fig3d_*` surfaces and
  vectors: the potential-flow reconstruction of `unistall/flowfield.py`. The separated region is not drawn (the
  model's separation point is marked on each map) and the vortex is a marker of the model's vortex clock that
  does not reproduce a suction footprint. Pressure, local Mach number and temperatures come from a linearised
  (Kármán-Tsien) correction of the incompressible pressure, which holds only where the local flow is
  subcritical: the region beyond the critical pressure is hatched, its size is printed on each figure and kept
  in `masked_region_record.csv`. Each colour scale spans the values it colours. The maps are in axes fixed to
  the section.
- `fig3d_response_surface.png`: the cycle damping in it is model output only; the figure prints how often the
  model gets the sign of the damping right on the held-out loops.

`figure_record.csv` holds the resolution and the smallest lettering of every figure, from which the report
works out the size of the lettering as printed.
