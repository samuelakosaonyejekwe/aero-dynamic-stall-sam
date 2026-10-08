# 04_solver

The solver of this case study is the Python package [`unistall`](../unistall): the load model is
`unistall.dsmodel.solve`, and its calibrated constants are read from `results/calibrated_constants.json`.
Nothing in this folder computes a load.

`unistall_solver.py` is a command-line entry for one condition: it has no physics of its own and only calls
the package (`unistall.dsmodel` or, with `--form state-space`, `unistall.statespace`) and writes the cycle to a
CSV file, for example `PYTHONPATH=. python3 04_solver/unistall_solver.py 10 10 0.10 0.30 --out cycle.csv`.

`run_case.py` only drives the package for the two cases set up in `03_model_setup/` and writes the
results to `05_solution/`: the loads and states over the last cycle, the scalar results, the convergence
studies, the error bands of the model on the held-out loops and a response surface. Run it from the
repository root: `PYTHONPATH=. python3 04_solver/run_case.py`.

The flow fields and surface pressures it writes come from `unistall.flowfield`. They are a potential-flow
reconstruction drawn around the predicted lift, not a flow solution: incompressible with a linearised
(Kármán-Tsien) correction of the pressure that holds only where the local flow is subcritical, quasi-steady, with a
transpiration through the suction surface (of either sign) standing for the separated region and a vortex
that only marks the model's vortex clock. They are compared with no measurement and no load depends on them.
The settings under `field_reconstruction.read_by_the_driver` in `03_model_setup/solver_config.json` are the
ones the driver passes on; the others in that file are a record of the module's constants. The fields are
stored as gzip-compressed CSV.
