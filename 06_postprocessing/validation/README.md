**Standing: the held-out results here are the same numbers as those of `results/`, not a second assessment.**

# 06_postprocessing/validation

Every file is written by `make_validation.py` from the files of this solver. Run it from the repository root:
`PYTHONPATH=. python3 06_postprocessing/validation/make_validation.py`.

| File | What it is |
|---|---|
| `validation_static.csv`, `fig_validation_static.png` | the digitised static points at Mach 0.30 beside what the model reads at the same incidence: its static normal force (not lift) and static moment, and for drag its slow-sweep value from `results/limit_static_drag.csv` |
| `validation_loops_summary.csv` | one row per scored loop of `results/validation_table.csv` for the tabulated model: conditions, loop errors, peak-lift error, stall-incidence errors, measured and predicted damping |
| `validation_targets.csv` | `results/validation_targets.csv`, passed through: each target with the measured value and whether it is met, as that file gives it |
| `fig_validation_loops.png` | measured against predicted loops for the loops listed below |
| `loop_<frame>_CL.csv`, `_CM.csv`, `_CD.csv` | for those loops, the model read at the measured incidences, each point on its own stroke |
| `calibration_constants.csv` | the constants of `results/calibrated_constants.json`, fitted and fixed, with literature values and bounds |
| `PROVENANCE.md` | where the measured frames come from |
| `figure_record.csv` | resolution and smallest lettering of the two figures |

The loops of the figure are chosen by a fixed rule (`chosen_loops`), not picked by eye:

- `frame_10309`: held-out loop with the lowest loop error in lift
- `frame_7117`: held-out loop with the median loop error in lift
- `frame_10117`: held-out loop with the highest loop error in lift
- `frame_9302`: the CALIBRATION loop set beside Case A (not an independent test)

A target that `validation_targets.csv` gives as not met is not met; nothing here changes that.
