# Record of when the targets and the data split were set

Author: Akosa Samuel Onyejekwe (independent)

The repository is published as a single commit, so its history cannot show
that the targets and the split into calibration and held-out loops were set
before the held-out loops were scored. This page gives the dated record. The
author keeps the working history it comes from and will show it to the editor
or a referee on request; the identifiers below let that history be checked
against this page.

| Item | First recorded | Commit in the working history | SHA-256 of the file as first recorded |
|---|---|---|---|
| Targets | 2026-10-06 10:46 (+03:00) | `98a3443986413b194d93adf2e431cc5e36dbc53a` | `bbdce806f6e6f022086098d600e150504a0b9657f4a5a5184465c39cb731dca0` |
| Calibration and held-out split | 2026-10-06 10:46 (+03:00) | `98a3443986413b194d93adf2e431cc5e36dbc53a` | `7a408f7ca22d2cfad14d28c0323e8e03cfe583840bcc15ecd4ecf02528a6e80a` |
| First scoring of this model on the held-out loops | 2026-10-06 13:22 (+03:00) | `6425f007798cb229b284204269eebba336ae8d8b` | |

**The split.** `data/split.json` is byte for byte the file first recorded
(SHA-256 `7a408f7ca22d2cfad14d28c0323e8e03cfe583840bcc15ecd4ecf02528a6e80a`).
It fixes five NACA 0012 frames by name, the ones in use before the rule on the
frame number was written, and applies the rule to all others. Two of the five
are not where the rule would put them: frame 9302 is a calibration loop
(the rule gives held-out) and frame 7113 is a held-out loop (the rule gives
calibration). `results/split_departures.csv` gives the held-out measures with
those two moved to the rule's sets and with them left out.

**The targets.** `data/targets.json` is not byte for byte the first record:
the groups were renamed and the file was reduced to the limits that apply.
The values are these.

| Target | First record | `data/targets.json` now |
|---|---|---|
| Held-out lift loop error, at most | 0.10 | 0.10 |
| Held-out moment loop error, at most | 0.15 | 0.15 |
| Held-out drag loop error, at most | 0.20 | 0.20 |
| Held-out incidence of maximum lift, at most | 0.7238 deg | 0.7238 deg |
| Held-out incidence of moment stall, at most | 1.0 deg | 1.0 deg |
| Held-out damping sign correct, at least | 0.80 | 0.80 |
| Attached flow against Theodorsen, amplitude and phase | 3 %, 3 deg | 3 %, 3 deg |
| Attached-flow loops, lift and moment loop error | 0.05, 0.10 | 0.05, 0.10 |
| Drag negative for at most, and no lower than | 0.05 of the cycle, -0.01 | 0.05 of the cycle, -0.01 |
| Shift of stall onset with a finer step, at most | 0.02 deg | 0.02 deg |
| Repeatability of the calibration, at most | 1e-6 | 1e-6 |

Four things changed after the first record, and none loosens a limit:

1. The limit on the incidence of maximum lift was first recorded as 0.7238
   deg, the stricter of two values. For part of 2026-10-06 the published file
   carried the looser one, 1.0 deg. It is 0.7238 deg again; the measured value
   meets both.
2. The group `numerical` (change in loop error with a finer step and with
   more cycles, at most 0.002) was added on 2026-10-07, after a convergence
   study showed that the step had to be set in semichords of travel. Its
   limit was written down before the study was repeated with the new step.
3. The structural and code-quality groups were added during 2026-10-06; they
   do not bear on the comparison with measurement. The target on the shift
   of stall onset with a finer step was first recorded under a name that gave
   step counts (720 to 5760); since the step is now set in semichords it is
   named for what is compared, the default march against four times as many
   steps. Its limit, 0.02 deg, is unchanged.
4. The group `equivalence_of_separation_laws` (a margin of one tenth of each
   accuracy target for calling the two separation laws equivalent) was added
   on 2026-10-07. The differences between the two models' mean loop errors
   had been seen by then (at most 0.0021); the direct differences between
   the predicted loads had not been computed. The lift difference turned out
   to lie outside the margin and is reported so.

**Where each limit comes from.** The first record formed each held-out limit
by one rule: the stricter of an absolute limit and three quarters of the
baseline, the value that the earlier form of the model had scored on the same
held-out group at that time.

| Measure | Absolute limit | Baseline | 75 % of baseline | Limit |
|---|---|---|---|---|
| Lift loop error | 0.10 | 0.1475 | 0.1106 | 0.10 |
| Moment loop error | 0.15 | 0.2581 | 0.1936 | 0.15 |
| Drag loop error | 0.20 | 0.4043 | 0.3032 | 0.20 |
| Incidence of maximum lift, deg | 1.0 | 0.965 | 0.7238 | 0.7238 |
| Incidence of moment stall, deg | 1.0 | 1.475 | 1.1063 | 1.0 |

So the four-figure limit on the incidence of maximum lift is 0.75 x 0.965
deg. The absolute limits, and the limit of 0.80 on the damping sign, are the
author's working limits for a section-load model to be useful; they are not
taken from a standard. `data/targets.json` carries the same table under
`held_out_validation.origin`.

What this record cannot show: the held-out loops had been scored with an
earlier form of the load model before the targets were written, so the limits
were set knowing roughly how hard they were. They were not adjusted afterwards
to fit a result.
