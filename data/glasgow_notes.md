# The Glasgow blind test: data, assumptions, protocol and what to run

Author: Akosa Samuel Onyejekwe (independent). Written 2026-10-07, before any
forecast existed and before any held-out run was scored.

Every number below is in a file written by a script; the file is named beside it.

## 1. What the data are

**Experiment.** University of Glasgow, 7 ft x 5 ft closed-return low-speed
tunnel; NACA 0012 ("model 11"), chord 0.55 m, spanning the 5 ft dimension,
pitched about the quarter chord by a hydraulic actuator; 30 surface pressure
transducers at mid-span; tests of March 1991. Source: Green, R. B. & Giuni, M.
(2017), *Dynamic stall database R and D 1570-AM-01: Final Report*, University
of Glasgow, doi:10.5525/gla.researchdata.464 (CC BY 4.0).

**Files used.** `unistall/fetch_glasgow.py` fetches, into the git-ignored
cache, (a) the 223 sinusoidal runs as redistributed in BL-DSM-JFS-2021,
directory "Glasgow Data", at commit 84e7945c5b16c6da4075c7fb467c01125dace36b,
three files a run, and (b) from the original deposit, by byte range, the final
report, the 28 static runs of the same model, the transducer positions, the
two run lists of model 11 and the two MATLAB sources that define the
coefficient files. Checksums: `data/glasgow_inventory.csv`,
`data/glasgow_static_inventory.csv`, `data/glasgow_sources.csv`.

**The sinusoidal runs** (`data/glasgow_inventory.csv`): 223 runs, all model 11,
clean (experiment type 0), motion type 1, each 128 samples that are the
average of 10 cycles. Mach 0.078 to 0.155 (191 of them at 0.1146 to 0.1192, 16
at 0.078, 16 at 0.155), chord Reynolds number 0.97 to 1.87 million, reduced
frequency 0.0097 to 0.175, measured mean incidence 2.5 to 20.3 deg, measured
amplitude 3.9 to 10.4 deg, peak incidence up to 29.9 deg, lowest incidence
-6.7 deg. The database holds 414 runs of this model (final report, Table 1);
the redistributing repository carries only the 223 clean sinusoidal ones. No
ramp run is used.

**What each run file holds** (final report sections 2 and 4, Tables 7 and 12;
`DSplot.m`, `export_coeffs_Callback`):

- `NNNNNNNN.dat`: a 32-number run information block (date, temperature,
  barometric pressure, nominal mean and amplitude, frequency, dynamic
  pressure, Reynolds number, Mach number, reduced frequency k = omega c/(2U),
  wind speed), then per sample the dynamic pressure [Pa], 30 pressure
  coefficients (transducer 1 at the upper-surface trailing edge, 30 at the
  lower-surface trailing edge) and the incidence [deg].
- `NNNNNNNN_coeffs.dat`: phase omega t [rad], incidence [deg], C_n, C_t, C_m.
  C_n is normal to the chord; C_t is along the chord, positive towards the
  leading edge; C_m is about the quarter chord, nose-up positive. All three
  are sums of pressure over the 30 transducers.
- `GUD_NNNNNNNN.mat`: the same table repacked by the repository's
  `generate_GUD.m`.

**Checks made on the files** (columns of the inventory, all true for all 223
runs): the incidence column is the same in the two text files; the `.mat`
arrays equal the coefficient file; each redistributed `NNNNNNNN.dat` has the
CRC-32 that the deposit's archive records for its own file of that name, so
the pressure files are the deposited ones unchanged. The coefficient files
are not in the deposit; they were exported by whoever built the repository.
`unistall/glasgow_static.py` recomputes them from the pressures by the
deposit's panel sum: on the 76 calibration-type runs the largest difference
is 0.0041 in C_n, 0.0008 in C_t and 0.0030 in C_m
(`results/glasgow_integration_check.csv`). The held-out coefficient files
were not compared; `score` checks their SHA-256 against the inventory.

**Derived here:** C_l = C_n cos a + C_t sin a; C_d = C_n sin a - C_t cos a;
measured mean = (max + min)/2 and amplitude = (max - min)/2 of the incidence
column.

## 2. Assumptions

1. No wall or blockage correction is in the files and none is applied. The
   final report names none. The deposit's folder name "CorrectedFiles" is
   taken to mean the test records corrected against the paper reports
   (section 5 of the report), not a wall correction. The model is 0.55 m in a
   2.13 m high section (chord/height 0.26), so the uncorrected interference
   is not small.
2. The incidence column is the geometric pitch angle. The measured motion is
   treated as a pure sine at the measured mean and amplitude; the nominal
   values differ from the measured ones by up to 0.8 deg in the mean and 0.5 deg in amplitude.
3. Mach number, reduced frequency, Reynolds number and wind speed are taken
   as stored. The speed of sound used for a run is its stored speed divided
   by its stored Mach number.
4. The section is the NACA 0012 of the report's formula, with the open
   trailing edge; the panel sum closes the outline at x = 1.
5. The measured zero-lift incidence of the static runs (+0.28 to +0.36 deg)
   is an offset of the incidence scale and is carried as such, as for the
   NASA data.
6. The attached-flow moment factor of `results/attached_moment_factor.json`
   (from NASA calibration loops) is kept. It acts in proportion to the
   compressible weight, so it has no effect at Mach 0.078 and little at 0.117.

## 3. What makes this test weaker than it looks

- **Drag is pressure drag.** No skin friction is in C_t, so the measured C_d
  is low by the friction drag everywhere and is near zero (even slightly
  negative) in attached flow. The model's C_D0 is therefore set to the static
  pressure drag at zero lift from the static runs (0.003 to 0.005), not to a
  handbook value. The drag target is the hardest to read: the normalising
  range of C_d is small on attached loops.
- **Thirty transducers.** The chord force comes from few points near the
  leading edge; the suction peak is under-resolved at high incidence.
- **No wall correction** (assumption 1), in a tunnel in which the model is
  large.
- **Second tunnel, but not a second calibration-free test.** The dynamic
  constants are refitted on Glasgow calibration-type runs. What is blind is
  the held-out two thirds of the runs, not the tunnel.
- **Calibration-type and held-out runs are neighbours.** The rule on the run
  number interleaves them: a held-out run usually has a calibration-type run
  at the same mean incidence and frequency with a different amplitude, or the
  next mean incidence. The blind error measures interpolation inside the test
  matrix.
- **Below Mach 0.20 the model is outside the range it was assessed in.** The
  attached-flow loads are a blend of the compressible indicial form and the
  incompressible one (`attached_flow.py`): compressible weight 0 at Mach
  0.078, 0.17 at 0.117, 0.55 at 0.155. At these Mach numbers the loads are
  therefore mostly Jones' approximation of Wagner's function with the
  apparent-mass terms, the impulsive time constants play almost no part, and
  the stall constants fitted at Mach 0.3 were never fitted here.
- **Static stations differ in Reynolds number as much as in Mach number**
  (1.0, 1.5 and 1.9 million); Mach number is only the label the interpolation
  uses. The outer two stations rest on one static run each.
- **Reynolds number 1 to 2 million, no trip.** Static stall and reattachment
  at this Reynolds number depend on transition; the static runs show it
  (section 5).
- **Peak incidence above 25 deg** (8 calibration-type and 16 held-out runs):
  beyond the static data; scored and reported, outside the targets.
- **The repository the files come from tuned its own model on these runs**,
  all of them. Nothing of that model or its constants is used here.

## 4. The split

`data/glasgow_split.json`, fixed 2026-10-07: run number divisible by 3 ->
calibration, otherwise held_out; no assignment by name. It is the rule
`data/split.json` gives for frames added later. Result (column `set` of the
inventory): **76 calibration-type, 147 held out**. In the primary group (Mach
inside the static range, peak incidence at or below 25 deg): **68
calibration-type** (4 at Mach 0.078, 59 near 0.117, 5 at 0.155) and **131
held out** (10, 111, 10). The 28 static runs are static inputs and are in
neither set.

## 5. Static inputs, and how far the two sources differ

`unistall/glasgow_static.py`; rule as `static_stations.py` and
`static_model.py` (slope and zero-lift incidence from a line through the up
stroke between -5 and 8 deg; static stall where the up-stroke C_l is largest;
C_N1 the normal force there; f by the Kirchhoff inverse with the sweep's own
slope; the "more separated" curve is the up stroke to stall and the down
stroke beyond; moment and drag of the same stroke).

**Source used by the model: the static runs of the same model**, the 9 dated
from 1 March 1991 on, the days of the sinusoidal tests
(`results/glasgow_static_stations.csv`):

| Mach | Re (million) | slope /deg | zero-lift deg | stall deg | C_l,max | C_N1 | pressure C_D0 | runs |
|---|---|---|---|---|---|---|---|---|
| 0.078 | 0.98 | 0.1018 | +0.34 | 13.33 | 1.144 | 1.126 | 0.0047 | 1 |
| 0.118 | 1.47 | 0.1046 | +0.36 | 14.48 | 1.337 | 1.310 | 0.0026 | 7 |
| 0.155 | 1.89 | 0.1081 | +0.28 | 14.63 | 1.427 | 1.397 | 0.0031 | 1 |

Scatter over the 7 runs at Mach 0.118: standard deviation 0.14 deg in stall
incidence and 0.010 in C_N1.

**Second source: the slowest calibration-type sinusoidal runs read as
sweeps** (k = 0.0098, runs 11012322 and 11013021, Mach 0.116). Less the
static source at Mach 0.118 (`results/glasgow_static_difference.csv`):
slope -0.4 %, zero-lift incidence -0.09 deg, **stall incidence +1.06 deg,
C_l,max +0.061, C_N1 +0.057 (+4.4 %)**, pressure C_D0 -0.003; separation
point f differs by up to 0.53 (at 15.5 deg, where one source has stalled and
the other has not), 0.12 RMS over 0 to 25 deg; static C_m by up to 0.12. A
run at k = 0.01 is not static: stall is delayed by about a degree. This
source is not used by the model.

**Third comparison: static runs of the earlier weeks** (17 runs, December
1990 to January 1991) less those used: C_N1 +2.1 % at Mach 0.118, +5.4 % at
0.078, +3.4 % at 0.16; stall incidence +0.16 to +0.54 deg; f up to 0.2 to 0.6
different in the stall break.

**Uncertainty to carry:** about 1 deg in static-stall incidence and 4 to 5 %
in C_N1, from either comparison. The forecast is repeated with the
quasi-steady source on the calibration-type runs to show what that does to
the loop errors; no held-out run is scored with it.

**A published polar.** No published static polar for the NACA 0012 at
Reynolds number 1 to 2 million was verified against its source for this
work, so none is used and no number from one appears here. The static runs
above are measurements on the same model in the same tunnel, which is the
better static input; what they cannot show is a bias common to this tunnel
(wall interference, three-dimensional effects at stall). The redistributing
repository holds a file labelled as Abbott and von Doenhoff static lift at
Reynolds number 3 million; it was not fetched or checked.

## 6. The protocol

In plain words, in the order it is carried out:

1. **Split.** Before any model result existed for these runs, each
   sinusoidal run was assigned by its run number alone: divisible by three,
   calibration; otherwise held out. The rule and its date are in
   `data/glasgow_split.json`.
2. **Targets.** The measures and limits are those already fixed for the NASA
   held-out loops (`data/targets.json`), copied unchanged into
   `data/glasgow_targets.json` together with the whole procedure, dated, with
   the statement that no held-out run had been scored.
3. **Forecast.** Using only calibration-type runs, the dynamic constants are
   fitted by the procedure of `unistall/calibrate.py` inside a five-fold
   cross-validation: each run is predicted by constants fitted without it.
   The mean errors of those left-out predictions are the forecast of the
   blind result.
4. **Go rule.** Written before the forecast existed and binding: go if the
   forecast mean loop errors in lift, moment and drag are each within 1.24
   times their targets (0.124, 0.186, 0.248); otherwise no-go. The factor
   comes from the NASA results as they stood: there the held-out error was
   0.94 of the cross-validated forecast, and the 95 % interval of a mean loop
   error was up to 16 % of the mean either side; 1.16/0.94 = 1.24.
5. **Seal.** If the verdict is go, the checksums of everything that fixes the
   test (constants, static inputs, targets, split, inventories, forecast,
   and the source files of the model, the metrics and the test itself) are
   written to `results/glasgow_seal.json` with the date and the forecast. No
   seal can be written on a no-go.
6. **Single scoring.** The held-out runs are scored once. The scoring
   refuses to run without the seal, if any sealed file has changed, or if a
   score already exists. The result is reported against the targets whatever
   it is. On a no-go the held-out runs are never scored, and the case is
   reported as a stated limit of the model, with the forecast numbers.

Throughout, the measured loads of a held-out run can be read only inside the
scoring: `fetch_glasgow.load_run` and `glasgow_blind.load` raise otherwise.

## 7. What to run next, in this order

Two cores only: prefix each with `taskset -c 0,1`, one at a time. The module
uses 2 worker processes.

1. `python3 -m unistall.glasgow_blind forecast`
   writes `results/glasgow_forecast.csv`, `glasgow_forecast_summary.csv`,
   `glasgow_verdict.json`. Cost: each fold fits at the search march (0.21
   million steps of the march per evaluation of all 68 runs) and the
   left-out runs are scored at the reporting resolution (1.9 million steps
   for the 68 runs, three candidates). The eight runs at k = 0.01 dominate
   the reporting cost.
2. `python3 -m unistall.glasgow_blind forecast --static=quasi_steady`
   (optional, calibration-type runs only; files with suffix `_quasi_steady`;
   takes no part in the rule).
3. Read `results/glasgow_verdict.json`.
   - **no-go:** stop. Do not run `constants`, `seal` or `score`. Report the
     forecast numbers and the limit.
   - **go:** continue.
4. `python3 -m unistall.glasgow_blind constants` -> `results/glasgow_constants.json`.
5. `python3 -m unistall.glasgow_static` must not be rerun between here and the
   scoring unless its outputs are unchanged (they are sealed).
6. `python3 -m unistall.glasgow_blind seal` -> `results/glasgow_seal.json`.
   Commit the seal before scoring, so its date is on record.
7. `python3 -m unistall.glasgow_blind score` -> `results/glasgow_blind_scores.csv`,
   `results/glasgow_blind_summary.csv`. Once.

Do not edit `unistall/dsmodel.py`, `attached_flow.py`, `static_model.py`,
`metrics.py`, `strokes.py`, `glasgow_static.py`, `fetch_glasgow.py` or
`glasgow_blind.py` between the seal and the scoring: the scoring will refuse.
A change of the model before the seal is allowed, but then the forecast must
be rerun, since the seal takes the verdict from the forecast file.

Not done by these modules, left for after the scoring: figures, the tables of
`make_tables.py`, the entries of `data_sources.csv` and `NOTICE` for the
Glasgow deposit, and listing the new modules in `code_quality.PACKAGE`.
