"""The stored held-out results must be the ones the code and data reproduce.
Author: Akosa Samuel Onyejekwe (independent)"""
import json
from pathlib import Path
import pytest
from unistall import metrics as mt

PKG = Path(__file__).resolve().parents[1]/"unistall"
from unistall.paths import RESULTS, ROOT
pytestmark = pytest.mark.skipif(not (mt.FRAME_CACHE/"frame_12102.mat").exists(),
                                reason="measured frames not fetched: run python3 -m unistall.fetch_frames")


@pytest.fixture(scope="module")
def primary_summary():
    from unistall import validate as v
    tab, _ = v.score_all()
    return v.targets_table(v.summaries(tab)).set_index("measure")


def test_held_out_measures_are_the_recorded_ones(primary_summary):
    """The held-out measures recomputed from the model and the fetched frames
    equal the recorded ones: the result files are current and reproducible.
    This is not a limit on the results; a deliberate change of the model is
    followed by a new record."""
    rec = json.load(open(RESULTS/"acceptance_thresholds.json"))
    for name, r in rec["measures"].items():
        got = float(primary_summary.measured[name])
        assert abs(got - r["recorded"]) <= rec["tolerance"], f"{name}: recomputed {got}, recorded {r['recorded']}"
        assert bool(primary_summary.met[name]) == r["met"], name


def test_separation_laws_with_identical_constants():
    """With the same constants the two separation laws score alike against
    measurement (every loop-error interval inside the margin of
    data/targets.json), and the verdict on the loads themselves is the one
    the stored differences give against that margin."""
    import json
    import pandas as pd
    margin = json.load(open(ROOT/"data"/"targets.json"))["equivalence_of_separation_laws"]
    e = pd.read_csv(RESULTS/"equivalence.csv")
    assert len(e) == 9
    for r in e.itertuples():
        m = margin[f"margin_{r.coefficient}"]
        assert r.margin == m
        assert -m < r.loop_error_ci95_lo and r.loop_error_ci95_hi < m
        assert bool(r.loads_within_margin) == (r.ci95_hi < m)
    per = pd.read_csv(RESULTS/"load_difference_per_frame.csv")
    assert per.groupby("constants").frame.nunique().eq(e.n_frames.iloc[0]).all()


def test_readme_states_the_load_difference_wherever_it_compares_the_laws():
    """The README may not say the two laws give the same loads; it must give
    the measured load difference and the verdict against the margin."""
    import pandas as pd
    text = " ".join((ROOT/"README.md").read_text().split())
    e = pd.read_csv(RESULTS/"equivalence.csv")
    finding = text[text.index("**1. Tabulated and fitted"):text.index("**2. What separate")]
    assert "same loads" not in finding
    cl = e[e.coefficient == "CL"]
    assert f"{cl.mean_rms_difference_over_range.max():.4f}" in finding
    assert ("not shown to be equivalent in lift" in finding) == (not cl.loads_within_margin.all())


def test_split_statement_names_every_frame_that_departs_from_the_rule():
    """Every NACA 0012 frame whose set differs from the frame-number rule is
    one fixed by name in data/split.json, and the README names each."""
    import pandas as pd
    man = pd.read_csv(ROOT/"data"/"data_manifest.csv")
    man = man[man.airfoil == "NACA 0012"]
    by_rule = man.frame.str.split("_").str[1].astype(int).mod(3).eq(0).map({True: "calibration", False: "held_out"})
    departing = set(man.frame[by_rule != man["set"]])
    fixed = pd.read_csv(RESULTS/"split_fixed_frames.csv")
    assert departing == set(fixed.frame[fixed.departs_from_rule])
    line = next(ln for ln in (ROOT/"README.md").read_text().splitlines() if ln.startswith("- The split into"))
    for frame in fixed.frame:
        assert frame.split("_")[1] in line
    for frame in departing:
        assert f"{frame.split('_')[1]} is in the" in line


def test_unsteady_moment_factor_is_reported_with_what_it_corrects():
    """The factor's result file carries, for every loop it was taken from, the
    part of the damping that comes from the lagged static table, the factor
    the loop would choose alone and whether the model stalls on it; with the
    factor the damping of each closed loop is within 0.02 of the measured;
    and the range is carried into the uncertainty summary."""
    import pandas as pd
    a = json.load(open(RESULTS/"attached_moment_factor.json"))
    lo, hi = a["basis"]["interval_95"]
    assert lo <= a["cm_unsteady"] <= hi and a["basis"]["n_loops"] >= 10
    for r in a["frames"]:
        for key in ("damping_static_part", "damping_unlagged_static", "best_factor_this_loop",
                    "peak_vortex_normal_force", "peak_CN_prime_over_CN1"):
            assert r[key] is not None
        if r["moment_loop_closed"]:
            assert abs(r["damping_with_factor"] - r["damping_measured"]) <= 0.02
    unc = pd.read_csv(RESULTS/"uncertainty_summary.csv")
    assert {"moment factor_min", "moment factor_max"} <= set(unc.columns)
    text = " ".join((ROOT/"README.md").read_text().split())
    assert "never reach static stall" not in text and "lagged reading of the static table" in text
    assert f"{lo:.2f} to {hi:.2f}" in text
