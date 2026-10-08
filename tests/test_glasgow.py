"""The Glasgow blind test: inventory, split, targets, the guard on held-out
runs and the seal. No network is used and no held-out run is scored; the
tests that read the cached files are skipped where the cache is absent.
Author: Akosa Samuel Onyejekwe (independent)"""
import io
import json

import pandas as pd
import pytest

from unistall import fetch_glasgow as fg
from unistall import glasgow_blind as gb
from unistall.paths import DATA

needs_cache = pytest.mark.skipif(not (fg.RUNS/"11011962.dat").exists() or not (fg.DEPOSIT/fg.CRC_NAME).exists(),
                                 reason="Glasgow files not fetched: python3 -m unistall.fetch_glasgow")


def _through_csv(df: pd.DataFrame) -> pd.DataFrame:
    return pd.read_csv(io.StringIO(df.to_csv(index=False)))


def _inventory() -> pd.DataFrame:
    return pd.read_csv(fg.INV)


# ---- inventory ---------------------------------------------------------------
@needs_cache
def test_inventory_is_what_the_fetch_script_builds_from_the_cache():
    pd.testing.assert_frame_equal(_inventory(), _through_csv(fg.inventory_from_cache()))


@needs_cache
def test_static_inventory_and_sources_are_what_the_fetch_script_builds():
    pd.testing.assert_frame_equal(pd.read_csv(fg.STATIC_INV), _through_csv(fg.static_inventory_from_cache()))
    pd.testing.assert_frame_equal(pd.read_csv(fg.SOURCES), _through_csv(fg.sources_from_cache()))


def test_inventory_holds_conditions_and_counts_only():
    inv = _inventory()
    assert not [c for c in inv.columns if c.lower().startswith(("cl_", "cm_", "cd_", "cn_", "ct_"))]
    assert (inv.model == fg.MODEL_NUMBER).all() and (inv.motion_type == fg.MOTION_SINUSOIDAL).all()
    assert inv.dat_crc32_equals_deposit.all() and inv.mat_equals_coeffs.all()
    assert (inv[["n_cn", "n_ct", "n_cm"]].values == inv[["n_pressure_rows"]].values).all()


# ---- split -------------------------------------------------------------------
def test_split_is_the_rule_on_the_run_number():
    rule = json.load(open(fg.SPLIT))
    assert rule["fixed_on"] == "2026-10-07"
    assert rule["rule"]["modulus"] == 3 and rule["rule"]["calibration_remainder"] == 0
    assert rule["rule"]["assignments_by_name"] == []
    inv = _inventory()
    expected = ["calibration" if run % 3 == 0 else "held_out" for run in inv.run]
    assert list(inv["set"]) == expected
    assert all(fg.split_of(run) == s for run, s in zip(inv.run, expected, strict=True))


def test_no_run_is_in_both_sets():
    inv = _inventory()
    assert inv.run.is_unique
    cal, held = set(inv[inv["set"] == "calibration"].run), set(inv[inv["set"] == "held_out"].run)
    assert not cal & held and len(cal) + len(held) == len(inv)
    assert set(inv["set"]) == {"calibration", "held_out"}


# ---- targets -----------------------------------------------------------------
def test_targets_are_those_of_the_nasa_assessment():
    ours = json.load(open(gb.TARGETS))
    nasa = json.load(open(DATA/"targets.json"))["held_out_validation"]
    mine = ours["held_out_validation"]
    # the limits themselves; the description of the group and the note on where
    # the NASA limits came from are not limits
    notes = ("group", "origin")
    assert {k: v for k, v in mine.items() if k not in notes} == {k: v for k, v in nasa.items() if k not in notes}
    assert ours["frozen_on"] == "2026-10-07" and ours["go_rule"]["status"] == "binding"
    assert "No held-out Glasgow run had been scored" in ours["statement"]
    assert ours["go_rule"]["factor"] == ours["go_rule"]["basis"]["factor"]


def test_go_rule_is_each_loop_error_within_the_factor_of_its_target():
    targets = json.load(open(gb.TARGETS))
    f, held = targets["go_rule"]["factor"], targets["held_out_validation"]
    at = {q: f*held[f"mean_{q}_max"] for q in gb.RULE_MEASURES}
    inside = pd.DataFrame({q: [0.98*v, 1.0*v] for q, v in at.items()})
    assert gb.verdict(inside, targets)["go"]
    for q in gb.RULE_MEASURES:
        outside = inside.copy()
        outside[q] = 1.02*at[q]
        v = gb.verdict(outside, targets)
        assert not v["go"] and not v["measures"][f"mean_{q}"]["within"]


# ---- the guard ---------------------------------------------------------------
def _one_of(set_name: str) -> int:
    inv = _inventory()
    return int(inv[inv["set"] == set_name].run.iloc[0])


def test_a_held_out_run_cannot_be_loaded_outside_the_scoring():
    run = _one_of("held_out")
    with pytest.raises(fg.HeldOutAccess):
        fg.load_run(run)
    for purpose in ("forecast", "constants", "seal", "score"):
        with pytest.raises(fg.HeldOutAccess):
            gb.load(run, purpose)
    assert fg._HELD_OUT_OPEN is False


def test_opening_the_scoring_does_not_admit_another_purpose():
    run = _one_of("held_out")
    with gb._scoring_open():
        with pytest.raises(fg.HeldOutAccess):
            gb.load(run, "forecast")
    assert fg._HELD_OUT_OPEN is False


@needs_cache
def test_a_calibration_run_loads_for_the_forecast():
    fr = gb.load(_one_of("calibration"), "forecast")
    assert len(fr["cl"]) == len(fr["cm"]) == len(fr["cd"]) == 128


def test_the_forecast_and_the_fit_take_calibration_runs_only():
    runs = gb.runs_table("calibration")
    assert (runs.run % 3 == 0).all() and set(runs.group) <= {"primary", "beyond_static_lift_range"}


# ---- the seal ----------------------------------------------------------------
def _fixture(tmp_path, scale: float) -> dict:
    """Stand-in files for a seal: three sealed files and a forecast whose
    chosen rows sit at `scale` times the limits of the go rule."""
    targets = json.load(open(gb.TARGETS))
    f, held = targets["go_rule"]["factor"], targets["held_out_validation"]
    rows = pd.DataFrame({q: [scale*f*held[f"mean_{q}_max"]]*4 for q in gb.RULE_MEASURES})
    rows["chosen"], rows["group"] = [True, True, False, False], "primary"
    rows.loc[~rows.chosen, list(gb.RULE_MEASURES)] = 9.0          # rows of a candidate not kept take no part
    forecast = tmp_path/"forecast.csv"
    rows.to_csv(forecast, index=False)
    files = {}
    for name in ("constants", "static", "source"):
        files[name] = tmp_path/f"{name}.txt"
        files[name].write_text(name)
    files["forecast"] = forecast
    return dict(files=files, forecast_file=forecast, targets_file=gb.TARGETS, seal_file=tmp_path/"seal.json")


def test_no_seal_is_written_on_a_no_go(tmp_path):
    fx = _fixture(tmp_path, 1.05)
    with pytest.raises(gb.NoGo):
        gb.seal(**fx)
    assert not fx["seal_file"].exists()
    with pytest.raises(gb.SealError):
        gb.check_seal(fx["files"], fx["seal_file"], tmp_path/"scores.csv")


def test_seal_records_checksums_and_is_written_once(tmp_path):
    fx = _fixture(tmp_path, 0.95)
    record = gb.seal(**fx)
    assert record["verdict"]["go"] and record["scored_on"] is None
    assert record["sha256"] == gb.checksums(fx["files"]) and set(record["sha256"]) == set(fx["files"])
    assert json.load(open(fx["seal_file"])) == record
    with pytest.raises(gb.SealError):
        gb.seal(**fx)
    assert gb.check_seal(fx["files"], fx["seal_file"], tmp_path/"scores.csv") == record


def test_scoring_is_refused_if_a_sealed_file_changed(tmp_path):
    fx = _fixture(tmp_path, 0.95)
    gb.seal(**fx)
    fx["files"]["constants"].write_text("changed")
    with pytest.raises(gb.SealError, match="constants"):
        gb.check_seal(fx["files"], fx["seal_file"], tmp_path/"scores.csv")


def test_scoring_is_refused_if_a_sealed_file_is_missing(tmp_path):
    fx = _fixture(tmp_path, 0.95)
    gb.seal(**fx)
    fx["files"]["static"].unlink()
    with pytest.raises(gb.SealError):
        gb.check_seal(fx["files"], fx["seal_file"], tmp_path/"scores.csv")


def test_scoring_is_refused_a_second_time(tmp_path):
    fx = _fixture(tmp_path, 0.95)
    record = gb.seal(**fx)
    scores = tmp_path/"scores.csv"
    scores.write_text("run\n")
    with pytest.raises(gb.SealError, match="once"):
        gb.check_seal(fx["files"], fx["seal_file"], scores)
    scores.unlink()
    record["scored_on"] = "2026-10-08"
    fx["seal_file"].write_text(json.dumps(record))
    with pytest.raises(gb.SealError, match="once"):
        gb.check_seal(fx["files"], fx["seal_file"], scores)


def test_the_real_seal_covers_every_file_that_fixes_the_test():
    files = gb.sealed_files()
    for name in ("constants", "targets", "split", "inventory", "static_inventory", "forecast"):
        assert name in files
    for src in ("dsmodel.py", "attached_flow.py", "static_model.py", "metrics.py", "glasgow_static.py",
                "glasgow_blind.py", "fetch_glasgow.py"):
        assert f"unistall/{src}" in files
    assert any(n.startswith("results/glasgow_static") for n in files)


def test_no_blind_score_exists_without_a_seal():
    assert gb.SEAL.exists() or not gb.SCORES.exists()


# ---- static inputs -----------------------------------------------------------
@needs_cache
def test_panel_sum_reproduces_the_stored_coefficients_of_a_calibration_run():
    import numpy as np

    from unistall import glasgow_static as gs
    run = _one_of("calibration")
    _, rows = fg.read_dat(fg.RUNS/f"{run}.dat")
    f, co = gs.forces(rows), fg.read_coeffs(fg.RUNS/f"{run}_coeffs.dat")
    for q, j in (("cn", 2), ("ct", 3), ("cm", 4)):
        assert np.abs(f[q] - co[:, j]).max() < 0.01


@needs_cache
def test_static_inputs_read_no_held_out_run():
    from unistall import glasgow_static as gs
    assert (gs.quasi_steady_runs()["set"] == "calibration").all()
    for source in gs.SOURCES:
        sm = gs.GlasgowStatic(0.117, source=source)
        assert 0.09 < sm.slope_per_deg < 0.12 and 12.0 < sm.alpha_stall_deg < 17.0 and 1.0 < sm.CN1 < 1.6
        assert sm.f_at(2.0) > 0.99 and sm.f_at(24.0) < 0.3


# ---- the summary of a scoring (made-up scores, no run solved) -----------------
def test_each_target_is_marked_met_or_not_met_from_the_group_means():
    import numpy as np

    from unistall import metrics as mt
    held = json.load(open(gb.TARGETS))["held_out_validation"]
    good = {"nRMS_CL": 0.05, "nRMS_CM": 0.10, "nRMS_CD": 0.30, "RMS_CL": 0.0, "RMS_CM": 0.0, "RMS_CD": 0.0,
            "CLmax_err_pct": 1.0, "dalpha_CLmax_deg": -0.5, "dalpha_Mstall_deg": 2.0, "Xi_sign_agree": 1.0,
            "CD_neg_frac": 0.0}
    table = pd.DataFrame([good]*5)
    summary = mt.summarise(table, {"primary": np.ones(len(table), bool)})
    met = gb.targets_table(summary, held).set_index("metric").met.to_dict()
    assert met == {"mean_nRMS_CL": True, "mean_nRMS_CM": True, "mean_nRMS_CD": False,
                   "mean_abs_dalpha_CLmax_deg": True, "mean_abs_dalpha_Mstall_deg": False, "mean_Xi_sign_agree": True}


def test_forecast_and_fit_use_two_processes():
    assert gb.PROCESSES == 2
