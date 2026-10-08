"""Every figure and table must come from the scripts, not from a hand.
Author: Akosa Samuel Onyejekwe (independent)"""
import pandas as pd

from unistall.paths import FIGURES, RESULTS, ROOT


def test_every_figure_has_a_function_that_draws_it():
    from unistall import make_figures
    drawn = {f.__name__ + ".png" for f in make_figures.FIGURES_ALL}
    present = {p.name for p in FIGURES.glob("*.png")}
    assert present == drawn, (sorted(present - drawn), sorted(drawn - present))


def test_tables_are_generated_from_the_result_files():
    from unistall import make_tables
    assert (RESULTS/"tables.md").read_text(encoding="utf-8") == make_tables.build()


def test_every_target_has_a_measured_value():
    import pandas as pd
    board = pd.read_csv(RESULTS/"targets_scoreboard.csv")
    assert set(board.met) <= {"yes", "no"}, board[~board.met.isin(["yes", "no"])].target.tolist()


def test_loads_reproduce_the_stored_reference():
    """The model's loads have not changed since the reference was stored.
    On the machine that stored it the difference is zero (code_quality.csv);
    another processor and maths library differ in the last digits, so the
    bound here is the one that holds across machines."""
    from unistall import code_quality
    assert code_quality.regression_difference() <= 1e-7


def test_default_resolution_meets_the_convergence_target():
    import json
    import pandas as pd
    from unistall.paths import DATA
    target = json.load(open(DATA/"targets.json"))["numerical"]["loop_error_change_max"]
    summary = pd.read_csv(RESULTS/"convergence_summary.csv")
    assert (summary.worst_change_in_loop_error <= target).all()


def test_case_study_outputs_pass_their_checks():
    """The numbered folders hold what the current solver produces and agree with one another."""
    import check_case_study
    failed = [text for text, ok in check_case_study.run() if not ok]
    assert not failed, failed


def test_bibliography_holds_every_requested_work():
    """Every key the bibliography builder asks for or enters by hand is an
    entry of docs/references.bib, none is unresolved, and the work on a
    separation point taken from the static lift curve is among them."""
    import re
    from unistall import build_bibliography as bb
    from unistall.paths import DOCS
    bib = (DOCS/"references.bib").read_text(encoding="utf-8")
    keys = set(re.findall(r"@\w+\{([^,]+),", bib))
    wanted = {c[0] for c in bb.C} | {m["key"] for m in bb.MANUAL}
    assert wanted == keys, sorted(wanted ^ keys)
    res = pd.read_csv(DOCS/"bibliography_resolution.csv")
    assert res.resolved.all() and set(res.key) == keys
    assert "hansen2004" in keys
    assert "Hansen, Gaunaa" in (ROOT/"README.md").read_text()


def test_manifest_is_what_its_script_builds():
    """data/data_manifest.csv, which decides the loops used, equals the output
    of unistall/build_manifest.py from the inventory and the report tables."""
    from unistall import build_manifest
    assert build_manifest.differences_from_committed() == []


def test_bibliography_entries_are_complete():
    """Every journal article has a year, a volume and pages or an article
    number; every entry has an author, a title and a year; and no surname is
    left in the registry's lower-case form."""
    import re
    from unistall.paths import DOCS
    bib = (DOCS/"references.bib").read_text(encoding="utf-8")
    for m in re.finditer(r"@(\w+)\{([^,]+),\n(.*?)\n\}", bib, re.S):
        f = dict(re.findall(r"^\s*(\w+) = \{(.*)\},?$", m.group(3), re.M))
        for field in ("author", "title", "year"):
            assert f.get(field), (m.group(2), field)
        if m.group(1) == "article":
            assert f.get("volume") and f.get("pages"), m.group(2)
        assert not re.search(r"\bMc[a-z]", f["author"]), m.group(2)
        assert not f["title"].isupper(), m.group(2)


def test_tables_are_well_formed_markdown():
    """Every table of results/tables.md has the same number of columns in
    every row as in its header (a bare bar in a cell would add one), and no
    header is a raw column name with underscores."""
    import re
    lines = (RESULTS/"tables.md").read_text(encoding="utf-8").split("\n")

    def columns(line):
        return len(re.findall(r"(?<!\\)\|", line)) - 1
    i = 0
    while i < len(lines) - 1:
        if lines[i].startswith("|") and re.fullmatch(r"\|(---\|)+", lines[i + 1]):
            n = columns(lines[i])
            assert "_" not in lines[i], lines[i]
            j = i + 2
            while j < len(lines) and lines[j].startswith("|"):
                assert columns(lines[j]) == n, lines[j]
                j += 1
            i = j
        else:
            i += 1


def test_source_log_lists_the_parts_the_data_files_use():
    """data/data_sources.csv names, for each NASA document, the tables,
    figures and pages that the data files cite, as unistall/source_log.py
    collects them."""
    from unistall import source_log
    old = pd.read_csv(source_log.LOG, dtype=str)
    assert list(old.parts_used.fillna("")) == list(source_log.build().parts_used.fillna(""))
    assert not old.status.str.contains("title page").any()


def test_published_reference_is_what_its_script_gives():
    """The comparison with the published Leishman-Beddoes curves is recomputed
    from the cached tracings and equals the stored file (skipped without them)."""
    import pytest
    from unistall import published_reference as pr
    if not all((pr.CACHE/pr._name(q, c, k)).exists() for c in pr.CASES for q in ("CN", "CM") for k in ("EXP", "BL")):
        pytest.skip("published curves not fetched: run python3 -m unistall.published_reference")
    new, old = pr.compare(), pd.read_csv(RESULTS/"published_reference.csv")
    for col in ("loop_error_published_model", "loop_error_this_model", "rms_difference_between_models_over_range"):
        assert (new[col] - old[col]).abs().max() < 1e-4, col
