"""README.md must be exactly what unistall/make_readme.py writes from the result
files, so none of its numbers can drift, and must not contain wording the
results do not support.
Author: Akosa Samuel Onyejekwe (independent)"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# wording the results do not support
BANNED = ["™", "(TM)", "universal", "novel", "patent", "competitive", "certification",
          "validated", "unconditionally", "real-time", "strong agreement", "CFD-like"]


def test_readme_is_generated_from_the_result_files():
    from unistall import make_readme
    assert (ROOT/"README.md").read_text(encoding="utf-8") == make_readme.build()


def test_readme_has_no_unsupported_wording():
    text = (ROOT/"README.md").read_text(encoding="utf-8").lower()
    hits = [b for b in BANNED if b.lower() in text]
    assert not hits, hits


def test_title_and_affiliation_come_from_one_place():
    """The software's title is the same in README.md and CITATION.cff, and the
    affiliation is the one string of project_meta wherever one is given."""
    import project_meta as meta
    from unistall.paths import ROOT
    readme = (ROOT/"README.md").read_text(encoding="utf-8")
    cff = (ROOT/"CITATION.cff").read_text(encoding="utf-8")
    assert readme.startswith(f"# {meta.SOFTWARE_TITLE}\n")
    assert f'title: "{meta.SOFTWARE_TITLE}"' in cff
    assert f"affiliation: {meta.AFFILIATION}" in cff and f"({meta.AFFILIATION})" in readme
