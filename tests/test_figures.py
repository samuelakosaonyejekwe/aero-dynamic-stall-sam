"""Every figure must be drawn from the result files as they are now, and must
keep to the rules of the figures: no black, nothing drawn outside what was
compared, every marked loop inside the mapped field.
Author: Akosa Samuel Onyejekwe (independent)"""
import numpy as np
import pytest
from PIL import Image

from unistall import dsmodel as dm
from unistall import make_figures as mf
from unistall import metrics as mt
from unistall.paths import FIGURES

NO_FRAMES = "measured frames not fetched: run python3 -m unistall.fetch_frames"
DARKEST_CHANNEL = 48        # a pixel whose red, green and blue are all below this counts as black; the ink is (31, 51, 80)


def _frames_present() -> bool:
    return (mt.FRAME_CACHE/"frame_12102.mat").exists()


def test_figure_records_the_files_it_was_drawn_from(tmp_path, monkeypatch):
    """The save helper stores the hash of every file read, refuses to write a
    figure that read none, and a later change of a file is found."""
    import matplotlib.pyplot as plt
    monkeypatch.setattr(mf, "FIGURES", tmp_path)
    monkeypatch.setattr(mf, "_READ", {})
    with pytest.raises(RuntimeError):
        mf._save(plt.figure(), "no_inputs")
    source = tmp_path/"source.csv"
    source.write_text("x\n1\n")
    mf._use(source)
    mf._save(plt.figure(), "one_input")
    png = tmp_path/"one_input.png"
    assert list(mf.figure_inputs(png).values()) == [mf.file_hash(source)]
    assert mf.stale_inputs(png) == []
    source.write_text("x\n2\n")
    assert len(mf.stale_inputs(png)) == 1


def test_figure_is_not_older_than_a_file_it_reads():
    """Each PNG carries the hashes of the result and data files it was drawn
    from; a figure whose file has changed since is stale. Files that are
    absent at the moment are not compared."""
    if not _frames_present():
        pytest.skip(NO_FRAMES)
    problems = {}
    for make in mf.FIGURES_ALL:
        png = FIGURES/f"{make.__name__}.png"
        if not png.exists():
            problems[make.__name__] = "not drawn"
        elif not mf.figure_inputs(png):
            problems[make.__name__] = "carries no record of its inputs"
        elif mf.stale_inputs(png):
            problems[make.__name__] = "older than " + ", ".join(mf.stale_inputs(png))
    assert not problems, "stale figures (python3 -m unistall.make_figures): " + "; ".join(
        f"{k}: {v}" for k, v in problems.items())


def test_figure_tables_come_from_the_present_inputs():
    """The two tables the figure script owns were computed from the present
    constants, static inputs and loops, on the present grid."""
    if not _frames_present():
        pytest.skip(NO_FRAMES)
    assert mf.stale_tables() == [], "computed from other inputs (python3 -m unistall.make_figures fig08 fig09 fig14)"


def test_figure_has_no_black():
    """No pixel of any figure is black or nearly so."""
    found = {}
    for png in sorted(FIGURES.glob("*.png")):
        with Image.open(png) as im:
            rgb = np.asarray(im.convert("RGB"))
        n = int((rgb.max(axis=2) < DARKEST_CHANNEL).sum())
        if n:
            found[png.name] = n
    assert not found, found


def test_figure_damping_map_field_covers_every_marked_loop():
    """The grid of the design-space maps reaches every measured loop that the damping map marks."""
    if not _frames_present():
        pytest.skip(NO_FRAMES)
    means, ks = mf.map_grid()
    loops = mf.map_loops()
    assert len(loops) > 0
    assert means.min() <= loops.alpha0_deg.min() and loops.alpha0_deg.max() <= means.max()
    assert ks.min() <= loops.k.min() and loops.k.max() <= ks.max()


def test_figure_mach_trends_draw_nothing_outside_the_compared_range():
    """Every point of the Mach-trend maps, and every loop marked on them, lies
    inside the span of the measured loops."""
    if not _frames_present():
        pytest.skip(NO_FRAMES)
    (k0, k1), (M0, M1) = mf.compared_ranges()
    t = mf.trend_field()
    assert len(t) > 0
    assert t.k.between(k0 - 1e-9, k1 + 1e-9).all() and t.M.between(M0 - 1e-9, M1 + 1e-9).all()
    for (a0, a1), _rows in t.groupby(["alpha_mean_deg", "alpha_amp_deg"]):
        loops = mf.trend_loops(a0, a1)
        assert loops.k.between(k0, k1).all() and loops.M.between(M0, M1).all()


def test_figure_cost_surface_grid_holds_the_calibrated_point():
    """The calibrated T_VL and T_p are themselves points of the cost surface."""
    tvl, tp = mf.surface_grid()
    c = dm.load_constants()
    assert np.isclose(tvl, c["Tvl"], atol=1e-6).any() and np.isclose(tp, c["Tp"], atol=1e-6).any()
