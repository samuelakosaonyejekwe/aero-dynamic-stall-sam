"""The reconstruction's helpers that have an answer known in advance.
Author: Akosa Samuel Onyejekwe (independent)"""
import numpy as np

from unistall import flowfield as ff


def test_mask_growth_does_not_wrap_round_the_grid():
    """A mask on one edge of the grid, grown by a few rings, must leave the
    opposite edge clear; a single node grows to a square of side 2 rings + 1."""
    mask = np.zeros((12, 15), bool)
    mask[:, 0] = True
    grown = ff._dilate(mask, 3)
    assert grown[:, :4].all() and not grown[:, 4:].any()
    one = np.zeros((11, 11), bool)
    one[5, 5] = True
    grown = ff._dilate(one, 2)
    assert grown.sum() == 25 and grown[3:8, 3:8].all()
    assert (ff._dilate(one, 0) == one).all()
