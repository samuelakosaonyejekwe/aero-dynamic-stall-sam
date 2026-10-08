# -*- coding: utf-8 -*-
"""
unistall / paths.py
---------------
Where things are kept.

Author: Akosa Samuel Onyejekwe (independent)

  data/             inputs: verified conditions of the measured loops, digitised
                    static data, calibration and held-out sets, targets
  results/          everything the scripts write: constants, scores, checks
  results/figures/  every figure, written by unistall/make_figures.py
  docs/             the formulation and the reference list
  cache             downloads (measured loops, source reports); not committed

ROOT, the folder that holds data/, results/ and docs/, is found in this order:
  1. the environment variable UNISTALL_ROOT, if set;
  2. the folder above the package, when the package is run from a clone of
     the repository (or installed from one with `pip install -e .`);
  3. <prefix>/share/unistall, where `pip install .` puts a copy of data/,
     results/ and docs/ so the model and its calibrated constants can be used
     from anywhere. The scripts that write results are meant for a clone.

The cache is UNISTALL_CACHE if set, unistall/cache/ in a clone, and otherwise
~/.cache/unistall.
"""
import os
import sys
from pathlib import Path

PACKAGE = Path(__file__).resolve().parent


def _root() -> Path:
    """The folder that holds data/, results/ and docs/ (see the module text)."""
    if os.environ.get("UNISTALL_ROOT"):
        return Path(os.environ["UNISTALL_ROOT"]).resolve()
    candidates = [PACKAGE.parent, Path(sys.prefix)/"share"/"unistall"]
    for c in candidates:
        if (c/"data"/"targets.json").exists():
            return c
    raise FileNotFoundError("unistall cannot find its data/ and results/ folders; looked in "
                            + ", ".join(str(c) for c in candidates) + ". Set UNISTALL_ROOT to the folder that holds them.")


def _cache(root: Path) -> Path:
    """The download cache (see the module text)."""
    if os.environ.get("UNISTALL_CACHE"):
        return Path(os.environ["UNISTALL_CACHE"]).resolve()
    return PACKAGE/"cache" if root == PACKAGE.parent else Path.home()/".cache"/"unistall"


ROOT = _root()
IN_CLONE = ROOT == PACKAGE.parent
DATA = ROOT/"data"
RESULTS = ROOT/"results"
FIGURES = RESULTS/"figures"
DOCS = ROOT/"docs"
CACHE = _cache(ROOT)
