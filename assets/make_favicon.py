# -*- coding: utf-8 -*-
"""
assets / make_favicon.py
------------------------
Draws favicon.ico at the repository root: a NACA 0012 outline at incidence on
a white tile with a navy border.

Run from the repository root:  PYTHONPATH=. python3 assets/make_favicon.py

Author: Akosa Samuel Onyejekwe (independent)
"""
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

import project_meta as pm
from unistall import naca4
from unistall.style import INK, PALETTE

ROOT = Path(__file__).resolve().parent.parent
SIZE = 256              # drawing size in pixels; the icon sizes are reductions of it
INCIDENCE_DEG = 12.0


def outline():
    """The section outline, pitched nose-up about mid-chord, in tile pixels."""
    X, Y, _, _ = naca4.outline(pm.SECTION, 80)
    px, py = X - 0.5, Y
    a = np.radians(INCIDENCE_DEG)
    rx, ry = px*np.cos(a) + py*np.sin(a), -px*np.sin(a) + py*np.cos(a)
    return [(SIZE*(0.5 + 0.8*u), SIZE*(0.5 - 0.8*v)) for u, v in zip(rx, ry, strict=True)]


def main():
    tile = Image.new("RGB", (SIZE, SIZE), "white")
    draw = ImageDraw.Draw(tile)
    draw.rectangle([0, 0, SIZE - 1, SIZE - 1], outline=INK, width=14)
    draw.polygon(outline(), fill=PALETTE[0], outline=INK)
    tile.save(ROOT/"favicon.ico", sizes=[(16, 16), (32, 32), (48, 48)])
    print("[favicon] written to %s" % (ROOT/"favicon.ico"))


if __name__ == "__main__":
    main()
