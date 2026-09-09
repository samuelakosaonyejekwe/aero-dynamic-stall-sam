"""
assets / make_social_preview.py
-------------------------------
Renders the repository's link-preview card, the 1280x640 image that {% seo %}
serves as og:image / twitter:image (see `image:` in _config.yml) and that
GitHub's own "Social preview" setting is uploaded from by hand.

  assets/social-preview.png   lossless master (the upload source)
  assets/social-preview.jpg   what the site links to

WHY THIS IS A SCRIPT.
The card used to be a hand-made image, and it was the one surface of this study
whose numbers traced to nothing. Under a heading reading "VALIDATION HEADLINE"
it showed a dynamic C_L,max of 2.07 and a moment break of -0.31. Neither is a
result this study reports: the matched validation point gives 1.912 and -0.236
(05_solution/metrics_A_validation.csv). 2.07 is the largest C_L,max anywhere in
the 36-point response-surface SWEEP, at a condition no reported case uses, and
also happens to be frame 9217's measured peak; -0.31 is that same held-out
frame's moment break. The third tile, "< 1 % static-polar error", is true but is
a CALIBRATION check -- the separation law is fitted to that polar -- which
README.md is careful to say is not the predictive claim, so it did not belong
under that heading either.

Every number on the card is now read from the artifacts the pipeline writes, in
the same way the report body and the README's quoted results are, so it cannot
drift again. run_all.py rebuilds it as a pipeline stage.

The design is unchanged: teal gradient, white section mark with an orange
leading-edge vortex, no black anywhere (the house rule the figures follow).
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
from project_meta import AUTHOR                      # single source of truth
import pandas as pd

W, H = 1280, 640

# Palette, sampled from the card this replaces so the branding is unchanged.
TEAL_L = np.array([12, 79, 79]) / 255.0             # gradient, left edge
TEAL_R = np.array([20, 119, 119]) / 255.0           # gradient, right edge
ORANGE = "#f4a13c"
WHITE = "#ffffff"
SUBTLE = "#d6e6e4"                                  # muted line, not grey-black

# ---------------------------------------------------------------- the numbers
mA = pd.read_csv(ROOT/"05_solution"/"metrics_A_validation.csv"
                 ).set_index("metric")["value"]
vs = pd.read_csv(ROOT/"06_postprocessing"/"validation"/"validation_realdata_summary.csv"
                 ).set_index("metric")["value"]
# The headline is the MATCHED VALIDATION POINT (Case A) and the held-out
# prediction error -- i.e. what the study validates -- not a sweep maximum and
# not the calibration check.
TILES = [
    (f"{float(mA['CL_max_dynamic']):.2f}",  "dynamic  C$_{l,max}$"),
    (f"−{abs(float(mA['CM_min(c/4)'])):.2f}", "C$_m$ c/4 break"),
    (f"{float(vs['mean |CLmax| error [%]']):.1f} %", "held-out peak-lift error"),
]

fig = plt.figure(figsize=(W/100.0, H/100.0), dpi=100)
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W); ax.set_ylim(H, 0)                # y down, so coords read as pixels
ax.axis("off")

# ---------------------------------------------------------------- background
grad = np.linspace(0, 1, W)[None, :, None]
ax.imshow((TEAL_L[None, None, :]*(1-grad) + TEAL_R[None, None, :]*grad),
          extent=(0, W, H, 0), aspect="auto", zorder=0)

# ---------------------------------------------------------------- section mark
CX, CY, R = 294.0, 320.0, 195.0
ax.add_patch(Circle((CX, CY), R, fill=False, ec=WHITE, lw=5.0, zorder=2))

def naca0012(n=200):
    """Upper and lower surface of the section, x/c in [0, 1]."""
    x = (1 - np.cos(np.linspace(0, np.pi, n)))/2.0
    yt = 5*0.12*(0.2969*np.sqrt(x) - 0.1260*x - 0.3516*x**2
                 + 0.2843*x**3 - 0.1015*x**4)
    return np.concatenate([x, x[::-1]]), np.concatenate([yt, -yt[::-1]])

CHORD_PX, INCID = 310.0, np.radians(14.0)           # chord length and nose-up set
xs, ys = naca0012()
xs, ys = xs*CHORD_PX, ys*CHORD_PX
xr = xs*np.cos(INCID) + ys*np.sin(INCID)
yr = -xs*np.sin(INCID) + ys*np.cos(INCID)
LEX, LEY = CX - 146.0, CY - 33.0                    # leading edge, inside the ring
ax.add_patch(Polygon(np.column_stack([LEX + xr, LEY - yr]), closed=True,
                     facecolor=WHITE, edgecolor=WHITE, lw=1.0, zorder=3))

# Leading-edge dynamic-stall vortex: the one orange element, an open spiral
# rolling up just above and ahead of the leading edge. It sits OFF the white
# section rather than on it -- drawn over the aerofoil the two shapes read as
# one smudge at thumbnail size, which is the size this card is usually seen at.
VX, VY, VR = LEX + 17.0, LEY - 47.0, 3.6            # centre and spiral pitch
th = np.linspace(0.30*np.pi, 3.10*np.pi, 240)
rr = 3.0 + VR*th
ax.plot(VX + rr*np.cos(th), VY - rr*np.sin(th),
        color=ORANGE, lw=6.5, solid_capstyle="round", zorder=4)
_a, _r = th[0], rr[0]
ax.annotate("", xytext=(VX + _r*np.cos(_a), VY - _r*np.sin(_a)),
            xy=(VX + _r*np.cos(_a) - 11, VY - _r*np.sin(_a) + 9),
            arrowprops=dict(arrowstyle="-|>", color=ORANGE, lw=0.1,
                            mutation_scale=24), zorder=5)

# ---------------------------------------------------------------- type
X0 = 560.0
ax.text(X0, 163, "UNISTALL", color=WHITE, fontsize=68, fontweight="bold",
        ha="left", va="center", zorder=6)
ax.text(X0, 238, "Dynamic-Stall Prediction", color=WHITE, fontsize=31,
        ha="left", va="center", zorder=6)
ax.plot([X0 - 2, X0 + 333], [277, 277], color=ORANGE, lw=5.0,
        solid_capstyle="butt", zorder=6)
ax.text(X0, 325, "Reduced-order UIBS solver   ·   NACA 0012", color=SUBTLE,
        fontsize=19.5, ha="left", va="center", zorder=6)
ax.text(X0, 363, "Validated vs. NASA dynamic-stall experiments", color=SUBTLE,
        fontsize=19.5, ha="left", va="center", zorder=6)
ax.text(X0, 424, AUTHOR, color=ORANGE, fontsize=24, fontweight="bold",
        ha="left", va="center", zorder=6)
ax.plot([X0 - 2, X0 + 48], [469, 469], color=ORANGE, lw=4.0,
        solid_capstyle="butt", zorder=6)
ax.text(X0 + 60, 469, "VALIDATION HEADLINE", color=ORANGE, fontsize=15,
        fontweight="bold", ha="left", va="center", zorder=6)

# The three tiles are laid out from their own captions rather than at fixed
# stations, so a longer number or caption cannot run into its neighbour.
_caps = [c for _, c in TILES]
_widths = [max(len(c), 8) for c in _caps]
_total = sum(_widths)
_avail = W - 40 - X0
_x = X0
for (val, cap), wgt in zip(TILES, _widths):
    _w = _avail*wgt/_total
    ax.text(_x + _w/2, 533, val, color=WHITE, fontsize=37, fontweight="bold",
            ha="center", va="center", zorder=6)
    ax.text(_x + _w/2, 581, cap, color=SUBTLE, fontsize=15,
            ha="center", va="center", zorder=6)
    _x += _w

png = HERE/"social-preview.png"
jpg = HERE/"social-preview.jpg"
fig.savefig(png, dpi=100, facecolor="none")
fig.savefig(jpg, dpi=100, facecolor="none", pil_kwargs={"quality": 92})
plt.close(fig)
print(f"[social] {png.name} + {jpg.name}: "
      + "  ".join(f"{v} ({c})" for v, c in
                  [(t[0], t[1].replace('$', '')) for t in TILES]))
