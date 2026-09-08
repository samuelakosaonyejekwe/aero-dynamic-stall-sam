"""
08_engineering_drawings / draw_engineering.py
---------------------------------------------
Professional, DIMENSIONED engineering drawings for the UNISTALL(TM) dynamic-stall
case study reference rotorcraft ("CS-MUH reference", a generic medium utility
helicopter).  Pure matplotlib (numpy allowed).  All output PNG @ dpi=150.

Sheets produced (landscape, 11 x 8.5 in / US-Letter proportion -- the working
"paper mm" grid below is 271.8 x 210, NOT the 297 x 210 of A4; the stated
drawing scales are exact on a sheet of that size):
  sheet1_general_arrangement_3view.png  Third-angle orthographic three-view.
  sheet2_isometric.png                  Isometric pictorial with envelope dims.
  sheet3_main_rotor_blade.png           Main-rotor blade plan + edge view, sec A-A.
  sheet4_section_AA_airfoil.png         NACA 0012 section A-A at r/R = 0.75.

HARD RULES honoured:
  * No black anywhere -- all line work / text in INK (#1f3350) or INK_SOFT.
  * Drafting conventions: border + title block, third-angle symbol, dimension
    lines with arrowheads + extension lines, value centred above the line,
    "ALL DIMENSIONS IN mm", general / tolerance notes.
  * No overlapping text: dims live in clear margins, leaders point outward.
"""
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import (Rectangle, Circle, Ellipse, Polygon, FancyBboxPatch,
                                Arc, PathPatch)
from matplotlib.path import Path as MplPath

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
from aero_style import apply_style, INK, INK_SOFT, PALETTE, GRID   # noqa: E402

apply_style()

# ----------------------------------------------------------------------------
# Sheet geometry (landscape, working in "paper mm": W x H proportional to fig)
# ----------------------------------------------------------------------------
FIG_W, FIG_H = 11.0, 8.5            # inches
H = 210.0                            # paper height units
W = H * FIG_W / FIG_H               # paper width units -> matches figure aspect
M = 8.0                             # outer margin -> border frame

from project_meta import (STUDY_DATE_ISO as DATE,      # single source of truth
                          AUTHOR_FULL as DRAWN_BY, AUTHOR_BAND)
PROJECTION = "THIRD ANGLE"

# ----------------------------------------------------------------------------
# CS-MUH reference aircraft — ONE definition, read by every sheet.
#
# Sheets 1 and 2 used to carry independent copies of these, and once the overall
# length was derived from each sheet's own drawn geometry rather than quoted,
# they disagreed: 20705 mm on sheet 1 against 20055 mm on sheet 2, because the
# two sheets put the hub and the tail rotor at different stations. Deriving
# everything from this one block is what makes the two sheets agree.
#
# The blade chord is read from 03_model_setup, the study's single source of
# truth for the case conditions, rather than restated as 527 here.
# ----------------------------------------------------------------------------
import pandas as pd
_flowB = pd.read_csv(ROOT/"03_model_setup"/"flow_conditions.csv"
                     ).set_index("parameter")["case_B_application"]

SPEC = dict(
    L_fus   = 15500.0,        # fuselage length, nose to tail        [mm]
    Wf      = 2360.0,         # fuselage width
    Hh      = 3760.0,         # rotor hub height above ground
    Hf      = 4000.0,         # fin tip height above ground
    Rdia    = 16360.0,        # main-rotor diameter  (= 2 * 8.18 m)
    trdia   = 3350.0,         # tail-rotor diameter
    track   = 2700.0,         # main-wheel track
    wbase   = 4830.0,         # wheelbase
    stab    = 4000.0,         # horizontal-stabiliser span
    n_blades = 4,
    hub_from_nose      = 4650.0,          # 0.30 of the fuselage length
    tailrotor_from_nose = 15500.0 - 700.0,  # on the fin, 700 forward of the tail
    blade_twist_deg    = -13.0,
    root_cutout_frac   = 0.20,
    section_station_frac = 0.75,
)
SPEC["R"] = SPEC["Rdia"]/2.0
SPEC["blade_chord"] = float(_flowB["chord_c"])*1000.0        # m -> mm
# overall length: forward-most point (rotor disc) to aft-most (tail-rotor disc)
SPEC["fwd_of_nose"] = max(0.0, SPEC["R"] - SPEC["hub_from_nose"])
SPEC["aft_of_tail"] = SPEC["tailrotor_from_nose"] + SPEC["trdia"]/2.0 - SPEC["L_fus"]
SPEC["L_ovl"] = SPEC["fwd_of_nose"] + SPEC["L_fus"] + max(0.0, SPEC["aft_of_tail"])

ARROW_KW = dict(arrowstyle="<->", color=INK, lw=0.8, mutation_scale=7,
                shrinkA=0, shrinkB=0)
LEAD_KW  = dict(arrowstyle="->",  color=INK_SOFT, lw=0.7, mutation_scale=8,
                shrinkA=0, shrinkB=1)

FS_DIM = 7.5      # dimension text
FS_NOTE = 7.0     # notes / leaders
FS_VIEW = 9.0     # view titles

# opaque halo so leader / centre / extension lines never strike through a value
TEXT_BG = dict(boxstyle="square,pad=0.18", fc="white", ec="none")


# ----------------------------------------------------------------------------
# Generic drafting helpers
# ----------------------------------------------------------------------------
def new_sheet():
    fig = plt.figure(figsize=(FIG_W, FIG_H))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.grid(False)
    return fig, ax


def draw_border(ax):
    # double border frame
    ax.add_patch(Rectangle((M, M), W - 2 * M, H - 2 * M,
                           fill=False, ec=INK, lw=1.6))
    ax.add_patch(Rectangle((M + 2, M + 2), W - 2 * M - 4, H - 2 * M - 4,
                           fill=False, ec=INK_SOFT, lw=0.6))


def third_angle_symbol(ax, x0, y0, h=8.0):
    """Truncated-cone two-view symbol, drawn inside box anchored bottom-left (x0,y0)."""
    d = h * 0.9
    r1, r2 = h * 0.45, h * 0.25         # big / small end radii
    cy = y0 + h / 2
    # left view: circle pair (front of cone) -> two concentric circles
    cx1 = x0 + h * 0.55
    ax.add_patch(Circle((cx1, cy), r1, fill=False, ec=INK, lw=0.9))
    ax.add_patch(Circle((cx1, cy), r2, fill=False, ec=INK, lw=0.9))
    # right view: trapezoid (side of truncated cone) with centre line
    cx2 = x0 + h * 1.45
    half_w = d * 0.7
    ax.add_patch(Polygon([(cx2 - half_w, cy - r1), (cx2 + half_w, cy - r2),
                          (cx2 + half_w, cy + r2), (cx2 - half_w, cy + r1)],
                         closed=True, fill=False, ec=INK, lw=0.9))
    # centre lines
    ax.plot([cx1 - r1 * 1.3, cx2 + half_w * 1.3], [cy, cy],
            color=INK_SOFT, lw=0.5, dashes=(6, 2, 1, 2))
    ax.text(x0 + h * 1.0, y0 - 2.2, "THIRD ANGLE", color=INK,
            fontsize=6.0, ha="center", va="top")


def title_block(ax, title, dwg_no, scale, material):
    """Title block in bottom-right corner inside border."""
    bw, bh = 96.0, 40.0
    x1 = W - M - 2 - bw
    y1 = M + 2
    x2, y2 = x1 + bw, y1 + bh
    ax.add_patch(Rectangle((x1, y1), bw, bh, fill=False, ec=INK, lw=1.0))

    # row lines
    r = [y1, y1 + 8, y1 + 16, y1 + 24, y2]   # 4 rows
    for yy in r[1:-1]:
        ax.plot([x1, x2], [yy, yy], color=INK_SOFT, lw=0.5)
    # column split
    xc = x1 + bw * 0.62
    ax.plot([xc, xc], [y1, r[3]], color=INK_SOFT, lw=0.5)

    def cell(xa, ya, key, val, kfs=5.6, vfs=7.4, vw="bold"):
        ax.text(xa + 1.4, ya + 5.8, key, color=INK_SOFT, fontsize=kfs,
                ha="left", va="center")
        ax.text(xa + 1.4, ya + 2.4, val, color=INK, fontsize=vfs,
                ha="left", va="center", fontweight=vw)

    # top band: author + title (spans full width, tall)
    ax.text(x1 + 2, r[3] + 11.5, AUTHOR_BAND, color=INK,
            fontsize=9.5, ha="left", va="center", fontweight="bold")
    ax.text(x1 + 2, r[3] + 5.0, title, color=INK,
            fontsize=7.6, ha="left", va="center", fontweight="bold")

    # left column cells
    cell(x1, r[2], "DRAWN BY", DRAWN_BY, vfs=5.3)
    cell(x1, r[1], "DATE", DATE)
    cell(x1, r[0], "MATERIAL / SECTION", material, vfs=6.4)
    # right column cells
    cell(xc, r[2], "DWG No.", dwg_no)
    cell(xc, r[1], "SCALE", scale)
    cell(xc, r[0], "UNITS / PROJ.", "mm  /  3rd ANGLE", vfs=6.0)

    return x1, y2   # for placing projection symbol just above


def notes_block(ax, x, y, lines, title="NOTES:"):
    ax.text(x, y, title, color=INK, fontsize=FS_NOTE + 0.5,
            ha="left", va="top", fontweight="bold")
    for i, ln in enumerate(lines):
        ax.text(x, y - 4.0 - i * 3.4, ln, color=INK, fontsize=FS_NOTE,
                ha="left", va="top")


def _ext_v(ax, x, y_from, y_to, gap=0.8, over=1.6):
    d = 1 if y_to > y_from else -1
    ax.plot([x, x], [y_from + d * gap, y_to + d * over],
            color=INK_SOFT, lw=0.5)


def _ext_h(ax, y, x_from, x_to, gap=0.8, over=1.6):
    d = 1 if x_to > x_from else -1
    ax.plot([x_from + d * gap, x_to + d * over], [y, y],
            color=INK_SOFT, lw=0.5)


def dim_h(ax, x1, x2, ydim, yref, text, fs=FS_DIM):
    """Horizontal linear dimension. Extension lines from yref up/down to ydim."""
    _ext_v(ax, x1, yref, ydim)
    _ext_v(ax, x2, yref, ydim)
    ax.annotate("", xy=(x2, ydim), xytext=(x1, ydim), arrowprops=ARROW_KW)
    ax.text((x1 + x2) / 2, ydim + 0.9, text, color=INK, fontsize=fs,
            ha="center", va="bottom", bbox=TEXT_BG, zorder=8)


def dim_v(ax, y1, y2, xdim, xref, text, fs=FS_DIM, side="left"):
    """Vertical linear dimension. Text rotated 90, placed on `side` of the line."""
    _ext_h(ax, y1, xref, xdim)
    _ext_h(ax, y2, xref, xdim)
    ax.annotate("", xy=(xdim, y2), xytext=(xdim, y1), arrowprops=ARROW_KW)
    dx = -0.9 if side == "left" else 0.9
    ha = "right" if side == "left" else "left"
    ax.text(xdim + dx, (y1 + y2) / 2, text, color=INK, fontsize=fs,
            ha="center", va="center", rotation=90, bbox=TEXT_BG, zorder=8)


def leader(ax, xy, xytext, text, ha="left", va="center", fs=FS_NOTE):
    ax.annotate(text, xy=xy, xytext=xytext, color=INK, fontsize=fs,
                ha=ha, va=va, arrowprops=LEAD_KW, bbox=TEXT_BG, zorder=8)


def view_title(ax, x, y, text):
    ax.text(x, y, text, color=INK, fontsize=FS_VIEW, ha="center",
            va="center", fontweight="bold", bbox=TEXT_BG, zorder=8)


def centre_mark(ax, cx, cy, r):
    """Centre cross-hair (long-dash) for circular features."""
    ax.plot([cx - r, cx + r], [cy, cy], color=INK_SOFT, lw=0.5,
            dashes=(7, 2, 1.5, 2))
    ax.plot([cx, cx], [cy - r, cy + r], color=INK_SOFT, lw=0.5,
            dashes=(7, 2, 1.5, 2))


def common_note(ax, x=M + 4, y=H - M - 4):
    ax.text(x, y, "ALL DIMENSIONS IN mm UNLESS STATED",
            color=INK, fontsize=FS_NOTE + 0.5, ha="left", va="top",
            fontweight="bold")


def finalise(ax, title, dwg_no, scale, material):
    draw_border(ax)
    bx, by = title_block(ax, title, dwg_no, scale, material)
    third_angle_symbol(ax, bx + 28, by + 5)
    common_note(ax)


# ============================================================================
# SHEET 1 -- General arrangement, third-angle three-view
# ============================================================================
def sheet1():
    fig, ax = new_sheet()
    S = 1.0 / 270.0          # paper units per real mm  (scale 1:270)

    # every dimension from SPEC; the overall length is DERIVED there, not quoted.
    # It used to be given as 19760, which matched nothing that is drawn: with the
    # hub 4650 mm aft of the nose the rotor disc reaches 3530 mm AHEAD of it, so
    # the "OVERALL" dimension started aft of the aircraft's forward-most point.
    L_fus, Wf = SPEC["L_fus"], SPEC["Wf"]
    Hh, Hf = SPEC["Hh"], SPEC["Hf"]
    Rdia, trdia = SPEC["Rdia"], SPEC["trdia"]
    track, wbase, stab = SPEC["track"], SPEC["wbase"], SPEC["stab"]
    HUB_FRAC = SPEC["hub_from_nose"]/L_fus
    fwd_of_nose, aft_of_tail = SPEC["fwd_of_nose"], SPEC["aft_of_tail"]
    L_ovl = SPEC["L_ovl"]
    rP = half = Rdia * S / 2          # rotor-disc radius in paper units

    # ---- FRONT VIEW (looking on the nose) -------------------------------
    cxF = 52.0                  # lateral centreline
    gyF = 98.0                  # ground line
    front_view(ax, cxF, gyF, S, Wf, Hh, Hf, Rdia, track, trdia, stab)
    view_title(ax, cxF, gyF - 25.5, "FRONT VIEW")

    # front dims below part (stacked: width, track, rotor span)
    dim_h(ax, cxF - Wf * S / 2, cxF + Wf * S / 2, gyF - 7.5, gyF, f"{int(Wf)}")
    dim_h(ax, cxF - track * S / 2, cxF + track * S / 2, gyF - 12.5, gyF,
          f"{int(track)}  TRACK")
    dim_h(ax, cxF - half, cxF + half, gyF - 18.5, gyF + Hh * S,
          f"Ø{int(Rdia)}  ROTOR DISC")
    # heights on the RIGHT of the front view (clear of border / side view)
    dim_v(ax, gyF, gyF + Hh * S, cxF + half + 7.0, cxF + half, f"{int(Hh)}",
          side="right")
    dim_v(ax, gyF, gyF + Hf * S, cxF + half + 15.0, cxF + half,
          f"{int(Hf)}", side="right")

    # ---- PLAN VIEW (top, placed ABOVE front -> third angle) -------------
    cyP = 160.0                 # vertical centre of plan
    plan_view(ax, cxF, cyP, S, L_fus, Wf, Rdia, trdia, stab)
    view_title(ax, cxF, 117.0, "PLAN VIEW")

    # rotor diameter (vertical dim on the left)
    dim_v(ax, cyP - rP, cyP + rP, cxF - rP - 5.0, cxF - rP, f"Ø{int(Rdia)}")
    # stabiliser span (horizontal, below plan)
    # keep the dimension line BELOW the tail-rotor disc so its value never sits
    # on top of the disc outline (tail rotor spans down to cyP - rP + ~1.5)
    ytail = cyP - L_fus * S * 0.55 + L_fus * S * 0.07
    dim_h(ax, cxF - stab * S / 2, cxF + stab * S / 2, cyP - rP - 8.5, ytail,
          f"{int(stab)}  STAB SPAN")
    # main rotor leader -> clear zone right of the disc
    leader(ax, (cxF + 0.62 * rP, cyP + 0.62 * rP),
           (cxF + rP + 6, cyP + rP - 6),
           f"4-BLADE MAIN ROTOR\nNACA 0012  Ø{int(Rdia)}", ha="left")

    # ---- SIDE VIEW (right side, placed to the RIGHT of front) -----------
    cxS = 150.0                 # nose-to-tail mid
    side_view(ax, cxS, gyF, S, L_fus, L_ovl, Hh, Hf, Rdia, trdia, wbase, HUB_FRAC)
    view_title(ax, cxS, gyF - 25.5, "SIDE VIEW (PORT)")

    # side dims: overall length, fuselage length, wheelbase below; height right
    _nose = cxS - L_fus * S / 2
    dim_h(ax, _nose - fwd_of_nose * S, _nose + (L_fus + aft_of_tail) * S,
          gyF - 18.5, gyF, f"{int(round(L_ovl))}  OVERALL")
    dim_h(ax, cxS - L_fus * S / 2, cxS + L_fus * S / 2, gyF - 12.5, gyF,
          f"{int(L_fus)}  FUSELAGE")
    nose_x = _nose
    dim_h(ax, nose_x, nose_x + wbase * S, gyF - 7.0, gyF,
          f"{int(wbase)}  W/BASE")
    dim_v(ax, gyF, gyF + Hh * S, cxS + L_ovl * S / 2 + 5.0,
          cxS + L_fus * S / 2, f"{int(Hh)}", side="right")
    leader(ax, (cxS + L_fus * S / 2 - 1.0, gyF + Hf * S - 1.0),
           (cxS + L_ovl * S / 2 + 8.0, gyF + Hf * S + 8.0),
           f"TAIL ROTOR\nØ{int(trdia)}", ha="left")

    notes_block(ax, M + 4, 70,
                ["1.  OUTLINE INDICATIVE; ROUNDS R VARY.",
                 "2.  GENERAL TOL. ±25 mm; ANGULAR ±0.5°.",
                 f"3.  MAIN ROTOR {SPEC['n_blades']} BLADES, Ø{int(Rdia)}.",
                 f"4.  TAIL ROTOR {SPEC['n_blades']} BLADES, Ø{int(trdia)}, PORT.",
                 f"5.  OVERALL LENGTH {int(round(L_ovl))} INCL ROTOR.",
                 "6.  PROJECTION: THIRD ANGLE."])

    finalise(ax, "GENERAL ARRANGEMENT - CS-MUH REFERENCE",
             "UN-CSMUH-001", "1:270", "GA / OML")
    fig.savefig(HERE / "sheet1_general_arrangement_3view.png", dpi=150,
                facecolor="white")
    plt.close(fig)


def front_view(ax, cx, gy, S, Wf, Hh, Hf, Rdia, track, trdia, stab):
    bw = Wf * S
    body_top = gy + 2500 * S
    clear = 700 * S                       # wheel ground clearance
    # fuselage cross-section (rounded body)
    ax.add_patch(FancyBboxPatch((cx - bw / 2, gy + clear), bw, body_top - gy - clear,
                                boxstyle="round,pad=0,rounding_size=" + str(bw * 0.32),
                                fill=True, fc=PALETTE[0], ec=INK, lw=1.1, alpha=0.10))
    ax.add_patch(FancyBboxPatch((cx - bw / 2, gy + clear), bw, body_top - gy - clear,
                                boxstyle="round,pad=0,rounding_size=" + str(bw * 0.32),
                                fill=False, ec=INK, lw=1.1))
    # cockpit windscreen line
    ax.plot([cx - bw * 0.30, cx + bw * 0.30],
            [body_top - 600 * S, body_top - 600 * S], color=INK_SOFT, lw=0.7)
    # mast + hub
    hub_y = gy + Hh * S
    ax.plot([cx, cx], [body_top, hub_y], color=INK, lw=1.4)
    ax.add_patch(Rectangle((cx - bw * 0.10, hub_y - 1.2), bw * 0.20, 1.4,
                           fill=True, fc=INK_SOFT, ec=INK, lw=0.8))
    # rotor disc edge-on (full diameter line) + slight droop tips
    half = Rdia * S / 2
    ax.plot([cx - half, cx + half], [hub_y, hub_y], color=INK, lw=1.3)
    ax.plot([cx - half, cx - half * 0.0, cx + half],
            [hub_y - 1.0, hub_y, hub_y - 1.0], color=INK_SOFT, lw=0.6)
    # vertical fin behind (dashed) up to fin tip
    fin_y = gy + Hf * S
    ax.plot([cx, cx], [hub_y, fin_y], color=INK_SOFT, lw=0.8, dashes=(4, 2))
    # horizontal stabiliser (line) lower on tail, behind
    sy = gy + 1900 * S
    ax.plot([cx - stab * S / 2, cx + stab * S / 2], [sy, sy],
            color=INK_SOFT, lw=0.9, dashes=(4, 2))
    # tail rotor edge-on (port side, small vertical line near fin top)
    ax.plot([cx - bw * 0.65, cx - bw * 0.65], [fin_y - trdia * S / 2, fin_y + trdia * S / 2],
            color=INK_SOFT, lw=0.9, dashes=(4, 2))
    # landing gear: two main wheels at track
    wr = 320 * S
    for sgn in (-1, 1):
        wx = cx + sgn * track / 2 * S
        ax.add_patch(Circle((wx, gy + wr), wr, fill=True, fc=PALETTE[4],
                            ec=INK, lw=0.9, alpha=0.20))
        ax.add_patch(Circle((wx, gy + wr), wr, fill=False, ec=INK, lw=0.9))
        ax.plot([wx, cx + sgn * bw * 0.30], [gy + 2 * wr, gy + clear],
                color=INK, lw=0.8)
    # ground line
    ax.plot([cx - half - 4, cx + half + 4], [gy, gy], color=INK, lw=1.0)
    for k in range(-6, 7):
        gx = cx + k * (half + 4) / 6.0
        ax.plot([gx, gx - 1.4], [gy, gy - 1.4], color=INK_SOFT, lw=0.5)


def plan_view(ax, cx, cy, S, L_fus, Wf, Rdia, trdia, stab):
    # rotor disc (true circle)
    rP = Rdia * S / 2
    ax.add_patch(Circle((cx, cy), rP, fill=False, ec=INK_SOFT, lw=0.8))
    centre_mark(ax, cx, cy, rP * 1.05)
    # fuselage (teardrop) -- nose up (+y), tail down (-y)
    hw = Wf * S / 2
    yN = cy + L_fus * S * 0.45     # nose
    yT = cy - L_fus * S * 0.55     # tail
    boom_hw = hw * 0.28
    pts = [(0, yN), (hw * 0.75, yN - L_fus * S * 0.12),
           (hw, cy + L_fus * S * 0.18), (hw, cy - L_fus * S * 0.05),
           (boom_hw, cy - L_fus * S * 0.32), (boom_hw, yT + L_fus * S * 0.02),
           (0, yT)]
    full = [(cx + dx, y) for dx, y in pts] + [(cx - dx, y) for dx, y in reversed(pts)]
    ax.add_patch(Polygon(full, closed=True, fill=True, fc=PALETTE[0],
                         ec=INK, lw=1.1, alpha=0.10))
    ax.add_patch(Polygon(full, closed=True, fill=False, ec=INK, lw=1.1))
    # cabin partition
    ax.plot([cx - hw * 0.9, cx + hw * 0.9],
            [cy + L_fus * S * 0.05, cy + L_fus * S * 0.05],
            color=INK_SOFT, lw=0.6)
    # 4 blades
    hubr = 700 * S
    for ang in (0, 90, 180, 270):
        a = np.deg2rad(ang)
        bx = cx + (rP - 2) * np.cos(a)
        by = cy + (rP - 2) * np.sin(a)
        ax.plot([cx + hubr * np.cos(a), bx], [cy + hubr * np.sin(a), by],
                color=INK, lw=2.0, solid_capstyle="round")
    ax.add_patch(Circle((cx, cy), hubr, fill=True, fc=INK_SOFT, ec=INK, lw=0.8))
    # horizontal stabiliser near tail
    sy = yT + L_fus * S * 0.07
    ax.add_patch(Rectangle((cx - stab * S / 2, sy - 1.0), stab * S, 2.0,
                           fill=True, fc=PALETTE[2], ec=INK, lw=0.8, alpha=0.18))
    ax.add_patch(Rectangle((cx - stab * S / 2, sy - 1.0), stab * S, 2.0,
                           fill=False, ec=INK, lw=0.8))
    # tail rotor disc (port side) at tail
    trx = cx - boom_hw - trdia * S / 2 - 0.5
    ax.add_patch(Circle((trx, yT + L_fus * S * 0.02), trdia * S / 2,
                        fill=False, ec=INK_SOFT, lw=0.8))


def side_view(ax, cx, gy, S, L_fus, L_ovl, Hh, Hf, Rdia, trdia, wbase, hub_frac=0.30):
    # profile polygon -- nose at left (-x), tail at right (+x)
    nose = cx - L_fus * S / 2
    L = L_fus * S
    def P(fx, fy):     # fractions of length / absolute height-mm
        return (nose + fx * L, gy + fy * S)
    top = [P(0.02, 1600), P(0.06, 2350), P(0.20, 2550), P(0.34, 2500),
           P(0.55, 2150), P(0.78, 1900), P(0.92, 1820), P(1.00, 1780)]
    belly = [P(1.00, 1500), P(0.85, 1450), P(0.62, 1250), P(0.40, 760),
             P(0.20, 700), P(0.08, 800), P(0.02, 1150)]
    poly = top + belly
    ax.add_patch(Polygon(poly, closed=True, fill=True, fc=PALETTE[0],
                         ec=INK, lw=1.1, alpha=0.10))
    ax.add_patch(Polygon(poly, closed=True, fill=False, ec=INK, lw=1.1))
    # cockpit windows
    ax.plot([nose + 0.05 * L, nose + 0.18 * L], [gy + 2200 * S, gy + 2350 * S],
            color=INK_SOFT, lw=0.7)
    # vertical fin at tail
    fin = [P(0.88, 1850), P(0.97, 4000), P(1.04, 4000), P(1.05, 1900)]
    ax.add_patch(Polygon(fin, closed=True, fill=True, fc=PALETTE[2],
                         ec=INK, lw=1.0, alpha=0.16))
    ax.add_patch(Polygon(fin, closed=True, fill=False, ec=INK, lw=1.0))
    # horizontal stabiliser (edge) on boom
    ax.plot([nose + 0.80 * L, nose + 0.95 * L], [gy + 1950 * S, gy + 1950 * S],
            color=INK, lw=1.6, solid_capstyle="round")
    # mast + hub + rotor disc edge-on (overall length = rotor extent)
    hub_x = nose + hub_frac * L
    hub_y = gy + Hh * S
    ax.plot([hub_x, hub_x], [gy + 2500 * S, hub_y], color=INK, lw=1.4)
    ax.add_patch(Rectangle((hub_x - 1.2, hub_y - 1.0), 2.4, 1.4,
                           fill=True, fc=INK_SOFT, ec=INK, lw=0.8))
    # rotor disc edge-on: this is the ROTOR DIAMETER, not the overall length
    half_r = Rdia * S / 2
    ax.plot([hub_x - half_r, hub_x + half_r], [hub_y, hub_y], color=INK, lw=1.3)
    half = L_ovl * S / 2                      # overall extent, used for the ground line
    # tail rotor (disc face on fin), at the SPEC station so sheets 1 and 2 agree
    trc = (P(SPEC["tailrotor_from_nose"]/SPEC["L_fus"], 3500))
    ax.add_patch(Circle(trc, trdia * S / 2, fill=False, ec=INK_SOFT, lw=0.9))
    centre_mark(ax, trc[0], trc[1], trdia * S / 2 * 1.1)
    # landing gear nose + main wheels
    wr = 320 * S
    nw_x = nose + 0.12 * L
    mw_x = nw_x + wbase * S
    for wx in (nw_x, mw_x):
        ax.add_patch(Circle((wx, gy + wr), wr, fill=True, fc=PALETTE[4],
                            ec=INK, lw=0.9, alpha=0.20))
        ax.add_patch(Circle((wx, gy + wr), wr, fill=False, ec=INK, lw=0.9))
    ax.plot([nw_x, nw_x], [gy + 2 * wr, gy + 760 * S], color=INK, lw=0.8)
    ax.plot([mw_x, mw_x], [gy + 2 * wr, gy + 760 * S], color=INK, lw=0.8)
    # ground line
    ax.plot([cx - half - 4, cx + half + 4], [gy, gy], color=INK, lw=1.0)
    for k in range(-7, 8):
        gx = cx + k * (half + 4) / 7.0
        ax.plot([gx, gx - 1.4], [gy, gy - 1.4], color=INK_SOFT, lw=0.5)


# ============================================================================
# SHEET 2 -- Isometric pictorial
# ============================================================================
COS30, SIN30 = np.cos(np.deg2rad(30)), np.sin(np.deg2rad(30))


def iso(x, y, z, S, ox, oy):
    """Isometric projection (z up). Returns paper coords."""
    X = (x - y) * COS30
    Y = (x + y) * SIN30 + z
    return ox + X * S, oy + Y * S


def iso_box(ax, c, dx, dy, dz, S, ox, oy, color, alpha=0.16, lw=1.0):
    """Axis-aligned box centred at c=(x,y,z) with full sizes dx,dy,dz."""
    x0, y0, z0 = c[0] - dx / 2, c[1] - dy / 2, c[2] - dz / 2
    x1, y1, z1 = c[0] + dx / 2, c[1] + dy / 2, c[2] + dz / 2
    V = {k: iso(*p, S, ox, oy) for k, p in {
        "A": (x0, y0, z0), "B": (x1, y0, z0), "C": (x1, y1, z0), "D": (x0, y1, z0),
        "E": (x0, y0, z1), "F": (x1, y0, z1), "G": (x1, y1, z1), "H": (x0, y1, z1),
    }.items()}
    faces = [("E", "F", "G", "H"),   # top
             ("B", "C", "G", "F"),   # right (+x)
             ("A", "B", "F", "E")]   # front (-y)
    for f in faces:
        ax.add_patch(Polygon([V[k] for k in f], closed=True, fill=True,
                             fc=color, ec=INK, lw=lw, alpha=alpha))
        ax.add_patch(Polygon([V[k] for k in f], closed=True, fill=False,
                             ec=INK, lw=lw))
    return V


def iso_disc(ax, c, R, S, ox, oy, plane="xy", ec=INK, lw=1.1, fc=None, alpha=0.12):
    t = np.linspace(0, 2 * np.pi, 80)
    if plane == "xy":      # horizontal disc (rotor) z const
        pts = [iso(c[0] + R * np.cos(a), c[1] + R * np.sin(a), c[2], S, ox, oy)
               for a in t]
    else:                  # plane == "xz": vertical disc facing y (tail rotor)
        pts = [iso(c[0] + R * np.cos(a), c[1], c[2] + R * np.sin(a), S, ox, oy)
               for a in t]
    if fc is not None:
        ax.add_patch(Polygon(pts, closed=True, fill=True, fc=fc, ec=ec,
                             lw=lw, alpha=alpha))
    ax.add_patch(Polygon(pts, closed=True, fill=False, ec=ec, lw=lw))
    return pts


# ----------------------------------------------------------------------------
# Isometric solid rendering (painter's algorithm + back-face culling)
# ----------------------------------------------------------------------------
# The projection iso(x,y,z) = ((x-y)cos30, (x+y)sin30 + z) is the standard
# isometric seen from (+x, -y, +z); a face is therefore nearer the viewer the
# larger its (x - y + z). Sorting on that and dropping back-faces is what turns
# a pile of translucent outlines into a solid that reads as one object.
_CAM = np.array([1.0, -1.0, 1.0]); _CAM /= np.linalg.norm(_CAM)
_LIGHT = np.array([0.35, -0.55, 0.76]); _LIGHT /= np.linalg.norm(_LIGHT)


def _tone(t, lo=(0.17, 0.29, 0.46), hi=(0.74, 0.82, 0.90)):
    """Shade ramp. Both ends are blue: the dark end never approaches black."""
    t = float(np.clip(t, 0.0, 1.0))
    return tuple(lo[k] + (hi[k] - lo[k])*t for k in range(3))


def _paint(ax, faces, S, ox, oy, edge="face", lw=0.35, zorder=2, alpha=1.0):
    """Draw 3-D quads/tris back-to-front, culling faces that point away."""
    drawn = []
    for P in faces:
        A = np.asarray(P, float)
        n = np.cross(A[1] - A[0], A[2] - A[0])
        ln = np.linalg.norm(n)
        if ln < 1e-9:
            continue
        n /= ln
        if n @ _CAM <= 0.02:                       # back-face
            continue
        depth = float((A[:, 0] - A[:, 1] + A[:, 2]).mean())
        drawn.append((depth, A, 0.06 + 0.88*max(0.0, float(n @ _LIGHT))))
    drawn.sort(key=lambda r: r[0])                 # far first
    for _, A, sh in drawn:
        pts = [iso(x, y, z, S, ox, oy) for x, y, z in A]
        fc = _tone(sh)
        # edge="face" hides the facet seams so the surface reads as one solid
        ec = fc if edge == "face" else edge
        ax.add_patch(Polygon(pts, closed=True, fill=True, fc=fc, ec=ec,
                             lw=lw, zorder=zorder, alpha=alpha))


def _ring(zc, hw, hh, n=26, e=2.7):
    """Super-ellipse fuselage cross-section: rounded, not a box."""
    t = np.linspace(0, 2*np.pi, n, endpoint=False)
    c, si = np.cos(t), np.sin(t)
    y = hw*np.sign(c)*np.abs(c)**(2.0/e)
    z = zc + hh*np.sign(si)*np.abs(si)**(2.0/e)
    return y, z


def _loft_faces(stations, n=26, cap_first=True, cap_last=True):
    """Quads between consecutive cross-sections of (x, z_centre, half_w, half_h)."""
    rings = [(x, ) + _ring(zc, hw, hh, n) for x, zc, hw, hh in stations]
    faces = []
    for i in range(len(rings) - 1):
        x0, y0, z0 = rings[i]; x1, y1, z1 = rings[i + 1]
        for j in range(n):
            k = (j + 1) % n
            faces.append([(x0, y0[j], z0[j]), (x0, y0[k], z0[k]),
                          (x1, y1[k], z1[k]), (x1, y1[j], z1[j])])
    for cap, idx in ((cap_first, 0), (cap_last, -1)):
        if cap:
            x, y, z = rings[idx]
            faces.append([(x, y[j], z[j]) for j in range(n)][::(1 if idx else -1)])
    return faces


def _slab_faces(c, dx, dy, dz, taper=1.0):
    """Box, optionally tapered in the +x half-span (fins, stabilisers)."""
    x0, x1 = c[0] - dx/2, c[0] + dx/2
    hy0, hy1 = dy/2, dy/2*taper
    hz0, hz1 = dz/2, dz/2*taper
    V = lambda x, sy, sz, hy, hz: (x, c[1] + sy*hy, c[2] + sz*hz)
    a = [V(x0, s, t, hy0, hz0) for s, t in ((-1, -1), (1, -1), (1, 1), (-1, 1))]
    b = [V(x1, s, t, hy1, hz1) for s, t in ((-1, -1), (1, -1), (1, 1), (-1, 1))]
    f = [a[::-1], b]
    for i in range(4):
        j = (i + 1) % 4
        f.append([a[i], a[j], b[j], b[i]])
    return f


def _iso_wheel(ax, centre, r, S, ox, oy, n=22):
    """Wheel in the x-z plane (rolling about y) plus its strut."""
    t = np.linspace(0, 2*np.pi, n, endpoint=False)
    face = [(centre[0] + r*np.cos(a), centre[1], centre[2] + r*np.sin(a)) for a in t]
    _paint(ax, [face], S, ox, oy, edge=INK, lw=0.7, zorder=4)


def sheet2():
    fig, ax = new_sheet()
    S = 0.0046
    ox, oy = 118.0, 96.0

    # model coords: x fore(+)/aft, y port(+)/stbd, z up. Hub at (0,0,Hh).
    Hh, Hf = SPEC["Hh"], SPEC["Hf"]
    Wf = SPEC["Wf"]
    Rdia, trdia = SPEC["Rdia"], SPEC["trdia"]
    track, wbase = SPEC["track"], SPEC["wbase"]
    L_fus = SPEC["L_fus"]
    # model coords put the hub at the origin, so the nose is hub_from_nose ahead
    nose_x = SPEC["hub_from_nose"]
    tail_x = nose_x - L_fus
    tr_x   = nose_x - SPEC["tailrotor_from_nose"]     # tail-rotor station
    L_ovl = SPEC["L_ovl"]

    # ---- fuselage: one lofted body, nose -> boom -> tail --------------------
    # (x, z of section centre, half-width, half-height)  [mm]
    FUS = [( 4600, 1780,  110,  170), ( 4150, 1810,  470,  520),
           ( 3300, 1900,  880,  830), ( 2000, 1950, 1140,  960),
           (  600, 1960, 1180, 1000), (-1200, 1930, 1080,  930),
           (-2600, 2030,  760,  660), (-4200, 2120,  520,  470),
           (-6400, 2210,  400,  380), (-8600, 2300,  340,  330),
           (-10200, 2400,  300,  300), (-10900, 2470,  270,  285)]
    _paint(ax, _loft_faces(FUS), S, ox, oy, zorder=2)

    # ---- empennage ----------------------------------------------------------
    fin_c = (tr_x, 0.0, (2470 + Hf)/2 + 120)
    _paint(ax, _slab_faces(fin_c, 1500, 190, Hf - 2470 + 240, taper=0.78),
           S, ox, oy, zorder=3)
    _paint(ax, _slab_faces((tail_x + 1900, 0.0, 2280), 1150, 4000, 190),
           S, ox, oy, zorder=3)

    # ---- landing gear: nose wheel + two mains, each on a strut --------------
    wr = 340.0
    for wx, wy, wz_top in ((nose_x - 1400, 0.0, 1250.0),
                           (nose_x - 1400 - wbase,  track/2, 1350.0),
                           (nose_x - 1400 - wbase, -track/2, 1350.0)):
        # run the strut up INTO the belly, not to a point short of it: stopping
        # at a nominal height left the wheels looking detached from the aircraft
        ax.plot(*zip(iso(wx, wy, wr, S, ox, oy), iso(wx, wy*0.45, wz_top, S, ox, oy)),
                color=INK, lw=1.3, zorder=4, solid_capstyle="round")
        _iso_wheel(ax, (wx, wy, wr), wr, S, ox, oy)

    # ---- mast, hub, rotor ---------------------------------------------------
    ax.plot(*zip(iso(0, 0, 2900, S, ox, oy), iso(0, 0, Hh, S, ox, oy)),
            color=INK, lw=2.0, zorder=5)
    # The rotor disc is the swept ENVELOPE: draw it as an outline, not a filled
    # translucent sheet. Filling it washed the whole aircraft out to a flat blur.
    disc = [iso(Rdia/2*np.cos(a), Rdia/2*np.sin(a), Hh, S, ox, oy)
            for a in np.linspace(0, 2*np.pi, 120)]
    ax.add_patch(Polygon(disc, closed=True, fill=False, ec=INK_SOFT, lw=0.9,
                         ls=(0, (7, 3)), zorder=6))
    for ang in (18, 108, 198, 288):                       # 4 blades, tapered
        a = np.deg2rad(ang); ca, sa = np.cos(a), np.sin(a)
        r0, r1, hw0, hw1 = 760.0, Rdia/2 - 150, 300.0, 190.0
        quad = [(r0*ca - hw0*sa, r0*sa + hw0*ca, Hh),
                (r1*ca - hw1*sa, r1*sa + hw1*ca, Hh),
                (r1*ca + hw1*sa, r1*sa - hw1*ca, Hh),
                (r0*ca + hw0*sa, r0*sa - hw0*ca, Hh)]
        ax.add_patch(Polygon([iso(*q, S, ox, oy) for q in quad], closed=True,
                             fill=True, fc=PALETTE[0], ec=INK, lw=0.8,
                             alpha=0.85, zorder=7))
    hub = iso(0, 0, Hh, S, ox, oy)
    ax.add_patch(Circle(hub, 1.5, fill=True, fc=INK_SOFT, ec=INK, lw=0.9, zorder=8))

    # ---- tail rotor on the port face of the fin -----------------------------
    tr_c = (tr_x, 260.0, fin_c[2] + 250)
    trd = [iso(tr_c[0] + trdia/2*np.cos(a), tr_c[1], tr_c[2] + trdia/2*np.sin(a),
               S, ox, oy) for a in np.linspace(0, 2*np.pi, 80)]
    ax.add_patch(Polygon(trd, closed=True, fill=False, ec=INK_SOFT, lw=0.8,
                         ls=(0, (5, 2.5)), zorder=6))
    for ang in (40, 130, 220, 310):
        a = np.deg2rad(ang)
        ax.plot(*zip(iso(tr_c[0] + 180*np.cos(a), tr_c[1], tr_c[2] + 180*np.sin(a), S, ox, oy),
                     iso(tr_c[0] + (trdia/2 - 120)*np.cos(a), tr_c[1],
                         tr_c[2] + (trdia/2 - 120)*np.sin(a), S, ox, oy)),
                color=INK, lw=1.5, zorder=7, solid_capstyle="round")

    # ---- envelope dimensions ------------------------------------------------
    # Short witness ticks at each end, not full-length extension lines: the old
    # version ran them from ground level down past the sheet, leaving two lines
    # dangling in space with nothing attached to them.
    def tick(x, y, z, dz=520.0):
        ax.plot(*zip(iso(x, y, z, S, ox, oy), iso(x, y, z + dz, S, ox, oy)),
                color=INK_SOFT, lw=0.5, zorder=9)

    def iso_dim(pA, pB, text, rot, voff=0.0, hoff=0.0):
        ax.annotate("", xy=pB, xytext=pA, arrowprops=ARROW_KW)
        ax.text((pA[0] + pB[0])/2 + hoff, (pA[1] + pB[1])/2 + voff, text,
                color=INK, fontsize=FS_DIM, ha="center", va="center",
                rotation=rot, bbox=TEXT_BG, zorder=10)

    yL, zL = -3050.0, -250.0                    # starboard side, just off the ground
    for xx in (nose_x, tail_x):
        tick(xx, yL, zL)
    iso_dim(iso(nose_x, yL, zL, S, ox, oy), iso(tail_x, yL, zL, S, ox, oy),
            f"FUSELAGE {int(round(L_fus))}", rot=30, voff=-2.2)

    hx, hy = nose_x + 2900, -3050.0             # height, clear of the disc
    ax.plot(*zip(iso(0, 0, Hh, S, ox, oy), iso(hx, hy, Hh, S, ox, oy)),
            color=INK_SOFT, lw=0.45, dashes=(5, 3), zorder=9)
    h1, h2 = iso(hx, hy, 0, S, ox, oy), iso(hx, hy, Hh, S, ox, oy)
    ax.annotate("", xy=h2, xytext=h1, arrowprops=ARROW_KW)
    ax.text(h1[0] + 1.5, (h1[1] + h2[1])/2, f"{int(Hh)}  HUB HT", color=INK,
            fontsize=FS_DIM, ha="left", va="center", rotation=90,
            bbox=TEXT_BG, zorder=10)

    # ---- callouts -----------------------------------------------------------
    leader(ax, iso(-Rdia/2*0.72, Rdia/2*0.72, Hh, S, ox, oy), (92, 172),
           f"MAIN ROTOR DISC\nØ{int(Rdia)}  ({SPEC['n_blades']} BLADES)", ha="left")
    leader(ax, hub, (176, 146), "MAIN ROTOR HUB", ha="left")
    leader(ax, iso(tr_c[0], tr_c[1], tr_c[2] + trdia/2, S, ox, oy), (20, 132),
           f"TAIL ROTOR Ø{int(trdia)}\n{SPEC['n_blades']} BLADES (PORT)", ha="left")
    leader(ax, iso(fin_c[0], 0, Hf, S, ox, oy), (20, 112),
           f"VERTICAL FIN\n{int(Hf)} TIP HT", ha="left")
    leader(ax, iso(nose_x - 1400, 0, 640, S, ox, oy), (185, 80),
           f"WHEELED GEAR\nTRACK {int(track)} / W-BASE {int(wbase)}", ha="left")

    ax.text(W/2, H - M - 16, "ISOMETRIC VIEW  -  CS-MUH REFERENCE",
            color=INK, fontsize=11, ha="center", va="center", fontweight="bold")

    notes_block(ax, M + 4, 182,
                ["1.  PICTORIAL VIEW; ENVELOPE DIMS",
                 "     NOMINAL, mm.",
                 f"2.  OVERALL LENGTH INCL ROTOR {int(round(L_ovl))}.",
                 f"3.  OVERALL HEIGHT {int(Hh)} (HUB),",
                 f"     FIN TIP {int(Hf)}.",
                 "4.  ROTOR DISC SHOWN AS SWEPT",
                 "     ENVELOPE (BROKEN LINE).",
                 "5.  GENERAL TOL. \u00b125 mm."])

    finalise(ax, "ISOMETRIC GENERAL ARRANGEMENT",
             "UN-CSMUH-002", "NTS", "PICTORIAL")
    fig.savefig(HERE / "sheet2_isometric.png", dpi=150, facecolor="white")
    plt.close(fig)


# ============================================================================
# SHEET 3 -- Main rotor blade
# ============================================================================
def sheet3():
    fig, ax = new_sheet()
    R = SPEC["R"]
    chord = SPEC["blade_chord"]
    cutout = SPEC["root_cutout_frac"] * R
    station = SPEC["section_station_frac"] * R      # -> section A-A
    S = 0.0205                # paper units per mm

    x0 = 48.0                 # hub centre x (paper)
    yP = 150.0                # plan centreline (pitch axis line)

    # ---- PLAN VIEW (planform) ------------------------------------------
    pa = 0.25 * chord         # pitch axis from LE
    yLE = yP + pa * S         # leading edge above axis
    yTE = yP - (chord - pa) * S
    xc = x0 + cutout * S      # root cutout station
    xt = x0 + R * S           # tip

    # blade planform (constant chord 527, square tip)
    plan = [(xc, yLE), (xt, yLE), (xt, yTE), (xc, yTE)]
    ax.add_patch(Polygon(plan, closed=True, fill=True, fc=PALETTE[0],
                         ec=INK, lw=1.2, alpha=0.10))
    ax.add_patch(Polygon(plan, closed=True, fill=False, ec=INK, lw=1.2))
    # root cutout / grip region (hub to cutout)
    ax.add_patch(Rectangle((x0, yP - 0.5 * chord * 0.4 * S),
                           cutout * S, chord * 0.4 * S,
                           fill=True, fc=INK_SOFT, ec=INK, lw=0.9, alpha=0.25))
    ax.plot([x0, x0], [yP - 0.5 * chord * 0.4 * S, yP + 0.5 * chord * 0.4 * S],
            color=INK, lw=1.2)
    # hub centre mark
    centre_mark(ax, x0, yP, 5)
    ax.text(x0 - 4, yP, "ROTOR\nAXIS", color=INK, fontsize=FS_NOTE,
            ha="right", va="center")
    # pitch axis line (0.25c) along span
    ax.plot([x0, xt + 4], [yP, yP], color=INK_SOFT, lw=0.7,
            dashes=(8, 2, 1.5, 2))
    ax.text(xt + 5, yP, "PITCH AXIS 0.25c", color=INK, fontsize=FS_NOTE,
            ha="left", va="center")

    # SECTION line A-A at 0.75R
    xs = x0 + station * S
    ax.plot([xs, xs], [yLE + 5, yTE - 5], color=PALETTE[1], lw=1.2,
            dashes=(9, 2, 1.5, 2))
    # BOTH section arrows must point the SAME way (the direction of sight for
    # section A-A). They were mirrored, which reads as two different sections.
    for yy in (yLE + 5, yTE - 5):
        ax.annotate("", xy=(xs - 3.2, yy), xytext=(xs, yy),
                    arrowprops=dict(arrowstyle="-|>", color=PALETTE[1], lw=1.2,
                                    mutation_scale=10))
    ax.text(xs, yLE + 6.5, "A", color=PALETTE[1], fontsize=10, ha="center",
            va="bottom", fontweight="bold")
    ax.text(xs, yTE - 6.5, "A", color=PALETTE[1], fontsize=10, ha="center",
            va="top", fontweight="bold")

    view_title(ax, x0 + R * S / 2, yP + 24, "BLADE PLANFORM")

    # plan dims
    dim_h(ax, x0, xt, yLE + 14, yLE, f"R = {int(round(R))}  (RADIUS)")
    dim_h(ax, x0, xc, yLE + 8, yLE,
          f"{int(round(cutout))}  ROOT CUTOUT ({SPEC['root_cutout_frac']:.2f}R)")
    dim_h(ax, x0, xs, yTE - 10, yTE,
          f"{int(round(station))}  ({SPEC['section_station_frac']:.2f}R) SEC A-A")
    dim_v(ax, yTE, yLE, xc - 4, xc, f"c = {int(round(chord))}")

    # ---- EDGE / SIDE VIEW (thickness) ----------------------------------
    yE = 92.0
    thick = 0.12 * chord
    tE = thick * S
    # edge silhouette: constant-thickness lens (NACA 0012, exaggerated)
    edge = [(xc, yE + tE / 2), (xt, yE + tE / 2),
            (xt, yE - tE / 2), (xc, yE - tE / 2)]
    ax.add_patch(Polygon(edge, closed=True, fill=True,
                         fc=PALETTE[0], ec=INK, lw=1.1, alpha=0.12))
    ax.add_patch(Polygon(edge, closed=True, fill=False, ec=INK, lw=1.1))
    ax.plot([x0, xc], [yE, yE], color=INK, lw=1.0)     # grip stub
    ax.plot([x0, xt + 4], [yE, yE], color=INK_SOFT, lw=0.6,
            dashes=(8, 2, 1.5, 2))
    view_title(ax, x0 + R * S / 2, yE - 12, "EDGE VIEW (THICKNESS EXAGGERATED)")
    dim_v(ax, yE - tE / 2, yE + tE / 2, xc - 4, xc, f"t={thick:.1f} (12%c)")
    leader(ax, (xt, yE), (xt + 6, yE + 9),
           f"LINEAR TWIST {SPEC['blade_twist_deg']:.0f}°\n(WASHOUT ROOT→TIP)", ha="left")

    notes_block(ax, M + 4, 60,
                ["1.  AIRFOIL: NACA 0012 CONSTANT.",
                 "2.  %d BLADES; Ø%d ROTOR." % (SPEC["n_blades"], int(SPEC["Rdia"])),
                 "3.  LINEAR TWIST %.0f° (WASHOUT)." % SPEC["blade_twist_deg"],
                 "4.  PITCH AXIS AT 0.25 CHORD.",
                 "5.  TOL: SPAN ±5, CHORD ±1.5,",
                 "     TWIST ±0.25°.",
                 "6.  SECTION A-A SHOWN ON SHT 4,",
                 "     VIEWED INBOARD (ARROWS)."])

    finalise(ax, "MAIN ROTOR BLADE", "UN-CSMUH-003", "1:50 (approx)",
             "AL / CFRP SPAR")
    fig.savefig(HERE / "sheet3_main_rotor_blade.png", dpi=150, facecolor="white")
    plt.close(fig)


# ============================================================================
# SHEET 4 -- Section A-A : NACA 0012 airfoil
# ============================================================================
def naca0012(chord, n=160):
    t = 0.12
    beta = np.linspace(0.0, np.pi, n)
    x = (1 - np.cos(beta)) / 2.0
    yt = 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x**2
                  + 0.2843 * x**3 - 0.1015 * x**4)   # open-TE coeff, matching
                                                     # 01_geometry/generate_geometry.py
    xu, yu = x, yt
    xl, yl = x[::-1], -yt[::-1]
    X = np.concatenate([xu, xl]) * chord
    Y = np.concatenate([yu, yl]) * chord
    return X, Y, x * chord, yt * chord


def sheet4():
    fig, ax = new_sheet()
    matplotlib.rcParams['hatch.color'] = INK_SOFT
    matplotlib.rcParams['hatch.linewidth'] = 0.4
    chord = SPEC["blade_chord"]
    S = 175.0 / chord          # fit chord ~175 paper units
    X, Y, xc_arr, yt_arr = naca0012(chord)

    ox = 48.0                  # LE x (paper)
    oy = 120.0                 # chord-line y (paper)
    px = ox + X * S
    py = oy + Y * S

    # hatched section fill (light, non-black) -- hatch colour from rcParams
    ax.add_patch(Polygon(np.column_stack([px, py]), closed=True, fill=True,
                         fc=PALETTE[0], ec=INK_SOFT, lw=0.0, alpha=0.10,
                         hatch="////"))
    # crisp outline on top
    ax.add_patch(Polygon(np.column_stack([px, py]), closed=True, fill=False,
                         ec=INK, lw=1.4))

    # chord line
    xTE = ox + chord * S
    ax.plot([ox, xTE], [oy, oy], color=INK_SOFT, lw=0.7,
            dashes=(8, 2, 1.5, 2))

    # max thickness location 0.30c
    xtmax = ox + 0.30 * chord * S
    tmax = 0.12 * chord
    yt_at = np.interp(0.30 * chord, xc_arr, yt_arr)
    ax.plot([xtmax, xtmax], [oy - yt_at * S, oy + yt_at * S],
            color=INK_SOFT, lw=0.6, dashes=(4, 2))

    # pitch axis 0.25c
    xpa = ox + 0.25 * chord * S
    ya = np.interp(0.25 * chord, xc_arr, yt_arr)
    ax.plot([xpa, xpa], [oy - ya * S - 6, oy + ya * S + 6],
            color=PALETTE[1], lw=1.0, dashes=(9, 2, 1.5, 2))
    ax.add_patch(Circle((xpa, oy), 1.4, fill=False, ec=PALETTE[1], lw=1.0))
    ax.plot([xpa - 2.1, xpa + 2.1], [oy, oy], color=PALETTE[1], lw=0.8)
    ax.plot([xpa, xpa], [oy - 2.1, oy + 2.1], color=PALETTE[1], lw=0.8)
    ax.text(xpa, oy + ya * S + 8, "PITCH AXIS 0.25c", color=PALETTE[1],
            fontsize=FS_NOTE, ha="center", va="bottom")

    # LE radius detail
    rle = 1.1019 * (0.12 ** 2) * chord     # ~ 8.36 mm
    leader(ax, (ox + rle * S, oy),
           (ox - 16, oy - 22),
           f"LE RADIUS\nr = {rle:.1f}", ha="left")
    ax.add_patch(Circle((ox + rle * S, oy), rle * S, fill=False,
                        ec=PALETTE[2], lw=0.8))

    view_title(ax, ox + chord * S / 2, oy + 40,
               "SECTION  A-A   NACA 0012   (r/R = %.2f, VIEWED INBOARD)"
               % SPEC["section_station_frac"])

    # ---- dimensions -----------------------------------------------------
    yt_max = tmax / 2
    # chord (below)
    dim_h(ax, ox, xTE, oy - 30, oy, f"CHORD c = {int(round(chord))}")
    # 0.30c station (below, nearer)
    dim_h(ax, ox, xtmax, oy - 22, oy, "0.30c = %.1f" % (0.30*chord))
    # 0.25c station
    dim_h(ax, ox, xpa, oy - 14, oy, "0.25c = %.1f" % (0.25*chord))
    # max thickness (vertical, right of tmax line)
    dim_v(ax, oy - yt_max * S, oy + yt_max * S, xtmax + 14, xtmax,
          f"t = {tmax:.1f}  (12%c)", side="right")

    notes_block(ax, M + 4, 70,
                ["1.  AIRFOIL NACA 0012 (SYMMETRIC).",
                 "2.  MAX THICKNESS 12%%c = %.1f @ 0.30c." % tmax,
                 "3.  LE RADIUS 1.1019 t²c = %.1f." % rle,
                 "4.  PITCH / FEATHER AXIS AT 0.25c.",
                 "5.  COORDS GENERATED ANALYTICALLY.",
                 "6.  SECTION TOL ±0.5; CONTOUR ±0.25.",
                 "7.  HATCHING = CUT SOLID (CONVENTION)."])

    finalise(ax, "SECTION A-A  -  NACA 0012 AIRFOIL",
             "UN-CSMUH-004", "1:3 (approx)", "NACA 0012")
    fig.savefig(HERE / "sheet4_section_AA_airfoil.png", dpi=150,
                facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    sheet1()
    sheet2()
    sheet3()
    sheet4()
    print("Done. PNGs written to:", HERE)
