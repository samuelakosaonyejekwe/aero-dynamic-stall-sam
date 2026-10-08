# -*- coding: utf-8 -*-
"""
08_engineering_drawings / draw_engineering.py
---------------------------------------------
Four illustrative drawing sheets for the two cases of the study, and the
table of the dimensions they show.

Standing: illustrations, not for manufacture. Case B is not a rotor
calculation; sheets 1 and 2 show the rotor whose blade station gives its
condition, and nothing in the study depends on them.

Every dimension printed on a sheet is read from DIMENSIONS below, where each
entry carries its source; a sheet prints no other dimension. Sheets 1 and 2
are schematic (a generic four-blade rotor and a generic body outline, both
sheets drawing the body from the one set of coordinates BODY), with only the
rotor radius, the blade chord and the analysis station dimensioned. Sheets 3
and 4 draw the blade planform and the section outline to scale; every sheet
states its scale, or that it has none, in its title block. The section is
the closed outline of unistall.naca4.outline, the one of 01_geometry.

Section A-A of sheet 3 is seen looking the way its cutting-plane arrows
point. The planform is seen from above with the leading edge uppermost, so
looking outboard puts the leading edge on the left; `section_side` works
that out, and the section is drawn on the side it gives.

Outputs (all in this folder)
  sheet1_general_arrangement_3view.png
  sheet2_isometric.png
  sheet3_main_rotor_blade.png
  sheet4_section_AA_airfoil.png
  dimensions.csv            name, value, unit, source
  sheet_record.csv          what each sheet states: its scale, the viewing
                            direction of section A-A, the body coordinates
  figure_record.csv         resolution and smallest lettering of each sheet

Run from the repository root:
  PYTHONPATH=. python3 08_engineering_drawings/draw_engineering.py

Author: Akosa Samuel Onyejekwe (independent)
"""
import csv
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch, Polygon, Rectangle
from scipy.spatial import ConvexHull

import project_meta as pm
from unistall import naca4
from unistall.style import INK, INK_SOFT, PALETTE, apply_style

HERE = Path(__file__).resolve().parent
GEOMETRY = HERE.parent/"01_geometry"/"section_geometry_summary.csv"

SRC_TUNNEL = "NASA TM-84245"
SRC_ROTOR = pm.ROTOR_SOURCE
SRC_ROOT = "%s (Table 3: root cutout %.2f in of %.2f in)" % (SRC_ROTOR.split(" (")[0], pm.ROOT_CUTOUT_IN,
                                                             pm.TIP_RADIUS_IN)
CASES = pm.read_setup()[0]
CHORD_A, CHORD_B = CASES[pm.CASE_A]["chord"], CASES[pm.CASE_B]["chord"]
THICKNESS_RATIO = naca4.thickness_ratio(pm.SECTION)
X_TMAX = float(pd.read_csv(GEOMETRY).set_index("property").loc["x_at_max_thickness_x_c", "value"])
PITCH_AXIS = 0.25                                  # pitch axis of the tunnel model, fraction of chord

# The only dimensions the sheets may print: key -> (name, value, unit, source).
DIMENSIONS = {
    "A_section": ("Case A section", "NACA " + pm.SECTION, "-", SRC_TUNNEL),
    "A_chord": ("Case A chord of the tunnel model", CHORD_A, "m", SRC_TUNNEL),
    "A_thickness_ratio": ("Case A maximum thickness, percent of chord", round(100.0*THICKNESS_RATIO), "% of chord",
                          "NACA %s designation" % pm.SECTION),
    "A_thickness": ("Case A maximum thickness", round(THICKNESS_RATIO*CHORD_A, 4), "m",
                    "derived: %d %% of the Case A chord" % round(100.0*THICKNESS_RATIO)),
    "A_tmax_station": ("Case A station of maximum thickness, fraction of chord", X_TMAX, "chord",
                       "four-digit thickness form (01_geometry)"),
    "A_pitch_axis": ("Case A pitch axis from the leading edge, fraction of chord", PITCH_AXIS, "chord", SRC_TUNNEL),
    "B_radius": ("Case B rotor radius", round(pm.ROTOR["radius_ft"]*pm.FT_TO_M, 2), "m", SRC_ROTOR),
    "B_radius_ft": ("Case B rotor radius, in feet", pm.ROTOR["radius_ft"], "ft", SRC_ROTOR),
    "B_chord": ("Case B blade chord", CHORD_B, "m", SRC_ROTOR),
    "B_chord_in": ("Case B blade chord, in inches", pm.ROTOR["chord_in"], "in", SRC_ROTOR),
    "B_blades": ("Case B number of blades", pm.ROTOR["blades"], "-", SRC_ROTOR),
    "B_speed": ("Case B rotor speed", pm.ROTOR["speed_rpm"], "rpm", SRC_ROTOR),
    "B_station": ("Case B analysis station, fraction of rotor radius", pm.STATION_OVER_R, "R",
                  "case definition of this study"),
    "B_root_end": ("Root end of the blade as drawn, fraction of rotor radius", pm.ROOT_END_OVER_R, "R", SRC_ROOT),
}

# The generic body of sheets 1 and 2, in rotor radii, hub at the origin, x aft, y to starboard, z up. It
# carries no dimension and stands for no aircraft.
BODY = dict(mast=0.15, cabin_x=-0.06, cabin_z=-0.30, cabin_a=0.42, cabin_b=0.13, cabin_c=0.15,
            boom_x=(0.28, 1.12), boom_z=(-0.27, -0.25), boom_half_width=(0.045, 0.015),
            boom_half_height=(0.06, 0.02),
            fin=((1.04, -0.26), (1.10, -0.06), (1.15, -0.06), (1.12, -0.26)),          # (x, z) corners
            tailplane_x=(1.00, 1.07), tailplane_half_span=0.11,
            skid_z=-0.52, skid_x=(-0.36, 0.24), strut_x=(-0.26, 0.14), strut_y=0.11)

SHEET_W, SHEET_H, MARGIN = 420.0, 297.0, 10.0      # millimetres on the sheet (A3)
SHEET_FORMAT = "A3"
N_SHEETS = 4
FS_SMALL, FS_TEXT, FS_VIEW, FS_BANNER = 11.0, 11.5, 13.5, 13.5
WHITE_BOX = dict(boxstyle="square,pad=0.25", fc="white", ec="none")
FILL = "#dbe6f2"
NOTE_ROTOR = ("Rotor radius, blade chord, number of blades and rotor speed are those of the UH-60A main rotor "
              "(%s). The real blade uses SC1095 and SC1094 R8 sections and a swept tip; the NACA %s section and "
              "the untwisted, unswept planform drawn are an illustration." % (SRC_ROTOR.split(" (")[0], pm.SECTION))
NOTE_ROOT = ("The blade is drawn from the root cutout of the source, %s R, on every sheet; inboard of it the root "
             "end is schematic." % pm.ROOT_END_OVER_R)
BANNER_SCHEMATIC = "SCHEMATIC; ILLUSTRATION ONLY"
SHEET_RECORD = []                                  # what each sheet states, written to sheet_record.csv
FIGURE_RECORD = []


def val(key):
    """The numerical value of a dimension."""
    return DIMENSIONS[key][1]


def label(key):
    """A dimension as printed on a sheet: value and unit, kept together on one line."""
    _, value, unit, _ = DIMENSIONS[key]
    return "%s" % value if unit == "-" else ("%s %s" % (value, unit)).replace(" ", "\u00a0")


def label_radius():
    return "R = %s (%s)" % (label("B_radius"), label("B_radius_ft"))


def label_chord_b():
    return "%s (%s)" % (label("B_chord"), label("B_chord_in"))


def label_station():
    return "r\u00a0=\u00a0%s" % label("B_station")


def write_dimension_table():
    """dimensions.csv from the dictionary."""
    with open(HERE/"dimensions.csv", "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["name", "value", "unit", "source"])
        for name, value, unit, source in DIMENSIONS.values():
            writer.writerow([name, value, unit, source])


# ---------------------------------------------------------------- sheet furniture
def line(ax, xs, ys, color=INK, lw=0.9, ls="-", zorder=3):
    ax.plot(xs, ys, color=color, lw=lw, ls=ls, zorder=zorder, solid_capstyle="round")


def title_block(ax, sheet_title, number, scale):
    """Title, study, author, date, sheet number, scale and status, bottom right."""
    x0, y0, w, row = 252.0, MARGIN, SHEET_W - MARGIN - 252.0, 10.0
    ax.add_patch(Rectangle((x0, y0), w, 5*row + 2.0, fc="white", ec=INK, lw=1.4, zorder=5))
    for n in range(1, 5):
        line(ax, [x0, x0 + w], [y0 + n*row, y0 + n*row], lw=0.8, zorder=6)
    line(ax, [x0 + 88.0, x0 + 88.0], [y0 + 2*row, y0 + 3*row], lw=0.8, zorder=6)
    line(ax, [x0 + 124.0, x0 + 124.0], [y0 + 2*row, y0 + 3*row], lw=0.8, zorder=6)
    kw = dict(color=INK, zorder=7, va="center")
    ax.text(x0 + w/2, y0 + 4.6*row, sheet_title, ha="center", fontsize=FS_VIEW, fontweight="bold", **kw)
    ax.text(x0 + w/2, y0 + 3.5*row, pm.TITLE, ha="center", fontsize=FS_SMALL, **kw)
    ax.text(x0 + 2.0, y0 + 2.5*row, "DRAWN BY  " + pm.AUTHOR_BAND, ha="left", fontsize=FS_SMALL - 0.4, **kw)
    ax.text(x0 + 106.0, y0 + 2.5*row, pm.STUDY_DATE_ISO, ha="center", fontsize=FS_SMALL, **kw)
    ax.text(x0 + 141.0, y0 + 2.5*row, "SHEET %d OF %d" % (number, N_SHEETS), ha="center", fontsize=FS_SMALL,
            fontweight="bold", **kw)
    ax.text(x0 + w/2, y0 + 1.5*row, "SCALE: " + scale, ha="center", fontsize=FS_SMALL, **kw)
    ax.text(x0 + w/2, y0 + 0.5*row, "ILLUSTRATIVE; NOT FOR MANUFACTURE", ha="center", fontsize=FS_TEXT,
            fontweight="bold", **kw)


def dimension_table(ax, keys):
    """The dictionary entries shown on this sheet, with their sources, bottom left."""
    x0, row, widths = MARGIN + 4.0, 6.3, (132.0, 30.0, 70.0)
    top = MARGIN + 4.0 + row*(len(keys) + 1)
    ax.text(x0, top + 4.0, "DIMENSIONS ON THIS SHEET AND THEIR SOURCES", fontsize=FS_SMALL, fontweight="bold",
            color=INK, va="center")
    xs = np.concatenate([[x0], x0 + np.cumsum(widths)])
    rows = [("dimension", "value", "source")] + [(DIMENSIONS[k][0], label(k), DIMENSIONS[k][3].split(" (")[0])
                                                 for k in keys]
    for n, cells in enumerate(rows):
        y = top - row*(n + 0.5)
        for x, cell in zip(xs, cells, strict=False):
            ax.text(x + 1.5, y, cell, fontsize=FS_SMALL, color=INK, va="center",
                    fontweight="bold" if n == 0 else "normal")
        line(ax, [xs[0], xs[-1]], [top - row*(n + 1)]*2, color=INK_SOFT, lw=0.5)
    line(ax, [xs[0], xs[-1]], [top, top], color=INK_SOFT, lw=0.5)
    for x in xs:
        line(ax, [x, x], [top - row*len(rows), top], color=INK_SOFT, lw=0.5)


def new_sheet(number, sheet_title, banner, keys, scale):
    """A blank sheet with its border, banner, title block and dimension table.
    `scale` is what the title block states."""
    fig = plt.figure(figsize=(SHEET_W/25.4, SHEET_H/25.4))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, SHEET_W)
    ax.set_ylim(0, SHEET_H)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.add_patch(Rectangle((MARGIN, MARGIN), SHEET_W - 2*MARGIN, SHEET_H - 2*MARGIN, fc="none", ec=INK, lw=1.6))
    ax.text(SHEET_W/2, SHEET_H - MARGIN - 8.0, banner, ha="center", va="center", fontsize=FS_BANNER,
            fontweight="bold", color=PALETTE[1])
    title_block(ax, sheet_title, number, scale)
    dimension_table(ax, keys)
    SHEET_RECORD.append(dict(sheet=number, item="scale stated in the title block", value=scale))
    return fig, ax


def save(fig, name):
    """Write a sheet and enter it in the figure record."""
    FIGURE_RECORD.append((name, 150.0, pm.smallest_lettering(fig)))
    with plt.rc_context({"savefig.bbox": "standard"}):
        fig.savefig(HERE/name, dpi=150, facecolor="white")
    plt.close(fig)


def note(ax, x, y, text, width=70):
    """A wrapped note with its top-left corner at (x, y)."""
    ax.text(x, y, textwrap.fill(text, width), fontsize=FS_SMALL, color=INK, va="top", ha="left", linespacing=1.35)


def view_title(ax, x, y, text):
    ax.text(x, y, text, fontsize=FS_VIEW, fontweight="bold", color=INK, ha="center", va="center")


def arrow(ax, p, q, style="<|-|>"):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle=style, mutation_scale=9, color=INK, lw=0.9, shrinkA=0,
                                 shrinkB=0, zorder=6))


def dim_h(ax, x1, x2, y, text, from_y, above=True):
    """A horizontal dimension at height y, with extension lines from from_y = (y at x1, y at x2)."""
    over = 2.0 if y > from_y[0] else -2.0
    for x, y_from in zip((x1, x2), from_y, strict=True):
        line(ax, [x, x], [y_from, y + over], lw=0.6)
    arrow(ax, (x1, y), (x2, y))
    ax.text(0.5*(x1 + x2), y + (2.0 if above else -2.0), text, fontsize=FS_TEXT, color=INK, ha="center",
            va="bottom" if above else "top")


def dim_v(ax, y1, y2, x, text, from_x):
    """A vertical dimension at position x, with extension lines from from_x; the text sits beside it."""
    over = 2.0 if x > from_x else -2.0
    for y in (y1, y2):
        line(ax, [from_x, x + over], [y, y], lw=0.6)
    arrow(ax, (x, y1), (x, y2))
    ax.text(x + (2.5 if over > 0 else -2.5), 0.5*(y1 + y2), text, fontsize=FS_TEXT, color=INK,
            ha="left" if over > 0 else "right", va="center")


def leader(ax, tip, anchor, text, ha="left"):
    """Text at anchor with an arrow to tip."""
    ax.annotate(text, xy=tip, xytext=anchor, fontsize=FS_TEXT, color=INK, ha=ha, va="center", zorder=7,
                arrowprops=dict(arrowstyle="-|>", color=INK, lw=0.8, shrinkA=3, shrinkB=0))


def section_outline(n=161):
    """The closed section outline for unit chord, leading edge at the origin
    (unistall.naca4.outline, the outline of 01_geometry)."""
    X, Y, _, _ = naca4.outline(pm.SECTION, n)
    return X, Y


def blade_corners(radial, normal, r_root, r_tip, chord):
    """The corners of a blade drawn about its pitch axis: PITCH_AXIS of the
    chord ahead of the axis, on the side `normal` points to, the rest behind."""
    lead, trail = PITCH_AXIS*chord, -(1.0 - PITCH_AXIS)*chord
    return [radial*r + normal*s for r, s in ((r_root, trail), (r_tip, trail), (r_tip, lead), (r_root, lead))]


def section_side(looking_outboard):
    """The side of the sheet the leading edge of section A-A belongs on, for
    a planform seen from above with the span to the right and the leading
    edge uppermost. To an observer looking along d with the blade's upper
    surface upwards, the right-hand side is d x up; with d outboard (+x on
    the sheet) that is the trailing edge, so the leading edge is on the
    left, and looking inboard it is on the right."""
    d = np.array([1.0 if looking_outboard else -1.0, 0.0, 0.0])
    right = np.cross(d, np.array([0.0, 0.0, 1.0]))            # sheet y is towards the leading edge
    return "right" if right[1] > 0.0 else "left"


# ---------------------------------------------------------------- sheet 1
def rotor_plan(ax, cx, cy, r_draw):
    """Plan of the rotor: disc, station circle, hub and the blades at true
    chord about their pitch axes, leading edges towards the sense of rotation
    drawn (counter-clockwise seen from above)."""
    chord = r_draw*val("B_chord")/val("B_radius")
    ax.add_patch(Circle((cx, cy), r_draw, fc="none", ec=INK_SOFT, lw=0.8, ls="-.", zorder=3))
    ax.add_patch(Circle((cx, cy), val("B_station")*r_draw, fc="none", ec=PALETTE[1], lw=1.0, ls="--", zorder=3))
    root = val("B_root_end")*r_draw
    for k in range(val("B_blades")):
        psi = np.radians(45.0) + 2.0*np.pi*k/val("B_blades")
        radial, normal = np.array([np.cos(psi), np.sin(psi)]), np.array([-np.sin(psi), np.cos(psi)])
        corners = blade_corners(radial, normal, root, r_draw, chord)
        ax.add_patch(Polygon(np.array(corners) + [cx, cy], closed=True, fc=PALETTE[0], ec=INK, lw=0.8, zorder=4))
        line(ax, [cx, cx + root*radial[0]], [cy, cy + root*radial[1]], lw=1.2, zorder=4)
    ax.add_patch(Circle((cx, cy), 2.6, fc="white", ec=INK, lw=1.0, zorder=5))
    t = np.radians(np.linspace(196.0, 226.0, 30))
    line(ax, cx + 1.07*r_draw*np.cos(t), cy + 1.07*r_draw*np.sin(t), lw=1.0)
    arrow(ax, (cx + 1.07*r_draw*np.cos(t[-2]), cy + 1.07*r_draw*np.sin(t[-2])),
          (cx + 1.07*r_draw*np.cos(t[-1]), cy + 1.07*r_draw*np.sin(t[-1])), style="-|>")


def body_plan(ax, cx, cy, r):
    """The generic body in plan, nose to the left, from BODY."""
    t = np.linspace(0.0, 2.0*np.pi, 80)
    ax.add_patch(Polygon(np.column_stack([cx + r*(BODY["cabin_x"] + BODY["cabin_a"]*np.cos(t)),
                                          cy + r*BODY["cabin_b"]*np.sin(t)]), closed=True,
                         fc=FILL, ec=INK, lw=0.9, zorder=2))
    (x0, x1), (w0, w1) = BODY["boom_x"], BODY["boom_half_width"]
    boom = [(cx + x0*r, cy + w0*r), (cx + x1*r, cy + w1*r), (cx + x1*r, cy - w1*r), (cx + x0*r, cy - w0*r)]
    ax.add_patch(Polygon(boom, closed=True, fc=FILL, ec=INK, lw=0.9, zorder=1))
    t0, t1 = BODY["tailplane_x"]
    half = BODY["tailplane_half_span"]
    ax.add_patch(Rectangle((cx + t0*r, cy - half*r), (t1 - t0)*r, 2.0*half*r, fc=FILL, ec=INK, lw=0.9, zorder=1))


def plan_view(ax, cx, cy, r):
    body_plan(ax, cx, cy, r)
    rotor_plan(ax, cx, cy, r)
    dim_h(ax, cx, cx + r, cy - r - 11.0, "rotor radius " + label_radius(), (cy - 4.0, cy - 14.0), above=False)
    leader(ax, (cx + 0.55*r*np.cos(np.pi/4), cy + 0.55*r*np.sin(np.pi/4) + 1.5), (cx + 30.0, cy + r + 11.0),
           "blade chord " + label_chord_b())
    ang = np.radians(118.0)
    leader(ax, (cx + val("B_station")*r*np.cos(ang), cy + val("B_station")*r*np.sin(ang)),
           (cx - 44.0, cy + r + 11.0), "analysis station " + label_station(), ha="right")
    ax.text(cx - 1.10*r, cy - 0.62*r, "sense of rotation\nas drawn", fontsize=FS_SMALL, color=INK, ha="right",
            va="center")
    view_title(ax, cx, cy + r + 24.0, "PLAN VIEW (FROM ABOVE)")


def body_elevation(ax, cx, cy, r):
    """The generic body from the side, nose to the left, from BODY; returns
    the heights it drew, in rotor radii."""
    line(ax, [cx, cx], [cy - BODY["mast"]*r, cy], lw=1.6)
    t = np.linspace(0.0, 2.0*np.pi, 80)
    ax.add_patch(Polygon(np.column_stack([cx + r*(BODY["cabin_x"] + BODY["cabin_a"]*np.cos(t)),
                                          cy + r*(BODY["cabin_z"] + BODY["cabin_c"]*np.sin(t))]),
                         closed=True, fc=FILL, ec=INK, lw=0.9, zorder=2))
    (x0, x1), (z0, z1), (h0, h1) = BODY["boom_x"], BODY["boom_z"], BODY["boom_half_height"]
    boom = [(cx + x0*r, cy + (z0 + h0)*r), (cx + x1*r, cy + (z1 + h1)*r), (cx + x1*r, cy + (z1 - h1)*r),
            (cx + x0*r, cy + (z0 - h0)*r)]
    ax.add_patch(Polygon(boom, closed=True, fc=FILL, ec=INK, lw=0.9, zorder=1))
    ax.add_patch(Polygon([(cx + x*r, cy + z*r) for x, z in BODY["fin"]], closed=True, fc=FILL, ec=INK, lw=0.9,
                         zorder=1))
    for x in BODY["strut_x"]:
        line(ax, [cx + x*r, cx + x*r], [cy + (BODY["cabin_z"] - 0.9*BODY["cabin_c"])*r, cy + BODY["skid_z"]*r])
    line(ax, [cx + BODY["skid_x"][0]*r, cx + BODY["skid_x"][1]*r], [cy + BODY["skid_z"]*r]*2, lw=1.6)
    return dict(mast_over_R=BODY["mast"], fin_top_over_R=max(z for _, z in BODY["fin"]))


def front_view(ax, cx, cy, r):
    """Rotor edge-on over the generic body, seen from the nose."""
    line(ax, [cx - r, cx + r], [cy, cy], color=PALETTE[0], lw=2.2, zorder=4)
    line(ax, [cx, cx], [cy - BODY["mast"]*r, cy], lw=1.6)
    ax.add_patch(Ellipse((cx, cy + BODY["cabin_z"]*r), 2.0*BODY["cabin_b"]*r, 2.0*BODY["cabin_c"]*r, fc=FILL,
                         ec=INK, lw=0.9, zorder=2))
    for s in (-1.0, 1.0):
        y = BODY["strut_y"]*r*s
        line(ax, [cx + 0.6*y, cx + y], [cy + (BODY["cabin_z"] - 0.85*BODY["cabin_c"])*r, cy + BODY["skid_z"]*r])
        line(ax, [cx + y - 2.0, cx + y + 2.0], [cy + BODY["skid_z"]*r]*2, lw=1.6)
    dim_h(ax, cx, cx + r, cy + 9.0, label_radius(), (cy + 2.0, cy + 2.0))
    view_title(ax, cx, cy + BODY["skid_z"]*r - 9.0, "FRONT VIEW")


def side_view(ax, cx, cy, r):
    """Rotor edge-on over the generic body, nose to the left; returns the
    body heights drawn."""
    line(ax, [cx - r, cx + r], [cy, cy], color=PALETTE[0], lw=2.2, zorder=4)
    drawn = body_elevation(ax, cx, cy, r)
    dim_h(ax, cx, cx + r, cy + 9.0, label_radius(), (cy + 2.0, cy + 2.0))
    view_title(ax, cx, cy + BODY["skid_z"]*r - 9.0, "SIDE VIEW")
    return drawn


def record_body(number, drawn):
    """Enter the body heights a sheet drew in the sheet record."""
    for key, value in drawn.items():
        SHEET_RECORD.append(dict(sheet=number, item="generic body: " + key, value=value))


def sheet1():
    keys = ["B_radius", "B_radius_ft", "B_chord", "B_chord_in", "B_blades", "B_speed", "B_station", "B_root_end"]
    fig, ax = new_sheet(1, "GENERAL ARRANGEMENT, THREE VIEWS", BANNER_SCHEMATIC, keys,
                        "NOT TO SCALE; VIEWS NOT ALIGNED IN PROJECTION")
    r = 62.0
    plan_view(ax, 124.0, 176.0, r)
    front_view(ax, 318.0, 246.0, r)
    record_body(1, side_view(ax, 312.0, 172.0, r))
    note(ax, 256.0, 118.0,
         "Case B illustration. Rotor: %s blades, %s. Dashed circle: analysis station at %s. "
         % (label("B_blades"), label("B_speed"), label_station()) + NOTE_ROTOR
         + " The body outline is generic and carries no dimension; the sense of rotation is the one assumed for "
         "the sketch. " + NOTE_ROOT)
    save(fig, "sheet1_general_arrangement_3view.png")


# ---------------------------------------------------------------- sheet 2
def iso(points, origin, scale):
    """Isometric projection of points (x, y, z) to sheet coordinates."""
    p = np.atleast_2d(np.asarray(points, float))
    u = (p[:, 0] - p[:, 1])*np.cos(np.radians(30.0))
    v = (p[:, 0] + p[:, 1])*0.5 + p[:, 2]
    return np.column_stack([origin[0] + scale*u, origin[1] + scale*v])


def hull_patch(ax, points, zorder):
    """Fill the silhouette of a convex body given the projections of points on it."""
    ax.add_patch(Polygon(points[ConvexHull(points).vertices], closed=True, fc=FILL, ec=INK, lw=0.9, zorder=zorder))


def ring(centre, ry, rz, n=40):
    """A cross-section ring in a plane of constant x."""
    t = np.linspace(0.0, 2.0*np.pi, n)
    return np.column_stack([np.full(n, centre[0]), centre[1] + ry*np.cos(t), centre[2] + rz*np.sin(t)])


def body_iso(ax, origin, scale):
    """The generic body below the rotor, from BODY: fin, tail boom and cabin,
    drawn far to near; returns the body heights drawn, in rotor radii."""
    fin = np.array([[x, 0.0, z] for x, z in BODY["fin"]])
    ax.add_patch(Polygon(iso(fin, origin, scale), closed=True, fc=FILL, ec=INK, lw=0.9, zorder=1))
    (x0, x1), (z0, z1) = BODY["boom_x"], BODY["boom_z"]
    (w0, w1), (h0, h1) = BODY["boom_half_width"], BODY["boom_half_height"]
    boom = np.vstack([ring((x0, 0, z0), w0, h0), ring((x1, 0, z1), w1, h1)])
    hull_patch(ax, iso(boom, origin, scale), 2)
    cabin = np.vstack([ring((BODY["cabin_x"] + BODY["cabin_a"]*np.cos(a), 0, BODY["cabin_z"]),
                            BODY["cabin_b"]*np.sin(a), BODY["cabin_c"]*np.sin(a))
                       for a in np.linspace(0.02, np.pi - 0.02, 30)])
    hull_patch(ax, iso(cabin, origin, scale), 3)
    line(ax, *iso([[0, 0, -BODY["mast"]], [0, 0, 0]], origin, scale).T, lw=1.6, zorder=4)
    return dict(mast_over_R=BODY["mast"], fin_top_over_R=float(fin[:, 2].max()))


def rotor_iso(ax, origin, scale):
    """The rotor disc, the station circle and the blades at true chord about
    their pitch axes, in the plane z = 0."""
    t = np.linspace(0.0, 2.0*np.pi, 200)
    circle = np.column_stack([np.cos(t), np.sin(t), np.zeros_like(t)])
    ax.add_patch(Polygon(iso(circle, origin, scale), closed=True, fc=PALETTE[0], alpha=0.06, ec="none", zorder=4))
    line(ax, *iso(circle, origin, scale).T, color=INK_SOFT, lw=0.8, ls="-.", zorder=5)
    line(ax, *iso(val("B_station")*circle, origin, scale).T, color=PALETTE[1], lw=1.0, ls="--", zorder=5)
    chord, root = val("B_chord")/val("B_radius"), val("B_root_end")
    for k in range(val("B_blades")):
        psi = np.radians(45.0) + 2.0*np.pi*k/val("B_blades")
        radial, normal = np.array([np.cos(psi), np.sin(psi), 0.0]), np.array([-np.sin(psi), np.cos(psi), 0.0])
        corners = blade_corners(radial, normal, root, 1.0, chord)
        ax.add_patch(Polygon(iso(corners, origin, scale), closed=True, fc=PALETTE[0], ec=INK, lw=0.8, zorder=6))
        line(ax, *iso([[0, 0, 0], root*radial], origin, scale).T, lw=1.2, zorder=6)
    ax.add_patch(Polygon(iso(0.035*circle, origin, scale), closed=True, fc="white", ec=INK, lw=1.0, zorder=7))


def sheet2():
    keys = ["B_radius", "B_radius_ft", "B_chord", "B_chord_in", "B_blades", "B_speed", "B_station", "B_root_end"]
    fig, ax = new_sheet(2, "ISOMETRIC VIEW", BANNER_SCHEMATIC, keys, "NOT TO SCALE (ISOMETRIC SKETCH)")
    origin, scale = (196.0, 184.0), 62.0
    record_body(2, body_iso(ax, origin, scale))
    rotor_iso(ax, origin, scale)
    hub, edge = iso([[0, 0, 0]], origin, scale)[0], iso([[0, 1, 0]], origin, scale)[0]
    arrow(ax, tuple(hub), tuple(edge))
    leader(ax, tuple(0.5*(hub + edge)), (20.0, 226.0), "rotor radius " + label_radius())
    leader(ax, tuple(iso([[0.5, -0.5, 0]], origin, scale)[0] + [0.0, -1.5]), (296.0, 138.0),
           "blade chord " + label_chord_b())
    leader(ax, tuple(iso([[-0.53, 0.53, 0]], origin, scale)[0]), (20.0, 124.0),
           "analysis station " + label_station())
    view_title(ax, origin[0], 264.0, "ISOMETRIC VIEW")
    note(ax, 256.0, 118.0,
         "Case B illustration. Rotor: %s blades, %s. Dashed circle: analysis station at %s. "
         % (label("B_blades"), label("B_speed"), label_station()) + NOTE_ROTOR
         + " The body outline is generic and carries no dimension. " + NOTE_ROOT)
    save(fig, "sheet2_isometric.png")


# ---------------------------------------------------------------- sheet 3
def blade_planform(ax, x_axis, y_axis, mm_per_m, looking_outboard):
    """The blade planform to scale, seen from above: rotation axis at x_axis,
    pitch axis at height y_axis, leading edge uppermost. The cutting-plane
    arrows of section A-A point the way the section is looked at."""
    radius, chord = mm_per_m*val("B_radius"), mm_per_m*val("B_chord")
    x_root, x_tip, x_st = x_axis + val("B_root_end")*radius, x_axis + radius, x_axis + val("B_station")*radius
    y_le, y_te = y_axis + PITCH_AXIS*chord, y_axis - (1.0 - PITCH_AXIS)*chord
    ax.add_patch(Rectangle((x_root, y_te), x_tip - x_root, chord, fc=FILL, ec=INK, lw=1.2, zorder=2))
    line(ax, [x_axis, x_root], [y_axis, y_axis], lw=2.0, ls=(0, (4, 2)))
    line(ax, [x_root, x_tip + 6.0], [y_axis, y_axis], color=INK_SOFT, lw=0.7, ls="-.", zorder=3)
    line(ax, [x_axis, x_axis], [y_axis - 30.0, y_axis + 30.0], color=INK_SOFT, lw=0.8, ls="-.")
    ax.text(x_axis + 2.0, y_axis + 27.0, "rotation axis", fontsize=FS_SMALL, color=INK, va="center")
    ax.text(x_axis + 2.0, y_axis - 6.0, "root end:\nschematic", fontsize=FS_SMALL, color=INK,
            ha="left", va="top")
    ax.text(x_root + 0.30*(x_st - x_root), y_le + 2.5, "leading edge", fontsize=FS_SMALL, color=INK, ha="center",
            va="bottom")
    ax.text(x_root + 0.42*(x_st - x_root), y_te - 2.5, "trailing edge", fontsize=FS_SMALL, color=INK, ha="center",
            va="top")
    ax.text(x_tip + 8.0, y_axis, "pitch axis", fontsize=FS_SMALL, color=INK, ha="left", va="center")
    dim_h(ax, x_axis, x_tip, y_axis - 44.0, "rotor radius " + label_radius(), (y_axis - 32.0, y_te - 2.0),
          above=False)
    dim_h(ax, x_axis, x_st, y_axis + 44.0, "analysis station " + label_station(), (y_axis + 32.0, y_le + 18.0))
    dim_v(ax, y_te, y_le, x_tip + 9.0 + 26.0, "blade chord\n" + label_chord_b().replace(" (", "\n("), x_tip + 1.5)
    sense = 1.0 if looking_outboard else -1.0
    for y, direction in ((y_le + 5.0, 1.0), (y_te - 5.0, -1.0)):
        line(ax, [x_st, x_st], [y, y + 9.0*direction], color=PALETTE[1], lw=2.0)
        arrow(ax, (x_st, y + 9.0*direction), (x_st + 10.0*sense, y + 9.0*direction), style="-|>")
        ax.text(x_st - 3.5*sense, y + 6.5*direction, "A", fontsize=FS_VIEW, fontweight="bold", color=PALETTE[1],
                ha="right" if sense > 0 else "left", va="center")
    line(ax, [x_st, x_st], [y_te, y_le], color=PALETTE[1], lw=1.0, ls="--", zorder=3)


def section_enlarged(ax, x_left, y_mid, chord_draw, looking_outboard):
    """Section A-A, enlarged, with its chord dimensioned, the leading edge on
    the side the viewing direction puts it; returns that side."""
    side = section_side(looking_outboard)
    xs, ys = section_outline()
    xs = xs if side == "left" else 1.0 - xs
    ax.add_patch(Polygon(np.column_stack([x_left + chord_draw*xs, y_mid + chord_draw*ys]), closed=True, fc=FILL,
                         ec=INK, lw=1.2, zorder=2))
    x_le, x_te = (x_left, x_left + chord_draw) if side == "left" else (x_left + chord_draw, x_left)
    ax.text(x_le, y_mid + 0.07*chord_draw, "leading edge", fontsize=FS_SMALL, color=INK,
            ha="left" if side == "left" else "right", va="bottom")
    ax.text(x_te, y_mid + 0.03*chord_draw, "trailing edge", fontsize=FS_SMALL, color=INK,
            ha="right" if side == "left" else "left", va="bottom")
    dim_h(ax, x_left, x_left + chord_draw, y_mid - 0.06*chord_draw - 9.0, "chord " + label_chord_b(),
          (y_mid - 2.0, y_mid - 2.0), above=False)
    view_title(ax, x_left + chord_draw/2, y_mid + 0.06*chord_draw + 12.0,
               "SECTION A-A, LOOKING %s (ENLARGED; %s)"
               % ("OUTBOARD" if looking_outboard else "INBOARD", label("A_section")))
    return side


def sheet3():
    keys = ["B_radius", "B_radius_ft", "B_chord", "B_chord_in", "B_blades", "B_speed", "B_station", "B_root_end",
            "A_section"]
    mm_per_m = 300.0/val("B_radius")
    fig, ax = new_sheet(3, "MAIN ROTOR BLADE", "PLANFORM DRAWN TO SCALE; ILLUSTRATION ONLY", keys,
                        "PLANFORM 1 : %.1f ON %s; SECTION A-A ENLARGED" % (1000.0/mm_per_m, SHEET_FORMAT))
    looking_outboard = True
    view_title(ax, 190.0, 270.0, "BLADE PLANFORM FROM ABOVE (ONE OF %s BLADES)" % label("B_blades"))
    blade_planform(ax, 30.0, 208.0, mm_per_m, looking_outboard)
    side = section_enlarged(ax, 36.0, 128.0, 150.0, looking_outboard)
    SHEET_RECORD.append(dict(sheet=3, item="section A-A: cutting-plane arrows point",
                             value="outboard" if looking_outboard else "inboard"))
    SHEET_RECORD.append(dict(sheet=3, item="section A-A: leading edge drawn on the", value=side))
    SHEET_RECORD.append(dict(sheet=3, item="section A-A: leading edge belongs on the",
                             value=section_side(looking_outboard)))
    note(ax, 256.0, 150.0,
         "Case B illustration: an aerofoil at the blade station %s of a rotor turning at %s. "
         % (label_station(), label("B_speed")) + NOTE_ROTOR
         + " The planform is drawn about its pitch axis, taken at %s of the chord as in Case A. " % PITCH_AXIS
         + NOTE_ROOT + " The attachment is not drawn.")
    save(fig, "sheet3_main_rotor_blade.png")


# ---------------------------------------------------------------- sheet 4
def hatched_section(ax, x_le, y_mid, chord_draw):
    """A cut section: light fill with coloured hatching and a navy outline."""
    xs, ys = section_outline()
    pts = np.column_stack([x_le + chord_draw*xs, y_mid + chord_draw*ys])
    ax.add_patch(Polygon(pts, closed=True, fc="#eef3f9", ec=PALETTE[0], hatch="///", lw=0.0, zorder=2))
    ax.add_patch(Polygon(pts, closed=True, fc="none", ec=INK, lw=1.3, zorder=3))
    line(ax, [x_le - 8.0, x_le + chord_draw + 8.0], [y_mid, y_mid], color=INK_SOFT, lw=0.7, ls="-.", zorder=4)


def case_a_section(ax, x_le, y_mid, mm_per_m):
    """The Case A tunnel-model section with chord, maximum thickness and pitch axis."""
    chord = mm_per_m*val("A_chord")
    half = 0.5*mm_per_m*val("A_thickness")
    hatched_section(ax, x_le, y_mid, chord)
    x_axis = x_le + val("A_pitch_axis")*chord
    ax.add_patch(Circle((x_axis, y_mid), 1.6, fc="white", ec=PALETTE[1], lw=1.4, zorder=5))
    line(ax, [x_axis, x_axis], [y_mid - half - 4.0, y_mid + half + 12.0], color=PALETTE[1], lw=0.8, ls="-.")
    dim_h(ax, x_le, x_axis, y_mid + half + 10.0, "pitch axis at %s" % label("A_pitch_axis"),
          (y_mid + 2.0, y_mid + half + 12.0))
    dim_h(ax, x_le, x_le + chord, y_mid - half - 12.0, "chord " + label("A_chord"), (y_mid - 2.0, y_mid - 2.0),
          above=False)
    x_max = x_le + val("A_tmax_station")*chord
    line(ax, [x_max, x_le + chord + 12.0], [y_mid + half, y_mid + half], lw=0.6)
    line(ax, [x_max, x_le + chord + 12.0], [y_mid - half, y_mid - half], lw=0.6)
    arrow(ax, (x_le + chord + 10.0, y_mid - half), (x_le + chord + 10.0, y_mid + half))
    ax.text(x_le + chord + 13.0, y_mid, "maximum\nthickness\n%s\n(%s)\nat %s" % (
        label("A_thickness"), label("A_thickness_ratio"), label("A_tmax_station")), fontsize=FS_TEXT, color=INK,
        ha="left", va="center")
    view_title(ax, x_le + chord/2, y_mid + half + 26.0,
               "CASE A: TUNNEL MODEL SECTION, %s" % label("A_section"))


def case_b_section(ax, x_le, y_mid, mm_per_m):
    """The Case B illustration section at the same scale, chord only."""
    chord = mm_per_m*val("B_chord")
    hatched_section(ax, x_le, y_mid, chord)
    dim_h(ax, x_le, x_le + chord, y_mid - 0.06*chord - 10.0, "chord " + label_chord_b(),
          (y_mid - 2.0, y_mid - 2.0), above=False)
    view_title(ax, x_le + chord/2, y_mid + 0.06*chord + 8.0,
               "CASE B: SECTION A-A OF SHEET 3 (%s AS ILLUSTRATION)" % label("A_section"))


def sheet4():
    keys = ["A_section", "A_chord", "A_thickness_ratio", "A_thickness", "A_tmax_station", "A_pitch_axis", "B_chord",
            "B_chord_in"]
    mm_per_m = 270.0/val("A_chord")
    fig, ax = new_sheet(4, "SECTION A-A, AEROFOIL", "BOTH SECTIONS DRAWN TO ONE SCALE; ILLUSTRATION ONLY", keys,
                        "1 : %.2f ON %s" % (1000.0/mm_per_m, SHEET_FORMAT))
    case_a_section(ax, 40.0, 214.0, mm_per_m)
    case_b_section(ax, 40.0, 130.0, mm_per_m)
    note(ax, 256.0, 104.0,
         "Outline: NACA four-digit thickness form with a closed trailing edge, symmetric. Case A is the model of "
         "the tunnel experiment (%s). Case B uses the same outline as an illustration; the real blade of the rotor "
         "named on sheet 3 uses SC1095 and SC1094 R8 sections." % SRC_TUNNEL)
    save(fig, "sheet4_section_AA_airfoil.png")


def main():
    apply_style()
    plt.rcParams["hatch.linewidth"] = 0.5
    write_dimension_table()
    sheet1()
    sheet2()
    sheet3()
    sheet4()
    pd.DataFrame(SHEET_RECORD).to_csv(HERE/"sheet_record.csv", index=False)
    pm.record_figures(HERE/"figure_record.csv", FIGURE_RECORD)
    print("[drawings] four sheets, dimensions.csv and sheet_record.csv written to %s" % HERE)


if __name__ == "__main__":
    main()
