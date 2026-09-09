# -*- coding: utf-8 -*-
"""
07_report / build_pdfs.py
-------------------------
Generates the two PDF deliverables that accompany case.docx:
  (1) UNISTALL_plots_album.pdf  — every figure / curve / contour / chart /
                                  engineering drawing, one per page, captioned.
  (2) UNISTALL_data_dossier.pdf — the input, metrics and validation tables in
                                  full, the large per-step outputs as captioned
                                  samples, and the solver config. Its second
                                  page is a contents list naming exactly what is
                                  in it and what is not (the reconstructed 2-D
                                  fields and the digitised experimental loops
                                  are too large to typeset and ship as CSV).
Author: Akosa Samuel Onyejekwe (independent).  No black is used anywhere.
"""
import sys, glob, json, textwrap
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from aero_style import apply_style, INK, INK_SOFT
apply_style()
# aero_style sets savefig.bbox="tight", which is right for the standalone PNGs
# but wrong here: PdfPages.savefig() honours it and crops every page to its own
# content, so these two documents had NO consistent page size -- 45 distinct
# geometries across the 104 album pages, the smallest an 80 x 36 pt sliver
# holding one section title. Every page here must be a full 11 x 8.5 sheet.
#
# It has to be cleared on the rcParam. Passing savefig(bbox_inches=None) does
# NOT work: matplotlib treats an explicit None as "unspecified" and falls back
# to this very rcParam, so the crop survives. Verified both ways.
matplotlib.rcParams["savefig.bbox"] = None
PAGE = {}
# Type 3 (matplotlib's default) draws each glyph as inline PDF operators. The
# text renders correctly but is not reliably selectable or searchable, which for
# a document whose dossier IS tables makes the data unsearchable. Type 42 embeds
# proper TrueType subsets instead.
matplotlib.rcParams["pdf.fonttype"] = 42
from project_meta import AUTHOR, STUDY_DATE        # single source of truth
ACC = "#2e6fb7"

def cover(pdf, title, subtitle):
    fig = plt.figure(figsize=(11, 8.5)); fig.patch.set_facecolor("white")
    fig.text(0.5, 0.70, "UNISTALL™ Dynamic-Stall Case Study", ha="center",
             fontsize=24, color=ACC, weight="bold")
    fig.text(0.5, 0.60, title, ha="center", fontsize=20, color=INK, weight="bold")
    fig.text(0.5, 0.53, subtitle, ha="center", fontsize=13, color=INK)
    fig.text(0.5, 0.40, "Prediction of Dynamic Stall on a Helicopter\n"
             "Main-Rotor Retreating Blade", ha="center", fontsize=14, color=INK)
    fig.text(0.5, 0.22, f"Author: {AUTHOR} (independent)", ha="center",
             fontsize=13, color=INK, weight="bold")
    fig.text(0.5, 0.17, f"Date: {STUDY_DATE}", ha="center", fontsize=11, color=INK)
    pdf.savefig(fig, **PAGE); plt.close(fig)

def section(pdf, text):
    fig = plt.figure(figsize=(11, 8.5))
    fig.text(0.5, 0.5, text, ha="center", va="center", fontsize=22,
             color=ACC, weight="bold")
    pdf.savefig(fig, **PAGE); plt.close(fig)

# ======================================================= 1. PLOTS ALBUM
ALBUM = [
 ("Geometry", [ROOT/"01_geometry"/"fig_geometry_profile.png",
               ROOT/"01_geometry"/"fig_geometry_thickness.png"]),
 ("Mesh", sorted((ROOT/"02_mesh").glob("fig_*.png"))),
 ("Calibration & static polar", [ROOT/"06_postprocessing"/"plots"/"static_polar_calibration.png"]),
 ("Hysteresis loops", sorted((ROOT/"06_postprocessing"/"plots").glob("hyst_*.png"))),
 ("Time histories & states", sorted((ROOT/"06_postprocessing"/"plots").glob("timehist_*.png"))
                            + sorted((ROOT/"06_postprocessing"/"plots").glob("states_*.png"))),
 ("Surface pressure & convergence",
   sorted((ROOT/"06_postprocessing"/"plots").glob("cp_distribution_*.png"))
   + [ROOT/"06_postprocessing"/"plots"/"convergence_residuals.png"]),
 ("Pressure contours", sorted((ROOT/"06_postprocessing"/"plots").glob("contour_Cp_*.png"))),
 ("Velocity: streamlines & vectors",
   sorted((ROOT/"06_postprocessing"/"plots").glob("contour_speed_stream_*.png"))
   + sorted((ROOT/"06_postprocessing"/"plots").glob("contour_vectors_*.png"))),
 ("Vorticity (dynamic-stall vortex)",
   sorted((ROOT/"06_postprocessing"/"plots").glob("contour_vorticity_*.png"))),
 ("Local Mach contours", sorted((ROOT/"06_postprocessing"/"plots").glob("contour_Mach_*.png"))),
 ("Temperature fields & profiles",
   sorted((ROOT/"06_postprocessing"/"plots").glob("contour_Tstatic_*.png"))
   + sorted((ROOT/"06_postprocessing"/"plots").glob("contour_Trecovery_*.png"))
   + sorted((ROOT/"06_postprocessing"/"plots").glob("temperature_profile_*.png"))),
 ("3-D outputs", sorted((ROOT/"06_postprocessing"/"plots").glob("fig3d_*.png"))),
 ("Validation", sorted((ROOT/"06_postprocessing"/"validation").glob("fig_*.png"))),
 ("Engineering drawings", sorted((ROOT/"08_engineering_drawings").glob("sheet*.png"))),
]
out1 = ROOT/"07_report"/"UNISTALL_plots_album.pdf"
with PdfPages(out1) as pdf:
    cover(pdf, "Figures & Contours Album", "All graphs, curves, contours, charts and drawings")
    n = 0
    for title, files in ALBUM:
        files = [f for f in files if Path(f).exists()]
        if not files: continue
        section(pdf, title)
        for f in files:
            img = plt.imread(str(f))
            fig = plt.figure(figsize=(11, 8.5))
            ax = fig.add_axes([0.04, 0.09, 0.92, 0.85]); ax.imshow(img); ax.axis("off")
            fig.text(0.5, 0.045, Path(f).stem.replace("_", " "), ha="center",
                     fontsize=10, color=ACC, style="italic")
            fig.text(0.5, 0.02, f"{title}", ha="center", fontsize=8, color=INK_SOFT)
            pdf.savefig(fig, **PAGE); plt.close(fig); n += 1
print(f"[pdf] plots album: {n} figure pages -> {out1.name}")

# ======================================================= 2. DATA DOSSIER
# max_rows was 34, which cut the 36-row response-surface tables to 34 and
# declared "first 34 of 36 rows" -- hiding two rows of a complete design
# space for no gain. 38 shows them whole; measured to add no new overflow
# (blocks outside the page rectangle unchanged) and no extra pages.
def table_page(pdf, df, title, max_rows=38, max_cols=10):
    """Render one table. A table WIDER than max_cols is split across pages, each
    carrying the first (identifier) column plus a block of the rest, rather than
    truncated: the old "first 10 of 16 cols" quietly dropped CLmax_exp,
    CMmin_model, CMmin_exp and both Xi_hat columns from validation_nasa_real.csv
    -- i.e. every measured value the model is being compared against."""
    if df.shape[1] > max_cols:
        key = df.columns[0]
        rest = list(df.columns[1:])
        blocks = [rest[i:i+max_cols-1] for i in range(0, len(rest), max_cols-1)]
        for bi, cols in enumerate(blocks, 1):
            _table_page_one(pdf, df[[key] + cols], title, max_rows,
                            f"columns {bi} of {len(blocks)}: {key} + "
                            f"{cols[0]}..{cols[-1]}")
        return
    _table_page_one(pdf, df, title, max_rows, None)


def _table_page_one(pdf, df, title, max_rows, block_note):
    d = df.copy()
    extra = []
    if block_note:
        extra.append(block_note)
    if len(d) > max_rows:
        extra.append(f"first {max_rows} of {len(df)} rows"); d = d.head(max_rows)
    # Missing entries as an em dash, not the string "nan" -- see build_docx.py.
    # This MUST happen element-wise before astype(str): on pandas 3.0
    # DataFrame.replace() treats a dict's keys as COLUMN names, and even
    # replace([...], value) leaves an already-stringified "nan" alone, so both
    # of those quietly did nothing and the cells still read "nan".
    na = bool(d.isna().to_numpy().any())
    d = d.map(lambda v: "\u2014" if (v is None or (isinstance(v, float) and v != v))
              else str(v))
    if na:
        extra.append("\u2014 = not applicable / no value defined")
    fig = plt.figure(figsize=(11, 8.5)); ax = fig.add_axes([0.03, 0.03, 0.94, 0.88]); ax.axis("off")
    # Fit the title to the page. At a fixed 13 pt the longest of these ran 5.8 pt
    # off the right edge and was clipped mid-word, so the page that says
    # "first 34 of 721 rows" actually read "first 34 of 721 rov". Measure the
    # rendered width and step the size down until it fits, which changes nothing
    # for the titles that already did.
    _title = title + ("   (" + "; ".join(extra) + ")" if extra else "")
    _fs = 13
    _r = fig.canvas.get_renderer()
    for _fs in (13, 12, 11, 10, 9, 8):
        _t = ax.set_title(_title, color=ACC, fontsize=_fs, weight="bold",
                          pad=14, loc="left")
        if _t.get_window_extent(renderer=_r).x1 <= fig.bbox.x1 - 6:
            break
    tbl = ax.table(cellText=d.values, colLabels=list(d.columns),
                   loc="upper center", cellLoc="center")
    tbl.auto_set_font_size(False)
    fs = 8 if d.shape[1] <= 7 else 6.5
    tbl.set_fontsize(fs); tbl.scale(1, 1.35)
    for (r, _c), cell in tbl.get_celld().items():
        cell.set_edgecolor(INK_SOFT)
        if r == 0:
            cell.set_facecolor("#dfe8f2"); cell.get_text().set_color(INK)
            cell.get_text().set_weight("bold")
        else:
            cell.set_facecolor("white" if r % 2 else "#f3f6fa")
            cell.get_text().set_color(INK)
    try: tbl.auto_set_column_width(col=list(range(len(d.columns))))
    except Exception: pass
    pdf.savefig(fig, **PAGE); plt.close(fig)

CSV_GROUPS = [
 ("Inputs — geometry & flow",
   [ROOT/"01_geometry"/"section_geometry_summary.csv",
    ROOT/"03_model_setup"/"flow_conditions.csv",
    ROOT/"03_model_setup"/"kinematics.csv",
    ROOT/"03_model_setup"/"material_thermo_properties.csv",
    ROOT/"03_model_setup"/"static_polar_reference.csv"]),
 ("Inputs — mesh",
   [ROOT/"02_mesh"/"mesh_quality_metrics.csv",
    ROOT/"02_mesh"/"mesh_radial_spacing.csv"]),
 ("Outputs — engineering metrics",
   [ROOT/"05_solution"/"metrics_A_validation.csv",
    ROOT/"05_solution"/"metrics_B_application.csv",
    ROOT/"05_solution"/"summary_all_cases.csv"]),
 ("Outputs — time histories & fields (samples)",
   [ROOT/"05_solution"/"time_history_A_validation.csv",
    ROOT/"05_solution"/"time_history_B_application.csv",
    ROOT/"05_solution"/"cp_distribution_A_validation.csv",
    ROOT/"05_solution"/"model_static_polar.csv",
    ROOT/"05_solution"/"response_surface.csv"]),
 ("Outputs — convergence",
   [ROOT/"05_solution"/"convergence"/"residuals_A_validation.csv",
    ROOT/"05_solution"/"convergence"/"residuals_B_application.csv"]),
 ("Validation & calibration",
   [ROOT/"06_postprocessing"/"validation"/"validation_static.csv",
    ROOT/"06_postprocessing"/"validation"/"calibration_constants.csv",
    ROOT/"06_postprocessing"/"validation"/"validation_nasa_real.csv",
    ROOT/"06_postprocessing"/"validation"/"validation_realdata_summary.csv"]),
]
out2 = ROOT/"07_report"/"UNISTALL_data_dossier.pdf"
with PdfPages(out2) as pdf:
    cover(pdf, "Data & Tables Dossier",
          "Inputs, metrics and validation tables in full; the large per-step\n"
          "outputs as captioned samples (see the contents page)")
    # Contents page. The cover used to promise "All CSVs, metrics and tables",
    # which was not true: the reconstructed 2-D fields, the digitised
    # experimental loops and the per-phase Cp of case B are not in here. Say so,
    # and say where they are.
    fig = plt.figure(figsize=(11, 8.5)); ax = fig.add_axes([0.06, 0.05, 0.88, 0.86])
    ax.axis("off")
    ax.set_title("What this dossier contains", color=ACC, fontsize=15,
                 weight="bold", loc="left", pad=16)
    _inc = "\n".join("    \u2022 " + t for t, _ in CSV_GROUPS)
    # Which tables are SAMPLED is derived from the row counts against
    # table_page's max_rows, not written by hand. The hand-written list had
    # already drifted: raising max_rows to 38 made response_surface.csv (36
    # rows) print in full while this page still called it a sample, and
    # model_static_polar.csv (89 rows) became a sample without being listed.
    _MAXROWS = table_page.__defaults__[0]
    _sampled = []
    for _t, _fs in CSV_GROUPS:
        for _f in _fs:
            if Path(_f).exists() and len(pd.read_csv(_f)) > _MAXROWS:
                _sampled.append((Path(_f).name, len(pd.read_csv(_f))))
    _sampled_txt = "".join("    \u2022 %s  \u2014 %d rows\n" % nr for nr in _sampled) \
                   or "    \u2022 (none: every table below prints in full)\n"
    _n_field = len(sorted((ROOT/"05_solution").glob("field_*.csv")))
    # row count DERIVED from the grid the config declares, not typed in. It was
    # the literal "37 400", which is nx*ny for the current grid but would have
    # gone on saying so after either was changed.
    _fc = json.load(open(ROOT/"03_model_setup"/"solver_config.json"))["field_reconstruction"]
    _field_rows = f"{_fc['grid_nx_solution']*_fc['grid_ny_solution']:,}".replace(",", "\u202f")
    _n_exp = len(sorted((ROOT/"06_postprocessing"/"validation").glob("exp_frame_*.csv")))
    ax.text(0, 0.98,
            "IN FULL, one page per table (wide tables continue over further pages):\n"
            + _inc
            + "\n\nAS CAPTIONED SAMPLES (row count given on each page):\n"
            + _sampled_txt
            + "\nNOT REPRODUCED HERE \u2014 too large to typeset, shipped as CSV in the "
              "repository:\n"
              f"    \u2022 05_solution/field_*.csv  \u2014 {_n_field} reconstructed 2-D fields, "
              f"{_field_rows} rows each\n"
              "    \u2022 05_solution/cp_distribution_B_application.csv  \u2014 the case-B "
              "counterpart of the sampled case-A table\n"
              f"    \u2022 06_postprocessing/validation/exp_frame_*.csv  \u2014 {_n_exp} "
              "digitised experimental loops (NASA TM-84245)\n"
              "\nEvery figure is in the companion volume UNISTALL_plots_album.pdf.",
            va="top", ha="left", fontsize=10.5, color=INK, linespacing=1.55)
    pdf.savefig(fig, **PAGE); plt.close(fig)
    # solver config page
    cfg = json.load(open(ROOT/"03_model_setup"/"solver_config.json"))
    fig = plt.figure(figsize=(11, 8.5)); ax = fig.add_axes([0.05, 0.04, 0.9, 0.9]); ax.axis("off")
    ax.set_title("Solver configuration (solver_config.json)", color=ACC, fontsize=13,
                 weight="bold", loc="left")
    # textwrap.wrap() collapses newlines, so wrapping the whole blob destroys the
    # JSON indentation. Wrap each line individually and keep its leading indent.
    lines = []
    for ln in json.dumps(cfg, indent=2).splitlines():
        indent = " " * (len(ln) - len(ln.lstrip()))
        lines += textwrap.wrap(ln, 108, initial_indent="", subsequent_indent=indent + "    ",
                               drop_whitespace=False, replace_whitespace=False) or [""]
    ax.text(0, 0.99, "\n".join(lines), va="top", ha="left", fontsize=6.2,
            family="monospace", color=INK, transform=ax.transAxes)
    pdf.savefig(fig, **PAGE); plt.close(fig)
    np_ = 0
    for title, files in CSV_GROUPS:
        section(pdf, title)
        for f in files:
            if not Path(f).exists(): continue
            try:
                table_page(pdf, pd.read_csv(f), Path(f).name); np_ += 1
            except Exception as e:
                section(pdf, f"{Path(f).name}\n[error: {e}]")
# np_ counts TABLES, not pages: a wide table now continues over several pages,
# so calling it "table pages" would understate the document.
print(f"[pdf] data dossier: {np_} tables -> {out2.name}")

# copies at project root
import shutil
shutil.copy(str(out1), str(ROOT/out1.name)); shutil.copy(str(out2), str(ROOT/out2.name))
print("[pdf] copied both PDFs to project root")
