# -*- coding: utf-8 -*-
# Run from the repository root:  PYTHONPATH=. python3 00_overview/write_folder_contents.py
"""
00_overview / write_folder_contents.py
--------------------------------------
Writes 00_overview/folder_contents.md: for each numbered folder, what each
kind of file is and which script writes it, and a short table of equivalent
names for a reader who knows a flow-solver layout of such a case.

Author: Akosa Samuel Onyejekwe (independent)

ENTRIES is the list. A file of the numbered folders that no entry covers is
an error: `uncovered` names it, the script stops, and check_case_study.py
makes the same test, so the list cannot fall behind the folders.
"""
import fnmatch
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
FOLDERS = ("00_overview", "01_geometry", "02_mesh", "03_model_setup", "04_solver", "05_solution",
           "06_postprocessing", "07_report", "08_engineering_drawings")
IGNORED = ("__pycache__",)

SETUP, RUN, MESH, MESHF = ("03_model_setup/generate_setup.py", "04_solver/run_case.py", "02_mesh/generate_mesh.py",
                           "02_mesh/field_on_mesh.py")
PLOTS, PLOTS3, VALID = ("06_postprocessing/make_all_plots.py", "06_postprocessing/make_3d_plots.py",
                        "06_postprocessing/validation/make_validation.py")
REPORT, DRAW, GEOM = "07_report/build_report.py", "08_engineering_drawings/draw_engineering.py", \
    "01_geometry/generate_geometry.py"
BY_HAND = "written by hand"

# (pattern relative to the repository root, what the files are, what writes them)
ENTRIES = (
    ("00_overview/case_definition.md", "the definition of the study in words, every number from 03_model_setup", SETUP),
    ("00_overview/folder_contents.md", "this list", "00_overview/write_folder_contents.py"),
    ("00_overview/write_folder_contents.py", "the script that writes this list", BY_HAND),
    ("01_geometry/generate_geometry.py", "script of the section geometry", BY_HAND),
    ("01_geometry/naca*_coordinates.csv", "coordinates of the closed section outline, fractions of chord", GEOM),
    ("01_geometry/section_geometry_summary.csv", "thickness, nose radius, form of the trailing edge, area", GEOM),
    ("01_geometry/fig_geometry_*.png", "the profile and its thickness distribution", GEOM),
    ("0[1268]_*/figure_record.csv", "resolution and smallest lettering of the folder's figures, from which the "
     "report works out the size of the lettering as printed", "the figure script(s) of the folder"),
    ("06_postprocessing/validation/figure_record.csv", "the same for the two validation figures", VALID),
    ("02_mesh/README.md", "standing of the grid: an illustration, no flow equation solved on it", BY_HAND),
    ("02_mesh/generate_mesh.py", "script of the O-grid", BY_HAND),
    ("02_mesh/field_on_mesh.py", "script that evaluates the reconstructed field at the grid nodes", BY_HAND),
    ("02_mesh/mesh_*.csv", "grid nodes, quality measures and the wall-normal spacing law", MESH),
    ("02_mesh/fig_mesh_field_*.png", "corrected pressure coefficient of the reconstruction at the grid nodes, at "
     "peak lift", MESHF),
    ("02_mesh/fig_mesh_*.png", "pictures of the grid", MESH),
    ("03_model_setup/generate_setup.py", "script of the conditions of the two cases", BY_HAND),
    ("03_model_setup/flow_conditions.csv", "chord, Mach number, stream speed, air state, Reynolds number and "
     "reduced frequency of each case, each row with its source", SETUP),
    ("03_model_setup/kinematics.csv", "the prescribed pitching motion of each case", SETUP),
    ("03_model_setup/station_condition.csv", "how the Mach number and reduced frequency of Case B follow from the "
     "rotor", SETUP),
    ("03_model_setup/air_properties.csv", "the gas properties the cases are run with", SETUP),
    ("03_model_setup/static_inputs.csv", "the static normal force, moment and separation point the model reads at "
     "each case's Mach number; they come from measurement", SETUP),
    ("03_model_setup/solver_config.json", "solver name, calibrated constants, march of each case, settings of the "
     "reconstruction", SETUP),
    ("04_solver/README.md", "what the folder holds: drivers only; the model is the package unistall/", BY_HAND),
    ("04_solver/run_case.py", "driver of the two cases: writes 05_solution", BY_HAND),
    ("04_solver/unistall_solver.py", "command-line entry to the packaged solver for one condition", BY_HAND),
    ("05_solution/README.md", "standing of the files: loads are the model's prediction, fields a reconstruction",
     BY_HAND),
    ("05_solution/time_history_*.csv", "loads and model states at every step of the reported cycle", RUN),
    ("05_solution/metrics_*.csv", "scalar results of a case", RUN),
    ("05_solution/summary_all_cases.csv", "one row of main results per case", RUN),
    ("05_solution/held_out_error_bands.csv", "error of the model in each headline quantity on the held-out loops",
     RUN),
    ("05_solution/static_station_spread.csv", "Case B with the static inputs of each Mach station", RUN),
    ("05_solution/model_static_polar_*.csv", "lift, drag, moment and separation point of the model in a slow "
     "sweep: its quasi-steady limit, not a measurement", RUN),
    ("05_solution/response_surface.csv", "peak lift, minimum moment and damping over mean incidence and reduced "
     "frequency", RUN),
    ("05_solution/runtime_environment.csv", "machine, library versions, CPU time of the march", RUN),
    ("05_solution/cp_distribution_*.csv", "reconstructed surface pressure, local Mach number and temperatures at "
     "six instants", RUN),
    ("05_solution/reconstruction_*.csv", "the reconstruction instant by instant as numbers", RUN),
    ("05_solution/cycle_scan_*.csv", "the reconstruction at instants spread over the cycle", RUN),
    ("05_solution/field_*_on_mesh.csv.gz", "the reconstruction at the nodes of the O-grid, gzip-compressed CSV", MESHF),
    ("05_solution/field_*.csv.gz", "the reconstructed field on a rectangular grid at one instant, gzip-compressed "
     "CSV", RUN),
    ("05_solution/convergence/residuals_*.csv", "change of the cycle with the number of cycles marched", RUN),
    ("05_solution/convergence/timestep_*.csv", "loads against the time step; observed order and discretisation "
     "uncertainty", RUN),
    ("05_solution/convergence/panel_convergence_*.csv", "the reconstruction against the number of panels", RUN),
    ("05_solution/convergence/reconstruction_timing.csv", "CPU time of the reconstruction", RUN),
    ("06_postprocessing/README.md", "standing of the figures", BY_HAND),
    ("06_postprocessing/make_all_plots.py", "script of the two-dimensional figures", BY_HAND),
    ("06_postprocessing/make_3d_plots.py", "script of the surfaces, maps over two variables and the vector view",
     BY_HAND),
    ("06_postprocessing/overlay_record.csv", "whether the measured loop was drawn with each case", PLOTS),
    ("06_postprocessing/masked_region_record.csv", "the region beyond the critical pressure printed on each figure "
     "of a corrected quantity", PLOTS),
    ("06_postprocessing/plots/fig3d_*.png", "response surface; surfaces of the reconstructed pressure and speed; "
     "section with velocity vectors", PLOTS3),
    ("06_postprocessing/plots/cp_phase_map_*.png", "reconstructed upper-surface pressure over chord and cycle",
     PLOTS3),
    ("06_postprocessing/plots/contour_*.png", "maps of the reconstructed field: pressure, local Mach number, "
     "static and recovery temperature, speed, vorticity, vectors", PLOTS),
    ("06_postprocessing/plots/*.png", "load loops, time histories, model states, static inputs, convergence, "
     "surface pressure and surface temperatures", PLOTS),
    ("06_postprocessing/validation/make_validation.py", "script of the validation folder", BY_HAND),
    ("06_postprocessing/validation/README.md", "what each file of the validation folder is", VALID),
    ("06_postprocessing/validation/PROVENANCE.md", "where the measured frames come from", VALID),
    ("06_postprocessing/validation/validation_*.csv", "static comparison, one row per scored loop, and the targets "
     "with their verdicts; the same numbers as results/", VALID),
    ("06_postprocessing/validation/loop_*.csv", "the model read at the measured points of the loops drawn", VALID),
    ("06_postprocessing/validation/calibration_constants.csv", "the calibrated and fixed constants", VALID),
    ("06_postprocessing/validation/fig_validation_*.png", "static comparison; measured against predicted loops",
     VALID),
    ("07_report/build_report.py", "the one builder of the Word file and the three PDFs", BY_HAND),
    ("07_report/rpt_*.py", "the parts of the builder: data, text, layout, album, dossier", BY_HAND),
    ("07_report/case.docx", "the report as a Word document", REPORT),
    ("07_report/UNISTALL_*.pdf", "the report, the plots album and the data dossier", REPORT),
    ("07_report/_equations/eq_*.png", "the equations of the report", REPORT),
    ("07_report/report_numbers.json", "every number quoted in the report with the file it was read from", REPORT),
    ("07_report/figure_lettering.csv", "smallest lettering of each figure of the report as printed", REPORT),
    ("07_report/page_fill.csv", "how full each page of the three PDFs is", REPORT),
    ("08_engineering_drawings/draw_engineering.py", "script of the drawing sheets", BY_HAND),
    ("08_engineering_drawings/sheet?_*.png", "four illustrative sheets, not for manufacture", DRAW),
    ("08_engineering_drawings/dimensions.csv", "every dimension printed on a sheet, with its source", DRAW),
    ("08_engineering_drawings/sheet_record.csv", "what each sheet states: scale, viewing direction of section A-A, "
     "body coordinates", DRAW),
)

# what a reader who knows a flow-solver layout would look for -> where it is here
EQUIVALENTS = (
    ("fluid or material properties", "03_model_setup/air_properties.csv holds the gas properties"),
    ("aerofoil polar or tabulated section data", "03_model_setup/static_inputs.csv holds the static inputs the "
     "model reads; they come from measurement. 05_solution/model_static_polar_<case>.csv is the model's own "
     "slow-sweep polar"),
    ("solver executable or solver script", "the model is the package unistall/; 04_solver/unistall_solver.py runs "
     "it for one condition and 04_solver/run_case.py for the two cases"),
    ("mesh and solution on the mesh", "02_mesh is an illustration only; no flow equation is solved. "
     "05_solution/field_*_on_mesh.csv.gz is the reconstruction evaluated at its nodes"),
    ("field or contour data", "05_solution/field_*.csv.gz, compressed (gzip CSV), read with pandas.read_csv"),
    ("validation data and comparison", "06_postprocessing/validation, the same numbers as results/"),
    ("separate report, album and dossier builders", "one builder, 07_report/build_report.py, writes the Word file "
     "and the three PDFs"),
)


def files() -> list:
    """Every file of the numbered folders, relative to the repository root."""
    found = []
    for folder in FOLDERS:
        found += [str(p.relative_to(ROOT)) for p in sorted((ROOT/folder).rglob("*"))
                  if p.is_file() and not any(part in IGNORED for part in p.parts)]
    return found


def covering(path: str) -> tuple | None:
    """The first entry whose pattern covers a file, or None."""
    return next((e for e in ENTRIES if fnmatch.fnmatchcase(path, e[0])), None)


def uncovered() -> list:
    """The files of the numbered folders that no entry covers."""
    return [f for f in files() if covering(f) is None]


def text() -> str:
    """The text of folder_contents.md."""
    counts = {e[0]: 0 for e in ENTRIES}
    for f in files():
        counts[covering(f)[0]] += 1
    lines = ["# Contents of the numbered folders", "",
             "Written by `00_overview/write_folder_contents.py`. Each row is a kind of file: the pattern of its "
             "name, what it is, and what writes it. n is the number of files of that kind present when this "
             "list was written.", ""]
    for folder in FOLDERS:
        lines += [f"## {folder}", "", "| Files | n | What they are | Written by |", "|---|---|---|---|"]
        for pattern, what, by in ENTRIES:
            inside = pattern.startswith(folder + "/") or (pattern.startswith("0[1268]_*/") and folder[:2] in
                                                         ("01", "02", "06", "08"))
            if inside:
                shown = pattern.split("/", 1)[1]
                n = counts[pattern] if not pattern.startswith("0[") else sum(
                    1 for f in files() if f.startswith(folder + "/") and covering(f)[0] == pattern)
                lines.append(f"| `{shown}` | {n} | {what} | {by if by == BY_HAND else '`' + by + '`' if '/' in by else by} |")
        lines.append("")
    lines += ["## Where to find what", "",
              "For a reader who expects the layout of a flow-solver case:", "",
              "| Looking for | Here |", "|---|---|"]
    lines += [f"| {a} | {b} |" for a, b in EQUIVALENTS]
    lines += ["", "## The two case names", "",
              "File names carry the keys `A_validation` and `B_application`. They are labels of the two conditions "
              "and nothing more. Case A is the tunnel condition of a measured loop that was used in the calibration, "
              "so the comparison drawn for it is a fit and not an independent test; Case B is a section in a steady "
              "stream at a blade station's Mach number and reduced frequency, not a rotor calculation. The comparison "
              "of the model with loops held out of the calibration is in `results/` and in "
              "`06_postprocessing/validation/`."]
    return "\n".join(lines) + "\n"


def main() -> int:
    """Write the list, or stop naming the files no entry covers."""
    target = HERE/"folder_contents.md"
    target.touch()
    missing = uncovered()
    if missing:
        print("[folder contents] no entry covers: " + ", ".join(missing))
        return 1
    target.write_text(text(), encoding="utf-8")
    print("[folder contents] %d files of %d kinds listed in %s" % (len(files()), len(ENTRIES), target.relative_to(ROOT)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
