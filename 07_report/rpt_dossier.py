# -*- coding: utf-8 -*-
"""
07_report / rpt_dossier.py
--------------------------
The data dossier: the data behind the report as tables, each with the file it
is read from, and an inventory of 05_solution with the size of every file and
what it holds.

Author: Akosa Samuel Onyejekwe (independent)

Units: as in the files shown; file sizes in kilobytes of 1000 bytes.
"""
import re
from pathlib import Path

import pandas as pd
from reportlab.lib.pagesizes import A4, landscape
from reportlab.platypus import CondPageBreak, Paragraph, Spacer

import project_meta as meta
from rpt_blocks import COMPUTER_NOTATION, escape
from rpt_data import cell, frame_rows, sibling

pdf = sibling("rpt_pdf")

DIGITS = 10                  # significant figures kept, so that values read as they do in the files
POLAR_EVERY = 4              # one row in this many of the model's static polar is shown
STATIC_EVERY = 8             # one row in this many of the static inputs is shown
HOLDS = (
    ("field_.*_on_mesh", "reconstructed pressure coefficient and speed at the nodes of the O-grid of 02_mesh at the "
                         "instant of peak lift, with the source of each value (reconstruction, not a flow "
                         "solution), gzip-compressed"),
    ("reconstruction_timing", "CPU time of the reconstruction of each case"),
    ("reconstruction_", "the reconstruction instant by instant: transpiration, circulations, the vortex marker and "
                        "the pressure beneath it, extremes of the corrected pressure, of the local Mach number and "
                        "of the temperatures, the size of the region beyond the critical pressure, closure of lift "
                        "(corrected and incompressible), moment and drag"),
    ("cycle_scan_", "the transpiration, the closure of the lift and the vortex marker at instants spread over the "
                    "cycle"),
    ("time_history_", "loads and model states at every step of the reported cycle"),
    ("cp_distribution_", "reconstructed surface pressure (corrected and incompressible), local Mach number, static "
                         "and recovery temperature, speed and transpiration at instants of the cycle, with the "
                         "flag of the region beyond the critical pressure (reconstruction, not a flow solution)"),
    ("field_", "reconstructed field on the rectangular grid at one instant: incompressible velocity, pressure "
               "coefficient (corrected and incompressible), local Mach number, static and recovery temperature, "
               "vorticity, the flag of the region beyond the critical pressure (values not physical there) and a "
               "flag for the nodes inside the section (reconstruction, not a flow solution), gzip-compressed"),
    ("metrics_", "scalar results of the case"),
    ("summary_all_cases", "one row of main results per case"),
    ("held_out_error_bands", "error of the load model in each headline quantity on the held-out loops"),
    ("static_station_spread", "Case B run with the static inputs of each Mach station in place of their blend"),
    ("model_static_polar_", "lift, drag, moment and separation point of the model in a slow sweep at the case's "
                            "Mach number: its quasi-steady limit, not a measurement"),
    ("response_surface", "peak lift, minimum moment, damping and onset incidence over mean incidence and "
                         "reduced frequency"),
    ("runtime_environment", "machine, library versions and CPU time of the march"),
    ("residuals_", "the reported loads and their largest change over the cycle against the number of cycles "
                   "marched"),
    ("timestep_refinement", "loads of each case against the number of steps per cycle"),
    ("timestep_order", "observed order of convergence, Richardson estimate and discretisation uncertainty of each "
                       "load"),
    ("panel_convergence_summary", "largest change of the quoted quantities of the reconstruction when the panels "
                                  "are doubled, against the stated tolerances"),
    ("panel_convergence_", "the reconstruction at the six instants against the number of panels"),
)


def flatten(node, prefix: str = "") -> list:
    """A nested dict as rows of (dotted key, value as text). A list of
    records is flattened record by record, each under its position in the
    list; a list of plain values is one row."""
    if isinstance(node, dict):
        rows = []
        for key, value in node.items():
            rows += flatten(value, f"{prefix}.{key}" if prefix else str(key))
        return rows
    if isinstance(node, list):
        if any(isinstance(v, (dict, list)) for v in node):
            rows = []
            for i, value in enumerate(node, start=1):
                rows += flatten(value, f"{prefix}.{i}")
            return rows
        return [[prefix, ", ".join(cell(v, DIGITS) for v in node)]]
    return [[prefix, cell(node, DIGITS)]]


class Dossier:
    """Collects the sections of the dossier."""

    def __init__(self, D, doc, st: dict) -> None:
        self.D, self.doc, self.st = D, doc, st
        self.story = []
        self.count = 0

    def add(self, title: str, source: str, header: list, rows: list, size: float = 8.0) -> None:
        """One table with its title and source."""
        self.count += 1
        for row in rows:
            for c in row:
                if "{" in str(c) or COMPUTER_NOTATION.search(str(c)):
                    raise ValueError(f"dossier table '{title}': cell '{c}' is a raw record or in computer notation")
        self.D.whole(source, f"dossier table '{title}'")
        self.story += [CondPageBreak(pdf.TABLE_START + 30), Paragraph(escape(f"{self.count}  {title}"), self.st["h2"]),
                       Paragraph(escape(f"Source: {source}. {len(rows)} rows."), self.st["table_caption"]),
                       pdf.table_flowable(header, rows, self.doc.avail_w, size), Spacer(1, 14)]

    def csv(self, title: str, rel: str, size: float = 8.0, df: pd.DataFrame | None = None) -> None:
        """A CSV file, or a part of it handed in as `df`, as a table."""
        header, rows = frame_rows(self.D.csv(rel) if df is None else df, DIGITS)
        self.add(title, rel, header, rows, size)

    def md(self, title: str, rel: str, heading: str) -> None:
        """A table of a Markdown file."""
        table = self.D.md(rel)[heading]
        self.add(title, f"{rel}, {heading}", list(table.columns), [list(r) for r in table.itertuples(index=False)])


def setup_tables(X: Dossier) -> None:
    """Conditions, kinematics, air, configuration, constants, static inputs."""
    D = X.D
    X.csv("Flow conditions", "03_model_setup/flow_conditions.csv")
    X.csv("Kinematics", "03_model_setup/kinematics.csv")
    X.csv("How the condition of Case B follows from the rotor (Case B is not a rotor calculation)",
          "03_model_setup/station_condition.csv")
    X.csv("Air properties", "03_model_setup/air_properties.csv")
    X.add("Solver configuration", "03_model_setup/solver_config.json", ["Key", "Value"],
          flatten(D.json("03_model_setup/solver_config.json")))
    cal = {k: v for k, v in D.json("results/calibrated_constants.json").items() if k != "first_stage_record"}
    X.add("Calibrated constants and the record of the fit", "results/calibrated_constants.json (without the copy "
          "kept under first_stage_record)", ["Key", "Value"], flatten(cal))
    X.add("Unsteady-moment factor (empirical; a correction of this model)", "results/attached_moment_factor.json",
          ["Key", "Value"], flatten(D.json("results/attached_moment_factor.json")))
    X.md("Constants, literature and calibrated", "results/tables.md", "Table 3. Constants")
    X.csv("Static inputs at the Mach stations", "results/static_by_mach_naca0012.csv", 7.5)
    static = D.csv("03_model_setup/static_inputs.csv")
    sampled = pd.concat([part.iloc[::STATIC_EVERY] for _, part in static.groupby("case", sort=False)])
    X.csv(f"Static inputs of the cases, sampled (one row in {STATIC_EVERY} of {len(static)})",
          "03_model_setup/static_inputs.csv", df=sampled)
    X.csv("Section geometry", "01_geometry/section_geometry_summary.csv")
    X.csv("Mesh quality (the load model does not read the grid)", "02_mesh/mesh_quality_metrics.csv")


def solution_tables(X: Dossier) -> None:
    """Metrics, convergence and the response surface."""
    D = X.D
    names = list(meta.CASES)
    merged = D.csv(f"05_solution/metrics_{names[0]}.csv").rename(columns={"value": names[0]})
    for name in names[1:]:
        merged = merged.merge(D.csv(f"05_solution/metrics_{name}.csv").rename(columns={"value": name}),
                              on=["metric", "units"])
    X.csv("Metrics of the cases (cycle damping is model output only)",
          " and ".join(f"05_solution/metrics_{n}.csv" for n in names), df=merged[["metric", *names, "units"]])
    summary = D.csv("05_solution/summary_all_cases.csv").set_index("case").T.reset_index(names="quantity")
    X.csv("Summary of the cases", "05_solution/summary_all_cases.csv", df=summary)
    for name in names:
        X.csv(f"Convergence with cycles marched, {meta.CASES[name]['label']}",
              f"05_solution/convergence/residuals_{name}.csv")
    X.csv("Convergence with the time step", "05_solution/convergence/timestep_refinement.csv", 7.0)
    X.csv("Observed order of convergence and discretisation uncertainty", "05_solution/convergence/timestep_order.csv",
          6.5)
    X.csv("The reconstruction against the number of panels: largest change on doubling",
          "05_solution/convergence/panel_convergence_summary.csv", 7.5)
    for name in names:
        X.csv(f"The reconstruction against the number of panels, {meta.CASES[name]['label']}",
              f"05_solution/convergence/panel_convergence_{name}.csv", 6.5)
    X.csv("CPU time of the reconstruction", "05_solution/convergence/reconstruction_timing.csv")
    X.csv("Error of the load model in each headline quantity on the held-out loops",
          "05_solution/held_out_error_bands.csv", 7.0)
    X.csv("Case B with the static inputs of each Mach station", "05_solution/static_station_spread.csv", 7.5)
    for name in names:
        polar = D.csv(f"05_solution/model_static_polar_{name}.csv")
        X.csv(f"Static polar of the model in a slow sweep, {meta.CASES[name]['label']}, sampled (one row in "
              f"{POLAR_EVERY} of {len(polar)}); the model's quasi-steady limit, not a measurement",
              f"05_solution/model_static_polar_{name}.csv", 7.5, df=polar.iloc[::POLAR_EVERY].drop(columns="note"))
    X.csv("Convergence of the calibration loops", "results/convergence_summary.csv", 7.5)
    X.csv("Response surface (a section in a steady stream with prescribed pitch; damping is model output only)",
          "05_solution/response_surface.csv", 7.5)
    for name in names:
        X.csv(f"The reconstruction instant by instant, {meta.CASES[name]['label']} (not a flow solution)",
              f"05_solution/reconstruction_{name}.csv", 6.5)
    X.csv("State-space form against the indicial march", "results/statespace_summary.csv", 7.0,
          df=D.csv("results/statespace_summary.csv").T.reset_index(names="quantity").rename(columns={0: "value"}))
    X.csv("A section free in torsion: the two forms of the model (model output only)",
          "results/structural_checks.csv", 6.5)
    X.csv("Calibration cost against the stall-onset level", "results/onset_level_check.csv")
    X.csv("Runtime environment", "05_solution/runtime_environment.csv")


def accuracy_tables(X: Dossier) -> None:
    """Held-out accuracy, the scoreboard and the dimensions."""
    X.md("Accuracy on the held-out loops", "results/tables.md", "Table 5. Accuracy on the held-out loops")
    X.csv("Held-out measures against their targets", "results/validation_targets.csv")
    X.csv("Held-out accuracy by stall depth", "results/validation_by_stall_depth.csv")
    X.csv("Cycle damping on the held-out loops (the model does not predict cycle damping through stall)",
          "results/validation_damping.csv", 7.0)
    X.csv("Tabulated against fitted separation law, same constants", "results/comparison_same_constants.csv", 7.5)
    X.csv("Every target against its measured value", "results/targets_scoreboard.csv", 7.5)
    X.csv("Scored loops, one row each (the same numbers as results/validation_table.csv)",
          "06_postprocessing/validation/validation_loops_summary.csv", 6.0)
    X.csv("Constants of the calibration, fitted and fixed", "06_postprocessing/validation/calibration_constants.csv")
    X.csv("Dimensions of the drawings (illustrative)", "08_engineering_drawings/dimensions.csv")
    X.csv("What each drawing sheet states", "08_engineering_drawings/sheet_record.csv")
    X.csv("Whether the measured loop was drawn with each case", "06_postprocessing/overlay_record.csv")
    X.csv("The region beyond the critical pressure printed on each figure of a corrected quantity",
          "06_postprocessing/masked_region_record.csv", 7.5)


def holds(name: str) -> str:
    """What a file of 05_solution holds."""
    for prefix, text in HOLDS:
        if re.match(prefix, name):
            return text
    return "not described in the dossier"


def inventory(X: Dossier) -> None:
    """Every file of 05_solution with its size and contents."""
    folder = X.D.root/"05_solution"
    rows = []
    for path in sorted(p for p in folder.rglob("*") if p.is_file()):
        shape = pd.read_csv(path).shape if path.name.endswith((".csv", ".csv.gz")) else ("", "")
        rows.append([str(path.relative_to(folder)), f"{path.stat().st_size/1000.0:.1f}", cell(shape[0]),
                     cell(shape[1]), holds(path.name)])
    X.add("File inventory of 05_solution", "05_solution (sizes read from the files at build time)",
          ["File", "Size, kB", "Rows", "Columns", "Holds"], rows)


def build(D, path: Path) -> tuple:
    """Write the dossier; returns (number of tables in it, rows of its page record)."""
    st = pdf.styles()
    doc = pdf.make_doc(path, f"{meta.SOLVER} data dossier", f"{meta.SOLVER} data dossier", landscape(A4))
    X = Dossier(D, doc, st)
    setup_tables(X)
    solution_tables(X)
    accuracy_tables(X)
    inventory(X)
    story = pdf.cover(st, "Data dossier", [
        "The data behind the report, as tables. Each table names the file it is read from; no number is typed "
        "by hand.",
        "The loads are the prediction of the load model. The flow fields listed in the inventory are a "
        "potential-flow reconstruction, incompressible with a linearised correction of the pressure, "
        "quasi-steady, not a flow solution.",
        "Case B is an aerofoil at a rotor-blade-station condition in a steady stream with prescribed pitch, not "
        "a rotor calculation."], doc.avail_h)
    doc.build(story + X.story)
    return X.count, pdf.short_pages(doc, path.name)
