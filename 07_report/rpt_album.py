# -*- coding: utf-8 -*-
"""
07_report / rpt_album.py
------------------------
The plots album: every image of the case study and of results/figures, one or
two to a page, each with a caption that says what it shows. Images of the
reconstruction carry the standing note on what it is.

Author: Akosa Samuel Onyejekwe (independent)

Units: lengths in points (1/72 inch).
"""
import re
from pathlib import Path

from reportlab.lib.units import cm
from reportlab.platypus import KeepTogether, PageBreak, Paragraph, Spacer

import project_meta as meta
from rpt_blocks import escape
from rpt_data import sibling

from unistall import flowfield

pdf = sibling("rpt_pdf")

TALL = 0.78                  # an image taller than this multiple of its width has a page to itself
CAPTION_ROOM = 4.6*cm        # height kept for a caption and the heading of a group
RECON = "Standing note: " + flowfield.STANDING_NOTE
MODEL_OUTPUT = "The cycle damping is model output only: the model does not predict cycle damping through stall."
OVERLAY = "06_postprocessing/overlay_record.csv"
CORRECTED = ("Cp", "Mach", "Tstatic", "Trecovery")
QUANTITY = {
    "Cp": "Pressure coefficient (Kármán-Tsien corrected)",
    "Mach": "Local Mach number, from the corrected pressure by the isentropic relations",
    "Tstatic": "Static temperature, an isentropic estimate from the corrected pressure",
    "Trecovery": "Recovery temperature (that of an adiabatic wall under a turbulent boundary layer, which the "
                 "reconstruction does not have; away from the wall, what such a wall would take at the local state)",
    "speed_stream": "Speed of the incompressible solution and streamlines",
    "vorticity": "Vorticity (zero outside the core of the vortex marker)", "vectors": "Velocity vectors",
}
MASKED = "06_postprocessing/masked_region_record.csv"
LOADS = {"cl": "lift", "cm": "moment about the quarter chord", "cd": "drag"}
SINGLE = {
    "static_inputs": "Static inputs of the load model at the Mach numbers of the cases: normal force, "
                     "separation point and moment against incidence.",
    "convergence_residuals": "Largest change of the loads over the cycle when one more cycle is marched, against "
                             "the number of cycles marched.",
    "timestep_refinement": "Refinement of the time step: difference of the peak loads from the finest march, "
                           "against the number of steps per cycle.",
    "fig3d_response_surface": "Peak lift and cycle damping given by the load model over mean incidence and "
                              "reduced frequency, at the Mach number and amplitude of Case B. " + MODEL_OUTPUT
                              + " The figure prints how often the model gets the sign of the damping right on the "
                              "held-out loops.",
    "fig3d_cp_phase_surface": "Reconstructed upper-surface pressure of Case A against chordwise position and cycle "
                              "phase as a surface, suction upwards, cut off at the critical value. " + RECON,
    "fig3d_field_surface_Cp": "The corrected pressure coefficient of Case A at the instant of greatest lift, drawn "
                              "as a surface over the plane, suction upwards, cut off at the critical value. " + RECON,
    "fig3d_field_surface_speed": "The speed of the incompressible solution of Case A at the instant of greatest "
                                 "lift, drawn as a surface over the plane. " + RECON,
    "fig3d_section_vectors": "The section drawn with span, with the reconstructed velocity on its near end plane. "
                             "The span is for the picture only; the calculation is two-dimensional. " + RECON,
    "fig_geometry_profile": "The NACA 0012 section, from the four-digit thickness form with a closed trailing edge.",
    "fig_geometry_thickness": "Thickness distribution of the NACA 0012 section.",
    "fig_mesh_full": "The O-grid round the section. An illustration: no equation of the flow is solved on it and "
                     "the load model does not read it.",
    "fig_mesh_le_zoom": "The O-grid at the leading edge. Not read by the load model.",
    "fig_mesh_te_zoom": "The O-grid at the trailing edge. Not read by the load model.",
    "fig_mesh_wall_spacing": "Spacing of the O-grid along the wall and away from it. Not read by the load model.",
}
VALIDATION = {
    "fig_validation_static": "Static data at Mach 0.30 beside what the model reads: its static normal force (set "
                             "beside the measured lift), its static moment, and its drag in a slow sweep. The same "
                             "data as results/ and data/; not a second assessment.",
    "fig_validation_loops": "Measured against predicted loops of lift, moment and drag for loops chosen by rule from "
                            "results/validation_table.csv: lowest, median and highest loop error in lift of the "
                            "held-out loops, and the calibration loop of Case A. The same numbers as results/; not a "
                            "second assessment.",
}
PER_CASE = {
    "timehist_loads": "lift, drag and moment through the reported cycle.",
    "states": "states of the load model through the reported cycle: static and lagged separation point, vortex "
              "normal force, lagged normal force over its onset level, and the vortex clock.",
    "cp_distribution": "reconstructed surface pressure at instants of the cycle, suction upwards, the band beyond "
                       "the critical value hatched. " + RECON,
    "cp_phase_map": "reconstructed upper-surface pressure over the chord and the cycle; " + flowfield.CRITICAL_NOTE
                    + " " + RECON,
    "temperature_profile": "static and recovery temperature along the surface at instants of the cycle, from the "
                           "corrected pressure by the isentropic relations; the curves stop where that pressure is "
                           "below the critical value. " + RECON,
    "fig_mesh_field": "pressure coefficient of the reconstructed field evaluated at the nodes of the O-grid, at "
                      "the instant of peak lift. The grid is only the set of points; no equation of the flow is "
                      "solved on it. " + RECON,
}


def case_of(stem: str) -> tuple:
    """(case key, its label) for a file name that carries a case key."""
    for key, case in meta.CASES.items():
        if key in stem:
            return key, case["label"]
    return "", ""


def contour_caption(stem: str, D) -> str:
    """Caption of a field image named contour_<quantity>_<case>_<phase>_a<deg>."""
    key, label = case_of(stem)
    m = re.match(rf"contour_(\w+)_{key}_(\w+)_a(\d+)$", stem)
    quantity, phase, deg = m.group(1), m.group(2), m.group(3)
    when = D.json("03_model_setup/solver_config.json")["field_reconstruction"]["field_phases"][phase]
    masked = ""
    if quantity in CORRECTED:
        row = D.csv(MASKED).set_index("file").loc[stem + ".png"]
        masked = (f" {flowfield.CRITICAL_NOTE} At this instant that region is {100.0*row['surface_fraction']:.1f} % "
                  f"of the surface and {100.0*row['field_fraction']:.3f} % of the stored field.")
    return (f"{QUANTITY[quantity]} of the reconstructed field. {label}; {when} (incidence {deg}° to the "
            f"nearest degree, from the file name).{masked} The separated region is not drawn; the model's "
            f"separation point is marked. {flowfield.AXES_NOTE} {RECON}")


def plot_caption(stem: str, D) -> str:
    """Caption of an image of 06_postprocessing/plots, 01_geometry or 02_mesh."""
    if stem in SINGLE:
        return SINGLE[stem]
    key, label = case_of(stem)
    if stem.startswith("contour_"):
        return contour_caption(stem, D)
    if stem.startswith("hyst_"):
        return hysteresis_caption(stem, key, label, D)
    for prefix, text in PER_CASE.items():
        if stem.startswith(prefix):
            return f"{label}: {text}"
    print(f"[album] no caption written for '{stem}': captioned by its file name")
    return f"{stem.replace('_', ' ')} (captioned by its file name; the album has no description of this image)."


def hysteresis_caption(stem: str, key: str, label: str, D) -> str:
    """Caption of a load loop, from the record of whether the measured loop
    was drawn (06_postprocessing/overlay_record.csv)."""
    load = LOADS[stem.split("_")[1]]
    rec = D.csv(OVERLAY).set_index("case").loc[key]
    frame = meta.CASES[key]["measured_frame"]
    if not frame:
        return f"{label}: loop of {load} against incidence. No measurement is compared with."
    if not bool(rec["drawn"]):
        return (f"{label}: loop of {load} against incidence. The measured points of {frame.replace('_', ' ')} are "
                f"NOT drawn: {rec['note']}.")
    return (f"{label}: loop of {load} against incidence at the nominal condition of the case, with the measured "
            f"points of {frame.replace('_', ' ')} of NASA TM-84245, a calibration loop, taken at its own "
            f"conditions (given in the legend).")


def sheet_caption(stem: str) -> str:
    """Caption of a drawing sheet."""
    number, words = stem.split("_", 1)
    return (f"Drawing {number.replace('sheet', 'sheet ')}: {words.replace('_', ' ').replace('AA', 'A-A').replace('3view', 'three views')}. "
            "Illustrative; every dimension from the table with its sources; not for manufacture.")


def result_captions(D) -> dict:
    """Captions of results/figures, from the figure list of README.md."""
    out = {}
    for table in D.md("README.md").values():
        if "Figure" not in table.columns or "Shows" not in table.columns:
            continue
        for link, shows in table.itertuples(index=False):
            m = re.search(r"\(results/figures/(\w+)\.png\)", link)
            if m:
                out[m.group(1)] = shows.rstrip(".") + "."
    return out


def groups(D) -> list:
    """(heading, [(path, caption)]) for every group of images, in order."""
    root = D.root
    res = result_captions(D)
    out = [(heading, [(p, plot_caption(p.stem, D)) for p in sorted((root/folder).glob("*.png"))])
           for heading, folder in (("Geometry (01_geometry)", "01_geometry"), ("The grid (02_mesh)", "02_mesh"))]
    plots = sorted((root/"06_postprocessing/plots").glob("*.png"))
    fields = [p for p in plots if p.stem.startswith("contour_")]
    out.append(("Loads, states and convergence (06_postprocessing)",
                [(p, plot_caption(p.stem, D)) for p in plots if p not in fields]))
    for key, case in meta.CASES.items():
        out.append((f"Reconstructed fields, {case['label']} (06_postprocessing)",
                    [(p, plot_caption(p.stem, D)) for p in fields if key in p.stem]))
    out.append(("Comparison with measurement (06_postprocessing/validation)",
                [(p, VALIDATION[p.stem]) for p in sorted((root/"06_postprocessing/validation").glob("*.png"))]))
    out.append(("Drawings (08_engineering_drawings)",
                [(p, sheet_caption(p.stem)) for p in sorted((root/"08_engineering_drawings").glob("*.png"))]))
    figs = sorted((root/"results/figures").glob("*.png"))
    missing = [p.name for p in figs if p.stem not in res]
    if missing:
        raise KeyError(f"README.md lists no caption for {missing}")
    out.append(("Calibration and accuracy of the load model (results/figures)",
                [(p, res[p.stem]) for p in figs]))
    return out


def is_tall(pics, path: Path) -> bool:
    """True for an image that is given a page to itself."""
    w, h = pics.size(path)
    return h/w > TALL


def image_block(pics, path: Path, caption: str, st: dict, doc, root: Path) -> KeepTogether:
    """An image with its caption, sized so that two wide images share a page
    and a tall one has a page to itself."""
    half = doc.avail_h/2 - CAPTION_ROOM
    box_h = doc.avail_h - CAPTION_ROOM if is_tall(pics, path) else half
    img = pdf.picture(pics, path, doc.avail_w, box_h)
    img.hAlign = "CENTER"
    name = f"<b>{escape(str(path.relative_to(root)))}</b>"
    return KeepTogether([img, Paragraph(f"{name}<br/>{escape(caption)}", st["caption"]), Spacer(1, 4)])


def build(D, pics, path: Path) -> tuple:
    """Write the album; returns (number of images, rows of its page record).
    Within a group the images that take a page to themselves come first and
    the ones that share a page after them, so that only the last page of a
    group can be left part empty."""
    st = pdf.styles()
    doc = pdf.make_doc(path, f"{meta.SOLVER} plots album", f"{meta.SOLVER} plots album")
    sets = groups(D)
    count = sum(len(items) for _, items in sets)
    D.note(count, "01_geometry, 02_mesh, 06_postprocessing/plots, 06_postprocessing/validation, "
           "08_engineering_drawings, results/figures",
           "number of *.png files", "{:d}")
    story = pdf.cover(st, "Plots album", [
        f"{count} images: every image of the case study and of results/figures.",
        "The flow-field images show a reconstruction, not a flow solution: " + flowfield.STANDING_NOTE,
        "Case B is an aerofoil at a rotor-blade-station condition in a steady stream with prescribed pitch, not "
        "a rotor calculation."], doc.avail_h)
    for heading, items in sets:
        ordered = sorted(items, key=lambda item: not is_tall(pics, item[0]))
        story.append(Paragraph(escape(heading), st["h1"]))
        story += [image_block(pics, p, cap, st, doc, D.root) for p, cap in ordered]
        story.append(PageBreak())
    doc.build(story[:-1])
    return count, pdf.short_pages(doc, path.name)
