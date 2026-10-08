# -*- coding: utf-8 -*-
# Run from the repository root:  PYTHONPATH=. python3 07_report/build_report.py
"""
07_report / build_report.py
---------------------------
Builds the documents of the case study from the files of the other folders:

  case.docx                    the case-study report as a Word document
  UNISTALL_report.pdf          the same report as a PDF
  UNISTALL_plots_album.pdf     every image of the case study and of results/figures
  UNISTALL_data_dossier.pdf    the data behind the report as tables
  _equations/eq_NNN.png        the equations of the report
  report_numbers.json          every number quoted in the report text, with the
                               file and the field it was read from
  figure_lettering.csv         each figure of the report: the width it is
                               placed at and the size of its smallest
                               lettering as printed
  page_fill.csv                how full each page of the three PDFs is

Author: Akosa Samuel Onyejekwe (independent)

Nothing here computes a result of the study: every number is read from a file
at build time. The build stops if the reduced-time step the load model
marches differs from the step the configuration gives
(rpt_text_mid.check_step).
Units: as in the files read; lettering in points.
"""
import importlib.util
import json
import sys
from pathlib import Path

import pandas as pd
from PIL import Image

import project_meta as meta

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
FIGURE_RECORDS = ("01_geometry", "02_mesh", "06_postprocessing", "06_postprocessing/validation",
                  "08_engineering_drawings")
REPORT_DPI = 190             # pixels per inch of page given to a figure of the report
ALBUM_DPI = 140              # the same for the plots album


def helper(name: str):
    """Load a helper module that lies beside this script."""
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, HERE/f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def write_numbers(D, equations: dict, counts: dict) -> None:
    """Write report_numbers.json."""
    eqs = [dict(number=e["number"], key=key, image=f"_equations/{e['path'].name}",
                formulation_lines=e["source"], formulation_text=e["source_text"])
           for key, e in equations.items()]
    record = dict(
        author=meta.AUTHOR_FULL,
        note="Every number quoted in the text of the report, in the order it is quoted, with the file and the "
             "field it was read from. Tables shown whole are listed under tables. Equations are listed with "
             "the numbered lines of docs/formulation.md they are taken from.",
        documents=counts, numbers=D.quoted, tables=D.tables_used, equations=eqs)
    (HERE/"report_numbers.json").write_text(json.dumps(record, indent=1, ensure_ascii=False, default=str) + "\n",
                                            encoding="utf-8")


def figure_lettering(placed: list) -> pd.DataFrame:
    """For each image placed in the report: the width it is placed at and
    the size of its smallest lettering as printed, from the figure record
    its script wrote (file, dots per inch, smallest lettering). An image
    whose script keeps no record (those of results/figures) is listed with
    no lettering size."""
    records = {}
    for folder in FIGURE_RECORDS:
        rec = ROOT/folder/"figure_record.csv"
        if rec.exists():
            for row in pd.read_csv(rec).itertuples():
                records[row.file] = (float(row.dpi), float(row.smallest_lettering_pt))
    rows = []
    for item in placed:
        path = Path(item["path"])
        with Image.open(path) as im:
            width_px = im.size[0]
        dpi, small = records.get(path.name, (float("nan"), float("nan")))
        native_in = width_px/dpi
        rows.append(dict(file=str(path.relative_to(ROOT)), placed_width_in=round(item["width_in"], 3),
                         native_width_in=round(native_in, 3), smallest_lettering_pt=small,
                         smallest_lettering_as_printed_pt=round(small*item["width_in"]/native_in, 2),
                         recorded=path.name in records))
    return pd.DataFrame(rows)


def main() -> None:
    """Build every document."""
    data = helper("rpt_data")
    blocks = helper("rpt_blocks")
    equations = helper("rpt_equations")
    D = data.Data()
    eqs = equations.render_all(HERE/"_equations", data.formulation_equations(data.ROOT/"docs/formulation.md"))
    R = blocks.Report(eqs, data.bib_entries(data.ROOT/"docs/references.bib"))
    for part in ("rpt_text_front", "rpt_text_mid", "rpt_text_back"):
        helper(part).build(R, D)
    R.resolve()
    pics = blocks.Pictures(REPORT_DPI)
    counts = dict(report_blocks=len(R.blocks), equations=len(eqs), references=len(R.cited))
    helper("rpt_docx").render_report(R, pics, HERE/"case.docx")
    pages, placed = helper("rpt_pdf").render_report(R, pics, HERE/f"{meta.SOLVER}_report.pdf")
    figure_lettering(placed).to_csv(HERE/"figure_lettering.csv", index=False)
    album_pics = blocks.Pictures(ALBUM_DPI)
    counts["album_images"], album_pages = helper("rpt_album").build(D, album_pics, HERE/f"{meta.SOLVER}_plots_album.pdf")
    counts["report_images_as"], counts["album_images_as"] = pics.used, album_pics.used
    counts["dossier_tables"], dossier_pages = helper("rpt_dossier").build(D, HERE/f"{meta.SOLVER}_data_dossier.pdf")
    pd.DataFrame(pages + album_pages + dossier_pages).to_csv(HERE/"page_fill.csv", index=False)
    counts["report_numbers_quoted"] = len(D.quoted)             # counted when every document has been built
    write_numbers(D, eqs, counts)
    for name in ("case.docx", f"{meta.SOLVER}_report.pdf", f"{meta.SOLVER}_plots_album.pdf",
                 f"{meta.SOLVER}_data_dossier.pdf", "report_numbers.json"):
        print(f"wrote 07_report/{name}  {(HERE/name).stat().st_size/1e6:.2f} MB")
    print(counts)


if __name__ == "__main__":
    main()
