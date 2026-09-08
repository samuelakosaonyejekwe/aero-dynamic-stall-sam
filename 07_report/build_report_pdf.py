# -*- coding: utf-8 -*-
"""
07_report / build_report_pdf.py
-------------------------------
Builds the consolidated deliverable  aero_dynamic_stall_report.pdf  as

    case.docx (rendered)  +  UNISTALL_data_dossier.pdf  +  UNISTALL_plots_album.pdf

case.docx is rendered here directly rather than exported by hand from Word.
The manual export was the reason the consolidated report drifted out of step
with the solver output; making it a pipeline stage removes that failure mode.

Author: Akosa Samuel Onyejekwe (independent).  No black is used anywhere.
"""
import html
import sys
from pathlib import Path

# This stage needs more than the solver pipeline does. Fail with an actionable
# message rather than a bare ImportError three frames deep.
_MISSING = []
for _mod, _pkg in (("fitz", "PyMuPDF"), ("docx", "python-docx"),
                   ("PIL", "pillow"), ("reportlab", "reportlab")):
    try:
        __import__(_mod)
    except ImportError:
        _MISSING.append(_pkg)
if _MISSING:
    sys.exit(
        "[report] missing report-build dependencies: " + ", ".join(_MISSING) + "\n"
        "         install them with:  pip install " + " ".join(_MISSING) + "\n"
        "         (on a PEP 668 / externally-managed Python, either use a venv --\n"
        "          python3 -m venv .venv && .venv/bin/pip install -r requirements.txt\n"
        "          -- or add --user to install into your own site-packages.)\n"
        "         The solver pipeline itself does not need any of these."
    )

import re
import fitz                                   # PyMuPDF, for the final merge
from docx import Document
from docx.oxml.ns import qn
from docx.enum.text import WD_ALIGN_PARAGRAPH
from PIL import Image as PILImage


from reportlab.lib.units import inch
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.enums import TA_LEFT, TA_CENTER, TA_RIGHT, TA_JUSTIFY
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (BaseDocTemplate, PageTemplate, Frame, Paragraph,
                                Spacer, Image, Table, TableStyle, PageBreak,
                                KeepTogether)

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from project_meta import AUTHOR, TITLE, SOLVER
HERE = ROOT/"07_report"
DOCX = HERE/"case.docx"
OUT  = ROOT/"aero_dynamic_stall_report.pdf"

INK      = colors.HexColor("#1f3350")         # deep navy — never black
INK_SOFT = colors.HexColor("#4a5d78")
ACC      = colors.HexColor("#2e6fb7")
HDR_BG   = colors.HexColor("#dfe8f2")
ROW_BG   = colors.HexColor("#f3f6fa")

# ---------------------------------------------------------------- fonts
# DejaVu Sans Condensed: the family with the full Unicode coverage the document
# needs (Greek alpha/beta/Xi/psi, plus the arrow, degree, +/- and (TM) glyphs).
# The Condensed cut specifically. This is not cosmetic: measured with
# pdfmetrics.stringWidth on a representative line of this report's body text,
# the regular cut sets 11.2% wider (1020.4 pt vs 917.9 pt at 10 pt), which is
# 91 rather than 101 characters on a 6.5 in line. Substituting it would
# repaginate the whole 73-page body, so it is not a drop-in replacement.
#
# The install location differs per distribution, so search rather than assume
# one. This used to hardcode the Debian path, which was harmless while this
# script was not distributed, but would fail on Fedora, Arch or macOS.
_FACES = {"Body": "DejaVuSansCondensed.ttf",
          "Body-Bold": "DejaVuSansCondensed-Bold.ttf",
          "Body-Italic": "DejaVuSansCondensed-Oblique.ttf",
          "Body-BoldItalic": "DejaVuSansCondensed-BoldOblique.ttf"}
_FONT_ROOTS = [Path("/usr/share/fonts"), Path("/usr/local/share/fonts"),
               Path("/opt/homebrew/share/fonts"), Path("/Library/Fonts"),
               Path.home()/".local"/"share"/"fonts", Path.home()/".fonts"]

def _find_face(fn):
    """First match for one face file under any known font root."""
    for root in _FONT_ROOTS:
        if not root.is_dir():
            continue
        direct = root/"truetype"/"dejavu"/fn          # the common case, no walk
        if direct.exists():
            return direct
        try:
            hit = next(root.rglob(fn), None)
        except (OSError, PermissionError):
            hit = None
        if hit is not None:
            return hit
    return None

_resolved, _missing_faces = {}, []
for _name, _fn in _FACES.items():
    _hit = _find_face(_fn)
    if _hit is None:
        _missing_faces.append(_fn)
    else:
        _resolved[_name] = _hit
if _missing_faces:
    sys.exit(
        "[report] cannot find these DejaVu Sans Condensed font files: "
        + ", ".join(_missing_faces) + "\n"
        "         Searched: " + ", ".join(str(r) for r in _FONT_ROOTS) + "\n"
        "         The report needs a family with full Unicode coverage (Greek\n"
        "         alpha/beta/Xi/psi, plus the arrow, degree and (TM) glyphs).\n"
        "         Debian/Ubuntu:  sudo apt install fonts-dejavu-core\n"
        "         Fedora:         sudo dnf install dejavu-sans-fonts\n"
        "         Arch:           sudo pacman -S ttf-dejavu\n"
        "         macOS:          brew install --cask font-dejavu\n"
        "         Note that matplotlib bundles DejaVu Sans but NOT the Condensed\n"
        "         cut this report is paginated for, so it is not a substitute.\n"
        "         Then re-run this stage. The solver pipeline does not need it."
    )
for name, _path in _resolved.items():
    pdfmetrics.registerFont(TTFont(name, str(_path)))
pdfmetrics.registerFontFamily("Body", normal="Body", bold="Body-Bold",
                              italic="Body-Italic", boldItalic="Body-BoldItalic")

# ------------------------------------------------- page geometry (matches docx)
doc_src = Document(str(DOCX))
_sec = doc_src.sections[0]
PAGE_W, PAGE_H = _sec.page_width.inches*inch, _sec.page_height.inches*inch
ML, MR = _sec.left_margin.inches*inch, _sec.right_margin.inches*inch
MT, MB = _sec.top_margin.inches*inch, _sec.bottom_margin.inches*inch
FRAME_W = PAGE_W - ML - MR
FRAME_H = PAGE_H - MT - MB

_ALIGN = {WD_ALIGN_PARAGRAPH.CENTER: TA_CENTER,
          WD_ALIGN_PARAGRAPH.RIGHT: TA_RIGHT,
          WD_ALIGN_PARAGRAPH.JUSTIFY: TA_JUSTIFY}

BODY_PT = 11.0
def _style(name, size, align=TA_LEFT, colour=INK, bold=False, space_after=6,
           space_before=0, leading_mult=1.24, left_indent=0):
    # bulletFontName MUST be set. ReportLab defaults it to Helvetica, which is a
    # base-14 font and is NOT embedded, so every bullet glyph in the document was
    # left to the reader's font substitution -- 16 spans of visible text in a
    # report that otherwise deliberately uses one embedded family throughout.
    return ParagraphStyle(name, fontName="Body-Bold" if bold else "Body",
                          fontSize=size, leading=size*leading_mult,
                          textColor=colour, alignment=align,
                          spaceAfter=space_after, spaceBefore=space_before,
                          leftIndent=left_indent, allowWidows=0, allowOrphans=0,
                          bulletFontName="Body", bulletFontSize=size)

S_BULLET = _style("bullet", BODY_PT, left_indent=16, space_after=3)
S_H1 = _style("h1", 16.0, colour=ACC, bold=True, space_before=16, space_after=8)
S_H2 = _style("h2", 13.5, colour=ACC, bold=True, space_before=12, space_after=6)
S_H3 = _style("h3", 12.0, colour=ACC, bold=True, space_before=10, space_after=5)
_HEADING = {"Heading 1": S_H1, "Heading 2": S_H2, "Heading 3": S_H3}


# ---------------------------------------------------------------- run markup
def _runs_markup(par, default_pt=BODY_PT):
    """Inline reportlab markup for one docx paragraph, preserving each run's
    bold / italic / size / colour."""
    out = []
    for r in par.runs:
        txt = html.escape(r.text).replace("\n", "<br/>").replace("\t", "&nbsp;&nbsp;")
        if not txt:
            continue
        rgb = getattr(r.font.color, "rgb", None) if r.font.color is not None else None
        size = r.font.size.pt if r.font.size is not None else None
        attrs = ""
        if size and abs(size - default_pt) > 0.01:
            attrs += ' size="%.2f"' % size
        if rgb is not None:
            attrs += ' color="#%s"' % str(rgb)
        if attrs:
            txt = "<font%s>%s</font>" % (attrs, txt)
        if r.italic:
            txt = "<i>%s</i>" % txt
        if r.bold:
            txt = "<b>%s</b>" % txt
        out.append(txt)
    return "".join(out)


def _dominant_pt(par, fallback):
    for r in par.runs:
        if r.text.strip() and r.font.size is not None:
            return r.font.size.pt
    return fallback


# ---------------------------------------------------------------- images
_IMG_CACHE = {}
def _image_flowables(par):
    """Every inline picture in this paragraph, as reportlab Image flowables.
    Widths are clamped to the text column: case.docx asks for up to 6.6 in on a
    6.0 in column, which in the Word export ran the figure off the page edge."""
    flows = []
    for blip in par._p.findall(".//" + qn("a:blip")):
        rid = blip.get(qn("r:embed"))
        if rid is None:
            continue
        part = doc_src.part.related_parts[rid]
        path = HERE/"_pdfimg"/Path(part.partname).name
        if rid not in _IMG_CACHE:
            path.parent.mkdir(exist_ok=True)
            path.write_bytes(part.blob)
            _IMG_CACHE[rid] = path
        path = _IMG_CACHE[rid]
        # requested display width (EMU -> pt), clamped to the frame.
        # wp:extent lives on the wp:inline / wp:anchor ancestor, five levels
        # above a:blip -- walk up by tag rather than counting parents, or every
        # picture (equations included) silently renders at full column width.
        node = blip
        while node is not None and node.tag not in (qn("wp:inline"), qn("wp:anchor")):
            node = node.getparent()
        ext = node.find(qn("wp:extent")) if node is not None else None
        req_w = (int(ext.get("cx"))/914400.0)*inch if ext is not None else FRAME_W
        w = min(req_w, FRAME_W)
        with PILImage.open(path) as im:
            iw, ih = im.size
        h = w*ih/float(iw)
        if h > FRAME_H*0.92:                     # never taller than a page
            h = FRAME_H*0.92; w = h*iw/float(ih)
        img = Image(str(path), width=w, height=h)
        img.hAlign = "CENTER"
        flows.append(img)
    return flows


def _is_caption(par):
    """build_docx.py captions: italic, 9.5 pt, accent colour, centred."""
    if par.alignment != WD_ALIGN_PARAGRAPH.CENTER or not par.runs:
        return False
    r = par.runs[0]
    return bool(r.italic) and r.font.size is not None and abs(r.font.size.pt - 9.5) < 0.2


# ---------------------------------------------------------------- tables
PAD = 6.0            # left+right cell padding, pt

def _col_floors(text_rows, ncol, hdr_pt, body_pt):
    """Minimum width each column needs so no cell breaks mid-word: the widest
    unbreakable token in it. Paths break at "/" and identifiers at "_", so those
    separators end a token; a filename or a bare word does not break at all."""
    floors = []
    for j in range(ncol):
        tokmax = 0.0
        for i, row in enumerate(text_rows):
            cell = row[j] if j < len(row) else ""
            fn, pt = ("Body-Bold", hdr_pt) if i == 0 else ("Body", body_pt)
            for tok in re.split(r"[\s/]+", cell):
                if tok:
                    tokmax = max(tokmax, pdfmetrics.stringWidth(tok, fn, pt))
        floors.append(min(tokmax + PAD, 0.42*FRAME_W))
    return floors


def _col_widths(text_rows, ncol, hdr_pt, body_pt):
    """Column widths from content: every column must at least fit its longest
    unbreakable word, the rest of the frame is shared out in proportion to the
    natural (unwrapped) width. Equal widths broke long headers mid-word."""
    nat, floor = [], []
    for j in range(ncol):
        wmax = tokmax = 0.0
        for i, row in enumerate(text_rows):
            cell = row[j] if j < len(row) else ""
            fn, pt = ("Body-Bold", hdr_pt) if i == 0 else ("Body", body_pt)
            wmax = max(wmax, pdfmetrics.stringWidth(cell, fn, pt))
            # A path is breakable at its separators, a filename is not. Splitting
            # on "/" as well as whitespace stops a long folder path reserving a
            # floor it does not need and starving the neighbouring column: adding
            # 06_postprocessing/validation/experimental to the data inventory
            # squeezed the file column until the three longest filenames broke
            # mid-extension, which is data-lossless but reads as a typo.
            for tok in re.split(r"[\s/]+", cell):
                if tok:
                    tokmax = max(tokmax, pdfmetrics.stringWidth(tok, fn, pt))
        nat.append(wmax + PAD)
        floor.append(min(tokmax + PAD, 0.42*FRAME_W))
    if sum(nat) <= FRAME_W:                       # fits: share out the slack
        extra = FRAME_W - sum(nat); tot = sum(nat) or 1.0
        return [w + extra*w/tot for w in nat]
    slack = FRAME_W - sum(floor)
    if slack <= 0:                                # even the floors do not fit
        k = FRAME_W/sum(floor)
        return [w*k for w in floor]
    over = [max(n - f, 0.0) for n, f in zip(nat, floor)]
    tot = sum(over) or 1.0
    return [f + slack*o/tot for f, o in zip(floor, over)]


def _table_flowable(tbl):
    text_rows = []
    for row in tbl.rows:
        text_rows.append([" ".join(p.text for p in c.paragraphs if p.text)
                          for c in row.cells])
    if not text_rows:
        return None
    ncol = max(len(r) for r in text_rows)

    # Shrink the type until every column's unbreakable minimum fits. This used to
    # test "abs(sum(widths) - FRAME_W) < 1.0", which every branch of _col_widths
    # satisfies by construction -- all three normalise to fill the frame -- so the
    # condition was true on the first iteration and the font NEVER shrank. The
    # visible cost was the summary table in section 11, whose eight columns were
    # squeezed until every header broke mid-word: "CL_m/ax", "CM_m/in",
    # "aero_dampin/g_Xi", "flutter_ri/sk", "A_validati/on".
    hdr_pt, body_pt = 9.0, 8.5
    while body_pt > 5.5 and sum(_col_floors(text_rows, ncol, hdr_pt, body_pt)) > FRAME_W:
        hdr_pt -= 0.5; body_pt -= 0.5
    widths = _col_widths(text_rows, ncol, hdr_pt, body_pt)
    s_h = _style("cellh%.1f" % hdr_pt, hdr_pt, bold=True, space_after=0, leading_mult=1.16)
    s_b = _style("cellb%.1f" % body_pt, body_pt, space_after=0, leading_mult=1.16)

    def _fit(txt, j, fn, pt):
        """Break an over-wide identifier at its own separators rather than let
        reportlab break it mid-word. A 16-column time-history table cannot give
        every column its full unbreakable width even at 5.5 pt, so headers such
        as alpha_mean_deg were rendered as "alpha_mean_de" + "g". Underscores and
        slashes are the natural break points; an explicit <br/> puts the break
        there instead of one character from the column edge."""
        esc = html.escape(txt)
        if not txt or pdfmetrics.stringWidth(txt, fn, pt) <= widths[j] - PAD:
            return esc
        parts = re.split(r"(?<=[_/])", txt)          # keep the separator on the left
        if len(parts) < 2:
            return esc
        out, cur = [], ""
        for part in parts:
            trial = cur + part
            if cur and pdfmetrics.stringWidth(trial, fn, pt) > widths[j] - PAD:
                out.append(cur); cur = part
            else:
                cur = trial
        if cur:
            out.append(cur)
        return "<br/>".join(html.escape(o) for o in out)

    rows = []
    for i, row in enumerate(text_rows):
        fn, pt = ("Body-Bold", hdr_pt) if i == 0 else ("Body", body_pt)
        cells = [Paragraph(_fit(row[j] if j < len(row) else "", j, fn, pt),
                           s_h if i == 0 else s_b) for j in range(ncol)]
        rows.append(cells)
    t = Table(rows, colWidths=widths, repeatRows=1, hAlign="CENTER")
    t.setStyle(TableStyle([
        ("GRID",       (0, 0), (-1, -1), 0.4, INK_SOFT),
        ("BACKGROUND", (0, 0), (-1, 0), HDR_BG),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, ROW_BG]),
        ("VALIGN",     (0, 0), (-1, -1), "MIDDLE"),
        ("LEFTPADDING", (0, 0), (-1, -1), PAD/2),
        ("RIGHTPADDING", (0, 0), (-1, -1), PAD/2),
        ("TOPPADDING", (0, 0), (-1, -1), 2),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 2),
    ]))
    return t


# ---------------------------------------------------------------- build story
def build_story():
    body = doc_src.element.body
    par_map = {p._p: p for p in doc_src.paragraphs}
    tbl_map = {t._tbl: t for t in doc_src.tables}
    story, pending = [], []          # pending: image(s) awaiting their caption

    def flush():
        if pending:
            story.append(KeepTogether(pending[:]) if len(pending) > 1 else pending[0])
            pending.clear()

    for child in body.iterchildren():
        if child.tag == qn("w:tbl"):
            flush()
            t = _table_flowable(tbl_map[child])
            if t is not None:
                story += [t, Spacer(1, 8)]
            continue
        if child.tag != qn("w:p"):
            continue
        par = par_map.get(child)
        if par is None:
            continue

        if par._p.findall(".//" + qn("w:br") + "[@" + qn("w:type") + "='page']"):
            flush(); story.append(PageBreak()); continue

        imgs = _image_flowables(par)
        if imgs:
            flush(); pending.extend(imgs); continue

        if pending and _is_caption(par):
            pending.append(Paragraph(_runs_markup(par, 9.5),
                                     _style("cap", 9.5, TA_CENTER, ACC, space_after=12)))
            flush(); continue
        flush()

        text = par.text
        if not text.strip():
            story.append(Spacer(1, 5)); continue

        sname = par.style.name
        if sname in _HEADING:
            story.append(Paragraph(html.escape(text), _HEADING[sname]))
        elif sname == "List Bullet":
            story.append(Paragraph(_runs_markup(par), S_BULLET,
                                   bulletText="•"))
        else:
            pt = _dominant_pt(par, BODY_PT)
            sa = par.paragraph_format.space_after
            st = _style("p%.1f" % pt, pt,
                        align=_ALIGN.get(par.alignment, TA_LEFT),
                        bold=False,
                        space_after=sa.pt if sa is not None else 6)
            story.append(Paragraph(_runs_markup(par, pt), st))
    flush()
    return story


# ---------------------------------------------------------------- page furniture
def _furniture(canvas, doc):
    canvas.saveState()
    canvas.setFont("Body", 8)
    canvas.setFillColor(INK_SOFT)
    canvas.drawCentredString(PAGE_W/2.0, MB*0.45, str(canvas.getPageNumber()))
    canvas.setStrokeColor(INK_SOFT); canvas.setLineWidth(0.4)
    canvas.line(ML, MB*0.72, PAGE_W-MR, MB*0.72)
    canvas.restoreState()


def render_docx(dest):
    d = BaseDocTemplate(str(dest), pagesize=(PAGE_W, PAGE_H),
                        leftMargin=ML, rightMargin=MR, topMargin=MT, bottomMargin=MB,
                        title=f"{TITLE} — {SOLVER}", author=AUTHOR,
                        subject="Dynamic stall / UIBS reduced-order solver")
    frame = Frame(ML, MB, FRAME_W, FRAME_H, id="body",
                  leftPadding=0, rightPadding=0, topPadding=0, bottomPadding=0)
    d.addPageTemplates([PageTemplate(id="main", frames=[frame], onPage=_furniture)])
    d.build(build_story())


if __name__ == "__main__":
    import tempfile, shutil
    parts = [HERE/"UNISTALL_data_dossier.pdf", HERE/"UNISTALL_plots_album.pdf"]
    missing = [p.name for p in parts if not p.exists()]
    if missing:
        raise SystemExit(f"[report] missing {missing} — run build_pdfs.py first")

    with tempfile.TemporaryDirectory() as td:
        main_pdf = Path(td)/"case_body.pdf"
        render_docx(main_pdf)
        merged = fitz.open()
        merged.insert_file(str(main_pdf)); n_body = merged.page_count
        counts = []
        for p in parts:
            with fitz.open(str(p)) as src:
                merged.insert_file(str(p)); counts.append((p.name, src.page_count))
        # ---- navigation. A 211-page report shipped with no contents page and no
        #      PDF outline at all, so a reader had no way to reach section 12 or
        #      the drawings except by scrolling. The outline is DERIVED from the
        #      headings already on the pages -- top-level "N. Title" in the body,
        #      then one entry per appendix volume and per album section divider --
        #      so it cannot disagree with the document it indexes.
        _toc = []
        _seen = set()
        for _i in range(n_body):
            # every heading on the page, not just the first: sections 2, 3 and 7
            # share a page with the one before them, and breaking after the first
            # match silently dropped all three from the outline.
            for _ln in merged[_i].get_text().splitlines():
                _m = re.match(r'^\s*(\d{1,2})\.\s+(\S.{2,70})$', _ln.strip())
                if _m and _m.group(1) not in _seen:
                    _seen.add(_m.group(1))
                    _toc.append([1, f"{_m.group(1)}. {_m.group(2).strip()}", _i + 1])
        _at = n_body
        for _name, _cnt in counts:
            _label = ("Appendix — Data & Tables Dossier" if "dossier" in _name.lower()
                      else "Appendix — Figures & Contours Album" if "album" in _name.lower()
                      else _name)
            _toc.append([1, _label, _at + 1])
            for _j in range(_at, _at + _cnt):          # album/dossier section dividers
                _lines = [l for l in merged[_j].get_text().splitlines() if l.strip()]
                if len(_lines) == 1 and len(_lines[0]) < 45:
                    _toc.append([2, _lines[0].strip(), _j + 1])
            _at += _cnt
        if _toc:
            merged.set_toc(_toc)
        merged.set_metadata({"title": f"{TITLE} — {SOLVER}", "author": AUTHOR,
                             "subject": "Prediction of dynamic stall on a helicopter "
                                        "main-rotor retreating blade (UIBS core)",
                             "creator": "07_report/build_report_pdf.py"})
        merged.save(str(OUT), garbage=4, deflate=True, clean=True)
        merged.close()
    shutil.rmtree(HERE/"_pdfimg", ignore_errors=True)
    with fitz.open(str(OUT)) as f:
        total = f.page_count
    print(f"[report] {OUT.name}: {n_body} body + "
          + " + ".join(f"{c} {n}" for n, c in counts) + f" = {total} pages "
          f"({OUT.stat().st_size/1e6:.1f} MB)")
