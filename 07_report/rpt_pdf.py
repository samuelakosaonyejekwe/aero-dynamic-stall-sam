# -*- coding: utf-8 -*-
"""
07_report / rpt_pdf.py
----------------------
Draws the report as a PDF with reportlab, and holds the page layout, styles
and table builder that the plots album and the data dossier share. No word
processor is needed.

Author: Akosa Samuel Onyejekwe (independent)

Units: lengths in points (1/72 inch).
"""
import io
import zlib
from pathlib import Path

import matplotlib
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_RIGHT
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import cm, inch
from reportlab.pdfbase import pdfdoc, pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (BaseDocTemplate, CondPageBreak, Flowable, Frame, Image, KeepTogether, NextPageTemplate,
                                PageBreak, PageTemplate, Paragraph, Spacer, Table, TableStyle)

import project_meta as meta
from rpt_blocks import column_widths, escape, strip_tags, words
from rpt_data import ACCENT, BAND, NAVY, RULE, SOFT

FONT = "ReportSans"
FACES = {"": "DejaVuSans.ttf", "-Bold": "DejaVuSans-Bold.ttf", "-Italic": "DejaVuSans-Oblique.ttf",
         "-BoldItalic": "DejaVuSans-BoldOblique.ttf"}
MARGIN_X, MARGIN_TOP, MARGIN_BOTTOM = 2.0*cm, 2.3*cm, 2.1*cm
SHEET_MARGIN = 1.1*cm        # margins of the landscape page a drawing sheet is given
SHEET_CAPTION = 1.5*cm       # height kept for the caption under a drawing sheet
TABLE_START = 3.6*cm       # room a table needs to begin on a page: its caption, its head and a few rows
HALF_FULL = 0.5              # a page filled to less than this fraction is entered as short
BODY_PT = 10.0
MIN_TABLE_PT = 5.5
INK, MUTED, LINE, TINT, MARK = (colors.HexColor(c) for c in (NAVY, SOFT, RULE, BAND, ACCENT))


def register_fonts() -> None:
    """Register the type family (DejaVu Sans, shipped with matplotlib)."""
    if FONT in pdfmetrics.getRegisteredFontNames():
        return
    folder = Path(matplotlib.get_data_path())/"fonts"/"ttf"
    for suffix, fname in FACES.items():
        pdfmetrics.registerFont(TTFont(FONT + suffix, str(folder/fname)))
    pdfmetrics.registerFontFamily(FONT, normal=FONT, bold=FONT + "-Bold", italic=FONT + "-Italic",
                                  boldItalic=FONT + "-BoldItalic")


def style(name: str, size: float, **kw) -> ParagraphStyle:
    """A paragraph style in the document's type and navy ink."""
    kw.setdefault("textColor", INK)
    kw.setdefault("leading", size*1.38)
    kw.setdefault("fontName", FONT)
    return ParagraphStyle(name, fontSize=size, **kw)


def styles() -> dict:
    """The paragraph styles of the documents."""
    register_fonts()
    bold = FONT + "-Bold"
    return dict(
        body=style("body", BODY_PT, spaceAfter=6, alignment=TA_LEFT),
        bullet=style("bullet", BODY_PT, spaceAfter=4, leftIndent=16, bulletIndent=4),
        note=style("note", BODY_PT, spaceAfter=0),
        h1=style("h1", 15, fontName=bold, spaceBefore=16, spaceAfter=8, keepWithNext=1),
        h2=style("h2", 11.5, fontName=bold, spaceBefore=10, spaceAfter=5, keepWithNext=1),
        caption=style("caption", 8.5, textColor=MUTED, alignment=TA_CENTER, spaceBefore=4, spaceAfter=12),
        table_caption=style("table_caption", 8.5, textColor=MUTED, alignment=TA_LEFT, spaceBefore=6,
                            spaceAfter=4),
        eq_number=style("eq_number", BODY_PT, alignment=TA_RIGHT),
        title=style("title", 24, fontName=bold, leading=31, alignment=TA_LEFT, spaceAfter=10),
        subtitle=style("subtitle", 13.5, leading=19, textColor=MUTED, spaceAfter=6),
        cover=style("cover", 11, leading=16, spaceAfter=4),
        cover_small=style("cover_small", 9.5, leading=14, textColor=MUTED, spaceAfter=3),
        kind=style("kind", 10, fontName=bold, textColor=MARK, spaceAfter=14),
    )


def markup(text: str) -> str:
    """The report's inline markup in reportlab's dialect."""
    text = text.replace("<sub>", '<sub rise="2" size="7.5">')
    return text.replace("<sup>", '<super rise="4" size="7.5">').replace("</sup>", "</super>")


def furniture(left: str, right: str):
    """A page-drawing function: running head, rule, author and page number."""
    def draw(canvas, doc) -> None:
        w, h = doc.pagesize
        canvas.saveState()
        canvas.setStrokeColor(LINE)
        canvas.setLineWidth(0.6)
        canvas.line(MARGIN_X, h - MARGIN_TOP + 10, w - MARGIN_X, h - MARGIN_TOP + 10)
        canvas.line(MARGIN_X, MARGIN_BOTTOM - 10, w - MARGIN_X, MARGIN_BOTTOM - 10)
        canvas.setFillColor(MUTED)
        canvas.setFont(FONT, 8)
        canvas.drawString(MARGIN_X, h - MARGIN_TOP + 15, left)
        canvas.drawRightString(w - MARGIN_X, h - MARGIN_TOP + 15, right)
        canvas.drawString(MARGIN_X, MARGIN_BOTTOM - 21, meta.AUTHOR_FULL)
        canvas.drawRightString(w - MARGIN_X, MARGIN_BOTTOM - 21, f"{doc.page}")
        canvas.restoreState()
    return draw


class Document(BaseDocTemplate):
    """A document that lets text flow past a figure that does not fit, and
    keeps a record of how full each page is.

    A flowable marked `floating` that does not fit in what is left of the
    page is held back and placed at the top of the next page, and the text
    after it moves up to fill the page. `page_fill` holds one row per page:
    its number, the fraction of the frame that is filled, and whether the
    page ends at a break asked for by the document (the end of a part)."""

    def __init__(self, *args: object, **kw: object) -> None:
        super().__init__(*args, **kw)
        self.page_fill, self._held, self._break_asked = [], [], False

    def handle_flowable(self, flowables: list) -> None:
        """Place the next flowable, holding back a float that does not fit.
        Held figures are placed at the top of the next ordinary page, and
        before any break the document asks for."""
        frame = getattr(self, "frame", None)
        f = flowables[0]
        if self._held and frame is not None:
            if frame._atTop and self.pageTemplate.id != "sheet":
                flowables[0:0] = self._held
                self._held = []
                f = flowables[0]
            elif isinstance(f, (PageBreak, NextPageTemplate)) or getattr(f, "last", False):
                flowables[0:0] = [PageBreak()]
                super().handle_flowable(flowables)
                return
        if isinstance(f, PageBreak):
            self._break_asked = True
        if getattr(f, "floating", False) and frame is not None and not frame._atTop:
            if f.wrapOn(self.canv, frame._aW, frame._aH)[1] > frame._aH:
                self._held.append(flowables.pop(0))
                return
        super().handle_flowable(flowables)

    def afterPage(self) -> None:
        """Enter the page that has just ended in the record."""
        frame = self.frame
        used = 1.0 - max(frame._y - frame._y1p, 0.0)/(frame._y2 - frame._y1)
        self.page_fill.append(dict(page=self.page, filled=round(used, 3), ends_a_part=self._break_asked))
        self._break_asked = False


class End(Spacer):
    """The last flowable of a story: it makes the document place any figure
    it is still holding back."""
    last = True

    def __init__(self) -> None:
        super().__init__(1, 0.1)


def short_pages(doc: Document, name: str) -> list:
    """Rows of the page record of a built document, with the document's name
    and whether the page counts as short: less than HALF_FULL filled, and
    neither the cover, nor the end of a part, nor the last page."""
    rows = []
    for row in doc.page_fill:
        last = row["page"] == len(doc.page_fill)
        rows.append(dict(document=name, **row, last_page=last,
                         short=bool(row["filled"] < HALF_FULL and row["page"] > 1 and not row["ends_a_part"]
                                    and not last)))
    return rows


def make_doc(path: Path, title: str, subject: str, pagesize: tuple = A4) -> Document:
    """A document with a bare cover page, furnished body pages and landscape
    pages with thin margins for drawing sheets."""
    register_fonts()
    doc = Document(str(path), pagesize=pagesize, leftMargin=MARGIN_X, rightMargin=MARGIN_X,
                   topMargin=MARGIN_TOP, bottomMargin=MARGIN_BOTTOM, title=title, author=meta.AUTHOR,
                   subject=subject, creator=meta.AUTHOR)
    w, h = pagesize
    frame = Frame(MARGIN_X, MARGIN_BOTTOM, w - 2*MARGIN_X, h - MARGIN_TOP - MARGIN_BOTTOM, id="f", leftPadding=0,
                  rightPadding=0, topPadding=0, bottomPadding=0)
    lw, lh = landscape(A4)
    sheet = Frame(SHEET_MARGIN, SHEET_MARGIN, lw - 2*SHEET_MARGIN, lh - 2*SHEET_MARGIN, id="s", leftPadding=0,
                  rightPadding=0, topPadding=0, bottomPadding=0)
    doc.addPageTemplates([PageTemplate(id="cover", frames=[frame]),
                          PageTemplate(id="body", frames=[frame], onPage=furniture(meta.TITLE, subject)),
                          PageTemplate(id="sheet", frames=[sheet], pagesize=(lw, lh), onPage=sheet_furniture)])
    doc.avail_w, doc.avail_h = frame._width, frame._height
    doc.sheet_w, doc.sheet_h = sheet._width, sheet._height
    return doc


def sheet_furniture(canvas, doc) -> None:
    """The page number of a drawing-sheet page, in its bottom margin."""
    w, _ = landscape(A4)
    canvas.saveState()
    canvas.setFillColor(MUTED)
    canvas.setFont(FONT, 8)
    canvas.drawString(SHEET_MARGIN, 0.45*cm, meta.AUTHOR_FULL)
    canvas.drawRightString(w - SHEET_MARGIN, 0.45*cm, f"{doc.page}")
    canvas.restoreState()


def cover(st: dict, kind_line: str, extra: list, avail_h: float) -> list:
    """The flowables of a cover page, spaced for a page `avail_h` high."""
    unit = avail_h/100.0
    rule = Table([[""]], colWidths=[6*cm], rowHeights=[3])
    rule.setStyle(TableStyle([("BACKGROUND", (0, 0), (-1, -1), MARK)]))
    rule.hAlign = "LEFT"
    out = [Spacer(1, 21*unit), Paragraph(escape(kind_line), st["kind"]), Paragraph(escape(meta.TITLE), st["title"]),
           Paragraph(escape(meta.SUBTITLE), st["subtitle"]), Spacer(1, 10), rule, Spacer(1, 18)]
    out += [Paragraph(escape(line), st["cover"]) for line in (meta.SOLVER_LONG, meta.METHOD)]
    out += [Spacer(1, 9*unit), Paragraph(escape(meta.AUTHOR_FULL), st["cover"]),
            Paragraph(escape(meta.STUDY_DATE), st["cover"]), Spacer(1, 5*unit)]
    out += [Paragraph(line, st["cover_small"]) for line in extra]
    return out + [NextPageTemplate("body"), PageBreak()]


def table_flowable(header: list, rows: list, avail: float, size: float = 8.5, rich: bool = False) -> Table:
    """A table that fits the width `avail`: wrapped cells, a tinted header row
    repeated on every page, muted rules. The type is made smaller until the
    longest words fit. With `rich` the cells carry the inline markup."""
    shown = rows
    if rich:
        rows = [[strip_tags(c) for c in row] for row in rows]
    widths = None
    while widths is None:
        def measure(text: str, pt: float = size) -> float:
            return pdfmetrics.stringWidth(text, FONT + "-Bold", pt)
        widths = column_widths(header, rows, avail - 2, measure)
        if widths is None and size <= MIN_TABLE_PT:
            floor = [max(measure(w) for c in col for w in words(c)) for col in zip(header, *rows, strict=True)]
            widths = [avail*f/sum(floor) for f in floor]
        elif widths is None:
            size -= 0.5
    head = style("th", size, fontName=FONT + "-Bold", leading=size*1.25)
    body = style("td", size, leading=size*1.25)
    data = [[Paragraph(escape(c), head) for c in header]]
    data += [[Paragraph(markup(c) if rich else escape(c), body) for c in row] for row in shown]
    table = Table(data, colWidths=widths, repeatRows=1)
    table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), TINT), ("LINEBELOW", (0, 0), (-1, 0), 0.8, INK),
        ("LINEABOVE", (0, 0), (-1, 0), 0.8, INK), ("LINEBELOW", (0, 1), (-1, -1), 0.3, LINE),
        ("LINEBELOW", (0, -1), (-1, -1), 0.8, INK), ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("LEFTPADDING", (0, 0), (-1, -1), 4), ("RIGHTPADDING", (0, 0), (-1, -1), 4),
        ("TOPPADDING", (0, 0), (-1, -1), 2.5), ("BOTTOMPADDING", (0, 0), (-1, -1), 3.5),
    ]))
    return table


def note_flowable(text: str, st: dict, avail: float) -> Table:
    """A statement in a tinted band with a coloured bar at its left."""
    table = Table([[Paragraph(markup(text), st["note"])]], colWidths=[avail])
    table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), TINT), ("LINEBEFORE", (0, 0), (0, -1), 3, MARK),
        ("LEFTPADDING", (0, 0), (-1, -1), 10), ("RIGHTPADDING", (0, 0), (-1, -1), 8),
        ("TOPPADDING", (0, 0), (-1, -1), 7), ("BOTTOMPADDING", (0, 0), (-1, -1), 8),
    ]))
    return table


class HexString(pdfdoc.PDFObject):
    """Bytes written into the PDF as a hexadecimal string."""

    def __init__(self, data: bytes) -> None:
        self.data = data

    def format(self, document) -> bytes:
        """The PDF text of the string."""
        return b"<" + self.data.hex().encode("ascii") + b">"


class PaletteXObject(pdfdoc.PDFObject):
    """A 256-colour image as a PDF image with an indexed colour space, one
    byte per pixel, deflated."""

    def __init__(self, image) -> None:
        self.size = image.size
        self.palette = bytes(image.getpalette()[:768])
        self.content = zlib.compress(image.tobytes(), 9)

    def format(self, document) -> bytes:
        """The PDF text of the image object."""
        stream = pdfdoc.PDFStream(content=self.content, filters=[])
        d = stream.dictionary
        d["Type"], d["Subtype"] = pdfdoc.PDFName("XObject"), pdfdoc.PDFName("Image")
        d["Width"], d["Height"], d["BitsPerComponent"] = self.size[0], self.size[1], 8
        d["Filter"] = pdfdoc.PDFName("FlateDecode")
        d["ColorSpace"] = pdfdoc.PDFArray([pdfdoc.PDFName("Indexed"), pdfdoc.PDFName("DeviceRGB"),
                                           len(self.palette)//3 - 1, HexString(self.palette)])
        return stream.format(document)


class PaletteImage(Flowable):
    """A 256-colour image placed at a given size."""

    def __init__(self, image, width: float, height: float) -> None:
        super().__init__()
        self.image, self.width, self.height = image, width, height
        self.name = f"Pal{id(image):x}"

    def wrap(self, avail_w: float, avail_h: float) -> tuple:
        """The space the image takes."""
        return self.width, self.height

    def draw(self) -> None:
        """Register the image with the document once, and paint it."""
        canvas = self.canv
        reg = canvas._doc.getXObjectName(self.name)
        if reg not in canvas._doc.idToObject:
            obj = PaletteXObject(self.image)
            canvas._doc.Reference(obj, reg)
            canvas._doc.addForm(self.name, obj)
        canvas.saveState()
        canvas.scale(self.width, self.height)
        canvas._code.append(f"/{reg} Do")
        canvas.restoreState()
        canvas._formsinuse.append(self.name)


def picture(pics, path: Path, max_w: float, max_h: float) -> Flowable:
    """An image scaled to fit a box: a JPEG as it is, or a palette image."""
    w, h = pics.fit(path, max_w, max_h)
    kind, data = pics.prepared(path, w/inch, "pdf")
    if kind == "jpeg":
        return Image(io.BytesIO(data), width=w, height=h)
    return PaletteImage(data, w, h)


def equation_flowable(block: dict, st: dict, avail: float) -> Table:
    """An equation centred in the column with its number at the right."""
    img = Image(str(block["path"]), width=block["width_in"]*inch, height=block["height_in"]*inch)
    number = Paragraph(f"({block['number']})", st["eq_number"])
    table = Table([["", img, number]], colWidths=[1.2*cm, avail - 2.4*cm, 1.2*cm])
    table.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "MIDDLE"), ("ALIGN", (1, 0), (1, 0), "CENTER"),
                               ("TOPPADDING", (0, 0), (-1, -1), 3), ("BOTTOMPADDING", (0, 0), (-1, -1), 7),
                               ("LEFTPADDING", (0, 0), (-1, -1), 0), ("RIGHTPADDING", (0, 0), (-1, -1), 0)]))
    return table


def figure_sizes(block: dict, pics, avail: float, sheet_box: tuple) -> list:
    """(path, width, height) in points of each image of a figure block as it
    is placed: within the column, or for a drawing sheet within `sheet_box`,
    the (width, height) of its landscape page."""
    paths = block["paths"]
    if block.get("sheet"):
        return [(paths[0], *pics.fit(paths[0], sheet_box[0], sheet_box[1] - SHEET_CAPTION))]
    gap = 8.0
    each = (avail*block["width"] - gap*(len(paths) - 1))/len(paths)
    return [(p, *pics.fit(p, each, block["max_h_in"]*inch)) for p in paths]


def figure_flowable(block: dict, st: dict, pics, avail: float, sheet_box: tuple) -> list:
    """One image, or a row of images, with its caption. An ordinary figure is
    kept on a page with its caption and floats: if it does not fit where it
    falls, it goes to the top of the next page and the text flows past it. A
    drawing sheet is given a landscape page to itself."""
    sizes = figure_sizes(block, pics, avail, sheet_box)
    images = []
    for path, w, h in sizes:
        kind, data = pics.prepared(path, w/inch, "pdf")
        images.append(Image(io.BytesIO(data), width=w, height=h) if kind == "jpeg" else PaletteImage(data, w, h))
    caption = Paragraph(markup(block["caption"]), st["caption"])
    if block.get("sheet"):
        images[0].hAlign = "CENTER"
        return [NextPageTemplate("sheet"), PageBreak(), images[0], caption, NextPageTemplate("body")]
    if len(images) == 1:
        body = images[0]
        body.hAlign = "CENTER"
    else:
        body = Table([images], colWidths=[w + 8.0 for _, w, _ in sizes])
        body.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "MIDDLE"), ("ALIGN", (0, 0), (-1, -1), "CENTER"),
                                  ("LEFTPADDING", (0, 0), (-1, -1), 0), ("RIGHTPADDING", (0, 0), (-1, -1), 0)]))
    kept = KeepTogether([Spacer(1, 4), body, caption])
    kept.floating = True
    return [kept]


def block_flowables(block: dict, st: dict, pics, avail: float, avail_h: float, sheet_box: tuple) -> list:
    """The flowables of one block of the report."""
    kind = block["kind"]
    simple = {"p": "body", "h2": "h2"}
    if kind in simple:
        return [Paragraph(markup(block["text"]), st[simple[kind]])]
    if kind == "h1":
        return [CondPageBreak(4*cm), Paragraph(markup(block["text"]), st["h1"])]
    if kind == "bullet":
        return [Paragraph(markup(block["text"]), st["bullet"], bulletText="•")]
    if kind == "note":
        return [Spacer(1, 3), note_flowable(block["text"], st, avail), Spacer(1, 9)]
    if kind == "eq":
        return [equation_flowable(block, st, avail)]
    if kind == "fig":
        return figure_flowable(block, st, pics, avail, sheet_box)
    if kind == "table":
        return [CondPageBreak(TABLE_START), Paragraph(markup(block["caption"]), st["table_caption"]),
                table_flowable(block["header"], block["rows"], avail, block["size"], block.get("rich", False)),
                Spacer(1, 12)]
    if kind == "break":
        return [PageBreak()]
    if kind == "title":
        return cover(st, block["kind_line"], [], avail_h)
    raise ValueError(f"unknown block '{kind}'")


def render_report(report, pics, path: Path) -> tuple:
    """Write the report as a PDF. Returns (rows of the page record, rows of
    the placed figures: block caption number, path, width placed [in])."""
    st = styles()
    doc = make_doc(path, meta.TITLE, f"{meta.SOLVER} case-study report")
    sheet_box = (doc.sheet_w, doc.sheet_h)
    story, placed = [], []
    after_sheet = False
    for block in report.blocks:
        sheet = bool(block.get("sheet"))
        if after_sheet and not sheet:
            story.append(PageBreak())               # the text goes on from a fresh page after a drawing sheet
        after_sheet = sheet
        story += block_flowables(block, st, pics, doc.avail_w, doc.avail_h, sheet_box)
        if block["kind"] == "fig":
            placed += [dict(path=p, width_in=w/inch) for p, w, _ in figure_sizes(block, pics, doc.avail_w, sheet_box)]
    doc.build(story + [End()])
    return short_pages(doc, path.name), placed
