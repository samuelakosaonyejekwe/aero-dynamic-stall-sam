# -*- coding: utf-8 -*-
"""
07_report / rpt_blocks.py
-------------------------
The report as a list of blocks (headings, paragraphs, equations, figures,
tables), written once and drawn twice: as a Word document and as a PDF. Also
the pieces both drawings share: the inline markup, the column widths of a
table and the images scaled for the page.

Author: Akosa Samuel Onyejekwe (independent)

Units: lengths in points (1/72 inch) unless a name says inches or pixels.
"""
import io
import re
import zlib
from pathlib import Path

from PIL import Image

MAX_PIXELS = 1400            # longest side of an embedded image, pixels
JPEG_QUALITY = 85
PALETTE_BIAS = 1.15          # a palette image is kept if it is no more than this many times the JPEG
TAG = re.compile(r"<(/?)(b|i|sub|sup)>")
COMPUTER_NOTATION = re.compile(r"\d[eE][+-]\d")      # 1.2e+06, 3e-05: not to appear in a table cell
ENTITIES = (("&lt;", "<"), ("&gt;", ">"), ("&amp;", "&"))


def runs(text: str) -> list:
    """Split text carrying <b>, <i>, <sub>, <sup> into (text, set of tags)."""
    out, on, pos = [], set(), 0
    for m in TAG.finditer(text):
        if m.start() > pos:
            out.append((plain(text[pos:m.start()]), frozenset(on)))
        if m.group(1):
            on.discard(m.group(2))
        else:
            on.add(m.group(2))
        pos = m.end()
    if pos < len(text):
        out.append((plain(text[pos:]), frozenset(on)))
    return out


def plain(text: str) -> str:
    """Text with its entities turned back into characters."""
    for ent, ch in ENTITIES:
        text = text.replace(ent, ch)
    return text


def strip_tags(text: str) -> str:
    """Text without its markup."""
    return plain(TAG.sub("", text))


def escape(text: str) -> str:
    """Text from a data file made safe to stand inside the markup."""
    return str(text).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def words(text: str) -> list:
    """The pieces of a cell between which a line may break (a no-break space
    is not such a place)."""
    return re.split(r"[ \n\t]+", text) or [""]


def column_widths(header: list, rows: list, avail: float, measure) -> list | None:
    """Widths that fit `avail`: every column at least as wide as its longest
    word, the rest shared in proportion to the longest cell. `measure(text)`
    returns a width in the units of `avail`. None if the words do not fit."""
    pad = measure("nn")
    cols = list(zip(header, *rows, strict=True))
    natural = [max(measure(c) for c in col) + pad for col in cols]
    floor = [max(measure(w) for c in col for w in words(c)) + pad for col in cols]
    if sum(natural) <= avail:
        grow = avail/sum(natural)
        return [n*min(grow, 1.6) for n in natural]
    if sum(floor) > avail:
        return None
    spare, want = avail - sum(floor), [n - f for n, f in zip(natural, floor, strict=True)]
    return [f + spare*w/sum(want) for f, w in zip(floor, want, strict=True)]


class Pictures:
    """Images prepared for embedding: flattened onto white, reduced to the
    size they are shown at, and stored as JPEG or with a 256-colour palette,
    whichever is smaller for the kind of file they go into. Each prepared
    image is kept, so an image shown twice at one size is encoded once."""

    def __init__(self, dpi: float) -> None:
        self.dpi = dpi                  # pixels per inch of the page an image is given
        self.used = dict(jpeg=0, palette=0)
        self._flat, self._kept = {}, {}

    def flat(self, path: Path) -> Image.Image:
        """The image on a white ground, at its full size."""
        key = str(path)
        if key not in self._flat:
            with Image.open(path) as im:
                im = im.convert("RGBA")
                flat = Image.new("RGB", im.size, "white")
                flat.paste(im, mask=im.split()[3])
            self._flat = {key: flat}    # only the last one is kept: the files are large
        return self._flat[key]

    def size(self, path: Path) -> tuple:
        """(width, height) of the source image in pixels."""
        with Image.open(path) as im:
            return im.size

    def fit(self, path: Path, max_w: float, max_h: float) -> tuple:
        """(width, height) that fits the box and keeps the proportions."""
        w, h = self.size(path)
        scale = min(max_w/w, max_h/h)
        return w*scale, h*scale

    def prepared(self, path: Path, shown_w_in: float, target: str) -> tuple:
        """(kind, data) for an image shown `shown_w_in` inches wide in a
        'pdf' or a 'docx'. kind 'jpeg': data is the bytes of a JPEG file.
        kind 'palette': data is a 256-colour image (pdf) or the bytes of
        its PNG file (docx)."""
        w, h = self.size(path)
        scale = min(1.0, shown_w_in*self.dpi/w, MAX_PIXELS/max(w, h))
        px = (max(1, round(w*scale)), max(1, round(h*scale)))
        key = (str(path), px, target)
        if key not in self._kept:
            small = self.flat(path).resize(px, Image.LANCZOS) if scale < 1.0 else self.flat(path)
            self._kept[key] = encode(small, target)
            self.used[self._kept[key][0]] += 1
        return self._kept[key]


def encode(image: Image.Image, target: str) -> tuple:
    """Encode an RGB image as JPEG or with a 256-colour palette, whichever
    is smaller in the target file; the palette, which keeps lines and
    lettering sharp, is preferred up to PALETTE_BIAS times the JPEG."""
    buf = io.BytesIO()
    image.save(buf, "JPEG", quality=JPEG_QUALITY)
    jpeg = buf.getvalue()
    paletted = image.quantize(colors=256, dither=Image.Dither.NONE)
    if target == "pdf":
        if len(zlib.compress(paletted.tobytes(), 6)) > PALETTE_BIAS*len(jpeg):
            return "jpeg", jpeg
        return "palette", paletted
    buf = io.BytesIO()
    paletted.save(buf, "PNG", compress_level=9)
    if len(buf.getvalue()) > PALETTE_BIAS*len(jpeg):
        return "jpeg", jpeg
    return "palette", buf.getvalue()


class Report:
    """The blocks of the report, in order, with figures, tables, sections and
    references numbered as they are added."""

    def __init__(self, equations: dict, bib: dict) -> None:
        self.blocks = []
        self.equations = equations
        self.bib = bib
        self.cited = []
        self._n = dict(fig=0, tab=0, h1=0, h2=0, eq=0)
        self._sections = {}

    def title(self, **fields) -> None:
        """The title page."""
        self.blocks.append(dict(kind="title", **fields))

    def h1(self, text: str, label: str = "") -> None:
        """A numbered section heading; `label` lets the text refer to it."""
        self._n["h1"] += 1
        self._n["h2"] = 0
        if label:
            self._sections[label] = str(self._n["h1"])
        self.blocks.append(dict(kind="h1", text=f"{self._n['h1']}  {text}"))

    def h2(self, text: str) -> None:
        """A numbered subsection heading."""
        self._n["h2"] += 1
        self.blocks.append(dict(kind="h2", text=f"{self._n['h1']}.{self._n['h2']}  {text}"))

    def p(self, text: str) -> None:
        """A paragraph."""
        self.blocks.append(dict(kind="p", text=text))

    def bullets(self, items: list) -> None:
        """A list of points."""
        for item in items:
            self.blocks.append(dict(kind="bullet", text=item))

    def note(self, text: str) -> None:
        """A statement set apart in a tinted band."""
        self.blocks.append(dict(kind="note", text=text))

    def eq(self, key: str) -> str:
        """An equation, by its key; returns its number as '(n)'."""
        e = self.equations[key]
        if e["number"] != self._n["eq"] + 1:
            raise ValueError(f"equation '{key}' is used out of the order it is numbered in")
        self._n["eq"] = e["number"]
        self.blocks.append(dict(kind="eq", **e))
        return f"({e['number']})"

    def eq_ref(self, key: str) -> str:
        """The number of an equation as '(n)'."""
        return f"({self.equations[key]['number']})"

    def fig(self, paths, caption: str, width: float = 1.0, max_h_in: float = 7.6, sheet: bool = False) -> str:
        """One image, or two side by side, with a numbered caption. `width` is
        the fraction of the column the figure takes. A `sheet` is given a
        landscape page to itself. Returns 'Figure n'."""
        self._n["fig"] += 1
        paths = [paths] if isinstance(paths, (str, Path)) else list(paths)
        self.blocks.append(dict(kind="fig", paths=[Path(p) for p in paths], width=width, max_h_in=max_h_in,
                                sheet=sheet, caption=f"<b>Figure {self._n['fig']}.</b> {caption}"))
        return f"Figure {self._n['fig']}"

    def next_fig(self, ahead: int = 1) -> str:
        """The name the next figure (or one further on) will have."""
        return f"Figure {self._n['fig'] + ahead}"

    def next_tab(self, ahead: int = 1) -> str:
        """The name the next table (or one further on) will have."""
        return f"Table {self._n['tab'] + ahead}"

    def table(self, header: list, rows: list, caption: str, size: float = 8.5, rich: bool = False) -> str:
        """A table with a numbered caption above it. With `rich` its cells
        carry the inline markup of the text (and must be escaped by the
        caller). Returns 'Table n'."""
        self._n["tab"] += 1
        for row in rows:
            for c in row:
                if COMPUTER_NOTATION.search(str(c)) or "{" in str(c):
                    raise ValueError(f"table '{caption[:60]}': cell '{c}' is in computer notation")
        self.blocks.append(dict(kind="table", header=[str(h) for h in header], size=size, rich=rich,
                                rows=[[str(c) for c in r] for r in rows],
                                caption=f"<b>Table {self._n['tab']}.</b> {caption}"))
        return f"Table {self._n['tab']}"

    def sec(self, label: str) -> str:
        """A mark that resolve() turns into the number of the labelled section."""
        return f"@@{label}@@"

    def resolve(self) -> None:
        """Put the section numbers in place of the marks left by sec()."""
        def fill(text: str) -> str:
            return re.sub(r"@@(\w+)@@", lambda m: self._sections[m.group(1)], text)
        for block in self.blocks:
            for key in ("text", "caption"):
                if key in block:
                    block[key] = fill(block[key])

    def page_break(self) -> None:
        """Start a new page."""
        self.blocks.append(dict(kind="break"))

    def cite(self, *keys: str) -> str:
        """Reference numbers for entries of the bibliography, as '[n]' or '[n, m]'."""
        nums = []
        for key in keys:
            if key not in self.bib:
                raise KeyError(f"docs/references.bib has no entry '{key}'")
            if key not in self.cited:
                self.cited.append(key)
            nums.append(str(self.cited.index(key) + 1))
        return "[" + ", ".join(nums) + "]"

    def reference_lines(self) -> list:
        """The cited entries, formatted, in the order they were first cited."""
        return [f"[{i}] {format_entry(self.bib[key])}" for i, key in enumerate(self.cited, start=1)]


def format_entry(e: dict) -> str:
    """One bibliography entry as a line of text."""
    parts = [f"{escape(e.get('author', ''))} ({e.get('year', '')}). {escape(e.get('title', '').strip())}."]
    if e.get("journal"):
        where = f"<i>{escape(e['journal'])}</i>"
        if e.get("volume"):
            where += f" {e['volume']}"
        if e.get("number"):
            where += f"({e['number']})"
        if e.get("pages"):
            where += f", {e['pages'].replace('--', '-')}"
        parts.append(where + ".")
    elif e.get("institution"):
        parts.append(f"{escape(e['institution'])}, {escape(e.get('number', ''))}.")
    if e.get("doi"):
        parts.append(f"doi:{escape(e['doi'])}")
    elif e.get("url"):
        parts.append(escape(e["url"]))
    return " ".join(parts)
