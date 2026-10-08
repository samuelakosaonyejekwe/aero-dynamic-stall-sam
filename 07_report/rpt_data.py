# -*- coding: utf-8 -*-
"""
07_report / rpt_data.py
-----------------------
Reads the files the report is written from and keeps a register of every
number that is quoted, with the file and the field it was read from. Nothing
here computes a load.

Author: Akosa Samuel Onyejekwe (independent)

Units: values are passed on in the units of the file they are read from.
"""
import importlib.util
import json
import re
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

# colours of the documents: deep navy for text and rules, muted tones elsewhere
NAVY = "#1F3350"
SOFT = "#4A5F7D"
RULE = "#9AA8BC"
BAND = "#E8EDF4"
ACCENT = "#2F6DB5"
WHITE = "#FFFFFF"

WORDS = ("zero", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten", "eleven", "twelve")


def sibling(name: str):
    """Load a helper module that lies beside this file, once."""
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, HERE/f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def as_text(value, fmt: str) -> str:
    """One value as the text that is printed. fmt is a format string, "{}" for
    text as it stands, "word" for a small whole number spelt out, or "pow10"
    for two significant figures with the power of ten written as such."""
    if fmt == "word":
        n = int(float(value))
        return WORDS[n] if 0 <= n < len(WORDS) else str(n)
    if fmt == "{}":
        return str(value)
    if fmt == "pow10":
        return power_of_ten(value, 2)
    if fmt.endswith("d}"):
        return fmt.format(int(round(float(value))))
    return fmt.format(float(value))


def md_tables(path: Path) -> dict:
    """Every pipe table of a Markdown file, keyed by the heading above it; a
    second table under one heading gets the key 'heading #2'."""
    tables, heading, rows, seen = {}, "", [], {}

    def flush() -> None:
        if len(rows) >= 2:
            seen[heading] = seen.get(heading, 0) + 1
            key = heading if seen[heading] == 1 else f"{heading} #{seen[heading]}"
            body = [r for r in rows[1:] if not set("".join(r)) <= set("-: ")]
            if all(len(r) == len(rows[0]) for r in body):      # a table with bare pipes in its cells is left out
                tables[key] = pd.DataFrame(body, columns=rows[0])
        rows.clear()

    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("|"):
            cells = [c.strip() for c in re.split(r"(?<!\\)\|", line.strip().strip("|"))]
            rows.append([c.replace("\\|", "|") for c in cells])
            continue
        flush()
        if line.startswith("#"):
            heading = line.lstrip("#").strip()
    flush()
    return tables


def formulation_equations(path: Path) -> dict:
    """The numbered equation lines of docs/formulation.md: number -> text."""
    out = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        m = re.match(r"\s{4}\((\d+[a-z]?)\)\s+(.*\S)", line)
        if m:
            out[m.group(1)] = m.group(2)
    return out


def bib_entries(path: Path) -> dict:
    """The entries of a BibTeX file: key -> dict of fields (plus 'type')."""
    text = path.read_text(encoding="utf-8")
    out = {}
    for m in re.finditer(r"@(\w+)\{([^,]+),(.*?)\n\}", text, re.S):
        fields = dict(type=m.group(1).lower())
        for f in re.finditer(r"^\s*(\w+)\s*=\s*\{(.*)\},?\s*$", m.group(3), re.M):
            fields[f.group(1).lower()] = f.group(2).replace("{", "").replace("}", "")
        out[m.group(2).strip()] = fields
    return out


class Data:
    """The files of the case study, read on demand, and the register of the
    numbers quoted from them."""

    def __init__(self) -> None:
        self.root = ROOT
        self.quoted = []            # one row per number quoted in the text
        self.tables_used = []       # one row per table shown whole
        self._cache = {}

    # ---- loading ---------------------------------------------------------
    def csv(self, rel: str) -> pd.DataFrame:
        """A CSV file of the repository as a table."""
        if rel not in self._cache:
            self._cache[rel] = pd.read_csv(ROOT/rel)
        return self._cache[rel]

    def json(self, rel: str) -> dict:
        """A JSON file of the repository."""
        if rel not in self._cache:
            self._cache[rel] = json.loads((ROOT/rel).read_text(encoding="utf-8"))
        return self._cache[rel]

    def md(self, rel: str) -> dict:
        """The pipe tables of a Markdown file of the repository."""
        if rel not in self._cache:
            self._cache[rel] = md_tables(ROOT/rel)
        return self._cache[rel]

    # ---- quoting ---------------------------------------------------------
    def note(self, value, rel: str, field: str, fmt: str = "{:.3f}") -> str:
        """Register a value read or derived from a file and return its text."""
        text = as_text(value, fmt)
        raw = value.item() if hasattr(value, "item") else value
        self.quoted.append(dict(text=text, value=raw, file=rel, field=field))
        return text

    def j(self, rel: str, dotted: str, fmt: str = "{:.3f}") -> str:
        """A value of a JSON file, addressed as 'a.b.0.c'."""
        node = self.json(rel)
        for part in dotted.split("."):
            node = node[int(part)] if isinstance(node, list) else node[part]
        return self.note(node, rel, dotted, fmt)

    def c(self, rel: str, col: str, fmt: str = "{:.3f}", **where) -> str:
        """The value of one column of a CSV file in the row that matches
        every 'column=value' given."""
        df = self.csv(rel)
        mask = pd.Series(True, index=df.index)
        for key, wanted in where.items():
            mask &= df[key].astype(str) == str(wanted)
        rows = df[mask]
        if len(rows) != 1:
            raise KeyError(f"{rel}: {where} matches {len(rows)} rows")
        field = col + "".join(f" [{k}={v}]" for k, v in where.items())
        return self.note(rows[col].iloc[0], rel, field, fmt)

    def whole(self, rel: str, what: str) -> None:
        """Register that a table of the document shows a file, or part of it."""
        self.tables_used.append(dict(file=rel, shown=what))


EXPONENT_TEXT = re.compile(r"[+-]?\d+(\.\d+)?[eE][+-]?\d+")    # a number a file holds as text with an exponent
SUPERSCRIPT = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def power_of_ten(value: float, digits: int = 3) -> str:
    """A number to `digits` significant figures, written with a power of ten
    where its size calls for one: 4.26 × 10⁶, 3.0 × 10⁻⁵."""
    text = f"{float(value):.{digits}g}"
    if "e" not in text:
        return text
    mantissa, exponent = text.split("e")
    return f"{mantissa} × 10{str(int(exponent)).translate(SUPERSCRIPT)}"


def cell(value, digits: int = 4) -> str:
    """A table cell as text: whole numbers without a decimal point, other
    numbers to `digits` significant figures with powers of ten written as
    such, booleans as yes or no."""
    if isinstance(value, bool) or type(value).__name__ == "bool_":
        return "yes" if value else "no"
    if value is None or (isinstance(value, float) and value != value):
        return ""
    if isinstance(value, (int, float)) or hasattr(value, "dtype"):
        v = float(value)
        if v == int(v) and abs(v) < 1e9:
            return str(int(v))
        return power_of_ten(v, digits)
    if isinstance(value, str) and EXPONENT_TEXT.fullmatch(value.strip()):
        return power_of_ten(float(value), digits)
    return str(value)


def frame_rows(df: pd.DataFrame, digits: int = 4) -> tuple:
    """A table as (header, rows) of text, headers with spaces for underscores."""
    header = [str(c).replace("_", " ") for c in df.columns]
    rows = [[cell(v, digits) for v in rec] for rec in df.itertuples(index=False, name=None)]
    return header, rows
