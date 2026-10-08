# -*- coding: utf-8 -*-
"""
unistall / code_quality.py
----------------------
Measures for the code-quality work package, written to results/code_quality.csv,
and the numerical regression reference every refactoring step is held to.

Author: Akosa Samuel Onyejekwe (independent)

  python3 -m unistall.code_quality --snapshot   store the load model's outputs for a
                                            fixed set of cases
                                            (results/regression_reference.npz)
  python3 -m unistall.code_quality              compare the current outputs with
                                            that reference and measure the code
"""
import ast
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
from unistall.paths import RESULTS
ROOT = HERE.parent
REFERENCE = RESULTS/"regression_reference.npz"
PACKAGE = ["attached_flow.py", "static_model.py", "dsmodel.py", "reference_lb.py", "metrics.py",
           "strokes.py", "structural.py", "flowfield.py", "naca4.py", "paths.py", "statespace.py"]
#  mean, amplitude (deg), k, M
CASES = [(10.0, 10.0, 0.10, 0.30), (15.0, 10.0, 0.05, 0.29), (9.0, 5.0, 0.20, 0.30),
         (5.0, 5.0, 0.10, 0.30), (15.0, 10.0, 0.10, 0.22), (12.0, 8.0, 0.074, 0.279)]
KEYS = ("CL", "CD", "CM", "CN", "CC", "f_sep", "CN_vortex", "tau_v")
HISTORY = r"used to|previously|earlier revision|until this audit|was wrong|no longer"


def _outputs():
    from unistall import dsmodel as dm, reference_lb as ref, attached_flow as af
    out = {}
    lit = dict(dm.DEFAULTS)
    for i, (a0, a1, k, M) in enumerate(CASES):
        o = dm.solve(a0, a1, k, M, consts=lit)
        r = dm.solve(a0, a1, k, M, consts=lit, static=ref.ExponentialStatic(M))
        for key in KEYS:
            out[f"tab{i}_{key}"] = o[key]
            out[f"ref{i}_{key}"] = r[key]
        out[f"tab{i}_onset"] = np.array([o["onset_alpha_deg"]])
    t = np.linspace(0.0, 0.5, 2001)
    a = af.solve_attached(np.radians(5 + 5*np.sin(60*t)), t, 102.0, 0.61, 0.30, 6.245, x_pitch=0.40)
    out["att_CN"], out["att_CM"] = a["CN"], a["CM"]
    return out


def regression_difference():
    ref = np.load(REFERENCE)
    cur = _outputs()
    return max(float(np.nanmax(np.abs(cur[k] - ref[k]))) if np.isfinite(ref[k]).any() else 0.0 for k in ref.files)


def measure():
    rows = {}
    rows["regression_max_abs_difference"] = regression_difference()
    files = [HERE/f for f in PACKAGE if (HERE/f).exists()]
    load_path = sorted(HERE.glob("*.py"))
    rows["import_path_edits_load_model_path"] = sum(len(re.findall(r"sys\.path\.(insert|append)", p.read_text())) for p in load_path)
    worst_exec, longest, unused, typed, public, hist = 0, 0, [], 0, 0, 0
    for p in load_path:
        tree = ast.parse(p.read_text())
        body = [s for s in tree.body if not isinstance(s, (ast.Import, ast.ImportFrom, ast.FunctionDef, ast.ClassDef))
                and not (isinstance(s, ast.Expr) and isinstance(s.value, ast.Constant))
                and not (isinstance(s, ast.If) and "__main__" in ast.unparse(s.test))]
        n_exec = sum((s.end_lineno - s.lineno + 1) for s in body if not isinstance(s, (ast.Assign, ast.AnnAssign)))
        worst_exec = max(worst_exec, n_exec)
    for p in files:
        src = p.read_text()
        tree = ast.parse(src)
        hist += len(re.findall(HISTORY, src, re.I))
        for f in [n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef)]:
            doc = ast.get_docstring(f)
            n_doc = len(doc.splitlines()) + 2 if doc else 0
            longest = max(longest, f.end_lineno - f.lineno + 1 - n_doc)
            names = {n.id for n in ast.walk(f) if isinstance(n, ast.Name)}
            unused += [f"{p.name}:{f.name}:{a.arg}" for a in f.args.args
                       if a.arg not in names and a.arg not in ("self", "cls")]
            if not f.name.startswith("_"):
                public += 1
                args = [a for a in f.args.args if a.arg not in ("self", "cls")]
                if f.returns is not None and all(a.annotation is not None for a in args) and doc:
                    typed += 1
    rows["module_level_executable_lines_worst_load_model_file"] = worst_exec
    rows["longest_package_function_lines_excluding_docstring"] = longest
    rows["unused_parameters_in_package"] = len(unused)
    rows["public_package_functions"] = public
    rows["public_package_functions_typed_and_documented_pct"] = round(100.0*typed/max(public, 1), 1)
    rows["revision_history_phrases_in_package"] = hist
    rule = ["F", "E7", "E9", "B", "C901", "PLR0912", "PLR0915"]
    py = subprocess.check_output(["git", "ls-files", "*.py"], text=True, cwd=ROOT).split()
    r = subprocess.run(["ruff", "check", "--select", ",".join(rule), "--output-format", "concise", *py],
                       capture_output=True, text=True, cwd=ROOT)
    rows["lint_findings_whole_repository"] = len([line for line in r.stdout.splitlines() if re.match(r"\S+:\d+:\d+:", line)])
    # checks passed over without a word: an except clause whose whole body is "pass"
    rows["silently_skipped_checks"] = sum(
        len(re.findall(r"except[^\n]*:\s*\n\s*pass\b", (ROOT/f).read_text())) for f in py)
    df = pd.DataFrame({"measure": list(rows), "value": list(rows.values())})
    df.to_csv(RESULTS/"code_quality.csv", index=False)
    if unused:
        print("unused:", unused)
    return df


if __name__ == "__main__":
    if "--snapshot" in sys.argv:
        np.savez_compressed(REFERENCE, **_outputs())
        print(f"[quality] reference written: {REFERENCE.name}")
    else:
        print(measure().to_string(index=False))
