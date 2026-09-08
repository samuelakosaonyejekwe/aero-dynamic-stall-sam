"""
run_all.py — reproduce the entire case study end-to-end.
Author: Akosa Samuel Onyejekwe (independent).
"""
import subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
# EXECUTION ORDER != folder order, deliberately. 03_model_setup is the single
# source of truth for the case conditions, and 01_geometry (chord) and 02_mesh
# (chord, rho, U, mu) both consume them, so the setup stage has to run first.
# They used to carry their own hardcoded copies of those values.
STEPS = [
    ("Model setup",     "03_model_setup/generate_setup.py"),
    ("Geometry",        "01_geometry/generate_geometry.py"),
    ("Mesh",            "02_mesh/generate_mesh.py"),
    ("Solver run",      "04_solver/run_case.py"),
    ("2-D postprocess", "06_postprocessing/make_all_plots.py"),
    ("3-D postprocess", "06_postprocessing/make_3d_plots.py"),
    ("Validation (static / calibration)", "06_postprocessing/validation/validate.py"),
    ("Validation (real NASA data)",  "06_postprocessing/validation/validate_nasa_real.py"),
    ("Validation (digitizer harness)", "06_postprocessing/validation/validate_digitized.py"),
    ("Engineering drawings", "08_engineering_drawings/draw_engineering.py"),
    ("Report (docx)",   "07_report/build_docx.py"),
    ("Report (PDFs)",   "07_report/build_pdfs.py"),
    ("Report (consolidated PDF)", "07_report/build_report_pdf.py"),
    # last, so it checks the artifacts this run actually produced
    ("Invariant check (physics & numerics)", "verify_invariants.py"),
    ("Claim check (README vs artifacts)", "check_claims.py"),
]
for name, rel in STEPS:
    script = ROOT/rel
    print(f"\n========== {name}: {rel} ==========")
    r = subprocess.run([sys.executable, script.name], cwd=script.parent)
    if r.returncode != 0:
        print(f"[run_all] FAILED at {name}"); sys.exit(r.returncode)
print("\n[run_all] complete — see aero_dynamic_stall_report.pdf, "
      "07_report/case.docx and the two component PDFs.")
