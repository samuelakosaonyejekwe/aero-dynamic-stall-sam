# -*- coding: utf-8 -*-
"""
run_all.py
----------
Reproduce the case study with one command, from the repository root:

    python3 run_all.py            every stage
    python3 run_all.py 4          from stage 4 on
    python3 run_all.py 2          offline: every stage but the download

Stage 1 downloads the measured loops, which only the overlay of Case A needs.
If it fails (no network), the run goes on: the loops of Case A are then drawn
without the measured points, the figures and the documents say so, and the
check on the overlay fails until the loops have been fetched.

Author: Akosa Samuel Onyejekwe (independent)

The stages use the calibrated constants and the result files already in
results/; the calibration and the comparison with the held-out measurements
are reproduced by the commands in README.md, which take several hours.
"""
import os
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
# (name, command, may fail): a stage that may fail does not stop the run
STAGES = [
    ("measured loops (download, checksum)", [sys.executable, "-m", "unistall.fetch_frames"], True),
    ("geometry", [sys.executable, "01_geometry/generate_geometry.py"], False),
    ("mesh", [sys.executable, "02_mesh/generate_mesh.py"], False),
    ("model setup", [sys.executable, "03_model_setup/generate_setup.py"], False),
    ("solution of the two cases", [sys.executable, "04_solver/run_case.py"], False),
    ("reconstructed field at the mesh nodes", [sys.executable, "02_mesh/field_on_mesh.py"], False),
    ("plots", [sys.executable, "06_postprocessing/make_all_plots.py"], False),
    ("surface and map plots", [sys.executable, "06_postprocessing/make_3d_plots.py"], False),
    ("comparison with measurement, set out", [sys.executable, "06_postprocessing/validation/make_validation.py"],
     False),
    ("engineering drawings", [sys.executable, "08_engineering_drawings/draw_engineering.py"], False),
    ("report documents", [sys.executable, "07_report/build_report.py"], False),
    ("list of the folders' contents", [sys.executable, "00_overview/write_folder_contents.py"], False),
    ("checks on the case study", [sys.executable, "check_case_study.py"], False),
]
NO_LOOPS = ("[run_all] stage 1 failed: the measured loops could not be fetched. The run goes on WITHOUT them: the "
            "loops of Case A will be drawn without the measured points and the documents will say so. Fetch them "
            "later with 'python3 -m unistall.fetch_frames' and run 'python3 run_all.py 7'. The stage that sets out the "
            "comparison with measurement needs them and stops the run.")


def main(first: int = 1) -> int:
    """Run the stages from number `first` (1-based); stop at the first
    failure of a stage that the later ones need."""
    env = dict(os.environ, PYTHONPATH=str(ROOT) + os.pathsep + os.environ.get("PYTHONPATH", ""))
    for n, (name, cmd, may_fail) in enumerate(STAGES, start=1):
        if n < first:
            continue
        print(f"\n========== stage {n} of {len(STAGES)}: {name} ==========", flush=True)
        t0 = time.time()
        if subprocess.run(cmd, cwd=ROOT, env=env).returncode != 0:
            if not may_fail:
                print(f"stage {n} ({name}) failed")
                return 1
            print(NO_LOOPS, flush=True)
            continue
        print(f"[run_all] stage {n} done in {time.time() - t0:.0f} s", flush=True)
    print("\n[run_all] all stages complete")
    return 0


if __name__ == "__main__":
    sys.exit(main(int(sys.argv[1]) if len(sys.argv) > 1 else 1))
