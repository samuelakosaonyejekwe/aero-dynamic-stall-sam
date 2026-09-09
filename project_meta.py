"""
project_meta.py — single source of truth for the case-study identity.

Author, title and the study date of record are stated ONCE here and imported by
every deliverable builder (report docx, report PDFs, engineering drawings), so
the README, the drawing title blocks and the report covers cannot disagree.
Each of those previously carried its own hardcoded copy.

Note the distinction: STUDY_DATE is the date of record for the REVISION, not a
build timestamp. It is stable across rebuilds — a cover date that moved every
time the pipeline ran would be worse than a stale one, and the date a file was
generated belongs in that file's metadata rather than on its cover. It is
advanced deliberately, by editing this file, whenever the study itself is
revised; it had been left at the date of the first release while the solver,
the validation and the report were all revised past it, which put a date on
every cover and every drawing title block that the content no longer matched.
check_claims.py asserts that README.md still quotes the value set here.
"""
AUTHOR       = "Akosa Samuel Onyejekwe"
AUTHOR_SUFFIX = "(independent)"
AUTHOR_BAND  = "AKOSA SAMUEL ONYEJEKWE"          # drawing title blocks (upper case)
AUTHOR_FULL  = f"{AUTHOR} {AUTHOR_SUFFIX}"

TITLE  = "Prediction of Dynamic Stall on a Helicopter Main-Rotor Retreating Blade"
SOLVER = "UNISTALL™ Universal Unsteady-Aerodynamics & Dynamic-Stall Solver"
METHOD = "Unified Indicial–Beddoes State-Space (UIBS)"   # used by 03_model_setup + 07_report

STUDY_DATE_ISO = "2026-09-09"        # ISO form, used in drawing title blocks
STUDY_DATE     = "9 September 2026"  # long form, used on report covers
