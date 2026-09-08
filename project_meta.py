"""
project_meta.py — single source of truth for the case-study identity.

Author, title and the study date of record are stated ONCE here and imported by
every deliverable builder (report docx, report PDFs, engineering drawings), so
the README, the drawing title blocks and the report covers cannot disagree.
Each of those previously carried its own hardcoded copy.

Note the distinction: STUDY_DATE is the date of record for the study and is
stable across rebuilds. The date a file was generated belongs in that file's
metadata, not on its cover — a cover date that moves every rebuild is worse
than a stale one.
"""
AUTHOR       = "Akosa Samuel Onyejekwe"
AUTHOR_SUFFIX = "(independent)"
AUTHOR_BAND  = "AKOSA SAMUEL ONYEJEKWE"          # drawing title blocks (upper case)
AUTHOR_FULL  = f"{AUTHOR} {AUTHOR_SUFFIX}"

TITLE  = "Prediction of Dynamic Stall on a Helicopter Main-Rotor Retreating Blade"
SOLVER = "UNISTALL™ Universal Unsteady-Aerodynamics & Dynamic-Stall Solver"
METHOD = "Unified Indicial–Beddoes State-Space (UIBS)"

STUDY_DATE_ISO = "2026-06-27"        # ISO form, used in drawing title blocks
STUDY_DATE     = "27 June 2026"      # long form, used on report covers
