# -*- coding: utf-8 -*-
"""
project_meta.py
---------------
The identity of the case study, stated once and imported by every builder
(report, drawings, plots), so covers, title blocks and the README agree; the
definition of the two cases; and the one reader of 03_model_setup, through
which every later stage takes the conditions it runs at.

Author: Akosa Samuel Onyejekwe (independent)
"""
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SETUP = ROOT/"03_model_setup"

AUTHOR = "Akosa Samuel Onyejekwe"
AUTHOR_SUFFIX = "(independent)"
AUTHOR_BAND = "AKOSA SAMUEL ONYEJEKWE"          # drawing title blocks
AUTHOR_FULL = f"{AUTHOR} {AUTHOR_SUFFIX}"
AFFILIATION = "Independent Researcher"            # where a form or a title page asks for an affiliation
# One title for each thing that carries one: the software (README, CITATION.cff),
# and the case study (TITLE below: report, site, drawings). The manuscript has its own.
SOFTWARE_TITLE = "UNISTALL: a dynamic-stall load model for a pitching aerofoil"
SOFTWARE_SUBTITLE = "A Leishman–Beddoes model with tabulated separation, assessed against NASA TM-84245"

TITLE = "Prediction of Dynamic Stall on a Pitching NACA 0012 Section"
SUBTITLE = "A case study at a tunnel condition and at a rotor-blade-station condition"
SOLVER = "UNISTALL"                              # the name of the load model in this case study
SOLVER_LONG = "UNISTALL, a Leishman-Beddoes dynamic-stall load model with tabulated separation (Python package unistall)"
METHOD = "Leishman-Beddoes indicial method with a tabulated separation point"

STUDY_DATE_ISO = "2026-10-07"
STUDY_DATE = "7 October 2026"

SECTION = "0012"                                 # NACA designation of the section of both cases
N_SURFACE_POINTS = 160                           # cosine-spaced stations on each surface of the coordinate file

# Sea-level standard air: ratio of specific heats, gas constant [J/kg/K],
# temperature [K], pressure [Pa]. The speed of sound follows from the first
# three (unistall.flowfield.Air.a_sound).
GAMMA, R_GAS, T_INF, P_INF = 1.4, 287.05, 288.15, 101325.0
FT_TO_M, IN_TO_M = 0.3048, 0.0254

# The rotor whose blade station gives Case B, as the source states it.
ROTOR_SOURCE = "Bousman, NASA/TP-2003-212265 (UH-60A main rotor)"
ROTOR = dict(radius_ft=26.83, chord_in=20.76, blades=4, speed_rpm=258)
ROOT_CUTOUT_IN, TIP_RADIUS_IN = 42.00, 322.00    # radial locations of the root cutout and the tip, Table 3 of the source
STATION_OVER_R = 0.75                            # the blade station of Case B
STATION_AZIMUTH_DEG = 270.0                      # retreating side, where the station sees its lowest speed
STATION_MACH = 0.28                              # Mach number of Case B: the station's at that azimuth
ROOT_END_OVER_R = round(ROOT_CUTOUT_IN/TIP_RADIUS_IN, 2)    # where the drawings start the blade

CASE_A, CASE_B = "A_validation", "B_application"


def speed_of_sound() -> float:
    """Speed of sound [m/s] of the standard air above, sqrt(gamma R T)."""
    return math.sqrt(GAMMA*R_GAS*T_INF)


def station_condition() -> dict:
    """The condition of Case B, derived from the rotor and the station.

    The rotor turns at Omega = 2 pi rpm / 60. A station at radius r sees the
    speed Omega r + V sin(psi) at azimuth psi in forward flight at speed V.
    Case B is defined by the Mach number the station sees at
    STATION_AZIMUTH_DEG; the forward speed that gives it follows, and so does
    the range of Mach number round the azimuth. The pitch of Case B is once
    per revolution, so its circular frequency is Omega and its reduced
    frequency is Omega c / (2 U). Lengths in metres, speeds in m/s, Omega in
    rad/s, frequency in Hz."""
    a = speed_of_sound()
    radius, chord = ROTOR["radius_ft"]*FT_TO_M, round(ROTOR["chord_in"]*IN_TO_M, 3)
    omega = 2.0*math.pi*ROTOR["speed_rpm"]/60.0
    rotational = omega*STATION_OVER_R*radius
    U = STATION_MACH*a
    forward = (U - rotational)/math.sin(math.radians(STATION_AZIMUTH_DEG))
    return dict(radius_m=radius, chord_m=chord, omega_rad_s=omega, rev_per_s=ROTOR["speed_rpm"]/60.0,
                station_radius_m=STATION_OVER_R*radius, rotational_speed_ms=rotational,
                tip_speed_ms=omega*radius, speed_of_sound_ms=a, station_speed_ms=U,
                forward_speed_ms=forward, advance_ratio=forward/(omega*radius),
                mach_min=(rotational - abs(forward))/a, mach_max=(rotational + abs(forward))/a,
                reduced_frequency=omega*chord/(2.0*U))


_STATION = station_condition()

# The two cases. A is a condition of the NASA TM-84245 experiment; B is the
# condition of an aerofoil at a rotor-blade station, used as an illustration.
# B is not a rotor calculation: steady stream, prescribed pitch. Its chord
# and reduced frequency follow from the rotor above (station_condition).
CASES = {
    CASE_A: dict(label="Case A: tunnel condition", alpha_mean_deg=10.0, alpha_amp_deg=10.0, k=0.10, M=0.30,
                 chord_m=0.61, measured_frame="frame_9302"),
    CASE_B: dict(label="Case B: rotor-blade-station condition", alpha_mean_deg=12.0, alpha_amp_deg=8.0,
                 k=_STATION["reduced_frequency"], M=STATION_MACH, chord_m=_STATION["chord_m"], measured_frame=None),
}


def short_label(name: str) -> str:
    """'Case A' or 'Case B' for a case key."""
    return CASES[name]["label"].split(":")[0]


def read_setup() -> tuple:
    """(cases, air, config): the conditions as 03_model_setup holds them.

    Every stage after the setup takes its numbers from here. Each case is the
    entry of solver_config.json with these added: a_mean, a_amp [deg], k, M,
    chord [m], U [m/s], Re, omega [rad/s], steps and cycles of the march, Tvl
    (semichords) and a heading line for figures. `air` is the
    unistall.flowfield.Air of air_properties.csv."""
    import pandas as pd

    from unistall import flowfield

    with open(SETUP/"solver_config.json", encoding="utf-8") as fh:
        cfg = json.load(fh)
    flow = pd.read_csv(SETUP/"flow_conditions.csv").set_index("parameter")
    kin = pd.read_csv(SETUP/"kinematics.csv").set_index("case")
    prop = pd.read_csv(SETUP/"air_properties.csv").set_index("property")["value"]
    air = flowfield.Air(**{f: float(prop[f]) for f in ("gamma", "R_gas", "prandtl", "T_inf", "p_inf")})
    cases = {}
    for name, c in cfg["cases"].items():
        f, kr = flow[name], kin.loc[name]
        line = "%s  (NACA %s, M = %.2f, k = %.3f, α = %.0f° ± %.0f°)" % (
            c["label"], cfg["section"], float(f["freestream_mach_M"]), float(kr["reduced_freq_k"]),
            float(kr["alpha_mean_deg"]), float(kr["alpha_amp_deg"]))
        if not c["measured_frame"]:
            line += "\n" + c["description"]
        cases[name] = dict(c, a_mean=float(kr["alpha_mean_deg"]), a_amp=float(kr["alpha_amp_deg"]),
                           k=float(kr["reduced_freq_k"]), M=float(f["freestream_mach_M"]),
                           chord=float(f["chord_c"]), U=float(f["freestream_velocity_U"]),
                           Re=float(f["reynolds_number_Re_c"]), omega=float(kr["omega_rad_s"]),
                           steps=int(c["march"]["steps_per_cycle"]), cycles=int(c["march"]["cycles"]),
                           Tvl=float(cfg["constants"]["Tvl"]), solver=cfg["solver_name"], line=line)
    return cases, air, cfg


def experiment_citation() -> str:
    """The two volumes of the experiment report as one citation string,
    authors in the order docs/references.bib gives them."""
    text = (ROOT/"docs"/"references.bib").read_text(encoding="utf-8")
    parts = []
    for key, volume in (("mccroskey1982v1", "Vol. 1"), ("mcalister1982v2", "Vol. 2")):
        entry = text[text.index("{" + key + ","):]
        authors = entry[entry.index("author = {") + 10:entry.index("},")]
        names = [a.split(",")[0].strip() for a in authors.split(" and ")]
        names = [n[:2] + n[2].upper() + n[3:] if n.startswith("Mc") else n for n in names]
        parts.append(", ".join(names[:-1]) + " & " + names[-1] + f" (1982), NASA TM-84245 {volume}")
    return "; ".join(parts)


# Shared by the figure scripts.
CP_TICKS = (-30.0, -20.0, -10.0, -5.0, -2.0, -1.0, -0.5, 0.0, 0.5, 1.0)   # marks of a colour bar of Cp
CMAP_TRIM = (0.07, 0.95)                         # part of a colour map used, leaving out its darkest ends
NOTE_PT = 8.5                                    # lettering of the notes under a figure


def trimmed_colormap(name: str) -> object:
    """The named matplotlib colour map without its darkest ends."""
    import numpy as np
    from matplotlib import colormaps
    from matplotlib.colors import LinearSegmentedColormap

    return LinearSegmentedColormap.from_list(name + "_trimmed", colormaps[name](np.linspace(*CMAP_TRIM, 256)))


def cp_colour_bar(bar: object, lo: float, hi: float) -> None:
    """Set a colour bar of the pressure coefficient the way every figure
    shows it: plain marks at CP_TICKS within the data, suction upwards."""
    ticks = [t for t in CP_TICKS if lo <= t <= hi]
    bar.set_ticks(ticks, labels=["%g" % t for t in ticks])
    bar.ax.minorticks_off()
    bar.ax.invert_yaxis()


def smallest_lettering(fig: object) -> float:
    """Size [pt] of the smallest visible lettering of a matplotlib figure."""
    from matplotlib.text import Text

    sizes = [t.get_fontsize() for t in fig.findobj(Text) if t.get_visible() and t.get_text().strip()]
    return float(min(sizes)) if sizes else float("nan")


def record_figures(record_csv: Path, rows: list) -> None:
    """Enter figures in a figure record: `rows` of (file name, dots per
    inch it was saved at, smallest lettering [pt]). Rows already in the
    record for other files are kept, so several scripts can share one
    record. The report uses it to work out the size of the lettering as
    printed."""
    import pandas as pd

    new = pd.DataFrame(rows, columns=["file", "dpi", "smallest_lettering_pt"])
    if Path(record_csv).exists():
        old = pd.read_csv(record_csv)
        new = pd.concat([old[~old["file"].isin(new["file"])], new], ignore_index=True)
    new.sort_values("file").to_csv(record_csv, index=False)
