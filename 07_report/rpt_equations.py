# -*- coding: utf-8 -*-
"""
07_report / rpt_equations.py
----------------------------
The equations shown in the report, typeset with matplotlib's mathtext and
written to _equations/eq_NNN.png. Each equation names the numbered lines of
docs/formulation.md it is taken from; the build stops if a line it names is
not in that file.

Author: Akosa Samuel Onyejekwe (independent)

Units: image sizes are returned in inches at the size they are to be placed.
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from PIL import Image  # noqa: E402

FONT_PT = 12.0           # size of the mathematics on the page
DPI = 300                # resolution of the images
MAX_WIDTH_IN = 5.7       # widest an equation may be placed
NAVY = "#1F3350"

HALF = r"\left(\dfrac{1+\sqrt{f^{\prime\prime}}}{2}\right)^{2}"

# key -> (lines of mathtext, numbered lines of docs/formulation.md), in the order the report uses them
EQUATIONS = {
    "motion": ([r"\alpha(t) = \alpha_m + \alpha_1 \sin(\omega t), \qquad k = \dfrac{\omega c}{2U}"], ()),
    "step": ([r"\Delta s = \dfrac{2U\,\Delta t}{c}, \qquad \beta^2 = 1 - M^2"], ("1", "2")),
    "rate": ([r"q_n = \dfrac{c}{U}\,\dfrac{\alpha_n-\alpha_{n-1}}{\Delta t}, \qquad "
              r"\alpha_{3/4,n} = \alpha_n + q_n\left(\dfrac{3}{4} - x_p\right)"], ("3", "4", "5")),
    "deficiency": ([r"X_{j,n} = X_{j,n-1}\,e^{-b_j\beta^2\Delta s} + A_j\,e^{-b_j\beta^2\Delta s/2}"
                    r"\left(\alpha_{3/4,n}-\alpha_{3/4,n-1}\right), \qquad j = 1,\,2"], ("6", "7")),
    "circulatory": ([r"\alpha_{e,n} = \alpha_{3/4,n} - \alpha_0 - X_{1,n} - X_{2,n}, \qquad "
                     r"C_N^{c} = C_{N\alpha}\,\alpha_e"], ("8", "9")),
    "impulsive": ([r"K^{\prime}_{\alpha,n} = K^{\prime}_{\alpha,n-1}\,e^{-\Delta t/T_\alpha} + "
                   r"\left(K_{\alpha,n}-K_{\alpha,n-1}\right)e^{-\Delta t/(2T_\alpha)}",
                   r"C_N^{nc,\alpha} = \dfrac{4T_\alpha}{M}\left(K_\alpha - K^{\prime}_\alpha\right), \qquad "
                   r"C_N^{pot} = C_N^{c} + C_N^{nc,\alpha} + C_N^{nc,q}"], ("14", "15", "18")),
    "pressure_lag": ([r"D_{p,n} = D_{p,n-1}\,e^{-\Delta s/T_p} + "
                      r"\left(C^{pot}_{N,n}-C^{pot}_{N,n-1}\right)e^{-\Delta s/(2T_p)}, \qquad "
                      r"C^{\prime}_N = C_N^{pot} - D_p"], ("26", "27")),
    "lookup": ([r"\alpha_f = \dfrac{C^{\prime}_N}{C_{N\alpha}} + \alpha_0, \qquad "
                r"f^{\prime} = f(\alpha_f;\,M)"], ("28", "29")),
    "static_f": ([r"f = \left(2\sqrt{\dfrac{C_{N,s}}{C_{N\alpha,s}\,(\alpha-\alpha_{0,s})}} - 1\right)^{2}"],
                 ("50",)),
    "boundary_layer_lag": ([r"D_{f,n} = D_{f,n-1}\,e^{-\Delta s/T_f} + "
                            r"\left(f^{\prime}_n-f^{\prime}_{n-1}\right)e^{-\Delta s/(2T_f)}",
                            r"f^{\prime\prime} = f^{\prime} - D_f, \qquad T_f = \dfrac{T_{f0}}{\sigma_1}"],
                           ("30", "31", "32")),
    "kirchhoff": ([r"C_N^{fs} = C_N^{c}" + HALF + r" + C_N^{nc}"], ("33",)),
    "vortex_feed": ([r"C_V = C_N^{c}\left[1 - " + HALF + r"\right], \qquad T_V = \dfrac{T_{V0}}{\sigma_3}"],
                    ("37", "38")),
    "vortex_lift": ([r"C^{v}_{N,n} = C^{v}_{N,n-1}\,e^{-\Delta s/T_V} + "
                     r"\left(C_{V,n}-C_{V,n-1}\right)e^{-\Delta s/(2T_V)} \qquad (0 < \tau_v \leq T_{VL})",
                     r"C^{v}_{N,n} = C^{v}_{N,n-1}\,e^{-\Delta s/T_V} \quad \mathrm{otherwise}, \qquad "
                     r"C_N = C_N^{fs} + C_N^{v}"], ("39", "40")),
    "chord_force": ([r"C_C = \eta\,C_N^{c}\,\tan\alpha\,\sqrt{f^{\prime\prime}}"], ("43",)),
    "lift_drag": ([r"C_L = C_N\cos\alpha + C_C\sin\alpha, \qquad "
                   r"C_D = C_N\sin\alpha - C_C\cos\alpha + C_{D0}"], ("45", "46")),
    "vortex_arm": ([r"x_v = x_{cp,v}\left[1-\cos\left(\pi\,\min\left(\tau_v/T_{VL},\,1\right)\right)\right]"],
                   ("48",)),
    "moment": ([r"C_M = C_{M,s}(\alpha^{\prime}_f;\,M) + g\left(K_f\,C_M^{c,q} + C_M^{nc,\alpha} + "
                r"C_M^{nc,q}\right) - x_v\,C_N^{v}",
                r"K_f = " + HALF + r", \qquad g = 1 - (1 - g_M)\,w"], ("49",)),
    "step_rule": ([r"N_{steps} \geq \dfrac{2\pi}{k\,\Delta s_{max}}"], ("54",)),
    "damping": ([r"\Xi = -\dfrac{1}{\pi\alpha_1^{2}}\oint C_M\,d\alpha"], ()),
    "ss_circulatory": ([r"\dfrac{dW_i}{ds} = b_i\beta^2\left(\alpha_{3/4} - W_i\right), \qquad "
                        r"X_i = A_i\left(\alpha_{3/4} - W_i\right), \qquad i = 1,\,2",
                        r"\dfrac{dW_5}{ds} = b_5\beta^2\left(q - W_5\right), \qquad "
                        r"C_M^{c,q} = -\dfrac{C_{N\alpha}}{16}\left[q - A_5\left(q - W_5\right)\right]"],
                       ("56", "57")),
    "ss_impulsive": ([r"\dfrac{dZ_a}{ds} = \dfrac{c}{2U}\,\dfrac{K_\alpha - Z_a}{T_\alpha}, \qquad "
                      r"\dfrac{dZ_q}{ds} = \dfrac{c}{2U}\,\dfrac{K_q - Z_q}{T_q}, \qquad "
                      r"C_N^{nc} = \dfrac{4T_\alpha}{M}Z_a + \dfrac{T_q}{M}Z_q"], ("58",)),
    "ss_moment": ([r"\dfrac{dZ_j}{ds} = \dfrac{c}{2U}\,\dfrac{K_\alpha - Z_j}{T_j} \quad (j = 3,\,4), \qquad "
                   r"\dfrac{dZ_m}{ds} = \dfrac{c}{2U}\,\dfrac{K_q - Z_m}{k_{mq}^2 T_I}",
                   r"C_M^{nc} = -\dfrac{A_3T_3}{M}Z_3 - \dfrac{A_4T_4}{M}Z_4 - "
                   r"\dfrac{7k_{mq}^2T_I}{12M}Z_m"], ("59",)),
    "ss_lags": ([r"\dfrac{dC^{\prime}_N}{ds} = \dfrac{C_N^{pot} - C^{\prime}_N}{T_p}, \qquad "
                 r"\dfrac{df^{\prime\prime}}{ds} = \dfrac{f(\alpha_f;\,M) - f^{\prime\prime}}{T_f}, \qquad "
                 r"\dfrac{d\alpha_m}{ds} = \dfrac{\alpha_f - \alpha_m}{c_m T_{f0}}"], ("60", "61", "62")),
    "ss_vortex": ([r"\dfrac{dC_N^{v}}{ds} = -\dfrac{C_N^{v}}{T_V} + \dfrac{d}{ds}\left[C_N^{c}\left(1 - K_f\right)"
                   r"\right] \quad \mathrm{(second\ term\ only\ while\ the\ vortex\ is\ on\ the\ chord)}"],
                  ("63",)),
}


def render_one(lines: list, path: Path) -> tuple:
    """Typeset one equation and return the (width, height) in inches at which
    it is placed: its true size, or less if it is wider than the column."""
    fig = plt.figure(figsize=(0.01, 0.01))
    fig.text(0, 0, "\n".join(f"${t}$" for t in lines), fontsize=FONT_PT, color=NAVY, linespacing=2.0)
    fig.savefig(path, dpi=DPI, bbox_inches="tight", pad_inches=0.03, facecolor="white")
    plt.close(fig)
    with Image.open(path) as im:
        w, h = im.size
    scale = min(1.0, MAX_WIDTH_IN/(w/DPI))
    return w/DPI*scale, h/DPI*scale


def render_all(outdir: Path, source_lines: dict) -> dict:
    """Write every equation to outdir and return key -> dict(path, width_in,
    height_in, number, source, source_text)."""
    outdir.mkdir(exist_ok=True)
    for old in outdir.glob("eq_*.png"):
        old.unlink()
    out = {}
    with plt.rc_context({"mathtext.fontset": "stix"}):
        for number, (key, (lines, src)) in enumerate(EQUATIONS.items(), start=1):
            missing = [s for s in src if s not in source_lines]
            if missing:
                raise KeyError(f"equation '{key}': docs/formulation.md has no line {missing}")
            path = outdir/f"eq_{number:03d}.png"
            w, h = render_one(lines, path)
            out[key] = dict(path=path, width_in=w, height_in=h, number=number, source=list(src),
                            source_text={s: source_lines[s] for s in src})
    return out
