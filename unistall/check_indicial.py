# -*- coding: utf-8 -*-
"""
unistall / check_indicial.py
------------------------
The compressible attached-flow branch against exact linear theory at the Mach
numbers of the study, where the comparison with Theodorsen's incompressible
solution does not reach.

Author: Akosa Samuel Onyejekwe (independent)

Two exact results of linearised subsonic theory are used (Lomax et al., NACA
TR 1077; piston theory for the starting values):

  STARTING VALUES  just after a step, per unit step:
                   C_N,alpha = 4/M      C_M,alpha(c/4) = -1/M
                   C_N,q     = 1/M      C_M,q(c/4)     = -7/(12 M)
  SHORT-TIME LIFT  for a step in incidence without pitch rate and
                   0 <= s <= 2M/(1+M) semichords:
                   C_N,alpha(s) = (4/M) [1 - (1 - M) s / (2M)]

The loads are marched by unistall/attached_flow.py with a small step; a step
in incidence is one step long, a step in pitch rate is a ramp in incidence,
and the response is read from the components the march returns. No
oscillatory compressible solution is used: at these Mach numbers the
frequency response of the model is checked only through these two results
and the measured loops.

Output: results/indicial_checks.csv
Usage:  python3 -m unistall.check_indicial
"""
import numpy as np
import pandas as pd

from unistall import attached_flow as af
from unistall.dsmodel import A_SOUND as dm_A_SOUND
from unistall.paths import RESULTS

MACH = (0.22, 0.30)
CHORD, A_SOUND = 0.61, dm_A_SOUND
DS = 1.0e-4                    # step of the march, semichords
STEP = 1.0e-3                  # size of the step in incidence [rad] and in pitch rate [-]


def _march(M: float, alpha: np.ndarray, dt: float, x_pitch: float = 0.25) -> list:
    m = af.AttachedFlow(M*A_SOUND, CHORD, M, 2.0*np.pi/np.sqrt(1.0 - M*M), dt, x_pitch=x_pitch, low_mach_blend=False)
    return [m.start(alpha[0], alpha[1])] + [m.step(a) for a in alpha[1:]]


def starting_values(M: float) -> list:
    """The four starting values at Mach number M: the model against theory."""
    U = M*A_SOUND
    dt = DS*CHORD/(2.0*U)
    n = 6
    a_step = np.concatenate([np.zeros(3), np.full(n, STEP)])                 # step in incidence
    r = _march(M, a_step, dt)[4]                                              # the first instant after the step
    a_ramp = np.concatenate([np.zeros(3), STEP*(U/CHORD)*dt*np.arange(1, n + 1)])   # q steps from 0 to STEP
    q = _march(M, a_ramp, dt)[4]
    rows = [("C_N per unit step in incidence", r["CN_imp_alpha"]/STEP, 4.0/M),
            ("C_M about c/4 per unit step in incidence", r["CM_imp_alpha"]/STEP, -1.0/M),
            ("C_N per unit step in pitch rate", q["CN_imp_q"]/STEP, 1.0/M),
            ("C_M about c/4 per unit step in pitch rate", (q["CM_imp_q"] + q["CM_circ_q"])/STEP, -7.0/(12.0*M))]
    return [dict(check="starting value", M=M, quantity=name, model=float(got), theory=float(exact),
                 error_pct=float(100.0*(got/exact - 1.0))) for name, got, exact in rows]


def short_time_lift(M: float) -> dict:
    """Largest difference between the model's lift after a step in incidence
    and the exact short-time result, over 0 <= s <= 2M/(1+M)."""
    U = M*A_SOUND
    dt = DS*CHORD/(2.0*U)
    s_end = 2.0*M/(1.0 + M)
    n = int(np.ceil(s_end/DS)) + 2
    # pitch axis at three-quarter chord, so that the step in incidence is seen
    # by the circulatory lags as incidence alone; the pitch-rate impulse of the
    # step is left out by summing the incidence terms only
    out = _march(M, np.concatenate([np.zeros(3), np.full(n, STEP)]), dt, x_pitch=0.75)
    s = DS*(np.arange(len(out)) - 3.0)                    # travel since the end of the step
    keep = (s >= DS) & (s <= s_end)
    model = np.array([o["CN_imp_alpha"] + o["CN_circ"] for o in out])[keep]/STEP
    exact = (4.0/M)*(1.0 - (1.0 - M)*s[keep]/(2.0*M))
    i = int(np.argmax(np.abs(model - exact)))
    return dict(check="short-time lift after a step in incidence", M=M,
                quantity=f"largest difference over 0 <= s <= {s_end:.3f} semichords (at s = {s[keep][i]:.3f})",
                model=float(model[i]), theory=float(exact[i]), error_pct=float(100.0*(model[i]/exact[i] - 1.0)))


def main() -> pd.DataFrame:
    rows = []
    for M in MACH:
        rows += starting_values(M) + [short_time_lift(M)]
    df = pd.DataFrame(rows)
    df.round(5).to_csv(RESULTS/"indicial_checks.csv", index=False)
    pd.set_option("display.width", 200)
    print(df.round(4).to_string(index=False))
    return df


if __name__ == "__main__":
    main()
