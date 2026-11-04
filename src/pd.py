"""Pharmacodynamics: from concentration to effect.

Sigmoid Emax (Hill) model:

    E = E0 + Emax * Ce^n / (EC50^n + Ce^n)

EC50 is the concentration giving half the maximum effect; the Hill coefficient n
sets how steep the switch is. Like a dimmer that saturates: past a point, more
drug barely brightens the room.

Effects often lag behind plasma levels because the drug acts somewhere else
(brain, muscle) that takes time to equilibrate. An effect compartment models
this with no mass of its own:

    dCe/dt = ke0 * (C - Ce)

so Ce follows C with a delay of half-life ln(2) / ke0. Plotting effect against
plasma concentration then traces a loop (hysteresis) rather than a line.
"""
import numpy as np


def effect_site(t, c, ke0):
    """Effect-site concentration Ce(t) for plasma concentration c sampled at times t.

    Between samples C is treated as a straight line, for which the equation above
    has an exact solution, so Ce is exact for piecewise-linear C:

        Ce(t+h) = C(t+h) - s/ke0 + (Ce(t) - C(t) + s/ke0) * exp(-ke0 h),  s = slope of C
    """
    t, c = np.asarray(t, dtype=float), np.asarray(c, dtype=float)
    if np.isinf(ke0):
        return c.copy()
    ce = np.zeros_like(c)
    for i in range(1, t.size):
        h = t[i] - t[i - 1]
        slope = (c[i] - c[i - 1]) / h
        decay = np.exp(-ke0 * h)
        ce[i] = c[i] - slope / ke0 + (ce[i - 1] - c[i - 1] + slope / ke0) * decay
    return ce


def emax_effect(ce, e0, emax, ec50, hill=1.0):
    """Sigmoid Emax model. emax < 0 for drugs that lower something (e.g. heart rate)."""
    ce = np.clip(np.asarray(ce, dtype=float), 0, None)
    return e0 + emax * ce**hill / (ec50**hill + ce**hill)
