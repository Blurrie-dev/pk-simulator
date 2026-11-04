"""Variability between patients: simulate many people instead of one.

Two people given the same dose don't get the same curve. Population PK describes
each parameter as a typical value times a random factor:

    CL_i = CL_typical * exp(eta_i),   eta_i ~ Normal(0, omega^2)

so parameters are log-normal: always positive, skewed to the right, and the
typical value is the median. omega is set from the coefficient of variation
(CV) people usually quote: omega^2 = ln(1 + CV^2).

It's like rolling a set of dice per patient for CL, V and ka, then running the
same model for each.
"""
from dataclasses import replace

import numpy as np

from src.metrics import steady_state_peak_trough

PARAMETERS = ("cl", "vd", "ka", "vmax")  # the ones that get a random factor


def omega(cv):
    """Standard deviation of eta for a log-normal parameter with this CV."""
    return np.sqrt(np.log(1 + cv**2))


def sample_patients(scenario, cv, n, seed=0):
    """n scenarios with random parameters around `scenario`'s values.

    cv: {parameter name: coefficient of variation as a fraction, e.g. 0.3 for 30%}
    Parameters are drawn independently. Same seed, same patients.
    """
    rng = np.random.default_rng(seed)
    etas = {name: rng.normal(0, omega(cv.get(name, 0.0)), n) for name in PARAMETERS}
    return [replace(scenario, **{name: getattr(scenario, name) * np.exp(etas[name][i]) for name in PARAMETERS})
            for i in range(n)]


def simulate(patients, t, n_doses, rtol=1e-5):
    """Repeated-dosing curve for every patient: one row per patient.

    rtol: solver accuracy for saturable elimination (numerical). 1e-5 is plenty
    when the spread between patients is tens of percent, and several times faster.
    """
    return np.array([p.repeated(t, n_doses, rtol) for p in patients])


def percentiles(curves, levels=(5, 25, 50, 75, 95)):
    """{level: curve} of concentration percentiles across patients at each time."""
    return dict(zip(levels, np.percentile(curves, levels, axis=0)))


def steady_state_summary(curves, t, tau):
    """Each patient's steady-state peak and trough (last dosing interval)."""
    peaks, troughs = np.array([steady_state_peak_trough(t, c, tau) for c in curves]).T
    return peaks, troughs


def target_attainment(peaks, troughs, window):
    """Fraction of patients whose whole steady-state interval sits in the window."""
    low, high = window
    return float(np.mean((troughs >= low) & (peaks <= high)))
