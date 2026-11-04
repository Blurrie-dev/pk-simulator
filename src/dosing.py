"""Dosing tools: find a regimen that hits a therapeutic window, and a loading dose.

Both lean on linearity. In a linear model doubling the dose doubles every
concentration, so the steady-state curve for any dose is the curve for 1 mg,
scaled. One simulation per dosing interval is enough to price every dose.
"""
from dataclasses import replace

import numpy as np

from src.metrics import steady_state_peak_trough
from src.simulate import MAX_DOSES


def steady_state_per_mg(scenario, tau):
    """Steady-state (peak, trough) in mg/L per mg of maintenance dose, at interval tau."""
    s = replace(scenario.regular, dose=1.0, tau=tau)
    if s.infusion:
        s = replace(s, duration=min(s.duration, tau))  # an infusion can't outlast its interval
    # 14 half-lives instead of 7: within 0.01% of steady state, not 1%
    n = min(2 * s.doses_needed(), MAX_DOSES)
    t = s.grid(n)
    return steady_state_peak_trough(t, s.repeated(t, n), tau)


def suggest_regimens(scenario, window, taus=(6, 8, 12, 24), dose_step=25.0):
    """For each interval, the dose that best centres steady state in the window.

    "Centred" means the geometric mean of peak and trough equals the geometric
    mean of the window's limits, so peak and trough sit equally far (as ratios)
    inside it. A regimen fits only if the curve's swing, peak / trough, is smaller
    than the window's, upper / lower: no dose can fit a wide swing into a narrow
    window, whatever its size. That's what shorter intervals or infusions fix.

    Returns one dict per interval: tau, dose (rounded to dose_step), daily dose,
    peak, trough and whether it fits.
    """
    low, high = window
    if scenario.saturable:
        raise ValueError("saturable elimination isn't linear: scale the dose by hand instead")
    if not 0 < low < high:
        raise ValueError("needs a window with 0 < lower < upper")
    rows = []
    for tau in taus:
        peak, trough = steady_state_per_mg(scenario, tau)
        ideal = np.sqrt(low * high / (peak * trough))
        dose = max(dose_step, round(ideal / dose_step) * dose_step)
        rows.append(dict(tau=tau, dose=dose, daily_dose=dose * 24 / tau,
                         peak=dose * peak, trough=dose * trough,
                         fits=bool(low <= dose * trough and dose * peak <= high)))
    return rows


def loading_dose(scenario):
    """Dose that makes the first peak equal the steady-state peak (mg).

    Repeated doses pile up until the peak is R times the first one, R being the
    accumulation ratio. Giving R times the dose first skips the wait, like filling
    the bath with a bucket before turning the tap down to a trickle. For an IV
    bolus with one compartment this is exactly dose / (1 - exp(-k * tau)).
    """
    s = scenario.regular
    peak_ss, _trough = steady_state_per_mg(s, s.tau)
    t = s.grid(1)
    first_peak = replace(s, dose=1.0).single(t).max()
    return s.dose * peak_ss / first_peak
