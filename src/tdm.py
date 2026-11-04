"""Therapeutic drug monitoring: estimate one patient's parameters from measured levels.

Population values say what a typical patient looks like; a couple of blood levels
say what *this* patient looks like. Bayesian (MAP) estimation combines them:
find the parameters that make the measured levels likely, without straying
further from the typical values than the levels justify.

Like a satnav fusing a rough map with a few GPS fixes: with no fixes you trust
the map, with many you trust the fixes.

Objective minimised (-2 log posterior, up to a constant):

    sum over levels  [ ln(var_j) + (obs_j - pred_j)^2 / var_j ]
  + sum over params  eta_p^2 / omega_p^2

where each parameter is typical * exp(eta), omega comes from its population CV
(see src/population.py), and var_j = sd_add^2 + (cv_prop * pred_j)^2 is the
assay plus model error.
"""
from dataclasses import dataclass, replace

import numpy as np
from scipy.optimize import minimize

from src.population import omega


@dataclass
class Fit:
    scenario: object    # Scenario with the estimated parameters
    factors: dict       # estimated / prior value, per fitted parameter
    objective: float
    predicted: np.ndarray  # fitted concentration at each observation time


def doses_before(scenario, last_time):
    """Number of scheduled doses up to the last observation."""
    return int(last_time // scenario.tau) + 1


def predict(scenario, times):
    """Concentrations at the observation times, under the scenario's dose schedule."""
    times = np.asarray(times, dtype=float)
    order = np.argsort(times)
    c = np.empty_like(times)
    c[order] = scenario.concentration(times[order], scenario.schedule(doses_before(scenario, times.max())))
    return c


def fit(scenario, times, levels, cv, parameters=("cl", "vd"), sd_add=0.1, cv_prop=0.15):
    """MAP estimate of `parameters` given measured levels (mg/L) at times (h after the first dose).

    cv: population variability per parameter (fraction), the prior's width.
    """
    times, levels = np.asarray(times, dtype=float), np.asarray(levels, dtype=float)
    parameters = [p for p in parameters if cv.get(p, 0) > 0]
    if times.size == 0 or not parameters:
        return Fit(scenario, {p: 1.0 for p in parameters}, 0.0, predict(scenario, times) if times.size else times)
    omegas = np.array([omega(cv[p]) for p in parameters])

    def with_etas(etas):
        return replace(scenario, **{p: getattr(scenario, p) * np.exp(e) for p, e in zip(parameters, etas)})

    def objective(etas):
        pred = predict(with_etas(etas), times)
        var = sd_add**2 + (cv_prop * pred) ** 2
        return float(np.sum(np.log(var) + (levels - pred) ** 2 / var) + np.sum((etas / omegas) ** 2))

    best = minimize(objective, np.zeros(len(parameters)), method="Nelder-Mead",
                    options=dict(xatol=1e-6, fatol=1e-8, maxiter=4000))
    fitted = with_etas(best.x)
    return Fit(fitted, {p: float(np.exp(e)) for p, e in zip(parameters, best.x)}, float(best.fun),
               predict(fitted, times))
