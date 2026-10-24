"""Run with `pytest` from the repo root."""
import numpy as np
import pytest

from src.metrics import auc, cmax, half_life, steady_state_peak_trough, time_to_steady_state, tmax
from src.models import (multiple_doses, one_compartment_iv, one_compartment_ode,
                        one_compartment_oral, time_grid)
from src.presets import PRESETS, WEIGHT

# Hypothetical drug: k = CL / Vd = 0.1 /h, half-life = ln(2) / k = 6.9 h
DOSE, CL, VD = 500, 5, 50
T = np.linspace(0, 72, 500)  # about 10 half-lives


def test_iv_numerical_matches_analytical():
    numerical = one_compartment_ode(T, DOSE, CL, VD)
    analytical = one_compartment_iv(T, DOSE, CL, VD)
    np.testing.assert_allclose(numerical, analytical, rtol=0.01)  # every point within 1%


@pytest.mark.parametrize("ka", [1.5, 0.1, 0.03])  # fast absorption, ka == k, flip-flop (ka < k)
def test_oral_numerical_matches_analytical(ka):
    numerical = one_compartment_ode(T, DOSE, CL, VD, ka=ka, f=0.7)
    analytical = one_compartment_oral(T, DOSE, CL, VD, ka=ka, f=0.7)
    np.testing.assert_allclose(numerical, analytical, rtol=0.01)


def test_auc_and_half_life_match_theory():
    t = np.linspace(0, 200, 4001)  # 29 half-lives: the curve is back to ~zero
    c = one_compartment_oral(t, DOSE, CL, VD, ka=1.5, f=0.7)
    assert auc(t, c) == pytest.approx(0.7 * DOSE / CL, rel=0.01)  # AUC = F * dose / CL
    assert half_life(t, c) == pytest.approx(np.log(2) * VD / CL, rel=0.01)


@pytest.mark.parametrize("dose", [50, 500, 5000])
def test_steady_state_after_4_to_5_half_lives(dose):
    cl, vd, tau = 2, 100, 12  # half-life 34.7 h, dosed every 12 h
    t_half = np.log(2) * vd / cl
    n_doses = 45              # 540 h = 15.6 half-lives, so the last trough is steady state
    t = time_grid(tau, n_doses)
    c = multiple_doses(one_compartment_oral, t, tau, n_doses, dose=dose, cl=cl, vd=vd, ka=1.0)
    assert 4 <= time_to_steady_state(t, c, tau) / t_half <= 5


# The preset tests check that parameters taken from the sources, scaled to 70 kg,
# reproduce what the same sources observed. Reference numbers [n] are in the README.

def test_theophylline_preset_reproduces_label():
    # [1]: a single 5 mg/kg dose peaks at about 10 mcg/mL (range 5-15) 1-2 h after
    # dosing; half-life in healthy non-smoking adults is 6.1-12.8 h.
    p = PRESETS["Theophylline (oral)"]
    t = np.linspace(0, 24, 2401)
    c = one_compartment_oral(t, 5 * WEIGHT, p["cl"], p["vd"], p["ka"], p["f"])
    assert 5 <= cmax(c) <= 15
    assert 1 <= tmax(t, c) <= 2
    assert 6.1 <= np.log(2) * p["vd"] / p["cl"] <= 12.8


def test_gentamicin_preset_reproduces_label():
    # [2]: at 3 mg/kg/day every 8 h, peaks are expected at 4-6 mcg/mL and troughs
    # should stay below 2 mcg/mL.
    p = PRESETS["Gentamicin (IV)"]
    t = time_grid(p["tau"], 10)
    c = multiple_doses(one_compartment_iv, t, p["tau"], 10, dose=p["dose"], cl=p["cl"], vd=p["vd"])
    peak, trough = steady_state_peak_trough(t, c, p["tau"])
    assert 4 <= peak <= 6
    assert trough < 2
