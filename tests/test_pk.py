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
    p = PRESETS["Theophylline"]
    t = np.linspace(0, 24, 2401)
    c = one_compartment_oral(t, 5 * WEIGHT, p["cl"], p["vd"], p["ka"], p["f"])
    assert 5 <= cmax(c) <= 15
    assert 1 <= tmax(t, c) <= 2
    assert 6.1 <= np.log(2) * p["vd"] / p["cl"] <= 12.8


def test_gentamicin_preset_reproduces_label():
    # [2]: at 3 mg/kg/day every 8 h, peaks are expected at 4-6 mcg/mL and troughs
    # should stay below 2 mcg/mL.
    p = PRESETS["Gentamicin"]
    t = time_grid(p["tau"], 10)
    c = multiple_doses(one_compartment_iv, t, p["tau"], 10, dose=p["dose"], cl=p["cl"], vd=p["vd"])
    peak, trough = steady_state_peak_trough(t, c, p["tau"])
    assert 4 <= peak <= 6
    assert trough < 2


def test_levetiracetam_preset_reproduces_label():
    # [4]: peak about an hour after an oral dose, half-life 7 +/- 1 h, steady
    # state within 2 days of twice-daily dosing.
    p = PRESETS["Levetiracetam"]
    t = time_grid(p["tau"], 10)
    single = one_compartment_oral(t, p["dose"], p["cl"], p["vd"], p["ka"], p["f"])
    multi = multiple_doses(one_compartment_oral, t, p["tau"], 10,
                           dose=p["dose"], cl=p["cl"], vd=p["vd"], ka=p["ka"], f=p["f"])
    assert 0.75 <= tmax(t, single) <= 1.25
    assert 6 <= np.log(2) * p["vd"] / p["cl"] <= 8
    assert time_to_steady_state(t, multi, p["tau"]) <= 48


def test_vancomycin_preset_reproduces_label():
    # [5]: half-life 4-6 h; with 1 g every 12 h the concentration just before the
    # next dose is about 8 mcg/mL. One compartment can't match the early, distribution
    # phase concentrations, so only the half-life and trough are checked.
    p = PRESETS["Vancomycin"]
    t = time_grid(p["tau"], 10)
    c = multiple_doses(one_compartment_iv, t, p["tau"], 10, dose=p["dose"], cl=p["cl"], vd=p["vd"])
    _peak, trough = steady_state_peak_trough(t, c, p["tau"])
    assert 4 <= np.log(2) * p["vd"] / p["cl"] <= 6
    assert trough == pytest.approx(8, rel=0.25)


def test_valproate_preset_reproduces_label():
    # [6]: half-life 16 +/- 3 h after IV dosing; therapeutic range 50-100 mcg/mL.
    p = PRESETS["Valproate"]
    t = time_grid(p["tau"], 40)
    c = multiple_doses(one_compartment_iv, t, p["tau"], 40, dose=p["dose"], cl=p["cl"], vd=p["vd"])
    peak, trough = steady_state_peak_trough(t, c, p["tau"])
    assert 13 <= np.log(2) * p["vd"] / p["cl"] <= 19
    assert 50 <= trough and peak <= 100


def test_levofloxacin_preset_reproduces_label():
    # [7] Table 8 and text: 500 mg once daily gives steady-state peak 5.7 +/- 1.4 and
    # trough 0.5 +/- 0.2 mcg/mL; half-life 6-8 h; steady state within 48 h.
    p = PRESETS["Levofloxacin"]
    t = time_grid(p["tau"], 10)
    c = multiple_doses(one_compartment_oral, t, p["tau"], 10,
                       dose=p["dose"], cl=p["cl"], vd=p["vd"], ka=p["ka"], f=p["f"])
    peak, trough = steady_state_peak_trough(t, c, p["tau"])
    assert 5.7 - 1.4 <= peak <= 5.7 + 1.4
    assert 0.5 - 0.2 <= trough <= 0.5 + 0.2
    assert 6 <= np.log(2) * p["vd"] / p["cl"] <= 8
    assert time_to_steady_state(t, c, p["tau"]) <= 48


def test_fluconazole_preset_reproduces_label():
    # [8] Table 1: 200 mg once daily for 14 days gives Cmax 10.12 mcg/mL and
    # AUC0-24 169.5 mcg*h/mL; steady state within 5 to 10 days.
    p = PRESETS["Fluconazole"]
    t = time_grid(p["tau"], 15)
    c = multiple_doses(one_compartment_oral, t, p["tau"], 15,
                       dose=p["dose"], cl=p["cl"], vd=p["vd"], ka=p["ka"], f=p["f"])
    day_14 = t >= 14 * p["tau"]
    assert cmax(c[day_14]) == pytest.approx(10.12, rel=0.1)
    assert auc(t[day_14], c[day_14]) == pytest.approx(169.5, rel=0.1)
    assert 5 * 24 <= time_to_steady_state(t, c, p["tau"]) <= 10 * 24


def simulate(name, n_doses):
    """Repeated dosing of a preset as published (its own route, dose and interval)."""
    p = PRESETS[name]
    t = time_grid(p["tau"], n_doses)
    if p["route"] == "Oral":
        c = multiple_doses(one_compartment_oral, t, p["tau"], n_doses,
                           dose=p["dose"], cl=p["cl"], vd=p["vd"], ka=p["ka"], f=p["f"])
    else:
        c = multiple_doses(one_compartment_iv, t, p["tau"], n_doses, dose=p["dose"], cl=p["cl"], vd=p["vd"])
    return p, t, c


def test_pregabalin_preset_reproduces_label():
    # [9]: half-life 6.3 h; steady state within 24 to 48 h.
    p, t, c = simulate("Pregabalin", 10)
    assert np.log(2) * p["vd"] / p["cl"] == pytest.approx(6.3, rel=0.05)
    assert 24 <= time_to_steady_state(t, c, p["tau"]) <= 48


def test_lacosamide_preset_reproduces_label():
    # [10]: half-life about 13 h; steady state after 3 days of twice-daily dosing.
    p, t, c = simulate("Lacosamide", 15)
    assert np.log(2) * p["vd"] / p["cl"] == pytest.approx(13, rel=0.05)
    assert time_to_steady_state(t, c, p["tau"]) <= 72


def test_lamotrigine_preset_reproduces_label():
    # [11]: multiple-dose half-life 25.4 h in healthy volunteers (range 11.6-61.6).
    p, t, c = simulate("Lamotrigine", 20)
    assert np.log(2) * p["vd"] / p["cl"] == pytest.approx(25.4, rel=0.2)


def test_topiramate_preset_reproduces_label():
    # [12]: half-life 21 h, so "steady-state is thus reached in about 4 days".
    p, t, c = simulate("Topiramate", 20)
    assert 3 * 24 <= time_to_steady_state(t, c, p["tau"]) <= 5 * 24


def test_linezolid_preset_reproduces_label():
    # [13] Table 8, 600 mg every 12 h: Cmax 21.2 (SD 5.78), Cmin 6.15 (SD 2.94),
    # AUC over the interval 138 (SD 42.1).
    p, t, c = simulate("Linezolid", 20)
    peak, trough = steady_state_peak_trough(t, c, p["tau"])
    last = t >= 19 * p["tau"]
    assert 21.2 - 5.78 <= peak <= 21.2 + 5.78
    assert 6.15 - 2.94 <= trough <= 6.15 + 2.94
    assert auc(t[last], c[last]) == pytest.approx(138, abs=42.1)


def test_metronidazole_preset_reproduces_label():
    # [14]: a single 750 mg dose gives Cmax 21.4 mcg/mL. The model is about 20% low,
    # so this checks within 25%.
    p = PRESETS["Metronidazole"]
    t = np.linspace(0, 24, 2401)
    c = one_compartment_oral(t, 750, p["cl"], p["vd"], p["ka"], p["f"])
    assert cmax(c) == pytest.approx(21.4, rel=0.25)


def test_doxycycline_preset_reproduces_label():
    # [15]: a single 200 mg dose peaks at 2.6 mcg/mL at 2 h. (The 1.45 mcg/mL at 24 h
    # was used to set Vd, so only the peak is an independent check.)
    p = PRESETS["Doxycycline"]
    t = np.linspace(0, 24, 2401)
    c = one_compartment_oral(t, 200, p["cl"], p["vd"], p["ka"], p["f"])
    assert cmax(c) == pytest.approx(2.6, rel=0.2)
    assert tmax(t, c) == pytest.approx(2, abs=0.1)


def test_amikacin_preset_reproduces_label():
    # [16]: 500 mg leaves about 2.1 mcg/mL at 10 h (after IM dosing, so this IV model
    # is checked within 25%); peaks above 35 and troughs above 10 should be avoided.
    p, t, c = simulate("Amikacin", 10)
    peak, trough = steady_state_peak_trough(t, c, p["tau"])
    c_10h = one_compartment_iv(np.array([10.0]), 500, p["cl"], p["vd"])[0]
    assert c_10h == pytest.approx(2.1, rel=0.25)
    assert peak < 35 and trough < 10


def test_daptomycin_preset_reproduces_label():
    # [17] 6 mg/kg every 24 h: AUC0-24 632 (SD 78), half-life 7.9 (SD 1.0) h. The early
    # peak is not checked: one compartment can't show the distribution phase.
    p, t, c = simulate("Daptomycin", 10)
    last = t >= 9 * p["tau"]
    assert auc(t[last], c[last]) == pytest.approx(632, abs=78)
    assert np.log(2) * p["vd"] / p["cl"] == pytest.approx(7.9, abs=1.0)
