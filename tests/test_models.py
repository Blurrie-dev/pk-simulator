"""Checks for infusion, two-compartment, dose schedules and saturable elimination."""
import numpy as np
import pytest

from src.metrics import auc, steady_state_peak_trough, time_to_steady_state
from src.models import (dosing_schedule, multiple_doses, one_compartment_infusion, one_compartment_iv,
                        one_compartment_oral, simulate_ode, superpose, time_grid,
                        two_compartment_constants, two_compartment_infusion, two_compartment_iv,
                        two_compartment_ode, two_compartment_oral)
from src.simulate import Scenario

DOSE, CL, VD = 500, 5, 50     # k = 0.1 /h
V1, Q, V2 = 20, 10, 60        # a clearly bi-exponential drug
T = np.linspace(0, 72, 721)


# --- infusion -----------------------------------------------------------------

def test_infusion_numerical_matches_analytical():
    numerical = simulate_ode(T, [(0, DOSE)], cl=CL, v1=VD, duration=2)
    analytical = one_compartment_infusion(T, DOSE, CL, VD, duration=2)
    np.testing.assert_allclose(numerical, analytical, rtol=1e-4, atol=1e-8)


def test_very_short_infusion_is_a_bolus():
    np.testing.assert_allclose(one_compartment_infusion(T[T > 1], DOSE, CL, VD, duration=1e-6),
                               one_compartment_iv(T[T > 1], DOSE, CL, VD), rtol=1e-5)


def test_long_infusion_plateaus_at_rate_over_clearance():
    # Css = R / CL: the tap and the drain balance. 100 h = 14 half-lives.
    rate = 10  # mg/h
    c = one_compartment_infusion(np.array([100.0]), rate * 200, CL, VD, duration=200)
    assert c[0] == pytest.approx(rate / CL, rel=1e-4)


def test_infusion_auc_is_dose_over_clearance():
    t = np.linspace(0, 200, 20001)
    assert auc(t, one_compartment_infusion(t, DOSE, CL, VD, 3)) == pytest.approx(DOSE / CL, rel=1e-3)


# --- two-compartment ------------------------------------------------------------

def test_two_compartment_rates_are_ordered():
    (a, alpha), (b, beta) = two_compartment_constants(CL, V1, Q, V2)
    assert alpha > beta > 0
    assert a + b == pytest.approx(1 / V1)  # C(0) = dose / V1: all of it starts in V1


@pytest.mark.parametrize("kind", ["iv", "oral", "infusion"])
def test_two_compartment_numerical_matches_analytical(kind):
    if kind == "iv":
        analytical = two_compartment_iv(T, DOSE, CL, V1, Q, V2)
        numerical = two_compartment_ode(T, DOSE, CL, V1, Q, V2)
    elif kind == "oral":
        analytical = two_compartment_oral(T, DOSE, CL, V1, Q, V2, ka=1.2, f=0.7)
        numerical = two_compartment_ode(T, DOSE, CL, V1, Q, V2, ka=1.2, f=0.7)
    else:
        analytical = two_compartment_infusion(T, DOSE, CL, V1, Q, V2, duration=1.5)
        numerical = two_compartment_ode(T, DOSE, CL, V1, Q, V2, duration=1.5)
    np.testing.assert_allclose(numerical, analytical, rtol=1e-4, atol=1e-8)


def test_two_compartment_oral_when_ka_equals_a_rate():
    (_, _), (_, beta) = two_compartment_constants(CL, V1, Q, V2)
    analytical = two_compartment_oral(T, DOSE, CL, V1, Q, V2, ka=beta)
    numerical = two_compartment_ode(T, DOSE, CL, V1, Q, V2, ka=beta)
    np.testing.assert_allclose(numerical, analytical, rtol=1e-4, atol=1e-8)


def test_two_compartment_auc_is_f_dose_over_clearance():
    # Exposure only depends on how fast drug leaves, not where it waits meanwhile.
    t = np.linspace(0, 600, 60001)
    c = two_compartment_oral(t, DOSE, CL, V1, Q, V2, ka=1.0, f=0.6)
    assert auc(t, c) == pytest.approx(0.6 * DOSE / CL, rel=1e-3)


def test_tiny_peripheral_compartment_is_one_compartment():
    c2 = two_compartment_iv(T, DOSE, CL, VD, q=1e-6, v2=1e-3)
    np.testing.assert_allclose(c2, one_compartment_iv(T, DOSE, CL, VD), rtol=1e-3)


def test_two_compartment_needs_a_peripheral_compartment():
    with pytest.raises(ValueError):
        two_compartment_constants(CL, V1, 0, V2)


# --- schedules and superposition --------------------------------------------------

def test_dosing_schedule_what_ifs():
    assert dosing_schedule(12, 3, 100) == [(0, 100), (12, 100), (24, 100)]
    assert dosing_schedule(12, 3, 100, loading_dose=300)[0] == (0, 300)
    assert dosing_schedule(12, 3, 100, missed=[2]) == [(0, 100), (24, 100)]
    assert dosing_schedule(12, 3, 100, missed=[2], double_next=True) == [(0, 100), (24, 200)]
    assert dosing_schedule(12, 3, 100, late={2: 5}) == [(0, 100), (17, 100), (24, 100)]


def test_multiple_doses_unchanged_by_refactor():
    t = time_grid(8, 5)
    expected = sum(np.where(t >= i * 8, one_compartment_oral(np.clip(t - i * 8, 0, None), 250, 2.5, 50, 1.0), 0)
                   for i in range(5))
    np.testing.assert_allclose(multiple_doses(one_compartment_oral, t, 8, 5, dose=250, cl=2.5, vd=50, ka=1.0),
                               expected)


@pytest.mark.parametrize("route", ["oral", "iv", "infusion"])
def test_superposition_matches_ode_for_irregular_schedule(route):
    # loading dose, a missed dose and a late dose, two compartments
    schedule = dosing_schedule(12, 6, 200, loading_dose=600, missed=[3], late={5: 4.5})
    t = time_grid(12, 6, extra=[time for time, _ in schedule] + [time + 1 for time, _ in schedule])
    if route == "oral":
        a = superpose(two_compartment_oral, t, schedule, cl=CL, v1=V1, q=Q, v2=V2, ka=0.8, f=0.9)
        n = simulate_ode(t, schedule, cl=CL, v1=V1, q=Q, v2=V2, ka=0.8, f=0.9)
    elif route == "iv":
        a = superpose(two_compartment_iv, t, schedule, cl=CL, v1=V1, q=Q, v2=V2)
        n = simulate_ode(t, schedule, cl=CL, v1=V1, q=Q, v2=V2)
    else:
        a = superpose(two_compartment_infusion, t, schedule, cl=CL, v1=V1, q=Q, v2=V2, duration=1)
        n = simulate_ode(t, schedule, cl=CL, v1=V1, q=Q, v2=V2, duration=1)
    np.testing.assert_allclose(n, a, rtol=1e-4, atol=1e-7)


def test_loading_dose_reaches_steady_state_at_once():
    # IV bolus, one compartment: loading dose = dose / (1 - exp(-k tau)) puts every
    # trough at its steady-state value from the first interval on.
    k, tau = CL / VD, 8
    loading = 100 / (1 - np.exp(-k * tau))
    s = Scenario(route="IV bolus", dose=100, tau=tau, cl=CL, vd=VD, loading_dose=loading)
    t = s.grid(10)
    c = s.repeated(t, 10)
    troughs = [c[np.isclose(t, (i + 1) * tau)][0] for i in range(9)]
    np.testing.assert_allclose(troughs, troughs[-1], rtol=1e-9)


# --- saturable (Michaelis-Menten) elimination ---------------------------------------

def test_saturable_with_huge_km_is_linear():
    # far below Km, vmax * C / (km + C) ~ (vmax / km) * C: plain clearance
    km = 1e6
    c = simulate_ode(T, [(0, DOSE)], cl=0, v1=VD, vmax=CL * km, km=km)
    np.testing.assert_allclose(c, one_compartment_iv(T, DOSE, CL, VD), rtol=1e-3)


def test_saturable_infusion_steady_state():
    # Constant input R balances vmax * Css / (km + Css): Css = km * R / (vmax - R)
    vmax, km, rate = 20.0, 4.0, 15.0
    t = np.linspace(0, 2000, 201)
    c = simulate_ode(t, [(0, rate * 2000)], cl=0, v1=45, vmax=vmax, km=km, duration=2000)
    assert c[-1] == pytest.approx(km * rate / (vmax - rate), rel=1e-3)


def test_saturable_exposure_rises_more_than_dose():
    t = np.linspace(0, 400, 4001)
    low = auc(t, simulate_ode(t, [(0, 300)], cl=0, v1=45, vmax=20, km=4))
    high = auc(t, simulate_ode(t, [(0, 600)], cl=0, v1=45, vmax=20, km=4))
    assert high / low > 2.2  # doubling the dose more than doubles AUC


def test_scenario_saturable_steady_state_matches_formula():
    s = Scenario(route="IV infusion", dose=8 * 15, tau=8, duration=8, saturable=True, vmax=20, km=4, vd=45)
    n = min(s.doses_needed(), 200)
    t = s.grid(n)
    c = s.repeated(t, n)
    _peak, trough = steady_state_peak_trough(t, c, s.tau)
    assert trough == pytest.approx(s.average_steady_state(), rel=0.01)


def test_saturated_scenario_has_no_steady_state():
    s = Scenario(route="IV bolus", dose=200, tau=8, saturable=True, vmax=20, km=4)
    assert s.saturated() and s.average_steady_state() == np.inf


# --- scenario plumbing ---------------------------------------------------------------

@pytest.mark.parametrize("route", ["Oral", "IV bolus", "IV infusion"])
@pytest.mark.parametrize("two", [False, True])
def test_scenario_linear_matches_ode(route, two):
    s = Scenario(route=route, two_compartment=two, q=Q, v2=V2, missed=(2,), late=((4, 3),))
    linear = s.repeated(s.grid(6), 6)
    ode = Scenario(**{**s.__dict__, "saturable": True, "vmax": s.cl * 1e7, "km": 1e7}).repeated(s.grid(6), 6)
    np.testing.assert_allclose(ode, linear, rtol=1e-3, atol=1e-6)


def test_scenario_steady_state_timing_uses_terminal_half_life():
    s = Scenario(route="IV bolus", two_compartment=True, vd=V1, cl=CL, q=Q, v2=V2, tau=12)
    (_, _), (_, beta) = two_compartment_constants(CL, V1, Q, V2)
    assert s.half_life() == pytest.approx(np.log(2) / beta)
    n = s.doses_needed()
    t = s.grid(n)
    tss = time_to_steady_state(t, s.repeated(t, n), s.tau)
    assert 3.5 <= tss / s.half_life() <= 5.5


def test_average_steady_state_and_time_in_window():
    from src.metrics import average_steady_state, time_in_window
    s = Scenario()
    n = 2 * s.doses_needed()
    t = s.grid(n)
    c = s.repeated(t, n)
    assert average_steady_state(t, c, s.tau) == pytest.approx(s.average_steady_state(), rel=1e-3)
    assert time_in_window(t, c, s.tau, (0, 1e9)) == pytest.approx(1)
    assert time_in_window(t, c, s.tau, (1e8, 1e9)) == 0
