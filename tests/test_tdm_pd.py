"""Checks for Bayesian level fitting and the effect model."""
import numpy as np
import pytest

from src import pd, tdm
from src.simulate import Scenario

CV = {"cl": 0.3, "vd": 0.2}


def test_fit_recovers_true_parameters_from_rich_data():
    prior = Scenario(route="IV bolus", dose=500, tau=12, cl=5, vd=50)
    truth = Scenario(route="IV bolus", dose=500, tau=12, cl=7, vd=40)
    times = np.array([1, 2, 4, 6, 8, 11, 13, 16, 20, 23.5])
    fit = tdm.fit(prior, times, tdm.predict(truth, times), CV, sd_add=0.01, cv_prop=0.01)
    assert fit.scenario.cl == pytest.approx(7, rel=0.01)
    assert fit.scenario.vd == pytest.approx(40, rel=0.01)


def test_one_level_is_shrunk_towards_the_prior():
    # A single trough can't pin down both CL and V: MAP lands between the prior and the truth.
    prior = Scenario(route="Oral", dose=250, tau=8, cl=2.5, vd=50, ka=1.0, f=0.8)
    truth = Scenario(route="Oral", dose=250, tau=8, cl=4.0, vd=50, ka=1.0, f=0.8)
    fit = tdm.fit(prior, [39.5], tdm.predict(truth, [39.5]), CV)
    assert 2.5 < fit.scenario.cl < 4.0


def test_no_levels_returns_the_prior():
    prior = Scenario()
    fit = tdm.fit(prior, [], [], CV)
    assert fit.scenario == prior


def test_fit_works_for_saturable_elimination():
    prior = Scenario(route="Oral", dose=100, tau=8, saturable=True, vmax=20.4, km=4, vd=45.5, ka=0.3, f=0.92)
    truth = Scenario(**{**prior.__dict__, "vmax": 17.0})
    times = np.array([100.0, 150.0, 190.0])
    fit = tdm.fit(prior, times, tdm.predict(truth, times), {"vmax": 0.3}, parameters=("vmax",),
                  sd_add=0.01, cv_prop=0.01)
    assert fit.scenario.vmax == pytest.approx(17.0, rel=0.02)


def test_predict_accepts_unsorted_times():
    s = Scenario()
    np.testing.assert_allclose(tdm.predict(s, [20, 5]), tdm.predict(s, [5, 20])[::-1])


def test_effect_site_without_delay_is_plasma():
    t = np.linspace(0, 24, 241)
    c = np.exp(-0.2 * t)
    np.testing.assert_allclose(pd.effect_site(t, c, np.inf), c)
    np.testing.assert_allclose(pd.effect_site(t, c, 1e4)[5:], c[5:], rtol=1e-3)


def test_effect_site_approaches_a_constant_level_with_ke0_half_life():
    t = np.linspace(0, 10, 101)
    ce = pd.effect_site(t, np.ones_like(t), ke0=np.log(2))  # half-life 1 h
    assert ce[np.isclose(t, 1)][0] == pytest.approx(0.5)
    assert ce[np.isclose(t, 3)][0] == pytest.approx(0.875)


def test_effect_site_matches_ode_for_a_real_curve():
    from scipy.integrate import solve_ivp
    s = Scenario()
    t = s.grid(3)
    c = s.repeated(t, 3)
    ce = pd.effect_site(t, c, 0.5)
    sol = solve_ivp(lambda tt, y: 0.5 * (np.interp(tt, t, c) - y), (0, t[-1]), [0.0], t_eval=t,
                    rtol=1e-9, atol=1e-12, max_step=0.05)
    np.testing.assert_allclose(ce, sol.y[0], atol=1e-4)


def test_emax_model():
    assert pd.emax_effect(10, 0, 100, 10) == pytest.approx(50)          # half effect at EC50
    assert pd.emax_effect(1e9, 5, 100, 10) == pytest.approx(105)       # saturates at E0 + Emax
    assert pd.emax_effect(20, 0, 100, 10, hill=4) > pd.emax_effect(20, 0, 100, 10)  # steeper
    assert pd.emax_effect(0, 70, -30, 5) == 70                          # inhibitory starts at E0


def test_settings_round_trip_and_validation():
    from src import settings
    state = {"preset": "Phenytoin", "dose": 130, "cl": 1.46, "adjust": True, "sex": "Female", "comparisons": []}
    values, rejected = settings.from_json(settings.to_json(state))
    assert values == {k: v for k, v in state.items() if k != "comparisons"} and not rejected
    values, rejected = settings.parse(settings.to_query(state))   # strings, as in a URL
    assert values["dose"] == 130 and values["adjust"] is True and values["cl"] == 1.46
    values, rejected = settings.parse({"dose": "99999", "route": "Rectal", "hack": 1, "tau": "12"})
    assert values == {"tau": 12} and set(rejected) == {"dose", "route", "hack"}
    assert settings.clamp("dose", 5000) == 2000
