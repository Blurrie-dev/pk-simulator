"""Checks for the regimen finder and loading dose."""
import numpy as np
import pytest

from src.dosing import loading_dose, steady_state_per_mg, suggest_regimens
from src.metrics import steady_state_peak_trough
from src.presets import PRESETS
from src.simulate import Scenario


def test_loading_dose_matches_formula_for_iv_bolus():
    s = Scenario(route="IV bolus", dose=100, tau=8, cl=5, vd=50)
    assert loading_dose(s) == pytest.approx(100 / (1 - np.exp(-0.1 * 8)), rel=1e-3)


def test_steady_state_per_mg_scales_linearly():
    s = Scenario()
    peak, trough = steady_state_per_mg(s, 8)
    n = 2 * s.doses_needed()
    t = s.grid(n)
    real = steady_state_peak_trough(t, s.repeated(t, n), 8)
    np.testing.assert_allclose(real, (250 * peak, 250 * trough), rtol=1e-9)


def test_suggested_regimens_land_where_they_say():
    s = Scenario.from_preset(PRESETS["Theophylline"])
    rows = suggest_regimens(s, (10, 20), taus=(6, 8, 12, 24), dose_step=10)
    assert [r["tau"] for r in rows] == [6, 8, 12, 24]
    for r in rows:
        check = Scenario.from_preset(PRESETS["Theophylline"], dose=r["dose"], tau=r["tau"])
        n = 2 * check.doses_needed()
        t = check.grid(n)
        peak, trough = steady_state_peak_trough(t, check.repeated(t, n), r["tau"])
        assert peak == pytest.approx(r["peak"], rel=1e-6)
        assert trough == pytest.approx(r["trough"], rel=1e-6)
        assert r["fits"] == (10 <= trough and peak <= 20)
    assert rows[0]["fits"] and not rows[-1]["fits"]  # 8 h half-life: every 6 h fits, daily can't


def test_suggest_regimens_rejects_bad_inputs():
    with pytest.raises(ValueError):
        suggest_regimens(Scenario(), (0, 10))
    with pytest.raises(ValueError):
        suggest_regimens(Scenario(saturable=True), (5, 10))
