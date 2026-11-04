"""Smoke tests: drive the real Streamlit app headless and check nothing errors."""
from pathlib import Path

import pytest
from streamlit.testing.v1 import AppTest

from src.presets import PRESETS

APP = Path(__file__).resolve().parent.parent / "app.py"


def run(**settings):
    at = AppTest.from_file(str(APP), default_timeout=60)
    at.run()
    for key, value in settings.items():
        if key == "preset":
            at.selectbox(key="preset").select(value).run()
        elif key in ("route", "compartments", "elimination"):
            at.radio(key=key).set_value(value).run()
        elif key in ("adjust", "log"):
            at.toggle(key=key).set_value(value).run()
        else:
            at.number_input(key=key).set_value(value).run()
    assert not at.exception, [e.message for e in at.exception]
    return at


@pytest.mark.parametrize("name", list(PRESETS))
def test_every_preset_runs(name):
    run(preset=name)


@pytest.mark.parametrize("route", ["Oral", "IV bolus", "IV infusion"])
@pytest.mark.parametrize("compartments", ["One", "Two"])
@pytest.mark.parametrize("elimination", ["Linear", "Saturable"])
def test_every_model_combination_runs(route, compartments, elimination):
    run(route=route, compartments=compartments, elimination=elimination)


def test_patient_adjustment_and_what_ifs_run():
    at = run(adjust=True, log=True, missed=3, late_dose=5, loading=500.0)
    assert any("Adjusted for this patient" in i.value for i in at.info)


def test_regimen_finder_applies_its_suggestion():
    at = run(preset="Theophylline")
    use = [b for b in at.button if b.label == "Use this regimen"]
    assert use
    use[0].click().run()
    assert not at.exception
    assert at.session_state.tau in (4, 6, 8, 12, 24, 48)


def test_saturated_phenytoin_warns():
    at = run(preset="Phenytoin", dose_box=500.0)
    assert any("no steady state" in e.value.lower() for e in at.error)
