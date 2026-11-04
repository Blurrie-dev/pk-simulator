"""Checks for patient adjustment and population variability."""
import numpy as np
import pytest

from src import patient, population
from src.simulate import Scenario


def test_cockcroft_gault_known_value():
    # 60 y, 72 kg, creatinine 1.0: (140 - 60) * 72 / 72 = 80 mL/min; women x 0.85
    assert patient.cockcroft_gault(60, 72, 1.0) == pytest.approx(80)
    assert patient.cockcroft_gault(60, 72, 1.0, female=True) == pytest.approx(68)


def test_reference_patient_changes_nothing():
    s = Scenario(two_compartment=True)
    adjusted = patient.adjust(s, patient.REFERENCE["weight"], patient.REFERENCE_CRCL, fe=0.7)
    for name in ("cl", "vd", "q", "v2", "vmax"):
        assert getattr(adjusted, name) == pytest.approx(getattr(s, name))


def test_kidney_function_only_affects_renal_part():
    s = Scenario(cl=10)
    half = patient.REFERENCE_CRCL / 2
    assert patient.adjust(s, 70, half, fe=1.0).cl == pytest.approx(5)    # all renal: halves
    assert patient.adjust(s, 70, half, fe=0.0).cl == pytest.approx(10)   # all hepatic: unchanged
    assert patient.adjust(s, 70, half, fe=0.6).cl == pytest.approx(7)    # 4 + 6 / 2


def test_weight_scaling_is_allometric():
    s = Scenario(cl=10, vd=50)
    big = patient.adjust(s, 140, patient.REFERENCE_CRCL * 2, fe=0.0)
    assert big.vd == pytest.approx(100)
    assert big.cl == pytest.approx(10 * 2**0.75)


def test_population_is_reproducible_and_log_normal():
    s = Scenario()
    a = population.sample_patients(s, {"cl": 0.3}, 4000, seed=1)
    b = population.sample_patients(s, {"cl": 0.3}, 4000, seed=1)
    cl = np.array([p.cl for p in a])
    assert [p.cl for p in b] == list(cl)
    assert np.median(cl) == pytest.approx(s.cl, rel=0.03)          # typical value = median
    assert np.std(cl) / np.mean(cl) == pytest.approx(0.3, rel=0.08)  # requested CV
    assert all(p.vd == s.vd for p in a)                              # no CV given: no variability


def test_percentiles_bracket_the_median_and_attainment_is_a_fraction():
    s = Scenario()
    patients = population.sample_patients(s, {"cl": 0.3, "vd": 0.2, "ka": 0.4}, 200, seed=2)
    n = s.doses_needed()
    t = s.grid(n)
    curves = population.simulate(patients, t, n)
    bands = population.percentiles(curves)
    assert np.all(bands[5] <= bands[50]) and np.all(bands[50] <= bands[95])
    peaks, troughs = population.steady_state_summary(curves, t, s.tau)
    pta = population.target_attainment(peaks, troughs, (5, 15))
    assert 0 < pta < 1
    assert population.target_attainment(peaks, troughs, (0, 1e9)) == 1
