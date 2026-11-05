"""Adjusting typical-adult parameters for one patient: body size and kidney function.

The presets describe a typical 70 kg adult. Two things change most from patient
to patient and are easy to measure: weight and kidney function (estimated as
creatinine clearance, CrCl, from a blood test). Bigger bathtubs and drains for
bigger people; a narrower drain when the kidneys work less well.
"""
from dataclasses import replace

# The patient the presets describe. With these inputs nothing is adjusted.
REFERENCE = dict(weight=70.0, age=40, female=False, scr=1.0)


def cockcroft_gault(age, weight, scr, female=False):
    """Creatinine clearance (mL/min), Cockcroft-Gault equation.

    CrCl = (140 - age) * weight / (72 * serum creatinine), x 0.85 for women.
    age in years, weight in kg, scr (serum creatinine) in mg/dL.
    Creatinine is made by muscle at a steady rate and cleared by the kidneys, so
    the higher it sits in blood, the slower the kidneys are clearing it.
    """
    crcl = (140 - age) * weight / (72 * scr)
    return crcl * 0.85 if female else crcl


REFERENCE_CRCL = cockcroft_gault(**REFERENCE)  # 97 mL/min


def size_factor(weight, exponent):
    """Allometric scaling: (weight / 70) ** exponent.

    Volumes scale with weight (exponent 1). Clearances, like metabolic rate,
    scale with weight to the power 0.75: twice the size doesn't mean twice the
    liver blood flow.
    """
    return (weight / REFERENCE["weight"]) ** exponent


def clearance_factor(weight, crcl, fe):
    """How much the patient's clearance differs from the reference patient's.

    The renal part, fraction fe of the total, follows creatinine clearance (which
    already includes weight). The rest, mostly liver metabolism, scales with size.
    """
    return (1 - fe) * size_factor(weight, 0.75) + fe * crcl / REFERENCE_CRCL


def adjust(scenario, weight, crcl, fe):
    """The scenario's PK parameters scaled to a patient. Doses are left alone."""
    return replace(
        scenario,
        cl=scenario.cl * clearance_factor(weight, crcl, fe),
        vmax=scenario.vmax * size_factor(weight, 0.75),  # saturable metabolism scales like CL
        vd=scenario.vd * size_factor(weight, 1),
        v2=scenario.v2 * size_factor(weight, 1),
        q=scenario.q * size_factor(weight, 0.75),
    )
