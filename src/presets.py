"""Parameter presets for the app. Every real-drug value is cited in the README
("Drug presets"); the [n] tags below match the reference list there.

Per-kg values are scaled to a 70 kg adult and rounded to 0.1 to fit the
app's slider steps.
"""
WEIGHT = 70  # kg


def per_kg_clearance(ml_min_kg):
    """mL/min/kg (how papers report it) -> L/h for a 70 kg adult (what the models use)."""
    return round(ml_min_kg * WEIGHT * 60 / 1000, 1)


PRESETS = {
    "Hypothetical drug": dict(
        route="Oral", dose=250, tau=8, cl=2.5, vd=50.0, ka=1.0, f=0.8, window=(5.0, 15.0),
        note="Made-up parameters, the same ones as the README figure.",
    ),
    "Theophylline (oral)": dict(
        route="Oral",
        dose=200, tau=8,                # 600 mg/day divided every 8 h: final adult titration step [1]
        cl=per_kg_clearance(0.65),      # healthy non-smoking adults [1]
        vd=round(0.45 * WEIGHT, 1),     # 0.45 L/kg [1]
        ka=2.0,                         # not in [1]: chosen so Tmax = 1.6 h, inside its "1-2 hr" [1]
        f=1.0,                          # "rapidly and completely absorbed" [1]
        window=(10.0, 20.0),            # therapeutic range 10-20 mcg/mL [1]
        note="Narrow therapeutic range: the whole curve should stay inside the band. "
             "The label titrates this dose using measured serum concentrations.",
    ),
    "Gentamicin (IV)": dict(
        route="IV bolus",
        dose=70, tau=8,                 # 3 mg/kg/day in three equal doses every 8 h [2]
        cl=per_kg_clearance(1.29),      # adults 20-39 y, normal renal function [3]
        vd=round(0.23 * WEIGHT, 1),     # 0.23 L/kg [3]
        window=(4.0, 12.0),             # expected peak 4-6, avoid prolonged peaks above 12 mcg/mL [2]
        note="Band = peak targets. Peaks should reach the band and troughs should fall below "
             "2 mg/L. Given as a short infusion in practice; modelled here as a bolus.",
    ),
}
