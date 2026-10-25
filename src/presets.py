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
    "Theophylline": dict(
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
    "Levetiracetam": dict(
        route="Oral",
        dose=1500, tau=12,              # recommended adult dose, 1500 mg twice daily [4]
        cl=per_kg_clearance(0.96),      # total body clearance [4]
        vd=round(per_kg_clearance(0.96) * 7 / 0.693, 1),  # not in [4]: CL x t1/2 / ln(2), t1/2 = 7 h [4]
        ka=3.7,                         # not in [4]: chosen so Tmax = 1.0 h, "about an hour" [4]
        f=1.0,                          # oral bioavailability 100% [4]
        window=(0.0, 0.0),              # no target range in [4], so no band
        note="No therapeutic range in the label, so no band is drawn. PK is linear over "
             "500-5000 mg, so adding up doses (superposition) is valid.",
    ),
    "Levofloxacin": dict(
        route="Oral",
        dose=500, tau=24,               # 500 mg once daily, the regimen measured in [7] Table 8
        cl=round(175 * 0.99 * 60 / 1000, 1),  # CL/F 175 mL/min x F 0.99 [7] Table 8
        vd=round(102 * 0.99, 1),        # Vd/F 102 L x F 0.99 [7] Table 8
        ka=3.25,                        # not in [7]: chosen so Tmax = 1.1 h, as in [7] Table 8
        f=0.99,                         # absolute bioavailability ~99% [7]
        window=(0.0, 0.0),              # no target range in [7], so no band
        note="No concentration range in the label, so no band is drawn. Linear PK; steady "
             "state within 48 h of once-daily dosing.",
    ),
    "Fluconazole": dict(
        route="Oral",
        dose=200, tau=24,               # 200 mg once daily, a regimen measured in [8] Table 1
        cl=per_kg_clearance(0.23),      # mean body clearance 0.23 mL/min/kg [8]
        vd=round(per_kg_clearance(0.23) * 30 / 0.693, 1),  # not in [8]: CL x t1/2 / ln(2), t1/2 = 30 h [8]
        ka=3.35,                        # not in [8]: chosen so Tmax = 1.5 h, inside its "1 and 2 hours" [8]
        f=0.9,                          # "over 90%" [8]; the lower bound
        window=(0.0, 0.0),              # no target range in [8], so no band
        note="No concentration range in the label, so no band is drawn. The 30 h half-life "
             "means steady state takes about 6 days, which is why the label uses a double "
             "loading dose on day 1.",
    ),
    "Valproate": dict(
        route="IV bolus",
        dose=250, tau=6,                # 1000 mg/day = 14 mg/kg/day, inside the initial 10-15 mg/kg/day,
                                        # divided as [6] requires above 250 mg/day; 6 h is a chosen interval
        cl=0.56,                        # 0.56 L/h per 1.73 m2, total valproate [6]
        vd=11.0,                        # 11 L per 1.73 m2 [6]
        window=(50.0, 100.0),           # therapeutic range 50-100 mcg/mL [6]
        note="Values are per 1.73 m2, taken as a 70 kg adult. Total valproate is only roughly "
             "linear: protein binding saturates at higher concentrations.",
    ),
    "Gentamicin": dict(
        route="IV bolus",
        dose=70, tau=8,                 # 3 mg/kg/day in three equal doses every 8 h [2]
        cl=per_kg_clearance(1.29),      # adults 20-39 y, normal renal function [3]
        vd=round(0.23 * WEIGHT, 1),     # 0.23 L/kg [3]
        window=(4.0, 12.0),             # expected peak 4-6, avoid prolonged peaks above 12 mcg/mL [2]
        note="Band = peak targets. Peaks should reach the band and troughs should fall below "
             "2 mg/L. Given as a short infusion in practice; modelled here as a bolus.",
    ),
    "Vancomycin": dict(
        route="IV bolus",
        dose=1000, tau=12,              # 2 g/day as 1 g every 12 h [5]
        cl=round(0.058 * WEIGHT, 1),    # 0.058 L/kg/h [5]
        vd=round(0.365 * WEIGHT, 1),    # midpoint of 0.3-0.43 L/kg [5]
        window=(0.0, 0.0),              # [5] says to measure troughs but gives no range, so no band
        note="No concentration range in the label, so no band is drawn. Vancomycin really has a "
             "fast distribution phase, so this one-compartment model gets the trough about right "
             "but underestimates the first hours after each dose.",
    ),
    "Pregabalin": dict(
        route="Oral",
        dose=150, tau=12,               # 300 mg/day in 2 doses, inside 150-600 mg/day in 2-3 doses [9]
        vd=round(0.5 * WEIGHT * 0.9, 1),  # Vd/F 0.5 L/kg [9] x F
        cl=round(round(0.5 * WEIGHT * 0.9, 1) * 0.693 / 6.3, 2),  # not in [9]: Vd x ln(2) / t1/2, t1/2 6.3 h [9]
        ka=3.6,                         # not in [9]: chosen so Tmax = 1.0 h, "within 1.5 hours" [9]
        f=0.9,                          # "greater than or equal to 90%" [9]; the lower bound
        window=(0.0, 0.0),              # no target range in [9]
        note="No concentration range in the label, so no band is drawn. Linear PK; cleared "
             "almost entirely by the kidneys.",
    ),
    "Lacosamide": dict(
        route="Oral",
        dose=200, tau=12,               # monotherapy 150-200 mg twice daily [10]
        vd=round(0.6 * WEIGHT, 1),      # 0.6 L/kg [10]
        cl=round(0.6 * WEIGHT * 0.693 / 13, 2),  # not in [10]: Vd x ln(2) / t1/2, t1/2 13 h [10]
        ka=1.34,                        # not in [10]: chosen so Tmax = 2.5 h, middle of "1-to-4-hour" [10]
        f=1.0,                          # bioavailability approximately 100% [10]
        window=(0.0, 0.0),              # no target range in [10]
        note="No concentration range in the label, so no band is drawn. Dose proportional "
             "from 100 to 800 mg.",
    ),
    "Lamotrigine": dict(
        route="Oral",
        dose=150, tau=12,               # 300 mg/day, inside the usual maintenance 225-375 mg/day [11]
        cl=round(0.58 * WEIGHT * 60 / 1000 * 0.98, 2),  # multiple-dose CL/F 0.58 mL/min/kg [11] x F
        vd=round(1.1 * WEIGHT * 0.98, 1),  # Vd/F 0.9-1.3 L/kg [11], midpoint, x F
        ka=2.63,                        # not in [11]: chosen so Tmax = 1.7 h, the multiple-dose mean [11]
        f=0.98,                         # absolute bioavailability 98% [11]
        window=(0.0, 0.0),              # no target range in [11]
        note="No concentration range in the label, so no band is drawn. Values are for adults "
             "taking no interacting drugs: valproate roughly halves lamotrigine clearance.",
    ),
    "Topiramate": dict(
        route="Oral",
        dose=200, tau=12,               # monotherapy 400 mg/day in two divided doses [12]
        cl=round(25 * 60 / 1000, 2),    # CL/F 20-30 mL/min [12], midpoint (apparent, so F = 1)
        vd=round(25 * 60 / 1000 * 21 / 0.693, 1),  # not in [12]: CL/F x t1/2 / ln(2), t1/2 21 h [12]
        ka=2.11,                        # not in [12]: chosen so Tmax = 2 h [12]
        f=1.0,                          # [12] gives only relative bioavailability; CL/F already includes F
        window=(0.0, 0.0),              # no target range in [12]
        note="No concentration range in the label, so no band is drawn. Uses apparent "
             "clearance (CL/F), so F is set to 1.",
    ),
    "Linezolid": dict(
        route="Oral",
        dose=600, tau=12,               # 600 mg every 12 h, measured in [13] Table 8
        cl=round(80 * 60 / 1000, 2),    # CL 80 mL/min at 600 mg every 12 h [13] Table 8
        vd=45.0,                        # 40-50 L at steady state [13], midpoint
        ka=3.49,                        # not in [13]: chosen so Tmax = 1.03 h, as in [13] Table 8
        f=1.0,                          # absolute bioavailability approximately 100% [13]
        window=(0.0, 0.0),              # no target range in [13]
        note="No concentration range in the label, so no band is drawn. Clearance falls slightly "
             "at higher doses, a small nonlinearity the label calls minor.",
    ),
    "Metronidazole": dict(
        route="Oral",
        dose=500, tau=6,                # 7.5 mg/kg every 6 h, "approximately 500 mg for a 70 kg adult" [14]
        cl=round(750 / 223, 2),         # not in [14]: CL/F = dose / AUC = 750 mg / 223 mg*h/L [14]
        vd=round(750 / 223 * 8 / 0.693, 1),  # not in [14]: CL/F x t1/2 / ln(2), t1/2 8 h [14]
        ka=2.07,                        # not in [14]: chosen so Tmax = 1.6 h, as measured in [14]
        f=1.0,                          # apparent values (CL/F, V/F) already include F
        window=(0.0, 0.0),              # no target range in [14]
        note="No concentration range in the label, so no band is drawn. Peaks run about 20% "
             "below the label's measured Cmax.",
    ),
    "Doxycycline": dict(
        route="Oral",
        dose=100, tau=12,               # 100 mg every 12 h for more severe infections [15]
        vd=61.1,                        # not in [15]: set so 200 mg gives 1.45 mcg/mL at 24 h, as in [15]
        cl=round(61.1 * 0.693 / 20, 2),  # not in [15]: Vd x ln(2) / t1/2, t1/2 18-22 h (20) [15]
        ka=2.08,                        # not in [15]: chosen so Tmax = 2 h [15]
        f=1.0,                          # "virtually completely absorbed" [15]
        window=(0.0, 0.0),              # no target range in [15]
        note="No concentration range in the label, so no band is drawn. The usual maintenance "
             "dose is 100 mg/day; 100 mg twice daily is for more severe infections.",
    ),
    "Amikacin": dict(
        route="IV bolus",
        dose=525, tau=12,               # 15 mg/kg/day as 7.5 mg/kg every 12 h [16]
        cl=round(100 * 60 / 1000, 2),   # mean serum clearance about 100 mL/min [16]
        vd=24.0,                        # mean apparent volume of distribution 24 L [16]
        window=(10.0, 35.0),            # avoid peaks above 35 and troughs above 10 mcg/mL [16]
        note="Band spans the label's limits: peaks should stay under 35 mg/L and troughs fall "
             "below 10 mg/L. Given IM or as an infusion in practice; modelled here as a bolus.",
    ),
    "Daptomycin": dict(
        route="IV bolus",
        dose=420, tau=24,               # 6 mg/kg once every 24 h (S. aureus bacteremia) [17]
        cl=round(9.1 * WEIGHT / 1000, 2),  # CL 9.1 mL/h/kg at 6 mg/kg [17]
        vd=round(0.101 * WEIGHT, 2),    # Vss 0.101 L/kg at 6 mg/kg [17]
        window=(0.0, 0.0),              # no target range in [17]
        note="No concentration range in the label, so no band is drawn. Exposure (AUC) matches "
             "the label, but like vancomycin the early peak is underestimated: daptomycin has a "
             "distribution phase this one-compartment model doesn't show.",
    ),
}
