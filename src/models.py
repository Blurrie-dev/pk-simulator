"""Pharmacokinetic models.

Units used throughout (keep them consistent and the maths works out):
    t       time, hours (h)
    dose    milligrams (mg)
    cl      clearance, litres per hour (L/h)
    vd      volume of distribution, litres (L)
    ka      absorption rate constant, per hour (1/h)
    f       bioavailability, fraction 0-1
    output  plasma concentration, mg/L

Every function takes an array of times `t` (starting at 0) and returns the plasma
concentration at those times.

Two ways of getting a concentration-time curve appear here:
  * Analytical: the rate equations have a closed-form solution, so we evaluate
    the formula directly.
  * Numerical: we give scipy.integrate.solve_ivp the rate equations ("how fast is
    each amount changing right now?") and it steps forward in time from the
    starting amounts, building the curve piece by piece. We use this when there
    is no simple formula (two-compartment) and to cross-check the formulas
    (one-compartment).
"""
import numpy as np
from scipy.integrate import solve_ivp

# solve_ivp accuracy settings. rtol = allowed relative error per step (1e-8 is about
# 8 significant digits); atol = absolute error floor, so the solver doesn't waste
# effort chasing amounts that are effectively zero. The defaults (1e-3, 1e-6) are
# only good to roughly 0.1%; tighter values make solver error negligible.
RTOL, ATOL = 1e-8, 1e-10


def _solve(rates, start, t):
    """Integrate `rates` from t=0 with solve_ivp; return the amount in each
    compartment at times t (one row per compartment, one column per time)."""
    # Integrate from 0 to the last requested time; t_eval asks for the solution
    # at exactly our time points. solve_ivp refuses a zero-length interval, so if
    # only t=0 is requested any end time will do.
    t_end = max(t.max(), 1.0)
    sol = solve_ivp(rates, (0, t_end), start, t_eval=t, rtol=RTOL, atol=ATOL)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.y


def one_compartment_iv(t, dose, cl, vd):
    """IV bolus into one well-mixed compartment, analytical solution.

    The whole dose is in plasma at t=0, so C0 = dose / Vd. Elimination is
    first-order: a constant fraction k = CL / Vd of the drug leaves per hour, so
    concentration decays exponentially: C(t) = C0 * exp(-k * t).
    """
    k = cl / vd
    return dose / vd * np.exp(-k * np.asarray(t))


def one_compartment_oral(t, dose, cl, vd, ka, f=1.0):
    """Oral dose with first-order absorption, analytical solution (Bateman function).

    Only a fraction F of the dose reaches the circulation. It moves from gut to
    plasma at ka * (amount in gut) and is eliminated at k * (amount in plasma).
    Solving those two equations gives

        C(t) = F * dose * ka / (Vd * (ka - k)) * (exp(-k*t) - exp(-ka*t))

    an elimination exponential minus an absorption exponential.
    """
    t = np.asarray(t)
    k = cl / vd
    if np.isclose(ka, k):
        # The formula divides by (ka - k), which is zero here. The limit as
        # ka -> k is this instead: same curve, no division by zero.
        return f * dose * k / vd * t * np.exp(-k * t)
    return f * dose * ka / (vd * (ka - k)) * (np.exp(-k * t) - np.exp(-ka * t))


def one_compartment_ode(t, dose, cl, vd, ka=None, f=1.0):
    """The same one-compartment model, solved numerically with solve_ivp.

    Exists to check the analytical formulas above. ka=None means IV bolus (dose
    starts in plasma); otherwise the bioavailable dose F * dose starts in the gut.
    """
    t = np.asarray(t, dtype=float)

    if ka is None:  # IV: everything starts in plasma, nothing in the gut
        ka, start = 0.0, [0.0, dose]
    else:           # oral: bioavailable fraction starts in the gut
        start = [f * dose, 0.0]

    def rates(_time, amounts):
        # solve_ivp calls this repeatedly with the current amounts (mg) and asks
        # how fast each is changing (mg/h). Time isn't used: the rates depend only
        # on how much drug is where, not on the clock.
        gut, plasma = amounts
        absorption = ka * gut           # mg/h moving gut -> plasma
        elimination = cl * plasma / vd  # CL x concentration = mg/h removed
        return [-absorption, absorption - elimination]

    gut, plasma = _solve(rates, start, t)
    return plasma / vd  # amount in plasma / volume = concentration


def two_compartment_iv(t, dose, cl, v1, q, v2):
    """IV bolus into a two-compartment model, solved numerically.

    Central compartment (V1): plasma and well-perfused organs. Dose goes in
        here and elimination (CL) happens from here.
    Peripheral compartment (V2): slowly equilibrating tissues (e.g. muscle, fat).
    Q (L/h): inter-compartmental clearance, how fast drug swaps between the two.

    The curve is bi-exponential: a fast early fall while drug distributes into
    tissue, then a slower terminal fall (elimination, with tissue drug leaking back).
    """
    t = np.asarray(t, dtype=float)

    def rates(_time, amounts):
        central, peripheral = amounts
        c1, c2 = central / v1, peripheral / v2
        elimination = cl * c1         # mg/h leaving the body
        distribution = q * (c1 - c2)  # net mg/h central -> peripheral (negative = flowing back)
        return [-elimination - distribution, distribution]

    central, peripheral = _solve(rates, [dose, 0.0], t)
    return central / v1


def time_grid(tau, n_doses, points_per_interval=100):
    """Time points from 0 to n_doses * tau with every dose time exactly on the grid.

    An IV bolus peaks at the exact moment of dosing, so a grid that stepped over
    the dose times would miss the peaks.
    """
    return np.arange(n_doses * points_per_interval + 1) / points_per_interval * tau


def multiple_doses(model, t, tau, n_doses, **params):
    """Repeat a single-dose model every tau hours, n_doses times.

    Every model above is linear: each rate is proportional to the amount of drug.
    That means doses don't interact, the body handles each one independently, so
    the total concentration is the sum of each dose's own curve shifted to start
    at its dosing time (the superposition principle).

    This breaks for non-linear PK, e.g. saturable (Michaelis-Menten) elimination
    as with phenytoin; that would need the ODE integrated from dose to dose.

    model:  any single-dose function above, e.g. one_compartment_oral
    params: that model's arguments, e.g. dose=500, cl=5, vd=50, ka=1.0
    """
    t = np.asarray(t, dtype=float)
    total = np.zeros_like(t)
    for i in range(n_doses):
        since_dose = t - i * tau  # time elapsed since dose i was given
        given = since_dose >= 0   # dose i only contributes after it is taken
        if given.any():
            total[given] += model(since_dose[given], **params)
    return total
