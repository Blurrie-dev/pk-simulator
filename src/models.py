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
    is no simple formula (saturable, Michaelis-Menten elimination) and to
    cross-check every formula (the tests compare the two).

Linear models (one- and two-compartment; bolus, infusion and oral) are combined
over many doses by superposition; saturable elimination is integrated dose by
dose with simulate_ode.
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


def one_compartment_infusion(t, dose, cl, vd, duration):
    """Constant-rate IV infusion (zero-order input) over `duration` hours, analytical.

    While the tap runs, drug goes in at R = dose / duration mg/h and concentration
    climbs towards the level where input equals output, R / CL, closing a fraction
    k of the remaining gap per hour. When the tap closes it decays exponentially
    from wherever it got to:

        during: C(t) = R / CL * (1 - exp(-k*t))
        after:  C(t) = C(T) * exp(-k*(t - T))

    Both are one formula if we track the time spent infusing so far, min(t, T).
    """
    t = np.asarray(t, dtype=float)
    k, rate = cl / vd, dose / duration
    infused_for = np.minimum(t, duration)
    return rate / cl * (1 - np.exp(-k * infused_for)) * np.exp(-k * (t - infused_for))


def two_compartment_constants(cl, v1, q, v2):
    """Exponential rates alpha (fast, distribution) and beta (slow, terminal) of the
    two-compartment model, and the per-mg coefficients of each exponential.

    With micro-constants k10 = CL/V1 (elimination), k12 = Q/V1 and k21 = Q/V2
    (swapping between compartments), the concentration after a bolus is

        C(t) = dose * (A * exp(-alpha*t) + B * exp(-beta*t))

    where alpha and beta are the roots of s^2 - (k10 + k12 + k21) s + k10 k21 = 0,
    A = (alpha - k21) / (V1 (alpha - beta)) and B = (k21 - beta) / (V1 (alpha - beta)).
    Needs q > 0 and v2 > 0: with no peripheral compartment use the one-compartment model.
    """
    if q <= 0 or v2 <= 0:
        raise ValueError("two-compartment model needs q > 0 and v2 > 0")
    k10, k12, k21 = cl / v1, q / v1, q / v2
    total = k10 + k12 + k21
    root = np.sqrt(total**2 - 4 * k10 * k21)
    alpha, beta = (total + root) / 2, (total - root) / 2
    a = (alpha - k21) / (v1 * (alpha - beta))
    b = (k21 - beta) / (v1 * (alpha - beta))
    return (a, alpha), (b, beta)


def _bateman(t, ka, lam):
    """Response to an oral dose of one exponential term exp(-lam*t): the same
    "elimination minus absorption" shape as one_compartment_oral, per unit coefficient."""
    if np.isclose(ka, lam):
        return ka * t * np.exp(-lam * t)  # the ka -> lam limit, as in one_compartment_oral
    return ka / (ka - lam) * (np.exp(-lam * t) - np.exp(-ka * t))


def two_compartment_iv(t, dose, cl, v1, q, v2):
    """IV bolus into a two-compartment model, analytical (bi-exponential).

    Central compartment (V1): plasma and well-perfused organs. Dose goes in
        here and elimination (CL) happens from here.
    Peripheral compartment (V2): slowly equilibrating tissues (e.g. muscle, fat).
    Q (L/h): inter-compartmental clearance, how fast drug swaps between the two.

    Think of two bathtubs joined by a pipe: the drain is in the first one. The curve
    falls fast at first while drug flows through the pipe into tissue (alpha phase),
    then slowly as elimination takes over with tissue drug leaking back (beta phase).
    """
    t = np.asarray(t, dtype=float)
    return dose * sum(coef * np.exp(-lam * t) for coef, lam in two_compartment_constants(cl, v1, q, v2))


def two_compartment_oral(t, dose, cl, v1, q, v2, ka, f=1.0):
    """Oral dose with first-order absorption into a two-compartment model, analytical.

    Linear model, so each exponential of the bolus response gets the same
    absorption treatment as the one-compartment Bateman function, and the
    results add up.
    """
    t = np.asarray(t, dtype=float)
    return f * dose * sum(coef * _bateman(t, ka, lam) for coef, lam in two_compartment_constants(cl, v1, q, v2))


def two_compartment_infusion(t, dose, cl, v1, q, v2, duration):
    """Constant-rate IV infusion into a two-compartment model, analytical.

    An infusion is a stream of tiny boluses, so each exponential term builds up
    as (1 - exp(-lam*t)) / lam while it runs and decays after, as in
    one_compartment_infusion.
    """
    t = np.asarray(t, dtype=float)
    rate = dose / duration
    infused_for = np.minimum(t, duration)
    return rate * sum(coef / lam * (1 - np.exp(-lam * infused_for)) * np.exp(-lam * (t - infused_for))
                      for coef, lam in two_compartment_constants(cl, v1, q, v2))


def two_compartment_ode(t, dose, cl, v1, q, v2, ka=None, f=1.0, duration=0.0):
    """The two-compartment model solved numerically, to check the formulas above.

    ka given: oral. duration > 0: infusion. Otherwise: IV bolus.
    """
    return simulate_ode(t, [(0.0, dose)], cl=cl, v1=v1, q=q, v2=v2, ka=ka, f=f, duration=duration)


def simulate_ode(t, schedule, cl, v1, ka=None, f=1.0, duration=0.0, q=0.0, v2=1.0,
                 vmax=None, km=None, rtol=RTOL, atol=ATOL):
    """General model solved numerically, dose by dose. Works when superposition doesn't.

    schedule: list of (time in h, amount in mg) dose events
    ka:       given -> oral (F * amount goes into the gut); None -> IV
    duration: > 0 -> each IV dose runs as an infusion over this many hours
    q, v2:    > 0 -> adds a peripheral compartment (two-compartment model)
    vmax, km: given -> saturable, Michaelis-Menten elimination: vmax (mg/h) *
              C / (km + C) instead of CL * C. At low C it behaves like CL = vmax / km;
              at high C it can't go faster than vmax, like a checkout with a fixed
              number of tills: double the queue doesn't double the throughput.

    Between dose events the rates are smooth, so solve_ivp integrates each stretch
    on its own; at a bolus we add the dose to the right compartment and carry on.
    Infusion starts and stops are also stretch boundaries, so the step in the
    input rate never falls inside a solver step.
    """
    t = np.asarray(t, dtype=float)
    t_end = t[-1]
    oral, infusion = ka is not None, duration > 0 and ka is None
    saturable = vmax is not None
    events = sorted((float(time), float(amount)) for time, amount in schedule if amount > 0)

    breaks = {0.0, t_end, *(time for time, _ in events)}
    if infusion:
        breaks |= {time + duration for time, _ in events}
    breaks = sorted(b for b in breaks if 0 <= b <= t_end)

    def rates(_time, amounts, infusion_rate):
        gut, central, peripheral = amounts
        c1 = central / v1
        absorption = ka * gut if oral else 0.0
        elimination = vmax * c1 / (km + c1) if saturable else cl * c1
        distribution = q * (c1 - peripheral / v2) if q > 0 else 0.0
        return [-absorption, infusion_rate + absorption - elimination - distribution, distribution]

    def give_boluses(at, state):
        if infusion:
            return
        for time, amount in events:
            if time == at:
                if oral:
                    state[0] += f * amount
                else:
                    state[1] += amount

    state = np.zeros(3)
    out = np.zeros((3, t.size))
    for start, stop in zip(breaks[:-1], breaks[1:]):
        give_boluses(start, state)
        rate = sum(amount / duration for time, amount in events
                   if infusion and time <= start < time + duration)
        inside = (t >= start) & (t < stop)
        t_eval = np.append(t[inside], stop)
        sol = solve_ivp(rates, (start, stop), state, t_eval=t_eval, args=(rate,), rtol=rtol, atol=atol)
        if not sol.success:
            raise RuntimeError(sol.message)
        out[:, inside] = sol.y[:, :-1]
        state = sol.y[:, -1].copy()
    give_boluses(t_end, state)  # a dose exactly at the last time point
    out[:, t == t_end] = state[:, None]
    return out[1] / v1


def dosing_schedule(tau, n_doses, dose, loading_dose=0.0, missed=(), late=None, double_next=False):
    """Dose events (time h, amount mg) for a regimen, with optional what-ifs.

    loading_dose: > 0 replaces the first dose with this amount
    missed:       dose numbers (1 = first) that are skipped
    double_next:  a skipped dose is made up by doubling the next one
    late:         {dose number: hours late}
    """
    times = [i * tau for i in range(n_doses)]
    amounts = [float(dose)] * n_doses
    if loading_dose > 0 and n_doses:
        amounts[0] = float(loading_dose)
    for number in missed:
        if 1 <= number <= n_doses:
            if double_next and number < n_doses:
                amounts[number] += amounts[number - 1]
            amounts[number - 1] = 0.0
    for number, hours in (late or {}).items():
        if 1 <= number <= n_doses:
            times[number - 1] += hours
    return [(time, amount) for time, amount in zip(times, amounts) if amount > 0]


def time_grid(tau, n_doses, points_per_interval=100, extra=()):
    """Time points from 0 to n_doses * tau with every dose time exactly on the grid.

    An IV bolus peaks at the exact moment of dosing, so a grid that stepped over
    the dose times would miss the peaks. `extra` adds more times that must be on
    the grid, such as late doses or the end of each infusion (where its peak is).
    """
    grid = np.arange(n_doses * points_per_interval + 1) / points_per_interval * tau
    extra = [e for e in extra if 0 <= e <= grid[-1]]
    return np.union1d(grid, extra) if extra else grid


def superpose(model, t, schedule, **params):
    """Add up single-dose curves, one per dose event, each shifted to its dosing time.

    Every linear model above has rates proportional to the amount of drug, so doses
    don't interact: the body handles each one independently, and the total is the
    sum of each dose's own curve (the superposition principle). That's why a dose
    can have any size or time: loading doses, missed or late doses all just work.

    This breaks for non-linear PK, e.g. saturable (Michaelis-Menten) elimination
    as with phenytoin; use simulate_ode for that.

    model:    any single-dose function above, e.g. one_compartment_oral
    schedule: list of (time h, amount mg)
    params:   that model's other arguments, e.g. cl=5, vd=50, ka=1.0
    """
    t = np.asarray(t, dtype=float)
    total = np.zeros_like(t)
    for time, amount in schedule:
        since_dose = t - time     # time elapsed since this dose
        given = since_dose >= 0   # a dose only contributes after it is taken
        if amount > 0 and given.any():
            total[given] += model(since_dose[given], dose=amount, **params)
    return total


def multiple_doses(model, t, tau, n_doses, dose, **params):
    """Repeat a single-dose linear model every tau hours, n_doses times, by superposition.

    model:  any single-dose function above, e.g. one_compartment_oral
    params: that model's other arguments, e.g. cl=5, vd=50, ka=1.0
    """
    return superpose(model, t, dosing_schedule(tau, n_doses, dose), **params)
