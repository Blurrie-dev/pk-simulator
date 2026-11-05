"""A scenario: one drug, one patient, one regimen. Turns the app's settings into curves.

Kept free of Streamlit so the tests (and notebooks) can drive exactly what the
app shows.
"""
from dataclasses import dataclass, fields, replace

import numpy as np

from src import models

ROUTES = ("Oral", "IV bolus", "IV infusion")
MAX_DOSES = 200  # keeps things responsive when the half-life is very long relative to tau


@dataclass(frozen=True)
class Scenario:
    route: str = "Oral"         # one of ROUTES
    dose: float = 250.0         # mg, maintenance dose
    tau: float = 8.0            # h, dosing interval
    cl: float = 2.5             # L/h (ignored when saturable)
    vd: float = 50.0            # L; V1, the central volume, with two compartments
    ka: float = 1.0             # 1/h, oral only
    f: float = 0.8              # bioavailability, oral only
    duration: float = 1.0       # h, infusion length, IV infusion only
    two_compartment: bool = False
    q: float = 5.0              # L/h, inter-compartmental clearance
    v2: float = 50.0            # L, peripheral volume
    saturable: bool = False     # Michaelis-Menten elimination instead of CL
    vmax: float = 20.0          # mg/h, maximum elimination rate
    km: float = 4.0             # mg/L, concentration at half of vmax
    loading_dose: float = 0.0   # mg, replaces the first dose when > 0
    missed: tuple = ()          # dose numbers skipped
    double_next: bool = False   # make up a missed dose by doubling the next
    late: tuple = ()            # ((dose number, hours late), ...)

    @classmethod
    def from_preset(cls, preset, **changes):
        """Scenario from a src.presets entry (ignoring its window, note and fe)."""
        names = {f.name for f in fields(cls)}
        return cls(**{**{k: v for k, v in preset.items() if k in names}, **changes})

    @property
    def oral(self):
        return self.route == "Oral"

    @property
    def infusion(self):
        return self.route == "IV infusion"

    @property
    def regular(self):
        """Same regimen without loading dose or what-ifs: what steady state means."""
        return replace(self, loading_dose=0.0, missed=(), late=(), double_next=False)

    def schedule(self, n_doses):
        return models.dosing_schedule(self.tau, n_doses, self.dose, self.loading_dose,
                                      self.missed, dict(self.late), self.double_next)

    def linear_model(self):
        """The single-dose analytical function and its parameters (dose excluded)."""
        if self.two_compartment:
            base = dict(cl=self.cl, v1=self.vd, q=self.q, v2=self.v2)
            if self.oral:
                return models.two_compartment_oral, dict(base, ka=self.ka, f=self.f)
            if self.infusion:
                return models.two_compartment_infusion, dict(base, duration=self.duration)
            return models.two_compartment_iv, base
        base = dict(cl=self.cl, vd=self.vd)
        if self.oral:
            return models.one_compartment_oral, dict(base, ka=self.ka, f=self.f)
        if self.infusion:
            return models.one_compartment_infusion, dict(base, duration=self.duration)
        return models.one_compartment_iv, base

    def concentration(self, t, schedule, rtol=models.RTOL):
        """Concentration at times t for any list of (time, amount) dose events.
        rtol only matters for saturable elimination, the one solved numerically."""
        if self.saturable:
            return models.simulate_ode(
                t, schedule, cl=self.cl, v1=self.vd, ka=self.ka if self.oral else None, f=self.f,
                duration=self.duration if self.infusion else 0.0,
                q=self.q if self.two_compartment else 0.0, v2=self.v2, vmax=self.vmax, km=self.km,
                rtol=rtol, atol=models.ATOL * 100)
        model, params = self.linear_model()
        return models.superpose(model, t, schedule, **params)

    def single(self, t):
        return self.concentration(t, [(0.0, self.dose)])

    def repeated(self, t, n_doses, rtol=models.RTOL):
        return self.concentration(t, self.schedule(n_doses), rtol)

    def grid(self, n_doses, points_per_interval=100):
        """Time grid with every dose, and the end of every infusion, on it."""
        times = [time for time, _ in self.schedule(n_doses)]
        if self.infusion:
            times += [time + self.duration for time in times]
        return models.time_grid(self.tau, n_doses, points_per_interval, extra=times)

    # --- half-lives and how long to simulate ---------------------------------

    def average_input_rate(self):
        """mg/h reaching the circulation on average."""
        return (self.f if self.oral else 1.0) * self.dose / self.tau

    def saturated(self):
        """Saturable elimination can't keep up: input >= vmax, so there is no steady
        state and the concentration keeps rising (the bathtub drain is at full flow)."""
        return self.saturable and self.average_input_rate() >= self.vmax

    def average_steady_state(self):
        """Average steady-state concentration (mg/L), or inf if there is none.

        Linear: input = CL * Css, so Css = F * dose / (CL * tau).
        Saturable: input = vmax * Css / (km + Css), so Css = km * R / (vmax - R).
        """
        rate = self.average_input_rate()
        if not self.saturable:
            return rate / self.cl
        return np.inf if rate >= self.vmax else self.km * rate / (self.vmax - rate)

    def elimination_rate(self):
        """Terminal elimination rate constant (1/h): k for one compartment, beta for
        two. With saturable elimination it depends on concentration, so this is the
        value at the average steady state (CL there is vmax / (km + Css))."""
        if self.saturable:
            css = self.average_steady_state()
            if not np.isfinite(css):
                return 0.0
            cl = self.vmax / (self.km + css)
        else:
            cl = self.cl
        if self.two_compartment:
            (_, _alpha), (_, beta) = models.two_compartment_constants(cl, self.vd, self.q, self.v2)
            return beta
        return cl / self.vd

    def half_life(self):
        k = self.elimination_rate()
        return np.log(2) / k if k > 0 else np.inf

    def approach_rate(self):
        """Rate constant (1/h) at which the gap to steady state closes.

        Linear: the elimination rate, so each half-life halves the gap. Saturable:
        slower. Near Css, a little extra drug raises the elimination rate by only
        d/dC [vmax C / (km + C)] = vmax km / (km + Css)^2, less than the average
        vmax / (km + Css). That's why phenytoin can take weeks to settle after a dose
        change though its half-life is about a day.
        """
        if not self.saturable:
            return self.elimination_rate()
        css = self.average_steady_state()
        if not np.isfinite(css):
            return 0.0
        volume = self.vd + (self.v2 if self.two_compartment else 0.0)
        return self.vmax * self.km / (self.km + css) ** 2 / volume

    def doses_needed(self):
        """Doses to simulate so the last is within ~1% of steady state: 7 half-lives
        of the slowest process. That's elimination, unless absorption is slower (ka < k,
        "flip-flop"). Can exceed MAX_DOSES; callers cap it."""
        slowest = self.approach_rate()
        if self.oral:
            slowest = min(slowest, self.ka)
        if slowest <= 0:
            return MAX_DOSES + 1
        return int(np.ceil(max(7 * np.log(2) / slowest, 4 * self.tau) / self.tau))
