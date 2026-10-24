"""Regenerate figures/single_vs_multiple.png, the figure in the README.

Run from the repo root: python make_figure.py
"""
from pathlib import Path

from src.metrics import steady_state_peak_trough
from src.models import multiple_doses, one_compartment_oral, time_grid
from src.plotting import plot_profile

# A hypothetical drug, not a real one: half-life = ln(2) * 50 / 2.5 = 13.9 h.
params = dict(dose=250, cl=2.5, vd=50, ka=1.0, f=0.8)
tau, n_doses = 8, 13  # every 8 h for 104 h, about 7.5 half-lives

t = time_grid(tau, n_doses)
single = one_compartment_oral(t, **params)
multi = multiple_doses(one_compartment_oral, t, tau, n_doses, **params)

fig = plot_profile(t, single, multi, window=(5, 15), steady_state=steady_state_peak_trough(t, multi, tau))
fig.axes[0].set_title("Oral 250 mg, single dose vs every 8 h (hypothetical drug, half-life 13.9 h)", loc="left")
Path("figures").mkdir(exist_ok=True)
fig.savefig("figures/single_vs_multiple.png", dpi=150)
