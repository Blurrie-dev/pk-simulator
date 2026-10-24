# pk-simulator

Interactive pharmacokinetic simulator: plasma concentration over time for different doses,
routes and dosing intervals.

Think of the body as a bathtub: dosing is the tap, clearance is the drain, and volume of
distribution is the tub size.

![Single vs repeated oral dosing](figures/single_vs_multiple.png)

A single dose never reaches the therapeutic window. Dosing every 8 h, each dose lands on what's
left of the previous ones, so concentration builds up until the amount eliminated per interval
equals the amount given. That's steady state, about 4 to 5 half-lives in, whatever the dose.

## Models

`src/models.py`, units mg, L, h (concentration in mg/L):

- One-compartment IV bolus: `C(t) = Dose/Vd * exp(-k*t)` with `k = CL/Vd`
- One-compartment oral, first-order absorption `ka` and bioavailability `F` (Bateman function)
- The same one-compartment model solved numerically with `scipy.integrate.solve_ivp`, as a check
- Two-compartment IV (central + peripheral), ODE system solved with `solve_ivp`
- Multiple dosing every `tau` hours, by superposition of single doses (valid for linear PK)

`src/metrics.py`, computed from the simulated curve: Cmax, Tmax, AUC (trapezoidal), terminal
half-life (log-linear fit), steady-state peak and trough, time to steady state (95% of the
steady-state trough).

## Setup

```bash
conda env create -f environment.yml
conda activate pk-simulator
```

or `pip install -r requirements.txt`.

## Run

```bash
streamlit run app.py      # interactive app
pytest                    # checks against analytical results
python make_figure.py     # regenerate the figure above
```

The tests check that:

- the numerical one-compartment solution matches the analytical formula within 1% (IV, oral,
  ka = k and flip-flop)
- AUC = F * Dose / CL
- steady state takes 4 to 5 half-lives for doses of 50, 500 and 5000 mg
- each drug preset reproduces what its sources observed (below)

## Deploy

On Streamlit Community Cloud, point a new app at `app.py`. When a repo has both files, Community
Cloud installs from `environment.yml` (conda) rather than `requirements.txt`.

## Structure

```
src/models.py      PK models
src/metrics.py     Cmax, Tmax, AUC, half-life, time to steady state
src/plotting.py    concentration-time plot used by the app and the figure
src/presets.py     drug presets for the app
tests/test_pk.py   pytest checks
app.py             Streamlit app
make_figure.py     builds figures/single_vs_multiple.png
notebooks/         exploration
```

## Drug presets

The app has a preset menu: a hypothetical drug (also used for the figure above) and two real
drugs with published values. Per-kg values are scaled to a 70 kg adult and rounded to 0.1 to fit
the slider steps. 1 mcg/mL = 1 mg/L. These are population averages for teaching, not for clinical
use.

### Theophylline, oral

| Parameter | In the app | Published value | Source |
|---|---|---|---|
| CL | 2.7 L/h | 0.65 mL/min/kg, otherwise healthy non-smoking adults 16-60 y (x 70 kg = 2.73 L/h) | [1] Table I |
| Vd | 31.5 L | 0.45 L/kg (range 0.3-0.7), based on ideal body weight | [1] Distribution |
| F | 1.0 | "rapidly and completely absorbed" from oral solution | [1] Absorption |
| ka | 2.0 /h | Not reported. Chosen so the model's Tmax (1.6 h) falls in the reported 1-2 h to peak | derived from [1] Absorption |
| Dose, tau | 200 mg every 8 h | 600 mg/day divided every 6-8 h, the final adult titration step | [1] Table V |
| Window | 10-20 mg/L | therapeutic range 10-20 mcg/mL | [1] |

Check (in `tests/test_pk.py`): the model's half-life is 8.1 h against 8.7 h (range 6.1-12.8) in
[1], and a single 5 mg/kg dose peaks at 9.7 mg/L at 1.6 h against "about 10 mcg/mL (range 5-15)
... 1-2 hr after the dose" in [1].

### Gentamicin, IV

| Parameter | In the app | Published value | Source |
|---|---|---|---|
| CL | 5.4 L/h | 1.29 mL/min/kg, patients 20-39 y with normal renal function | [3] |
| Vd | 16.1 L | 0.23 L/kg | [3] |
| Dose, tau | 70 mg every 8 h | 3 mg/kg/day in three equal doses every 8 h | [2] Dosage and Administration |
| Window | 4-12 mg/L | peak expected at 4-6 mcg/mL; avoid prolonged peaks above 12 mcg/mL and troughs above 2 mcg/mL | [2] |

For gentamicin the band shows the peak targets: peaks should reach it and troughs should fall below
2 mg/L, unlike theophylline where the whole curve should stay inside. In practice it is given IM
or as an IV infusion, so the bolus model slightly overstates the peak.

Check (in `tests/test_pk.py`): the model's half-life is 2.1 h against 2.2 h in [3], and at steady
state the peak is 4.7 mg/L and the trough 0.3 mg/L, inside the label's expectations.

### References

1. Elixophyllin (theophylline anhydrous) oral solution, prescribing information. Nostrum
   Laboratories, Inc. DailyMed, version of 14 Dec 2018.
   https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=b488f996-2f80-438a-8333-b837f9acd4fc
2. Gentamicin sulfate injection, prescribing information. Fresenius Kabi USA, LLC. DailyMed,
   version of 17 Oct 2024.
   https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=84f5c763-1cd3-4d85-9afb-934db8666fbf
3. Bauer LA, Blouin RA. Gentamicin pharmacokinetics: effect of aging in patients with normal renal
   function. J Am Geriatr Soc. 1982;30(5):309-311. https://doi.org/10.1111/j.1532-5415.1982.tb05618.x

## What I found

_To be written._
