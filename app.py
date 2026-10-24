"""Interactive PK simulator. Run from the repo root: streamlit run app.py

Streamlit reruns this whole script top to bottom every time a slider moves.
"""
import numpy as np
import streamlit as st

from src.metrics import auc, cmax, half_life, steady_state_peak_trough, time_to_steady_state, tmax
from src.models import multiple_doses, one_compartment_iv, one_compartment_oral, time_grid
from src.plotting import plot_profile
from src.presets import PRESETS

MAX_DOSES = 200  # keeps the app responsive when the half-life is very long relative to tau


def load_preset():
    """Copy the chosen preset into session_state. Each slider has a key and reads
    its value from session_state under that key, so this moves all the sliders."""
    st.session_state.update(PRESETS[st.session_state.preset])


if "preset" not in st.session_state:  # first run in this browser tab
    st.session_state.preset = next(iter(PRESETS))
    load_preset()

st.set_page_config(page_title="PK simulator", layout="wide")
st.title("PK simulator")
st.caption("One-compartment model. The body is a bathtub: dosing is the tap, "
           "clearance (CL) the drain, volume of distribution (Vd) the tub size.")

with st.sidebar:
    preset = st.selectbox("Preset", list(PRESETS), key="preset", on_change=load_preset)
    st.caption(PRESETS[preset]["note"] + " Sources for real-drug values are in the README.")
    route = st.radio("Route", ["Oral", "IV bolus"], key="route", horizontal=True)
    oral = route == "Oral"
    dose = st.slider("Dose (mg)", 10, 2000, step=10, key="dose")
    cl = st.slider("Clearance, CL (L/h)", 0.1, 50.0, step=0.1, key="cl")
    vd = st.slider("Volume of distribution, Vd (L)", 5.0, 500.0, step=0.1, key="vd")
    ka = st.slider("Absorption rate constant, ka (1/h)", 0.05, 5.0, step=0.05, key="ka", disabled=not oral)
    f = st.slider("Bioavailability, F", 0.05, 1.0, step=0.05, key="f", disabled=not oral)
    tau = st.slider("Dosing interval, tau (h)", 1, 48, key="tau")
    window = st.slider("Therapeutic window (mg/L)", 0.0, 50.0, step=0.5, key="window")

if oral:
    model, params = one_compartment_oral, dict(dose=dose, cl=cl, vd=vd, ka=ka, f=f)
else:
    model, params = one_compartment_iv, dict(dose=dose, cl=cl, vd=vd)

# Simulate 7 half-lives of the slowest process, so the last dose is within 1% of
# steady state. That's elimination, unless absorption is slower (ka < k, "flip-flop").
k = cl / vd
slowest = min(k, ka) if oral else k
n_doses = int(np.ceil(max(7 * np.log(2) / slowest, 4 * tau) / tau))
reached = n_doses <= MAX_DOSES
n_doses = min(n_doses, MAX_DOSES)

t = time_grid(tau, n_doses)
single = model(t, **params)
multi = multiple_doses(model, t, tau, n_doses, **params)
peak, trough = steady_state_peak_trough(t, multi, tau)

st.pyplot(plot_profile(t, single, multi, window, (peak, trough) if reached else None))

st.subheader("Single dose")
col1, col2, col3, col4, col5 = st.columns(5)
col1.metric("Cmax", f"{cmax(single):.2f} mg/L")
col2.metric("Tmax", f"{tmax(t, single):.1f} h")
col3.metric("AUC 0-t", f"{auc(t, single):.0f} mg·h/L")
col4.metric("Elimination half-life", f"{np.log(2) / k:.1f} h", help="ln(2) x Vd / CL, from the parameters.")
col5.metric("Terminal half-life", f"{half_life(t, single):.1f} h",
            help="Fitted to the slope of ln(C) at the end of the curve, as non-compartmental analysis "
                 "does with real samples. It differs from the elimination half-life when absorption is "
                 "the slower process (ka < k, flip-flop) or when ka is close to k.")

st.subheader("Repeated dosing")
if reached:
    col1, col2, col3 = st.columns(3)
    col1.metric("Steady-state peak", f"{peak:.2f} mg/L")
    col2.metric("Steady-state trough", f"{trough:.2f} mg/L")
    col3.metric("Time to steady state (95%)", f"{time_to_steady_state(t, multi, tau):.0f} h")
else:
    st.warning(f"The half-life is very long compared with tau, so only the first {MAX_DOSES} "
               "doses are shown and steady state isn't reached yet.")
