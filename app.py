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


SLIDERS = ["dose", "cl", "vd", "ka", "f", "tau"]


def load_preset():
    """Copy the chosen preset into session_state. Each widget has a key and reads
    its value from session_state under that key, so this moves all the inputs."""
    p = PRESETS[st.session_state.preset]
    st.session_state.route = p["route"]
    st.session_state.window_lo, st.session_state.window_hi = p["window"]
    for key in SLIDERS:
        if key in p:  # IV presets have no ka or F
            st.session_state[key] = st.session_state[key + "_box"] = p[key]


def copy_value(source, target):
    st.session_state[target] = st.session_state[source]


def slider_with_box(label, key, lo, hi, step, disabled=False):
    """A slider plus a box for typing an exact value. Each widget keeps its own
    session_state key; whichever one changes copies its value into the other."""
    slider_col, box_col = st.columns([3, 2], vertical_alignment="bottom")
    slider_col.slider(label, lo, hi, step=step, key=key, disabled=disabled,
                      on_change=copy_value, args=(key, key + "_box"))
    box_col.number_input(label, lo, hi, step=step, key=key + "_box", disabled=disabled,
                         label_visibility="collapsed", on_change=copy_value, args=(key + "_box", key))
    return st.session_state[key]


if "preset" not in st.session_state:  # first run in this browser tab
    st.session_state.preset = next(iter(PRESETS))
    load_preset()

st.set_page_config(page_title="PK simulator", layout="wide")
st.title("PK simulator")
st.caption("One-compartment model. The body is a bathtub: dosing is the tap, "
           "clearance (CL) the drain, volume of distribution (Vd) the tub size.")

with st.sidebar:
    st.subheader("Drug", divider="gray")
    options = sorted(PRESETS, key=lambda name: (name != "Hypothetical drug", name))  # A-Z after the default
    preset = st.selectbox("Preset", options, key="preset", on_change=load_preset)
    st.caption(PRESETS[preset]["note"] + " Sources for real-drug values are in the README.")
    route = st.radio("Route", ["Oral", "IV bolus"], key="route", horizontal=True)
    oral = route == "Oral"

    st.subheader("Dosing", divider="gray")
    dose = slider_with_box("Dose (mg)", "dose", 10, 2000, step=10)
    tau = slider_with_box("Dosing interval, tau (h)", "tau", 1, 48, step=1)

    st.subheader("Pharmacokinetics", divider="gray")
    cl = slider_with_box("Clearance, CL (L/h)", "cl", 0.1, 50.0, step=0.1)
    vd = slider_with_box("Volume of distribution, Vd (L)", "vd", 5.0, 500.0, step=0.1)
    ka = slider_with_box("Absorption rate constant, ka (1/h)", "ka", 0.05, 5.0, step=0.05, disabled=not oral)
    f = slider_with_box("Bioavailability, F", "f", 0.05, 1.0, step=0.05, disabled=not oral)

    st.subheader("Therapeutic window (mg/L)", divider="gray")
    low_col, high_col = st.columns(2)
    window = (low_col.number_input("Lower", 0.0, 1000.0, step=0.5, key="window_lo"),
              high_col.number_input("Upper", 0.0, 1000.0, step=0.5, key="window_hi"))

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

with st.container(border=True):
    st.pyplot(plot_profile(t, single, multi, window, (peak, trough) if reached else None))

# border=True draws each metric as a card; units sit in the labels so values fit
st.subheader("Single dose")
col1, col2, col3, col4, col5 = st.columns(5)
col1.metric("Cmax (mg/L)", f"{cmax(single):.2f}", border=True)
col2.metric("Tmax (h)", f"{tmax(t, single):.1f}", border=True)
col3.metric("AUC 0-t (mg·h/L)", f"{auc(t, single):.0f}", border=True)
col4.metric("t½ elimination (h)", f"{np.log(2) / k:.1f}", border=True,
            help="ln(2) x Vd / CL, from the parameters.")
col5.metric("t½ terminal (h)", f"{half_life(t, single):.1f}", border=True,
            help="Fitted to the slope of ln(C) at the end of the curve, as non-compartmental analysis "
                 "does with real samples. It differs from the elimination half-life when absorption is "
                 "the slower process (ka < k, flip-flop) or when ka is close to k.")

st.subheader("Repeated dosing")
if reached:
    col1, col2, col3 = st.columns(3)
    col1.metric("Steady-state peak (mg/L)", f"{peak:.2f}", border=True)
    col2.metric("Steady-state trough (mg/L)", f"{trough:.2f}", border=True)
    col3.metric("Time to steady state, 95% (h)", f"{time_to_steady_state(t, multi, tau):.0f}", border=True)
else:
    st.warning(f"The half-life is very long compared with tau, so only the first {MAX_DOSES} "
               "doses are shown and steady state isn't reached yet.")
