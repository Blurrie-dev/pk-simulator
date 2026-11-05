"""Interactive PK simulator. Run from the repo root: streamlit run app.py

Streamlit reruns this whole script top to bottom every time an input changes.
Inputs live in the sidebar; the main area has one tab per tool.
"""
from dataclasses import asdict

import numpy as np
import pandas as pd
import streamlit as st

from src import charts, dosing, patient, population, settings, tdm
from src import pd as pharmacodynamics
from src.metrics import (auc, average_steady_state, cmax, half_life, steady_state_peak_trough,
                         time_in_window, time_to_steady_state, tmax)
from src.models import two_compartment_constants
from src.presets import PRESETS
from src.simulate import MAX_DOSES, ROUTES, Scenario

PAIRED = ["dose", "tau", "cl", "vd", "ka", "f", "duration", "q", "v2", "vmax", "km"]  # slider + box
DEFAULTS = Scenario()  # generic values for inputs a preset doesn't set
MAX_COMPARISONS = 4


# --- session state helpers ----------------------------------------------------------

def set_value(key, value):
    """Set an input and, for slider + box pairs, its box too."""
    st.session_state[key] = value
    if key in PAIRED:
        st.session_state[key + "_box"] = value


def load_preset():
    """Copy the chosen preset into session_state. Each widget has a key and reads
    its value from session_state under that key, so this moves all the inputs."""
    p = PRESETS[st.session_state.preset]
    st.session_state.route = p["route"]
    st.session_state.window_lo, st.session_state.window_hi = p["window"]
    st.session_state.elimination = "Saturable" if p.get("saturable") else "Linear"
    st.session_state.compartments = "Two" if p.get("two_compartment") else "One"
    st.session_state.fe = p.get("fe", 0.5)
    st.session_state.loading = 0
    st.session_state.missed = st.session_state.late_dose = 0
    for key in PAIRED:
        value = p.get(key, getattr(DEFAULTS, key))
        set_value(key, round(value * 24, 1) if key == "vmax" else value)  # vmax shown per day


def copy_value(source, target):
    st.session_state[target] = st.session_state[source]


def slider_with_box(label, key, lo, hi, step, disabled=False, help=None):
    """A slider plus a box for typing an exact value. Each widget keeps its own
    session_state key; whichever one changes copies its value into the other."""
    slider_col, box_col = st.columns([3, 2], vertical_alignment="bottom")
    slider_col.slider(label, lo, hi, step=step, key=key, disabled=disabled, help=help,
                      on_change=copy_value, args=(key, key + "_box"))
    box_col.number_input(label, lo, hi, step=step, key=key + "_box", disabled=disabled,
                         label_visibility="collapsed", on_change=copy_value, args=(key + "_box", key))
    return st.session_state[key]


def apply_regimen(dose, tau):
    set_value("dose", int(min(max(round(dose), 10), 2000)))
    set_value("tau", int(tau))


def apply_fit(factors):
    """Scale the sidebar values by fitted / prior. A ratio works whether or not the
    patient panel is adjusting them, and vmax's per-day units cancel."""
    for key, factor in factors.items():
        set_value(key, settings.clamp(key, st.session_state[key] * factor))


def apply_settings(values):
    for key, value in values.items():
        set_value(key, value)


def load_file():
    upload = st.session_state.get("upload")
    if upload is None:
        return
    try:
        values, rejected = settings.from_json(upload.getvalue().decode())
    except ValueError:
        st.session_state.load_message = ("error", "That file isn't a saved scenario (JSON).")
        return
    apply_settings(values)
    st.session_state.load_message = ("success", f"Loaded {len(values)} settings"
                                     + (f"; ignored {', '.join(rejected)}" if rejected else "") + ".")


def share_link():
    st.query_params.clear()
    st.query_params.update(settings.to_query(st.session_state))


def fill_example_levels(scenario):
    """Two levels from a made-up patient who clears the drug 40% faster than typical,
    measured around the 4th dose: a peak-ish level and a trough."""
    from dataclasses import replace as _replace
    truth = _replace(scenario, cl=scenario.cl * 1.4, vmax=scenario.vmax * 1.4)
    times = [3 * scenario.tau + min(2.0, scenario.tau / 2), 4 * scenario.tau - 0.25]
    levels = tdm.predict(truth, times)
    st.session_state.levels_data = pd.DataFrame({"time_h": times, "level_mg_L": np.round(levels, 2)})
    st.session_state.levels_version += 1


def pin(label, t, c):
    st.session_state.comparisons = (st.session_state.comparisons + [(label, t, c)])[-MAX_COMPARISONS:]


# --- cached simulations (keyed on plain dicts so Streamlit can hash them) -----------

@st.cache_data(max_entries=64)
def simulate_profile(scenario_dict, n_doses):
    s = Scenario(**scenario_dict)
    t = s.grid(n_doses)
    rtol = 1e-6  # only used by saturable (numerical) models
    return t, s.concentration(t, [(0.0, s.dose)], rtol), s.repeated(t, n_doses, rtol)


@st.cache_data(max_entries=16, show_spinner="Simulating patients...")
def simulate_population(scenario_dict, cv, n_patients, seed, n_doses):
    s = Scenario(**scenario_dict)
    t = s.grid(n_doses, points_per_interval=50)
    patients = population.sample_patients(s, dict(cv), n_patients, seed)
    curves = population.simulate(patients, t, n_doses)
    return t, curves, s.repeated(t, n_doses)


if "preset" not in st.session_state:  # first run in this browser tab
    st.session_state.preset = next(iter(PRESETS))
    st.session_state.comparisons = []
    st.session_state.update(weight=70.0, age=40, sex="Male", scr=1.0, late_hours=4.0, adjust=False,
                            double_next=False, levels_version=0,
                            levels_data=pd.DataFrame({"time_h": [], "level_mg_L": []}))
    load_preset()
    shared, _rejected = settings.parse(st.query_params.to_dict())  # a shared link
    if "preset" in shared:
        st.session_state.preset = shared["preset"]
        load_preset()
    apply_settings(shared)

st.set_page_config(page_title="PK simulator", layout="wide")


# --- sidebar: every input -----------------------------------------------------------

with st.sidebar:
    st.subheader("Drug", divider="gray")
    options = sorted(PRESETS, key=lambda name: (name != "Hypothetical drug", name))  # A-Z after the default
    preset = st.selectbox("Preset", options, key="preset", on_change=load_preset)
    st.caption(PRESETS[preset]["note"] + " Sources for real-drug values are in the README.")

    st.subheader("Model", divider="gray")
    route = st.radio("Route", ROUTES, key="route", horizontal=True)
    compartments = st.radio("Compartments", ["One", "Two"], key="compartments", horizontal=True,
                            help="Two adds a slowly equilibrating tissue compartment: a second bathtub "
                                 "joined to the first by a pipe.")
    elimination = st.radio("Elimination", ["Linear", "Saturable"], key="elimination", horizontal=True,
                           help="Saturable = Michaelis-Menten: like a checkout with a fixed number of tills, "
                                "elimination can't go faster than Vmax however much drug is waiting.")
    oral, infusion = route == "Oral", route == "IV infusion"
    two, saturable = compartments == "Two", elimination == "Saturable"

    st.subheader("Dosing", divider="gray")
    dose = slider_with_box("Dose (mg)", "dose", 10, 2000, step=10)
    tau = slider_with_box("Dosing interval, tau (h)", "tau", 1, 48, step=1)
    if infusion:
        duration = slider_with_box("Infusion duration (h)", "duration", 0.1, 24.0, step=0.1)
        if duration > tau:
            st.warning("The infusion is longer than the interval; it is cut to tau.")
            duration = float(tau)
    else:
        duration = st.session_state.duration
    loading = st.number_input("Loading dose (mg, 0 = none)", 0, 20000, step=10, key="loading",
                              help="Replaces the first dose. The Regimen finder tab suggests one.")

    with st.expander("What if a dose is missed or late?"):
        missed = st.number_input("Miss dose number (0 = none)", 0, MAX_DOSES, key="missed")
        double_next = st.checkbox("...and double the next one to make up", key="double_next")
        late_dose = st.number_input("Take dose number late (0 = none)", 0, MAX_DOSES, key="late_dose")
        late_hours = st.slider("...by this many hours", 0.5, 24.0, step=0.5, key="late_hours")

    st.subheader("Pharmacokinetics", divider="gray")
    st.caption("Values for a typical 70 kg adult. Adjust for a patient below.")
    cl = slider_with_box("Clearance, CL (L/h)", "cl", 0.1, 50.0, step=0.1, disabled=saturable)
    vd = slider_with_box("Central volume, V1 (L)" if two else "Volume of distribution, Vd (L)",
                         "vd", 1.0, 500.0, step=0.1)
    ka = slider_with_box("Absorption rate constant, ka (1/h)", "ka", 0.05, 5.0, step=0.05, disabled=not oral)
    f = slider_with_box("Bioavailability, F", "f", 0.05, 1.0, step=0.01, disabled=not oral)
    if two:
        q = slider_with_box("Inter-compartmental clearance, Q (L/h)", "q", 0.1, 100.0, step=0.1)
        v2 = slider_with_box("Peripheral volume, V2 (L)", "v2", 1.0, 1000.0, step=1.0)
    else:
        q, v2 = st.session_state.q, st.session_state.v2
    if saturable:
        vmax_day = slider_with_box("Vmax (mg/day)", "vmax", 10.0, 3000.0, step=10.0,
                                   help="Maximum elimination rate, when the enzymes are saturated.")
        km = slider_with_box("Km (mg/L)", "km", 0.1, 100.0, step=0.1,
                             help="Concentration at which elimination runs at half of Vmax.")
    else:
        vmax_day, km = st.session_state.vmax, st.session_state.km

    st.subheader("Patient", divider="gray")
    adjust = st.toggle("Adjust for a specific patient", key="adjust")
    if adjust:
        weight = st.slider("Weight (kg)", 3.0, 250.0, step=0.5, key="weight")
        age = st.slider("Age (years)", 18, 100, key="age")
        female = st.radio("Sex", ["Male", "Female"], horizontal=True, key="sex") == "Female"
        scr = st.slider("Serum creatinine (mg/dL)", 0.2, 10.0, step=0.1, key="scr")
        fe = st.slider("Fraction cleared by the kidneys, fe", 0.0, 1.0, step=0.01, key="fe",
                       help="Share of CL that follows kidney function; set by the preset.")
        crcl = patient.cockcroft_gault(age, weight, scr, female)
        st.caption(f"Creatinine clearance (Cockcroft-Gault): **{crcl:.0f} mL/min** "
                   f"(reference patient {patient.REFERENCE_CRCL:.0f}).")

    st.subheader("Therapeutic window (mg/L)", divider="gray")
    low_col, high_col = st.columns(2)
    window = (low_col.number_input("Lower", 0.0, 1000.0, step=0.5, key="window_lo"),
              high_col.number_input("Upper", 0.0, 1000.0, step=0.5, key="window_hi"))
    has_window = window[1] > window[0]

    st.subheader("Save and share", divider="gray")
    st.download_button("Save scenario (JSON)", settings.to_json(st.session_state), "pk_scenario.json",
                       "application/json", icon=":material/save:", width="stretch")
    st.file_uploader("Load a saved scenario", type="json", key="upload", on_change=load_file)
    if "load_message" in st.session_state:
        kind, text = st.session_state.pop("load_message")
        getattr(st, kind)(text)
    st.button("Put this scenario in the link", on_click=share_link, icon=":material/link:", width="stretch",
              help="Updates the address bar; copy it to share exactly these settings.")


# --- build the scenario -------------------------------------------------------------

typical = Scenario(
    route=route, dose=float(dose), tau=float(tau), cl=cl, vd=vd, ka=ka, f=f, duration=duration,
    two_compartment=two, q=q, v2=v2, saturable=saturable, vmax=vmax_day / 24, km=km,
    loading_dose=float(loading), missed=(missed,) if missed else (), double_next=double_next,
    late=((late_dose, late_hours),) if late_dose else (),
)
scenario = patient.adjust(typical, weight, crcl, fe) if adjust else typical
what_if = bool(scenario.missed or scenario.late or scenario.loading_dose)

n_needed = scenario.doses_needed()
reached = n_needed <= MAX_DOSES and not scenario.saturated()
n_doses = min(n_needed, MAX_DOSES if not saturable else 60)  # the numerical model is slower
t, single, multi = simulate_profile(asdict(scenario), n_doses)
regular = simulate_profile(asdict(scenario.regular), n_doses)[2] if what_if else multi
peak, trough = steady_state_peak_trough(t, regular, scenario.tau)


# --- main area ----------------------------------------------------------------------

st.title("PK simulator")
st.caption("The body is a bathtub: dosing is the tap, clearance (CL) the drain, volume of "
           "distribution (Vd) the tub size.")

if adjust:
    st.info(f"Adjusted for this patient: CL {typical.cl:.2f} -> **{scenario.cl:.2f} L/h**, "
            f"V {typical.vd:.1f} -> **{scenario.vd:.1f} L**"
            + (f", Vmax {typical.vmax * 24:.0f} -> **{scenario.vmax * 24:.0f} mg/day**" if saturable else "")
            + ".", icon=":material/person:")
if scenario.saturated():
    st.error(f"Input ({scenario.average_input_rate() * 24:.0f} mg/day reaching the blood) is at or above "
             f"Vmax ({scenario.vmax * 24:.0f} mg/day): elimination is maxed out, so concentrations keep "
             "rising with every dose. There is no steady state.")

profile_tab, population_tab, finder_tab, fit_tab, effect_tab, about_tab = st.tabs(
    ["Concentration-time", "Population variability", "Regimen finder", "Fit to levels", "Effect (PD)",
     "About the models"])

with profile_tab:
    opt1, opt2, opt3 = st.columns(3)
    log = opt1.toggle("Log scale", key="log", help="On a log scale, first-order elimination is a straight "
                                                   "line; its slope gives the half-life.")
    show_single = opt2.toggle("Show single dose", value=True, key="show_single")
    show_doses = opt3.toggle("Show dose times", value=False, key="show_doses")

    comparisons = list(st.session_state.comparisons)
    if what_if:
        comparisons.append(("Regular schedule", t, regular))
    dose_times = [time for time, _ in scenario.schedule(n_doses)] if show_doses else ()
    with st.container(border=True):
        st.plotly_chart(charts.profile(t, single, multi, window, (peak, trough) if reached else None,
                                       comparisons, dose_times, log, show_single),
                        width="stretch", key="profile_chart")

    b1, b2, b3 = st.columns(3)
    label = (f"{preset}: {dose:g} mg every {tau} h, {route.lower()}"
             + (", 2-cpt" if two else "") + (", saturable" if saturable else "") + (", patient" if adjust else ""))
    b1.button("Pin this curve for comparison", on_click=pin, args=(label, t, multi),
              icon=":material/push_pin:", width="stretch")
    b2.button("Clear pinned curves", on_click=lambda: st.session_state.update(comparisons=[]),
              disabled=not st.session_state.comparisons, width="stretch")
    table = pd.DataFrame({"time_h": t, "single_dose_mg_per_L": single, "repeated_mg_per_L": multi})
    b3.download_button("Download curve (CSV)", table.to_csv(index=False).encode(), "pk_profile.csv",
                       "text/csv", icon=":material/download:", width="stretch")

    st.subheader("Single dose")
    col1, col2, col3, col4, col5 = st.columns(5)
    col1.metric("Cmax (mg/L)", f"{cmax(single):.3g}", border=True)
    col2.metric("Tmax (h)", f"{tmax(t, single):.1f}", border=True)
    col3.metric("AUC 0-t (mg·h/L)", f"{auc(t, single):.4g}", border=True)
    t_half_help = ("ln(2) / beta, the terminal (slow) phase of the two-compartment model." if two else
                   "ln(2) x Vd / CL, from the parameters.")
    if saturable:
        t_half_help = ("Depends on concentration with saturable elimination; this is the value at the "
                       "average steady state, ln(2) x V x (Km + Css) / Vmax.")
    col4.metric("t½ elimination (h)", f"{scenario.half_life():.1f}" if np.isfinite(scenario.half_life()) else "∞",
                border=True, help=t_half_help)
    try:
        fitted = f"{half_life(t, single):.1f}"
    except (TypeError, ValueError, np.linalg.LinAlgError):
        fitted = "n/a"
    col5.metric("t½ terminal, fitted (h)", fitted, border=True,
                help="Fitted to the slope of ln(C) at the end of the curve, as non-compartmental analysis "
                     "does with real samples. It differs from the elimination half-life when absorption is "
                     "the slower process (ka < k, flip-flop) or when ka is close to k.")
    if two:
        (_, alpha), (_, beta) = two_compartment_constants(scenario.cl, scenario.vd, scenario.q, scenario.v2)
        st.caption(f"Two compartments: distribution half-life {np.log(2) / alpha:.2f} h, "
                   f"terminal half-life {np.log(2) / beta:.1f} h.")

    st.subheader("Repeated dosing, at steady state" + (" (regular schedule)" if what_if else ""))
    if reached:
        cols = st.columns(5 if has_window else 4)
        cols[0].metric("Peak (mg/L)", f"{peak:.3g}", border=True)
        cols[1].metric("Trough (mg/L)", f"{trough:.3g}", border=True)
        cols[2].metric("Average (mg/L)", f"{average_steady_state(t, regular, scenario.tau):.3g}", border=True,
                       help="AUC over one interval at steady state, divided by tau.")
        cols[3].metric("Time to reach, 95% (h)", f"{time_to_steady_state(t, regular, scenario.tau):.0f}",
                       border=True)
        if has_window:
            cols[4].metric("Time in window", f"{time_in_window(t, regular, scenario.tau, window):.0%}",
                           border=True, help="Share of a steady-state dosing interval spent inside the window.")
    elif not scenario.saturated():
        st.warning(f"The half-life is very long compared with tau, so only the first {n_doses} "
                   "doses are shown and steady state isn't reached yet.")

with population_tab:
    st.markdown("No two patients are alike. Peaks and troughs below are at steady state. Each simulated patient gets their own CL, V and ka, drawn "
                "around the values in the sidebar (log-normal, so the sidebar value is the median "
                "patient). The bands show where 50% and 90% of patients' concentrations fall.")
    c1, c2, c3, c4 = st.columns(4)
    cv = {"cl": c1.slider("CL variability (CV %)", 0, 100, 30, 5, key="cv_cl") / 100,
          "vd": c2.slider("V variability (CV %)", 0, 100, 20, 5, key="cv_vd") / 100}
    if oral:
        cv["ka"] = c3.slider("ka variability (CV %)", 0, 150, 40, 5, key="cv_ka") / 100
    if saturable:
        cv["vmax"] = c4.slider("Vmax variability (CV %)", 0, 100, 30, 5, key="cv_vmax") / 100
    sizes = [20, 50] if saturable else [100, 250, 500, 1000]
    p1, p2 = st.columns(2)
    n_patients = p1.select_slider("Patients", sizes, value=sizes[-1] if saturable else 250, key="n_patients")
    seed = p2.number_input("Random seed", 0, 10_000, 0, key="seed", help="Same seed, same patients.")
    pop_log = st.toggle("Log scale", key="pop_log")

    if scenario.saturated():
        st.info("No steady state for the typical patient; reduce the dose first.")
    else:
        pt, curves, typical_curve = simulate_population(asdict(scenario.regular), tuple(sorted(cv.items())),
                                                        n_patients, seed, n_doses)
        bands = population.percentiles(curves)
        with st.container(border=True):
            st.plotly_chart(charts.population(pt, bands, typical_curve, window, pop_log),
                            width="stretch", key="population_chart")
        peaks, troughs = population.steady_state_summary(curves, pt, scenario.tau)
        m = st.columns(4 if has_window else 3)
        m[0].metric("Median peak (mg/L)", f"{np.median(peaks):.3g}",
                    help=f"90% of patients between {np.percentile(peaks, 5):.3g} and {np.percentile(peaks, 95):.3g}",
                    border=True)
        m[1].metric("Median trough (mg/L)", f"{np.median(troughs):.3g}",
                    help=f"90% of patients between {np.percentile(troughs, 5):.3g} and "
                         f"{np.percentile(troughs, 95):.3g}", border=True)
        m[2].metric("Peaks, 90% range (mg/L)",
                    f"{np.percentile(peaks, 5):.3g} - {np.percentile(peaks, 95):.3g}", border=True)
        if has_window:
            m[3].metric("Target attainment", f"{population.target_attainment(peaks, troughs, window):.0%}",
                        border=True, help="Patients whose whole steady-state interval is inside the window.")
        h1, h2 = st.columns(2)
        h1.plotly_chart(charts.histogram(peaks, "Steady-state peaks", window[1] if has_window else None),
                        width="stretch", key="peak_hist")
        h2.plotly_chart(charts.histogram(troughs, "Steady-state troughs", window[0] if has_window else None),
                        width="stretch", key="trough_hist")

with finder_tab:
    st.markdown("Finds, for each dosing interval, the dose that puts steady state in the middle of the "
                "therapeutic window, using the current drug and patient. A short half-life with a long "
                "interval swings too far between peak and trough for any dose to fit.")
    if saturable:
        st.info("The finder relies on linearity (double the dose, double the concentration), which "
                "saturable elimination breaks. Adjust the dose by hand and watch the steady state.")
    elif not (0 < window[0] < window[1]):
        st.info("Set a therapeutic window with a lower limit above 0 in the sidebar to use the finder.")
    else:
        f1, f2 = st.columns(2)
        taus = f1.multiselect("Intervals to try (h)", [4, 6, 8, 12, 24, 48], [6, 8, 12, 24], key="finder_taus")
        step = f2.select_slider("Round doses to (mg)", [1, 5, 10, 25, 50, 100, 250], 25, key="finder_step")
        rows = dosing.suggest_regimens(scenario.regular, window, sorted(taus), step) if taus else []
        if rows:
            df = pd.DataFrame(rows).rename(columns={
                "tau": "Every (h)", "dose": "Dose (mg)", "daily_dose": "Daily dose (mg)",
                "peak": "Peak (mg/L)", "trough": "Trough (mg/L)", "fits": "Fits window"})
            st.dataframe(df.style.format({"Peak (mg/L)": "{:.3g}", "Trough (mg/L)": "{:.3g}",
                                          "Daily dose (mg)": "{:.0f}", "Dose (mg)": "{:g}"}),
                         hide_index=True, width="stretch")
            fitting = [r for r in rows if r["fits"]]
            if fitting:
                best = max(fitting, key=lambda r: r["tau"])  # fewest doses a day that still fits
                st.success(f"Longest interval that fits: **{best['dose']:g} mg every {best['tau']} h**.")
                st.button("Use this regimen", on_click=apply_regimen, args=(best["dose"], best["tau"]),
                          icon=":material/check:")
            else:
                st.warning("No interval fits: the swing between peak and trough is wider than the window. "
                           "Try shorter intervals or an IV infusion.")
    st.divider()
    st.subheader("Loading dose")
    ld = dosing.loading_dose(scenario) if not saturable else None
    if ld is None:
        st.caption("With saturable elimination the loading dose isn't a simple multiple; "
                   "use Css x V as a starting point.")
    else:
        st.markdown(f"Steady state takes several half-lives. A first dose of **{ld:.0f} mg** gives a first "
                    f"peak at the steady-state peak straight away, like filling the bath with a bucket "
                    f"before turning the tap down to a trickle.")
        st.button("Use this loading dose", on_click=lambda v=ld: st.session_state.update(loading=int(round(v, -1))),
                  icon=":material/bolt:")

with fit_tab:
    st.markdown("Measured a patient's blood levels? Enter them with the time since the **first** dose "
                "(the schedule in the sidebar is assumed to have been followed). The fit finds this "
                "patient's parameters, pulled towards the typical values as much as the population "
                "variability (Population tab) says is plausible: Bayesian (MAP) estimation, as in TDM software.")
    e1, e2 = st.columns([3, 2])
    with e1:
        levels_df = st.data_editor(
            st.session_state.levels_data, num_rows="dynamic", width="stretch",
            key=f"levels_{st.session_state.levels_version}",
            column_config={"time_h": st.column_config.NumberColumn("Time after first dose (h)", min_value=0.0),
                           "level_mg_L": st.column_config.NumberColumn("Level (mg/L)", min_value=0.0)})
        st.button("Fill in example levels", on_click=fill_example_levels, args=(scenario,))
        st.caption("The example: two levels from a made-up patient who clears the drug 40% faster than typical.")
    with e2:
        options = ["vmax", "km", "vd", "ka"] if saturable else ["cl", "vd", "ka"]
        if not oral:
            options = [o for o in options if o != "ka"]
        names = {"cl": "CL", "vd": "V", "ka": "ka", "vmax": "Vmax", "km": "Km"}
        to_fit = st.multiselect("Parameters to estimate", options, [o for o in options if o in ("cl", "vmax")] +
                                (["vd"] if "vd" in options else []), format_func=names.get, key="fit_params",
                                help="With one or two levels, estimate one or two parameters at most.")
        sd_add = st.number_input("Assay error, SD (mg/L)", 0.0, 10.0, 0.1, 0.05, key="sd_add")
        cv_prop = st.slider("Proportional error (CV %)", 1, 50, 15, key="cv_prop") / 100
    data = levels_df.dropna()
    data = data[(data.time_h > 0) & (data.level_mg_L > 0)]
    if data.empty:
        st.info("Add at least one level to fit.")
    else:
        prior_cv = {**cv, "km": cv.get("km", 0.3), "vmax": cv.get("vmax", 0.3), "ka": cv.get("ka", 0.4)}
        result = tdm.fit(scenario, data.time_h.values, data.level_mg_L.values, prior_cv, to_fit, sd_add, cv_prop)
        n_fit = max(n_doses, tdm.doses_before(scenario, data.time_h.max()) + 1)
        ft = scenario.grid(n_fit)
        with st.container(border=True):
            st.plotly_chart(charts.fitted_levels(ft, scenario.repeated(ft, n_fit), result.scenario.repeated(ft, n_fit),
                                                 data.time_h, data.level_mg_L, window, st.session_state.get("log")),
                            width="stretch", key="fit_chart")
        rows = [{"Parameter": names[p], "Typical / patient-adjusted": getattr(scenario, p) * (24 if p == "vmax" else 1),
                 "Fitted": getattr(result.scenario, p) * (24 if p == "vmax" else 1),
                 "Change": f"{result.factors[p] - 1:+.0%}"} for p in result.factors]
        st.dataframe(pd.DataFrame(rows).style.format({"Typical / patient-adjusted": "{:.3g}", "Fitted": "{:.3g}"}),
                     hide_index=True, width="stretch")
        st.caption("Vmax in mg/day. " + ("Fitted half-life: "
                   f"{result.scenario.half_life():.1f} h (typical {scenario.half_life():.1f} h)."
                   if np.isfinite(result.scenario.half_life()) else ""))
        st.button("Use the fitted parameters", on_click=apply_fit, args=(result.factors,), icon=":material/check:",
                  help="Moves the sidebar values, so the other tabs (and the Regimen finder) use this patient.")

with effect_tab:
    st.markdown("From concentration to effect, with the sigmoid Emax (Hill) model. A delay `ke0` "
                "lets the effect lag behind plasma, as when a drug acts in the brain: the effect then "
                "traces a loop against concentration (hysteresis) instead of a line.")
    d1, d2, d3, d4, d5 = st.columns(5)
    default_ec50 = float(round(np.sqrt(window[0] * window[1]), 2)) if window[0] > 0 and has_window else \
        float(round(max(average_steady_state(t, regular, scenario.tau), 0.01), 2))
    e0 = d1.number_input("Baseline, E0", -1000.0, 1000.0, 0.0, key="e0")
    emax_value = d2.number_input("Emax", -1000.0, 1000.0, 100.0, key="emax",
                                 help="Maximum change from baseline; negative for a drug that lowers something.")
    ec50 = d3.number_input("EC50 (mg/L)", 0.001, 10000.0, default_ec50, key="ec50")
    hill = d4.number_input("Hill coefficient", 0.1, 10.0, 1.0, 0.1, key="hill")
    delay = d5.number_input("Delay half-life (h)", 0.0, 48.0, 1.0, 0.25, key="ke0_half",
                            help="ln(2) / ke0. 0 = effect follows plasma instantly.")
    label = st.text_input("Effect name", "Effect (%)", key="effect_label")
    ke0 = np.inf if delay == 0 else np.log(2) / delay
    effect = pharmacodynamics.emax_effect(pharmacodynamics.effect_site(t, multi, ke0), e0, emax_value, ec50, hill)
    target = st.number_input("Target effect (optional, 0 = none)", -1000.0, 1000.0, 0.0, key="target_effect")
    st.plotly_chart(charts.effect_profile(t, effect, label, target or None), width="stretch", key="effect_chart")
    if reached:
        last = t >= t[-1] - scenario.tau
        k1, k2, k3 = st.columns(3)
        k1.metric("Highest effect at steady state", f"{effect[last].max():.3g}", border=True)
        k2.metric("Lowest effect at steady state", f"{effect[last].min():.3g}", border=True)
        if target:
            above = effect[last] >= target if emax_value >= 0 else effect[last] <= target
            k3.metric("Time at target", f"{above.mean():.0%}", border=True,
                      help="Share of a steady-state interval at or beyond the target effect.")
    t_half = scenario.half_life()
    span = min(t[-1], max(scenario.tau, 5 * t_half if np.isfinite(t_half) else scenario.tau))
    single_t = np.union1d(np.linspace(0, span, 400), [scenario.duration] if infusion else [])
    single_c = scenario.single(single_t)
    single_effect = pharmacodynamics.emax_effect(pharmacodynamics.effect_site(single_t, single_c, ke0),
                                                 e0, emax_value, ec50, hill)
    st.plotly_chart(charts.hysteresis(single_c, single_effect, single_t, label), width="stretch",
                    key="hysteresis_chart")

with about_tab:
    st.markdown(r"""
**One compartment.** The body is one well-mixed bathtub of volume $V_d$. Drug drains at a rate
proportional to how much is there: $k = CL / V_d$ per hour, so concentration falls exponentially and
$t_{1/2} = \ln 2 / k$.

**Two compartments.** A second bathtub (tissue, $V_2$) joined to the first (plasma, $V_1$) by a
pipe of size $Q$. After a dose the level falls fast while drug flows into tissue (distribution),
then slowly once the tubs balance (terminal phase). Same AUC, different shape.

**Routes.** An IV bolus pours the dose in at once; an infusion runs the tap at a steady rate for a
set time, which trims the peak; an oral dose drips in from the gut at rate $k_a$, and only a
fraction $F$ arrives.

**Saturable elimination (Michaelis-Menten).** Elimination rate $= V_{max} C / (K_m + C)$. Below
$K_m$ it behaves like ordinary clearance; above it the enzymes are busy and the rate flattens out.
Small dose increases can then cause big concentration jumps, and steady state takes much longer.
Phenytoin is the classic example.

**Patient adjustment.** Volumes scale with weight, clearances with weight$^{0.75}$ (allometry).
The renal part of clearance (fraction $f_e$) follows creatinine clearance from Cockcroft-Gault.

**Fitting to levels (TDM).** Bayesian MAP estimation: the parameters that best explain the
measured levels, penalised for straying from the typical values by the population variability.
One level mostly moves one parameter; several levels across an interval pin down CL and V.

**Effect (PD).** Sigmoid Emax: $E = E_0 + E_{max} C_e^n / (EC_{50}^n + C_e^n)$, with an effect
compartment $dC_e/dt = k_{e0} (C - C_e)$ for effects that lag behind plasma.

**Population variability.** Each parameter is the typical value times $e^{\eta}$, with $\eta$
normally distributed: a log-normal spread, as in population PK software.

Teaching tool: population averages, not for clinical decisions.
""")
