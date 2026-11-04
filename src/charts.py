"""Interactive (Plotly) charts for the Streamlit app. Hover for values, drag to
zoom, double-click to reset. src/plotting.py keeps the static matplotlib figure
for the README."""
import numpy as np
import plotly.graph_objects as go

from src.plotting import AQUA, BLUE, GREY, GRID, ORANGE, SURFACE

COMPARISON_COLOURS = ["#8a5cd6", "#c2417a", "#3b8c3b", "#b8860b"]
FONT = "Inter, sans-serif"


def _base(log, y_max, title=None):
    fig = go.Figure()
    fig.update_layout(
        template="plotly_white", paper_bgcolor=SURFACE, plot_bgcolor=SURFACE, height=440,
        font=dict(family=FONT, color="#0b0b0b"), margin=dict(l=10, r=10, t=40 if title else 10, b=10),
        title=dict(text=title, x=0, font=dict(size=15)) if title else None,
        legend=dict(orientation="h", yanchor="top", y=-0.15, xanchor="center", x=0.5),
        hovermode="x unified",
    )
    fig.update_xaxes(title="Time (h)", showgrid=False, linecolor=GREY, ticks="outside", zeroline=False)
    fig.update_yaxes(title="Plasma concentration (mg/L)", gridcolor=GRID, linecolor=GREY, zeroline=False)
    if log:
        # log axes can't show zero: start three decades below the top
        fig.update_yaxes(type="log", range=[np.log10(y_max) - 3, np.log10(y_max * 1.3)])
    else:
        fig.update_yaxes(range=[0, y_max * 1.1])
    return fig


def _window(fig, window):
    low, high = window
    if high > low:  # (0, 0) means no window to draw
        fig.add_hrect(y0=max(low, 1e-12), y1=high, fillcolor=AQUA, opacity=0.15, line_width=0,
                      annotation_text="Therapeutic window", annotation_position="top left",
                      annotation_font=dict(color=GREY, size=11))


def _line(fig, t, c, name, colour, dash=None, width=2):
    fig.add_trace(go.Scatter(x=t, y=c, name=name, mode="lines",
                             line=dict(color=colour, width=width, dash=dash),
                             hovertemplate="%{y:.3g} mg/L"))


def profile(t, single, multi, window, steady_state=None, comparisons=(), dose_times=(), log=False,
            show_single=True):
    """Single dose vs repeated dosing, with the therapeutic window, steady-state
    peak and trough, pinned comparison curves and dose markers."""
    tops = [multi.max(), window[1]] + [c.max() for _, _, c in comparisons]
    fig = _base(log, max(tops) or 1.0)
    _window(fig, window)
    for i, (label, tc, c) in enumerate(comparisons):
        _line(fig, tc, c, label, COMPARISON_COLOURS[i % len(COMPARISON_COLOURS)], width=1.5)
    if show_single:
        _line(fig, t, single, "Single dose", ORANGE, dash="dash")
    _line(fig, t, multi, "Repeated dosing", BLUE)
    if steady_state:
        for i, value in enumerate(steady_state):
            fig.add_trace(go.Scatter(x=[t[0], t[-1]], y=[value, value], mode="lines", name="Steady-state peak / trough",
                                     line=dict(color=GREY, width=1, dash="dot"), legendgroup="ss",
                                     showlegend=i == 0, hoverinfo="skip"))
    if len(dose_times):
        # markers on a hidden second axis so they sit at the bottom on linear and log scales
        fig.add_trace(go.Scatter(x=list(dose_times), y=[0.02] * len(dose_times), yaxis="y2", name="Doses",
                                 mode="markers", marker=dict(symbol="triangle-up", size=9, color=GREY),
                                 hovertemplate="dose at %{x:.1f} h"))
        fig.update_layout(yaxis2=dict(overlaying="y", range=[0, 1], visible=False, fixedrange=True))
    fig.update_xaxes(range=[0, t[-1]])
    return fig


def population(t, bands, typical, window, log=False):
    """Prediction bands across simulated patients: 90% (5th-95th percentile) and
    50% (25th-75th) of people fall inside the outer and inner bands."""
    fig = _base(log, max(bands[95].max(), window[1]) or 1.0)
    _window(fig, window)
    for low, high, opacity, name in ((5, 95, 0.18, "90% of patients"), (25, 75, 0.32, "50% of patients")):
        fig.add_trace(go.Scatter(x=np.concatenate([t, t[::-1]]), y=np.concatenate([bands[high], bands[low][::-1]]),
                                 fill="toself", fillcolor=f"rgba(42,120,214,{opacity})", line=dict(width=0),
                                 name=name, hoverinfo="skip"))
    _line(fig, t, bands[50], "Median patient", BLUE)
    _line(fig, t, typical, "Typical parameters", GREY, dash="dot", width=1.5)
    fig.update_xaxes(range=[0, t[-1]])
    return fig


def histogram(values, name, window_edge=None):
    fig = go.Figure(go.Histogram(x=values, nbinsx=30, marker_color=BLUE, opacity=0.8, name=name))
    fig.update_layout(template="plotly_white", paper_bgcolor=SURFACE, plot_bgcolor=SURFACE, height=260,
                      font=dict(family=FONT), margin=dict(l=10, r=10, t=30, b=10), showlegend=False,
                      title=dict(text=name, x=0, font=dict(size=13)), bargap=0.05)
    fig.update_xaxes(title="mg/L")
    fig.update_yaxes(title="Patients", gridcolor=GRID)
    if window_edge:
        fig.add_vline(x=window_edge, line=dict(color=AQUA, width=2, dash="dash"))
    return fig


def fitted_levels(t, prior, fitted, times, levels, window, log=False):
    """Measured levels against the population (prior) and fitted (this patient) curves."""
    fig = _base(log, max(prior.max(), fitted.max(), max(levels, default=0), window[1]) or 1.0)
    _window(fig, window)
    _line(fig, t, prior, "Population prediction", GREY, dash="dot", width=1.5)
    _line(fig, t, fitted, "Fitted to this patient", BLUE)
    fig.add_trace(go.Scatter(x=list(times), y=list(levels), mode="markers", name="Measured levels",
                             marker=dict(color=ORANGE, size=11, symbol="circle", line=dict(color="white", width=1.5)),
                             hovertemplate="%{y:.3g} mg/L at %{x:.1f} h"))
    fig.update_xaxes(range=[0, t[-1]])
    return fig


def effect_profile(t, effect, label, target=None):
    """Effect over time, for the repeated-dosing curve."""
    fig = go.Figure()
    fig.update_layout(template="plotly_white", paper_bgcolor=SURFACE, plot_bgcolor=SURFACE, height=380,
                      font=dict(family=FONT), margin=dict(l=10, r=10, t=10, b=10), hovermode="x unified",
                      legend=dict(orientation="h", yanchor="top", y=-0.18, xanchor="center", x=0.5))
    fig.update_xaxes(title="Time (h)", showgrid=False, linecolor=GREY, range=[0, t[-1]])
    fig.update_yaxes(title=label, gridcolor=GRID, linecolor=GREY)
    if target is not None:
        fig.add_hline(y=target, line=dict(color=AQUA, width=2, dash="dash"), annotation_text="Target effect",
                      annotation_position="top left", annotation_font=dict(color=GREY, size=11))
    _line(fig, t, effect, "Effect, repeated dosing", "#8a5cd6")
    return fig


def hysteresis(c, effect, t, label):
    """Effect against plasma concentration after one dose, coloured by time. With an
    effect-site delay the rising and falling limbs differ, tracing a loop."""
    fig = go.Figure(go.Scatter(
        x=c, y=effect, mode="lines+markers", name="Single dose",
        marker=dict(size=4, color=t, colorscale="Viridis", showscale=True, colorbar=dict(title="h")),
        line=dict(color=GRID, width=1), hovertemplate="C %{x:.3g} mg/L, effect %{y:.3g}<extra></extra>"))
    fig.update_layout(template="plotly_white", paper_bgcolor=SURFACE, plot_bgcolor=SURFACE, height=380,
                      font=dict(family=FONT), margin=dict(l=10, r=10, t=30, b=10), showlegend=False,
                      title=dict(text="Effect vs plasma concentration, one dose", x=0, font=dict(size=13)))
    fig.update_xaxes(title="Plasma concentration (mg/L)", gridcolor=GRID)
    fig.update_yaxes(title=label, gridcolor=GRID)
    return fig
