"""Concentration-time plot shared by the Streamlit app and the README figure."""
from matplotlib.figure import Figure

BLUE, ORANGE, AQUA, GREY, GRID = "#2a78d6", "#eb6834", "#1baf7a", "#52514e", "#e6e6e6"
SURFACE = "#fcfcfb"  # matches backgroundColor in .streamlit/config.toml


def plot_profile(t, single, multi, window, steady_state=None):
    """Single dose vs repeated dosing, with the therapeutic window shaded.

    window:       (lower, upper) therapeutic range in mg/L
    steady_state: optional (peak, trough) at steady state, drawn as dotted lines
    """
    # Figure() rather than plt.subplots(): no global pyplot state, so a web app
    # that redraws on every slider move doesn't pile up open figures.
    fig = Figure(figsize=(9, 4.5), layout="constrained", facecolor=SURFACE)
    ax = fig.subplots()
    ax.set_facecolor(SURFACE)

    if window[1] > window[0]:  # (0, 0) means no window to draw
        ax.axhspan(*window, color=AQUA, alpha=0.15, lw=0, label="Therapeutic window")
    ax.plot(t, single, color=ORANGE, lw=2, ls="--", label="Single dose")
    ax.plot(t, multi, color=BLUE, lw=2, label="Repeated dosing")
    if steady_state:
        peak, trough = steady_state
        ax.axhline(peak, color=GREY, lw=1, ls=":", label="Steady-state peak / trough")
        ax.axhline(trough, color=GREY, lw=1, ls=":")

    ax.set_xlim(0, t[-1])
    ax.set_ylim(0, max(multi.max(), window[1]) * 1.1)
    ax.set_xlabel("Time (h)")
    ax.set_ylabel("Plasma concentration (mg/L)")
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(GREY)
    ax.tick_params(colors=GREY)
    ax.grid(axis="y", color=GRID, lw=0.8)
    fig.legend(loc="outside lower center", ncol=4, frameon=False)
    return fig
