"""PK metrics computed from a simulated concentration-time curve.

These work on the curve itself (arrays of times t and concentrations c), the way
non-compartmental analysis works on measured samples. They don't need to know
which model produced the curve.
"""
import numpy as np


def cmax(c):
    """Peak concentration (mg/L)."""
    return np.max(c)


def tmax(t, c):
    """Time of the peak (h). With repeated doses: the first time the overall maximum is hit."""
    return t[np.argmax(c)]


def auc(t, c):
    """Area under the concentration-time curve by the trapezoidal rule (mg*h/L).

    Between neighbouring samples the curve is treated as a straight line, so each
    slice is a trapezoid: width (t[i+1] - t[i]) times mean height (c[i] + c[i+1]) / 2.
    """
    t, c = np.asarray(t), np.asarray(c)
    return np.sum(np.diff(t) * (c[1:] + c[:-1]) / 2)


def half_life(t, c, tail=0.3):
    """Terminal half-life (h), estimated from a single-dose curve.

    In the terminal phase C = C0 * exp(-lambda_z * t), so ln(C) plotted against t
    is a straight line with slope -lambda_z. We fit that line to the last `tail`
    fraction of the time span, then t1/2 = ln(2) / lambda_z.
    """
    t, c = np.asarray(t), np.asarray(c)
    late = (t >= t[-1] * (1 - tail)) & (c > 0)  # c > 0 because ln(0) is undefined
    slope, _intercept = np.polyfit(t[late], np.log(c[late]), 1)
    return np.log(2) / -slope


def _intervals(t, tau):
    """Indices of the time points in each complete dosing interval [i*tau, (i+1)*tau].

    Both ends are included so the last interval keeps its final point: for an IV
    bolus that is the true trough, the instant before a next dose would go in.
    """
    return [np.flatnonzero((t >= i * tau) & (t <= (i + 1) * tau)) for i in range(int(t[-1] // tau))]


def steady_state_peak_trough(t, c, tau):
    """Highest and lowest concentration in the last dosing interval: the
    steady-state peak and trough, if the simulation ran long enough to get there."""
    last = _intervals(t, tau)[-1]
    return c[last].max(), c[last].min()


def time_to_steady_state(t, c, tau, fraction=0.95):
    """Time (h) at which the trough first reaches `fraction` of its steady-state value.

    The trough is the lowest concentration in each dosing interval. The last
    interval's trough is taken as the steady-state value, so the simulation must
    run long enough to get there (7 half-lives puts it within 1%).

    Each half-life closes half of the remaining gap to steady state, so after n
    half-lives you are at 1 - 0.5**n of it: 94% after 4, 97% after 5. 95% takes
    log2(1 / 0.05) = 4.3 half-lives. The dose never appears in that, only k does.
    """
    lowest = [i[np.argmin(c[i])] for i in _intervals(t, tau)]  # index of each interval's trough
    troughs = c[lowest]
    first = np.argmax(troughs >= fraction * troughs[-1])  # argmax of booleans = first True
    return t[lowest[first]]
