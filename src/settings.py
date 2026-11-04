"""Saving, loading and sharing the app's inputs.

A scenario can leave the app two ways: as a JSON file, or squeezed into the URL
(?preset=Phenytoin&dose=130...) so a link reproduces exactly what you see. Both
come back in through `parse`, which only accepts known inputs of the right type
and inside their ranges, so a hand-edited file or link can't break the app.
"""
import json

from src.presets import PRESETS
from src.simulate import ROUTES

# input key -> (type, allowed values or (min, max))
SCHEMA = {
    "preset": (str, tuple(PRESETS)),
    "route": (str, ROUTES),
    "compartments": (str, ("One", "Two")),
    "elimination": (str, ("Linear", "Saturable")),
    "dose": (int, (10, 2000)),
    "tau": (int, (1, 48)),
    "cl": (float, (0.1, 50.0)),
    "vd": (float, (1.0, 500.0)),
    "ka": (float, (0.05, 5.0)),
    "f": (float, (0.05, 1.0)),
    "duration": (float, (0.1, 24.0)),
    "q": (float, (0.1, 100.0)),
    "v2": (float, (1.0, 1000.0)),
    "vmax": (float, (10.0, 3000.0)),   # mg/day, as shown in the app
    "km": (float, (0.1, 100.0)),
    "loading": (int, (0, 20000)),
    "missed": (int, (0, 200)),
    "double_next": (bool, None),
    "late_dose": (int, (0, 200)),
    "late_hours": (float, (0.5, 24.0)),
    "adjust": (bool, None),
    "weight": (float, (3.0, 250.0)),
    "age": (int, (18, 100)),
    "sex": (str, ("Male", "Female")),
    "scr": (float, (0.2, 10.0)),
    "fe": (float, (0.0, 1.0)),
    "window_lo": (float, (0.0, 1000.0)),
    "window_hi": (float, (0.0, 1000.0)),
}
RANGES = {key: limits for key, (kind, limits) in SCHEMA.items() if kind in (int, float)}


def clamp(key, value):
    lo, hi = RANGES[key]
    kind = SCHEMA[key][0]
    return kind(min(max(value, lo), hi))


def _convert(key, raw):
    kind, allowed = SCHEMA[key]
    if kind is bool:
        if isinstance(raw, str):
            if raw.lower() not in ("true", "false", "1", "0"):
                raise ValueError
            return raw.lower() in ("true", "1")
        return bool(raw)
    if kind is str:
        if raw not in allowed:
            raise ValueError
        return raw
    value = kind(float(raw))
    lo, hi = allowed
    if not lo <= value <= hi:
        raise ValueError
    return value


def parse(data):
    """Known, valid inputs from a dict (JSON or query params); returns (values, rejected keys)."""
    values, rejected = {}, []
    for key, raw in data.items():
        if key not in SCHEMA:
            rejected.append(key)
            continue
        try:
            values[key] = _convert(key, raw)
        except (TypeError, ValueError):
            rejected.append(key)
    return values, rejected


def snapshot(state):
    """The saveable inputs currently in session_state."""
    return {key: state[key] for key in SCHEMA if key in state}


def to_json(state):
    return json.dumps(snapshot(state), indent=2)


def from_json(text):
    data = json.loads(text)
    if not isinstance(data, dict):
        raise ValueError("expected a JSON object")
    return parse(data)


def to_query(state):
    """Strings for st.query_params: a shareable link."""
    return {key: str(value) for key, value in snapshot(state).items()}
