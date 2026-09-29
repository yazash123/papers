"""Gillespie simulation of a competing-exits motor under a loading protocol.

The load depends only on position and position changes only when a step
fires, so within one race all rates are constant and the direct Gillespie
method is exact: draw the race time from Exp(K) and the winning exit with
probability k_i/K.  Entry waits and post-exit costs are exponential by default
(``fixed_delays=True`` makes them deterministic; means are unchanged).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .model import ExitKind, Motor


@dataclass
class Run:
    time: float            # attached time (s) until termination or t_max
    x: float               # final position (nm)
    load: float            # load at the end (pN): the load at detachment if terminated
    max_load: float        # largest load reached during the run
    terminated: bool       # False if stopped by t_max
    counts: dict           # number of firings of each exit (state.exit -> n)
    fuel: float


class _RateCache:
    """Cumulative exit rates per (state, load); loads repeat on a lattice."""

    def __init__(self, motor: Motor):
        self.motor = motor
        self.cache = {}

    def get(self, state_name: str, F: float):
        key = (state_name, round(float(F), 9))
        hit = self.cache.get(key)
        if hit is None:
            s = self.motor.state(state_name)
            k = np.array([float(e.k(F, s.kT)) for e in s.exits])
            hit = (np.cumsum(k), float(k.sum()), s.exits)
            self.cache[key] = hit
        return hit


def _delay(rng, mean, fixed):
    if mean <= 0:
        return 0.0
    return mean if fixed else rng.exponential(mean)


def _pick_destination(rng, exit_):
    dest = exit_.destinations()
    if len(dest) == 1:
        return next(iter(dest))
    keys = list(dest)
    p = np.array([dest[k] for k in keys])
    return keys[rng.choice(len(keys), p=p)]


def simulate_run(motor: Motor, load, rng, t_max=np.inf, x0=0.0,
                 fixed_delays=False, cache=None, trace=False):
    """One run from position x0 until termination or t_max.

    ``load`` maps position (nm) to hindering load (pN), e.g. ``Trap(0.05)``.
    Returns a ``Run`` (and, if ``trace``, a list of (t, x) after each step).
    """
    cache = cache or _RateCache(motor)
    t, x, fuel = 0.0, float(x0), 0.0
    max_load = float(load(x))
    counts = {}
    path = [(0.0, x)] if trace else None
    state = None  # None: at the start of a new attempt
    while True:
        if state is None:
            t += _delay(rng, motor.home_state.entry_time, fixed_delays)
            state = motor.home
        if t >= t_max:
            break
        F = float(load(x))
        cum, K, exits = cache.get(state, F)
        t += rng.exponential(1.0 / K)
        if t >= t_max:
            break
        i = int(np.searchsorted(cum, rng.random() * K, side="right"))
        e = exits[min(i, len(exits) - 1)]
        key = f"{state}.{e.name}"
        counts[key] = counts.get(key, 0) + 1
        fuel += e.fuel
        if e.kind == ExitKind.TERMINATING:
            run = Run(t, x, F, max_load, True, counts, fuel)
            return (run, path) if trace else run
        x += e.step
        max_load = max(max_load, float(load(x)))
        if trace and e.step != 0:
            path.append((t, x))
        t += _delay(rng, e.cost_time, fixed_delays)
        state = _pick_destination(rng, e)
    run = Run(min(t, t_max), x, float(load(x)), max_load, False, counts, fuel)
    return (run, path) if trace else run


def simulate_runs(motor: Motor, load, n_runs, seed=0, t_max=np.inf, x0=0.0,
                  fixed_delays=False) -> dict:
    """Many independent runs; returns arrays of per-run results."""
    rng = np.random.default_rng(seed)
    cache = _RateCache(motor)
    runs = [simulate_run(motor, load, rng, t_max, x0, fixed_delays, cache) for _ in range(n_runs)]
    keys = sorted({k for r in runs for k in r.counts})
    return {
        "time": np.array([r.time for r in runs]),
        "x": np.array([r.x for r in runs]),
        "load": np.array([r.load for r in runs]),
        "max_load": np.array([r.max_load for r in runs]),
        "terminated": np.array([r.terminated for r in runs]),
        "fuel": np.array([r.fuel for r in runs]),
        "counts": {k: np.array([r.counts.get(k, 0) for r in runs]) for k in keys},
    }
