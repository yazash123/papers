"""Competing-exits model: a waiting state left by the first of several exits.

A *waiting state* is left by the first of several independent Poisson exits
(competing exponential clocks).  Exit i has a load-dependent rate k_i(F), so
the probability that exit i wins is P_i = k_i / sum_j k_j, and the time to
leave is exponential with rate K = sum_j k_j, independent of which exit wins.

Exits come in four kinds:

* PRODUCTIVE  - the outcome the system is after (forward step, correct
                product, capture);
* WRONG       - the competing outcome (backstep, wrong product);
* TERMINATING - ends the run (detachment);
* RESTART     - the attempt is abandoned and the system goes back to the
                start, paying a time (``cost_time``) and/or fuel cost (e.g.
                ATP leaves before the head commits; a proofreading discard).

An *attempt* starts with an optional entry wait (mean ``entry_time``, e.g. the
wait for ATP to bind), then the race.  After a PRODUCTIVE or WRONG exit the
system spends ``cost_time`` (e.g. the rest of the mechanochemical cycle) and
then starts a new attempt.

Several waiting states can be chained into a small network (``Motor``): an exit
may lead to another waiting state (``to``) instead of back to the start.  This
is used for kinesin's fast events after a backstep, and later for multi-step
proofreading.  The single-state formulas are on ``WaitingState``; the network
versions (first-step analysis) are on ``Motor``.

Units: rates 1/s, times s, loads pN (F > 0 hindering), distances nm,
kT pN*nm.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace
from enum import Enum
from typing import Mapping, Optional, Union

import numpy as np
from scipy.optimize import brentq

from .rates import Scaled

Destination = Union[None, str, Mapping[Optional[str], float]]


class ExitKind(str, Enum):
    PRODUCTIVE = "productive"
    WRONG = "wrong"
    TERMINATING = "terminating"
    RESTART = "restart"


STEP_KINDS = (ExitKind.PRODUCTIVE, ExitKind.WRONG)


@dataclass(frozen=True)
class Exit:
    """One way out of a waiting state.

    name       label used to refer to the exit
    kind       an ExitKind
    rate       rate law, called as rate(F, kT) -> 1/s (see ``rates``)
    step       displacement (nm) when the exit fires (+d forward, -d back)
    cost_time  mean time (s) spent after the exit fires, before the system
               reaches its destination (post-step rest of cycle, or the
               delay of a restart)
    fuel       fuel units consumed when the exit fires (bookkeeping only)
    to         where the system goes next: None = back to the start of a new
               attempt (entry wait first); a state name = straight into that
               waiting state; or a mapping {destination: probability}.
               Ignored for TERMINATING exits.
    """

    name: str
    kind: ExitKind
    rate: object
    step: float = 0.0
    cost_time: float = 0.0
    fuel: float = 0.0
    to: Destination = None

    def k(self, F, kT):
        return self.rate(F, kT)

    def destinations(self) -> dict:
        """Destination probabilities as a dict {name or None: p}."""
        if self.to is None or isinstance(self.to, str):
            return {self.to: 1.0}
        dest = dict(self.to)
        total = sum(dest.values())
        if not np.isclose(total, 1.0):
            raise ValueError(f"destination probabilities of {self.name!r} sum to {total}")
        return dest


@dataclass(frozen=True)
class WaitingState:
    """A waiting state with competing exits (all formulas at constant load F)."""

    name: str
    exits: tuple
    kT: float
    entry_time: float = 0.0

    def __post_init__(self):
        names = [e.name for e in self.exits]
        if len(set(names)) != len(names):
            raise ValueError("exit names must be unique")
        object.__setattr__(self, "exits", tuple(self.exits))

    # ----- access -----------------------------------------------------------
    def exit(self, name: str) -> Exit:
        for e in self.exits:
            if e.name == name:
                return e
        raise KeyError(name)

    def of_kind(self, kind: ExitKind) -> list:
        return [e for e in self.exits if e.kind == kind]

    def scaled(self, name: str, factor: float) -> "WaitingState":
        """Copy with one exit's rate multiplied by ``factor`` (e.g. a gate)."""
        new = tuple(replace(e, rate=Scaled(e.rate, factor)) if e.name == name else e
                    for e in self.exits)
        self.exit(name)  # raises if missing
        return replace(self, exits=new)

    # ----- splitting --------------------------------------------------------
    def rate(self, name: str, F=0.0):
        return self.exit(name).k(F, self.kT)

    def total_rate(self, F=0.0):
        return sum(e.k(F, self.kT) for e in self.exits)

    def splitting(self, F=0.0) -> dict:
        """P_i = k_i / sum_j k_j for every exit."""
        K = self.total_rate(F)
        return {e.name: e.k(F, self.kT) / K for e in self.exits}

    def prob(self, name: str, F=0.0):
        return self.rate(name, F) / self.total_rate(F)

    def prob_kind(self, kind: ExitKind, F=0.0):
        K = self.total_rate(F)
        return sum(e.k(F, self.kT) for e in self.of_kind(kind)) / K

    # ----- the bet ----------------------------------------------------------
    def log_odds(self, a: str, b: str, F=0.0):
        """ln(k_a / k_b), in nats."""
        return np.log(self.rate(a, F) / self.rate(b, F))

    def lever(self, a: str, b: str, F=0.0, dF=1e-4):
        """Load distance over which F erodes the odds: -kT d ln(k_a/k_b)/dF (nm).

        For two Bell rates this is delta_a - delta_b at every load.
        """
        return -self.kT * (self.log_odds(a, b, F + dF) - self.log_odds(a, b, F - dF)) / (2 * dF)

    def balance_point(self, a: str, b: str, bracket=(-100.0, 200.0)):
        """Load at which the odds between exits a and b are 1:1.

        For Bell rates this is F = kT ln(k_a0/k_b0) / (delta_a - delta_b).
        """
        f = lambda F: float(self.log_odds(a, b, F))
        return brentq(f, *bracket, xtol=1e-12)

    def price(self, name: str, F=0.0):
        """Price (nats) of conditioning one attempt on exit ``name``: ln(1/P)."""
        return -np.log(self.prob(name, F))

    def attempts_per(self, name: str, F=0.0):
        """Mean number of races run per firing of exit ``name``: 1/P."""
        return 1.0 / self.prob(name, F)

    # ----- time and motion (single state; every exit returns to the start) --
    def _check_single(self):
        for e in self.exits:
            if e.kind != ExitKind.TERMINATING and e.to is not None:
                raise ValueError("use Motor for exits that lead to other states")

    def attempt_time(self, F=0.0):
        """Mean duration of one attempt, excluding post-step costs:
        entry wait + race + (restart delay if the attempt restarts)."""
        P = self.splitting(F)
        restart = sum(P[e.name] * e.cost_time for e in self.of_kind(ExitKind.RESTART))
        return self.entry_time + 1.0 / self.total_rate(F) + restart

    def mean_dwell(self, F=0.0):
        """Mean time from one step (forward or back) to the next, conditional on
        the next event being a step:

            attempt_time / (1 - P_restart) + mean post-step cost.

        (The race time does not depend on which exit wins, so conditioning on a
        step rather than a termination leaves the time unchanged.)"""
        self._check_single()
        P = self.splitting(F)
        P_r = sum(P[e.name] for e in self.of_kind(ExitKind.RESTART))
        steps = [e for e in self.exits if e.kind in STEP_KINDS]
        P_s = sum(P[e.name] for e in steps)
        post = sum(P[e.name] * e.cost_time for e in steps) / P_s
        return self.attempt_time(F) / (1.0 - P_r) + post

    def velocity(self, F=0.0):
        """Mean drift per unit attached time (nm/s): displacement per attempt
        over time per attempt (renewal-reward).  Without terminating exits this
        is (mean step displacement) / mean_dwell."""
        self._check_single()
        P = self.splitting(F)
        steps = [e for e in self.exits if e.kind in STEP_KINDS]
        disp = sum(P[e.name] * e.step for e in steps)
        time = self.attempt_time(F) + sum(P[e.name] * e.cost_time for e in steps)
        return disp / time

    def run_stats(self, F=0.0) -> dict:
        """Mean steps, displacement (nm) and attached time (s) per run at a
        constant load, for a state with at least one terminating exit."""
        self._check_single()
        P = self.splitting(F)
        P_t = sum(P[e.name] for e in self.of_kind(ExitKind.TERMINATING))
        if P_t <= 0:
            raise ValueError("no terminating exit: runs never end")
        steps = [e for e in self.exits if e.kind in STEP_KINDS]
        n_steps = sum(P[e.name] for e in steps) / P_t
        disp = sum(P[e.name] * e.step for e in steps) / P_t
        time = (self.attempt_time(F) + sum(P[e.name] * e.cost_time for e in steps)) / P_t
        return {"steps": n_steps, "displacement": disp, "time": time}

    def fuel_per_attempt(self, F=0.0):
        P = self.splitting(F)
        return sum(P[e.name] * e.fuel for e in self.exits)


@dataclass(frozen=True)
class Motor:
    """A network of waiting states.  Attempts start at ``home``.

    Constant-load results use first-step analysis over one *cycle*: from the
    start of an attempt at home until the system is back at the start of a new
    attempt, or the run terminates.
    """

    states: tuple
    home: str
    step_size: float = 8.2

    def __post_init__(self):
        object.__setattr__(self, "states", tuple(self.states))
        names = [s.name for s in self.states]
        if len(set(names)) != len(names):
            raise ValueError("state names must be unique")
        if self.home not in names:
            raise ValueError(f"home state {self.home!r} not in network")
        kTs = {s.kT for s in self.states}
        if len(kTs) != 1:
            raise ValueError("all states must share one kT")
        for s in self.states:
            for e in s.exits:
                for d in e.destinations():
                    if d is not None and d not in names:
                        raise ValueError(f"exit {e.name!r} leads to unknown state {d!r}")

    @classmethod
    def single(cls, state: WaitingState, step_size: float = 8.2) -> "Motor":
        return cls((state,), state.name, step_size)

    @property
    def kT(self):
        return self.states[0].kT

    def state(self, name: str) -> WaitingState:
        for s in self.states:
            if s.name == name:
                return s
        raise KeyError(name)

    @property
    def home_state(self) -> WaitingState:
        return self.state(self.home)

    def replace_state(self, new: WaitingState) -> "Motor":
        return replace(self, states=tuple(new if s.name == new.name else s for s in self.states))

    def scaled(self, state: str, exit_name: str, factor: float) -> "Motor":
        return self.replace_state(self.state(state).scaled(exit_name, factor))

    # ----- first-step analysis at constant load ------------------------------
    def cycle_stats(self, F=0.0) -> dict:
        """Expected time, displacement, steps, fuel and termination probability
        over one cycle started at home (entry wait included)."""
        names = [s.name for s in self.states]
        idx = {n: i for i, n in enumerate(names)}
        n = len(names)
        A = np.eye(n)
        b = np.zeros((n, 5))  # columns: time, displacement, steps, fuel, P(term)
        for s in self.states:
            i = idx[s.name]
            P = s.splitting(F)
            b[i, 0] += 1.0 / s.total_rate(F)
            for e in s.exits:
                p = P[e.name]
                if e.kind != ExitKind.TERMINATING:  # time after termination is not attached time
                    b[i, 0] += p * e.cost_time
                b[i, 1] += p * e.step
                b[i, 2] += p * (e.kind in STEP_KINDS)
                b[i, 3] += p * e.fuel
                if e.kind == ExitKind.TERMINATING:
                    b[i, 4] += p
                    continue
                for dest, q in e.destinations().items():
                    if dest is not None:  # None: cycle ends
                        A[i, idx[dest]] -= p * q
        x = np.linalg.solve(A, b)[idx[self.home]]
        return {"time": self.home_state.entry_time + x[0], "displacement": x[1],
                "steps": x[2], "fuel": x[3], "p_term": x[4]}

    def velocity(self, F=0.0):
        """Mean drift per unit attached time (nm/s)."""
        c = self.cycle_stats(F)
        return c["displacement"] / c["time"]

    def run_stats(self, F=0.0) -> dict:
        """Mean steps, displacement and attached time per run at constant load."""
        c = self.cycle_stats(F)
        if c["p_term"] <= 0:
            raise ValueError("no terminating exit reachable: runs never end")
        return {"steps": c["steps"] / c["p_term"],
                "displacement": c["displacement"] / c["p_term"],
                "time": c["time"] / c["p_term"]}

    def stall_force(self, bracket=(-20.0, 60.0)):
        """Load at which the drift velocity is zero."""
        return brentq(lambda F: self.velocity(F), *bracket, xtol=1e-10)
