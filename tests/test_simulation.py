"""Gillespie simulation against the analytic results.

Tolerances are statistical: each check asserts |sim - exact| < 4 standard
errors, with the standard error estimated from the simulated runs.
"""
import numpy as np
import pytest

from competing_exits import Clamp, Motor, Trap, simulate_run, simulate_runs, trap_statistics
from competing_exits.kinesin import head_race_v2, kondo_kif5a
from competing_exits.simulate import _RateCache


def within(sim_values, exact, z=4.0):
    sim_values = np.asarray(sim_values, dtype=float)
    se = sim_values.std(ddof=1) / np.sqrt(len(sim_values))
    assert abs(sim_values.mean() - exact) < z * se, (sim_values.mean(), exact, se)


@pytest.mark.parametrize("F,atp", [(0.0, 1000.0), (5.0, 1000.0), (3.0, 10.0)])
def test_v2_constant_load_velocity_and_splitting(F, atp):
    """No terminating exit: simulate many fixed-length stretches of time."""
    state = head_race_v2(atp)
    m = Motor.single(state)
    rng = np.random.default_rng(1)
    cache = _RateCache(m)
    t_max = 2.0
    v, fwd_frac = [], []
    for _ in range(300):
        run = simulate_run(m, Clamp(F), rng, t_max=t_max, cache=cache)
        v.append(run.x / t_max)
        nf = run.counts.get("atp_bound.forward", 0)
        nb = run.counts.get("atp_bound.back", 0)
        if nf + nb:
            fwd_frac.append(nf / (nf + nb))
    # each 2-s stretch starts afresh, so the mean is biased by less than one
    # cycle per run (dwell/t_max < 3%); allow for that on top of the noise
    exact = state.velocity(F)
    se = np.std(v, ddof=1) / np.sqrt(len(v))
    bias = 8.2 / t_max
    assert abs(np.mean(v) - exact) < 4 * se + bias
    P = state.splitting(F)
    within(fwd_frac, P["forward"] / (P["forward"] + P["back"]), z=5)


def test_v2_dwell_distribution_mean():
    """Mean time between steps equals mean_dwell (with ATP wait, restarts, T)."""
    state = head_race_v2(10.0)
    m = Motor.single(state)
    rng = np.random.default_rng(2)
    run, path = simulate_run(m, Clamp(4.0), rng, t_max=200.0, trace=True)
    t = np.array([p[0] for p in path[1:]])
    # the recorded time is when the step fires; the post-step cost comes after,
    # so successive differences are exactly the step-to-step dwells
    dwells = np.diff(t)
    within(dwells, state.mean_dwell(4.0))


@pytest.mark.parametrize("F", [2.0, 8.0])
def test_constant_load_run_statistics_kondo(F):
    m = kondo_kif5a()
    sims = simulate_runs(m, Clamp(F), n_runs=4000, seed=3)
    exact = m.run_stats(F)
    within(sims["time"], exact["time"])
    within(sims["x"], exact["displacement"])


def test_constant_load_run_statistics_with_fast_events():
    m = kondo_kif5a(fast=True, q=0.3)
    sims = simulate_runs(m, Clamp(12.0), n_runs=6000, seed=4)
    exact = m.run_stats(12.0)
    within(sims["time"], exact["time"])
    within(sims["x"], exact["displacement"])


@pytest.mark.parametrize("gate,fast", [(1.0, False), (0.1, False), (1.0, True)])
def test_trap_simulation_matches_exact_lattice(gate, fast):
    m = kondo_kif5a(gate=gate, fast=fast)
    trap = Trap(0.05)
    exact = trap_statistics(m, trap)
    sims = simulate_runs(m, trap, n_runs=5000, seed=5)
    assert sims["terminated"].all()
    within(sims["load"], exact["mean_load_at_termination"])
    within(sims["time"], exact["mean_time"])


def test_trap_with_entry_and_delays_matches_exact_lattice():
    """A motor with an entry wait, restarts, post-step delays and a detachment
    exit, in a trap: exercises every phase type of the lattice solver."""
    from competing_exits import Bell, Exit, ExitKind, WaitingState, constant
    s = WaitingState("w", (
        Exit("f", ExitKind.PRODUCTIVE, Bell(400.0, 3.0), step=8.2, cost_time=0.01),
        Exit("b", ExitKind.WRONG, Bell(3.0, 0.0), step=-8.2, cost_time=0.01),
        Exit("c", ExitKind.RESTART, constant(60.0), cost_time=0.004),
        Exit("d", ExitKind.TERMINATING, Bell(1.0, -0.6)),
    ), kT=4.1, entry_time=0.02)
    m = Motor.single(s)
    exact = trap_statistics(m, Trap(0.08))
    for fixed in (False, True):  # only the means of delays matter for these
        sims = simulate_runs(m, Trap(0.08), n_runs=4000, seed=6, fixed_delays=fixed)
        within(sims["load"], exact["mean_load_at_termination"])
        within(sims["time"], exact["mean_time"])


def test_simulate_run_stops_at_t_max():
    m = Motor.single(head_race_v2(1000.0))  # no terminating exit
    run = simulate_run(m, Clamp(0.0), np.random.default_rng(7), t_max=0.5)
    assert not run.terminated
    assert run.time == pytest.approx(0.5)
    assert 0 < run.x < 0.5 * 600  # moved, but no faster than ~v2 allows
