"""WaitingState.transport_stats: randomness and dwell statistics."""
import numpy as np
import pytest

from competing_exits import Clamp, Exit, ExitKind, Motor, WaitingState, constant, simulate_run
from competing_exits.kinesin import head_race_v2

P, W, R = ExitKind.PRODUCTIVE, ExitKind.WRONG, ExitKind.RESTART


def test_poisson_stepper_has_randomness_one():
    s = WaitingState("s", kT=4.1, exits=(Exit("f", P, constant(50.0), step=8.0),))
    t = s.transport_stats()
    assert t["v"] == pytest.approx(400.0)
    assert t["randomness"] == pytest.approx(1.0)
    assert t["dwell_cv"] == pytest.approx(1.0)


def test_two_equal_sequential_steps_give_randomness_one_half():
    s = WaitingState("s", kT=4.1, entry_time=0.02, exits=(Exit("f", P, constant(50.0), step=8.0),))
    t = s.transport_stats()
    assert t["randomness"] == pytest.approx(0.5)
    assert t["dwell_cv"] == pytest.approx(np.sqrt(0.5))


def test_fixed_rest_time_lowers_randomness():
    s = WaitingState("s", kT=4.1, exits=(Exit("f", P, constant(1e9), step=8.0, cost_time=0.02, cost_cv2=0.0),))
    assert s.transport_stats()["randomness"] == pytest.approx(0.0, abs=1e-6)


def test_biased_random_walk():
    a, b, d = 30.0, 10.0, 8.0
    s = WaitingState("s", kT=4.1, exits=(Exit("f", P, constant(a), step=d), Exit("b", W, constant(b), step=-d)))
    t = s.transport_stats()
    assert t["v"] == pytest.approx(d * (a - b))
    assert t["D"] == pytest.approx(d * d * (a + b) / 2)
    assert t["randomness"] == pytest.approx((a + b) / (a - b))


def test_transport_matches_gillespie_for_v2():
    """Head race v2 at 4 pN (restarts, backsteps, exponential entry and rest times)."""
    s = head_race_v2(1000.0)
    F = 4.0
    t = s.transport_stats(F)
    m = Motor.single(s, step_size=8.2)
    rng = np.random.default_rng(3)
    tmax, n = 6.0, 2500
    xs = np.array([simulate_run(m, Clamp(F), rng, t_max=tmax).x for _ in range(n)])
    v_sim = xs.mean() / tmax
    D_sim = xs.var() / (2 * tmax)
    assert v_sim == pytest.approx(t["v"], rel=0.03)
    assert D_sim == pytest.approx(t["D"], rel=0.08)


def test_dwell_statistics_match_gillespie():
    s = head_race_v2(1000.0)
    F = 5.0
    t = s.transport_stats(F)
    m = Motor.single(s, step_size=8.2)
    rng = np.random.default_rng(5)
    run, path = simulate_run(m, Clamp(F), rng, t_max=400.0, trace=True)
    times = np.array([p[0] for p in path[1:]])
    dw = np.diff(times)
    assert dw.mean() == pytest.approx(t["dwell_mean"], rel=0.05)
    assert dw.std() / dw.mean() == pytest.approx(t["dwell_cv"], rel=0.05)
