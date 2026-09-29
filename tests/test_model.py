"""Analytic identities of the competing-exits model (no kinesin numbers)."""
import numpy as np
import pytest

from competing_exits import (Bell, Clamp, Exit, ExitKind, Motor, PiecewiseBell, Scaled,
                             WaitingState, constant, kT_at, trap_statistics)

P, W, T, R = ExitKind.PRODUCTIVE, ExitKind.WRONG, ExitKind.TERMINATING, ExitKind.RESTART
KT = 4.1


def toy(entry=0.01, restart_cost=0.005, post=0.02, with_term=True):
    exits = [
        Exit("a", P, Bell(300.0, 3.0), step=+8.0, cost_time=post),
        Exit("b", W, Bell(4.0, -0.5), step=-8.0, cost_time=post),
        Exit("r", R, constant(50.0), cost_time=restart_cost),
    ]
    if with_term:
        exits.append(Exit("t", T, Bell(1.0, -0.6)))
    return WaitingState("s", tuple(exits), kT=KT, entry_time=entry)


def test_rate_laws():
    assert kT_at(25.0) == pytest.approx(4.1164, abs=1e-4)
    law = Bell(10.0, 2.0)
    assert law(0.0, KT) == 10.0
    assert law(KT, KT) == pytest.approx(10.0 * np.exp(-2.0))
    pw = PiecewiseBell(5.0, 3.0, delta_below=0.2, delta_above=2.0)
    assert pw(3.0 - 1e-9, KT) == pytest.approx(pw(3.0, KT))  # continuous at the knee
    assert pw(0.0, KT) == pytest.approx(5.0 * np.exp(3.0 * 0.2 / KT))
    assert pw(5.0, KT) == pytest.approx(5.0 * np.exp(-2.0 * 2.0 / KT))
    assert Scaled(law, 0.1)(1.0, KT) == pytest.approx(0.1 * law(1.0, KT))


@pytest.mark.parametrize("F", [0.0, 2.5, 7.0])
def test_splitting_and_prices(F):
    s = toy()
    P_ = s.splitting(F)
    assert sum(P_.values()) == pytest.approx(1.0)
    k = {e.name: e.k(F, KT) for e in s.exits}
    for name in k:
        assert P_[name] == pytest.approx(k[name] / sum(k.values()))
    # log-odds are a difference of prices (bookkeeping identity)
    assert s.log_odds("a", "b", F) == pytest.approx(s.price("b", F) - s.price("a", F))
    assert s.attempts_per("a", F) == pytest.approx(np.exp(s.price("a", F)))


def test_balance_point_and_lever():
    s = toy()
    Fs = KT * np.log(300.0 / 4.0) / (3.0 - (-0.5))
    assert s.balance_point("a", "b") == pytest.approx(Fs, rel=1e-10)
    assert s.lever("a", "b", F=1.0) == pytest.approx(3.5, rel=1e-6)
    assert s.log_odds("a", "b", Fs) == pytest.approx(0.0, abs=1e-10)


def test_dwell_and_velocity_closed_form():
    s = toy(with_term=False)
    F = 2.0
    ka, kb, kr = 300 * np.exp(-F * 3 / KT), 4 * np.exp(F * 0.5 / KT), 50.0
    K = ka + kb + kr
    Pr, Pst = kr / K, (ka + kb) / K
    attempt = 0.01 + 1 / K + Pr * 0.005
    dwell = attempt / Pst + 0.02
    assert s.mean_dwell(F) == pytest.approx(dwell)
    assert s.velocity(F) == pytest.approx(8.0 * (ka - kb) / (ka + kb) / dwell)


def test_single_state_motor_matches_state_formulas():
    s = toy()
    m = Motor.single(s, step_size=8.0)
    for F in (0.0, 3.0, 8.0):
        assert m.velocity(F) == pytest.approx(s.velocity(F))
        a, b = m.run_stats(F), s.run_stats(F)
        for key in a:
            assert a[key] == pytest.approx(b[key])


def test_run_stats_closed_form():
    s = WaitingState("s", (
        Exit("f", P, constant(80.0), step=8.0),
        Exit("b", W, constant(2.0), step=-8.0),
        Exit("d", T, constant(1.0)),
    ), kT=KT)
    r = s.run_stats(0.0)
    assert r["steps"] == pytest.approx(82.0)
    assert r["displacement"] == pytest.approx(8.0 * 78.0)
    assert r["time"] == pytest.approx(1.0)       # 1/k_d
    assert s.velocity() == pytest.approx(8.0 * 78.0)   # d (k_f - k_b)


def test_chained_state_first_step_analysis():
    """After a 'b' exit the system enters a fast state with probability q; from
    there it takes extra 'b2' steps (again re-entering with probability q) or
    terminates.  Expected extra wrong steps per 'b' = q r / (1 - q r)."""
    q = 0.3
    s0 = WaitingState("0", (
        Exit("f", P, constant(50.0), step=8.0),
        Exit("b", W, constant(5.0), step=-8.0, to={"x": q, None: 1 - q}),
        Exit("d", T, constant(0.5)),
    ), kT=KT)
    sx = WaitingState("x", (
        Exit("b2", W, constant(600.0), step=-8.0, to={"x": q, None: 1 - q}),
        Exit("d2", T, constant(400.0)),
    ), kT=KT)
    m = Motor((s0, sx), home="0", step_size=8.0)
    c = m.cycle_stats(0.0)
    r = 0.6
    extra = q * r / (1 - q * r)
    Pf, Pb, Pd = 50 / 55.5, 5 / 55.5, 0.5 / 55.5
    assert c["steps"] == pytest.approx(Pf + Pb * (1 + extra))
    assert c["displacement"] == pytest.approx(8.0 * (Pf - Pb * (1 + extra)))
    p_fast_term = q * (1 - r) / (1 - q * r)
    assert c["p_term"] == pytest.approx(Pd + Pb * p_fast_term)
    t_fast = 1 / 1000.0 / (1 - q * r)  # expected time in the fast state per entry
    assert c["time"] == pytest.approx(1 / 55.5 + Pb * q * t_fast)


def test_lattice_solver_under_clamp_matches_constant_load_formulas():
    """With a constant force the lattice solver must reproduce run_stats."""
    s = toy()
    m = Motor.single(s, step_size=8.0)
    for F in (1.0, 5.0):
        r = trap_statistics(m, Clamp(F), n_max=2500)
        a = m.run_stats(F)
        assert r["absorbed"] == pytest.approx(1.0, abs=1e-6)  # truncation tail
        assert r["mean_time"] == pytest.approx(a["time"], rel=1e-6)
        assert r["mean_load_at_termination"] == pytest.approx(F)
        n_f, n_b = r["mean_counts"]["s.a"], r["mean_counts"]["s.b"]
        assert 8.0 * (n_f - n_b) == pytest.approx(a["displacement"], rel=1e-6)


def test_gate_scaling():
    s = toy()
    s10 = s.scaled("b", 0.1)
    assert s10.rate("b", 2.0) == pytest.approx(0.1 * s.rate("b", 2.0))
    # a 10x stronger gate moves the balance point out by kT ln 10 / lever
    shift = s10.balance_point("a", "b") - s.balance_point("a", "b")
    assert shift == pytest.approx(KT * np.log(10) / 3.5)
    with pytest.raises(KeyError):
        s.scaled("nope", 2.0)
