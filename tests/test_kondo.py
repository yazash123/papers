"""Part C checks: mouse KIF5A rates of Kondo et al. (2023), Table 1."""
import numpy as np
import pytest

from competing_exits import Trap, trap_statistics
from competing_exits.kinesin import kondo_kif5a


def test_unloaded_odds_44_to_1():
    s = kondo_kif5a().home_state
    assert np.exp(s.log_odds("forward", "back", 0.0)) == pytest.approx(44.0, rel=0.01)


def test_backstep_rate_barely_feels_load():
    s = kondo_kif5a().home_state
    assert s.rate("back", 0.0) == pytest.approx(1.99)
    assert s.rate("back", 6.0) / s.rate("back", 0.0) == pytest.approx(1.9, rel=0.01)


def test_forward_rate_flat_to_knee_then_falls():
    s = kondo_kif5a().home_state
    k = [float(s.rate("forward", F)) for F in (0.0, 3.19, 6.0)]
    assert k[0] / k[1] == pytest.approx(np.exp(0.17 * 3.19 / s.kT))  # gentle below
    assert k[1] == pytest.approx(77.3)
    assert k[1] / k[2] == pytest.approx(np.exp(1.78 * (6.0 - 3.19) / s.kT))  # steep above


def test_fitted_rates_cross_near_9_3_pN():
    s = kondo_kif5a().home_state
    assert s.balance_point("forward", "back") == pytest.approx(9.3, abs=0.06)


def test_unloaded_velocity_707():
    """Kondo et al. p.471: 707 nm/s at no load = d (k_f - k_b)."""
    assert kondo_kif5a().velocity(0.0) == pytest.approx(707.0, rel=0.002)


def test_full_model_zero_velocity_near_9_1_pN():
    """With fast backsteps at ~15% of slow ones (Kondo p.470) the zero-velocity
    load moves from the 9.34-pN crossing to ~9.1 pN (Kondo p.471: 9.1)."""
    m = kondo_kif5a(fast=True, q=0.22)
    c = m.cycle_stats(5.0)
    s = m.home_state
    n_fast = m.state("fast")
    r = n_fast.prob("fast_back", 5.0)
    assert 0.22 * r / (1 - 0.22 * r) == pytest.approx(0.15, abs=0.01)
    assert m.stall_force() == pytest.approx(9.1, abs=0.1)
    assert kondo_kif5a().stall_force() == pytest.approx(s.balance_point("forward", "back"))


# Earlier, simpler calculation (0.05 pN/nm trap, starting unloaded):
# gate: (1:1 load, mean load at letting go, mean time attached)
EARLIER = {100: (17.9, 6.7, 0.52), 10: (13.6, 6.7, 0.52), 1: (9.3, 6.4, 0.53),
           1 / 3: (7.3, 5.8, 0.56), 0.1: (5.1, 4.5, 0.64)}
# Exact values from the lattice solver (this code), pinned as regression targets.
EXACT = {100: (17.883, 6.6557, 0.50896), 10: (13.613, 6.6205, 0.51054), 1: (9.343, 6.2943, 0.52576),
         1 / 3: (7.306, 5.7051, 0.55570), 0.1: (5.074, 4.3997, 0.63116)}


@pytest.mark.parametrize("gate", list(EARLIER))
def test_gate_table(gate):
    m = kondo_kif5a(gate=gate)
    F11 = m.home_state.balance_point("forward", "back")
    r = trap_statistics(m, Trap(0.05))
    assert r["absorbed"] == pytest.approx(1.0, abs=1e-9)
    ex = EXACT[gate]
    assert F11 == pytest.approx(ex[0], abs=0.001)
    assert r["mean_load_at_termination"] == pytest.approx(ex[1], abs=0.001)
    assert r["mean_time"] == pytest.approx(ex[2], abs=0.0005)
    # agreement with the earlier calculation: 1:1 loads to rounding; loads and
    # times within the Monte Carlo error of ~1000 runs (SD of the load at
    # detachment is ~2.4 pN, so SE ~0.08 pN; SD of time ~0.5 s, SE ~0.016 s)
    old = EARLIER[gate]
    assert F11 == pytest.approx(old[0], abs=0.05)
    assert r["mean_load_at_termination"] == pytest.approx(old[1], abs=0.15)
    assert r["mean_time"] == pytest.approx(old[2], abs=0.02)
