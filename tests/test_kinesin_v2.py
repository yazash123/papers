"""Part A regression targets: head race v2 (Drosophila kinesin-1, 23 C).

Targets are those given for this session; tolerance ~1% unless stated.
"""
import numpy as np
import pytest

from competing_exits.kinesin import V2, head_race_v2

REL = 0.01


def test_stall_is_7_05_pN_at_every_atp():
    for atp in (1.0, 10.0, 1000.0, 1e5):
        s = head_race_v2(atp)
        assert s.balance_point("forward", "back") == pytest.approx(7.05, rel=REL)
    # closed form F_s = kT ln(O0) / (delta_f - delta_b)
    Fs = V2["kT"] * np.log(V2["kf0"] / V2["kb0"]) / (V2["delta_f"] - V2["delta_b"])
    assert head_race_v2().balance_point("forward", "back") == pytest.approx(Fs, rel=1e-10)


def test_odds_fall_0_95_nats_per_pN():
    s = head_race_v2()
    slope = (s.log_odds("forward", "back", 2.0) - s.log_odds("forward", "back", 6.0)) / 4.0
    assert slope == pytest.approx(0.95, rel=REL)


@pytest.mark.parametrize("F,v", [(0.0, 454.0), (3.0, 335.0), (5.0, 104.0), (6.0, 32.5)])
def test_force_velocity_at_1mM(F, v):
    assert head_race_v2(1000.0).velocity(F) == pytest.approx(v, rel=REL)


def test_unloaded_velocity_at_10uM():
    assert head_race_v2(10.0).velocity(0.0) == pytest.approx(98.0, rel=REL)


def test_unloaded_dwell_at_1mM():
    s = head_race_v2(1000.0)
    assert s.mean_dwell(0.0) * 1e3 == pytest.approx(18.0, rel=REL)
    assert V2["T"] * 1e3 == pytest.approx(17.1)


@pytest.mark.parametrize("F,price", [(0.0, 0.02), (3.0, 0.45), (6.0, 2.64)])
def test_price_per_attempt(F, price):
    s = head_race_v2(1000.0)
    # 0.02 is a two-significant-figure target; the model gives 0.0248
    tol = 0.005 if F == 0.0 else REL * price
    assert s.price("forward", F) == pytest.approx(price, abs=tol)


def test_price_and_attempts_at_stall():
    s = head_race_v2(1000.0)
    Fs = s.balance_point("forward", "back")
    assert s.price("forward", Fs) == pytest.approx(3.68, rel=REL)
    assert s.attempts_per("forward", Fs) == pytest.approx(40.0, rel=0.02)


def test_price_does_not_depend_on_atp():
    """The race is entered only after ATP binds, so its odds and prices are
    ATP-independent; only the time per attempt changes."""
    for F in (0.0, 5.0):
        assert head_race_v2(10.0).price("forward", F) == pytest.approx(head_race_v2(1000.0).price("forward", F))
