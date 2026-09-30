"""The bet: relative entropies, the I-projection family, and the v3 race."""
import numpy as np
import pytest

from competing_exits import bet
from competing_exits.diffusion import kl_gaussian
from competing_exits.headrace import V3_FIXED, head_race_v3, v3_search


def test_commitment_and_mismatch_match_closed_forms_when_walls_are_far():
    x = np.linspace(-60, 60, 120001)
    p0 = bet.tilted(bet.tether_logdensity(0.21, 0.0), x)
    q = bet.tilted(bet.tether_logdensity(0.4, 0.9), x)
    assert bet.commitment(q, p0, x) == pytest.approx(kl_gaussian(0.9, 1 / 0.4, 0.0, 1 / 0.21), rel=1e-6)
    assert bet.mismatch(p0, q, x) == pytest.approx(kl_gaussian(0.0, 1 / 0.21, 0.9, 1 / 0.4), rel=1e-6)


def test_load_tilts_a_tether_by_f_over_kappa():
    x = np.linspace(-60, 60, 120001)
    kT, F, k = 4.087, 6.0, 0.3
    q = bet.tilted(bet.tether_logdensity(k, 1.0), x, F, kT)
    m, v = bet.moments(q, x)
    assert m == pytest.approx(1.0 - 0.5 * F / kT / k, rel=1e-6)
    assert v == pytest.approx(1 / k, rel=1e-6)


def test_restriction_costs_at_least_the_commitment():
    """For a docking potential V >= 0 (restriction only), the free-energy cost
    -ln <exp(-V)>_p0 of imposing q ~ p0 exp(-V) is >= D(q || p0), with equality for a
    hard wall."""
    x = bet.grid()
    p0 = bet.tilted(bet.tether_logdensity(0.21, 0.0), x)
    rng = np.random.default_rng(1)
    for _ in range(20):
        V = np.abs(rng.normal(0, 1.5)) * np.exp(-rng.uniform(0.2, 1.2) * x) * (rng.random() < 0.5) \
            + np.clip(rng.normal(0, 1) * (x - rng.uniform(-6, 2)), 0, None)
        w = p0 * np.exp(-V)
        cost = -np.log(np.trapezoid(w, x))
        q = bet.normalise(w, x)
        assert cost >= bet.commitment(q, p0, x) - 1e-9
    wall = (x > -2.0).astype(float)
    q = bet.normalise(p0 * wall, x)
    cost = -np.log(np.trapezoid(p0 * wall, x))
    assert cost == pytest.approx(bet.commitment(q, p0, x), rel=1e-3)


def test_restriction_cost_is_the_exact_least_cost_and_bounds_the_commitment():
    """For q ~ p0 exp(-V), V >= 0, the cost -ln<e^-V>_p0 is smallest when min V = 0,
    where it equals ln max(q/p0); that is >= D(q || p0), with equality for a hard wall."""
    x = bet.grid()
    p0 = bet.tilted(bet.tether_logdensity(0.21, 0.0), x)
    rng = np.random.default_rng(2)
    for _ in range(20):
        V = np.abs(rng.normal(0, 2.0)) * np.exp(-rng.uniform(0.2, 1.2) * x) + rng.uniform(0, 1) * (x - rng.uniform(-4, 4)) ** 2 / 8
        V = V - V.min()
        w = p0 * np.exp(-V)
        cost = -np.log(np.trapezoid(w, x))
        q = bet.normalise(w, x)
        assert bet.restriction_cost(q, p0, x) == pytest.approx(cost, rel=1e-9, abs=1e-12)
        assert bet.restriction_cost(q, p0, x) >= bet.commitment(q, p0, x) - 1e-12
    for k, xe in ((0.41, 0.91), (0.10, 4.8), (0.23, 0.4)):
        q = bet.tilted(bet.tether_logdensity(k, xe), x)
        assert bet.restriction_cost(q, p0, x) >= bet.commitment(q, p0, x)
    wall = (x > -2.0).astype(float)
    q = bet.normalise(p0 * wall, x)
    assert bet.restriction_cost(q, p0, x) == pytest.approx(bet.commitment(q, p0, x), rel=1e-3)


def test_iprojection_family_contains_p0_and_reduces_to_it():
    x = bet.grid()
    ld = bet.iprojection_logdensity(0.21, (0.0, 0.0), (0.37, 1.1))
    q = bet.tilted(ld, x)
    p0 = bet.tilted(bet.tether_logdensity(0.21, 0.0), x)
    assert bet.commitment(q, p0, x) == pytest.approx(0.0, abs=1e-12)


def test_v3_exponential_race_matches_exact_race_with_tilted_start():
    """The competing-exits treatment of v3 (capture as a Poisson exit) against the full
    diffusive race with the gate and clock as uniform killing, through and past stall."""
    kT = V3_FIXED["kT"]
    for kappa, x_eq, B in ((0.41, 0.91, 0.0), (0.102, 4.8, 11.26)):
        hs = v3_search(kappa, x_eq, B)
        s = head_race_v3(0, 0, B, 6.4, 0.58, 57.6, 1.07, 14.2e-3, search=hs)
        for F in (0.0, 4.0, 7.0, 10.0, 12.0):
            kb = 6.4 * np.exp(-F * 0.58 / kT)
            ex = hs.exact_race(F, kT, {"gate": kb, "clock": 57.6})
            assert ex.prob["right"] == pytest.approx(s.prob("forward", F), rel=2e-3)
            assert ex.prob["gate"] == pytest.approx(s.prob("back", F), rel=2e-3)


def test_v3_exponential_race_matches_exact_race_at_assisting_and_superstall_loads():
    """The same check from -15 pN (assisting) to +15 pN (superstall), for the stiff
    best-fit tether and the floppy, forward-parked edge of the valley."""
    kT = V3_FIXED["kT"]
    for kappa, x_eq, B in ((0.41, 0.91, 0.0), (0.102, 4.8, 11.26), (0.23, 0.4, 5.37)):
        hs = v3_search(kappa, x_eq, B)
        s = head_race_v3(0, 0, B, 6.4, 0.58, 57.6, 1.07, 14.2e-3, search=hs)
        for F in (-15.0, -10.0, -5.0, -2.0, 13.0, 15.0):
            kb = 6.4 * np.exp(-F * 0.58 / kT)
            ex = hs.exact_race(F, kT, {"gate": kb, "clock": 57.6})
            for exact, approx in ((ex.prob["right"], s.prob("forward", F)), (ex.prob["gate"], s.prob("back", F)),
                                  (ex.prob["clock"], s.prob("clock", F))):
                assert exact == pytest.approx(approx, rel=2e-3, abs=1e-9)
