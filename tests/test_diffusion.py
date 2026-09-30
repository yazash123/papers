"""The 1D diffusive first-passage solver against continuum formulas."""
import numpy as np
import pytest

from competing_exits.diffusion import (
    Search1D, bernoulli, exact_mfpt_reflect_left, exact_splitting_mfpt, gaussian_density,
    kl_gaussian, relative_entropy, stokes_einstein_D,
)

D = 87.3e6  # nm^2/s


def v1_like(F=0.0, kT=4.114):
    f = 0.5 * F / kT
    return lambda x: (0.15 * (x - 2.5) ** 2 + f * x
                      + 9.144 * np.exp(-(x - 6) ** 2 / 0.98) + 9.144 * np.exp(-(x + 6) ** 2 / 0.98))


def test_stokes_einstein_gives_87_3_nm2_per_us():
    assert stokes_einstein_D(4.114, 1.0, 2.5) / 1e6 == pytest.approx(87.3, rel=1e-3)


def test_bernoulli_function():
    z = np.array([-30.0, -1.0, -1e-8, 0.0, 1e-8, 1.0, 30.0])
    ref = np.where(z == 0, 1.0, z / np.expm1(np.where(z == 0, 1.0, z)))
    assert np.allclose(bernoulli(z), ref, rtol=1e-12)
    # detailed balance: B(z)/B(-z) = exp(-z)
    assert np.allclose(bernoulli(z) / bernoulli(-z), np.exp(-z))


@pytest.mark.parametrize("F", [0.0, 3.0, 6.0])
def test_chain_matches_exact_splitting_and_mfpt(F):
    x = np.linspace(-8, 8, 400001)
    pb, T = exact_splitting_mfpt(x, v1_like(F)(x), D, -4.0)
    s = Search1D(v1_like(F), D, N=6400)
    c = s.capture(-4.0)                     # exact birth-death formulas
    assert c["right"] == pytest.approx(pb, rel=1e-7)
    assert c["mean_time"] == pytest.approx(T, rel=2e-5)
    r = s.race(-4.0)                        # linear solves (stiff: ~1e-6 accuracy)
    assert r.prob["right"] == pytest.approx(c["right"], rel=1e-5)
    assert r.mean_time == pytest.approx(c["mean_time"], rel=1e-5)
    assert r.prob["left"] + r.prob["right"] == pytest.approx(1.0, abs=1e-5)


def test_chain_converges_second_order():
    x = np.linspace(-8, 8, 400001)
    _, T = exact_splitting_mfpt(x, v1_like()(x), D, -4.0)
    e1 = abs(Search1D(v1_like(), D, N=1600).race(-4.0).mean_time / T - 1)
    e2 = abs(Search1D(v1_like(), D, N=3200).race(-4.0).mean_time / T - 1)
    assert e1 / e2 == pytest.approx(4.0, rel=0.1)


def test_reflecting_left_end():
    x = np.linspace(-8, 8, 400001)
    T = exact_mfpt_reflect_left(x, v1_like()(x), D, 0.0)
    s = Search1D(v1_like(), D, left="reflect", N=6400)
    assert s.capture(0.0)["mean_time"] == pytest.approx(T, rel=2e-5)
    assert s.race(0.0).mean_time == pytest.approx(T, rel=2e-5)
    assert s.race(0.0).prob["right"] == pytest.approx(1.0, abs=1e-6)


def test_reflecting_wall_is_exact_for_free_diffusion():
    """Half-cell treatment of the wall: E[T] = L^2/(2D) to machine precision."""
    L = 10.0
    for N in (100, 1000):
        s = Search1D(lambda x: 0 * x, D, a=0.0, b=L, left="reflect", N=N)
        assert s.capture(0.0)["mean_time"] == pytest.approx(L ** 2 / (2 * D), rel=1e-12)


def test_free_diffusion_with_a_clock_matches_closed_form():
    """U = 0 on [0, L], both ends absorbing, a Poisson clock at rate s:
    P(right before clock | x) = sinh(q x)/sinh(q L), q = sqrt(s/D)."""
    L, s, x0 = 10.0, 3e5, 4.0
    q = np.sqrt(s / D)
    r = Search1D(lambda x: 0 * x, D, a=0.0, b=L, N=4000).race(x0, {"clock": s})
    assert r.prob["right"] == pytest.approx(np.sinh(q * x0) / np.sinh(q * L), rel=1e-4)
    assert sum(r.prob.values()) == pytest.approx(1.0, abs=1e-10)
    # E[T; any exit] = (1 - P_absorbed_by_ends ... ) : E[min(tau, T_clock)] = P_clock / s
    assert r.mean_time == pytest.approx(r.prob["clock"] / s, rel=1e-6)


def test_time_moments_of_free_diffusion():
    """Reflecting at 0, absorbing at L, start at 0: E[T] = L^2/(2D), E[T^2] = 5 L^4/(12 D^2)."""
    L = 10.0
    r = Search1D(lambda x: 0 * x, D, a=0.0, b=L, left="reflect", N=4000).race(0.0)
    assert r.m1["right"] == pytest.approx(L ** 2 / (2 * D), rel=1e-9)
    assert r.m2["right"] == pytest.approx(5 * L ** 4 / (12 * D ** 2), rel=1e-5)


def test_start_distribution_averages_point_starts():
    s = Search1D(v1_like(), D, N=3200)
    dens = gaussian_density(0.0, 2.0)
    avg = s.race(dens).mean_time
    xs = s.x[s.nodes]
    w = dens(xs) / dens(xs).sum()
    pointwise = sum(wi * s.race(float(xi)).mean_time for wi, xi in zip(w[::200], xs[::200])) / w[::200].sum()
    assert avg == pytest.approx(pointwise, rel=1e-3)


def test_relative_entropy_of_gaussians():
    x = np.linspace(-40, 40, 80001)
    p = gaussian_density(1.0, 2.0)(x)
    q = gaussian_density(-0.5, 5.0)(x)
    assert relative_entropy(p, q, x) == pytest.approx(kl_gaussian(1.0, 2.0, -0.5, 5.0), rel=1e-6)
    assert relative_entropy(p, p, x) == pytest.approx(0.0, abs=1e-12)
