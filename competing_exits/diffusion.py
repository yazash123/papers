"""First passage of an overdamped particle in a 1D landscape (the diffusive head search).

The particle diffuses with coefficient D (nm^2/s) in a potential U(x) given in
units of kT on an interval [a, b].  Each end is either absorbing (a binding
site: the search ends there) or reflecting.  In addition the search can be
ended at any moment by Poisson *clocks* (uniform killing rates: ATP leaving,
a load-insensitive backstep gate, ...).  This is the continuous analogue of a
``WaitingState`` race: the exits are the absorbing ends and the clocks.

Discretisation.  The interval is cut into N equal cells and the diffusion is
replaced by a birth-death chain on the N + 1 nodes with the Scharfetter-Gummel
rates

    k(i -> i+1) = (D/h^2) B(U[i+1] - U[i]),   k(i+1 -> i) = (D/h^2) B(U[i] - U[i+1]),
    B(z) = z / (e^z - 1),

which obey detailed balance with respect to exp(-U) exactly and converge to the
Fokker-Planck operator with O(h^2) error.  Every quantity is then a quantity of
an absorbing Markov chain (as in ``lattice.py``): with G = (-Q)^-1 on the
transient nodes and an initial distribution p,

    P_j        = p G r_j            (probability that exit j ends the search)
    E[T 1_j]   = p G^2 r_j          (first time moment restricted to exit j)
    E[T^2 1_j] = 2 p G^3 r_j

where r_j is the vector of rates into exit j.  Q is tridiagonal, so each of
these is a banded solve.  ``exact_*`` functions give the continuum formulas
(nested integrals) used by the tests to check the chain.

Units: x in nm, D in nm^2/s, rates in 1/s, U in kT.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Mapping, Optional

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.linalg import solve_banded

ABSORB, REFLECT = "absorb", "reflect"


def bernoulli(z):
    """B(z) = z / (e^z - 1), evaluated stably (B(0) = 1)."""
    z = np.asarray(z, dtype=float)
    out = np.empty_like(z)
    small = np.abs(z) < 1e-6
    out[small] = 1.0 - z[small] / 2.0 + z[small] ** 2 / 12.0
    zz = z[~small]
    out[~small] = zz / np.expm1(zz)
    return out


def log_bernoulli(z):
    """ln B(z) without overflow: ln z - z for large z, ln(-z) for very negative z."""
    z = np.asarray(z, dtype=float)
    out = np.empty_like(z)
    big, neg = z > 30.0, z < -30.0
    mid = ~(big | neg)
    out[big] = np.log(z[big]) - z[big] + np.log1p(np.exp(-z[big]) / (1.0 - np.exp(-z[big])))
    out[neg] = np.log(-z[neg]) + np.log1p(-np.exp(z[neg]))
    out[mid] = np.log(bernoulli(z[mid]))
    return out


def stokes_einstein_D(kT: float, eta_mPa_s: float, radius_nm: float) -> float:
    """D = kT / (6 pi eta r) in nm^2/s (kT in pN nm, eta in mPa s = 1e-3 pN s/nm^2 * 1e-6)."""
    eta = eta_mPa_s * 1e-9  # pN s / nm^2   (1 mPa s = 1e-3 Pa s = 1e-3 N s/m^2 = 1e-9 pN s/nm^2)
    return kT / (6.0 * np.pi * eta * radius_nm)


@dataclass
class RaceResult:
    """Outcome of a diffusive race: exit probabilities and time moments.

    prob[j]   P(exit j ends the search)
    m1[j]     E[T; exit j]   (so E[T | j] = m1[j] / prob[j])
    m2[j]     E[T^2; exit j]
    """

    prob: dict
    m1: dict
    m2: dict

    @property
    def mean_time(self):
        return sum(self.m1.values())

    def conditional_mean(self, j):
        return self.m1[j] / self.prob[j]

    def conditional_cv2(self, j):
        m = self.m1[j] / self.prob[j]
        return (self.m2[j] / self.prob[j]) / m ** 2 - 1.0

    def rate(self, j):
        """Effective exponential rate of exit j: P_j / E[T] (exact for a Poisson race)."""
        return self.prob[j] / self.mean_time


@dataclass
class Search1D:
    """Overdamped search on [a, b] in potential U(x) (kT), discretised on N cells.

    left, right : 'absorb' or 'reflect'
    U           : callable U(x) -> array (kT)
    D           : diffusion coefficient (nm^2/s)
    """

    U: Callable
    D: float
    a: float = -8.0
    b: float = 8.0
    left: str = ABSORB
    right: str = ABSORB
    N: int = 3200
    x: np.ndarray = field(init=False, repr=False)
    Ux: np.ndarray = field(init=False, repr=False)

    def __post_init__(self):
        self.x = np.linspace(self.a, self.b, self.N + 1)
        self.Ux = np.asarray(self.U(self.x), dtype=float)
        h = self.x[1] - self.x[0]
        dU = np.diff(self.Ux)
        c = self.D / h ** 2
        self.k_up = c * bernoulli(dU)      # i -> i+1
        self.k_dn = c * bernoulli(-dU)     # i+1 -> i
        self.log_k_up = np.log(c) + log_bernoulli(dU)
        # transient nodes
        lo = 1 if self.left == ABSORB else 0
        hi = self.N - 1 if self.right == ABSORB else self.N
        self.nodes = np.arange(lo, hi + 1)

    @property
    def h(self):
        return self.x[1] - self.x[0]

    # ----- the chain -------------------------------------------------------
    def _banded(self, kill_total: float):
        """(-Q) on the transient nodes as a banded matrix (1 sub, 1 super).

        A node on a reflecting end owns half a cell (finite-volume view), so its
        outgoing rate is doubled; this makes the wall second-order accurate."""
        idx = self.nodes
        n = len(idx)
        up = np.zeros(n)   # rate i -> i+1 for transient node i
        dn = np.zeros(n)   # rate i -> i-1
        for m, i in enumerate(idx):
            if i < self.N:
                up[m] = self.k_up[i]
            if i > 0:
                dn[m] = self.k_dn[i - 1]
        if self.left == REFLECT:
            up[0] *= 2.0
        if self.right == REFLECT:
            dn[-1] *= 2.0
        diag = up + dn + kill_total
        ab = np.zeros((3, n))
        ab[1] = diag
        # (-Q)[m, m+1] = -rate(m -> m+1) if m+1 is transient
        ab[0, 1:] = -up[:-1]
        # (-Q)[m, m-1] = -rate(m -> m-1)
        ab[2, :-1] = -dn[1:]
        return ab, up, dn

    def exit_vectors(self, clocks: Mapping[str, float]):
        idx = self.nodes
        n = len(idx)
        out = {}
        if self.left == ABSORB:
            v = np.zeros(n)
            v[0] = self.k_dn[0]          # node 1 -> node 0
            out["left"] = v
        if self.right == ABSORB:
            v = np.zeros(n)
            v[-1] = self.k_up[self.N - 1]  # node N-1 -> node N
            out["right"] = v
        for name, k in clocks.items():
            out[name] = np.full(n, float(k))
        return out

    def start_vector(self, start):
        """Initial distribution over the transient nodes.

        start: a position (float), or a callable/array giving a density on the grid
        (it is sampled at the nodes and renormalised over the transient nodes)."""
        n = len(self.nodes)
        p = np.zeros(n)
        if np.isscalar(start):
            i = int(np.argmin(np.abs(self.x[self.nodes] - start)))
            p[i] = 1.0
            return p
        dens = start(self.x) if callable(start) else np.asarray(start, dtype=float)
        p = np.clip(dens[self.nodes], 0.0, None)
        # trapezoid weights: a node on a reflecting wall owns half a cell
        if self.left == REFLECT:
            p[0] *= 0.5
        if self.right == REFLECT:
            p[-1] *= 0.5
        s = p.sum()
        if s <= 0:
            raise ValueError("start distribution has no mass on the transient nodes")
        return p / s

    def capture(self, start) -> dict:
        """Exit probabilities and mean exit time with no clocks, from the exact
        birth-death formulas (machine precision; no linear solve).

        With stationary weights pi_i (exp(-U), halved on a reflecting wall node) and
        edge resistances R_j = 1 / (pi_j k(j -> j+1)), the committor to the right end is
        h(i) = sum_{j<i} R_j / sum_j R_j and, with I_j = sum_{k<=j, transient} pi_k,
        T(i) = h(i) sum_j R_j I_j - sum_{j<i} R_j I_j  (both ends absorbing), or
        T(i) = sum_{j>=i} R_j I_j                        (left end reflecting).
        These are the discrete forms of the continuum integrals."""
        if self.right != ABSORB:
            raise ValueError("capture() needs an absorbing right end")
        U = self.Ux
        # log stationary weights (up to a constant) and log resistances of edges j -> j+1
        lpi = -(U - U.min())
        if self.left == REFLECT:
            lpi = lpi.copy()
            lpi[0] += np.log(0.5)
        lup = self.log_k_up.copy()
        if self.left == REFLECT:
            lup[0] += np.log(2.0)
        lR = -(lpi[:-1] + lup)                   # edges 0..N-1 (log space: no underflow)
        # I_j = sum of pi over transient nodes k <= j (edge j spans nodes j, j+1)
        first = 1 if self.left == ABSORB else 0
        pi = np.exp(lpi)
        cum = np.cumsum(pi)
        I = cum - (cum[0] if first == 1 else 0.0)  # I[0] = 0 when node 0 is absorbing
        I = I[:-1]                                 # edges 0..N-1
        shift = lR.max()
        Rs = np.exp(lR - shift)                    # scaled resistances
        RI = Rs * I
        cR = np.concatenate([[0.0], np.cumsum(Rs)])     # sum_{j<i} R_j, i = 0..N
        cRI = np.concatenate([[0.0], np.cumsum(RI)])
        nodes = self.nodes
        p = self.start_vector(start)
        if self.left == ABSORB:
            h = cR[nodes] / cR[-1]
            T = (h * cRI[-1] - cRI[nodes]) * np.exp(shift)
            pr = float(p @ h)
            return {"right": pr, "left": 1.0 - pr, "mean_time": float(p @ T)}
        T = (cRI[-1] - cRI[nodes]) * np.exp(shift)
        return {"right": 1.0, "mean_time": float(p @ T)}

    def race(self, start, clocks: Optional[Mapping[str, float]] = None) -> RaceResult:
        """Exit probabilities and time moments for a search started from ``start``,
        ended by an absorbing end ('left'/'right') or one of the Poisson ``clocks``."""
        clocks = dict(clocks or {})
        ab, _, _ = self._banded(sum(clocks.values()))
        p = self.start_vector(start)
        # Backward (committor) solves: y_j = G r_j is the probability of leaving by
        # exit j from each node, z_j = G y_j and w_j = G z_j give the time moments.
        # The chain is stiff (fast diffusion, slow escape: condition number ~1e9), so
        # probabilities carry absolute errors of ~1e-8; this is far below anything
        # used in the analysis.
        prob, m1, m2 = {}, {}, {}
        for name, r in self.exit_vectors(clocks).items():
            y = solve_banded((1, 1), ab, r)
            z = solve_banded((1, 1), ab, y)
            w = solve_banded((1, 1), ab, z)
            prob[name] = float(p @ y)
            m1[name] = float(p @ z)
            m2[name] = 2.0 * float(p @ w)
        return RaceResult(prob, m1, m2)

    # ----- equilibrium and relative entropy --------------------------------
    def boltzmann(self, U: Optional[Callable] = None):
        """Normalised equilibrium density exp(-U) on the grid (per nm)."""
        Uv = self.Ux if U is None else np.asarray(U(self.x), dtype=float)
        w = np.exp(-(Uv - Uv.min()))
        return w / np.trapezoid(w, self.x)


# ---------------------------------------------------------------------------
# Continuum formulas (for tests and checks)
# ---------------------------------------------------------------------------
def exact_splitting_mfpt(x, Ux, D, x0):
    """Both ends absorbing.  Returns (P_right, mean exit time) from x0, using the
    nested-integral solution of D e^U (e^-U T')' = -1 (trapezoid quadrature)."""
    Umax, Umin = Ux.max(), Ux.min()
    eU = np.exp(Ux - Umax)
    emU = np.exp(-(Ux - Umin))
    Iu = cumulative_trapezoid(eU, x, initial=0.0)
    I = cumulative_trapezoid(emU, x, initial=0.0)
    J = cumulative_trapezoid(eU * I, x, initial=0.0)
    iu = np.interp(x0, x, Iu)
    jj = np.interp(x0, x, J)
    pb = iu / Iu[-1]
    T = (pb * J[-1] - jj) * np.exp(Umax - Umin) / D
    return pb, T


def exact_mfpt_reflect_left(x, Ux, D, x0):
    """Reflecting at x[0], absorbing at x[-1]: T(x0) = (1/D) int_x0^b e^U(y) int_a^y e^-U."""
    Umax, Umin = Ux.max(), Ux.min()
    eU = np.exp(Ux - Umax)
    emU = np.exp(-(Ux - Umin))
    I = cumulative_trapezoid(emU, x, initial=0.0)
    J = cumulative_trapezoid(eU * I, x, initial=0.0)
    return (J[-1] - np.interp(x0, x, J)) * np.exp(Umax - Umin) / D


# ---------------------------------------------------------------------------
# Distributions and relative entropy on a grid
# ---------------------------------------------------------------------------
def gaussian_density(mean, var):
    def dens(x):
        return np.exp(-(np.asarray(x) - mean) ** 2 / (2.0 * var)) / np.sqrt(2 * np.pi * var)
    return dens


def relative_entropy(p, q, x):
    """D(p || q) in nats for densities sampled on the grid x (trapezoid)."""
    p = np.asarray(p, dtype=float)
    q = np.asarray(q, dtype=float)
    p = p / np.trapezoid(p, x)
    q = q / np.trapezoid(q, x)
    mask = p > 0
    integrand = np.zeros_like(p)
    integrand[mask] = p[mask] * (np.log(p[mask]) - np.log(np.maximum(q[mask], 1e-300)))
    return float(np.trapezoid(integrand, x))


def kl_gaussian(m1, v1, m0, v0):
    """D(N(m1, v1) || N(m0, v0)) in nats (closed form)."""
    return 0.5 * (np.log(v0 / v1) + (v1 + (m1 - m0) ** 2) / v0 - 1.0)
