"""The docking 'bet' on where the free head will be, and what it costs.

Three distributions of the free head's position x (nm, relative to the bound
head), all on the search interval [-8, 8] nm with reflecting walls:

    p0  the reference, "no bet": the head on the undocked tether, a floppy
        harmonic tether centred on the bound head, p0 ~ exp(-kappa0 x^2 / 2);
    q   the bet: the head's equilibrium distribution in the docked landscape;
    p   where the head actually is when docking happens (from tracking).

Under a hindering load F every distribution is tilted by exp(-f x), f = F/(2 kT)
(the head feels half the stalk load).

Costs (in nats = kT):
    commitment  D(q || p0): the least free energy that imposing q on a head that
                would otherwise follow p0 can cost, when docking acts only by
                restricting the head (a repulsive or excluding potential V >= 0):
                -ln <e^-V>_p0 >= D(q || p0), with equality for hard walls.
    mismatch    D(p || q): the free energy dissipated when a head distributed as
                p relaxes in the docked landscape q after the quench (the
                "mismatch cost" of betting q when the truth is p).
    reverse     D(q || p): dissipated if the bet is abandoned (undocking, e.g. ATP
                leaves) after the head has relaxed into q.

The least-committed docked landscape consistent with a set of expectations
<phi_k>_q = c_k is the I-projection of p0: q ~ p0 exp(-sum_k lambda_k phi_k(x)).
For forward capture over a high barrier the load dependence of the capture rate
measures <exp(-f x)>_q (see REPORT2.md), so the natural family is
phi_k(x) = exp(-f_k x) at the loads the motor meets.
"""
from __future__ import annotations

import numpy as np

from .diffusion import relative_entropy

X_LO, X_HI = -8.0, 8.0


def grid(n: int = 16001):
    return np.linspace(X_LO, X_HI, n)


def normalise(w, x):
    w = np.asarray(w, dtype=float)
    return w / np.trapezoid(w, x)


def tilted(log_density, x, F=0.0, kT=4.087):
    """exp(log_density(x) - f x) normalised on the grid (f = F / (2 kT))."""
    f = 0.5 * F / kT
    lw = log_density(x) - f * x
    return normalise(np.exp(lw - lw.max()), x)


def tether_logdensity(kappa, x_eq):
    return lambda x: -0.5 * kappa * (np.asarray(x) - x_eq) ** 2


def iprojection_logdensity(kappa0, lambdas, fs, x0=0.0):
    """log of p0(x) exp(-sum_k lambda_k exp(-f_k x)), p0 = N(x0, 1/kappa0) (unnormalised)."""
    lambdas = np.atleast_1d(lambdas)
    fs = np.atleast_1d(fs)

    def ld(x):
        x = np.asarray(x, dtype=float)
        out = -0.5 * kappa0 * (x - x0) ** 2
        for lam, f in zip(lambdas, fs):
            out = out - lam * np.exp(-f * x)
        return out
    return ld


def commitment(q, p0, x):
    """D(q || p0) in nats."""
    return relative_entropy(q, p0, x)


def mismatch(p, q, x):
    """D(p || q) in nats."""
    return relative_entropy(p, q, x)


def moments(dens, x):
    m = np.trapezoid(x * dens, x)
    v = np.trapezoid((x - m) ** 2 * dens, x)
    return float(m), float(v)
