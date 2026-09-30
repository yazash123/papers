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
    commitment  D(q || p0): the relative entropy of the bet; the quantity the
                MaxEnt (minimum-relative-entropy) member minimises.  It is a LOWER
                bound on the free energy that docking must supply to impose q when
                docking acts only by restricting the head (V >= 0):
                -ln <e^-V>_p0 = D(q || p0) + <V>_q >= D(q || p0), equality for
                hard walls.
    restriction ln max_x q(x)/p0(x) (the max-divergence): the EXACT least free
                energy -ln <e^-V>_p0 over all V >= 0 that produce q (V is fixed by q
                up to a constant; the constant is smallest when min V = 0).  Always
                >= D(q || p0).  This, not D, is the cost to compare with a docking
                budget (added after the session-2 review).
    mismatch    D(p || q): the free energy dissipated when a head distributed as
                p relaxes in the docked landscape q after the quench (the
                "mismatch cost" of betting q when the truth is p).
    reverse     D(q || p0): dissipated if the bet is abandoned (undocking, e.g. ATP
                leaves) after the head has relaxed into q, as it relaxes back to p0.

The least-committed docked landscape consistent with a set of expectations
<phi_k>_q = c_k is the I-projection of p0: q ~ p0 exp(-sum_k lambda_k phi_k(x)).
For forward capture over a high barrier, with a sharp load-independent transition
state at x_TS and the same F/2 lever in well and transition state,
ln k_f(F) = const - f x_TS - ln <exp(-f x)>_q, so the capture rates constrain
<exp(-f x)>_q and the natural family is phi_k(x) = exp(-f_k x) (REPORT2.md, Part B).
The fitted x_TS moves with the binding barrier (6.1-7.7 nm across the valley), so
this two-parameter family is a heuristic MaxEnt family, not the exact one.
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


def restriction_cost(q, p0, x=None):
    """ln max q/p0 over the support of q (nats): the least free energy -ln<e^-V>_p0
    over docking potentials V >= 0 that turn p0 into q.  `x` is unused (kept for a
    signature parallel to commitment)."""
    q = np.asarray(q, dtype=float)
    p0 = np.asarray(p0, dtype=float)
    m = q > 0
    return float(np.max(np.log(q[m]) - np.log(p0[m])))


def mismatch(p, q, x):
    """D(p || q) in nats."""
    return relative_entropy(p, q, x)


def moments(dens, x):
    m = np.trapezoid(x * dens, x)
    v = np.trapezoid((x - m) ** 2 * dens, x)
    return float(m), float(v)
