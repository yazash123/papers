"""Load-dependent rate laws.

Sign convention used throughout the package: F > 0 is a *hindering* load
(pN), and a rate with load distance ``delta`` (nm) obeys

    k(F) = k0 * exp(-F * delta / kT)

so delta > 0 means the rate is slowed by hindering load, delta < 0 means it
is sped up.  kT is in pN*nm (4.116 at 25 C, 4.087 at 23 C).

A rate law is any object with ``__call__(F, kT) -> rate``.  Everything else in
the package only relies on that, so new shapes (e.g. a saturating law) can be
added without touching the rest of the code.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

K_B = 1.380649e-2  # Boltzmann constant in pN*nm/K


def kT_at(celsius: float) -> float:
    """Thermal energy in pN*nm at a temperature in degrees Celsius."""
    return K_B * (celsius + 273.15)


@dataclass(frozen=True)
class Bell:
    """k(F) = k0 exp(-F delta / kT)."""

    k0: float
    delta: float = 0.0

    def __call__(self, F, kT):
        return self.k0 * np.exp(-np.asarray(F, dtype=float) * self.delta / kT)


@dataclass(frozen=True)
class PiecewiseBell:
    """Two exponential branches joined continuously at F_ref.

    k(F) = k_ref exp(-(F - F_ref) delta_below / kT)   for F <  F_ref
    k(F) = k_ref exp(-(F - F_ref) delta_above / kT)   for F >= F_ref

    This is the form of Kondo et al. (2023) eq. 5, whose distance d_j has the
    opposite sign to ``delta`` here (their k = lambda exp(d (L - L_j)/kT)).
    """

    k_ref: float
    F_ref: float
    delta_below: float
    delta_above: float

    def __call__(self, F, kT):
        F = np.asarray(F, dtype=float)
        delta = np.where(F < self.F_ref, self.delta_below, self.delta_above)
        return self.k_ref * np.exp(-(F - self.F_ref) * delta / kT)


@dataclass(frozen=True)
class Scaled:
    """A rate law multiplied by a constant factor (e.g. a gate factor)."""

    law: object
    factor: float

    def __call__(self, F, kT):
        return self.factor * self.law(F, kT)


def constant(k: float) -> Bell:
    """A load-independent rate."""
    return Bell(k, 0.0)
