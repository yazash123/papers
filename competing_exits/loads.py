"""Loading protocols: the hindering load felt at motor position x (nm)."""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Trap:
    """Stationary optical trap: F = stiffness * (x - x_trap).  pN/nm, nm."""

    stiffness: float
    x_trap: float = 0.0

    def __call__(self, x):
        return self.stiffness * (x - self.x_trap)


@dataclass(frozen=True)
class Clamp:
    """Force clamp: the same load F whatever the position."""

    F: float

    def __call__(self, x):
        return self.F + 0.0 * x
