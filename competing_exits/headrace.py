"""Diffusive head races: the tethered head searches for a binding site in a 1D
landscape, and the capture rates feed a competing-exits ``WaitingState``.

Geometry (head race v1, as specified for session 2):
    x (nm) = position of the free head relative to the bound head; rear site at
    -8 nm, front site at +8 nm (absorbing when a site is a capture exit);
    overdamped, D = kT / (6 pi eta r) with r = 2.5 nm (87.3 nm^2/us at 1 mPa s);
    a tether (1/2) kappa (x - x_eq)^2; binding barriers B exp(-(x -+ 6)^2/(2 w^2))
    at +-6 nm with w = 0.7 nm (the Gaussian standard deviation; see REPORT2.md for
    why this shape); a hindering load F on the stalk acts on the head as F/2, i.e.
    adds (F/2) x / kT to U.

A capture rate is the exponential-race rate of the diffusive search,
k_site(F) = P(site first) / E[search time], computed without the other clocks
(the competing-exits approximation).  ``exact_race`` gives the full non-
exponential race with the clocks, to check that approximation.

Units: nm, s, pN, kT in pN nm, U in kT.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace
from functools import lru_cache
from typing import Optional

import numpy as np

from .diffusion import ABSORB, REFLECT, Search1D, gaussian_density, stokes_einstein_D
from .model import Exit, ExitKind, WaitingState
from .rates import Bell, constant

P_, W_, R_ = ExitKind.PRODUCTIVE, ExitKind.WRONG, ExitKind.RESTART

HEAD_RADIUS_NM = 2.5
SITE = 8.0
BARRIER_X = 6.0
BARRIER_W = 0.7


# ---------------------------------------------------------------------------
# Landscapes
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Landscape:
    """U(x)/kT = (1/2) kappa (x - x_eq)^2 + (F/2) x / kT + barriers.

    B_front, B_rear: barrier heights (kT) at +x_b and -x_b; None = no barrier.
    """

    kappa: float
    x_eq: float
    B_front: Optional[float] = None
    B_rear: Optional[float] = None
    x_b: float = BARRIER_X
    w: float = BARRIER_W

    def U(self, F: float, kT: float):
        f = 0.5 * F / kT
        k, xe, xb, w = self.kappa, self.x_eq, self.x_b, self.w
        Bf, Br = self.B_front, self.B_rear

        def u(x):
            x = np.asarray(x, dtype=float)
            out = 0.5 * k * (x - xe) ** 2 + f * x
            if Bf:
                out = out + Bf * np.exp(-(x - xb) ** 2 / (2 * w * w))
            if Br:
                out = out + Br * np.exp(-(x + xb) ** 2 / (2 * w * w))
            return out

        return u

    def tether_only(self) -> "Landscape":
        return replace(self, B_front=None, B_rear=None)

    def barrier_heights(self, F: float, kT: float, x_grid=None):
        """(front, rear) barrier heights in kT: U(+-x_b) - min U (the definition that
        reproduces the v1 targets; see REPORT2.md), and the true maxima."""
        x = np.linspace(-SITE, SITE, 32001) if x_grid is None else x_grid
        u = self.U(F, kT)(x)
        umin = u.min()
        at = (float(self.U(F, kT)(self.x_b) - umin), float(self.U(F, kT)(-self.x_b) - umin))
        front = u[(x > 0.5 * self.x_b) & (x < SITE)].max() - umin
        rear = u[(x < -0.5 * self.x_b) & (x > -SITE)].max() - umin
        return {"at_xb": at, "max": (float(front), float(rear))}


# ---------------------------------------------------------------------------
# The search and its capture rates
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class HeadSearch:
    """A diffusive search in ``landscape`` from ``start`` (a position or a density
    callable), with the rear site absorbing (v1) or reflecting (v3)."""

    landscape: Landscape
    D: float                      # nm^2/s at the reference conditions
    start: object = -4.0
    rear: str = ABSORB
    N: int = 3200

    def search(self, F: float, kT: float, D: Optional[float] = None) -> Search1D:
        return Search1D(self.landscape.U(F, kT), self.D if D is None else D,
                        a=-SITE, b=SITE, left=self.rear, right=ABSORB, N=self.N)

    def rates(self, F: float, kT: float, D: Optional[float] = None) -> dict:
        """Capture rates {'front': k, 'rear': k} (1/s) of the exponential race."""
        return _rates_cached(self, float(F), float(kT), None if D is None else float(D))

    def exact_race(self, F: float, kT: float, clocks: dict, D: Optional[float] = None):
        """The full race: search + Poisson clocks, no exponential approximation."""
        return self.search(F, kT, D).race(self.start, clocks)

    def mean_capture_time(self, F: float, kT: float) -> float:
        return self.search(F, kT).capture(self.start)["mean_time"]


@lru_cache(maxsize=200000)
def _rates_cached(hs: HeadSearch, F: float, kT: float, D: Optional[float]):
    c = hs.search(F, kT, D).capture(hs.start)
    T = c["mean_time"]
    return {"front": c["right"] / T, "rear": c.get("left", 0.0) / T}


@dataclass(frozen=True)
class DiffusiveCapture:
    """Rate law for an exit: capture at 'front' or 'rear' by a HeadSearch.

    ``eta_rel`` scales the viscosity (D -> D / eta_rel); ``D_scale`` multiplies D
    (used for temperature)."""

    search: HeadSearch
    site: str = "front"
    eta_rel: float = 1.0
    D_scale: float = 1.0

    def __call__(self, F, kT):
        F = np.asarray(F, dtype=float)
        D = self.search.D * self.D_scale / self.eta_rel
        if F.ndim == 0:
            return self.search.rates(float(F), kT, D)[self.site]
        return np.array([self.search.rates(float(f), kT, D)[self.site] for f in F.ravel()]).reshape(F.shape)


# ---------------------------------------------------------------------------
# Head race v1 (regression model)
# ---------------------------------------------------------------------------
V1 = dict(
    kT=4.114,          # pN nm
    eta=1.0,           # mPa s
    kappa=0.30,        # kT/nm^2, docked tether
    x_eq=2.5,          # nm
    B=9.144,           # kT, both binding barriers
    x0=-4.0,           # nm, search start
    kc=100.0,          # 1/s, clock (ATP leaves)
    t_atp=0.5e-3,      # s, ATP wait
    T=8e-3,            # s, rest of cycle after a forward capture
    d=8.0,             # nm, step (reproduces the v1 speed targets; see REPORT2.md)
)


def v1_search(**over) -> HeadSearch:
    p = {**V1, **over}
    land = Landscape(p["kappa"], p["x_eq"], B_front=p["B"], B_rear=p["B"])
    D = stokes_einstein_D(p["kT"], p["eta"], HEAD_RADIUS_NM)
    return HeadSearch(land, D, start=p["x0"], rear=ABSORB)


def head_race_v1(**over) -> WaitingState:
    """v1: after ATP binds, the head searches both sites; front capture is a step
    (+d, then the rest of the cycle T); rear capture and the clock are futile
    restarts with no displacement and no rest time."""
    p = {**V1, **over}
    hs = v1_search(**over)
    return WaitingState("atp_bound", kT=p["kT"], entry_time=p["t_atp"], exits=(
        Exit("forward", P_, DiffusiveCapture(hs, "front"), step=+p["d"], cost_time=p["T"], fuel=1),
        Exit("rear", R_, DiffusiveCapture(hs, "rear"), fuel=1),
        Exit("clock", R_, constant(p["kc"])),
    ))


# ---------------------------------------------------------------------------
# Head race v3: docked landscape for the front search, a gate for the rear exit
# ---------------------------------------------------------------------------
V3_FIXED = dict(
    kT=4.087,          # pN nm, 23 C (Carter & Cross)
    eta=0.932,         # mPa s, water at 23 C
    d=8.2,             # nm, step
    start_mean=0.2,    # nm, tracked tethered head (Mickolajczyk 2015: 8.4 - 8.2 nm)
    start_kappa=0.21,  # kT/nm^2, width of p (the undocked tether; bracketed in Part B)
)


def start_density(mean: float, kappa: float):
    return gaussian_density(mean, 1.0 / kappa)


def v3_search(kappa, x_eq, B, start_mean=None, start_kappa=None, kT=None, eta=None, N=3200) -> HeadSearch:
    f = V3_FIXED
    kT = f["kT"] if kT is None else kT
    eta = f["eta"] if eta is None else eta
    sm = f["start_mean"] if start_mean is None else start_mean
    sk = f["start_kappa"] if start_kappa is None else start_kappa
    land = Landscape(kappa, x_eq, B_front=B, B_rear=None)
    D = stokes_einstein_D(kT, eta, HEAD_RADIUS_NM)
    return HeadSearch(land, D, start=start_density(sm, sk), rear=REFLECT, N=N)


def head_race_v3(kappa, x_eq, B, kb0, delta_b, kc, kon, T, atp_uM=1000.0,
                 eta_rel=1.0, T_cv2=1.0, search: Optional[HeadSearch] = None, **fixed) -> WaitingState:
    """v3: ATP binding (entry wait 1/(kon [ATP])) docks the neck linker; the head then
    searches for the FRONT site in the docked landscape (rear reflecting: no diffusive
    rear capture).  It competes with a load-insensitive backstep gate
    kb0 exp(-F delta_b / kT) and with the clock kc (ATP leaves: restart).  Both steps
    are followed by the rest of the cycle T."""
    f = {**V3_FIXED, **fixed}
    hs = search or v3_search(kappa, x_eq, B, kT=f["kT"], eta=f["eta"])
    d = f["d"]
    return WaitingState("atp_bound", kT=f["kT"], entry_time=1.0 / (kon * atp_uM), exits=(
        Exit("forward", P_, DiffusiveCapture(hs, "front", eta_rel=eta_rel), step=+d, cost_time=T, fuel=1, cost_cv2=T_cv2),
        Exit("back", W_, Bell(kb0, delta_b), step=-d, cost_time=T, fuel=1, cost_cv2=T_cv2),
        Exit("clock", R_, constant(kc)),
    ))
