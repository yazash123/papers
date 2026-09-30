"""Head race v3 fitted to Carter & Cross (2005): shared code for Parts A and B.

Data (Carter & Cross 2005, Nature 435:308, Fig. 2 legend, p.310; verified in
session 1):
    forward:back ratio  = 802 exp(-0.95 F)
    forward dwell 1 mM  = 3.6 exp(0.57 F) ms      (fits above 3 pN)
    forward dwell 10 uM = 25.6 exp(0.55 F) ms
sampled at 1-pN spacing over 3-9 pN, as in session 1 (analysis/part_a_v2.py).

v3 parameters
    kappa, x_eq, B   docked landscape (kT/nm^2, nm, kT); B >= 0 (a binding barrier)
    kb0, delta_b     backstep gate kb0 exp(-F delta_b / kT)
    kc               clock (ATP leaves; restart)
    kon              ATP binding (1/(uM s))
    T                rest of the cycle after a step (s)
Fixed: geometry of v1 (sites +-8 nm, barrier at +6 nm, sd 0.7 nm), D from
Stokes-Einstein at 23 C (water 0.932 mPa s, r = 2.5 nm), load F/2 on the head,
start distribution p = N(0.2 nm, 1/kappa_p) (tracking; kappa_p = 0.21 kT/nm^2),
d = 8.2 nm, kT = 4.087 pN nm.
"""
from __future__ import annotations

from functools import lru_cache

import numpy as np
from scipy.optimize import least_squares

from competing_exits.headrace import V3_FIXED, head_race_v3, v3_search

KT = V3_FIXED["kT"]
LOADS = np.arange(3.0, 9.01, 1.0)
RATIO = lambda F: 802.0 * np.exp(-0.95 * F)
DWELL_1MM = lambda F: 3.6e-3 * np.exp(0.57 * F)
DWELL_10UM = lambda F: 25.6e-3 * np.exp(0.55 * F)
SIG_DWELL, SIG_RATIO = 0.10, 0.20   # session 1's assumed per-bin scatter (see REPORT2 for data-driven values)

# v2 refit (session 1) as the starting point for the non-landscape parameters
V2_REFIT = dict(kb0=8.1, delta_b=0.69, kc=53.0, kon=0.98, T=15.7e-3)


@lru_cache(maxsize=400000)
def lnkf_cached(kappa, x_eq, B, loads=tuple(LOADS), start_mean=None, start_kappa=None, eta=None, N=3200):
    hs = v3_search(kappa, x_eq, B, start_mean=start_mean, start_kappa=start_kappa, eta=eta, N=N)
    return np.array([np.log(hs.rates(float(F), KT)["front"]) for F in loads])


def lnkf(kappa, x_eq, B, loads=LOADS, **kw):
    return lnkf_cached(float(kappa), float(x_eq), float(B), tuple(float(f) for f in loads), **kw)


def unpack(theta):
    """theta = (B, ln kb0, delta_b, ln kc, ln kon, T_ms)."""
    B, lkb, db, lkc, lkon, Tms = theta
    return dict(B=B, kb0=np.exp(lkb), delta_b=db, kc=np.exp(lkc), kon=np.exp(lkon), T=Tms * 1e-3)


def pack(p):
    return np.array([p["B"], np.log(p["kb0"]), p["delta_b"], np.log(p["kc"]), np.log(p["kon"]), p["T"] * 1e3])


def curves(kappa, x_eq, p, loads=LOADS):
    """Model ratio and mean forward dwells at 1 mM and 10 uM (as in session 1: the mean
    dwell per step, WaitingState.mean_dwell)."""
    kf = np.exp(lnkf(kappa, x_eq, p["B"], loads))
    kb = p["kb0"] * np.exp(-loads * p["delta_b"] / KT)
    kc = p["kc"]
    out = {"ratio": kf / kb}
    for atp, key in ((1000.0, "d1"), (10.0, "d2")):
        K = kf + kb + kc
        attempt = 1.0 / (p["kon"] * atp) + 1.0 / K
        Pr = kc / K
        out[key] = attempt / (1.0 - Pr) + p["T"]
    return out


def residuals(kappa, x_eq, theta, sig_dwell=SIG_DWELL, sig_ratio=SIG_RATIO, loads=LOADS):
    p = unpack(theta)
    c = curves(kappa, x_eq, p, loads)
    return np.concatenate([
        np.log(c["ratio"] / RATIO(loads)) / sig_ratio,
        np.log(c["d1"] / DWELL_1MM(loads)) / sig_dwell,
        np.log(c["d2"] / DWELL_10UM(loads)) / sig_dwell,
    ])


LOWER = np.array([0.0, np.log(1e-3), -5.0, np.log(0.1), np.log(1e-3), 0.0])
UPPER = np.array([40.0, np.log(1e4), 5.0, np.log(1e5), np.log(1e3), 200.0])


def fit_at(kappa, x_eq, theta0=None, **kw):
    """Best (B, gate, clock, kon, T) for a fixed docked tether (kappa, x_eq)."""
    if theta0 is None:
        theta0 = pack({**V2_REFIT, "B": guess_B(kappa, x_eq)})
    x0 = np.clip(theta0, LOWER + 1e-9, UPPER - 1e-9)
    r = least_squares(lambda th: residuals(kappa, x_eq, th, **kw), x0, bounds=(LOWER, UPPER),
                      x_scale="jac", xtol=1e-10, ftol=1e-10, gtol=1e-10, max_nfev=4000, diff_step=1e-6)
    return r.x, float(np.sum(r.fun ** 2))


def guess_B(kappa, x_eq, target_kf0=5000.0):
    """B such that the unloaded capture rate is near the v2 refit's k_f0 (for a start)."""
    from scipy.optimize import brentq
    f = lambda B: lnkf(kappa, x_eq, B, loads=(0.0,))[0] - np.log(target_kf0)
    lo, hi = 0.0, 40.0
    try:
        if f(lo) < 0:
            return 0.0
        return brentq(f, lo, hi, xtol=1e-4)
    except ValueError:
        return 5.0


def motor(kappa, x_eq, p, atp_uM=1000.0, **kw):
    return head_race_v3(kappa, x_eq, p["B"], p["kb0"], p["delta_b"], p["kc"], p["kon"], p["T"], atp_uM=atp_uM, **kw)
