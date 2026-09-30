"""Shared code for Part B (valley, MaxEnt member, predictions near limits).

Everything computed with this module is a PREDICTION: model parameters come only
from the Carter & Cross fit (Part A, v3) and from the papers named below; no
limit dataset (Visscher, Block, Sozanski, Taniguchi's backstep counts, Hong,
Budaitis, ...) enters any parameter.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np

from competing_exits import bet
from competing_exits.diffusion import REFLECT, stokes_einstein_D
from competing_exits.headrace import (HEAD_RADIUS_NM, HeadSearch, IProjLandscape, Landscape, V3_FIXED,
                                      head_race_v3)
from competing_exits.model import ExitKind
from v3fit import KT, LOWER, UPPER, V2_REFIT, lnkf, pack, residuals, unpack

ROOT = Path(__file__).resolve().parents[1]

# ---------------------------------------------------------------------------
# Fixed inputs (with sources)
# ---------------------------------------------------------------------------
KAPPA0_BRACKET = [0.01, 0.03, 0.10, 0.21, 0.30]   # undocked tether stiffness, kT/nm^2 (REPORT2 §B.2)
KAPPA0_NAMED = {"floppy (Guydosh 2009 pulling)": 0.03, "WLC (Kutys 2010: Lp 0.7 nm, 28 res.)": 0.21}
P_MEAN = 0.2            # nm: tethered intermediate 8.4 nm from the rear site (Mickolajczyk 2015) minus 8.2
P_MEAN_RANGE = (-1.8, 2.2)   # +-2 nm label offset (brief)
BUDGET = 1.2            # kT, docking free energy (Rice 2003 via Block 2007, Xu 2021); range 1-2
X = bet.grid(16001)


# ---------------------------------------------------------------------------
# Members of the landscape family
# ---------------------------------------------------------------------------
@dataclass
class Member:
    """A docked landscape plus its refitted v3 parameters."""

    name: str
    landscape: object          # Landscape (Gaussian tether) or IProjLandscape
    theta: np.ndarray          # (B, ln kb0, delta_b, ln kc, ln kon, T_ms)
    chi2: float
    sig: tuple

    @property
    def p(self):
        return unpack(self.theta)

    def logq(self):
        """Log docked density (zero load, tether part only: the binding barrier is a
        property of the site and is common to docked and undocked states)."""
        L = self.landscape
        if isinstance(L, Landscape):
            return bet.tether_logdensity(L.kappa, L.x_eq)
        return bet.iprojection_logdensity(L.kappa0, L.lambdas, L.fs, L.x0)

    def search(self, eta_rel=1.0, kT=None, D_scale=1.0, B=None, start=None):
        kT = KT if kT is None else kT
        L = self.landscape
        from dataclasses import replace
        L = replace(L, B_front=self.p["B"] if B is None else B)
        D = stokes_einstein_D(V3_FIXED["kT"], V3_FIXED["eta"], HEAD_RADIUS_NM) * D_scale / eta_rel
        st = start if start is not None else bet_start()
        return HeadSearch(L, D, start=st, rear=REFLECT, start_tilt=True)


def bet_start(mean=P_MEAN, kappa=V3_FIXED["start_kappa"]):
    from competing_exits.headrace import start_density
    return start_density(mean, kappa)


def gaussian_member(name, kappa, x_eq, theta, chi2, sig):
    return Member(name, Landscape(kappa, x_eq, B_front=unpack(theta)["B"]), np.asarray(theta), chi2, sig)


# ---------------------------------------------------------------------------
# Fitting the I-projection family (exact capture rates)
# ---------------------------------------------------------------------------
def iproj_lnkf(kappa0, lambdas, fs, B, loads):
    return _iproj_lnkf(float(kappa0), tuple(float(l) for l in lambdas), tuple(float(f) for f in fs), float(B),
                       tuple(float(F) for F in loads))


@lru_cache(maxsize=200000)
def _iproj_lnkf(kappa0, lambdas, fs, B, loads):
    L = IProjLandscape(kappa0, lambdas, fs, B_front=B)
    D = stokes_einstein_D(V3_FIXED["kT"], V3_FIXED["eta"], HEAD_RADIUS_NM)
    hs = HeadSearch(L, D, start=bet_start(), rear=REFLECT, start_tilt=True)
    return np.array([np.log(hs.rates(float(F), KT)["front"]) for F in loads])


# ---------------------------------------------------------------------------
# Densities and costs
# ---------------------------------------------------------------------------
P_KAPPA = V3_FIXED["start_kappa"]   # kT/nm^2: width of p = the search's start distribution (Part A)


def densities(member: Member, kappa0: float, F: float = 0.0, p_mean: float = P_MEAN, p_kappa=None, kT=None):
    """p0, q and p on the grid at load F (all tilted by the load).

    p, the head's distribution when docking starts, is the SAME density the front
    search starts from (tracked mean P_MEAN, width 1/P_KAPPA); it does not depend on
    the reference kappa0.  p_kappa = kappa0 is the variant in which p is the undocked
    equilibrium itself (shifted to the tracked mean)."""
    kT = KT if kT is None else kT
    pk = P_KAPPA if p_kappa is None else p_kappa
    p0 = bet.tilted(bet.tether_logdensity(kappa0, 0.0), X, F, kT)
    q = bet.tilted(member.logq(), X, F, kT)
    p = bet.tilted(bet.tether_logdensity(pk, p_mean), X, F, kT)
    return p0, q, p


def costs(member: Member, kappa0: float, F: float = 0.0, **kw):
    p0, q, p = densities(member, kappa0, F, **kw)
    # after undocking (clock: ATP leaves) the head relaxes from q back to the undocked
    # equilibrium p0, dissipating D(q || p0)
    return {"commitment": bet.commitment(q, p0, X), "mismatch": bet.mismatch(p, q, X),
            "reverse": bet.commitment(q, p0, X)}


# ---------------------------------------------------------------------------
# The motor for a member under a condition
# ---------------------------------------------------------------------------
def motor(member: Member, atp_uM=1000.0, eta_rel=1.0, gate_viscous=False, T_cv2=1.0):
    p = member.p
    kb0 = p["kb0"] / (eta_rel if gate_viscous else 1.0)
    return head_race_v3(0.0, 0.0, p["B"], kb0, p["delta_b"], p["kc"], p["kon"], p["T"], atp_uM=atp_uM,
                        eta_rel=eta_rel, T_cv2=T_cv2, search=member.search())


def bookkeeping(member, kappa0, F, atp_uM=1000.0, eta_rel=1.0, gate_viscous=False, T_cv2=1.0, **kw):
    """Kinetics and dissipation accounting at one condition."""
    m = motor(member, atp_uM, eta_rel, gate_viscous, T_cv2)
    ts = m.transport_stats(F, step_size=V3_FIXED["d"])
    P = m.splitting(F)
    c = costs(member, kappa0, F, **kw)
    attempts_per_net = 1.0 / (P["forward"] - P["back"]) if P["forward"] > P["back"] else np.inf
    steps_per_net = (P["forward"] + P["back"]) / (P["forward"] - P["back"]) if P["forward"] > P["back"] else np.inf
    attempt_rate = 1.0 / ts["cycle_time"]
    step_rate = (P["forward"] + P["back"]) * attempt_rate
    D = c["mismatch"]
    return {
        "F": F, "v": ts["v"], "randomness": ts["randomness"], "dwell_mean": ts["dwell_mean"],
        "dwell_cv": ts["dwell_cv"], "odds": P["forward"] / P["back"],
        "P_forward": P["forward"], "P_back": P["back"], "P_clock": P["clock"],
        "attempts_per_step": ts["attempts_per_step"], "attempts_per_net_step": attempts_per_net,
        "steps_per_net_step": steps_per_net,
        "mismatch_per_quench": D,
        # primary accounting: one quench per committed step (= per ATP hydrolysed)
        "mismatch_per_step": D,
        "mismatch_per_net_step": D * steps_per_net,
        "mismatch_per_s": D * step_rate,
        # upper-bound accounting: a quench at every ATP binding (attempt), incl. futile ones
        "mismatch_per_step_all_attempts": D * ts["attempts_per_step"],
        "mismatch_per_net_step_all_attempts": D * attempts_per_net,
        "mismatch_per_s_all_attempts": D * attempt_rate,
        "reverse_per_step_all_attempts": c["reverse"] * P["clock"] * ts["attempts_per_step"],
        "commitment": c["commitment"],
        "attempt_rate": attempt_rate, "step_rate": step_rate,
        # total dissipation per net forward step for Delta mu = 20.5 kT (one ATP per step)
        "total_dissipation_per_net_step_dmu20.5": 20.5 * steps_per_net - F * V3_FIXED["d"] / KT,
        "tur_bound_per_net_step_model": 2.0 / ts["randomness"] if ts["randomness"] > 0 else np.inf,
    }


def save(name, data):
    from common import save_json
    return save_json(name, data)
