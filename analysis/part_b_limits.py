"""Part B.3-B.4: predictions near the motor's limits for the MaxEnt members and
for representative other members of the valley.  PREDICTIONS ONLY: nothing here
reads a limit dataset.

Members (session-1 weights, 95% valley, refitted exactly):
    maxent_k0=0.03   MaxEnt member for the floppy undocked tether (Guydosh 2009)
    maxent_k0=0.21   MaxEnt member for the worm-like-chain undocked tether (Kutys 2010)
    best             best fit (stiffest part of the valley, binding barrier B = 0)
    floppy_edge      the floppiest docked tether in the valley
    middle           a member half-way along the valley
(a) load 0-12 pN at 1 mM ATP, and randomness against [ATP] at 1.05, 3.59, 5.69 pN;
(b) viscosity eta/eta0 = 1-10: the search slows as 1/eta; the gate is
    viscosity-independent (default) or slows as 1/eta (variant); clock, ATP binding
    and rest of cycle unchanged;
(c) temperature 5-40 C: see TEMPERATURE below.
Rest-of-cycle shape: two equal exponential stages (T_cv2 = 1/2) by default
(Mickolajczyk 2015: one rate-limiting transition in each of the 1HB and 2HB
states); exponential (T_cv2 = 1) as a variant.

Output: results/part_b_limits.json.  Run: cd analysis && PYTHONPATH=.. python part_b_limits.py
"""
from __future__ import annotations

import json
import time
from dataclasses import replace

import numpy as np
from scipy.optimize import brentq

from competing_exits import bet
from competing_exits.diffusion import REFLECT, stokes_einstein_D
from competing_exits.headrace import HEAD_RADIUS_NM, HeadSearch, IProjLandscape, V3_FIXED, head_race_v3
from competing_exits.rates import K_B
from partb_common import (KAPPA0_NAMED, P_MEAN, ROOT, X, Member, bookkeeping, bet_start, costs, gaussian_member,
                          motor, save)
from v3fit import KT, V2_REFIT, fit_at, pack, unpack

LOADS = np.round(np.arange(0.0, 12.01, 0.25), 2)
# the conditions (not the values) of the limit datasets: Carter & Cross 2005 (1 mM, 10 uM,
# -15..15 pN), Visscher 1999 (2 mM, 5 uM), Block 2003 (1.6 mM, 4.2 uM; assisting loads)
LOADS_WIDE = np.round(np.arange(-15.0, 15.01, 0.25), 2)
ATP_CONDITIONS = (1000.0, 10.0, 2000.0, 5.0, 1600.0, 4.2)
ETAS = np.array([1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0])
TEMPS_C = np.array([5.0, 10.0, 15.0, 20.0, 23.0, 25.0, 30.0, 35.0, 40.0])
ATPS = np.geomspace(1.0, 5000.0, 25)
T_CV2 = 0.5
SIG = (0.10, 0.20)
T0 = 296.15  # K, Carter & Cross 23 C


# ---------------------------------------------------------------------------
# Temperature model (default = Taniguchi's enthalpies; see PREDICTIONS.md)
# ---------------------------------------------------------------------------
def eta_water(T):
    """Water viscosity (mPa s), Vogel form 0.02414 * 10^(247.8/(T - 140)) (standard
    fit; 1.002 at 20 C, 0.890 at 25 C)."""
    return 0.02414 * 10 ** (247.8 / (T - 140.0))


def dlneta_dT(T):
    return -247.8 * np.log(10) / (T - 140.0) ** 2


# NOTE (found in Part C, after PREDICTIONS.md was committed): Taniguchi's Delta H
# excludes their diffusion-limited prefactor A_T ~ T/eta(T) (Methods, p.346); the raw
# Arrhenius enthalpy of their k_f0 (Table 1) is 26.7 +- 1.1 k_B T0.  The default below
# treats 18.3 as the TOTAL enthalpy of k_f(0); it is kept so that PREDICTIONS.md stays
# reproducible.  It does not affect the 1:1-load prediction (enthalpy difference 0.1 in
# both readings) or anything at 23 C.  See REPORT2.md, Part C.6.
H_FORWARD = 18.3   # k_B T0 (Taniguchi 2005 Table 2): enthalpy of the forward capture rate at F = 0
H_GATE = 18.2      # k_B T0 (Taniguchi 2005 Table 2): backstep
H_OTHER = 18.3     # k_B T0: clock, ATP binding, rest of cycle (ASSUMPTION; variants below)
H_D = 1.0 - T0 * dlneta_dT(T0)     # enthalpy of D = kT/(6 pi eta r), in k_B T0 (~7.9)
H_B = H_FORWARD - H_D              # what the binding barrier must carry


def arrhenius(H, T):
    """exp(-(H/k_B T0) (T0/T - 1)): rate factor at T relative to T0 for enthalpy H (in k_B T0)."""
    return np.exp(-H * (T0 / T - 1.0))


def member_at_temperature(member: Member, Tc, H_rest=H_OTHER, tether="entropic"):
    """Parameters at temperature Tc (C).

    Default (entropic landscape): the docked tether, the binding barrier and the
    start distribution keep their shape in units of kT; D follows Stokes-Einstein
    with water's viscosity; the rest of the forward capture enthalpy, H_B =
    H_FORWARD - H_D, is a load-independent prefactor (applied as a factor on D, which
    multiplies the capture rate exactly).  So k_f(F = 0) has Taniguchi's 18.3 k_B T0.
    Gate: Taniguchi's 18.2 k_B T0.  Clock, ATP binding: H_OTHER (assumption).  Rest of
    cycle: H_rest.  Variant tether='enthalpic': every tether is a spring of fixed
    stiffness in pN/nm, so kappa', kappa_p (start and p) and kappa0 in kT/nm^2 all
    scale as T0/T; x_eq fixed (Gaussian members only)."""
    T = Tc + 273.15
    kT = K_B * T
    p = member.p
    L = member.landscape
    s = T0 / T if tether == "enthalpic" else 1.0
    if s != 1.0 and hasattr(L, "kappa"):
        L = replace(L, kappa=L.kappa * s)
    L = replace(L, B_front=p["B"])
    D = stokes_einstein_D(kT, eta_water(T), HEAD_RADIUS_NM) * arrhenius(H_B, T)
    start = bet_start(kappa=V3_FIXED["start_kappa"] * s)
    hs = HeadSearch(L, D, start=start, rear=REFLECT, start_tilt=True)
    return dict(search=hs, kT=kT, kb0=p["kb0"] * arrhenius(H_GATE, T), delta_b=p["delta_b"],
                kc=p["kc"] * arrhenius(H_OTHER, T), kon=p["kon"] * arrhenius(H_OTHER, T),
                T_rest=p["T"] / arrhenius(H_rest, T), member=replace(member, landscape=L), scale=s)


def temperature_point(member, Tc, F, kappa0, H_rest=H_OTHER, tether="entropic"):
    m = member_at_temperature(member, Tc, H_rest, tether)
    s = head_race_v3(0, 0, 0, m["kb0"], m["delta_b"], m["kc"], m["kon"], m["T_rest"], search=m["search"],
                     T_cv2=T_CV2, kT=m["kT"])
    ts = s.transport_stats(F, step_size=V3_FIXED["d"])
    P = s.splitting(F)
    try:
        stall = brentq(lambda f: np.log(s.rate("forward", f) / s.rate("back", f)), 0.0, 25.0)
    except ValueError:
        stall = np.nan
    c = costs(m["member"], kappa0 * m["scale"], F, kT=m["kT"], p_kappa=V3_FIXED["start_kappa"] * m["scale"])
    D = c["mismatch"]
    step_rate = (P["forward"] + P["back"]) / ts["cycle_time"]
    return {"T_C": Tc, "F": F, "v": ts["v"], "randomness": ts["randomness"], "odds": P["forward"] / P["back"],
            "stall_1to1": stall, "dwell_cv": ts["dwell_cv"], "dwell_mean": ts["dwell_mean"],
            # same accounting as partb_common.bookkeeping: one quench per committed step
            "mismatch_per_step": D, "mismatch_per_s": D * step_rate,
            "mismatch_per_step_all_attempts": D * ts["attempts_per_step"],
            "mismatch_per_s_all_attempts": D / ts["cycle_time"],
            "commitment": c["commitment"], "kf0": s.rate("forward", 0.0)}


# ---------------------------------------------------------------------------
def pick_members():
    V = json.loads((ROOT / "results" / "part_b_valley.json").read_text())
    A = json.loads((ROOT / "results" / "part_a_v3.json").read_text())
    W = V["weights"]["session1"]
    kappas = np.array(A["grid"]["kappa"])
    xeqs = np.array(A["grid"]["x_eq"])
    chi = np.array(A["session1"]["chi2_grid"], dtype=float)
    th = np.array(A["session1"]["theta_grid"], dtype=float)
    chimin = W["chi2_min"]
    members = {}
    for label, k0 in (("maxent_k0=0.03", 0.03), ("maxent_k0=0.21", 0.21)):
        m = W["per_kappa0"][f"{k0:g}"]["maxent_95%"]
        members[label] = gaussian_member(label, m["kappa"], m["x_eq"], np.array(m["theta"]), m["chi2"], SIG)
    b = A["session1"]["best"]
    members["best"] = gaussian_member("best", b["kappa"], b["x_eq"], pack({k: b[k] for k in ("B", "kb0", "delta_b", "kc", "kon", "T")}), b["chi2"], SIG)
    valley = np.isfinite(chi) & (chi - chimin <= 5.99)
    ii, jj = np.where(valley)
    # floppy edge: smallest kappa in the valley (best x_eq there)
    i_min = ii.min()
    j = jj[ii == i_min][np.argmin(chi[i_min, jj[ii == i_min]])]
    members["floppy_edge"] = gaussian_member("floppy_edge", kappas[i_min], xeqs[j], th[i_min, j], chi[i_min, j], SIG)
    # middle: geometric middle kappa of the valley, best x_eq there
    i_mid = int(round((ii.min() + ii.max()) / 2))
    cols = jj[ii == i_mid]
    j = cols[np.argmin(chi[i_mid, cols])]
    members["middle"] = gaussian_member("middle", kappas[i_mid], xeqs[j], th[i_mid, j], chi[i_mid, j], SIG)
    # re-polish every member's inner parameters exactly
    for k, m in members.items():
        L = m.landscape
        thb, c = fit_at(L.kappa, L.x_eq, m.theta, sig_dwell=SIG[0], sig_ratio=SIG[1])
        members[k] = gaussian_member(k, L.kappa, L.x_eq, thb, c, SIG)
    # the I-projection (MaxEnt-family) members
    from part_b_maxent import F_IPROJ, fit_iproj
    fs = tuple(0.5 * F / KT for F in F_IPROJ)
    for k0 in (0.03, 0.21):
        ip = W["iprojection"].get(f"{k0:g}", {}).get("member")
        if not ip:
            continue
        lam = np.array(ip["lambdas"], dtype=float)
        v, c = fit_iproj(k0, SIG, fixed_lambda=lam, theta0=np.array(ip["theta"], dtype=float))
        th = v[2:]
        L = IProjLandscape(k0, tuple(float(l) for l in lam), fs, B_front=unpack(th)["B"])
        members[f"iproj_k0={k0}"] = Member(f"iproj_k0={k0}", L, th, c, SIG)
    return members, chimin


def describe(m: Member, chimin):
    p = m.p
    L = m.landscape
    if hasattr(L, "kappa"):
        shape = {"family": "gaussian", "kappa": L.kappa, "x_eq": L.x_eq}
    else:
        shape = {"family": "iprojection", "kappa0": L.kappa0, "lambdas": list(L.lambdas), "fs": list(L.fs)}
    out = {**shape, "chi2": m.chi2, "dchi2": m.chi2 - chimin, **{k: float(v) for k, v in p.items()}}
    for name, k0 in KAPPA0_NAMED.items():
        c = costs(m, k0, 0.0)
        out[f"commitment_k0={k0}"] = c["commitment"]
        out[f"mismatch_k0={k0}"] = c["mismatch"]
        # p as the undocked equilibrium (width kappa0) instead of the start distribution
        out[f"mismatch_p=undocked_k0={k0}"] = costs(m, k0, 0.0, p_kappa=k0)["mismatch"]
    # sensitivity of the imprint to p (q held fixed): tracked mean +-2 nm (label offset),
    # and a wider p (SD 3.0 nm, the spread of Mickolajczyk's substep sizes)
    out["mismatch_p_mean-1.8"] = costs(m, 0.21, 0.0, p_mean=P_MEAN - 2.0)["mismatch"]
    out["mismatch_p_mean+2.2"] = costs(m, 0.21, 0.0, p_mean=P_MEAN + 2.0)["mismatch"]
    out["mismatch_p_sd3.0"] = costs(m, 0.21, 0.0, p_kappa=1.0 / 9.0)["mismatch"]
    q = bet.tilted(m.logq(), X)
    mq, vq = bet.moments(q, X)
    out["q_mean_sd"] = (mq, float(np.sqrt(vq)))
    return out


def main():
    t0 = time.time()
    members, chimin = pick_members()
    out = {"conditions": {"loads": LOADS, "loads_wide": LOADS_WIDE, "atp_conditions": ATP_CONDITIONS, "etas": ETAS, "temps_C": TEMPS_C, "atps": ATPS, "T_cv2": T_CV2,
                          "sigma": SIG, "H": {"forward": H_FORWARD, "gate": H_GATE, "other": H_OTHER, "D": H_D, "B": H_B}},
           "members": {}}
    for name, m in members.items():
        rec = {"params": describe(m, chimin)}
        # (a) load
        for k0name, k0 in KAPPA0_NAMED.items():
            rec[f"load_k0={k0}"] = [bookkeeping(m, k0, float(F), T_cv2=T_CV2) for F in LOADS]
        rec["load_exponential_rest"] = [bookkeeping(m, 0.21, float(F), T_cv2=1.0) for F in LOADS]
        rec["load_10uM"] = [bookkeeping(m, 0.21, float(F), atp_uM=10.0, T_cv2=T_CV2) for F in LOADS]
        rec["load_wide"] = {f"{a:g}": [bookkeeping(m, 0.21, float(F), atp_uM=a, T_cv2=T_CV2) for F in LOADS_WIDE]
                            for a in ATP_CONDITIONS}
        rec["randomness_vs_atp"] = {f"{F}": [motor(m, float(a), T_cv2=T_CV2).transport_stats(F, step_size=8.2)["randomness"] for a in ATPS]
                                    for F in (1.05, 3.59, 5.69)}
        # (b) viscosity
        for gv in (False, True):
            key = "viscosity_gate_viscous" if gv else "viscosity"
            rec[key] = {f"F={F}": [bookkeeping(m, 0.21, F, eta_rel=float(e), gate_viscous=gv, T_cv2=T_CV2) for e in ETAS]
                        for F in (0.0, 3.0)}
        # (c) temperature
        rec["temperature"] = {}
        for variant, kw in (("default", {}), ("rest_H10", {"H_rest": 10.0}), ("rest_H26", {"H_rest": 26.2}),
                            ("enthalpic_tether", {"tether": "enthalpic"})):
            rec["temperature"][variant] = {f"F={F}": [temperature_point(m, float(Tc), F, 0.21, **kw) for Tc in TEMPS_C]
                                           for F in (0.0, 5.0)}
        out["members"][name] = rec
        print(name, rec["params"], "(%.0fs)" % (time.time() - t0))
    save("part_b_limits.json", out)


if __name__ == "__main__":
    main()
