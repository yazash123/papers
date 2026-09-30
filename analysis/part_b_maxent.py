"""Part B.1-B.2: the valley of docked landscapes, the minimum-relative-entropy
(MaxEnt) member relative to the undocked tether p0, and the docking-budget cut.

Valley: docked tethers (kappa', x') whose Carter & Cross chi^2, with every other
v3 parameter refitted (B, gate, clock, k_on, T), is within Delta chi^2 <= 5.99
(95%, 2 parameters) of the best fit; 2.30 (68%) is also reported.  Two weightings:
session 1's assumed per-bin scatter (sigma ln dwell 0.10, ln ratio 0.20) and the
scatter of Carter & Cross's own plotted bins about their fits (0.20, 0.33).

MaxEnt member: argmin over the valley of D(q || p0), p0 = N(0, 1/kappa0) on
[-8, 8] (reflecting), for kappa0 across the bracket.  Also fitted: the
I-projection family q ~ p0 exp(-l1 e^{-f1 x} - l2 e^{-f2 x}) (f at 3 and 9 pN),
which is the MaxEnt family when the capture rates constrain <e^{-f x}>_q.

Output: results/part_b_valley.json.  Run: cd analysis && PYTHONPATH=.. python part_b_maxent.py
"""
from __future__ import annotations

import json
import time

import numpy as np
from scipy.optimize import least_squares

from competing_exits import bet
from partb_common import (BUDGET, KAPPA0_BRACKET, P_MEAN, ROOT, X, Member, gaussian_member, iproj_lnkf, save)
from v3fit import (DWELL_1MM, DWELL_10UM, KT, LOADS, LOWER, RATIO, UPPER, V2_REFIT, fit_at, lnkf, pack, unpack)

LEVELS = {"68%": 2.30, "95%": 5.99}
F_IPROJ = (3.0, 9.0)   # loads whose tilts define the two I-projection basis functions


def load_part_a():
    return json.loads((ROOT / "results" / "part_a_v3.json").read_text())


def commitment_grid(kappas, xeqs, kappa0):
    p0 = bet.tilted(bet.tether_logdensity(kappa0, 0.0), X)
    out = np.empty((len(kappas), len(xeqs)))
    for i, k in enumerate(kappas):
        for j, x in enumerate(xeqs):
            q = bet.tilted(bet.tether_logdensity(k, x), X)
            out[i, j] = bet.commitment(q, p0, X)
    return out


def inner(kappa, x, sig, starts):
    best = (None, np.inf)
    for th0 in starts:
        try:
            th, c = fit_at(kappa, x, th0, sig_dwell=sig[0], sig_ratio=sig[1])
        except Exception:
            continue
        if c < best[1]:
            best = (th, c)
    return best


def refine_member(kappas, xeqs, chi, thetas, D, chimin, sig, level, kappa0):
    """Minimum of D over {Delta chi^2 <= level}, refined on a fine local grid."""
    mask = np.isfinite(chi) & (chi - chimin <= level)
    if not mask.any():
        return None
    Dm = np.where(mask, D, np.inf)
    i, j = np.unravel_index(np.argmin(Dm), Dm.shape)
    ki = np.log(kappas)
    lo_k, hi_k = ki[max(i - 1, 0)], ki[min(i + 1, len(ki) - 1)]
    lo_x, hi_x = xeqs[max(j - 1, 0)], xeqs[min(j + 1, len(xeqs) - 1)]
    p0 = bet.tilted(bet.tether_logdensity(kappa0, 0.0), X)
    best = None
    for lk in np.linspace(lo_k, hi_k, 9):
        for x in np.linspace(lo_x, hi_x, 9):
            k = float(np.exp(lk))
            th, c = inner(k, x, sig, [thetas[i, j]])
            if th is None or c - chimin > level:
                continue
            d = bet.commitment(bet.tilted(bet.tether_logdensity(k, x), X), p0, X)
            if best is None or d < best["D"]:
                best = {"kappa": k, "x_eq": float(x), "D": d, "chi2": c, "dchi2": c - chimin, "theta": th}
    return best


def iproj_residuals(v, kappa0, sig):
    l1, l2 = v[0], v[1]
    th = v[2:]
    p = unpack(th)
    fs = tuple(0.5 * F / KT for F in F_IPROJ)
    kf = np.exp(iproj_lnkf(kappa0, (l1, l2), fs, p["B"], LOADS))
    kb = p["kb0"] * np.exp(-LOADS * p["delta_b"] / KT)
    out = [np.log(kf / kb / RATIO(LOADS)) / sig[1]]
    for atp, data in ((1000.0, DWELL_1MM), (10.0, DWELL_10UM)):
        K = kf + kb + p["kc"]
        d = (1.0 / (p["kon"] * atp) + 1.0 / K) / (1.0 - p["kc"] / K) + p["T"]
        out.append(np.log(d / data(LOADS)) / sig[0])
    r = np.concatenate(out)
    r[~np.isfinite(r)] = 1e3      # capture impossible (e.g. an extreme wall): heavy penalty
    return r


def fit_iproj(kappa0, sig, lam0=(0.0, 0.0), theta0=None, fixed_lambda=None):
    th0 = pack({**V2_REFIT, "B": 5.0}) if theta0 is None else theta0
    if fixed_lambda is not None:
        f = lambda th: iproj_residuals(np.concatenate([fixed_lambda, th]), kappa0, sig)
        r = least_squares(f, np.clip(th0, LOWER + 1e-9, UPPER - 1e-9), bounds=(LOWER, UPPER), x_scale="jac",
                          xtol=1e-10, ftol=1e-10, gtol=1e-10, max_nfev=3000, diff_step=1e-6)
        return np.concatenate([fixed_lambda, r.x]), float(np.sum(r.fun ** 2))
    lo = np.concatenate([[-50.0, -50.0], LOWER])
    hi = np.concatenate([[200.0, 200.0], UPPER])
    x0 = np.clip(np.concatenate([lam0, th0]), lo + 1e-9, hi - 1e-9)
    r = least_squares(lambda v: iproj_residuals(v, kappa0, sig), x0, bounds=(lo, hi), x_scale="jac",
                      xtol=1e-10, ftol=1e-10, gtol=1e-10, max_nfev=4000, diff_step=1e-6)
    return r.x, float(np.sum(r.fun ** 2))


def iproj_density(kappa0, lambdas):
    fs = tuple(0.5 * F / KT for F in F_IPROJ)
    return bet.tilted(bet.iprojection_logdensity(kappa0, lambdas, fs), X)


def iproj_scan(kappa0, sig, chimin, level=5.99, n=24):
    """Min D(q_lambda || p0) over the I-projection valley.  Along each ray from
    lambda = 0 (p0 itself) the commitment grows, so on each ray the candidate is the
    first point entering the valley: a geometric scan in the ray length, then
    bisection to the boundary.  The member is the smallest D over the rays."""
    p0 = bet.tilted(bet.tether_logdensity(kappa0, 0.0), X)
    th_p0, chi_p0 = fit_iproj(kappa0, sig, fixed_lambda=np.zeros(2))
    rec = {"p0_in_valley": bool(chi_p0 - chimin <= level), "chi2_p0": chi_p0}
    if rec["p0_in_valley"]:
        rec["member"] = {"lambdas": [0.0, 0.0], "D": 0.0, "chi2": chi_p0, "dchi2": chi_p0 - chimin,
                         "theta": th_p0[2:].tolist()}
        return rec
    best = None
    for ang in np.linspace(-np.pi, np.pi, n, endpoint=False):
        u = np.array([np.cos(ang), np.sin(ang)])
        th = th_p0[2:]
        lo_t, lo_th, hit = 0.0, th, None
        for t in np.geomspace(1e-5, 100.0, 36):
            thr, c = fit_iproj(kappa0, sig, fixed_lambda=t * u, theta0=th)
            th = thr[2:]
            if c - chimin <= level:
                hit = (t, th, c)
                break
            lo_t, lo_th = t, th
        if hit is None:
            continue
        hi_t, hi_th, hi_c = hit
        for _ in range(10):  # bisection to the boundary
            mid = 0.5 * (lo_t + hi_t)
            thr, c = fit_iproj(kappa0, sig, fixed_lambda=mid * u, theta0=hi_th)
            if c - chimin <= level:
                hi_t, hi_th, hi_c = mid, thr[2:], c
            else:
                lo_t = mid
        lam = hi_t * u
        d = bet.commitment(iproj_density(kappa0, lam), p0, X)
        if best is None or d < best["D"]:
            best = {"lambdas": lam.tolist(), "D": d, "chi2": hi_c, "dchi2": hi_c - chimin, "theta": hi_th.tolist(),
                    "angle": float(ang)}
    rec["member"] = best
    return rec


def iproj_best(kappa0, sig):
    """Best chi^2 within the I-projection family (lambdas free)."""
    best = (None, np.inf)
    for lam0 in ((0.0, 0.0), (-0.01, 0.005), (0.05, 0.0), (0.0, 0.02), (-0.05, 0.02)):
        try:
            v, c = fit_iproj(kappa0, sig, lam0=np.array(lam0))
        except Exception:
            continue
        if c < best[1]:
            best = (v, c)
    return best


def main():
    t0 = time.time()
    A = load_part_a()
    kappas = np.array(A["grid"]["kappa"])
    xeqs = np.array(A["grid"]["x_eq"])
    out = {"levels": LEVELS, "kappa0_bracket": KAPPA0_BRACKET, "budget_kT": BUDGET, "grid": A["grid"], "weights": {}}
    for signame in ("session1", "data"):
        S = A[signame]
        sig = (S["sigma_ln_dwell"], S["sigma_ln_ratio"])
        chi = np.array(S["chi2_grid"], dtype=float)
        thetas = np.array(S["theta_grid"], dtype=float)
        chimin = min(S["best"]["chi2"], np.nanmin(chi))
        W = {"sigma": sig, "chi2_min": chimin, "best": S["best"], "per_kappa0": {}}
        # valley extent (range of kappa', x_eq with the refitted B and gate along it)
        for lname, lev in LEVELS.items():
            m = np.isfinite(chi) & (chi - chimin <= lev)
            ii, jj = np.where(m)
            W[f"valley_{lname}"] = {"n_points": int(m.sum()),
                                   "kappa_range": [float(kappas[ii].min()), float(kappas[ii].max())] if m.any() else None,
                                   "x_eq_range": [float(xeqs[jj].min()), float(xeqs[jj].max())] if m.any() else None}
        for kappa0 in KAPPA0_BRACKET:
            D = commitment_grid(kappas, xeqs, kappa0)
            rec = {"D_grid": D}
            for lname, lev in LEVELS.items():
                mem = refine_member(kappas, xeqs, chi, thetas, D, chimin, sig, lev, kappa0)
                rec[f"maxent_{lname}"] = mem
                # budget: fraction of the valley affordable
                m = np.isfinite(chi) & (chi - chimin <= lev)
                rec[f"affordable_fraction_{lname}"] = float((D[m] <= BUDGET).mean()) if m.any() else None
                rec[f"D_range_{lname}"] = [float(D[m].min()), float(D[m].max())] if m.any() else None
            # is the undocked tether itself (no bet) in the valley?
            th, c = inner(kappa0, 0.0, sig, [pack({**V2_REFIT, "B": b}) for b in (0.0, 3.0, 8.0)])
            rec["no_bet_chi2"] = c
            rec["no_bet_dchi2"] = c - chimin
            W["per_kappa0"][f"{kappa0:g}"] = rec
            print(signame, "kappa0", kappa0, "maxent95", {k: v for k, v in (rec["maxent_95%"] or {}).items() if k != "theta"},
                  "no-bet dchi2 %.2f" % rec["no_bet_dchi2"], "(%.0fs)" % (time.time() - t0))
        # I-projection family (the MaxEnt family) across the kappa0 bracket
        W["iprojection"] = {}
        for kappa0 in KAPPA0_BRACKET:
            vb, cb = iproj_best(kappa0, sig)
            ref = min(chimin, cb)
            ip = iproj_scan(kappa0, sig, ref)
            ip["family_best_chi2"] = cb
            ip["family_best_lambdas"] = vb[:2].tolist() if vb is not None else None
            ip["reference_chi2"] = ref
            W["iprojection"][f"{kappa0:g}"] = ip
            print(signame, "iproj kappa0", kappa0, ip["member"] and {k: v for k, v in ip["member"].items() if k != "theta"},
                  "chi2(p0) %.2f" % ip["chi2_p0"], "(%.0fs)" % (time.time() - t0))
        out["weights"][signame] = W
    save("part_b_valley.json", out)


if __name__ == "__main__":
    main()
