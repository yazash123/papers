"""Part A (session 2): head race v3 fitted to Carter & Cross, the chi^2 landscape
over the docked tether (kappa, x_eq), and checks against Kondo et al. 2023.

Output: results/part_a_v3.json (used by Part B), figure fig5 (panel a/b).
Run: cd analysis && PYTHONPATH=.. python part_a_v3.py      (~5 min)
"""
from __future__ import annotations

import csv
import time

import numpy as np

from common import save_json
from competing_exits.kinesin import head_race_v2, kondo_kif5a, kondo_values
from v3fit import (KT, LOADS, RATIO, DWELL_1MM, DWELL_10UM, V2_REFIT, curves, fit_at, guess_B,
                   lnkf, motor, pack, unpack, residuals)

KAPPAS = np.geomspace(0.02, 0.6, 26)
XEQS = np.linspace(-4.0, 8.0, 31)
SIGMAS = {"session1": (0.10, 0.20), "data": (0.20, 0.33)}   # (sigma ln dwell, sigma ln ratio)


def digitised_scatter():
    """RMS log deviation of Carter & Cross's plotted bins from their own fits
    (data/digitized, 1 mM, loads 2.5-9.5 pN)."""
    import pathlib
    root = pathlib.Path(__file__).resolve().parents[1] / "data" / "digitized"

    def rows(fn):
        with open(root / fn) as fh:
            return list(csv.DictReader([l for l in fh if not l.startswith("#")]))
    dw = [(float(r["load_pN"]), float(r["mean_dwell_ms"])) for r in rows("carter2005_fig2b_dwell.csv")
          if r["atp"] == "1mM" and r["step"] == "forward" and r["marker"] == "bin" and 2.5 <= float(r["load_pN"]) <= 9.6]
    ra = [(float(r["load_pN"]), float(r["forward_to_backward_ratio"])) for r in rows("carter2005_fig2a_ratio.csv")
          if r["atp"] == "1mM"]
    ddw = np.array([np.log(d * 1e-3 / DWELL_1MM(F)) for F, d in dw])
    dra = np.array([np.log(q / RATIO(F)) for F, q in ra])
    return {"dwell_1mM_rms": float(np.sqrt(np.mean(ddw ** 2))), "dwell_1mM_rms_without_8.4pN": float(np.sqrt(np.mean(np.delete(ddw, np.argmin(ddw)) ** 2))),
            "ratio_1mM_rms": float(np.sqrt(np.mean(dra ** 2))), "n_dwell": len(ddw), "n_ratio": len(dra)}


def best_inner(kappa, x_eq, sig, starts):
    best = (None, np.inf)
    for th0 in starts:
        try:
            th, chi = fit_at(kappa, x_eq, th0, sig_dwell=sig[0], sig_ratio=sig[1])
        except Exception:
            continue
        if chi < best[1]:
            best = (th, chi)
    return best


def scan(sig):
    """Profile chi^2 over (kappa, x_eq); inner fit over (B, gate, clock, kon, T)."""
    chi = np.full((len(KAPPAS), len(XEQS)), np.nan)
    thetas = np.full((len(KAPPAS), len(XEQS), 6), np.nan)
    base = pack({**V2_REFIT, "B": 5.0})
    for i, k in enumerate(KAPPAS):
        prev = None
        for j, x in enumerate(XEQS):
            starts = [pack({**V2_REFIT, "B": guess_B(k, x)}), base]
            if prev is not None:
                starts.append(prev)
            if i > 0 and np.isfinite(thetas[i - 1, j, 0]):
                starts.append(thetas[i - 1, j])
            th, c = best_inner(k, x, sig, starts)
            if th is not None:
                chi[i, j], thetas[i, j] = c, th
                prev = th
    return chi, thetas


def refine_global(chi, thetas, sig):
    """Polish the global minimum with (kappa, x_eq) free."""
    from scipy.optimize import minimize
    i, j = np.unravel_index(np.nanargmin(chi), chi.shape)
    th0 = thetas[i, j]

    def obj(v):
        k, x = np.exp(v[0]), v[1]
        th, c = best_inner(k, x, sig, [th0])
        return c
    r = minimize(obj, [np.log(KAPPAS[i]), XEQS[j]], method="Nelder-Mead", options=dict(xatol=1e-4, fatol=1e-6, maxiter=400))
    k, x = np.exp(r.x[0]), r.x[1]
    th, c = best_inner(k, x, sig, [th0])
    return k, x, th, c


def summary(k, x, th):
    p = unpack(th)
    m1, m10 = motor(k, x, p, 1000.0), motor(k, x, p, 10.0)
    out = {"kappa": k, "x_eq": x, **{kk: float(v) for kk, v in p.items()},
           "kf0": float(np.exp(lnkf(k, x, p["B"], loads=(0.0,))[0])),
           "stall_1to1_pN": None, "v0_1mM": m1.velocity(0.0), "v0_10uM": m10.velocity(0.0)}
    from scipy.optimize import brentq
    try:
        out["stall_1to1_pN"] = brentq(lambda F: np.log(m1.rate("forward", F) / m1.rate("back", F)), 0.0, 20.0)
    except ValueError:
        pass
    out["v_1mM"] = {f"{F:.0f}": m1.velocity(F) for F in (0, 2, 4, 6, 8)}
    return out


def kondo_check(best):
    """Compare the v3 structure with what Kondo et al. measured for KIF5A (not a fit)."""
    k, x = best["kappa"], best["x_eq"]
    p = {kk: best[kk] for kk in ("B", "kb0", "delta_b", "kc", "kon", "T")}
    m = motor(k, x, p, 1000.0)
    F = np.array([0.0, 2.0, 3.19, 4.0, 6.0, 8.0])
    kf = np.array([m.rate("forward", f) for f in F])
    kb = np.array([m.rate("back", f) for f in F])
    kv = kondo_values()
    kondo = kondo_kif5a().home_state
    kf_k = np.array([kondo.rate("forward", f) for f in F])
    kb_k = np.array([kondo.rate("back", f) for f in F])
    # local load distances d ln k / dF * (-kT), from 0-3 and 4-8 pN
    dist = lambda r, a, b: -KT * np.log(r[b] / r[a]) / (F[b] - F[a])
    ts = [m.transport_stats(f) for f in (3.0, 6.0, 7.5)]
    return {
        "loads": F, "v3_kf": kf, "v3_kb": kb, "kondo_kf": kf_k, "kondo_kb": kb_k,
        "v3_forward_distance_0_to_2pN": dist(kf, 0, 1), "v3_forward_distance_4_to_8pN": dist(kf, 3, 5),
        "kondo_forward_distance_below_knee": -kv["d_f_minus"], "kondo_forward_distance_above_knee": -kv["d_f_plus"],
        "v3_gate_rise_0_to_6pN": kb[4] / kb[0], "kondo_gate_rise_0_to_6pN": kb_k[4] / kb_k[0],
        "v3_backward_over_forward_dwell": 1.0,   # single waiting state: identical dwell distributions
        "kondo_backward_over_forward_dwell": [1.22, 0.23],
        "v3_dwell_cv_1mM": {"3 pN": ts[0]["dwell_cv"], "6 pN": ts[1]["dwell_cv"], "7.5 pN": ts[2]["dwell_cv"]},
    }


def main():
    t0 = time.time()
    out = {"data_driven_scatter": digitised_scatter(), "grid": {"kappa": KAPPAS, "x_eq": XEQS}}
    for name, sig in SIGMAS.items():
        chi, th = scan(sig)
        k, x, thb, cb = refine_global(chi, th, sig)
        out[name] = {"sigma_ln_dwell": sig[0], "sigma_ln_ratio": sig[1], "chi2_grid": chi,
                     "theta_grid": th, "best": {**summary(k, x, thb), "chi2": cb}}
        # v2 refit chi2 under the same weights, for comparison (7 parameters)
        print(name, "best:", out[name]["best"], "chi2 %.2f" % cb, "(%.0fs)" % (time.time() - t0))
    out["kondo_check"] = kondo_check(out["session1"]["best"])
    save_json("part_a_v3.json", out)
    print("kondo:", {k: v for k, v in out["kondo_check"].items() if not isinstance(v, np.ndarray)})


if __name__ == "__main__":
    main()
