"""Part A: head race v2 against Carter & Cross (2005), refit and identifiability.

Data used (Carter & Cross 2005, Nature 435:308, Fig. 2 legend, p.310):
  forward:back ratio  = 802 exp(-0.95 F)             (Fig. 2a inset)
  forward dwell 1 mM  = 0.0036 exp(0.57 F) s          (Fig. 2b, fits above 3 pN)
  forward dwell 10 uM = 0.0256 exp(0.55 F) s          (Fig. 2b, fits above 3 pN)
The paper gives these fitted curves, not tables of binned means, so the refit
targets the curves sampled at 1-pN spacing over 3-9 pN (C&C bin their data at
1-pN intervals, p.310).  Weights are an ASSUMPTION: sigma(ln dwell) = 0.10 per
bin (typical s.e.m. of a bin mean in Fig. 2b) and sigma(ln ratio) = 0.20 (the
ratio rests on few backsteps at low load).  Intervals below scale with these.

Unloaded velocities read off Fig. 2c (trap-off bead velocities, squares):
~829 nm/s at 1 mM, ~142 nm/s at 10 uM -- READ OFF A FIGURE (pixel positions).
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares

from common import COLORS, MUTED, FIGURES, plot_style, save_json

from competing_exits.kinesin import V2, head_race_v2

LOADS = np.arange(3.0, 9.01, 1.0)
RATIO = lambda F: 802.0 * np.exp(-0.95 * F)
DWELL_1MM = lambda F: 3.6e-3 * np.exp(0.57 * F)
DWELL_10UM = lambda F: 25.6e-3 * np.exp(0.55 * F)
SIG_DWELL, SIG_RATIO = 0.10, 0.20
UNLOADED_READ = {"1mM": 829.0, "10uM": 142.0}  # nm/s, read off Fig. 2c

NAMES = ["ln_kf0", "delta_f", "ln_kb0", "delta_b", "ln_kc", "ln_kon", "T_ms"]


def to_params(theta):
    lkf, df, lkb, db, lkc, lkon, Tms = theta
    return dict(kf0=np.exp(lkf), delta_f=df, kb0=np.exp(lkb), delta_b=db,
                kc=np.exp(lkc), kon=np.exp(lkon), T=Tms * 1e-3)


def from_params(p):
    return np.array([np.log(p["kf0"]), p["delta_f"], np.log(p["kb0"]), p["delta_b"],
                     np.log(p["kc"]), np.log(p["kon"]), p["T"] * 1e3])


def model_curves(p, F=LOADS):
    s1, s2 = head_race_v2(1000.0, **p), head_race_v2(10.0, **p)
    ratio = s1.rate("forward", F) / s1.rate("back", F)
    return ratio, np.array([s1.mean_dwell(f) for f in F]), np.array([s2.mean_dwell(f) for f in F])


def residuals(theta):
    ratio, d1, d2 = model_curves(to_params(theta))
    return np.concatenate([
        np.log(ratio / RATIO(LOADS)) / SIG_RATIO,
        np.log(d1 / DWELL_1MM(LOADS)) / SIG_DWELL,
        np.log(d2 / DWELL_10UM(LOADS)) / SIG_DWELL,
    ])


def chi2(theta):
    return float(np.sum(residuals(theta) ** 2))


def fit(theta0, fixed=None):
    """Least squares with parameter j fixed at a value if fixed=(j, value)."""
    free = [i for i in range(len(theta0)) if not fixed or i != fixed[0]]

    def full(x):
        th = np.array(theta0, dtype=float)
        th[free] = x
        if fixed:
            th[fixed[0]] = fixed[1]
        return th

    lower = np.array([-np.inf, 0.0, -np.inf, -5.0, np.log(0.1), -np.inf, 0.0])[free]
    upper = np.array([np.inf, 10.0, np.inf, 5.0, np.log(1e5), np.inf, 200.0])[free]
    x0 = np.clip(np.asarray(theta0, dtype=float)[free], lower + 1e-9, upper - 1e-9)
    r = least_squares(lambda x: residuals(full(x)), x0, bounds=(lower, upper),
                      x_scale="jac", xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=20000)
    th = full(r.x)
    return th, chi2(th)


def profile(j, grid, theta_best):
    out, start = [], theta_best.copy()
    for val in grid:
        th, c = fit(start, fixed=(j, val))
        out.append((val, c, th))
        start = th
    return out


def interval(prof, chi_min, level):
    ok = [v for v, c, _ in prof if c - chi_min <= level]
    return (min(ok), max(ok)) if ok else (np.nan, np.nan)


def main():
    v2 = {k: V2[k] for k in ("kf0", "delta_f", "kb0", "delta_b", "kc", "kon", "T")}
    s = head_race_v2(1000.0)

    # ---- 1. how well v2 reproduces the three relations over 3-9 pN ---------
    ratio, d1, d2 = model_curves(v2)
    dev = {
        "ratio": np.log(ratio / RATIO(LOADS)),
        "dwell_1mM": np.log(d1 / DWELL_1MM(LOADS)),
        "dwell_10uM": np.log(d2 / DWELL_10UM(LOADS)),
    }
    theta_v2 = from_params(v2)

    # ---- 2. refit -----------------------------------------------------------
    best, chi_best = theta_v2, chi2(theta_v2)
    rng = np.random.default_rng(0)
    for trial in range(40):  # multistart around v2
        th0 = theta_v2 + rng.normal(0, [0.5, 0.5, 0.5, 0.3, 1.0, 0.3, 5.0])
        th0[4] = np.log(rng.uniform(5, 2000))
        th, c = fit(th0)
        if c < chi_best - 1e-9:
            best, chi_best = th, c
    best, chi_best = fit(best)
    p_best = to_params(best)

    # ---- 3. profiles ----------------------------------------------------------
    grids = {
        "ln_kc": np.log(np.geomspace(3, 3000, 41)),
        "ln_kf0": np.linspace(best[0] - 2.5, best[0] + 2.5, 41),
        "delta_f": np.linspace(max(0.5, best[1] - 2.5), best[1] + 2.5, 41),
        "ln_kb0": np.linspace(best[2] - 2.0, best[2] + 2.0, 41),
        "delta_b": np.linspace(best[3] - 2.0, best[3] + 2.0, 41),
        "ln_kon": np.linspace(best[5] - 1.0, best[5] + 1.0, 41),
        "T_ms": np.linspace(0.0, 40.0, 41),
    }
    profiles, intervals = {}, {}
    for name, grid in grids.items():
        j = NAMES.index(name)
        # profile both ways from the optimum so the path stays continuous
        lo = [g for g in grid if g <= best[j]][::-1]
        hi = [g for g in grid if g > best[j]]
        prof = sorted(profile(j, lo, best) + profile(j, hi, best), key=lambda t: t[0])
        profiles[name] = prof
        i1, i4 = interval(prof, chi_best, 1.0), interval(prof, chi_best, 4.0)
        conv = (lambda v: np.exp(v)) if name.startswith("ln_") else (lambda v: v)
        intervals[name] = {"best": conv(best[j]), "68%": [conv(i1[0]), conv(i1[1])],
                           "95%": [conv(i4[0]), conv(i4[1])],
                           "at_grid_edge": bool(i4[0] <= grid[0] + 1e-12 or i4[1] >= grid[-1] - 1e-12)}

    # derived combinations along the kc profile (what stays fixed as kc slides)
    kc_track = []
    for val, c, th in profiles["ln_kc"]:
        p = to_params(th)
        sp = head_race_v2(1000.0, **p)
        Fs = sp.balance_point("forward", "back")
        kc_track.append({
            "kc": np.exp(val), "dchi2": c - chi_best, "kf0": p["kf0"], "delta_f": p["delta_f"],
            "kb0": p["kb0"], "delta_b": p["delta_b"], "kon": p["kon"], "T_ms": p["T"] * 1e3,
            "lnO0": np.log(p["kf0"] / p["kb0"]), "lever": p["delta_f"] - p["delta_b"],
            "stall": Fs, "price_at_stall": sp.price("forward", Fs),
            "v0_1mM": sp.velocity(0.0), "v0_10uM": head_race_v2(10.0, **p).velocity(0.0),
        })

    # ---- 4. unloaded speed ----------------------------------------------------
    t_unl_1mM = 8.2 / UNLOADED_READ["1mM"]
    t_unl_10uM = 8.2 / UNLOADED_READ["10uM"]
    unloaded = {
        "v2_v0_1mM": s.velocity(0.0), "v2_v0_10uM": head_race_v2(10.0).velocity(0.0),
        "measured_trap_off_1mM": UNLOADED_READ["1mM"], "measured_trap_off_10uM": UNLOADED_READ["10uM"],
        "implied_dwell_1mM_ms": t_unl_1mM * 1e3, "implied_dwell_10uM_ms": t_unl_10uM * 1e3,
        "v2_T_ms": V2["T"] * 1e3,
        "v2_atp_wait_10uM_ms": 1e3 / (V2["kon"] * 10.0),
        # the unloaded 10-uM dwell must contain the ATP wait: kon >= 1/(10 uM * dwell)
        "min_kon_for_10uM_unloaded": 1.0 / (10.0 * t_unl_10uM),
        # if the rest of the cycle equals the unloaded 1-mM dwell minus the race and ATP wait
        "max_T_consistent_with_1mM_ms": (t_unl_1mM - 1 / (V2["kon"] * 1000) - 1 / s.total_rate(0.0)) * 1e3,
        # offset that would turn the true unloaded dwell into v2's 18.0 ms
        "offset_needed_ms": (s.mean_dwell(0.0) - t_unl_1mM) * 1e3,
        "loaded_bins_read_1mM": {"1.6 pN": 372.0, "2.6 pN": 341.0, "3.6 pN": 210.0, "4.5 pN": 113.0},
        "v2_at_bins": {f"{F} pN": s.velocity(F) for F in (1.6, 2.6, 3.6, 4.5)},
    }

    # the width of "fits" is set by the assumed per-bin scatter: chi2 scales as
    # (0.10/sigma)^2, so recompute the clock interval for other sigmas
    kc_vals = np.array([t["kc"] for t in kc_track])
    dchi = np.array([t["dchi2"] for t in kc_track])
    kc_by_sigma = {}
    for sig in (0.10, 0.15, 0.20, 0.30):
        ok = kc_vals[dchi * (SIG_DWELL / sig) ** 2 <= 4.0]
        kc_by_sigma[f"{sig:.2f}"] = [float(ok.min()), float(ok.max())]
    joint = {"v2_dchi2": chi2(theta_v2) - chi_best,
             "chi2_7dof_95pct": 14.07}

    out = {
        "v2_parameters": v2,
        "kc_95pct_by_assumed_sigma": kc_by_sigma,
        "v2_joint": joint,
        "weights": {"sigma_ln_dwell": SIG_DWELL, "sigma_ln_ratio": SIG_RATIO, "loads": LOADS},
        "v2_log_deviation_3_to_9pN": dev,
        "v2_rms_log_deviation": {k: float(np.sqrt(np.mean(v ** 2))) for k, v in dev.items()},
        "v2_chi2": chi2(theta_v2),
        "refit_best": {**p_best, "T_ms": p_best["T"] * 1e3, "chi2": chi_best, "n_points": 3 * len(LOADS)},
        "profile_intervals": intervals,
        "kc_profile": kc_track,
        "ratio_1to1_load_from_CC_fit": np.log(802.0) / 0.95,
        "lever_from_CC_slope_nm": 0.95 * V2["kT"],
        "unloaded": unloaded,
    }
    save_json("part_a_v2.json", out)

    # ---- figure ---------------------------------------------------------------
    plt = plot_style()
    fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.3))
    F = np.linspace(0, 10, 201)
    r, a1, a2 = model_curves(v2, F)
    ax[0].semilogy(F, RATIO(F), color=MUTED, lw=2, ls="--", label="C&C fit 802e$^{-0.95F}$")
    ax[0].semilogy(F, r, color=COLORS[0], lw=3, label="v2 $k_f/k_b$ (coincides)")
    ax[0].semilogy(F, RATIO(F), color="#0b0b0b", lw=1, ls="--")
    ax[0].axhline(1, color=MUTED, lw=0.8)
    ax[0].set(xlabel="hindering load F (pN)", ylabel="forward : back", title="Odds")
    ax[0].legend(loc="lower left")
    for curve, model, c, lab in [(DWELL_1MM, a1, COLORS[0], "1 mM"), (DWELL_10UM, a2, COLORS[1], "10 µM")]:
        sel = F >= 3
        ax[1].semilogy(F[sel], 1e3 * curve(F[sel]), color=c, lw=2, ls="--")
        ax[1].semilogy(F, 1e3 * model, color=c, lw=2)
        ax[1].text(0.2, 1e3 * model[0] * 1.35, lab, color="#0b0b0b", fontsize=8.5, va="bottom")
    ax[1].plot([], [], color=MUTED, ls="--", label="C&C fits (3–9 pN)")
    ax[1].plot([], [], color=MUTED, label="v2")
    ax[1].legend(loc="upper left")
    ax[1].set(xlabel="hindering load F (pN)", ylabel="mean dwell (ms)", title="Dwell per step")
    kcs = np.array([t["kc"] for t in kc_track])
    dch = np.array([t["dchi2"] for t in kc_track])
    ax[2].semilogx(kcs, dch, color=COLORS[0], lw=2, marker="o", ms=3, label="σ = 10% per bin")
    ax[2].semilogx(kcs, dch / 4, color=COLORS[1], lw=2, marker="s", ms=3, label="σ = 20% per bin")
    ax[2].axvline(V2["kc"], color=MUTED, lw=0.8, ls="--")
    ax[2].text(V2["kc"] * 1.08, 0.4, "v2", fontsize=8, color=MUTED)
    ax[2].legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, fontsize=8)
    ax[2].axhline(1, color=MUTED, lw=0.8, ls=":")
    ax[2].axhline(4, color=MUTED, lw=0.8, ls=":")
    ax[2].text(3.5, 1.15, "Δχ² = 1", fontsize=8, color=MUTED)
    ax[2].text(3.5, 4.15, "Δχ² = 4", fontsize=8, color=MUTED)
    ax[2].set(xlabel="clock $k_c$ (1/s), others refitted", ylabel="Δχ² (clock profile)",
              ylim=(-0.3, 12))
    fig.tight_layout()
    fig.savefig(FIGURES / "fig1_v2_vs_carter_cross.png")

    # console summary
    print("v2 rms log deviation:", out["v2_rms_log_deviation"], "chi2", round(out["v2_chi2"], 2))
    print("refit:", {k: round(float(v), 4) for k, v in out["refit_best"].items()})
    for k, v in intervals.items():
        print(f"  {k:8s} best {v['best']:.4g}  68% {v['68%'][0]:.4g}–{v['68%'][1]:.4g}"
              f"  95% {v['95%'][0]:.4g}–{v['95%'][1]:.4g}  edge={v['at_grid_edge']}")
    print("kc 95% by sigma:", kc_by_sigma, "v2 joint dchi2", joint)
    print("unloaded:", {k: v for k, v in unloaded.items() if not isinstance(v, dict)})


if __name__ == "__main__":
    main()
