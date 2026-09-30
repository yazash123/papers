"""Part C: confront the Part B predictions (committed in PREDICTIONS.md, commit
84a3d19, before this script was first run) with the limit datasets, feature by
feature; the thermodynamic uncertainty relation (TUR) bound from measured
randomness; the difference predictions for hidden dissipation.

Output: results/part_c_confront.json, results/tables_partC.md, figures/fig7_confront.png.
Run: cd analysis && PYTHONPATH=.. python part_c_confront.py
"""
from __future__ import annotations

import csv
import json
import sys

import numpy as np

from common import COLORS, FIGURES, MUTED, RESULTS, plot_style, save_json
from partb_common import ROOT

DIG = ROOT / "data" / "digitized"
NAMES = ["iproj_k0=0.03", "iproj_k0=0.21", "maxent_k0=0.03", "maxent_k0=0.21", "best", "middle", "floppy_edge"]
MAXENT = NAMES[:4]
KT = 4.087            # pN nm at 23 C (v3)
D_STEP = 8.2          # nm
DMU = 20.5            # kT per ATP (Takaki's reference line; REPORT2 0.3)
ARIGA_HIDDEN = 16.0   # kT per step at 2 pN, second-hand (REPORT2 0.2), uncertainty unknown
CC_RATIO = (802.0, 0.95)   # Carter & Cross Fig. 2a inset fit: ratio = 802 exp(-0.95 F) (p.310)

# Taniguchi et al. 2005, Table 1 (p.345), text-extracted: T (K), kf0, kb0 (1/s, mean +- s.e.m.), df, db (nm)
TANIGUCHI_T1 = [
    (280.0, 100.0, 11.0, 0.28, 0.04, 2.4, 0.1, 0.0, 0.1),
    (287.0, 209.0, 17.0, 0.61, 0.06, 2.3, 0.1, 0.1, 0.1),
    (298.0, 544.0, 54.0, 1.6, 0.2, 2.6, 0.1, 0.0, 0.1),
    (308.0, 1353.0, 78.0, 3.8, 0.3, 2.8, 0.1, 0.0, 0.1),
]
K_B = 1.380649e-2     # pN nm / K


def rows(fn):
    with open(DIG / fn) as fh:
        return list(csv.DictReader([l for l in fh if not l.startswith("#")]))


def load_predictions():
    return json.loads((ROOT / "results" / "part_b_limits.json").read_text())


def curve(L, name, atp_key, qty):
    rec = L["members"][name]["load_wide"][atp_key]
    F = np.array([r["F"] for r in rec])
    y = np.array([r[qty] for r in rec], dtype=float)
    return F, y


def model_at(L, atp_key, qty, F, names=NAMES):
    """Model values at loads F for each member: array (len(names), len(F))."""
    out = []
    for n in names:
        Fm, y = curve(L, n, atp_key, qty)
        out.append(np.interp(F, Fm, y))
    return np.array(out)


# ---------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------
def data_sets():
    d = {}
    cc = rows("carter2005_fig2c_velocity.csv")
    for atp, key in (("1mM", "1000"), ("10uM", "10")):
        d[f"cc_v_{key}"] = np.array([(float(r["load_pN"]), float(r["velocity_nm_s"])) for r in cc
                                     if r["atp"] == atp and r["marker"] == "bin"])
        d[f"cc_v0_{key}"] = [float(r["velocity_nm_s"]) for r in cc if r["atp"] == atp and r["marker"] != "bin"][0]
    dw = rows("carter2005_fig2b_dwell.csv")
    for atp, key in (("1mM", "1000"), ("10uM", "10")):
        for step in ("forward", "backward"):
            d[f"cc_dwell_{step}_{key}"] = np.array([(float(r["load_pN"]), float(r["mean_dwell_ms"])) for r in dw
                                                   if r["atp"] == atp and r["step"] == step and r["marker"] == "bin"])
    ra = rows("carter2005_fig2a_ratio.csv")
    for atp, key in (("1mM", "1000"), ("10uM", "10")):
        a = np.array([(float(r["load_pN"]), float(r["forward_to_backward_ratio"])) for r in ra if r["atp"] == atp])
        d[f"cc_ratio_{key}"] = np.unique(a, axis=0)       # one duplicated vector point removed
    v4b = rows("visscher1999_fig4b_randomness_vs_load.csv")
    d["vis_r_2000"] = np.array([(float(r["load_pN"]), float(r["r"]), float(r["sem_approx"])) for r in v4b])
    v3a = rows("visscher1999_fig3a_force_velocity.csv")
    d["vis_v_2000"] = np.array([(float(r["load_pN"]), float(r["velocity_nm_s"]), float(r["sem_approx"])) for r in v3a if r["atp"] == "2mM"])
    d["vis_v_5"] = np.array([(float(r["load_pN"]), float(r["velocity_nm_s"]), float(r["sem_approx"])) for r in v3a if r["atp"] == "5uM"])
    v4a = rows("visscher1999_fig4a_randomness_vs_atp_raw.csv")
    d["vis_r_atp_5.69"] = np.array([(float(r["atp_uM"]), float(r["r"]), abs(float(r["err_approx"]))) for r in v4a
                                    if r["marker"] == "open_5.69pN"])
    d["vis_r_atp_low_loads"] = np.array([(float(r["atp_uM"]), float(r["r"]), abs(float(r["err_approx"]))) for r in v4a
                                         if r["marker"] == "filled"])
    b4 = rows("block2003_fig4_velocity_randomness.csv")
    for atp, key in (("1.6mM", "1600"), ("4.2uM", "4.2")):
        d[f"blk_r_{key}"] = np.array([(float(r["load_pN_hindering"]), float(r["value"]), float(r["sem_approx"])) for r in b4
                                      if r["quantity"] == "randomness" and r["atp"] == atp])
        d[f"blk_v_{key}"] = np.array([(float(r["load_pN_hindering"]), float(r["value"]), float(r["sem_approx"])) for r in b4
                                      if r["quantity"] == "velocity_nm_s" and r["atp"] == atp])
    sz = rows("sozanski2015_fig2b_velocity_vs_viscosity.csv")
    d["soz"] = [(r["crowder"], float(r["eta_eff_over_eta0"]), float(r["velocity_um_s"])) for r in sz]
    return d


# ---------------------------------------------------------------------------
# comparisons
# ---------------------------------------------------------------------------
def compare_sem(L, data, atp_key, qty):
    """Per data point (F, y, sem): model range over all members and over MaxEnt
    members, and z = (y - model)/sem for every member."""
    F, y, s = data[:, 0], data[:, 1], data[:, 2]
    M = model_at(L, atp_key, qty, F)
    z = (y[None, :] - M) / s[None, :]
    return {"F": F, "y": y, "sem": s, "model": M, "z": z,
            "chi2_per_member": dict(zip(NAMES, np.sum(z ** 2, axis=1))), "n": len(F)}


def compare_log(L, data, atp_key, qty, sigma_ln):
    """For data without error bars: residual in ln units, and z with the
    per-bin scatter sigma_ln (C&C's own scatter about their fits, Part A)."""
    F, y = data[:, 0], data[:, 1]
    M = model_at(L, atp_key, qty, F)
    with np.errstate(invalid="ignore", divide="ignore"):
        res = np.log(y[None, :]) - np.log(M)
    return {"F": F, "y": y, "model": M, "ln_residual": res, "z": res / sigma_ln}


def taniguchi_one_to_one():
    """1:1 load kf = kb from Taniguchi's Table 1 fits: F = kT ln(kf0/kb0)/(df - db),
    with first-order error propagation."""
    out = []
    for T, kf, skf, kb, skb, df, sdf, db, sdb in TANIGUCHI_T1:
        kT = K_B * T
        lr = np.log(kf / kb)
        slr = np.hypot(skf / kf, skb / kb)
        dd = df - db
        sdd = np.hypot(sdf, sdb)
        F = kT * lr / dd
        sF = F * np.hypot(slr / lr, sdd / dd)
        out.append({"T_K": T, "ln_ratio0": lr, "F_1to1": F, "sem": sF})
    Ts = np.array([o["T_K"] for o in out])
    Fs = np.array([o["F_1to1"] for o in out])
    ss = np.array([o["sem"] for o in out])
    w = 1 / ss ** 2
    A = np.vstack([np.ones_like(Ts), Ts - 294.0]).T
    cov = np.linalg.inv(A.T @ (w[:, None] * A))
    beta = cov @ (A.T @ (w * Fs))
    slope, sslope = beta[1], np.sqrt(cov[1, 1])
    rel = slope / beta[0]            # fractional change per K
    return {"points": out, "fit_intercept_294K": beta[0], "slope_pN_per_K": slope, "slope_sem": sslope,
            "relative_slope_per_K": rel, "relative_slope_sem": sslope / beta[0],
            "entropic_prediction_relative_slope_per_K": 1.0 / 294.0}


def tur(d):
    """Lower bound on the entropy produced per net forward step, 2/r (k_B), from each
    measured randomness point, against the total dissipation per net step
    Delta mu x (steps per net step) - F d / kT, with the steps per net step from
    Carter & Cross's fitted ratio (802 e^-0.95F)."""
    out = {}
    for key in ("vis_r_2000", "blk_r_1600"):
        a = d[key]
        F, r, s = a[:, 0], a[:, 1], a[:, 2]
        R = CC_RATIO[0] * np.exp(-CC_RATIO[1] * np.clip(F, 0, None))
        steps_per_net = (R + 1) / (R - 1)
        total = DMU * steps_per_net - F * D_STEP / KT
        out[key] = {"F": F, "r": r, "sem": s, "bound": 2.0 / r, "bound_lo": 2.0 / (r + 2 * s),
                    "bound_hi": 2.0 / np.maximum(r - 2 * s, 1e-3), "total_per_net_step": total,
                    "fraction_captured": (2.0 / r) / total}
    return out


def sozanski(d, v0=0.80, v0_range=(0.75, 0.84)):
    """Velocity relative to the uncrowded control: "In the absence of crowders, at
    ATP concentration of 1 mM, we observed average kinesin-1 velocities of about
    800 nm/s" (Sozanski 2015, p.218102-2).  The range 0.75-0.84 um/s is my allowance
    for "about"."""
    out = []
    for crowder, eta, v in d["soz"]:
        out.append({"crowder": crowder, "eta": eta, "v_um_s": v, "rel": v / v0,
                    "rel_lo": v / v0_range[1], "rel_hi": v / v0_range[0]})
    return out


def main(write=True):
    L = load_predictions()
    d = data_sets()
    out = {"source_commit_of_predictions": "84a3d19"}
    # (a) load
    out["cc_v_1000"] = compare_log(L, d["cc_v_1000"][d["cc_v_1000"][:, 1] > 0], "1000", "v", 0.31)
    out["cc_v_10"] = compare_log(L, d["cc_v_10"][d["cc_v_10"][:, 1] > 0], "10", "v", 0.31)
    for key in ("cc_v_1000", "cc_v_10"):
        a = d[key]
        out[key + "_all"] = {"F": a[:, 0], "y": a[:, 1], "model": model_at(L, key.split("_")[-1], "v", a[:, 0])}
    out["vis_v_2000"] = compare_sem(L, d["vis_v_2000"], "2000", "v")
    out["vis_v_5"] = compare_sem(L, d["vis_v_5"], "5", "v")
    out["blk_v_1600"] = compare_sem(L, d["blk_v_1600"], "1600", "v")
    out["blk_v_4.2"] = compare_sem(L, d["blk_v_4.2"], "4.2", "v")
    # (b) randomness
    out["vis_r_2000"] = compare_sem(L, d["vis_r_2000"], "2000", "randomness")
    out["blk_r_1600"] = compare_sem(L, d["blk_r_1600"], "1600", "randomness")
    out["blk_r_4.2"] = compare_sem(L, d["blk_r_4.2"], "4.2", "randomness")
    atps = np.array(L["conditions"]["atps"])
    ra = d["vis_r_atp_5.69"]
    M = np.array([np.interp(np.log(ra[:, 0]), np.log(atps), L["members"][n]["randomness_vs_atp"]["5.69"]) for n in NAMES])
    out["vis_r_atp_5.69"] = {"atp": ra[:, 0], "y": ra[:, 1], "sem": ra[:, 2], "model": M, "z": (ra[:, 1][None, :] - M) / ra[:, 2][None, :]}
    # (c) odds: Carter & Cross ratio at both ATP; the 10 uM points are outside the fit
    for key in ("1000", "10"):
        a = d[f"cc_ratio_{key}"]
        out[f"cc_ratio_{key}"] = compare_log(L, a, key, "odds", 0.33)
    # 10 uM vs 1 mM, point by point at matched loads (both vector-extracted at the same bin centres)
    a1, a10 = d["cc_ratio_1000"], d["cc_ratio_10"]
    pairs = []
    for F, y in a10:
        j = np.argmin(np.abs(a1[:, 0] - F))
        if abs(a1[j, 0] - F) < 0.1:
            pairs.append((F, y, a1[j, 1], np.log(y / a1[j, 1])))
    out["cc_ratio_10_over_1000"] = np.array(pairs)
    # (d) dwells: model mean dwell is the same before forward and backward steps
    for key in ("1000", "10"):
        for step in ("forward", "backward"):
            a = d[f"cc_dwell_{step}_{key}"]
            out[f"cc_dwell_{step}_{key}"] = compare_log(L, np.column_stack([a[:, 0], a[:, 1] / 1e3]), key, "dwell_mean", 0.31)
    # (e) viscosity
    out["sozanski"] = sozanski(d)
    etas = np.array(L["conditions"]["etas"])
    out["model_viscosity_rel"] = {n: np.array([r["v"] for r in L["members"][n]["viscosity"]["F=0.0"]]) /
                                  L["members"][n]["viscosity"]["F=0.0"][0]["v"] for n in NAMES}
    out["model_viscosity_rel_gate"] = {n: np.array([r["v"] for r in L["members"][n]["viscosity_gate_viscous"]["F=0.0"]]) /
                                       L["members"][n]["viscosity_gate_viscous"]["F=0.0"][0]["v"] for n in NAMES}
    out["etas"] = etas
    # (f) temperature
    out["taniguchi_1to1"] = taniguchi_one_to_one()
    # TUR
    out["tur"] = tur(d)
    if write:
        save_json("part_c_confront.json", out)
    return out, d, L


# ---------------------------------------------------------------------------
# tables
# ---------------------------------------------------------------------------
def fmt(x, n=2):
    return "–" if x is None or not np.isfinite(x) else f"{x:.{n}f}"


def tables(out, d, L, fh=sys.stdout):
    p = lambda *a: print(*a, file=fh)
    rng = lambda M, j, idx=None: (np.nanmin(M[idx, j] if idx is not None else M[:, j]),
                                  np.nanmax(M[idx, j] if idx is not None else M[:, j]))
    me = np.array([NAMES.index(n) for n in MAXENT])

    def sem_table(title, key, nd=2):
        c = out[key]
        p(f"\n#### {title}\n")
        p("| load (pN) | data ± s.e.m. | model, all members | model, MaxEnt members | z (MaxEnt members) |")
        p("|---|---|---|---|---|")
        for j in range(c["n"]):
            lo, hi = rng(c["model"], j)
            mlo, mhi = rng(c["model"], j, me)
            zs = c["z"][me, j]
            p(f"| {c['F'][j]:.2f} | {fmt(c['y'][j], nd)} ± {fmt(c['sem'][j], nd)} | {fmt(lo, nd)}–{fmt(hi, nd)} | "
              f"{fmt(mlo, nd)}–{fmt(mhi, nd)} | {fmt(zs.min(), 1)} to {fmt(zs.max(), 1)} |")
        chi = c["chi2_per_member"]
        p(f"\nχ² over {c['n']} points: " + ", ".join(f"{n} {chi[n]:.0f}" for n in NAMES))

    def log_table(title, key, scale=1.0, nd=1, unit=""):
        c = out[key]
        p(f"\n#### {title}\n")
        p(f"| load (pN) | data{unit} | model, all members | model, MaxEnt | ln(data/model), MaxEnt |")
        p("|---|---|---|---|---|")
        for j in range(len(c["F"])):
            lo, hi = rng(c["model"], j)
            mlo, mhi = rng(c["model"], j, me)
            res = c["ln_residual"][me, j]
            p(f"| {c['F'][j]:.2f} | {fmt(c['y'][j] * scale, nd)} | {fmt(lo * scale, nd)}–{fmt(hi * scale, nd)} | "
              f"{fmt(mlo * scale, nd)}–{fmt(mhi * scale, nd)} | {fmt(np.nanmin(res), 2)} to {fmt(np.nanmax(res), 2)} |")

    p("### C.1 Load: velocity through stall (Carter & Cross Fig. 2c; no error bars extracted)\n")
    for key, title in (("cc_v_1000_all", "1 mM"), ("cc_v_10_all", "10 µM")):
        c = out[key]
        p(f"\n#### {title} (all bins; the 3–9 pN window is not independent of the fit)\n")
        p("| load (pN) | data (nm/s) | model, all members | model, MaxEnt |")
        p("|---|---|---|---|")
        for j in range(len(c["F"])):
            lo, hi = rng(c["model"], j)
            mlo, mhi = rng(c["model"], j, me)
            p(f"| {c['F'][j]:.2f} | {c['y'][j]:.1f} | {lo:.1f}–{hi:.1f} | {mlo:.1f}–{mhi:.1f} |")
    sem_table("Visscher 1999 Fig. 3a, 2 mM ATP: velocity (nm/s)", "vis_v_2000", 0)
    sem_table("Visscher 1999 Fig. 3a, 5 µM ATP: velocity (nm/s)", "vis_v_5", 1)
    sem_table("Block 2003 Fig. 4A, 1.6 mM ATP: velocity (nm/s)", "blk_v_1600", 0)
    sem_table("Block 2003 Fig. 4A, 4.2 µM ATP: velocity (nm/s)", "blk_v_4.2", 1)
    p("\n### C.2 Randomness\n")
    sem_table("Visscher 1999 Fig. 4b, 2 mM ATP", "vis_r_2000")
    sem_table("Block 2003 Fig. 4C, 1.6 mM ATP", "blk_r_1600")
    sem_table("Block 2003 Fig. 4C, 4.2 µM ATP", "blk_r_4.2")
    c = out["vis_r_atp_5.69"]
    p("\n#### Visscher 1999 Fig. 4a, 5.69 pN\n")
    p("| [ATP] (µM) | data ± err | model, all | model, MaxEnt |")
    p("|---|---|---|---|")
    for j in range(len(c["atp"])):
        p(f"| {c['atp'][j]:.0f} | {c['y'][j]:.2f} ± {c['sem'][j]:.2f} | {np.min(c['model'][:, j]):.2f}–{np.max(c['model'][:, j]):.2f} | "
          f"{np.min(c['model'][me, j]):.2f}–{np.max(c['model'][me, j]):.2f} |")
    p("\n### C.3 Forward:back odds (Carter & Cross Fig. 2a inset)\n")
    log_table("1 mM (in the fit window 3–9 pN)", "cc_ratio_1000", nd=3)
    log_table("10 µM (NOT fitted: a test of P2)", "cc_ratio_10", nd=3)
    pr = out["cc_ratio_10_over_1000"]
    p("\n#### 10 µM over 1 mM at the same load (P2 predicts ln ratio = 0)\n")
    p("| load (pN) | odds 10 µM | odds 1 mM | ln(10 µM / 1 mM) |")
    p("|---|---|---|---|")
    for F, a, b, l in pr:
        p(f"| {F:.2f} | {a:.3f} | {b:.3f} | {l:+.2f} |")
    hi = pr[pr[:, 0] > 6.0]
    p(f"\nAbove 6 pN: mean ln(10 µM/1 mM) = {hi[:, 3].mean():+.2f} over {len(hi)} load bins "
      f"(s.d. {hi[:, 3].std(ddof=1):.2f}; s.e. of the mean {hi[:, 3].std(ddof=1) / np.sqrt(len(hi)):.2f}); "
      f"all bins: {pr[:, 3].mean():+.2f} ± {pr[:, 3].std(ddof=1) / np.sqrt(len(pr)):.2f}")
    p("\n### C.4 Dwells (Carter & Cross Fig. 2b; model: one dwell distribution for both step directions)\n")
    for key in ("1000", "10"):
        for step in ("forward", "backward"):
            log_table(f"{step} steps, {'1 mM' if key == '1000' else '10 µM'} (ms)", f"cc_dwell_{step}_{key}", scale=1e3, nd=0)
    p("\n### C.5 Viscosity (Sozański 2015 Fig. 2b; v₀ = 0.80 µm/s assumed, range 0.75–0.84)\n")
    etas = out["etas"]
    p("| crowder | η/η₀ | v (µm/s) | v/v₀ | model v/v₀ at that η (all members) | MaxEnt members |")
    p("|---|---|---|---|---|---|")
    for s in out["sozanski"]:
        mods = np.array([np.interp(s["eta"], etas, out["model_viscosity_rel"][n]) for n in NAMES])
        p(f"| {s['crowder']} | {s['eta']:.2f} | {s['v_um_s']:.3f} | {s['rel']:.2f} ({s['rel_lo']:.2f}–{s['rel_hi']:.2f}) | "
          f"{mods.min():.2f}–{mods.max():.2f} | {mods[me].min():.2f}–{mods[me].max():.2f} |")
    t = out["taniguchi_1to1"]
    p("\n### C.6 Temperature: Taniguchi 2005 Table 1 → 1:1 load\n")
    p("| T (K) | ln(k_f0/k_b0) | 1:1 load (pN) |")
    p("|---|---|---|")
    for o in t["points"]:
        p(f"| {o['T_K']:.0f} | {o['ln_ratio0']:.2f} | {o['F_1to1']:.2f} ± {o['sem']:.2f} |")
    p(f"\nWeighted linear fit: relative slope {t['relative_slope_per_K'] * 100:+.2f} ± {t['relative_slope_sem'] * 100:.2f} % per K; "
      f"the entropic default predicts {t['entropic_prediction_relative_slope_per_K'] * 100:+.2f} % per K "
      f"(difference {(t['relative_slope_per_K'] - t['entropic_prediction_relative_slope_per_K']) / t['relative_slope_sem']:+.1f} s.e.m.)")
    p("\n### C.7 TUR: entropy per net forward step ≥ 2/r (k_B)\n")
    p("| data | load (pN) | r ± s.e.m. | bound 2/r (k_BT) | 95% range | total per net step, Δμ 20.5 (k_BT) | fraction the bound captures |")
    p("|---|---|---|---|---|---|---|")
    for key, lab in (("vis_r_2000", "Visscher 2 mM"), ("blk_r_1600", "Block 1.6 mM")):
        c = out["tur"][key]
        for j in range(len(c["F"])):
            p(f"| {lab} | {c['F'][j]:.2f} | {c['r'][j]:.2f} ± {c['sem'][j]:.2f} | {c['bound'][j]:.2f} | "
              f"{c['bound_lo'][j]:.2f}–{c['bound_hi'][j]:.2f} | {c['total_per_net_step'][j]:.1f} | {c['fraction_captured'][j]:.2f} |")


# ---------------------------------------------------------------------------
# figure
# ---------------------------------------------------------------------------
def figure(out, d, L):
    plt = plot_style()
    fig, axs = plt.subplots(2, 3, figsize=(12.8, 7.4))
    band = dict(alpha=0.28, lw=0)
    me = np.array([NAMES.index(n) for n in MAXENT])

    def model_band(ax, atp_key, qty, color, Fmin, Fmax, scale=1.0, label=None):
        F, _ = curve(L, NAMES[0], atp_key, qty)
        sel = (F >= Fmin) & (F <= Fmax)
        M = np.array([curve(L, n, atp_key, qty)[1] for n in NAMES])[:, sel] * scale
        ax.fill_between(F[sel], M.min(0), M.max(0), color=color, **band)
        ax.plot(F[sel], M[me].mean(0), color=color, lw=1.4, label=label)

    # (a) velocity: C&C 1 mM and 10 uM, Visscher 2 mM
    ax = axs[0, 0]
    model_band(ax, "1000", "v", COLORS[0], -15, 15, label="model 1 mM")
    model_band(ax, "10", "v", COLORS[1], -15, 15, label="model 10 µM")
    a = d["cc_v_1000"]
    ax.plot(a[:, 0], a[:, 1], "o", ms=4, color=COLORS[0], mec="white", mew=0.6, label="C&C 1 mM")
    a = d["cc_v_10"]
    ax.plot(a[:, 0], a[:, 1], "s", ms=4, color=COLORS[1], mec="white", mew=0.6, label="C&C 10 µM")
    ax.axhline(0, color=MUTED, lw=0.8)
    ax.axvspan(3, 9, color="#f2f2f2", zorder=0)
    ax.set(xlabel="hindering load F (pN)", ylabel="velocity (nm/s)", ylim=(-60, 900))
    ax.set_title("a  force–velocity through stall (grey: fit window)", fontsize=8.5)
    ax.legend(fontsize=6.8, loc="upper right")
    # inset-like second axis for superstall is avoided; show superstall in panel b
    # (b) superstall: velocity and dwell
    ax = axs[0, 1]
    model_band(ax, "1000", "dwell_mean", COLORS[0], 0, 15, scale=1e3, label="model 1 mM")
    model_band(ax, "10", "dwell_mean", COLORS[1], 0, 15, scale=1e3, label="model 10 µM")
    for key, col, mk in (("1000", COLORS[0], "o"), ("10", COLORS[1], "s")):
        a = d[f"cc_dwell_backward_{key}"]
        ax.plot(a[:, 0], a[:, 1], mk, ms=4.5, color=col, mfc="white", mew=1.1)
        a = d[f"cc_dwell_forward_{key}"]
        a = a[a[:, 0] > 0]
        ax.plot(a[:, 0], a[:, 1], mk, ms=4, color=col, mec="white", mew=0.6)
    ax.set(xlabel="hindering load F (pN)", ylabel="mean dwell (ms)", yscale="log", ylim=(8, 2e4))
    ax.set_title("b  dwells (filled before forward, open before back steps)", fontsize=8.5)
    ax.legend(fontsize=6.8, loc="upper left")
    # (c) odds
    ax = axs[0, 2]
    model_band(ax, "1000", "odds", COLORS[2], 0, 11, label="model (1 mM = 10 µM)")
    for key, col, mk, lab in (("1000", COLORS[0], "o", "C&C 1 mM (fitted)"), ("10", COLORS[1], "s", "C&C 10 µM (not fitted)")):
        a = d[f"cc_ratio_{key}"]
        ax.plot(a[:, 0], a[:, 1], mk, ms=4.5, color=col, mec="white", mew=0.6, label=lab)
    ax.axhline(1, color=MUTED, lw=0.8)
    ax.set(xlabel="hindering load F (pN)", ylabel="forward : back steps", yscale="log", ylim=(0.02, 500))
    ax.set_title("c  step odds: [ATP] moves the 1:1 load", fontsize=8.5)
    ax.legend(fontsize=6.8, loc="upper right")
    # (d) randomness vs load
    ax = axs[1, 0]
    model_band(ax, "2000", "randomness", COLORS[3], 0, 6.6, label="model 2 mM")
    model_band(ax, "1600", "randomness", COLORS[4], -8, 5.0, label="model 1.6 mM")
    a = d["vis_r_2000"]
    ax.errorbar(a[:, 0], a[:, 1], a[:, 2], fmt="o", ms=4, color=COLORS[3], mec="white", mew=0.6, lw=1, label="Visscher 2 mM")
    a = d["blk_r_1600"]
    ax.errorbar(a[:, 0], a[:, 1], a[:, 2], fmt="s", ms=4, color=COLORS[4], mec="white", mew=0.6, lw=1, label="Block 1.6 mM")
    ax.set(xlabel="hindering load F (pN)", ylabel="randomness r", ylim=(0, 2.2))
    ax.set_title("d  randomness against load", fontsize=8.5)
    ax.legend(fontsize=6.8, loc="upper left")
    # (e) viscosity
    ax = axs[1, 1]
    etas = out["etas"]
    Mv = np.array([out["model_viscosity_rel"][n] for n in NAMES])
    ax.fill_between(etas, Mv.min(0), Mv.max(0), color=COLORS[0], **band)
    ax.plot(etas, Mv[me].mean(0), color=COLORS[0], lw=1.4, label="model (all members: band)")
    marks = {"BSA": "o", "Dextran 10 kg/mol": "s", "Dextran 500 kg/mol": "D", "PEG 1000 kg/mol": "^",
             "PEG 18 kg/mol": "v", "PEG 6 kg/mol": "P", "Sucrose": "X", "TetraEG": "*"}
    for cr, mk in marks.items():
        pts = [s for s in out["sozanski"] if s["crowder"] == cr]
        ax.plot([s["eta"] for s in pts], [s["rel"] for s in pts], mk, ms=5, color=MUTED, mfc="white", mew=1.0, label=cr)
    ax.set(xlabel="η_eff/η₀", ylabel="v / v₀ (unloaded, 1 mM)", xlim=(0.9, 10), ylim=(-0.05, 1.75), xscale="log")
    ax.set_title("e  viscosity: Sozański 2015 (v₀ = 0.80 µm/s, their text)", fontsize=8.5)
    ax.legend(fontsize=5.8, loc="upper right", ncol=3, columnspacing=0.8, handletextpad=0.3)
    # (f) TUR and the imprint on one per-net-step scale (kT)
    ax = axs[1, 2]
    for key, col, mk, lab in (("vis_r_2000", COLORS[3], "o", "TUR bound 2/r, Visscher"), ("blk_r_1600", COLORS[4], "s", "TUR bound 2/r, Block")):
        c = out["tur"][key]
        ax.errorbar(c["F"], c["bound"], [c["bound"] - c["bound_lo"], c["bound_hi"] - c["bound"]], fmt=mk, ms=4,
                    color=col, mec="white", mew=0.6, lw=1, label=lab)
    Fg = np.linspace(-5, 6.6, 300)
    R = CC_RATIO[0] * np.exp(-CC_RATIO[1] * np.clip(Fg, 0, None))
    ax.plot(Fg, DMU * (R + 1) / (R - 1) - Fg * D_STEP / KT, color=INK, lw=1.4, label="total per net step (Δμ 20.5, C&C odds)")
    ax.plot([2.0], [ARIGA_HIDDEN], "D", ms=6, color=COLORS[7], label="Ariga hidden (~16, 2 pN)")
    F1, _ = curve(L, NAMES[0], "1000", "mismatch_per_net_step")
    sel = (F1 <= 6.6) & (F1 >= -5)
    Mm = np.array([curve(L, n, "1000", "mismatch_per_net_step")[1] for n in NAMES])[:, sel]
    ax.fill_between(F1[sel], Mm[me].min(0), Mm[me].max(0), color=COLORS[2], **band)
    ax.plot(F1[sel], Mm[me].max(0), color=COLORS[2], lw=1.2, label="predicted imprint (MaxEnt members)")
    ax.set(xlabel="hindering load F (pN)", ylabel="k_BT per net forward step", yscale="log", ylim=(3e-3, 60), xlim=(-5.5, 7))
    ax.set_title("f  TUR bound, total, hidden, and the imprint", fontsize=8.5)
    ax.legend(fontsize=6.3, loc="lower left")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig7_confront.png")
    plt.close(fig)


INK = "#0b0b0b"

if __name__ == "__main__":
    out, d, L = main()
    with open(RESULTS / "tables_partC.md", "w") as fh:
        tables(out, d, L, fh)
    figure(out, d, L)
