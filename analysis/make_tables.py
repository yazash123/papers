"""Print the markdown tables of PREDICTIONS.md (and REPORT2.md) straight from
results/*.json, so that every number in the text is traceable to a script output.

Run: cd analysis && PYTHONPATH=.. python make_tables.py > ../results/tables_partB.md
"""
from __future__ import annotations

import json
import sys

import numpy as np

from partb_common import ROOT

NAMES = ["iproj_k0=0.03", "iproj_k0=0.21", "maxent_k0=0.03", "maxent_k0=0.21", "best", "middle", "floppy_edge"]
LABEL = {"iproj_k0=0.03": "I-projection, κ₀ 0.03", "iproj_k0=0.21": "I-projection, κ₀ 0.21",
         "maxent_k0=0.03": "Gaussian MaxEnt, κ₀ 0.03", "maxent_k0=0.21": "Gaussian MaxEnt, κ₀ 0.21",
         "best": "best fit", "middle": "mid-valley", "floppy_edge": "floppy edge"}


def f(x, n=2):
    if x is None:
        return "–"
    if isinstance(x, (list, tuple)):
        return "(" + ", ".join(f(v, n) for v in x) + ")"
    x = float(x)
    if not np.isfinite(x):
        return "∞" if x > 0 else "−∞"
    if abs(x) >= 1000:
        return f"{x:.0f}"
    return f"{x:.{n}f}"


def at(rec, F):
    for r in rec:
        if abs(r["F"] - F) < 1e-9:
            return r
    raise KeyError(F)


def main(out=sys.stdout):
    V = json.loads((ROOT / "results" / "part_b_valley.json").read_text())
    L = json.loads((ROOT / "results" / "part_b_limits.json").read_text())
    A = json.loads((ROOT / "results" / "part_a_v3.json").read_text())
    p = lambda *a: print(*a, file=out)
    names = [n for n in NAMES if n in L["members"]]

    p("### Valley extent\n")
    p("| weights | χ²_min | 68% region: κ′ range, x′ range | 95% region: κ′ range, x′ range |")
    p("|---|---|---|---|")
    for w in ("session1", "data"):
        W = V["weights"][w]
        a, b = W["valley_68%"], W["valley_95%"]
        p(f"| {w} σ = {tuple(W['sigma'])} | {f(W['chi2_min'])} | {f(a['kappa_range'], 3)}, {f(a['x_eq_range'], 1)} nm | "
          f"{f(b['kappa_range'], 3)}, {f(b['x_eq_range'], 1)} nm |")

    p("\n### MaxEnt members across the κ₀ bracket (95% valley)\n")
    p("| weights | κ₀ | no bet (q = p₀): Δχ² | Gaussian member (κ′, x′) | D(q‖p₀) | I-projection member (λ₁, λ₂) | D(q‖p₀) | valley D range | affordable fraction (D ≤ 1.2) |")
    p("|---|---|---|---|---|---|---|---|---|")
    for w in ("session1", "data"):
        W = V["weights"][w]
        for k0, rec in W["per_kappa0"].items():
            g = rec["maxent_95%"] or {}
            ip = (W["iprojection"].get(k0) or {}).get("member") or {}
            p(f"| {w} | {k0} | {f(rec['no_bet_dchi2'], 1)} | ({f(g.get('kappa'), 3)}, {f(g.get('x_eq'), 2)}) | {f(g.get('D'), 3)} | "
              f"{f(ip.get('lambdas'), 4)} | {f(ip.get('D'), 3)} | {f(rec['D_range_95%'], 2)} | {f(rec['affordable_fraction_95%'], 2)} |")

    p("\n### Members used for the predictions (session-1 weights)\n")
    p("| member | landscape | q mean, SD (nm) | Δχ² | B (kT) | k_b0 (s⁻¹) | δ_b (nm) | k_c (s⁻¹) | k_on (µM⁻¹s⁻¹) | T (ms) | commitment D(q‖p₀), κ₀ 0.21 / 0.03 |")
    p("|---|---|---|---|---|---|---|---|---|---|---|")
    for n in names:
        q = L["members"][n]["params"]
        shape = (f"κ′ {f(q['kappa'], 3)}, x′ {f(q['x_eq'], 2)} nm" if q["family"] == "gaussian"
                 else f"p₀(κ₀ {q['kappa0']}) × soft rear wall λ = {f(q['lambdas'], 4)}")
        p(f"| {LABEL[n]} | {shape} | {f(q['q_mean_sd'], 2)} | {f(q['dchi2'])} | {f(q['B'])} | {f(q['kb0'])} | {f(q['delta_b'])} | "
          f"{f(q['kc'], 1)} | {f(q['kon'])} | {f(q['T'] * 1e3, 1)} | {f(q['commitment_k0=0.21'], 3)} / {f(q['commitment_k0=0.03'], 3)} |")

    p("\n### The imprint at F = 0 and its sensitivity to p (q held fixed), kT per committed step\n")
    p("| member | D(p‖q), p = start N(0.2, 1/0.21) | p mean −1.8 nm | p mean +2.2 nm | p SD 3.0 nm | p = undocked equilibrium, κ₀ 0.21 | p = undocked equilibrium, κ₀ 0.03 |")
    p("|---|---|---|---|---|---|---|")
    for n in names:
        q = L["members"][n]["params"]
        p(f"| {LABEL[n]} | {f(q['mismatch_k0=0.21'], 3)} | {f(q['mismatch_p_mean-1.8'], 3)} | {f(q['mismatch_p_mean+2.2'], 3)} | "
          f"{f(q['mismatch_p_sd3.0'], 3)} | {f(q['mismatch_p=undocked_k0=0.21'], 3)} | {f(q['mismatch_p=undocked_k0=0.03'], 3)} |")

    p("\n### (a) Load, 1 mM ATP (κ₀ = 0.21 for p and p₀; rest of cycle two equal stages)\n")
    for F in (0.0, 2.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0):
        pass
    p("| member | v(0) | v(4) | v(6) | v(7) | v(10) nm/s | r(2) | r(5) | r(5.75) | r(6.5) | odds(4) | odds(7) | odds(10) | 1:1 load (pN) | dwell CV(7) |")
    p("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for n in names:
        rec = L["members"][n]["load_k0=0.21"]
        Fs = np.array([r["F"] for r in rec])
        lo = np.array([np.log(r["odds"]) for r in rec])
        i = np.where(np.diff(np.sign(lo)) != 0)[0]
        stall = Fs[i[0]] - lo[i[0]] * (Fs[i[0] + 1] - Fs[i[0]]) / (lo[i[0] + 1] - lo[i[0]]) if len(i) else np.nan
        g = lambda F, k: at(rec, F)[k]
        p(f"| {LABEL[n]} | {f(g(0,'v'),0)} | {f(g(4,'v'),0)} | {f(g(6,'v'),1)} | {f(g(7,'v'),1)} | {f(g(10,'v'),1)} | "
          f"{f(g(2,'randomness'))} | {f(g(5,'randomness'))} | {f(g(5.75,'randomness'))} | {f(g(6.5,'randomness'))} | "
          f"{f(g(4,'odds'),1)} | {f(g(7,'odds'))} | {f(g(10,'odds'),3)} | {f(stall)} | {f(g(7,'dwell_cv'))} |")

    p("\n### (a′) Imprint kT·D(p‖q) vs load, 1 mM ATP (p does not depend on κ₀)\n")
    p("| member | per committed step at F = 0 / 2 / 4 / 6 / 7 / 8 / 10 / 12 pN | per net forward step at 2 / 6 pN | per second at 0 / 2 / 6 pN | upper bound (every ATP binding a quench): per step at 2 / 6 pN |")
    p("|---|---|---|---|---|")
    for n in names:
        rec = L["members"][n]["load_k0=0.21"]
        g = lambda F, k: at(rec, F)[k]
        p(f"| {LABEL[n]} | " + " / ".join(f(g(F, 'mismatch_per_step'), 3) for F in (0, 2, 4, 6, 7, 8, 10, 12)) + " | "
          f"{f(g(2,'mismatch_per_net_step'),3)} / {f(g(6,'mismatch_per_net_step'),2)} | "
          f"{f(g(0,'mismatch_per_s'),1)} / {f(g(2,'mismatch_per_s'),1)} / {f(g(6,'mismatch_per_s'),1)} | "
          f"{f(g(2,'mismatch_per_step_all_attempts'),3)} / {f(g(6,'mismatch_per_step_all_attempts'),2)} |")

    p("\n### (b) Viscosity (F = 0, 1 mM): v(η)/v(η₀)\n")
    etas = L["conditions"]["etas"]
    p("| member | gate | " + " | ".join(f"η/η₀ = {e:g}" for e in etas) + " |")
    p("|---|---|" + "---|" * len(etas))
    for n in names:
        for key, lab in (("viscosity", "η-independent"), ("viscosity_gate_viscous", "∝ 1/η")):
            rec = L["members"][n][key]["F=0.0"]
            v0 = rec[0]["v"]
            p(f"| {LABEL[n]} | {lab} | " + " | ".join(f(r["v"] / v0) for r in rec) + " |")

    p("\n### (c) Temperature, 1 mM ATP\n")
    temps = L["conditions"]["temps_C"]
    p("| member | quantity (variant) | " + " | ".join(f"{t:g} °C" for t in temps) + " |")
    p("|---|---|" + "---|" * len(temps))
    rec = L["members"]["best"]["temperature"]["default"]["F=5.0"]
    p("| all members (spread < 0.01 pN) | 1:1 load, pN (default: entropic landscape) | " + " | ".join(f(r["stall_1to1"]) for r in rec) + " |")
    for n in names:
        rec = L["members"][n]["temperature"]["enthalpic_tether"]["F=5.0"]
        p(f"| {LABEL[n]} | 1:1 load, pN (enthalpic tether) | " + " | ".join(f(r["stall_1to1"]) for r in rec) + " |")
    for n in names:
        for var in ("default", "rest_H10", "rest_H26"):
            rec = L["members"][n]["temperature"][var]["F=5.0"]
            p(f"| {LABEL[n]} | v at 5 pN, nm/s ({var}) | " + " | ".join(f(r["v"], 0) for r in rec) + " |")
        rec0 = L["members"][n]["temperature"]["default"]["F=0.0"]
        p(f"| {LABEL[n]} | v at 0 pN, nm/s (default) | " + " | ".join(f(r["v"], 0) for r in rec0) + " |")
        rec = L["members"][n]["temperature"]["default"]["F=5.0"]
        p(f"| {LABEL[n]} | randomness at 5 pN (default) | " + " | ".join(f(r["randomness"]) for r in rec) + " |")
        p(f"| {LABEL[n]} | imprint per s at 5 pN, kT/s (default) | " + " | ".join(f(r["mismatch_per_s"], 2) for r in rec) + " |")

    p("\n### (d) At the conditions of the limit datasets (values from `load_wide`)\n")
    W0 = L["members"]
    cond = [("1000", (-10.0, -5.0, 0.0, 6.0, 8.0, 10.0, 12.0, 14.0), "v", 1, "v (nm/s), 1 mM (Carter & Cross)"),
            ("10", (-10.0, -5.0, 0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0), "v", 1, "v (nm/s), 10 µM (Carter & Cross)"),
            ("2000", (0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 5.5, 6.0, 6.5), "v", 0, "v (nm/s), 2 mM (Visscher)"),
            ("2000", (0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 5.5, 5.75, 6.0), "randomness", 2, "r, 2 mM (Visscher)"),
            ("5", (0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 5.5), "v", 1, "v (nm/s), 5 µM (Visscher)"),
            ("1600", (-8.0, -6.0, -4.0, -2.0, 0.0, 2.0, 3.0, 4.0, 4.75), "randomness", 2, "r, 1.6 mM (Block)"),
            ("1600", (-8.0, -6.0, -4.0, -2.0, 0.0, 2.0, 3.0, 4.0, 4.75), "v", 0, "v (nm/s), 1.6 mM (Block)"),
            ("4.2", (-8.0, -6.0, -4.0, -2.0, 0.0, 2.0, 3.0, 4.0, 4.75), "randomness", 2, "r, 4.2 µM (Block)"),
            ("1000", (0.0, 2.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0), "odds", 3, "forward:back odds, 1 mM"),
            ("10", (0.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0), "odds", 3, "forward:back odds, 10 µM"),
            ("1000", (0.0, 2.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0), "dwell_cv", 2, "dwell CV, 1 mM"),
            ("1000", (0.0, 2.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0), "dwell_mean", 1, "mean dwell (ms), 1 mM")]
    for atp, Fs, key, nd, title in cond:
        p(f"\n{title}\n")
        p("| member | " + " | ".join(f"{F:g} pN" for F in Fs) + " |")
        p("|---|" + "---|" * len(Fs))
        for n in names:
            rec = W0[n]["load_wide"][atp]
            scale = 1e3 if key == "dwell_mean" else 1.0
            p(f"| {LABEL[n]} | " + " | ".join(f(at(rec, F)[key] * scale, nd) for F in Fs) + " |")
    p("\nRandomness vs [ATP] at fixed load (Visscher Fig. 4a conditions)\n")
    atps = L["conditions"]["atps"]
    pick = [0, 4, 8, 12, 14, 16, 18, 20, 24]
    p("| member | load | " + " | ".join(f"{atps[i]:.3g} µM" for i in pick) + " |")
    p("|---|---|" + "---|" * len(pick))
    for n in names:
        for Fk, vals in W0[n]["randomness_vs_atp"].items():
            p(f"| {LABEL[n]} | {Fk} pN | " + " | ".join(f(vals[i]) for i in pick) + " |")


if __name__ == "__main__" and len(sys.argv) == 1:
    main()


# ---------------------------------------------------------------------------
# B.4: which features near a limit tell the MaxEnt members apart?
# ---------------------------------------------------------------------------
# realistic precision for each observable (per condition; my estimates from the
# error bars in the folder, see PREDICTIONS.md): used only to judge detectability
PRECISION = {
    "v at 6 pN (nm/s)": 15.0,              # Visscher Fig. 3a s.e.m. near stall ~15-20 nm/s
    "1:1 load (pN)": 0.3,                  # stall-force s.e.m. ~0.2-0.6 pN (Visscher Fig. 3b, Taniguchi Table 1)
    "r at 5.75 pN": 0.10,                  # Visscher Fig. 4b s.e.m. at 5.76 pN ~0.10
    "r at 5 pN": 0.06,
    "ln odds at 10 pN": 0.3,               # counting error of ln(Nf/Nb) with ~50 events per bin
    "dwell CV at 7 pN": 0.05,              # ~400 dwells per bin
    "v(5 eta0)/v(eta0)": 0.05,             # Sozanski s.e.m. ~0.01-0.02 um/s on ~0.4-0.8 um/s
    "1:1 load, 35 C / 15 C": 0.05,         # ~0.3 pN on 7 pN
    "mismatch per step at 2 pN (kT)": 1.0,  # guessed Harada-Sasa precision (no folder value), see Part D
    "mismatch difference 6 pN - 2 pN (kT)": 1.0,
}


def features(L, name):
    rec = L["members"][name]["load_k0=0.21"]
    g = lambda F, k: at(rec, F)[k]
    Fs = np.array([r["F"] for r in rec])
    lo = np.array([np.log(r["odds"]) for r in rec])
    i = np.where(np.diff(np.sign(lo)) != 0)[0]
    stall = Fs[i[0]] - lo[i[0]] * (Fs[i[0] + 1] - Fs[i[0]]) / (lo[i[0] + 1] - lo[i[0]]) if len(i) else np.nan
    visc = L["members"][name]["viscosity"]["F=0.0"]
    eta = L["conditions"]["etas"]
    j5 = int(np.argmin(np.abs(np.array(eta) - 5.0)))
    tem = L["members"][name]["temperature"]["default"]["F=5.0"]
    temps = L["conditions"]["temps_C"]
    t15, t35 = temps.index(15.0), temps.index(35.0)
    return {
        "v at 6 pN (nm/s)": g(6, "v"),
        "1:1 load (pN)": stall,
        "r at 5 pN": g(5, "randomness"),
        "r at 5.75 pN": g(5.75, "randomness"),
        "ln odds at 10 pN": np.log(g(10, "odds")),
        "dwell CV at 7 pN": g(7, "dwell_cv"),
        "v(5 eta0)/v(eta0)": visc[j5]["v"] / visc[0]["v"],
        "1:1 load, 35 C / 15 C": tem[t35]["stall_1to1"] / tem[t15]["stall_1to1"],
        "mismatch per step at 2 pN (kT)": g(2, "mismatch_per_step"),
        "mismatch difference 6 pN - 2 pN (kT)": g(6, "mismatch_per_step") - g(2, "mismatch_per_step"),
    }


def discrimination(out=sys.stdout):
    L = json.loads((ROOT / "results" / "part_b_limits.json").read_text())
    p = lambda *a: print(*a, file=out)
    names = [n for n in NAMES if n in L["members"]]
    F = {n: features(L, n) for n in names}
    maxent = [n for n in names if n.startswith(("iproj", "maxent"))]
    others = [n for n in names if n not in maxent]
    p("\n### B.4 Discrimination: spread across members vs achievable precision\n")
    p("| feature | " + " | ".join(LABEL[n] for n in names) + " | spread (max−min) | MaxEnt vs others (max |Δ|) | precision | spread / precision |")
    p("|---|" + "---|" * len(names) + "---|---|---|---|")
    for key, prec in PRECISION.items():
        vals = np.array([F[n][key] for n in names], dtype=float)
        sp = np.nanmax(vals) - np.nanmin(vals)
        dm = max(abs(F[a][key] - F[b][key]) for a in maxent for b in others)
        p(f"| {key} | " + " | ".join(f(v, 3) for v in vals) + f" | {f(sp, 3)} | {f(dm, 3)} | {prec:g} | {f(sp / prec, 2)} |")


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "--disc":
    discrimination()
