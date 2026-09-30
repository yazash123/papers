"""Re-run of the Part B predictions at the exact conditions of Ariga, Tomishige &
Mizuno 2018 (PRL 121:218101), after the paper itself became available.

POST HOC with respect to PREDICTIONS.md (commit 84a3d19): the model, the seven
members and every parameter are unchanged; only the conditions are Ariga's:
    F = 2 pN hindering; 25 +- 1 C;
    high ATP: 1 mM ATP, 0.1 mM ADP, 1 mM Pi;  low ATP: 10 uM ATP, 1 uM ADP, 1 mM Pi;
    Delta mu = 84.5 +- 2.5 pN nm in both (p.2).
v3 has no ADP/Pi dependence and was fitted to Carter & Cross (Drosophila full-length
KHC, 23 C); Ariga used human cysteine-light kinesin-1 truncated at residue 490
(Suppl. p.7).  Temperature: the native 23 C model, and 25 C under two readings of
Taniguchi's enthalpies (as committed, and corrected; REPORT2 C.6).

Ariga's measured energy flows (Table I, p.5; mean +- s.d.; pN nm/s):
                          high ATP (n = 8)   low ATP (n = 11)
    power F0 <v>          1150 +- 120        410 +- 60
    gamma <v>^2           10.6 +- 1.9        1.35 +- 0.37
    Harada-Sasa integral  53.4 +- 41.4       2.74 +- 1.52   (to 300 / 50 Hz)
    J_x                   63.9 +- 41.5       4.09 +- 1.56
    Delta mu / tau        6160 +- 560        2190 +- 310    (tau = d/<v>, d = 8 nm)
Head-motion dissipation estimated by the authors (Suppl. p.14): ~500 (high ATP),
~5 pN nm/s (low ATP).

Output: results/ariga_rerun.json, results/tables_ariga.md.
Run: cd analysis && PYTHONPATH=.. python ariga_rerun.py
"""
from __future__ import annotations

import json

import numpy as np

import part_b_limits as PL
from common import RESULTS, save_json
from competing_exits.headrace import IProjLandscape, V3_FIXED, head_race_v3
from competing_exits.rates import K_B
from partb_common import ROOT, Member, bookkeeping, costs, gaussian_member
from v3fit import pack

NAMES = ["iproj_k0=0.03", "iproj_k0=0.21", "maxent_k0=0.03", "maxent_k0=0.21", "best", "middle", "floppy_edge"]
MAXENT = NAMES[:4]
F = 2.0
DMU = 84.5            # pN nm (Ariga p.2)
DMU_ERR = 2.5
D_ARIGA = 8.0         # nm, Ariga's step size for tau = d/<v>
ARIGA = {  # Table I, experiment columns; pN nm/s
    "1000": dict(power=(1150, 120), gv2=(10.6, 1.9), hs=(53.4, 41.4), Jx=(63.9, 41.5), dmu_rate=(6160, 560),
                 head=500.0, n=8),
    "10": dict(power=(410, 60), gv2=(1.35, 0.37), hs=(2.74, 1.52), Jx=(4.09, 1.56), dmu_rate=(2190, 310),
               head=5.0, n=11),
}
ARIGA_FIT = {"1000": dict(kf0=981, kb0=22.8, kc=129, df=3.3, db=0.47),
             "10": dict(kf0=889, kb0=0.61, kc=32.5, df=4.0, db=-0.83)}   # Fig. 3 legend, p.3


def members():
    L = json.loads((ROOT / "results" / "part_b_limits.json").read_text())
    out = {}
    for n in NAMES:
        p = L["members"][n]["params"]
        th = pack({k: p[k] for k in ("B", "kb0", "delta_b", "kc", "kon", "T")})
        if p["family"] == "gaussian":
            out[n] = gaussian_member(n, p["kappa"], p["x_eq"], th, p["chi2"], (0.1, 0.2))
        else:
            Lp = IProjLandscape(p["kappa0"], tuple(p["lambdas"]), tuple(p["fs"]), B_front=p["B"])
            out[n] = Member(n, Lp, th, p["chi2"], (0.1, 0.2))
    return out


def ariga_hidden(atp):
    """Hidden dissipation per step from Table I: Delta mu - F d - J_x d/<v> (pN nm)."""
    a = ARIGA[atp]
    v = a["power"][0] / F
    jx_step = a["Jx"][0] * D_ARIGA / v
    jx_err = a["Jx"][1] * D_ARIGA / v
    h = DMU - F * D_ARIGA - jx_step
    return {"v_nm_s": v, "Jx_per_step": jx_step, "hidden_pNnm": h, "hidden_err": float(np.hypot(DMU_ERR, jx_err)),
            "hidden_fraction": h / DMU, "work_per_step": F * D_ARIGA}


def at_25C(member, atp, variant):
    """The member at 25 C (298.15 K) under a temperature reading, at load F and [ATP]."""
    saved = (PL.H_B, PL.H_GATE, PL.H_OTHER)
    try:
        if variant == "corrected":
            # Taniguchi's Delta H excludes the T/eta prefactor (REPORT2 C.6): barrier 18.3 on top
            # of D; gate raw 26.4; rest of cycle 13.8 (their k_c); clock, k_on as committed
            PL.H_B, PL.H_GATE = PL.H_FORWARD, 26.4
            m = PL.member_at_temperature(member, 25.0, H_rest=13.8)
        else:
            m = PL.member_at_temperature(member, 25.0)
    finally:
        PL.H_B, PL.H_GATE, PL.H_OTHER = saved
    s = head_race_v3(0, 0, 0, m["kb0"], m["delta_b"], m["kc"], m["kon"], m["T_rest"], search=m["search"],
                     atp_uM=float(atp), T_cv2=PL.T_CV2, kT=m["kT"])
    ts = s.transport_stats(F, step_size=V3_FIXED["d"])
    P = s.splitting(F)
    c = costs(m["member"], 0.21, F, kT=m["kT"])
    step_rate = (P["forward"] + P["back"]) / ts["cycle_time"]
    return {"v": ts["v"], "step_rate": step_rate, "back_fraction": P["back"] / (P["forward"] + P["back"]),
            "randomness": ts["randomness"], "imprint_kT": c["mismatch"], "kT": m["kT"]}


def summarise(v, step_rate, back_fraction, randomness, imprint_kT, kT):
    """Energy flows in Ariga's units (pN nm/s) and per step."""
    net = 1 - 2 * back_fraction                  # net forward steps per committed step
    work_step = F * V3_FIXED["d"] * net         # pN nm per committed step
    total_step = DMU - work_step                 # dissipated per committed step (one ATP each)
    imprint = imprint_kT * kT
    return {"v": v, "power": F * v, "dmu_rate": DMU * step_rate, "step_rate": step_rate, "back_fraction": back_fraction,
            "total_dissipation_per_step_pNnm": total_step, "total_dissipation_per_step_kT": total_step / kT,
            "imprint_per_step_pNnm": imprint, "imprint_per_step_kT": imprint_kT, "imprint_rate_pNnm_s": imprint * step_rate,
            "imprint_fraction_of_total": imprint / total_step,
            "tur_bound_per_net_step_kT": 2.0 / randomness, "randomness": randomness}


def main():
    M = members()
    out = {"conditions": {"F_pN": F, "dmu_pNnm": [DMU, DMU_ERR], "T_C": 25.0, "atp_uM": [1000, 10]},
           "ariga_table1": ARIGA, "ariga_fit_fig3": ARIGA_FIT, "ariga_hidden": {a: ariga_hidden(a) for a in ARIGA},
           "members": {}}
    for n, m in M.items():
        rec = {}
        for atp in ("1000", "10"):
            b = bookkeeping(m, 0.21, F, atp_uM=float(atp), T_cv2=PL.T_CV2)
            rec[f"23C_{atp}"] = summarise(b["v"], b["step_rate"], b["P_back"] / (b["P_forward"] + b["P_back"]),
                                          b["randomness"], b["mismatch_per_step"], V3_FIXED["kT"])
            for var in ("committed", "corrected"):
                t = at_25C(m, atp, var)
                rec[f"25C_{var}_{atp}"] = summarise(t["v"], t["step_rate"], t["back_fraction"], t["randomness"],
                                                    t["imprint_kT"], t["kT"])
        out["members"][n] = rec
    # Ariga's own two-state fits: implied odds and 1:1 loads (Bell, kb rising with load when db > 0)
    kT25 = K_B * 298.15
    out["ariga_fit_implied"] = {a: {"odds_at_2pN": p["kf0"] * np.exp(-F * p["df"] / kT25) / (p["kb0"] * np.exp(F * p["db"] / kT25)),
                                    "one_to_one_pN": kT25 * np.log(p["kf0"] / p["kb0"]) / (p["df"] + p["db"])}
                                for a, p in ARIGA_FIT.items()}
    save_json("ariga_rerun.json", out)
    tables(out)
    return out


def tables(out):
    lines = []
    p = lines.append
    p("### Ariga 2018 re-run (post hoc; conditions: 2 pN, 25 °C, Δμ = 84.5 pN·nm)\n")
    p("Ariga's hidden dissipation per step, from Table I: Δμ − F·d − J_x·d/⟨v⟩ (d = 8 nm, as they use)\n")
    p("| condition | ⟨v⟩ (nm/s) | work per step (pN·nm) | J_x per step (pN·nm) | hidden per step (pN·nm) | fraction of Δμ |")
    p("|---|---|---|---|---|---|")
    for a, lab in (("1000", "1 mM ATP"), ("10", "10 µM ATP")):
        h = out["ariga_hidden"][a]
        p(f"| {lab} | {h['v_nm_s']:.0f} | {h['work_per_step']:.1f} | {h['Jx_per_step']:.2f} | "
          f"{h['hidden_pNnm']:.1f} ± {h['hidden_err']:.1f} | {h['hidden_fraction']:.3f} |")
    for key, title in (("23C", "native model, 23 °C"), ("25C_committed", "25 °C, temperature model as committed"),
                       ("25C_corrected", "25 °C, Taniguchi's enthalpies read correctly")):
        for a, lab, ref in (("1000", "1 mM", ARIGA["1000"]), ("10", "10 µM", ARIGA["10"])):
            p(f"\n#### {lab} ATP, {title}\n")
            p(f"Ariga (experiment): power {ref['power'][0]} ± {ref['power'][1]}, Δμ/τ {ref['dmu_rate'][0]} ± {ref['dmu_rate'][1]} pN·nm/s; "
              f"authors' head-motion estimate ≈ {ref['head']:g} pN·nm/s\n")
            p("| member | v (nm/s) | power F·v (pN·nm/s) | Δμ × step rate (pN·nm/s) | backstep fraction | total dissipated per step (kT) | imprint per step (kT) | imprint rate (pN·nm/s) | imprint / total | TUR 2/r (kT per net step) |")
            p("|---|---|---|---|---|---|---|---|---|---|")
            for n in NAMES:
                r = out["members"][n][f"{key}_{a}"]
                p(f"| {n} | {r['v']:.0f} | {r['power']:.0f} | {r['dmu_rate']:.0f} | {r['back_fraction']:.4f} | "
                  f"{r['total_dissipation_per_step_kT']:.2f} | {r['imprint_per_step_kT']:.3f} | {r['imprint_rate_pNnm_s']:.1f} | "
                  f"{r['imprint_fraction_of_total']:.4f} | {r['tur_bound_per_net_step_kT']:.2f} |")
    p("\n#### Ariga's two-state fits (Fig. 3 legend), implied at 25 °C\n")
    p("| condition | odds f:b at 2 pN | 1:1 load (pN) |")
    p("|---|---|---|")
    for a, lab in (("1000", "1 mM"), ("10", "10 µM")):
        r = out["ariga_fit_implied"][a]
        p(f"| {lab} | {r['odds_at_2pN']:.1f} | {r['one_to_one_pN']:.2f} |")
    (RESULTS / "tables_ariga.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
