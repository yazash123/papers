"""Part E: are kinesin's direction odds entropic?

Three readings of the odds O0 = k_f0/k_b0 and the lever l (stall F = kT ln O0 / l):
  R1  entropic odds, fixed lever:  ln O0 independent of T, l fixed  -> F proportional to T
  R2  odds and lever both entropic: ln O0 independent of T, l ~ T   -> F flat
  R3  enthalpic odds, fixed lever:  kT ln O0 independent of T       -> F flat

Data (all from the PDFs unless stated):
  Taniguchi et al. 2005, Table 1 (p.345): k_f0, k_b0, d_f, d_b at 280/287/298/308 K
      (bovine brain kinesin, 1 mM ATP, trap 0.035-0.038 pN/nm);
    Table 2 (p.345): dH_f = 18.3 +/- 1.1, dH_b = 18.2 +/- 1.4 kBT0, T0 = 298 K;
    p.343: "At low temperature (7 C) ... the maximum force was unchanged".
  Hong et al. 2016, Fig. 2 caption (p.1289): KIF5A force production
    5.3 +/- 0.2 pN at 295 K vs 5.2 +/- 0.2 pN at 280.5 K (detachment force in a
    fixed trap: grip-limited for KIF5A, see Part C, so not a clean stall).
  Kawaguchi & Ishiwata 2000 (BBRC 272:895): flat stall 7.34 +/- 0.33 pN over
    15-35 C -- AS QUOTED IN THE SESSION BRIEF; the paper is not in papers/.
"""
from __future__ import annotations

import numpy as np

from common import COLORS, MUTED, FIGURES, plot_style, save_json

from competing_exits.rates import K_B

T = np.array([280.0, 287.0, 298.0, 308.0])
KF0 = np.array([100.0, 209.0, 544.0, 1353.0]); SKF0 = np.array([11.0, 17.0, 54.0, 78.0])
KB0 = np.array([0.28, 0.61, 1.6, 3.8]); SKB0 = np.array([0.04, 0.06, 0.2, 0.3])
DF = np.array([2.4, 2.3, 2.6, 2.8]); SDF = np.full(4, 0.1)
DB = np.array([0.0, 0.1, 0.0, 0.0]); SDB = np.full(4, 0.1)
T0 = 298.0


def wls(x, y, s):
    """Weighted straight-line fit y = a + b x; returns a, b, cov, chi2."""
    w = 1 / s ** 2
    A = np.vstack([np.ones_like(x), x]).T
    cov = np.linalg.inv(A.T @ (A * w[:, None]))
    a, b = cov @ (A.T @ (w * y))
    chi2 = float(np.sum(w * (y - a - b * x) ** 2))
    return a, b, cov, chi2


def main():
    kT = K_B * T
    lnO = np.log(KF0 / KB0)
    s_lnO = np.sqrt((SKF0 / KF0) ** 2 + (SKB0 / KB0) ** 2)
    lev = DF - DB
    s_lev = np.sqrt(SDF ** 2 + SDB ** 2)
    F11 = kT * lnO / lev
    s_F11 = F11 * np.sqrt((s_lnO / lnO) ** 2 + (s_lev / lev) ** 2)

    # ln O0 = a + b (T0/T - 1): b = -ddH/(kB T0) (enthalpic part), a = ln O0 at T0
    x = T0 / T - 1
    a, b, cov, chi_lnO = wls(x, lnO, s_lnO)
    ddH_fit = -b  # in kB T0 units: ln O0 = ddS/kB - ddH/(kB T) = a - (-b)...
    # enthalpic odds (R3) would need all of ln O0(T0) to be enthalpic: b = ln O0(T0)
    z_R3 = (a - b) / np.sqrt(cov[0, 0] + cov[1, 1] - 2 * cov[0, 1])  # b vs a
    # lever: l = l0 (T/T0)^alpha  ->  ln l = ln l0 + alpha ln(T/T0)
    al0, alpha, cov_l, chi_l = wls(np.log(T / T0), np.log(lev), s_lev / lev)
    s_alpha = np.sqrt(cov_l[1, 1])
    chi_fixed = float(np.sum(((lev - np.average(lev, weights=1 / s_lev ** 2)) / s_lev) ** 2))
    chi_prop = float(np.sum(((lev - np.average(lev / T, weights=(T / s_lev) ** 2) * T) / s_lev) ** 2))

    # predictions for the stall between 15 and 35 C (288 -> 308 K)
    ratio_R1 = 308.0 / 288.0
    pred = {"R1_entropic_odds_fixed_lever": ratio_R1, "R2_both_entropic": 1.0, "R3_enthalpic_odds": 1.0}
    stall_R1 = {"15C": 6.9, "35C": 6.9 * ratio_R1}

    # the discriminating measurement: ln(forward/back) at a fixed load at 15 and 35 C
    def ln_ratio(F, Tk, reading, lnO0, l0):
        kTk = K_B * Tk
        if reading == "R1":
            return lnO0 - F * l0 / kTk
        if reading == "R2":
            return lnO0 - F * l0 * (Tk / T0) / kTk
        if reading == "R3":
            return lnO0 * (K_B * T0) / kTk - F * l0 / kTk
    disc = {}
    for label, lnO0, l0 in (("Drosophila (C&C: ln O0 6.69, l 3.88 nm)", np.log(802.0), 3.88),
                            ("bovine (Taniguchi: ln O0 5.85, l 2.6 nm)", 5.85, 2.6)):
        rows = {}
        for F in (0.0, 3.0, 5.0, 7.0):
            rows[f"{F:.0f} pN"] = {r: ln_ratio(F, 308.0, r, lnO0, l0) - ln_ratio(F, 288.0, r, lnO0, l0)
                                   for r in ("R1", "R2", "R3")}
        disc[label] = rows
    # counting error: sd(ln(Nf/Nb)) ~ sqrt(1/Nf + 1/Nb); backsteps needed per
    # temperature for a 3-sigma separation of R1 from R2 at 5 pN (Drosophila)
    delta = disc["Drosophila (C&C: ln O0 6.69, l 3.88 nm)"]["5 pN"]["R1"]
    odds5 = np.exp(np.log(802.0) - 5 * 3.88 / (K_B * 298.0))
    # need sqrt(2) * sqrt((1 + 1/odds)/Nb) <= delta / 3
    Nb_needed = 2 * (1 + 1 / odds5) * (3 / delta) ** 2

    out = {
        "taniguchi": {"T_K": T, "ln_O0": lnO, "s_ln_O0": s_lnO, "lever_nm": lev, "s_lever": s_lev,
                      "F_1to1_pN": F11, "s_F_1to1": s_F11},
        "ln_O0_fit": {"ln_O0_at_298K": a, "slope_b": b, "s_b": np.sqrt(cov[1, 1]), "chi2": chi_lnO,
                      "enthalpic_part_kBT0": -b, "z_enthalpic_reading_R3": z_R3,
                      "R3_needs_b_equal": a},
        "taniguchi_table2_ddH_kBT0": 0.1, "taniguchi_table2_s_ddH": float(np.hypot(1.1, 1.4)),
        "lever_fit": {"alpha": alpha, "s_alpha": s_alpha, "chi2_power_law": chi_l,
                      "chi2_fixed_lever_R1": chi_fixed, "chi2_lever_prop_T_R2": chi_prop, "dof": 3},
        "stall_ratio_35C_over_15C": pred, "stall_R1_example": stall_R1,
        "hong_KIF5A_force": {"295K": [5.3, 0.2], "280.5K": [5.2, 0.2],
                             "R1_prediction_280.5K": 5.3 * 280.5 / 295.0},
        "kawaguchi_ishiwata_as_quoted": [7.34, 0.33],
        "discriminating_measurement_dlnratio_35C_minus_15C": disc,
        "backsteps_per_temperature_for_3sigma_R1_vs_R2_at_5pN": Nb_needed,
        "odds_at_5pN_drosophila_298K": odds5,
    }
    save_json("part_e_temperature.json", out)
    print("ln O0:", np.round(lnO, 3), "+/-", np.round(s_lnO, 3))
    print("lever:", lev, " F11:", np.round(F11, 2), "+/-", np.round(s_F11, 2))
    print("ln O0 fit: a=%.3f b=%.3f+/-%.3f chi2=%.2f  z(R3)=%.1f" % (a, b, np.sqrt(cov[1, 1]), chi_lnO, z_R3))
    print("lever alpha=%.2f+/-%.2f chi2 power %.2f fixed %.2f propT %.2f" % (alpha, s_alpha, chi_l, chi_fixed, chi_prop))
    print("disc", {k: {f: {r: round(v, 3) for r, v in d.items()} for f, d in rows.items()} for k, rows in disc.items()})
    print("Nb needed", Nb_needed, "odds5", odds5)
    make_figure(out)


def make_figure(out):
    plt = plot_style()
    t = out["taniguchi"]
    Tk = np.array(t["T_K"])
    fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.3))
    TT = np.linspace(278, 310, 50)
    # ln O0
    a = out["ln_O0_fit"]["ln_O0_at_298K"]
    ax[0].errorbar(Tk, t["ln_O0"], yerr=t["s_ln_O0"], fmt="o", color="#0b0b0b", ms=5, capsize=2, label="Taniguchi Table 1")
    ax[0].plot(TT, np.full_like(TT, a), color=COLORS[0], lw=2, label="entropic odds (R1, R2)")
    ax[0].plot(TT, a * 298 / TT, color=COLORS[1], lw=2, ls="--", label="enthalpic odds (R3)")
    ax[0].set(xlabel="temperature (K)", ylabel="ln O₀ = ln(k_f0/k_b0)", ylim=(5.0, 6.6))
    ax[0].legend(fontsize=7.5, loc="lower left")
    # lever
    lev = np.array(t["lever_nm"])
    ax[1].errorbar(Tk, lev, yerr=t["s_lever"], fmt="o", color="#0b0b0b", ms=5, capsize=2, label="Taniguchi $d_f-d_b$")
    w = 1 / np.array(t["s_lever"]) ** 2
    lbar = np.average(lev, weights=w)
    ax[1].plot(TT, np.full_like(TT, lbar), color=COLORS[1], lw=2, ls="--", label="fixed lever (R1, R3)")
    c = np.average(lev / Tk, weights=w * Tk ** 2)
    ax[1].plot(TT, c * TT, color=COLORS[0], lw=2, label="lever ∝ T (R2)")
    ax[1].set(xlabel="temperature (K)", ylabel="lever ℓ (nm)", ylim=(1.8, 3.2))
    ax[1].legend(fontsize=7.5, loc="upper left")
    # stall
    ax[2].fill_between([288, 308], 7.34 - 0.33, 7.34 + 0.33, color="#e5e5e5", lw=0)
    ax[2].text(288.5, 7.72, "Kawaguchi & Ishiwata 2000:\n7.34 ± 0.33 pN, 15–35 °C\n(as quoted; not in folder)", fontsize=7, color=MUTED)
    ax[2].plot([288, 308], [6.9, 6.9 * 308 / 288], color=COLORS[1], lw=2, ls="--", label="R1: F ∝ T")
    ax[2].plot([288, 308], [7.1, 7.1], color=COLORS[0], lw=2, label="R2, R3: flat")
    ax[2].errorbar([295.0, 280.5], [5.3, 5.2], yerr=[0.2, 0.2], fmt="s", color=COLORS[2], ms=5, capsize=2,
                   label="Hong 2016 KIF5A (detachment force)")
    ax[2].set(xlabel="temperature (K)", ylabel="force (pN)", ylim=(4.5, 8.5), xlim=(278, 310))
    ax[2].legend(fontsize=7, loc="lower right")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig4_temperature.png")


if __name__ == "__main__":
    main()
