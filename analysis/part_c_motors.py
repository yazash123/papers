"""Part C(d): other motors on one gate-vs-grip picture.

For every motor we compute, with one protocol (0.05 pN/nm trap, starting
unloaded):

  grip  = mean load at detachment for an infinitely strong gate, from the
          motor's own gate-free velocity v(F) and detachment rate k_off(F).
          Represented as one waiting state with a forward exit of rate v(F)/d
          and a terminating exit k_off(F) (exact lattice solution);
  gate  = the load at which forward and backward stepping odds are 1:1, only
          where a backstep rate or backstep fraction was actually measured.

Sources and what was guessed are listed per motor in MOTORS below and in
REPORT.md.  kT = 4.1 pN nm wherever the paper gives no temperature.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

from common import COLORS, MUTED, FIGURES, plot_style, save_json

from competing_exits import Exit, ExitKind, Motor, Trap, WaitingState, trap_statistics
from competing_exits.kinesin import head_race_v2, kondo_kif5a

D = 8.2
KAPPA = 0.05
KT_RT = 4.1  # assumed where a paper does not state its temperature


class FromFunction:
    """Rate law wrapping a plain function of load (kT ignored)."""

    def __init__(self, f):
        self.f, self.cache = f, {}

    def __call__(self, F, kT):
        key = round(float(F), 9)
        if key not in self.cache:
            self.cache[key] = max(float(self.f(float(F))), 1e-12)
        return self.cache[key]


def grip(v_inf, k_off, kT=KT_RT, kappa=KAPPA):
    """Mean load at detachment in a trap for a motor with no backsteps."""
    s = WaitingState("w", (
        Exit("forward", ExitKind.PRODUCTIVE, FromFunction(lambda F: v_inf(F) / D), step=D),
        Exit("detach", ExitKind.TERMINATING, FromFunction(k_off)),
    ), kT=kT)
    r = trap_statistics(Motor.single(s, step_size=D), Trap(kappa), n_max=int(45 / (kappa * D)) + 5)
    return r["mean_load_at_termination"], r["mean_time"]


# ---- Drosophila kinesin-1 (truncated DmK, Andreasson et al. 2015 eLife) ----
# 3-state velocity model, Table 2 (p.10); hindering F > 0 here (their F < 0).
def dmk_velocity(F, delta1=4.6, Fi=26.0, kT=KT_RT):
    k1 = 4900.0 * np.exp(-F * delta1 / kT)
    k2 = 95.0
    k3 = 260.0 * np.exp((Fi - F) * 0.35 / kT)
    return D / (1 / k1 + 1 / k2 + 1 / k3)


def dmk_koff_fig6(F, kT=KT_RT):
    """Unbinding under hindering load, Fig. 6 fit (p.11): 1.11 exp(F 0.60/kT)."""
    return 1.11 * np.exp(F * 0.60 / kT)


NL_TABLE1 = {  # Table 1 (p.7): L0- (nm), deltaL- (nm), hindering loads -6..0 pN
    "DmK-WT": (1120.0, 2.0), "DmK-1AA": (360.0, 1.6), "DmK-2AA": (410.0, 1.8),
    "DmK-3AA": (320.0, 1.6), "DmK-5AA": (440.0, 2.4), "DmK-6AA": (270.0, 1.9),
}


def nl_koff(name, kT=KT_RT):
    L0, dL = NL_TABLE1[name]
    mut = name != "DmK-WT"
    v = (lambda F: dmk_velocity(F, 4.0, 0.0, kT)) if mut else (lambda F: dmk_velocity(F, kT=kT))
    return v, (lambda F: v(F) / (L0 * np.exp(-F * dL / kT)))


# ---- Kinesin-2 KIF3 (Andreasson et al. 2015 Curr Biol) ---------------------
# Fig. 3C (p.1168): per-head rates; 2 mM ATP as in the run-length data.
KIF3_HEAD = {
    "A": dict(k1=3.06, km1=12.0, k2=530.0, k3=67.9, k4=np.inf, k5=23.5, delta=1.09),
    "B": dict(k1=2.06, km1=124.0, k2=1730.0, k3=81.2, k4=439.0, k5=17.0, delta=2.28),
}
KIF3_RUN = {"KIF3A/B": (182.0, 1.7), "KIF3A/A": (102.0, 1.6), "KIF3B/B": (177.0, 1.3)}  # Fig. 4E L0, delta
KIF3_BACK_4PN = {"KIF3A/B": 0.06, "KIF3A/A": 0.08, "KIF3B/B": 0.03, "Kinesin-1": 0.03}  # p.1169, 5 uM ATP


def kif3_head(h, atp_uM=2000.0, gate=1.0, kT=KT_RT):
    p = KIF3_HEAD[h]
    post = 1 / p["k3"] + (0 if np.isinf(p["k4"]) else 1 / p["k4"])
    from competing_exits import Bell, constant
    return WaitingState(h, (
        Exit("forward", ExitKind.PRODUCTIVE, Bell(p["k2"], p["delta"]), step=D, cost_time=post),
        Exit("back", ExitKind.WRONG, constant(p["k5"] / gate), step=-D),
        Exit("unbind", ExitKind.RESTART, constant(p["km1"])),
    ), kT=kT, entry_time=1 / (p["k1"] * atp_uM))


def kif3_velocity(heads, F, gate=1.0, atp_uM=2000.0):
    """Heads alternate after every step (forward or back), Fig. 3A."""
    disp, time = 0.0, 0.0
    for h in heads:
        s = kif3_head(h, atp_uM, gate)
        P = s.splitting(F)
        Ps = P["forward"] + P["back"]
        disp += D * (P["forward"] - P["back"]) / Ps
        time += s.mean_dwell(F)
    return disp / time


def main():
    rows = []

    # KIF5A, Kondo (exact own model; single state with all three exits)
    m = kondo_kif5a(gate=1e6)
    r = trap_statistics(m, Trap(KAPPA))
    s = kondo_kif5a().home_state
    rows.append(dict(
        motor="KIF5A (mouse)", source="Kondo 2023 Table 1", gate_F11=s.balance_point("forward", "back"),
        gate_basis="fitted k_f, k_b (direct)", grip=r["mean_load_at_termination"],
        grip_basis="Kondo k_f, k_d", ln_O0=float(s.log_odds("forward", "back", 0.0)),
        p_back_4pN=float(s.prob("back", 4.0) / (s.prob("back", 4.0) + s.prob("forward", 4.0))),
    ))
    # check: fluid-free single-state representation gives the same grip
    g_check, _ = grip(lambda F: D * float(s.rate("forward", F)), lambda F: float(s.rate("detach", F)), kT=s.kT)

    # Drosophila kinesin-1: gate from C&C (v2), grip from Andreasson eLife
    v2 = head_race_v2(2000.0)
    Fs_v2 = v2.balance_point("forward", "back")
    g_dmk, t_dmk = grip(lambda F: dmk_velocity(F), dmk_koff_fig6)
    v_wt, koff_wt_vL = nl_koff("DmK-WT")
    g_dmk_vL, _ = grip(v_wt, koff_wt_vL)
    # Block-lab backstep fraction at 4 pN (Andreasson 2015b p.1169): 3%
    ln_odds4 = np.log(0.97 / 0.03)
    rows.append(dict(
        motor="Drosophila KHC", source="C&C 2005 (gate); Andreasson 2015 eLife (grip)",
        gate_F11=Fs_v2, gate_basis="v2 fit to C&C ratio (full-length, 23 C)",
        gate_F11_alt=[4.0 + 4.087 * ln_odds4 / 3.88, 4.0 + KT_RT * ln_odds4 / 2.22],
        gate_alt_basis="3% backsteps at 4 pN (truncated DmK, Block lab), lever 3.88 or 2.22 nm",
        grip=g_dmk, grip_alt_vL=g_dmk_vL, grip_basis="3-state v(F) + k_off Fig. 6 (1.11 s^-1, 0.60 nm)",
        ln_O0=float(np.log(v2.rate("forward", 0) / v2.rate("back", 0))),
        p_back_4pN=float(v2.prob("back", 4.0) / (v2.prob("back", 4.0) + v2.prob("forward", 4.0))),
        p_back_4pN_measured=0.03,
    ))

    # neck-linker mutants: detachment only
    for name in ("DmK-1AA", "DmK-2AA", "DmK-3AA", "DmK-5AA", "DmK-6AA"):
        v, koff = nl_koff(name)
        g, t = grip(v, koff)
        rows.append(dict(motor=name, source="Andreasson 2015 eLife Tables 1-2", gate_F11=None,
                         gate_basis="no backstep data; 'no processive backstepping' (p.13)",
                         grip=g, grip_basis="v/L from Table 1 and 3-state model (mutant)"))

    # kinesin-2 KIF3
    for name, heads in (("KIF3A/B", "AB"), ("KIF3A/A", "AA"), ("KIF3B/B", "BB")):
        L0, dL = KIF3_RUN[name]
        v_inf = lambda F, h=heads: kif3_velocity(h, F, gate=1e9)
        v_g = lambda F, h=heads: kif3_velocity(h, F)
        koff = lambda F, v=v_g, L0=L0, dL=dL: max(v(F), 1e-6) / (L0 * np.exp(-F * dL / KT_RT))
        g, t = grip(v_inf, koff)
        F0 = brentq(v_g, 0.0, 40.0)
        per_head = {h: KT_RT * np.log(KIF3_HEAD[h]["k2"] / KIF3_HEAD[h]["k5"]) / KIF3_HEAD[h]["delta"] for h in sorted(set(heads))}
        pb = []
        for h in heads:
            st = kif3_head(h)
            pb.append(st.prob("back", 4.0) / (st.prob("back", 4.0) + st.prob("forward", 4.0)))
        p_meas = KIF3_BACK_4PN[name]
        lever = np.mean([KIF3_HEAD[h]["delta"] for h in heads])
        rows.append(dict(
            motor=name, source="Andreasson 2015 Curr Biol Fig. 3C, 4E", gate_F11=F0,
            gate_F11_alt=[F0, 4.0 + KT_RT * np.log((1 - p_meas) / p_meas) / lever],
            gate_alt_basis="measured backstep fraction at 4 pN (p.1169) with the fitted lever",
            gate_basis="zero velocity of fitted cycle (k5 backsteps)", gate_per_head=per_head,
            grip=g, grip_basis="v/L, L0 and delta from Fig. 4E (loaded regime)",
            p_back_4pN=float(np.mean(pb)), p_back_4pN_measured=KIF3_BACK_4PN[name],
            unloaded_velocity=kif3_velocity(heads, 0.0),
        ))

    # Gicking 2022 simulation inputs (Table 1, p.12): linear F-V to Fs = 6 pN,
    # k_off = k0 exp(F/Fdetach); backstep 3/s for all -> not a gate measurement
    for name, V0, k0, Fd in (("Kin1 (Gicking)", 586.0, 0.96, 6.8), ("Kin2 KIF3A/A (Gicking)", 307.0, 0.76, 3.0),
                             ("Kin3 KIF1A (Gicking)", 910.0, 0.16, 1.3)):
        g, t = grip(lambda F, V0=V0: V0 * max(1 - F / 6.0, 1e-9), lambda F, k0=k0, Fd=Fd: k0 * np.exp(F / Fd))
        rows.append(dict(motor=name, source="Gicking 2022 Table 1 (model inputs)", gate_F11=None,
                         gate_basis="k_backstep = 3/s assumed for all three", grip=g,
                         grip_basis="linear F-V (Fs = 6 pN assumed), Bell k_off"))

    # Budaitis 2019: mean +/- SD detachment forces (text p.6, Fig. 3B p.9), trap stiffness not given
    budaitis = {"KIF5C WT": [4.6, 0.8], "CNB": [0.91, 0.6], "Latch": [0.84, 0.4], "CNB+Latch": [0.81, 0.5]}

    out = {"rows": rows, "kondo_grip_single_state_check": g_check, "budaitis_observed": budaitis,
           "kappa": KAPPA, "kT_assumed": KT_RT}
    save_json("part_c_motors.json", out)
    for r_ in rows:
        print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r_.items()
               if k in ("motor", "gate_F11", "grip", "p_back_4pN", "p_back_4pN_measured", "gate_per_head", "gate_F11_alt", "grip_alt_vL", "unloaded_velocity")})
    print("kondo grip check", g_check)
    make_figure(out)


def make_figure(out):
    plt = plot_style()
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    lim = 16
    ax.plot([0, lim], [0, lim], color=MUTED, lw=0.8, ls="--")
    ax.text(12.3, 13.3, "1:1 load = grip", color=MUTED, fontsize=8, rotation=35)
    marks = {"KIF5A (mouse)": (COLORS[0], "o"), "Drosophila KHC": (COLORS[1], "s"),
             "KIF3A/B": (COLORS[2], "^"), "KIF3A/A": (COLORS[2], "v"), "KIF3B/B": (COLORS[2], "D")}
    for r in out["rows"]:
        if r["gate_F11"] is None:
            continue
        c, mk = marks[r["motor"]]
        ax.plot(r["grip"], r["gate_F11"], marker=mk, color=c, ms=8, ls="none", mec="white", mew=1.2)
        if "gate_F11_alt" in r:
            ax.plot([r["grip"]] * 2, r["gate_F11_alt"], color=c, lw=1.5)
        dx, dy, ha = {"KIF3A/A": (0.25, -0.3, "left"), "KIF3B/B": (0.3, -0.2, "left"),
                      "KIF3A/B": (-0.3, -0.2, "right"), "Drosophila KHC": (-0.3, -0.6, "right")}.get(
            r["motor"], (0.25, 0.35, "left"))
        ax.text(r["grip"] + dx, r["gate_F11"] + dy, r["motor"], fontsize=8.5, color="#0b0b0b", ha=ha)
    # detachment-only motors as ticks along the bottom
    y0 = 0.6
    for r in out["rows"]:
        if r["gate_F11"] is not None:
            continue
        ax.plot(r["grip"], y0, marker="|", color=MUTED, ms=12, mew=1.5)
    names = [r for r in out["rows"] if r["gate_F11"] is None]
    ax.text(0.2, 1.3, "grip only (no gate data): NL mutants, Gicking Kin1-3", fontsize=7.5, color=MUTED)
    ax.text(0.2, 15.2, "bars: 1:1 load implied by measured backstep\nfractions at 4 pN (Block lab)", fontsize=7.5, color=MUTED, va="top")
    ax.set(xlim=(0, lim), ylim=(0, lim), xlabel="grip: mean load at detachment, infinite gate,\n0.05 pN/nm trap (pN)",
           ylabel="gate: load at 1:1 forward:back odds (pN)")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig3_gate_vs_grip_motors.png")


if __name__ == "__main__":
    main()
