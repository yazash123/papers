"""Part C: is KIF5A's backstep gate sized to its grip?  (Kondo et al. 2023)

Model: one waiting state (Kondo's state-0) with forward step, slow backstep and
slow detachment (Table 1, p.470; 25 C).  Optionally the fast events: after a
slow backstep the motor enters state-3 with probability q and takes fast
backsteps / fast detachments (Table 1; scheme Fig. 4A, p.468).

The gate factor g divides the backstep rate (g > 1: stronger gate).

Shoulder definitions (the weakest gate for which detachment, not backsteps,
limits force):
  S95 / S90 : smallest g whose mean load at detachment is 95% / 90% of its
              value for an infinitely strong gate (g -> infinity);
  Scross    : g at which the 1:1 load equals the mean detachment load of an
              infinitely strong gate.
Where the measured gate (g = 1) sits is reported as log10(1/g_shoulder):
positive = stronger than the shoulder (in decades).

q (fast-event branching) is NOT reported by Kondo et al.; it is reconstructed:
  q = 0.22 from 'fast/slow backsteps ~15% at 3-8 pN' (p.470) -- default;
  q = 0.37 from the high-load fractions of Fig. 3B (fast back ~0.15, fast
           detach ~0.09 vs slow back ~0.5) -- READ OFF A FIGURE.
"""
from __future__ import annotations

import time

import numpy as np
from scipy.optimize import brentq

from common import COLORS, MUTED, FIGURES, plot_style, save_json

from competing_exits import Clamp, Trap, simulate_runs, trap_statistics
from competing_exits.kinesin import KONDO, kondo_kif5a, kondo_values

G_INF = 1e6
EARLIER = {100: (17.9, 6.7, 0.52), 10: (13.6, 6.7, 0.52), 1: (9.3, 6.4, 0.53),
           1 / 3: (7.3, 5.8, 0.56), 0.1: (5.1, 4.5, 0.64)}


def n_max_for(kappa, F_top=45.0, d=8.2):
    return int(np.ceil(F_top / (kappa * d))) + 5


def trap_point(g, kappa=0.05, **kw):
    m = kondo_kif5a(gate=g, **kw)
    r = trap_statistics(m, Trap(kappa), n_max=n_max_for(kappa))
    return m, r


def F11(m):
    return m.home_state.balance_point("forward", "back")


def detach_load(g, kappa=0.05, **kw):
    return trap_point(g, kappa, **kw)[1]["mean_load_at_termination"]


def shoulders(kappa=0.05, **kw):
    top = detach_load(G_INF, kappa, **kw)
    out = {"F_det_inf": top}
    for name, frac in (("S95", 0.95), ("S90", 0.90)):
        f = lambda lg: detach_load(10 ** lg, kappa, **kw) - frac * top
        out[name] = 10 ** brentq(f, -4, 6, xtol=1e-4)
    # gate whose 1:1 load equals the grip: ln(k_f g / k_b) = 0 at F = top
    s = kondo_kif5a(gate=1.0, **kw).home_state
    out["Scross"] = float(s.rate("back", top) / s.rate("forward", top))
    for name in ("S95", "S90", "Scross"):
        out["decades_above_" + name] = -np.log10(out[name])
    return out


# ----- constant-force clamp -------------------------------------------------
def clamp_max_force(g, **kw):
    """Largest clamp force at which a run still makes one net step on average
    (run displacement = d).  Above it the motor, on average, goes nowhere
    before letting go."""
    m = kondo_kif5a(gate=g, **kw)
    d = m.step_size
    f = lambda F: m.run_stats(F)["displacement"] - d
    if f(0.0) <= 0:
        return 0.0  # gate so weak that even unloaded runs make no net step
    return brentq(f, 0.0, 60.0, xtol=1e-8)


def clamp_shoulders(**kw):
    top = clamp_max_force(G_INF, **kw)
    out = {"F_max_inf": top}
    for name, frac in (("S95", 0.95), ("S90", 0.90)):
        f = lambda lg: clamp_max_force(10 ** lg, **kw) - frac * top
        out[name] = 10 ** brentq(f, -4, 6, xtol=1e-4)
        out["decades_above_" + name] = -np.log10(out[name])
    return out


# ----- parametric bootstrap over Table 1 ------------------------------------
def draw_values(rng, inflate=1.0):
    v = {}
    for k, val in KONDO.items():
        if isinstance(val, tuple):
            mu, se = val
            if k.startswith("lam"):
                # rates: lognormal with the same relative SE (keeps them positive
                # without clipping; matters for the poorly known fast rates)
                v[k] = mu * np.exp(rng.normal(0.0, inflate * se / mu) - 0.5 * (inflate * se / mu) ** 2)
            else:
                v[k] = rng.normal(mu, inflate * se)
        else:
            v[k] = val
    return v


def bootstrap(n, inflate=1.0, fast=False, seed=0):
    rng = np.random.default_rng(seed)
    rows = []
    for _ in range(n):
        vals = draw_values(rng, inflate)
        kw = {"values": vals}
        if fast:
            kw.update(fast=True, q=rng.uniform(0.15, 0.40), fast_sign=rng.choice([1.0, -1.0]))
        m, r = trap_point(1.0, **kw)
        sh = shoulders(**kw)
        rows.append({"F11": F11(m), "F_det": r["mean_load_at_termination"],
                     "time": r["mean_time"], **sh})
    keys = rows[0].keys()
    summary = {}
    for k in keys:
        x = np.array([row[k] for row in rows])
        summary[k] = {"median": np.median(x), "p2.5": np.percentile(x, 2.5),
                      "p97.5": np.percentile(x, 97.5), "p16": np.percentile(x, 16),
                      "p84": np.percentile(x, 84)}
    return summary, rows


def main():
    t0 = time.time()
    out = {}
    base = kondo_kif5a()
    s0 = base.home_state
    out["kT"] = s0.kT
    out["checks"] = {
        "odds_unloaded": float(np.exp(s0.log_odds("forward", "back", 0.0))),
        "ln_odds_unloaded": float(s0.log_odds("forward", "back", 0.0)),
        "kb0": float(s0.rate("back", 0.0)), "kb6_over_kb0": float(s0.rate("back", 6.0) / s0.rate("back", 0.0)),
        "kf0": float(s0.rate("forward", 0.0)), "kd0": float(s0.rate("detach", 0.0)),
        "crossing_kf_kb": F11(base), "v0": base.velocity(0.0),
        "lever_above_knee": float(s0.lever("forward", "back", 6.0)),
        "lever_below_knee": float(s0.lever("forward", "back", 1.0)),
        "zero_velocity_with_fast_q0.22": kondo_kif5a(fast=True, q=0.22).stall_force(),
        "zero_velocity_with_fast_q0.37": kondo_kif5a(fast=True, q=0.37).stall_force(),
        "zero_velocity_fast_q0.22_figure_sign": kondo_kif5a(fast=True, q=0.22, fast_sign=-1).stall_force(),
    }
    fm = kondo_kif5a(fast=True, q=0.22)
    r3 = fm.state("fast").prob("fast_back", 5.0)
    out["checks"]["fast_per_slow_backstep_q0.22"] = 0.22 * r3 / (1 - 0.22 * r3)
    out["checks"]["fast_detach_per_slow_backstep_q0.22"] = 0.22 * (1 - r3) / (1 - 0.22 * r3)

    # ----- the gate table: exact and Gillespie -------------------------------
    table = []
    for g, old in EARLIER.items():
        m, r = trap_point(g)
        sims = simulate_runs(m, Trap(0.05), n_runs=4000, seed=11)
        mf, rf = trap_point(g, fast=True, q=0.22)
        table.append({
            "gate": g, "earlier": old,
            "F11": F11(m), "F_det_exact": r["mean_load_at_termination"], "time_exact": r["mean_time"],
            "sd_F_det": r["sd_load_at_termination"],
            "F_det_sim": sims["load"].mean(), "F_det_sim_se": sims["load"].std(ddof=1) / np.sqrt(4000),
            "time_sim": sims["time"].mean(), "time_sim_se": sims["time"].std(ddof=1) / np.sqrt(4000),
            "fast_F_det": rf["mean_load_at_termination"], "fast_time": rf["mean_time"],
            "fast_zero_velocity": mf.stall_force(),
        })
    out["gate_table"] = table

    # ----- sweep over gate, several protocols --------------------------------
    gates = np.geomspace(0.01, 1e4, 49)
    sweep = {}
    variants = {
        "slow_k0.05": dict(kappa=0.05),
        "fast_q0.22_k0.05": dict(kappa=0.05, fast=True, q=0.22),
        "fast_q0.37_k0.05": dict(kappa=0.05, fast=True, q=0.37),
    }
    for kappa in (0.02, 0.09, 0.2, 0.5):
        variants[f"slow_k{kappa}"] = dict(kappa=kappa)
        variants[f"fast_q0.22_k{kappa}"] = dict(kappa=kappa, fast=True, q=0.22)
    for name, kw in variants.items():
        kappa = kw.pop("kappa")
        F_det = [detach_load(g, kappa, **kw) for g in gates]
        T_att = [trap_point(g, kappa, **kw)[1]["mean_time"] for g in gates]
        sweep[name] = {"kappa": kappa, "F_det": F_det, "time": T_att,
                       "shoulders": shoulders(kappa, **kw), **{k: v for k, v in kw.items()}}
    out["sweep_gates"] = gates
    out["fast_q0.22_gain_by_gate"] = {g: detach_load(g, 0.05, fast=True, q=0.22) for g in (1.0, 2.0, 3.0, 10.0, G_INF)}
    out["F11_vs_gate"] = [F11(kondo_kif5a(gate=g)) for g in gates]
    out["sweep"] = sweep

    # clamp protocol
    out["clamp"] = {
        "slow": {"F_max": [clamp_max_force(g) for g in gates], "shoulders": clamp_shoulders()},
        "fast_q0.22": {"F_max": [clamp_max_force(g, fast=True, q=0.22) for g in gates],
                       "shoulders": clamp_shoulders(fast=True, q=0.22)},
    }

    # ----- uncertainty -------------------------------------------------------
    boot = {}
    for label, kw in {"table_SE": dict(inflate=1.0), "SE_x5": dict(inflate=5.0),
                      "SE_x5_fast": dict(inflate=5.0, fast=True)}.items():
        boot[label], _ = bootstrap(300, seed=21, **kw)
    out["bootstrap"] = boot

    save_json("part_c_kondo.json", out)
    print(f"done in {time.time() - t0:.0f} s")
    print("checks:", {k: round(float(v), 4) for k, v in out["checks"].items()})
    for row in table:
        print({k: (round(float(v), 4) if not isinstance(v, tuple) else v) for k, v in row.items()})
    for name, sv in sweep.items():
        print(name, {k: round(float(v), 3) for k, v in sv["shoulders"].items()})
    for name, sv in out["clamp"].items():
        print("clamp", name, {k: round(float(v), 3) for k, v in sv["shoulders"].items()})
    for label, b in boot.items():
        print("boot", label, {k: [round(float(v["median"]), 3), round(float(v["p2.5"]), 3), round(float(v["p97.5"]), 3)]
                              for k, v in b.items()})
    make_figure(out)


def make_figure(out):
    plt = plot_style()
    gates = np.array(out["sweep_gates"])
    fig, ax = plt.subplots(1, 2, figsize=(10.5, 3.8), gridspec_kw={"width_ratios": [1.25, 1]})
    a = ax[0]
    a.semilogx(gates, out["F11_vs_gate"], color=MUTED, lw=1.5, ls="--")
    a.text(600, 16.2, "1:1 load\n(odds)", color=MUTED, fontsize=8.5, ha="left")
    styles = [("slow_k0.05", COLORS[0], "-", "o", "0.05 pN/nm, slow events"),
              ("fast_q0.22_k0.05", COLORS[1], "-", "s", "0.05 pN/nm, + fast events"),
              ("slow_k0.5", COLORS[2], "-", "^", "0.5 pN/nm, slow events")]
    for name, c, ls, mk, lab in styles:
        sv = out["sweep"][name]
        a.semilogx(gates, sv["F_det"], color=c, lw=2, ls=ls, marker=mk, ms=3.5, markevery=4, label=lab)
        s95 = sv["shoulders"]["S95"]
        a.plot([s95], [0.95 * sv["shoulders"]["F_det_inf"]], marker="|", ms=14, color=c, mew=2)
    a.axvline(1.0, color="#0b0b0b", lw=0.8)
    a.text(1.08, 1.0, "measured\ngate", fontsize=8.5, va="bottom")
    a.set(xlabel="gate strength g (backstep rate ÷ g)", ylabel="mean load at detachment (pN)",
          ylim=(0, 20), xlim=(gates[0], gates[-1]))
    a.legend(loc="upper left", fontsize=8)
    a.set_title("Trap, starting unloaded (| = 95% shoulder)", fontsize=9.5)

    b = ax[1]
    boot = out["bootstrap"]
    labels = [("table_SE", "Table 1 SEs"), ("SE_x5", "SEs × 5"), ("SE_x5_fast", "SEs × 5, + fast events")]
    y = np.arange(len(labels))[::-1]
    for yi, (key, lab) in zip(y, labels):
        for j, (sh, c, mk) in enumerate([("decades_above_S95", COLORS[0], "o"), ("decades_above_S90", COLORS[1], "s"),
                                         ("decades_above_Scross", COLORS[2], "^")]):
            st = boot[key][sh]
            yy = yi + 0.22 * (1 - j)
            b.plot([st["p2.5"], st["p97.5"]], [yy, yy], color=c, lw=2)
            b.plot([st["median"]], [yy], marker=mk, color=c, ms=6, ls="none",
                   label={0: "95% shoulder", 1: "90% shoulder", 2: "1:1 load = grip"}[j] if yi == y[0] else None)
    b.axvline(0, color="#0b0b0b", lw=0.8)
    b.set_yticks(y, [lab for _, lab in labels])
    b.set(xlabel="decades the measured gate sits above the shoulder")
    b.set_title("Uncertainty (median, 95% range)", fontsize=9.5)
    b.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=3, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig2_gate_vs_grip_kif5a.png")


if __name__ == "__main__":
    import json
    import sys
    from common import RESULTS
    if "--figure-only" in sys.argv:
        make_figure(json.loads((RESULTS / "part_c_kondo.json").read_text()))
    else:
        main()
