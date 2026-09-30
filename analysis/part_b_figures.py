"""Figures for Part B (predictions; no data plotted).

fig5_valley.png   the valley of docked landscapes, the commitment D(q||p0) of the
                  MaxEnt member across the kappa0 bracket against the docking budget,
                  and where each member puts the head
fig6_limits.png   predictions near the limits (load, viscosity, temperature) and the
                  mismatch cost per committed step
Run after part_b_maxent.py and part_b_limits.py:
    cd analysis && PYTHONPATH=.. python part_b_figures.py
"""
from __future__ import annotations

import json

import numpy as np

from common import COLORS, FIGURES, MUTED, plot_style
from competing_exits import bet
from partb_common import BUDGET, P_KAPPA, P_MEAN, ROOT, X

ORDER = ["iproj_k0=0.03", "iproj_k0=0.21", "maxent_k0=0.03", "maxent_k0=0.21", "best", "middle", "floppy_edge"]
LABEL = {"iproj_k0=0.03": "I-projection κ₀ 0.03", "iproj_k0=0.21": "I-projection κ₀ 0.21",
         "maxent_k0=0.03": "Gaussian MaxEnt κ₀ 0.03", "maxent_k0=0.21": "Gaussian MaxEnt κ₀ 0.21",
         "best": "best fit", "middle": "mid-valley", "floppy_edge": "floppy edge"}
# MaxEnt members solid, the other members dashed (identity never by colour alone)
STYLE = {n: dict(color=COLORS[i], lw=1.8, ls="-" if n.startswith(("iproj", "maxent")) else "--")
         for i, n in enumerate(ORDER)}
INK = "#0b0b0b"


def load():
    V = json.loads((ROOT / "results" / "part_b_valley.json").read_text())
    L = json.loads((ROOT / "results" / "part_b_limits.json").read_text())
    A = json.loads((ROOT / "results" / "part_a_v3.json").read_text())
    return V, L, A


def q_density(params):
    if params["family"] == "gaussian":
        return bet.tilted(bet.tether_logdensity(params["kappa"], params["x_eq"]), X)
    return bet.tilted(bet.iprojection_logdensity(params["kappa0"], params["lambdas"], params["fs"]), X)


def legend(fig, plt, ncol=4, y=-0.01):
    handles = [plt.Line2D([], [], **STYLE[n]) for n in ORDER]
    fig.legend(handles, [LABEL[n] for n in ORDER], loc="lower center", ncol=ncol, fontsize=7.5,
               bbox_to_anchor=(0.5, y))


def fig_valley(V, L, A):
    plt = plot_style()
    fig, axs = plt.subplots(1, 3, figsize=(12.8, 4.3))
    # (a) the valley (session-1 weights) and commitment contours
    ax = axs[0]
    k = np.array(A["grid"]["kappa"])
    x = np.array(A["grid"]["x_eq"])
    W = V["weights"]["session1"]
    d = np.array(A["session1"]["chi2_grid"], dtype=float) - W["chi2_min"]
    XX, KK = np.meshgrid(x, k)
    dd = np.where(np.isfinite(d), d, 99)
    ax.contourf(XX, KK, dd, levels=[0, 2.30, 5.99], colors=["#bcd4f1", "#dce8f7"])
    ax.contour(XX, KK, dd, levels=[2.30, 5.99], colors=[COLORS[0]], linewidths=[0.9, 1.5])
    for kappa0, ls in ((0.21, "--"), (0.03, ":")):
        D = np.array(W["per_kappa0"][f"{kappa0:g}"]["D_grid"], dtype=float)
        c2 = ax.contour(XX, KK, D, levels=[0.1, 0.3, BUDGET], colors=[MUTED], linewidths=0.9, linestyles=ls)
        ax.clabel(c2, fmt=lambda v: f"{v:g}", fontsize=6.5)
    ax.plot([0.0], [0.21], "s", ms=6, color=MUTED, mec="white")
    ax.annotate("p₀ (κ₀ 0.21)", (0.0, 0.21), xytext=(-8, -12), textcoords="offset points", fontsize=7, color=INK)
    for n in ORDER:
        p = L["members"][n]["params"]
        if p["family"] != "gaussian":
            continue
        ax.plot([p["x_eq"]], [p["kappa"]], "o", ms=6.5, color=STYLE[n]["color"], mec="white", mew=1.0, zorder=6)
    ax.set_yscale("log")
    ax.set(xlabel="docked tether centre x′ (nm)", ylabel="docked stiffness κ′ (kT/nm²)", xlim=(-2, 8), ylim=(0.04, 0.65))
    ax.set_title("a  valley (Δχ² ≤ 2.30 dark, ≤ 5.99 light) with\nD(q‖p₀) contours in kT (dashed κ₀ 0.21, dotted 0.03)",
                 fontsize=8.5)
    # (b) commitment of the MaxEnt member vs kappa0, against the budget
    ax = axs[1]
    ax.axhspan(1.0, 2.0, color="#eeeeee", zorder=0)
    ax.axhline(BUDGET, color=MUTED, lw=1.0)
    ax.annotate("docking budget ~1.2 kT (range 1–2)", (0.0105, BUDGET * 1.12), fontsize=7, color=INK)
    floor = 1e-3
    for w, ls, mk in (("session1", "-", "o"), ("data", "--", "s")):
        Wd = V["weights"][w]
        k0s = np.array([float(s) for s in Wd["per_kappa0"]])
        g = np.array([Wd["per_kappa0"][s]["maxent_95%"]["D"] for s in Wd["per_kappa0"]])
        ip = np.array([Wd["iprojection"][s]["member"]["D"] for s in Wd["iprojection"]])
        ax.plot(k0s, np.maximum(g, floor), ls=ls, marker=mk, ms=4, color=COLORS[2], lw=1.4)
        ax.plot(k0s, np.maximum(ip, floor), ls=ls, marker=mk, ms=4, color=COLORS[0], lw=1.4)
    ax.annotate("p₀ itself inside the valley\n(data weights, κ₀ 0.21)", (0.21, floor), xytext=(-60, 22),
                textcoords="offset points", fontsize=6.8, color=INK, arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.7))
    ax.set(xscale="log", yscale="log", ylim=(floor * 0.8, 4), xlabel="undocked tether stiffness κ₀ (kT/nm²)",
           ylabel="least commitment D(q‖p₀) in the valley (kT)")
    ax.set_title("b  MaxEnt member: green Gaussian family, blue I-projection;\nsolid session-1 weights, dashed data weights",
                 fontsize=8.5)
    # (c) where the bet puts the head
    ax = axs[2]
    p0 = bet.tilted(bet.tether_logdensity(0.21, 0.0), X)
    pp = bet.tilted(bet.tether_logdensity(P_KAPPA, P_MEAN), X)
    ax.plot(X, p0, color=MUTED, lw=1.3, ls=":")
    ax.plot(X, pp, color=INK, lw=1.0, ls=(0, (1, 1)))
    for n in ORDER:
        ax.plot(X, q_density(L["members"][n]["params"]), **STYLE[n])
    ax.annotate("p₀, κ₀ 0.21 (dotted grey)\np, the start (dotted black)", (-7.6, 0.165), fontsize=7, color=INK)
    ax.axvline(8.0, color=MUTED, lw=0.8)
    ax.annotate("front site", (7.9, 0.02), ha="right", fontsize=7, color=INK)
    ax.set(xlabel="free-head position x (nm)", ylabel="density (1/nm)", xlim=(-8, 8.3), ylim=(0, 0.3))
    ax.set_title("c  where each member's docked landscape\nputs the head (q, F = 0)", fontsize=8.5)
    legend(fig, plt)
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    fig.savefig(FIGURES / "fig5_valley.png")
    plt.close(fig)


def fig_limits(L):
    plt = plot_style()
    fig, axs = plt.subplots(2, 3, figsize=(12.8, 7.2))
    for n in ORDER:
        M = L["members"][n]
        wide = M["load_wide"]
        F = np.array([r["F"] for r in wide["1000"]])
        sel = (F >= 0) & (F <= 12)
        axs[0, 0].plot(F[sel], [r["v"] for r, s in zip(wide["1000"], sel) if s], **STYLE[n])
        F2 = np.array([r["F"] for r in wide["2000"]])
        r2 = np.array([r["randomness"] for r in wide["2000"]])
        s2 = (F2 >= 0) & (F2 <= 6.8)
        axs[0, 1].plot(F2[s2], r2[s2], **STYLE[n])
        odds = np.array([r["odds"] for r in wide["1000"]])
        s3 = (F >= 0) & (F <= 15)
        axs[0, 2].plot(F[s3], odds[s3], **STYLE[n])
        rec = M["load_k0=0.21"]
        Fm = np.array([r["F"] for r in rec])
        axs[1, 0].plot(Fm, [r["mismatch_per_step"] for r in rec], **STYLE[n])
        vis = M["viscosity"]["F=0.0"]
        eta = L["conditions"]["etas"]
        axs[1, 1].plot(eta, [r["v"] / vis[0]["v"] for r in vis], **STYLE[n])
        tem = M["temperature"]["default"]["F=5.0"]
        tem_e = M["temperature"]["enthalpic_tether"]["F=5.0"]
        T = [r["T_C"] for r in tem]
        axs[1, 2].plot(T, [r["stall_1to1"] for r in tem_e], **STYLE[n])
    # the entropic default is one curve for all members (1:1 load proportional to T)
    tem = L["members"]["best"]["temperature"]["default"]["F=5.0"]
    axs[1, 2].plot([r["T_C"] for r in tem], [r["stall_1to1"] for r in tem], color=INK, lw=2.6, alpha=0.35, zorder=0)
    axs[1, 2].annotate("default (entropic landscape):\nall members, ∝ T", (6, 7.55), fontsize=7, color=INK)
    axs[0, 0].axhline(0, color=MUTED, lw=0.8)
    axs[0, 0].set(xlabel="hindering load F (pN)", ylabel="velocity (nm/s)")
    axs[0, 0].set_title("a  velocity, 1 mM ATP", fontsize=8.5)
    axs[0, 1].set(xlabel="hindering load F (pN)", ylabel="randomness r = 2D/(v d)", ylim=(0, 5))
    axs[0, 1].set_title("b  randomness, 2 mM ATP (diverges at stall)", fontsize=8.5)
    axs[0, 2].axhline(1, color=MUTED, lw=0.8)
    axs[0, 2].set(xlabel="hindering load F (pN)", ylabel="forward : back step odds", yscale="log")
    axs[0, 2].set_title("c  step odds, 1 mM ATP", fontsize=8.5)
    axs[1, 0].set(xlabel="hindering load F (pN)", ylabel="mismatch kT·D(p‖q) per committed step (kT)", yscale="log")
    axs[1, 0].set_title("d  predicted imprint per step, 1 mM ATP", fontsize=8.5)
    axs[1, 1].set(xlabel="relative viscosity η/η₀", ylabel="v / v(η₀), F = 0, 1 mM ATP", ylim=(0.5, 1.02))
    axs[1, 1].set_title("e  viscosity (search ∝ 1/η, gate η-independent)", fontsize=8.5)
    axs[1, 2].set(xlabel="temperature (°C)", ylabel="1:1 load (pN)")
    axs[1, 2].set_title("f  1:1 load vs temperature\n(coloured: enthalpic-tether variant)", fontsize=8.5)
    legend(fig, plt, y=-0.005)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(FIGURES / "fig6_limits.png")
    plt.close(fig)


def main():
    V, L, A = load()
    fig_valley(V, L, A)
    fig_limits(L)


if __name__ == "__main__":
    main()
