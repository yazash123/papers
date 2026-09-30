"""Figures for Part B (predictions; no data plotted).

fig5_valley.png       the valley of docked tethers, commitment contours, MaxEnt members,
                      and the members' arrival rates
fig6_limits.png       predictions near the limits (load, viscosity, temperature)
Run after part_b_maxent.py and part_b_limits.py:
    cd analysis && PYTHONPATH=.. python part_b_figures.py
"""
from __future__ import annotations

import json

import numpy as np

from common import COLORS, FIGURES, MUTED, plot_style
from competing_exits import bet
from partb_common import ROOT, X

ORDER = ["maxent_k0=0.03", "maxent_k0=0.21", "best", "middle", "floppy_edge"]
LABEL = {"maxent_k0=0.03": "MaxEnt (κ₀ 0.03)", "maxent_k0=0.21": "MaxEnt (κ₀ 0.21)", "best": "best fit",
         "middle": "mid-valley", "floppy_edge": "floppy edge"}
STYLE = {name: dict(color=COLORS[i], lw=2.0) for i, name in enumerate(ORDER)}
INK = "#0b0b0b"


def load():
    V = json.loads((ROOT / "results" / "part_b_valley.json").read_text())
    L = json.loads((ROOT / "results" / "part_b_limits.json").read_text())
    A = json.loads((ROOT / "results" / "part_a_v3.json").read_text())
    return V, L, A


def direct_label(ax, x, y, text, color, dx=0.0, dy=0.0):
    ax.annotate(text, (x, y), xytext=(4 + dx, dy), textcoords="offset points", fontsize=7.5, color=INK,
                va="center")
    ax.plot([x], [y], "o", ms=4, color=color, mec="white", mew=0.8, zorder=5)


def fig_valley(V, L, A):
    plt = plot_style()
    fig, axs = plt.subplots(1, 3, figsize=(12.6, 3.9))
    k = np.array(A["grid"]["kappa"])
    x = np.array(A["grid"]["x_eq"])
    W = V["weights"]["session1"]
    chi = np.array(A["session1"]["chi2_grid"], dtype=float)
    d = chi - W["chi2_min"]
    ax = axs[0]
    XX, KK = np.meshgrid(x, k)
    ax.contourf(XX, KK, np.where(np.isfinite(d), d, 99), levels=[0, 2.30, 5.99], colors=["#bcd4f1", "#dce8f7"])
    cs = ax.contour(XX, KK, np.where(np.isfinite(d), d, 99), levels=[2.30, 5.99], colors=[COLORS[0]], linewidths=[1.0, 1.6])
    for kappa0, ls in ((0.21, "--"), (0.03, ":")):
        D = np.array(W["per_kappa0"][f"{kappa0:g}"]["D_grid"], dtype=float)
        c2 = ax.contour(XX, KK, D, levels=[0.1, 0.3, 1.2], colors=[MUTED], linewidths=0.9, linestyles=ls)
        ax.clabel(c2, fmt=lambda v: f"{v:g}", fontsize=6.5)
    for name in ORDER:
        p = L["members"][name]["params"]
        ax.plot([p["x_eq"]], [p["kappa"]], "o", ms=7, color=STYLE[name]["color"], mec="white", mew=1.0, zorder=6)
        ax.annotate(LABEL[name], (p["x_eq"], p["kappa"]), xytext=(5, 3), textcoords="offset points", fontsize=7, color=INK)
    ax.set_yscale("log")
    ax.set(xlabel="docked tether centre x′ (nm)", ylabel="docked stiffness κ′ (kT/nm²)", xlim=(-2, 8), ylim=(0.04, 0.65),
           title="a  valley (Δχ² ≤ 2.30, 5.99) and D(q‖p₀)\n(dashed κ₀ = 0.21, dotted κ₀ = 0.03; kT)")
    ax.title.set_fontsize(8.5)
    # b: arrival rates k_f(F) of the members
    ax = axs[1]
    for name in ORDER:
        rec = L["members"][name]["load_k0=0.21"]
        F = np.array([r["F"] for r in rec])
        # capture rate = forward exit rate: from P_forward and the race (K); store via odds*kb isn't available,
        # so plot the per-attempt forward probability instead (what the bet delivers)
        ax.plot(F, [r["P_forward"] for r in rec], **STYLE[name])
    ax.set_yscale("log")
    ax.set(xlabel="hindering load F (pN)", ylabel="P(forward capture wins an attempt)", ylim=(1e-4, 1.2),
           title="b  per-attempt forward capture, 1 mM ATP")
    ax.title.set_fontsize(8.5)
    for name in ORDER:
        rec = L["members"][name]["load_k0=0.21"]
        r = rec[-1]
        ax.annotate(LABEL[name], (r["F"], r["P_forward"]), xytext=(3, 0), textcoords="offset points", fontsize=6.5,
                    color=INK, va="center")
    ax.set_xlim(0, 14.5)
    # c: densities
    ax = axs[2]
    p0 = bet.tilted(bet.tether_logdensity(0.21, 0.0), X)
    ax.plot(X, p0, color=MUTED, lw=1.6, ls="--")
    ax.annotate("p₀ (undocked, κ₀ 0.21)", (-6.0, p0[np.searchsorted(X, -6.0)]), xytext=(0, 8), textcoords="offset points",
                fontsize=7, color=INK)
    for name in ORDER:
        p = L["members"][name]["params"]
        q = bet.tilted(bet.tether_logdensity(p["kappa"], p["x_eq"]), X)
        ax.plot(X, q, **STYLE[name])
    ip = V["weights"]["session1"]["iprojection"].get("0.21", {}).get("member")
    if ip and ip.get("lambdas") is not None:
        fs = tuple(0.5 * F / 4.087 for F in (3.0, 9.0))
        q = bet.tilted(bet.iprojection_logdensity(0.21, ip["lambdas"], fs), X)
        ax.plot(X, q, color=COLORS[6], lw=2.0)
        ax.annotate("I-projection member (κ₀ 0.21)", (-7.5, 0.02), fontsize=7, color=INK)
    ax.set(xlabel="free-head position x (nm)", ylabel="density (1/nm)", title="c  where the bet puts the head (F = 0)",
           xlim=(-8, 8))
    ax.title.set_fontsize(8.5)
    handles = [plt.Line2D([], [], **STYLE[n]) for n in ORDER]
    fig.legend(handles, [LABEL[n] for n in ORDER], loc="lower center", ncol=5, fontsize=7.5, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(FIGURES / "fig5_valley.png")


def fig_limits(V, L):
    plt = plot_style()
    fig, axs = plt.subplots(2, 3, figsize=(12.6, 6.8))
    # (a) load: velocity, randomness, mismatch per step
    for name in ORDER:
        rec = L["members"][name]["load_k0=0.21"]
        F = np.array([r["F"] for r in rec])
        axs[0, 0].plot(F, [r["v"] for r in rec], **STYLE[name])
        r_ = np.array([r["randomness"] for r in rec])
        axs[0, 1].plot(F, np.where(np.abs(r_) < 20, r_, np.nan), **STYLE[name])
        axs[0, 2].plot(F, [r["mismatch_per_step"] for r in rec], **STYLE[name])
        rec3 = L["members"][name]["load_k0=0.03"]
        axs[0, 2].plot(F, [r["mismatch_per_step"] for r in rec3], color=STYLE[name]["color"], lw=1.2, ls=":")
    axs[0, 0].axhline(0, color=MUTED, lw=0.8)
    axs[0, 0].set(xlabel="load F (pN)", ylabel="velocity (nm/s)", title="a  velocity, 1 mM ATP")
    axs[0, 1].set(xlabel="load F (pN)", ylabel="randomness r = 2D/(v d)", ylim=(0, 6), title="b  randomness (diverges where v → 0)")
    axs[0, 2].set(xlabel="load F (pN)", ylabel="mismatch cost per step (kT)", yscale="log",
                  title="c  kT·D(p‖q) per step (solid κ₀ 0.21, dotted 0.03)")
    # (b) viscosity
    for name in ORDER:
        rec = L["members"][name]["viscosity"]["F=0.0"]
        eta = L["conditions"]["etas"]
        v0 = rec[0]["v"]
        axs[1, 0].plot(eta, [r["v"] / v0 for r in rec], **STYLE[name])
        recg = L["members"][name]["viscosity_gate_viscous"]["F=0.0"]
        axs[1, 0].plot(eta, [r["v"] / v0 for r in recg], color=STYLE[name]["color"], lw=1.2, ls=":")
    axs[1, 0].set(xlabel="relative viscosity η/η₀", ylabel="v / v(η₀), F = 0", ylim=(0, 1.05),
                  title="d  viscosity (solid: gate η-independent; dotted: gate ∝ 1/η)")
    # (c) temperature: stall and velocity at 5 pN
    for name in ORDER:
        rec = L["members"][name]["temperature"]["default"]["F=5.0"]
        T = [r["T_C"] for r in rec]
        axs[1, 1].plot(T, [r["stall_1to1"] for r in rec], **STYLE[name])
        rec_e = L["members"][name]["temperature"]["enthalpic_tether"]["F=5.0"]
        axs[1, 1].plot(T, [r["stall_1to1"] for r in rec_e], color=STYLE[name]["color"], lw=1.2, ls=":")
        axs[1, 2].plot(T, [r["mismatch_per_s"] for r in rec], **STYLE[name])
    axs[1, 1].set(xlabel="temperature (°C)", ylabel="1:1 load (pN)", title="e  forward:back 1:1 load vs T\n(solid entropic tether, dotted enthalpic)")
    axs[1, 2].set(xlabel="temperature (°C)", ylabel="mismatch cost rate at 5 pN (kT/s)", yscale="log",
                  title="f  mismatch cost per second, 5 pN")
    for ax in axs.ravel():
        ax.title.set_fontsize(8.5)
    handles = [plt.Line2D([], [], **STYLE[n]) for n in ORDER]
    fig.legend(handles, [LABEL[n] for n in ORDER], loc="lower center", ncol=5, fontsize=8, bbox_to_anchor=(0.5, -0.01))
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(FIGURES / "fig6_limits.png")


def main():
    V, L, A = load()
    fig_valley(V, L, A)
    fig_limits(V, L)


if __name__ == "__main__":
    main()
