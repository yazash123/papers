"""Kinesin instances of the competing-exits model (data, not core code).

Every number here is tied to a source; see papers/INDEX.md and REPORT.md.
"""
from __future__ import annotations

from .model import Exit, ExitKind, Motor, WaitingState
from .rates import Bell, PiecewiseBell, constant, kT_at

P, W, T, R = ExitKind.PRODUCTIVE, ExitKind.WRONG, ExitKind.TERMINATING, ExitKind.RESTART

# ---------------------------------------------------------------------------
# Head race v2: Drosophila kinesin-1, 23 C, fitted to Carter & Cross (2005).
# Parameters as specified for this session (not re-derived here; see
# analysis/part_a_v2.py for the refit).
# ---------------------------------------------------------------------------
V2 = dict(
    kT=4.087,      # pN nm (23 C)
    d=8.2,         # nm, step size
    kf0=4182.0,    # 1/s, forward capture at zero load
    delta_f=4.27,  # nm
    kb0=5.2,       # 1/s, backstep (gate) at zero load
    delta_b=0.39,  # nm
    kc=100.0,      # 1/s, clock: ATP leaves (restart), load-independent
    kon=1.54,      # 1/(uM s), ATP binding
    T=17.1e-3,     # s, rest of the cycle after a step (load/ATP independent)
)


def head_race_v2(atp_uM: float = 1000.0, **overrides) -> WaitingState:
    """One bound ATP buys one attempt, ended by forward capture, backstep or
    the clock (ATP leaves; the head waits for a new ATP)."""
    p = {**V2, **overrides}
    d = p["d"]
    return WaitingState(
        name="atp_bound",
        kT=p["kT"],
        entry_time=1.0 / (p["kon"] * atp_uM),
        exits=(
            Exit("forward", P, Bell(p["kf0"], p["delta_f"]), step=+d, cost_time=p["T"], fuel=1),
            Exit("back", W, Bell(p["kb0"], p["delta_b"]), step=-d, cost_time=p["T"], fuel=1),
            Exit("clock", R, constant(p["kc"])),
        ),
    )


# ---------------------------------------------------------------------------
# Mouse KIF5A, Kondo, Sasaki & Higuchi (2023) Traffic 24:463, Table 1 (p.470),
# 25 +/- 1 C (p.473).  Their eq. 5: k = lambda exp(d (L - L_j)/kT); here
# delta = -d.  Values are mean and SE as printed.
# ---------------------------------------------------------------------------
KONDO = dict(
    kT=kT_at(25.0),
    d=8.2,                 # nm (their measured mean 8.17 +/- 0.15, p.464)
    lam_f=(77.3, 0.9), L_f=(3.19, 0.05), d_f_plus=(-1.78, 0.04), d_f_minus=(-0.17, 0.01),
    lam_b=(1.99, 0.01), d_b=(0.44, 0.01),
    lam_d=(1.12, 0.07), L_d=(1.32, 0.55), d_d_plus=(0.45, 0.01), d_d_minus=(-0.78, 0.05),
    lam_bf=(2200.0, 800.0), d_bf=(0.11, 0.04),
    lam_df=(1510.0, 690.0), d_df=(0.13, 0.06),
)


def kondo_values(overrides=None) -> dict:
    """Central values of Table 1 (optionally overridden, e.g. by a bootstrap draw)."""
    v = {k: (val[0] if isinstance(val, tuple) else val) for k, val in KONDO.items()}
    v.update(overrides or {})
    return v


def kondo_kif5a(gate: float = 1.0, fast: bool = False, q: float = 0.22,
                fast_sign: float = +1.0, values: dict | None = None) -> Motor:
    """KIF5A as one waiting state (state-0 of Kondo et al.) with forward step,
    slow backstep and slow detachment.

    gate      the backstep rate is divided by ``gate`` (gate > 1 is stronger)
    fast      add Kondo's fast events: after a slow backstep the motor enters
              a transient state '3' with probability q, from which it takes a
              fast backstep (and again enters '3' with probability q) or
              detaches fast.  q is not given in the paper (see REPORT.md).
    fast_sign +1 uses Table 1's signs for d_bf, d_df (fast rates rise with
              load); -1 flips them, as Fig. 3C appears to show.
    """
    v = values or kondo_values()
    kT, d = v["kT"], v["d"]
    forward = PiecewiseBell(v["lam_f"], v["L_f"], delta_below=-v["d_f_minus"], delta_above=-v["d_f_plus"])
    back = Bell(v["lam_b"] / gate, -v["d_b"])
    detach = PiecewiseBell(v["lam_d"], v["L_d"], delta_below=-v["d_d_minus"], delta_above=-v["d_d_plus"])
    to_back = {"fast": q, None: 1.0 - q} if fast else None
    s0 = WaitingState("0", kT=kT, exits=(
        Exit("forward", P, forward, step=+d, fuel=1),
        Exit("back", W, back, step=-d, fuel=1, to=to_back),
        Exit("detach", T, detach),
    ))
    if not fast:
        return Motor.single(s0, step_size=d)
    s3 = WaitingState("fast", kT=kT, exits=(
        Exit("fast_back", W, Bell(v["lam_bf"], -fast_sign * v["d_bf"]), step=-d, to={"fast": q, None: 1.0 - q}),
        Exit("fast_detach", T, Bell(v["lam_df"], -fast_sign * v["d_df"])),
    ))
    return Motor((s0, s3), home="0", step_size=d)
