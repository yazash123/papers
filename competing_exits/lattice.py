"""Exact run statistics under a position-dependent load (e.g. an optical trap).

The motor moves on a lattice x = n * d.  Every phase of the motor -- the entry
wait, each waiting state, each post-exit delay -- at every lattice site is a
transient state of an absorbing Markov chain; termination is absorption.  The
expected number of visits v to each transient state solves

    (I - P^T) v = e_start,

where P is the jump-chain transition matrix.  From v follow, exactly, the
distribution of the position (hence load) at termination, the mean attached
time (visits times mean holding times) and the mean number of each exit.

Only the means of the entry and delay times matter for these quantities, so
exponential and fixed delays give identical results here.

The lattice is truncated at |n| <= n_max; ``absorbed`` (the total termination
probability) should be 1 to within rounding, and is returned so this can be
checked.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from .model import ExitKind, Motor


def trap_statistics(motor: Motor, load, n_max: int = 300, x0: float = 0.0) -> dict:
    d = motor.step_size
    for s in motor.states:
        for e in s.exits:
            if e.kind != ExitKind.TERMINATING and not np.isclose(abs(e.step), d) and e.step != 0:
                raise ValueError(f"exit {e.name!r}: steps must be 0 or +/- step_size")
    sites = np.arange(-n_max, n_max + 1)
    n_sites = len(sites)
    home = motor.home_state
    has_entry = home.entry_time > 0

    # phase list per site: ('entry',), ('state', s), ('cost', s, e)
    phases = []
    if has_entry:
        phases.append(("entry",))
    for s in motor.states:
        phases.append(("state", s.name))
    for s in motor.states:
        for e in s.exits:
            if e.kind != ExitKind.TERMINATING and e.cost_time > 0:
                phases.append(("cost", s.name, e.name))
    ph_index = {p: i for i, p in enumerate(phases)}
    n_ph = len(phases)
    N = n_sites * n_ph

    def idx(site_i, phase):
        return site_i * n_ph + ph_index[phase]

    def start_phase():
        return ("entry",) if has_entry else ("state", motor.home)

    def arrive(site_i, dest):
        """Index of the transient state entered when going to ``dest``."""
        if dest is None:
            return idx(site_i, start_phase())
        return idx(site_i, ("state", dest))

    rows, cols, vals = [], [], []
    hold = np.zeros(N)
    p_term = np.zeros(N)
    exit_prob = {}  # (state, exit) -> array over transient states
    for si, n in enumerate(sites):
        F = float(load(n * d))
        if has_entry:
            i = idx(si, ("entry",))
            hold[i] = home.entry_time
            rows.append(i); cols.append(idx(si, ("state", motor.home))); vals.append(1.0)
        for s in motor.states:
            i = idx(si, ("state", s.name))
            k = np.array([float(e.k(F, s.kT)) for e in s.exits])
            K = k.sum()
            hold[i] = 1.0 / K
            for e, ke in zip(s.exits, k):
                p = ke / K
                exit_prob.setdefault((s.name, e.name), np.zeros(N))[i] = p
                if e.kind == ExitKind.TERMINATING:
                    p_term[i] += p
                    continue
                sj = si + int(round(e.step / d))
                if sj < 0 or sj >= n_sites:
                    continue  # leaves the truncated lattice (checked via 'absorbed')
                if e.cost_time > 0:
                    rows.append(i); cols.append(idx(sj, ("cost", s.name, e.name))); vals.append(p)
                else:
                    for dest, q in e.destinations().items():
                        rows.append(i); cols.append(arrive(sj, dest)); vals.append(p * q)
        for s in motor.states:
            for e in s.exits:
                if e.kind == ExitKind.TERMINATING or e.cost_time <= 0:
                    continue
                i = idx(si, ("cost", s.name, e.name))
                hold[i] = e.cost_time
                for dest, q in e.destinations().items():
                    rows.append(i); cols.append(arrive(si, dest)); vals.append(q)

    P = sp.csr_matrix((vals, (rows, cols)), shape=(N, N))
    A = (sp.identity(N, format="csr") - P.T).tocsc()
    start = np.zeros(N)
    start[idx(int(np.searchsorted(sites, round(x0 / d))), start_phase())] = 1.0
    v = spla.spsolve(A, start)

    absorbed_at = v * p_term                       # termination probability per transient state
    per_site = absorbed_at.reshape(n_sites, n_ph).sum(axis=1)
    loads = np.array([float(load(n * d)) for n in sites])
    absorbed = per_site.sum()
    mean_load = (per_site * loads).sum() / absorbed
    return {
        "absorbed": absorbed,
        "mean_load_at_termination": mean_load,
        "sd_load_at_termination": np.sqrt(max(0.0, (per_site * loads**2).sum() / absorbed - mean_load**2)),
        "mean_time": float(v @ hold),
        "site_loads": loads,
        "termination_distribution": per_site,
        "mean_counts": {f"{s}.{e}": float(v @ p) for (s, e), p in exit_prob.items()},
    }
