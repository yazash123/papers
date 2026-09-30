# Competing exits

A small, tested Python package for one kinetic idea: a system waits in one
state and leaves by the **first of several independent Poisson exits**. It is
written for kinesin here, but nothing in the core is kinesin-specific; kinetic
proofreading and microtubule search-and-capture fit the same template.

Results and their caveats are in [REPORT.md](REPORT.md) (session 1) and
[REPORT2.md](REPORT2.md) (session 2: is neck-linker docking a "state bet"?).
Session 2's predictions were fixed in [PREDICTIONS.md](PREDICTIONS.md) before
any comparison with the limit datasets. The papers used, and what each was used
for, are in [papers/INDEX.md](papers/INDEX.md).

## The model

Exit *i* has rate k_i(F) = k_i⁰ exp(−F δ_i / kT). Here F > 0 is a hindering
load in pN, δ_i a load distance in nm, and kT is in pN·nm (4.087 at 23 °C,
4.116 at 25 °C). Because the exits are independent exponential clocks:

| quantity | formula |
|---|---|
| splitting probability | P_i = k_i / Σ_j k_j |
| time to leave | exponential, rate K = Σ_j k_j, independent of which exit wins |
| bet (log-odds of a vs b) | ln(k_a/k_b) = ln O₀ − F(δ_a − δ_b)/kT |
| balance point (odds 1:1) | F* = kT ln O₀ / ℓ, with lever ℓ = δ_a − δ_b |
| price of conditioning an attempt on exit i | ln(1/P_i) nats (1/P_i attempts per success) |

Exits come in four kinds (`ExitKind`):

* **PRODUCTIVE**: the wanted outcome (forward step, correct product, capture);
* **WRONG**: the competing outcome (backstep, wrong product);
* **TERMINATING**: ends the run (detachment);
* **RESTART**: abandons the attempt and returns to the start, at a time cost
  (`cost_time`) and/or fuel cost (`fuel`). Examples: ATP leaving before commitment,
  or a proofreading discard.

An attempt may begin with an entry wait (`entry_time`, e.g. waiting for ATP to
bind). A PRODUCTIVE or WRONG exit may carry a post-exit time (`cost_time`, e.g. the
rest of the cycle after a step). With these, `WaitingState` gives in closed form
the mean dwell per step, the velocity (drift per unit attached time), and the
mean steps, displacement and attached time per run at constant load.

An exit can also lead into another waiting state (`to=`, optionally with
probabilities). `Motor` chains states into a small network and solves it at
constant load by first-step analysis. This is used for kinesin's fast events
after a backstep, and later for multi-step proofreading.

Under a position-dependent load (an optical trap) there is no closed form:

* `simulate_runs` runs a direct Gillespie simulation. It is exact because rates
  change only when a step changes the load.
* `trap_statistics` solves the same problem *exactly* as an absorbing Markov
  chain on position × phase (a sparse linear solve). It gives the mean and
  distribution of the load at detachment and the mean attached time. The tests
  check the simulation against it.

## Layout

```
competing_exits/     the package
  rates.py           rate laws: Bell, PiecewiseBell (Kondo eq. 5), Scaled, constant
  model.py           ExitKind, Exit, WaitingState (analytic), Motor (networks)
  loads.py           Trap(stiffness), Clamp(F)
  simulate.py        Gillespie simulation
  lattice.py         exact trap statistics (absorbing Markov chain)
  kinesin.py         kinesin instances as data: head race v2, Kondo KIF5A
  diffusion.py       1D diffusive first passage (Scharfetter-Gummel chain): exact
                     capture probabilities/times, races against Poisson clocks
  headrace.py        diffusive head races: landscapes, capture rates as rate laws,
                     head race v1 (regression) and v3 (docked front search + gate)
  bet.py             the docking "bet": commitment D(q||p0), mismatch D(p||q),
                     the I-projection family
tests/               pytest suite (model identities, simulation vs exact, regression targets)
analysis/            scripts that produce results/*.json and figures/*.png
data/digitized/      every value read off a figure (analysis/digitize_figures.py)
papers/              the PDFs and INDEX.md
```

## Running

```bash
pip install numpy scipy matplotlib pytest     # or: pip install -e .[analysis,test]
python -m pytest -q                           # ~30 s
cd analysis
PYTHONPATH=.. python part_a_v2.py             # Part A refit             (~15 s)
PYTHONPATH=.. python part_c_kondo.py          # Part C, KIF5A sizing     (~4 min)
PYTHONPATH=.. python part_c_motors.py         # Part C(d), other motors  (~2 s)
PYTHONPATH=.. python part_e_temperature.py    # Part E, temperature      (~1 s)
# session 2
PYTHONPATH=.. python digitize_figures.py      # figure readings -> data/digitized (needs pymupdf)
PYTHONPATH=.. python part_a_v3.py             # v3 fit and chi^2 map     (~10 min)
PYTHONPATH=.. python part_b_maxent.py         # valley, MaxEnt members   (~5 min)
PYTHONPATH=.. python part_b_limits.py         # predictions near limits  (~2 min)
PYTHONPATH=.. python part_b_figures.py        # fig5, fig6
PYTHONPATH=.. python make_tables.py > ../results/tables_partB.md && PYTHONPATH=.. python make_tables.py --disc >> ../results/tables_partB.md
PYTHONPATH=.. python part_c_confront.py       # data comparison, TUR     (~1 min)
```

## Example

```python
from competing_exits import Exit, ExitKind, WaitingState, Bell, constant, Motor, Trap, trap_statistics
from competing_exits.kinesin import head_race_v2, kondo_kif5a

s = head_race_v2(atp_uM=1000)
s.balance_point("forward", "back")   # 7.05 pN
s.velocity(3.0)                      # 335 nm/s
s.price("forward", 6.0)              # 2.64 nats per attempt

m = kondo_kif5a(gate=10)             # backstep rate / 10
trap_statistics(m, Trap(0.05))["mean_load_at_termination"]   # 6.62 pN

# a new system: two outcomes, a discard (restart) and no load dependence
proofreader = WaitingState("bound", kT=4.1, exits=(
    Exit("correct", ExitKind.PRODUCTIVE, constant(10.0)),
    Exit("wrong", ExitKind.WRONG, constant(0.1)),
    Exit("discard", ExitKind.RESTART, constant(5.0), cost_time=0.02, fuel=1),
))
proofreader.log_odds("correct", "wrong")      # ln 100
```
