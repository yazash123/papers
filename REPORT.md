# Competing exits for kinesin: report

**Summary.** I rebuilt head race v2 as a small, general competing-exits package
with 46 tests. Every v2 regression target is reproduced (stall 7.05 pN,
0.95 nats/pN, speeds, dwell, prices). Refitting v2 to Carter & Cross's three
relations pins the odds and the stall, but not the ATP-release clock. The
clock's allowed range depends on how noisy one assumes the binned dwells are:
about 40–80 s⁻¹ at 10% per bin, 28–134 s⁻¹ at 20%. v2's unloaded speed
(454 nm/s) cannot be raised to the measured ~830 nm/s by any refit to the
3–9 pN data. So below ~3 pN the motor is faster than v2's load-independent
"rest of cycle" T and ATP binding allow.

For mouse KIF5A (Kondo et al. 2023) I reproduce the stated checks and the
earlier gate table, to within the earlier calculation's Monte Carlo error. In
the protocol the question was posed in (0.05 pN/nm trap, starting unloaded),
the measured gate sits **at** the shoulder of the gate–grip curve:
- a stronger gate would add at most 0.37 pN (6%) to the mean load at
  detachment;
- the 95% shoulder is at gate factor 1.1 (bootstrap 0.7–2.2, with Table 1's
  errors inflated five-fold).

Two things weaken "sized to the grip":
- Kondo's fast detachments after backsteps move the shoulder to a
  ~2.5× stronger gate (range 1–7.5).
- In stiffer traps or a force clamp the shoulder moves to 2.6–7.5× stronger
  gates.

So the measured gate is "not much stronger than needed" in a soft trap, and
2–7× weaker than the shoulder under stiffer loading. Across motors the claim
that gate strength tracks detachment under load is **not supported**. Kinesin-2
(KIF3) detaches at 2.5–5 pN in the same model trap, yet its fitted and measured
backstep odds reach 1:1 only at 8–13 pN, as late as or later than kinesin-1's.
Only five motors have gate data at all, three of them from one fit.

On the literature (Part D): the stall law, the load levers and the "gate"
reading are all already present in Fisher & Kolomeisky (2001), Liepelt &
Lipowsky (2007) and Carter & Cross (2005) themselves. What is new here is the
sizing question and its (mixed) answer. Part E: the published temperature data
cannot separate the three readings of the odds; §6 names the measurement that
would.

Figures: [fig1](figures/fig1_v2_vs_carter_cross.png) (v2 vs Carter & Cross, clock profile),
[fig2](figures/fig2_gate_vs_grip_kif5a.png) (KIF5A gate vs grip, uncertainty),
[fig3](figures/fig3_gate_vs_grip_motors.png) (all motors),
[fig4](figures/fig4_temperature.png) (temperature).
Numbers behind every table are in `results/*.json`, produced by `analysis/*.py`.

---

## 1. Part A: head race v2 (Drosophila kinesin-1, 23 °C)

### 1.1 The data v2 was fitted to, checked against the PDF

All three relations are in the legend of Fig. 2, p.310 of Carter & Cross 2005
(*Nature* 435:308):

* Ratio: "Inset, ratio of forward to backward steps plotted against the load
  … The fit is: ratio = 802e^−0.95×load, with the load in piconewtons."
* Dwells: "Mean dwell data for forward steps above 3 pN were fitted to single
  exponentials (solid lines). The fits are dwell = 0.0036e^0.57×load at 1 mM ATP
  and dwell = 0.0256e^0.55×load at 10 µM ATP, with the load in piconewtons."
  Dwells are in seconds, so 3.6 and 25.6 ms. The PDF text layer renders "µM" as
  "mM"; the rendered page shows µM.
* Conditions: "All recordings were made at 23 °C" (Methods, p.312), with
  full-length Drosophila kinesin.

Two details matter below:
* Carter & Cross write "At the stall force of 7.2 pN (ref. 17) this ratio is 1"
  (p.309). But their own fit reaches 1 at ln 802/0.95 = **7.04 pN**; the 7.2 pN
  is Nishiyama et al.'s (2002) stall.
* Their step finder compared windows of W = 8 ms either side of each point, and
  excluded steps larger than 12 nm (Methods, p.312). This matters in §1.4.

### 1.2 Regression targets (all reproduced; `tests/test_kinesin_v2.py`)

Parameters as specified: kT = 4.087 pN·nm, d = 8.2 nm, k_f0 = 4182 s⁻¹,
δ_f = 4.27 nm, k_b0 = 5.2 s⁻¹, δ_b = 0.39 nm, k_c = 100 s⁻¹,
k_on = 1.54 µM⁻¹s⁻¹, T = 17.1 ms.

| quantity | target | model |
|---|---|---|
| stall (odds 1:1), any [ATP] | 7.05 pN | 7.047 pN (kT ln(4182/5.2)/3.88) |
| odds fall | 0.95 nats/pN | 0.949 nats/pN |
| speed at 1 mM, F = 0, 3, 5, 6 pN | 454, 335, 104, 32.5 nm/s | 454.3, 334.6, 104.0, 32.5 |
| speed at 10 µM, F = 0 | 98 nm/s | 97.6 |
| dwell per step, 1 mM, F = 0 | 18.0 ms (17.1 is T) | 18.00 ms |
| price ln(1/P_f) at 0, 3, 6 pN | 0.02, 0.45, 2.64 nats | 0.025, 0.452, 2.638 |
| price at stall; attempts per forward step | 3.68 nats; ~40 | 3.681; 39.7 |

The only loose match is the zero-load price. The target 0.02 is given to one
significant figure; the model gives 0.0248, and the test uses an absolute
tolerance of 0.005 there.

### 1.3 Refit and identifiability (`analysis/part_a_v2.py`, fig1)

The paper gives fitted curves, not tables of bin means. So I refitted v2's
seven parameters to the three curves, sampled at 1-pN spacing over 3–9 pN
(21 points). The weights are **my assumption**: σ(ln dwell) = 0.10 per bin
(typical s.e.m. of a bin in Fig. 2b) and σ(ln ratio) = 0.20. Profile
likelihoods refit all other parameters at each value.

**How well v2 matches the curves:**
* Odds: exactly.
* Dwells: 15% RMS (1 mM) and 14% RMS (10 µM). The model cannot be a pure
  exponential in F, so it sags below the 1-mM curve at 9 pN (−33%).
* v2 is not the least-squares optimum under these weights: Δχ² = 12.8 above the
  best fit. It is still inside the joint 95% region for 7 parameters
  (χ²₇ = 14.1).

| parameter | v2 | best refit | 68% (profile) | 95% (profile) |
|---|---|---|---|---|
| k_f0 (s⁻¹) | 4182 | 5063 | 3940–6500 | 2710–9460 |
| δ_f (nm) | 4.27 | 4.41 | 4.28–4.53 | 4.16–4.78 |
| k_b0 (s⁻¹) | 5.2 | 8.1 | 6.6–9.8 | 4.9–14.7 |
| δ_b (nm) | 0.39 | 0.69 | 0.59–0.79 | 0.49–0.99 |
| k_c (s⁻¹) | 100 | 53 | 48–57 | 40–80 |
| k_on (µM⁻¹s⁻¹) | 1.54 | 0.98 | 0.89–1.14 | 0.77–1.40 |
| T (ms) | 17.1 | 15.7 | ~14–18 | ~11–20 |

**What is and is not pinned.**
* The odds pin ln O₀ and the lever ℓ = δ_f − δ_b = 0.95 kT = 3.88 nm, so the
  stall is fixed. Along the whole clock profile the stall stays at
  7.06–7.17 pN.
* The individual δ's are set by the dwells.
* The clock is the loosest parameter, because it trades off against k_on. At
  high load the dwells constrain mainly k_c/k_on and k_f, so k_on rises from
  0.58 to 2.4 µM⁻¹s⁻¹ as k_c goes from 20 to 160 s⁻¹. Its 95% range depends
  directly on the assumed scatter:

| assumed σ per bin | k_c 95% range (s⁻¹) |
|---|---|
| 10% | 40–80 |
| 15% | 34–95 |
| 20% | 28–134 |
| 30% | 20–630 |

The brief's "about 20 to 160 s⁻¹ fits" corresponds to about 20–25% scatter per
bin. The error bars in Fig. 2b look comparable, but I could not quantify them
without reading them off the figure.

**What moves with the clock:** the price at stall rises from 2.4 to 4.2 nats
(11–68 attempts per forward step) as k_c goes from 20 to 160 s⁻¹. So the price
of the timed inference is not pinned by these data. The stall is.

### 1.4 The unloaded-speed miss: what it implies

v2 gives 454 nm/s unloaded at 1 mM and 98 nm/s at 10 µM. The trap-off bead
velocities in Carter & Cross's Fig. 2c (squares) are ≈ 829 and ≈ 142 nm/s.
These are **read off the figure** from pixel positions. Carter & Cross's *own*
loaded bins at 1 mM are ≈ 372 and 341 nm/s at 1.6 and 2.6 pN (also read off),
close to v2's 425 and 370. So v2 agrees with the loaded data; the factor of two
lies between the trap-off beads and the first loaded bins.

What the unloaded speeds imply:

1. **T cannot be a load-independent biochemical time.** At 1 mM the unloaded
   dwell is 8.2/829 = 9.9 ms, less than T = 17.1 ms. After the race and ATP
   wait, at most ~9 ms is left for the rest of the cycle.
2. **ATP binding must be faster at zero load.** At 10 µM the unloaded dwell is
   57.7 ms, shorter than v2's ATP wait alone (65 ms). That requires
   k_on ≥ 1.73 µM⁻¹s⁻¹ at zero load, against v2's 1.54 and the refit's 0.98.
3. **No refit to the 3–9 pN relations removes the miss.** Across the clock
   profile (k_c = 20–160 s⁻¹) the refitted unloaded 1-mM speed stays at
   478–563 nm/s. The load dependence below ~3 pN is simply outside the data v2
   was fitted to.
4. **The knee has precedent.** Kondo's KIF5A forward rate has exactly such a
   knee (flat below 3.2 pN, then falling), and Fisher & Kolomeisky's fit puts
   ~13% of the load dependence on ATP binding.

**A hypothesis, not tested here.** The offset needed to turn the true unloaded
dwell (9.9 ms) into v2's 18.0 ms is 8.1 ms, which equals the 8-ms step-finder
window. If dwells shorter than roughly the window were missed and merged, then
for exponential dwells the surviving mean is inflated by roughly the dead time
(memorylessness). Kondo et al. raise exactly this problem for earlier dwell data:
"dwell times could have been overestimated because mean dwell time was
calculated from data where undetectably short dwell times <4 ms were not
included" (p.463). If so, part of T is a detection offset, not biochemistry.

The stall and the odds would be unaffected, since they come from step counts.
Checking this needs the raw traces, or a simulation of their step finder on
synthetic data.

---

## 2. Part B: the code

`competing_exits/` (see README). The core knows nothing about kinesin:
* rate laws;
* four exit kinds;
* one waiting state with closed forms: splitting, log-odds, lever, balance
  point, price, attempts per success, dwell, velocity and run statistics;
* small networks of states (first-step analysis);
* a Gillespie simulator;
* an **exact** solver for position-dependent loads (absorbing Markov chain on
  lattice × phase).

The kinesin instances live in `kinesin.py` as data.

Analytic checks in the 46 tests (`python -m pytest -q`, ~10 s):
* Identities: splitting sums to 1; log-odds = difference of prices; closed-form
  balance point; lever = δ_a − δ_b; gate scaling moves the balance point by
  kT ln g/ℓ.
* Single-state formulas agree with network first-step analysis. A chained fast
  state agrees with hand-summed geometric series.
* Under a constant-force clamp, the lattice solver reproduces the closed-form
  run statistics.
* Gillespie agrees with the analytic results within 4 standard errors for:
  - velocity and splitting at constant load (v2 with ATP wait, restarts and
    post-step time);
  - the dwell mean;
  - run time and displacement at constant load, with and without fast events;
  - trap statistics against the exact lattice solution, including entry waits,
    restarts and fixed delays.

---

## 3. Part C: is KIF5A's gate sized to its grip? (Kondo et al. 2023)

### 3.1 Rates and checks (`tests/test_kondo.py`)

Table 1 (p.470) was read from the rendered page; the text layer drops the
minus signs. The rate law is k = λ exp(d (L − L_j)/kT). Parameters:
* forward: λ = 77.3 s⁻¹, L_f = 3.19 pN, d = −1.78 nm above the knee and −0.17 nm
  below;
* backstep: λ = 1.99 s⁻¹, d = +0.44 nm;
* detachment: λ = 1.12 s⁻¹, L_d = 1.32 pN, d = +0.45 nm above and −0.78 nm below;
* fast backstep: 2200 s⁻¹, d = +0.11 nm; fast detachment: 1510 s⁻¹, d = +0.13 nm.

Conditions: 25 ± 1 °C (p.473), so kT = 4.116 pN·nm; trap 0.090 pN/nm for
200-nm beads.

| check | stated | reproduced |
|---|---|---|
| unloaded odds | 44:1 | 44.3:1 (3.79 nats) |
| backstep rate unloaded; rise by 6 pN | 1.99 s⁻¹; ~1.9× | 1.99; 1.90× |
| fitted k_f, k_b cross | ~9.3 pN | 9.34 pN |
| unloaded velocity (p.471) | 707 nm/s | 706.8 = d (k_f − k_b) |
| full model zero velocity (p.471) | ~9.1 pN | 9.09 pN, with fast events at q = 0.22 (see 3.3) |

The last check needs one parameter the paper does not give: the probability
q that a backstep enters the fast state-3. I reconstructed it from the text,
"The number fraction of fast to slow backward steps … (~15% at 3–8 pN,
Figure 3B)" (p.470). In the model, fast/slow = qr/(1 − qr) with
r = k_bf/(k_bf + k_df) ≈ 0.59, so q = 0.22. Their supplementary equation S53
is not in the folder, so this is a reconstruction; with the slow events only
the zero-velocity point is 9.34 pN.

The lever is 2.22 nm above the knee (0.61 nm below it), against 3.88 nm for
Drosophila kinesin.

### 3.2 Reproducing the earlier gate table (0.05 pN/nm trap, start unloaded, slow events only)

"Earlier" is the table given in the brief. "Exact" is the lattice solver.
"Gillespie" is 4000 runs (± SE).

| gate | earlier (1:1 / mean load at letting go / time) | exact | Gillespie | SD of load at letting go |
|---|---|---|---|---|
| 100× stronger | 17.9 / 6.7 / 0.52 s | 17.88 / 6.66 / 0.509 | 6.64 ± 0.04 / 0.502 ± 0.006 | 2.6 pN |
| 10× stronger | 13.6 / 6.7 / 0.52 | 13.61 / 6.62 / 0.511 | 6.60 ± 0.04 / 0.504 ± 0.006 | 2.6 |
| as measured | 9.3 / 6.4 / 0.53 | 9.34 / 6.29 / 0.526 | 6.29 ± 0.04 / 0.529 ± 0.007 | 2.4 |
| 3× weaker | 7.3 / 5.8 / 0.56 | 7.31 / 5.71 / 0.556 | 5.70 ± 0.03 / 0.551 ± 0.008 | 2.1 |
| 10× weaker | 5.1 / 4.5 / 0.64 | 5.07 / 4.40 / 0.631 | 4.39 ± 0.03 / 0.638 ± 0.009 | 1.6 |

The 1:1 loads agree exactly, using kT at 25 °C. The earlier mean loads are
0.05–0.1 pN higher and the times ~0.01 s longer, all in the same direction.
With the load's SD of ~2.4 pN, a ~1000-run simulation has an SE of ~0.08 pN,
so the gap is within the earlier calculation's Monte Carlo error. The
same-sign pattern suggests common random numbers across gates. I found no
model difference that explains it: d = 8.0 nm or a softer trap move it the
wrong way. **Reproduced.**

### 3.3 (a) Kondo's fast events

Kondo's scheme (Fig. 4A, p.468) has two fast events. After a slow backstep, the
motor passes through state 2 and enters state-3 with probability q; from
state-3 it takes a fast backstep (then again state-3 with probability q, else
state-0) or detaches fast. Both happen at ~1–2 ms⁻¹, independent of ATP (p.469).
In the model the slow backstep exit leads to a network state `fast`.

Two gaps in the paper:
* **q** is not reported. I use q = 0.22 (from the 15% statement above) and
  q = 0.37 (from the high-load fractions in Fig. 3B: fast back ~0.15, fast
  detach ~0.09, slow back ~0.5; **read off a figure**).
* **The fast rates' load signs.** Table 1 gives d_bf, d_df > 0 (rates rising
  with load), but Fig. 3C draws both fast rates *falling* with load, and
  Fig. 3A's fast dwell rises from ~0.3 to ~1 ms. This is an inconsistency in
  the paper. Both signs were run; the zero-velocity point moves by only
  0.01 pN.

| gate | mean load at letting go: slow only → with fast (q = 0.22) | time attached | zero velocity |
|---|---|---|---|
| 100× | 6.66 → 6.65 pN | 0.509 → 0.508 s | 17.63 pN |
| 10× | 6.62 → 6.57 | 0.511 → 0.502 | 13.36 |
| 1 | 6.29 → 5.90 | 0.526 → 0.450 | 9.09 |
| 1/3 | 5.71 → 4.88 | 0.556 → 0.367 | 7.05 |
| 1/10 | 4.40 → 3.01 | 0.631 → 0.230 | 4.82 |

Every backstep now carries ~10% risk of a fast detachment (q(1−r)/(1−qr)). As
the context document anticipated, this makes weak gates look worse still. At
the measured gate it costs 0.39 pN (0.69 pN at q = 0.37) of mean detachment
load and 15% of the attached time.

### 3.4 Where the shoulder is

Definitions, all computed exactly:
* **S95 / S90:** the weakest gate whose mean load at detachment is 95% / 90% of
  the value for an infinitely strong gate.
* **S×:** the gate whose 1:1 load equals that infinite-gate detachment load.

The position of the measured gate is given as decades above the shoulder,
log₁₀(1/g_shoulder).

**(b) Uncertainty.** Parametric bootstrap, 300 draws per row. Table 1's
parameters are drawn as independent normals from their printed SEs; the
correlations are unknown. Those SEs are fit errors and look small next to the
scatter of Fig. 3C (e.g. λ_b = 1.99 ± 0.01 s⁻¹), so I also inflate them 5×. The
last row also draws q ∈ U(0.15, 0.40) and the fast-rate sign at random.
Values are median [95% range].

| variant | 1:1 load (pN) | mean load at letting go (pN) | g at S95 | decades above S95 | above S90 | above S× |
|---|---|---|---|---|---|---|
| Table 1 SEs | 9.33 [9.11, 9.61] | 6.27 [6.01, 6.55] | 1.10 [0.99, 1.25] | −0.04 [−0.10, 0.01] | 0.29 [0.24, 0.34] | 0.64 [0.57, 0.69] |
| SEs × 5 | 9.27 [8.29, 10.78] | 6.04 [4.74, 7.84] | 1.08 [0.69, 2.23] | −0.03 [−0.35, 0.16] | 0.30 [0.00, 0.49] | 0.69 [0.20, 0.95] |
| SEs × 5 + fast events | 9.29 [8.26, 10.85] | 5.65 [4.21, 6.94] | 2.85 [1.05, 7.46] | −0.46 [−0.87, −0.02] | −0.12 [−0.53, 0.30] | 0.67 [0.21, 1.07] |

**(c) Loading protocol.** Shoulders for the central parameter values. "Gain"
is the rise in mean load at letting go from the measured gate to an
infinitely strong gate.

| protocol | infinite-gate load at letting go | measured gate | gain | g at S95 (slow / fast q = 0.22) |
|---|---|---|---|---|
| trap 0.02 pN/nm | 4.82 pN | 4.66 pN | 0.16 pN (3%) | 0.66 / 2.4 |
| trap 0.05 pN/nm | 6.66 | 6.29 | 0.37 (6%) | 1.1 / 2.5 |
| trap 0.09 pN/nm (Kondo's own) | 7.87 | 7.27 | 0.60 (8%) | 1.6 / 3.0 |
| trap 0.2 pN/nm | 9.62 | 8.56 | 1.06 (11%) | 2.6 / 4.1 |
| trap 0.5 pN/nm | 11.97 | 10.31 | 1.67 (14%) | 3.8 / 5.9 |
| force clamp: largest force with ≥ 1 net step per run | 10.63 | 8.60 | 2.03 (19%) | 6.0 / 7.5 |

A stiffer trap loads the motor faster than it can detach, so it gets closer
to its 1:1 point. There, backsteps rather than detachment limit it, and the
shoulder moves to stronger gates. The "shoulder" is not a property of the motor
alone; it depends on the loading protocol, by a factor ~6 in gate strength over
these protocols.

### 3.5 (d) Other motors: does gate strength track detachment under load?

One protocol for every motor: the 0.05 pN/nm trap, starting unloaded.
* **Grip:** the mean load at detachment for an infinitely strong gate, computed
  from each motor's own gate-free v(F) and k_off(F).
* **Gate:** the 1:1 load, used only where a backstep rate or fraction was
  measured.

kT = 4.1 pN·nm where a paper gives no temperature.

| motor | gate: 1:1 load (pN) | grip (pN) | source of gate / grip |
|---|---|---|---|
| mouse KIF5A | 9.34 | 6.66 | Kondo Table 1 (own rates for both) |
| Drosophila KHC | 7.05 (7.7–10.4 from Block-lab backsteps) | 5.10 (5.97 via v/L) | C&C via v2; alternative: 3% backsteps at 4 pN (Andreasson 2015b p.1169), lever 3.88 or 2.22 nm / Andreasson eLife 3-state v(F) (Table 2) and k_off = 1.11 e^(0.60F/kT) (Fig. 6) |
| KIF3A/B | 9.4 (10.7 from measured 6%) | 3.41 | Andreasson 2015b: zero velocity of fitted cycle (Fig. 3C) / v/L, L0 = 182 nm, δ = 1.7 nm (Fig. 4E) |
| KIF3A/A | 11.7 (13.2 from measured 8%) | 2.51 | same; L0 = 102 nm, δ = 1.6 nm |
| KIF3B/B | 8.3 (10.3 from measured 3%) | 5.19 | same; L0 = 177 nm, δ = 1.3 nm |
| DmK +1 to +6 AA neck-linker inserts | not measured | 3.6–4.5 | Andreasson eLife Tables 1–2 (v/L) |
| Kin1 / KIF3A/A / KIF1A (Gicking) | not measured (k_back = 3 s⁻¹ assumed for all) | 4.7 / 3.7 / 4.8 | Gicking Table 1 model inputs |
| KIF5C CNB, Latch, CNB+Latch (Budaitis) | not measured | 0.6, 0.4, 0.5 (observed, trap stiffness not given; WT 4.8) | Budaitis Fig. 3B, detachment check only |

Checks behind this table:
* The single-state representation reproduces KIF5A's exact grip (6.660 vs 6.660).
* The KIF3 cycle reproduces its paper's unloaded velocities (467–526 nm/s;
  Fig. 2A–B shows ~470–520).
* The KIF3 model backstep fractions at 4 pN (9.9%, 11.4%, 8.3%) are higher than
  measured (6%, 8%, 3%). So the fitted cycle, if anything, *under*states
  kinesin-2's gate.

**Answer: not supported.** Only five motors have any gate data, and three of
them come from one kinetic fit. Among them, the 1:1 loads cluster at 7–13 pN
while grips span 2.5–6.7 pN, with no positive relation.
* Kinesin-2 detaches within a few steps under load: runs of ~25 nm at 4 pN,
  Fig. 4D.
* Yet its backstep fraction at 4 pN (3–8%) is like kinesin-1's (3%, same
  assay, same paper).
* So kinesin-2's gate sits far above its shoulder, not at it. The same holds in
  the model trap, where the 1:1 point is 3–5× its grip. The KIF5A and
  Drosophila KHC ratios are 1.4.

Other motor data sets only add points on the grip axis:
* The **neck-linker inserts** lose ~20–35% of their grip. Andreasson et al.
  report "no processive backstepping" under superstall loads for all
  non-cysteine-light constructs (p.13) and an intact stepping gate up to +3 AA
  (Fig. 2). That is qualitative, and consistent with a gate that did not
  weaken along with the grip, the opposite of "tracking".
* Andreasson's **cysteine-light caveat** (pp.11–15): CL constructs are slower
  and more load-sensitive under hindering load than wild type. HsK-CL-6AA
  backsteps processively where DmK-6AA does not. So load data from CL motors
  (e.g. earlier NL-extension backstepping) should not be used as gate data;
  none were used here.
* The **Budaitis** docking mutants let go at ~0.5 pN, far below any gate-limited
  load. That is a detachment check: detachment alone limits them.

### 3.6 (e) Verdict: does the measured gate sit at the shoulder?

**In the protocol the question was posed in (0.05 pN/nm trap, starting
unloaded, Kondo's slow events), yes, to within its uncertainty.**
* The measured gate is 0.04 decades below the 95% shoulder (95% range −0.10 to
  +0.01 with Table 1's SEs; −0.35 to +0.16 with SEs ×5) and 0.3 decades above
  the 90% shoulder.
* A stronger gate buys at most 0.37 pN (6%) of mean load at detachment; a
  10× weaker gate loses 1.9 pN.

Adding Kondo's fast detachments (with a reconstructed q) moves the measured
gate to 0.46 decades *below* the 95% shoulder [0.87 below to 0.02 below]. So
KIF5A would gain ~13% (0.8 pN) from a gate ~3× stronger.

In stiffer traps and in a force clamp the measured gate is 2.6–7.5× weaker
than the 95% shoulder.

The plain statement the numbers support: KIF5A's gate is **not oversized**. It
is at the shoulder, or up to a factor of a few below it, depending on the
fast events and the loading protocol. Without a cost of gate strength, "right
size" can mean only "not much stronger than a soft-trap grip requires". Across
motors, gate strength does **not** track detachment under load (kinesin-2).

---

## 4. Part D: where this sits in the literature

**Fisher & Kolomeisky 2001** (*PNAS* 98:7748) already contain the stall law and
the levers:
* They put load into each rate through load-distribution factors,
  u_j(F) = u_j⁰e^(−θ_j⁺Fd/kT) and w_j(F) = w_j⁰e^(+θ_j⁻Fd/kT) (eq. 3). These
  *are* this model's δ_i = θ_i d.
* Their stall force F_S = (kT/d) ln Π_j(u_j⁰/w_j⁰) (eq. 6) is the balance-point
  formula. For one state with one forward and one backward exit it reads
  F_S = kT ln(u⁰/w⁰)/ℓ.
* They also add detachment rates to get run lengths (eqs. 16–17).

What they do not allow is ℓ ≠ d. They assume Σ(θ⁺ + θ⁻) = 1, i.e. a
backstep is the microscopic reverse of the forward cycle. Applied to kinesin's
odds, that fixed lever would make the odds fall at 8.2/kT ≈ 2.0 nats/pN. The
measured slopes are 0.95 (Carter & Cross) and 0.54 nats/pN above the knee
(Kondo): levers of 3.9 and 2.2 nm, i.e. 0.47 d and 0.27 d. That the odds fall
at half or less of the reversal slope is the quantitative form of Carter &
Cross's point that backsteps are ATP-triggered, not reversals.

**Liepelt & Lipowsky 2007** (*PRL* 98:258102) go further:
* Forward and backward steps are two distinct ATP-hydrolysing chemomechanical
  cycles, F and B, and the stall is where their fluxes balance (their Fig. 2b
  fits Carter & Cross's ratio). That *is* "stall is where the forward:back odds
  reach 1:1, with both exits consuming ATP".
* They also include unbinding.
* But they route both steps through one mechanical transition whose forward and
  reverse rates split the full step's work, Φ₂₅ = e^(−θF̃), Φ₅₂ = e^((1−θ)F̃),
  with θ = 0.65 for Carter & Cross's data. So their stepping transition
  backwards *speeds up* strongly with load, as e^(0.35·8.0·F/kT).
* That is not the load-insensitive gate. Kondo later measured the gate directly:
  backsteps rise only 1.9× by 6 pN.

**The "gate" reading is not new either.** Carter & Cross themselves write that
"the probability of a forward step decreases exponentially with increasing
load, whereas the probability of a backward step is constant" (pp.311–312).
Kondo et al. fit exactly this competing-exits form.

**Bluntly:**
* The stall law, the levers and the gate reading are all in the literature. The
  stall law is FK's eq. 6 without the ℓ = d constraint, and is implicit in C&C's
  fit and in LL's F/B flux balance.
* The ln(1/P) pricing is a bookkeeping identity.
* What is new here is (i) posing the gate's *size* against a terminating exit
  under a realistic loading protocol and (ii) the answer. KIF5A's gate is at, or
  a factor of a few below, the soft-trap shoulder. The shoulder moves with the
  protocol. Across families the gate does not track the grip.
* A short paper could honestly claim a quantitative sizing test with a mixed,
  partly negative result, plus the code. It could not claim a new stall law or
  a new reading of the backstep.

---

## 5. Reproduced and not reproduced

Reproduced:
* All Part A targets.
* The Carter & Cross relations, verified in the PDF.
* Kondo's odds (44:1), backstep rate and its 1.9× rise, the 9.3-pN crossing and
  707 nm/s.
* Their 9.1-pN zero velocity, but only with a reconstructed q.
* The earlier gate table, within the earlier Monte Carlo error.

Not reproduced, or not reproducible from what is in the folder:
* v2's unloaded 454 nm/s vs the measured ~830 nm/s. This is a model miss below
  3 pN (§1.4).
* v2 is not the least-squares optimum of my (assumed-weight) refit, although it
  lies in its 95% region.
* Kondo's supplementary equations (S16–S18, S53) are not in the folder, so the
  full-model check uses a reconstructed q.
* Andreasson 2015b's kinesin-1 fit (Table S1) and KIF3 detachment tables
  (S3–S4) are in a supplement that is not in the folder.

---

## 6. Part E: temperature (see end of file)

---

## 7. Every guess and every value read off a figure

| item | where | what was done |
|---|---|---|
| trap-off velocities 829 / 142 nm/s; loaded bins 372, 341, 210, 113 nm/s | C&C Fig. 2c | pixel positions on the rendered figure |
| σ(ln dwell) = 0.10, σ(ln ratio) = 0.20 per 1-pN bin | Part A refit | assumed; k_c interval reported for 10–30% |
| fast-state branching q = 0.22 | Kondo p.470 ("~15%") | reconstructed via qr/(1−qr) = 0.15 |
| q = 0.37 | Kondo Fig. 3B high-load fractions | read off a figure |
| sign of fast-rate load dependence | Kondo Table 1 vs Fig. 3C | both signs run; negligible effect |
| Table 1 parameter correlations | Kondo | unknown; drawn independently; SEs also inflated ×5 |
| kT = 4.1 pN·nm | Andreasson (both), Gicking (21 °C stated there) | no temperature stated in the Andreasson papers |
| Drosophila KHC gate from a full-length motor (C&C) combined with a grip from a truncated one (Andreasson) | §3.5 | different constructs and labs; alternative gate from Block-lab backsteps shown |
| lever used to turn a 4-pN backstep fraction into a 1:1 load | §3.5 | fitted levers (v2 3.88 nm, Kondo 2.22 nm, KIF3 per-head δ) |
| KIF3 run lengths at loads < 1 pN | Andreasson 2015b Fig. 4E | loaded-regime exponential used from 0 pN; the unloaded runs are longer (their "two regimes") |
| Gicking's F_s = 6 pN, k_back = 3 s⁻¹, linear force–velocity | Gicking Table 1 | model inputs, not measurements; grip only |
| Budaitis trap stiffness | Budaitis | not given; observed forces used only as a qualitative check |

## 8. Open issues

* **The loaded–unloaded gap in Carter & Cross.** Test the step-finder
  dead-time hypothesis (§1.4) on synthetic traces. If confirmed, T shrinks and
  the unloaded speed is recovered without changing the stall.
* **Kondo's q and supplementary equations.** The full-model and fast-event
  results depend on q; the paper's own S-equations would fix it.
* **Parameter correlations.** Kondo's Table 1 errors are fit SEs without
  covariances. Refitting their Fig. 3C points would give honest intervals
  (this needs digitising).
* **No cost of gate strength has been identified.** Without one, "sized to the
  grip" can only mean "not too weak". The fast-event and stiff-trap results
  make even that protocol-dependent.
* **Gate data for more motors.** Backstep rates vs load exist in the folder for
  only KIF5A and the KIF3 family (one fit). The "tracking" claim needs
  backsteps measured under load for NL mutants (non-CL) and kinesin-3.
* **Detachment in v2.** v2 has no terminating exit. The Drosophila point
  combines a full-length gate with a truncated-construct grip.

## 9. Changes after independent review

(See end of file.)
