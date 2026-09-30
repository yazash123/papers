# Competing exits, session 2: is neck-linker docking a state bet?

**Plain summary.**

I tested whether kinesin-1 neck-linker docking behaves as a "state bet". In
that picture the docked landscape is the least-committed change from the
undocked tether that still reproduces the measured arrival times. The bet's
cost would appear as a mismatch dissipation, kT·D(p‖q), inside the motor's
hidden dissipation, and would show most clearly near the motor's limits.

I built the model from a diffusive first-passage solver. That solver
reproduces the earlier "v1" model's targets, and a new v3 model with it fits
Carter & Cross (C&C) as well as session 1 did. I then fixed the predictions in
PREDICTIONS.md and committed that file before comparing them with the limit
datasets. Part A had already fitted C&C and checked Kondo, and an unrun
comparison-script skeleton existed untracked (Part B).

The outcome is negative, and it is specific:
* **Within v3, the data need only a small bet.**
  - C&C's arrival times require at most 0.4 kT of commitment, and at most
    0.17 kT in the two-parameter I-projection family.
  - These are upper bounds at a nominal 95% level; at 68%, up to 0.47 kT.
  - They need no bet at all if the undocked neck linker is as stiff as a
    worm-like chain and C&C's own scatter is used, though marginally
    (Δχ² 5.05 against 5.99).
  - With the exact restriction cost (corrected after review), the ~1.2 kT
    docking budget leaves the MaxEnt members affordable. It excludes most of
    the valley if the undocked tether is that stiff.
* **The bet's imprint is too small to isolate.**
  - It is 0.01–0.17 kT per step at 2 pN for the MaxEnt members, about 1% of
    ≈16 kT of hidden dissipation. The 16 kT is my 0.8 × 20.5 kT, from
    second-hand statements of Ariga et al.'s 80%.
  - Under assisting load it reaches ≈1 kT for one MaxEnt member.
  - An Ariga-type Harada–Sasa measurement cannot separate the imprint from
    the chemistry. In Takaki's fitted cycle (my reconstruction from their
    code), the chemical dissipation alone changes by 5.1 kT between 2 and
    6 pN, 11–320 times the predicted contrast.
  - Resolving 0.1 kT would also need 10⁶–10⁷ steps.
* **Near the limits, the landscapes the data allow cannot be told apart**,
  except in their unloaded viscosity response. The bet is therefore not
  visible in the character of the curves.
* **The limit data contradict v3 at several points:**
  - past stall the motor walks backward 1.4–5.7 times faster than predicted,
    with a flat ≈ 300-ms dwell where the model's grows to 0.7–1.8 s;
  - backward-step dwells are 1.3–4 times longer than forward-step dwells at
    intermediate loads;
  - small crowders slow the unloaded motor far more than any member
    predicts: 3–13 times the predicted slowdown for the MaxEnt members, and
    1.3–2.6 times for the floppiest. This is a different construct,
    truncated GFP-kinesin.
  - In Taniguchi's bovine-brain data the balance load does not rise with
    temperature (1.8–2.5 s.e.m.).

  C&C's own 10 µM odds lie above the 1 mM odds near stall, which v3 cannot
  produce. But that selection was made after seeing the data, and the
  authors read their data as ATP-independent.

  Every member fails alike, so these results falsify v3's specific choices:
  - load-independent ATP binding;
  - a backstep gate confined to the ATP-bound race and slowing under load;
  - the search as the only viscous step;
  - an entropic temperature scaling.

  They do not falsify a search-plus-gate scaffold in general. And since the
  bet is sized inside v3, they neither convict nor clear the bet.

The thermodynamic uncertainty relation (TUR), applied to measured randomness,
bounds the total dissipation. The bound captures a quarter to a third of the
total at low load and 5–13% at Visscher's highest load (depending on the
assumed odds and Δμ). It goes slack toward stall and says nothing about the
imprint.

An independent review found errors. They are corrected here and listed under
"Review fixes" at the end. PREDICTIONS.md is left as committed, and its
errors are listed in Part B.

Deliverables: [PREDICTIONS.md](PREDICTIONS.md) (commit `84a3d19`, before any
comparison with the limit datasets); figures 5–7 in `figures/`; tables in
`results/tables_partB.md` and `results/tables_partC.md`.

---

## Part 0. What exists

### 0.1 What is and is not in the folder

**Not in `papers/`**, and not reachable from this session (arxiv.org, PubMed
Central, journals.aps.org, pnas.org and similar hosts are blocked by the
environment's network policy):
* Ariga, Tomishige & Mizuno 2018, *PRL* 121:218101 (arXiv:1704.05302);
* Takaki, Mugnai & Thirumalai 2022, *PNAS* 119:e2208083119;
* Rice et al. 2003, *Biophys. J.* 84:1844 (neck-linker docking thermodynamics);
* Kawaguchi & Ishiwata 2000 (temperature series; noted already in session 1);
* Clancy et al. 2011, *NSMB* 18:1020 (source of the 600 s⁻¹ / 100 s⁻¹ rates
  used by Sozanski et al. and by head race v1);
* the supplements of Carter & Cross 2005 (Fig. S2, dwell histograms) and
  Guydosh & Block 2009 (Suppl. Fig. 7, 0.4-pN load).

There is also no `CONTEXT.md` in the repository on either branch, although the
brief asks for it; session 1's work was on branch
`claude/adoring-goldberg-17q637` and was fast-forwarded into this branch.

So the three quantities the brief calls "exact" (Ariga's hidden dissipation,
Takaki's split, Rice's docking free energy) can only be recorded **second-hand**
below, from papers in the folder that cite them, and, for Takaki, from the
authors' deposited code. Each is labelled with its source.

### 0.2 Ariga, Tomishige & Mizuno 2018 (hidden dissipation)

What the folder says, verbatim:

| source (in folder) | statement |
|---|---|
| Hwang & Karplus 2019, *PNAS*, p.7 of PDF | "For kinesin, under a 2-pN hindering load, 80% of the input energy by ATP was dissipated in sources other than work against the load or viscous drag (126). The authors concluded this 'hidden' dissipation to be internal, which should be due to the conformational fluctuation of the motor head itself." |
| Xie et al. 2019, *IJMS* 20:4911, p.19 of PDF | model efficiency "of about 16% for HsK is in agreement with the available single-molecule value determined recently at F = 2 pN [62]" (ref. 62 = Ariga 2018) |
| Brown & Sivak 2020, *Chem. Rev.* 120:434, p.437 | "under physiologically relevant conditions, kinesin dissipates most of its input free energy.¹⁰¹" |
| S.-W. Wang 2018 (arXiv:1710.10531), p.4 | Ariga et al. set the temperature from the high-frequency ratio C̃/2R̃′ "assuming that FDT is satisfied there"; with perturbation asymmetry the FDT is violated even at high frequency, which biases this calibration. |

Search-engine summaries of the abstract (not verifiable here, flagged): the
sum of probe dissipation and work "did not amount to the input free energy
change"; "kinesin loses 80% of input energy to heat"; the hidden part was
attributed to internal dissipation using a Langevin model of the probe
coupled to a two-state Markov stepper.

**What can be confirmed.**
* **Method.** Harada–Sasa: the probe's velocity correlation and response
  spectra, measured with optical tweezers, give the heat dissipated through
  the probe's degree of freedom. This is by construction an **ensemble-averaged,
  steady-state rate**; a per-step number is that rate divided by the mean
  stepping rate. It is not a per-step measurement.
* **Condition.** A **single** hindering load of **2 pN** (two independent
  folder sources). The ATP concentration, temperature, construct and probe
  details are **not recoverable** from the folder.
* **Value.** ~80% of the ATP free energy is "hidden". The absolute value and
  its uncertainty are not given in any folder source. For use below I take
  hidden ≈ 0.8 × Δμ with Δμ ≈ 20.5 k_BT (the reference line in Takaki's
  code, §0.3), i.e. **≈ 16 k_BT per step, uncertainty unknown** (a range of
  15–17 k_BT follows from Δμ = 19–21 k_BT alone).
* **A consistency check (mine).** At 2 pN the work per 8.2-nm step is
  16.4 pN·nm = 4.0 k_BT ≈ 20% of Δμ. So "80% hidden" means that almost all of
  the non-work free energy escaped the probe: the probe-visible dissipation
  was a small fraction of Δμ.

So the brief's belief is confirmed as far as the folder allows: one condition,
inferred from ensemble-averaged steady-state spectra. It is not confirmed
against the paper itself.

### 0.3 Takaki, Mugnai & Thirumalai 2022 (split of ~20.5 k_BT)

The paper is not available. Its deposited code
(`github.com/kibidanngo/information-flow-pnas`, `Kin1_PNAS.nb`) was cloned
read-only. It contains the fitted 5-state kinesin cycle (`param3`: k12 = 384.4,
k23 = 371.2, k34 = 459.5, k45 = 29.8 µM⁻¹s⁻¹·[ATP], k51 = 3×10⁵; reverse
k21 = 1.5·[ADP], k32 = 33.2, k43 = 2.7·[Pi], k54 = 94.7, k15 = 0.7 s⁻¹;
χ₁₂ = 0.4, θ = 1, d = 8 nm, k_BT = 4.114 pN·nm), fitted to force–velocity
data at 1 mM ATP/100 µM ADP/1 mM Pi and 10 µM/1 µM/1 mM. It defines
Q_chem = ln Π(k₊/k₋) over the chemical transitions 1→2→3→4→5 and
Q_mec = ln(k51/k15) for the mechanical step, and plots 20.5 k_BT as the input
free energy. The chemical and mechanical entropy production per cycle are
σ_chem = Q_chem + ln(p₁/p₅) and σ_mec = Q_mec + ln(p₅/p₁).

My evaluation of that code (a reconstruction, **not a number printed in the
paper**), at 1 mM ATP:

| load | Δμ (model) | σ_chem | σ_mec | work | v (model) |
|---|---|---|---|---|---|
| 0 pN | 20.3 k_BT | 14.4 | 6.0 | 0 | 884 nm/s |
| 2 pN | 20.3 | 11.3 | 5.1 | 3.9 | 568 nm/s |
| 4 pN | 20.3 | 8.7 | 3.8 | 7.8 | 177 nm/s |

A search-engine summary of the abstract (unverified) says "only a fraction of
the energy from ATP hydrolysis is used to advance kinesin motors against load,
with the rest associated with chemical transitions in the two heads", which is
consistent with σ_chem > σ_mec above.

### 0.4 The docking free energy (Rice et al. 2003)

Not in the folder. Second-hand values:

| source | statement |
|---|---|
| Block 2007, *Biophys. J.* 92:2986, p.2991 | Rice et al. "estimated the free energy associated with neck-linker docking, and found it to be just ~3 kJ/mol (note that k_BT is 2.6 kJ/mol)", i.e. ≈1.15 k_BT; Hackney argued the value is an underestimate because Rice used AMP-PNP |
| Xu et al. 2021, *J. Phys. Chem. B* 125:2627, p.2627 | docking "only has a free energy gain of ∼1.2 k_BT" (ref. 4 = Rice 2003) |
| Taniguchi et al. 2005, *Nat. Chem. Biol.* 1:342, p.342 and p.345 | "only 1–2 k_BT" |
| Sumi 2017, *Sci. Rep.* 7:1163, p.2 | "about 3 kJ/mol" |

I use **1.2 k_BT** as the docking budget, with 1–2 k_BT as its range.

### 0.5 Datasets that approach a limit

Every value marked *(fig)* was read off a figure by `analysis/digitize_figures.py`
(vector extraction for Carter & Cross, pixel digitisation otherwise); the CSVs
are in `data/digitized/`.

**(a) Load: force–velocity near and past stall**

| dataset | motor, conditions | how close to the limit |
|---|---|---|
| Visscher et al. 1999, *Nature* 400:184, Fig. 3a (p.186) | squid kinesin, force clamp, 2 mM and 5 µM ATP | 2 mM: 854 → 148 nm/s at 6.3 pN *(fig)*; "stall force of ~7.0 pN" (p.186). 5 µM: 56 → 5.7 nm/s at 5.6 pN *(fig)*; stall ~5.5 pN. Stops short of stall. |
| Carter & Cross 2005, *Nature* 435:308, Fig. 2c (p.310) | full-length *Drosophila*, 23 °C, 1 mM and 10 µM | **through and past stall**: 1 mM from 824 nm/s (trap off) to −31 nm/s at 12.6–13.5 pN; 10 µM from 142 to −22 nm/s *(fig, vector)* |
| Block et al. 2003, *PNAS* 100:2351, Fig. 4A (p.2354) | 1.6 mM, 4.2 µM ATP; 2-D force clamp | 1.6 mM: ~650–700 → 128 nm/s at 4.7 pN hindering *(fig)*; "near zero between −4 and −6 pN" (p.2354), in Block's sign convention (backward, i.e. hindering, loads negative): 4–6 pN hindering in this report's convention |
| Nishiyama et al. 2002, *NCB* 4:790, Table 1 (p.790) | bovine brain kinesin, 298 K | v = 930, 230, 0 nm/s at 0, 3.8, 7.6 pN (model rates from dwells) |
| Kondo et al. 2023, *Traffic* 24:463 | mouse KIF5A, 25 °C | zero velocity ~9.1 pN; loads to >40 pN (session 1) |
| Khalil et al. 2008, *PNAS* 105:19247, Table 1 (p.19250) | *Drosophila* K401 WT and cover-strand mutants | stall forces 4.96 ± 0.05 (WT), 3.02 ± 0.03 (2G), 1.37 ± 0.04 pN (DEL, biased high) |

**(b) Randomness against load and ATP**

| dataset | conditions | how close to the limit |
|---|---|---|
| Visscher 1999, Fig. 4b (p.187) | 2 mM ATP | r = 0.41–0.48 from 0.5 to 4.2 pN; **0.68 ± 0.06 (5.0 pN), 0.77 ± 0.07 (5.4 pN), 1.13 ± 0.10 (5.76 pN)** *(fig)*. Last point ~1.2 pN below stall. |
| Visscher 1999, Fig. 4a (p.187) | r vs [ATP] at 1.05, 3.59, 5.69 pN | 5.69 pN: r = 0.74, 0.79, 0.83 at 0.4, 2, 5 mM *(fig)*; "unable to measure r for 5.69 pN at more than a few ATP levels" (p.187) |
| Block 2003, Fig. 4C (p.2354) | 1.6 mM, 4.2 µM | 1.6 mM: r = 0.36–0.42 from −4.7 (assisting) to 2.7 pN hindering; **0.55 ± 0.04 at 3.7 pN, 0.86 ± 0.06 at 4.8 pN** *(fig)*; text: r = 0.38 ± 0.01 over "−3 to 5 pN" in Block's convention, i.e. from 3 pN hindering to 5 pN assisting. 4.2 µM: 0.8–1.2 |
| Yildiz et al. 2008, *Cell* 134:1030, p.1033 | one head labelled, low ATP | r = 0.57 (WT); 0.97–1.01 for neck-linker-extended mutants (futile cycles) |

**(c) Forward:back ratio against load**

| dataset | how close |
|---|---|
| Carter & Cross 2005, Fig. 2a inset (p.310) | ratio = 802 e^(−0.95F) through 1:1 at 7.04 pN (their fit); 1 mM points 2.2–10.2 pN, down to ~0.08 *(fig, vector)*; the 10 µM points reach 1:1 only near 8.7 pN *(fig)* |
| Nishiyama 2002, Fig. 4b and Table 1 | R0 = 221, d = 2.9 nm (d_f = 3.0, d_b = 0.1 nm); 1:1 at 7.6 pN |
| Taniguchi 2005, p.343 and Table 1 | P(forward) = 0.5 at ~8 pN at every temperature; ln(ratio) linear in load ("data not shown") |
| Kondo 2023 (session 1) | 44:1 unloaded; lever 2.22 nm above 3.19 pN; 1:1 at 9.3 pN |

**(d) Dwell distributions near stall**

| dataset | finding |
|---|---|
| Carter & Cross 2005, p.309 | "At any particular backward load, dwell times are exponentially distributed (Supplementary Fig. S2)" — supplement not in folder. Mean backward-step dwells (Fig. 2b, *fig*): at 1 mM, 59, 138, 113, 196 ms at 3.6–6.6 pN against forward 30, 41, 66, 152 ms; similar near stall (235 vs 246 ms at 7.5 pN); **load-independent superstall** (~250–460 ms, 7.5–14.5 pN); at 10 µM ~0.4–2.3 s. |
| Nishiyama 2002, p.793 | dwells before back steps and detachments "could be essentially described by the curve fit used for the 8-nm steps" |
| Taniguchi 2005, p.343 | forward dwell histograms need **two** rate-limiting transitions even at 35 °C and 5.5–7.5 pN (χ²(29) = 30 vs 110 for one); backward and detachment time constants equal to forward |
| Kondo 2023, pp.467–468 | slow-event dwells single-exponential; backward/forward dwell ratio 1.22 ± 0.23, detachment/forward 1.21 ± 0.32 (3–12 pN); fast backsteps with 0.4-ms dwells |

**(e) Viscosity**

Sozański et al. 2015, *PRL* 115:218102 (297 ± 2 K, 1 mM ATP, TIRF at
10 frames/s, truncated GFP-kinesin-1):
* Fig. 2b *(fig)*: speed falls from ~0.75–0.84 µm/s to ~0.4 µm/s by
  η_eff/η₀ ≈ 3 (small crowders). "No motion observed" at η_eff/η₀ = 3.5–6 for
  PEG 6k, PEG 18k and TetraEG, but sucrose still moves (0.15 µm/s) at 8.8.
* Text (p.218102-2): "stalling at about 5 mPa·s"; "movement events become rare
  … up to the point where no motion is observed" (p.218102-3). So the "stop" is
  partly a loss of runs, not necessarily a stalled stepping cycle.
* Michaelis–Menten (p.218102-3): V_max = 103.4 ± 0.2 s⁻¹, K_M = 311 ± 2 µM
  (control); 52.4 ± 1.4 s⁻¹, 474 ± 54 µM (30% sucrose).
* Their fit maps one unit of η_eff/η₀ to ~1.5 pN of external load (p.218102-4).
* Their model quotes Clancy 2011 for 600 s⁻¹ (tethered-head step, "2→3") and
  ~100 s⁻¹ (ATP dissociation); these are v1's search and clock rates.

The limit is approached closely in the sense that motion ceases, but the
velocities of the last moving points are noisy and the crowders disagree.

**(f) Temperature**

| dataset | range | finding |
|---|---|---|
| Taniguchi 2005, Tables 1–2 (p.345) | 280–308 K, bovine brain, 1 mM ATP | k_f0 = 100 → 1353 s⁻¹, k_b0 = 0.28 → 3.8 s⁻¹; ΔH‡ = 18.3 ± 1.1 (forward), 18.2 ± 1.4 k_BT₀ (back); ΔS‡ = 9.9 ± 1.1, 3.9 ± 1.4 k_B; slow dwell phase doubles per 10 °C drop, fast phase less (p.344). *Added in Part C:* these ΔH‡ exclude the diffusion-limited prefactor A_T = k_BT/(3πηr x₀²) (Methods, p.347). The raw Arrhenius enthalpies of Table 1's k_f0 and k_b0 are 26.7 ± 1.1 and 26.4 ± 1.4 k_BT₀ (my weighted fits); k_c (their load-independent step) has 13.8 ± 0.7 k_BT₀. |
| Hong et al. 2016, *Biophys. J.* 111:1287, p.1289 | ~5–27 °C (Arrhenius fit to 8 °C), KIF5A | velocity Ea ≈ 65 kJ/mol (≈ 26 k_BT at 298 K); force 5.3 ± 0.2 pN (295 K) vs 5.2 ± 0.2 pN (280.5 K) |
| Kawaguchi 2008 (review) | — | no new temperature data |

Neither end of 5–40 °C is reached with backstep counts: Taniguchi covers
7–35 °C, Hong's velocities go down to ~5–8 °C.

**(g) Docking and neck-linker mutants**

| dataset | finding |
|---|---|
| Budaitis et al. 2019, *eLife* 8:e44146, p.6 | KIF5C(1-560): cover-neck-bundle (CNB), Latch and CNB+Latch mutants are **faster** unloaded (0.771, 0.761, 0.788 vs 0.617 µm/s WT) and more processive (2.07, 4.27, 5.33 vs 0.99 µm), but detach at 0.91, 0.84, 0.81 pN without stalling (WT 4.6 ± 0.8 pN). No ATPase measured; MD (p.8) suggests enhanced catalytic-site closure. Room temperature. |
| Khalil et al. 2008, Table 1 (p.19250) | 2G mutant: unloaded v and run length "at least those of WT", "suggesting that its ATPase machinery was unaffected" (an inference, p.19249); stall 61% of WT |
| Yildiz 2008, p.1033 | neck-linker-extended mutants: k_cat similar to WT but 5–10× slower motion (futile hydrolysis) |
| Andreasson et al. 2015, *eLife* | neck-linker length series; cysteine-light caveat (session 1) |

**(h) Where the waiting head is (for p) and how floppy the undocked tether is (for p₀)**

| dataset | finding |
|---|---|
| Mickolajczyk et al. 2015, *PNAS* E7186, p.E7187 | labelled head's substep 8.4 ± 3.0 nm (mean ± SD, 78% component) from its rear site: the tethered intermediate sits over the bound head; 1HB state 8.0 ± 0.5 ms at 1 mM; the motor **waits for ATP with two heads bound** (≥ 10 µM), and the tethered state follows ATP binding and precedes hydrolysis; "ATP binding only partially docks the NL and hydrolysis completes docking" |
| Guydosh & Block 2009, *Nature* 461:125, pp.126–127 | with a bead on one head, ±1.7 pN alternating loads move the unbound head through ~23 nm; "similar transitions" at ±0.4 pN (p.126); "the unbound head is mobile, and can be readily pulled about its bound partner and even rotated" (p.127) |
| Kutys et al. 2010, *PLoS Comput. Biol.*, p.4 | WLC neck linker: L_p = 0.7 nm fits MD; 0.364 nm per residue (14 residues ≈ 5.1 nm) |
| Xu et al. 2021, p.2629 | WLC with l_c = 10 nm, l_p = 0.5 nm |

---

## Part A. The diffusive head race: v1 rebuilt, v3 built

### A.1 The first-passage tools (`competing_exits/diffusion.py`)

The overdamped search in a landscape U(x) on [−8, 8] nm is a
Scharfetter–Gummel birth–death chain with rates (D/h²)·B(ΔU), where
B(z) = z/(eᶻ − 1). It obeys detailed balance with respect to e^(−U) exactly,
and converges as h².

Two solvers use the same chain:
* **Capture.** Capture probabilities and mean capture times come from the
  exact birth–death nested sums, evaluated in log space. They are
  machine-precise even behind 40-kT walls.
* **Races.** Races against Poisson clocks (the gate, ATP release) are solved
  as absorbing chains. They give each exit's probability and the first two
  moments of its time.

Reflecting walls use half cells, which makes free diffusion with a wall exact.

Checks (`tests/test_diffusion.py`):
* the continuum nested-integral formulas: splitting to 10⁻⁷, mean time to
  2×10⁻⁵ at N = 6400;
* the free-diffusion closed forms with a clock and with a reflecting wall;
* h² convergence.

### A.2 v1 reproduced (`tests/test_headrace_v1.py`)

| target | brief | this code |
|---|---|---|
| unloaded search rate | 600 s⁻¹ (sets B = 9.144 kT) | 600.2 s⁻¹ |
| front barrier, 0 → 6 pN | 10.98 → 14.42 kT | 10.98 → 14.42 |
| rear barrier, 0 → 6 pN | 20.06 → 14.67 kT | **19.98** → 14.67 |
| speed at 0, 3, 6 pN | 780, 454, 104 nm/s | 780.5, 454.4, 103.6 |
| P_win at 0, 3, 6 pN | 0.857, 0.533, 0.121 | 0.857, 0.533, 0.121 |
| front = rear | 6.17 pN | 6.171 pN |

**The choices the brief left open:**
* **Barrier shape.** A Gaussian with **standard deviation** 0.7 nm. It is the
  only one of nine shapes tried that gives 600 s⁻¹ at B = 9.144 kT; the
  others give 470–4350 s⁻¹. The nine were: a Gaussian with 0.7 nm as the SD,
  the e-folding length or the FWHM; and cos², tent, parabola and box shapes
  with 0.7 nm as the half- or full width.
* **Sensitivity to that choice.** If B is re-tuned to 600 s⁻¹ for each shape,
  the other targets move to:
  - P_win(3 pN) 0.533–0.561 and P_win(6 pN) 0.121–0.139;
  - v(3 pN) 454–482 and v(6 pN) 104–121 nm/s.

  The 6.17-pN balance point does not move, because it is where the two
  barrier tops are level.
* **Barrier height.** "Barrier" means U(±6 nm) − min U.
  - This reproduces three of the four barrier targets to 0.01 kT.
  - The rear barrier at zero load comes out 19.98 kT, not 20.06. No shape or
    definition gives all four. The true maxima of the potential are
    11.01 → 14.51 and 20.16 → 14.76.
  - The 0.08-kT discrepancy has no effect on any dynamic target.
* **The rest of v1.** The speed targets fix it:
  - the step is d = 8.0 nm (8.2 gives 800 nm/s at zero load);
  - rear capture is a futile restart, with no displacement and no rest time,
    like the clock.

  The alternatives fail: rear capture as a −8-nm backstep gives 22 nm/s at
  6 pN, and as a futile event followed by the 8-ms rest, 96 nm/s.
* **The race approximation.** P_win is k_f/(k_f + k_r + k_c), with
  k = P(site)/E[search time]. The full diffusive race, with the clock as a
  uniform killing rate, agrees with it to four digits.

### A.3 v3: the docked landscape for the front search, a gate for the rear exit

**The model** (`competing_exits/headrace.py: head_race_v3`):
* **Front search.** After ATP binds, the head searches for the FRONT site
  (+8 nm) in a docked landscape: a harmonic tether (κ′, x′) plus a Gaussian
  binding barrier B at +6 nm. The geometry, D and the F/2 load rule are v1's.
* **Start.** The search starts from p = N(0.2 nm, 1/0.21 nm²), tilted by the
  load.
  - The mean is Mickolajczyk 2015's 8.4-nm substep minus 8.2 nm.
  - The width is the worm-like-chain neck linker of Kutys 2010:
    κ = 3/(2 L_p L) with L_p = 0.7 nm and L = 28 × 0.364 nm.
* **Rear.** The rear end reflects: there is no diffusive rear capture. A
  Poisson gate k_b(F) = k_b0·e^(−Fδ_b/kT) competes with the search instead, as
  the session-1 lesson requires.
* **Other exits.** The clock k_c (ATP leaves, restart). After either step
  comes the rest of the cycle, T (in Part B, two equal exponential stages).
  The ATP wait is 1/(k_on[ATP]).
* **Parameters.** Eight are fitted: κ′, x′, B ≥ 0, k_b0, δ_b, k_c, k_on, T.
  The others are fixed: kT = 4.087 pN·nm (23 °C), water 0.932 mPa·s,
  d = 8.2 nm, D = 93.1 nm²/µs (Stokes–Einstein, 2.5-nm head).
* **Race approximation.** The competing-exits treatment (capture as a Poisson
  exit) agrees with the full diffusive race with the gate and clock:
  - within 2×10⁻³ from 0 to 12 pN (`tests/test_bet.py`);
  - within 2×10⁻³ from −15 to +15 pN (`tests/test_bet.py`, added after the
    review; the largest deviation, 1.7×10⁻³, is for the floppy edge at
    −10 pN).

**The fit** is to C&C's three relations (ratio, dwell at 1 mM, dwell at
10 µM), each sampled at 3, 4, …, 9 pN as in session 1 (`analysis/part_a_v3.py`,
`results/part_a_v3.json`). There are two weightings:
* session 1's assumed per-bin scatter: σ(ln dwell) 0.10, σ(ln ratio) 0.20;
* the scatter of C&C's own plotted bins about their fits, which I digitised
  from the vector figure. The RMS is 0.31 for the 1-mM dwells (0.15 without
  one outlier) and 0.33 for the 1-mM ratio; rounded, I use 0.20 and 0.33.

| | session-1 weights | data weights | v2 refit (session 1) |
|---|---|---|---|
| χ² (21 points) | 16.54 | 4.33 | 16.14 |
| κ′ (kT/nm²), x′ (nm) | 0.410, 0.91 | 0.403, 0.87 | — |
| B (kT) | 0 (at its bound) | 0 (at its bound) | — |
| k_f(F = 0) (s⁻¹) | 2170 | 2260 | 5063 |
| k_b0 (s⁻¹), δ_b (nm) | 6.42, 0.58 | 6.30, 0.57 | 8.06, 0.69 |
| k_c (s⁻¹), k_on (µM⁻¹s⁻¹), T (ms) | 57.6, 1.07, 14.2 | 56.9, 1.06, 14.4 | 52.9, 0.98, 15.7 |
| 1:1 load (pN) | 7.15 | 7.13 | 7.06 |
| v at 0 pN, 1 mM (nm/s) | 522 | 516 | 454 (v2) |

**What the fit says.**
* **Fit quality.** The diffusive front search fits C&C as well as v2's
  exponential forward rate does.
* **No barrier needed.** The best docked tether needs no binding barrier
  (B = 0). The capture barrier is then the tether's own rise to the site,
  (κ′/2)(8 − x′)² ≈ 10.3 kT.
* **Effective load distance.** Under 3–9 pN the head's mean position in the
  tilted docked landscape is 0.0 to −1.8 nm, which is 8.0–9.8 nm below the
  capture point (+8 nm). Half of that (the head feels F/2) is 4.0–4.9 nm,
  close to v2's δ_f ≈ 4.3–4.4 nm.
* **Unloaded speed still missed.** Session 1's miss remains: 522 against
  824 nm/s (the trap-off bead velocity, vector-extracted from C&C Fig. 2c;
  session 1 read 829).
* **The gate slows under load.** δ_b = 0.58 nm > 0: the fitted gate becomes
  slower under hindering load. This is inherited from C&C's dwells at
  3–9 pN, and Part C shows it fails past stall.

**Check against Kondo et al. 2023** (KIF5A, 25 °C), with v3 expressed in
Kondo's definitions (steps of a direction per unit time: that direction's
fraction divided by the mean dwell):
* **Forward rate: the same shape.**
  - v3 is nearly flat below ~2 pN (effective distance 0.29 nm) and falls
    above (3.7 nm at 4–8 pN).
  - Kondo's forward rate has a knee at 3.19 pN, with 0.17 nm below and
    1.78 nm above. Drosophila KHC is steeper under load.
* **Backstep rate: disagrees in structure.**
  - v3's effective backstep rate rises ~12× from 0 to 6 pN (0.19 → 2.2 s⁻¹)
    and then falls.
  - Kondo's rises 1.9× (1.99 → 3.78 s⁻¹) and keeps rising.
  - In v3 the gate fires only during the ATP-bound race, which is a small part
    of the cycle at low load. Kondo's gate acts from a state that fills the
    cycle.
* **Backward/forward dwell ratio: consistent.** v3 gives 1 by construction;
  Kondo reports 1.22 ± 0.23.

---

## Part B. The valley, the MaxEnt member, and the predictions

Part B is written up in full in **[PREDICTIONS.md](PREDICTIONS.md)**, committed as
`84a3d19`. PREDICTIONS.md is left as committed. The errors found in it since
are listed at the end of this Part, and the corrections are applied here and
in Part C.

**What git does and does not show about the order.**
* **What it shows.** Commit `84a3d19` (02:16) precedes the first commit of any
  comparison script or output (`1b4549b`, 02:35). `84a3d19` contains no
  comparison script or result, its scripts read nothing from
  `data/digitized/`, and PREDICTIONS.md has not changed since.
* **What it cannot show.**
  - An unrun skeleton of `analysis/part_c_confront.py` (data loaders and a TUR
    function) existed as an untracked file before `84a3d19`. It was first run
    after that commit.
  - The README had listed the command `part_c_confront.py` since `c34b200`.

  Git records neither an untracked file nor whether a script was run, so both
  facts rest on my word.
* **"Before comparing" is limited to the limit datasets.** Part A had already
  compared v3 with C&C (the fit) and with Kondo, and PREDICTIONS §0 lists
  everything I had seen.

**B.0 Why the I-projection family, and when it is the MaxEnt family.**
The minimiser of D(q‖p₀) under linear constraints ⟨φ_k⟩_q = c_k is
q ∝ p₀·exp(−Σ_k λ_k φ_k). For any feasible q, the Pythagorean identity gives
D(q‖p₀) = D(q‖q\*) + D(q\*‖p₀).

The capture rates give such constraints under three conditions. For capture
over a high barrier from a start inside the well (the Kramers regime, where
the rate does not depend on the start):

  k_f(F) ≈ D·[∫_well e^(−U−fx) dx · ∫_barrier e^(U+fx) dx]⁻¹.

If the transition state is sharp, load-independent and at x_TS, and the well
and transition state feel the same lever F/2, this becomes

  ln k_f(F) = const − f·x_TS − ln⟨e^(−fx)⟩_q, with f = F/2kT.

With the constant absorbed by the fitted barrier B, the rates at loads F_k
fix ⟨e^(−f_k x)⟩_q, which are linear constraints, so φ_k = e^(−f_k x).

**How well the three conditions hold.** A post-hoc check across the seven
members over 3–9 pN reproduces ln k_f to within 0.003–0.04, but only with a
member-dependent x_TS: 7.70 nm for the barrier-free best fit, 6.10–6.60 nm
when B > 0. With x_TS fixed at 8 nm, ln k_f drifts by 0.22–1.39 between 3 and
9 pN. Since x_TS moves with the refitted barrier, the two-parameter
I-projection family (φ at 3 and 9 pN) is a well-motivated heuristic, not the
exact MaxEnt family. Its members are valid members of the valley, so their D
values are **upper bounds** on the least commitment.

**B.1 The valley.**
* Under session-1 weights, the docked tethers that fit C&C within
  Δχ² ≤ 5.99 form a curved ridge: κ′ = 0.10–0.46 kT/nm², x′ = 0.4–5.2 nm.
  - One end is stiff and centred near the bound head, with no barrier.
  - The other end is floppy and parked forward, with a barrier of ~11 kT.
* Under data weights the ridge is wider: κ′ 0.06–0.52, x′ 0–8 nm.
* Two landscape parameters describe the ridge at the data's precision, but
  the binding barrier is refitted along it.
* **The confidence levels are nominal.** The 21 "points" are samples of three
  fitted two-parameter curves, not independent bins, so Δχ² ≤ 5.99 (and
  2.30) does not have its nominal coverage.
* The MaxEnt members sit at the valley's edge (Δχ² 5.2–5.99), so they depend
  on the level chosen. At the 68% level the least Gaussian-family commitment
  is 0.04–0.47 kT (session-1 weights), against up to 0.38 kT at 95%.

**B.2 The MaxEnt member and the budget.**
* **No bet (q = p₀).**
  - Under session-1 weights it misses C&C at every κ₀ in the bracket
    0.01–0.30 kT/nm² (Δχ² 18–89).
  - Under C&C's own scatter it fits at κ₀ = 0.21, the worm-like-chain value.
    But Δχ² is 5.05 against a nominal 5.99, so this is marginal.
* **Least commitment D(q‖p₀).** At most 0.38 kT in the Gaussian family and at
  most 0.17 kT in the two-parameter I-projection family. These are upper
  bounds on the least commitment over all landscapes that fit, within v3.
* **The MaxEnt bet's shape.** In the I-projection family it is a soft rear
  wall, not a forward shift: it suppresses excursions behind x ≈ −4 to −5 nm.
* **The budget, corrected (post hoc).** PREDICTIONS.md compared D(q‖p₀) with
  the ~1.2 kT docking budget.
  - For docking that only restricts the head (V ≥ 0), D is a lower bound on
    the cost. The exact least cost of producing q is ln max(q/p₀), the
    max-divergence (`bet.restriction_cost`; results/tables_partC.md, C.8).
  - With the exact cost, **the MaxEnt members remain affordable at their own
    κ₀**: I-projection 0.19 and 0.01 kT, Gaussian 0.71 and 0.23 kT (κ₀ 0.03
    and 0.21).
  - The budget cuts much more of the valley than PREDICTIONS says.
    - At κ₀ = 0.21 only 31% of it is affordable (session-1 weights; data
      weights 17%), against 88% and 57% committed.
    - The excluded points span κ′ 0.10–0.23 and x′ 0.8–5.2 nm: every
      landscape broader than the undocked tether (κ′ < κ₀), plus the
      forward-shifted ones.
    - For a floppy p₀ (κ₀ = 0.01) the stiffest landscapes are excluded
      instead.
  - A restriction-only docking potential cannot make the head floppier than
    it is undocked without paying ln max(q/p₀).

**Choice of p.** p, the head's distribution when docking starts, is the same
density the search starts from, for every κ₀. The reading "p = the undocked
equilibrium" is carried as a variant. It raises the imprint to 0.58–0.98 kT
when κ₀ = 0.03. For κ₀ = 0.03 that variant is internally inconsistent,
because the fit's start distribution is the narrower worm-like-chain one.

**B.3–B.4 Predictions and discrimination.** Seven members were carried to the
limits: two I-projection and two Gaussian MaxEnt members (κ₀ 0.03 and 0.21),
the best fit, a mid-valley member and the floppy edge. The conditions were:
* load from −15 to +15 pN, at the ATP concentrations of C&C, Visscher and
  Block;
* viscosity from 1 to 10 η₀;
* temperature from 5 to 40 °C.

Their predictions for the near-stall kinetics (v, the 1:1 load, randomness,
odds, dwell CV) differ by less than about 1.4 times a realistic measurement
precision. The only strongly discriminating feature is the unloaded
viscosity response, v(5η₀)/v(η₀) = 0.63–0.91. It follows the unloaded capture
rate, which the 3–9 pN constraints do not fix.

The imprint per committed step:
* at F = 0 to 2 pN: 0.006–0.23 kT for the MaxEnt members, and up to 1.02 kT
  for the floppy edge;
* under assisting load it is larger for some members (correction 2 below).

Figures 5 and 6 are in PREDICTIONS.md.

**Errors in PREDICTIONS.md found after the commit** (not patched there):
1. **Budget (R3, P8, the "affordable fraction" column).** It used D(q‖p₀),
   a lower bound on the cost. The exact least cost is ln max(q/p₀); see B.2.
   The MaxEnt members stay affordable, but the rest of R3 is wrong.
2. **B.4, "everywhere at ≤ 2 pN it stays below 1 kT, and below ~0.2 kT for
   the MaxEnt members".** This is contradicted by the committed JSON:
   - the floppy edge has 1.02 kT at F = 0;
   - under assisting load the imprint rises for some members. The I-projection
     κ₀ 0.03 member has 0.96 kT per step at −5 pN (51 kT/s); the Gaussian
     MaxEnt κ₀ 0.03 member has 0.33; the floppy edge has 1.52 at −5 pN and
     1.34 at −2 pN.

   So for the I-projection κ₀ 0.03 member the imprint is largest at *both*
   ends of the load range (P10). Caveat: under assisting load the search
   takes microseconds, so capture can outrun relaxation, and the
   full-relaxation D(p‖q) may overstate what is dissipated.
3. **R4, "~9 kT at −8 nm".** This holds for κ₀ 0.21; for κ₀ 0.03 the wall
   reaches 34 kT at −8 nm. (It also dips to −0.014 kT at positive x, a
   negligible breach of V ≥ 0.)
4. **Temperature.** Taniguchi's ΔH‡ excludes a T/η prefactor, so the
   committed rate scalings are too weak (C.6).
5. **"Git history records the order."** True for the committed files only;
   see the disclosure above.
6. **Accounting.** "One quench per committed step" was justified by saying
   x-independent binding loops cannot dissipate. That holds only if unbinding
   is x-dependent, i.e. detailed balance with a docking potential. In v3 the
   clock is x-independent, so bind/unbind loops are a driven (flashing)
   cycle. The primary accounting is therefore an **assumption** about the
   real motor, and the "every ATP binding a quench" upper bound is what v3
   literally implies.
7. **The κ₀ values "from Guydosh".** 0.01 and 0.03 are guesses. Guydosh
   report ~23 nm only at ±1.7 pN; at ±0.4 pN they report "similar
   transitions" (supplement not in the folder). The same displacement at
   loads 4× apart is itself non-Hookean.
8. **"The MaxEnt family" should read "the two-parameter I-projection
   family"** (B.0).

---

## Part C. Confrontation with the limit datasets

Everything in this Part was computed after the predictions were committed.
The script is `analysis/part_c_confront.py`, and it writes
`results/part_c_confront.json`, `results/tables_partC.md` and
`figures/fig7_confront.png`. Unless marked "post hoc", nothing in PREDICTIONS.md
was changed.

![confront](figures/fig7_confront.png)

**Figure 7.** Predictions (bands: all seven members; lines: the mean of the
MaxEnt members) against the data:
* (a) C&C force–velocity. The grey band is the constraint window.
* (b) C&C mean dwells before forward (filled) and backward (open) steps.
* (c) C&C forward:back odds at 1 mM (fitted) and 10 µM (not fitted).
* (d) Randomness: Visscher at 2 mM, Block at 1.6 mM.
* (e) Sozański's velocity against the effective viscosity at the head.
* (f) The TUR bound 2/r from the measured randomness. Also shown, per net
  forward step: the total dissipation (Δμ = 20.5 kT with C&C's odds),
  Ariga's hidden dissipation, and the predicted imprint.

### C.1 Scorecard

[v3] marks a test of the kinetic scaffold shared by all members; [bet] marks a
test specific to the hypothesis. "z" is (data − model)/s.e.m.

| # | prediction | data | outcome |
|---|---|---|---|
| P1 [v3] | v through stall. Past stall: most negative −10 to −12 nm/s at 9–10 pN, then back toward 0. Saturation under assisting load. | C&C Fig. 2c past stall: **−16 to −31 nm/s at 8.5–14.5 pN, most negative at 12.6–13.5 pN**. Assisting, 1 mM: 381–739 nm/s. Assisting, 10 µM: 130–309 nm/s against a predicted 63–93. | **Fails** past stall: the backward speed is 1.4–5.7× too slow (model −4.5 to −12.1 nm/s), with the wrong trend. Fails at 10 µM under assisting load (2–4×). At 1 mM under assisting load, only the fast members reach the data. |
| P2 [v3] | Odds at 10 µM equal those at 1 mM | C&C Fig. 2a: over all bins, ln(odds 10 µM / odds 1 mM) = +0.25 ± 0.28 (from the bins' scatter; counting errors not extractable). Above 6 pN, a **post-hoc** selection, +0.95 ± 0.29; the 10 µM points reach 1:1 near 8.7 pN. But C&C conclude the opposite: "The stall force does not seem to depend on the ATP concentration" (p.309). Nishiyama fit 1 mM and 10 µM with one relation (Fig. 4b). Visscher find a **lower** stall force at 5 µM (abstract; ~5.5 vs ~7 pN). | **Not decided.** Suggestive near stall, contrary to the authors' reading; the literature disagrees on the sign. |
| P3 [v3] | r against load | Block, 1.6 mM: r = 0.55 ± 0.04 and 0.86 ± 0.06 at 3.7 and 4.8 pN (model 0.47–0.58 and 0.80–0.92). Block at ≤ 2.7 pN: 0.36–0.42 (model 0.37–0.46). Visscher, 2 mM: fits to 3.6 pN, but 0.44–1.13 at 4.2–5.8 pN (model 0.58–1.62). Block, 4.2 µM: 1.20 ± 0.05 at 2 pN (model 0.91–0.93). | **Mixed.** Block near stall passes (z ≤ 2.1). Block at low load is 1–4 s.e.m. high. Visscher at 4.2–5.8 pN is z −4 to −7 (squid kinesin, whose force–velocity curve differs from C&C's; Visscher's Figs. 4a and 4b differ by ~0.3 in r near 5.7 pN, so these z overstate). Block at 4.2 µM and 2 pN is z +5.5. |
| P4 [v3] | r against [ATP] | Visscher Fig. 4a at 5.69 pN: 0.74–0.83 (model 1.44–1.53). At 1.05 and 3.59 pN (markers merged): 1.09–1.35 at ≤ 10 µM (model 0.77–1.07); 0.78–1.23 at 40–70 µM (model 0.37–0.70); 0.37–0.70 at 100–400 µM (model 0.33–0.58). | **Fails** at 5.69 pN (z ≈ −6 nominally; a different construct, and this panel's 0.79 ± 0.12 at ~2 mM disagrees with Fig. 4b's 1.13 ± 0.10 at 5.76 pN). The data sit above the model at ≤ 100 µM. |
| P5 [v3] | Backward and forward dwells equal; mean dwell past stall grows with load | C&C Fig. 2b: backward/forward mean dwell 1.3–3.4 at 3.6–6.6 pN (1 mM) and 1.3–4.3 at 3.4–8.5 pN (10 µM; the 4.3 rests on an anomalously short forward bin, and C&C note that at 5–8 pN each direction's dwell distribution is distorted by the other, p.310, p.312); 0.7–1.5 near and past stall (1 mM). Kondo 1.22 ± 0.23. Backward dwell past stall **252–348 ms, flat, at 10.5–14.5 pN** (model 0.68–1.81 s). At 10 µM: 0.38–1.42 s at 9.7–14.5 pN (model 3.7–11.2 s). | **Fails** at intermediate loads and past stall; consistent near stall and with Kondo. |
| P6 [v3, bet] | v(5η₀)/v(η₀) = 0.63–0.91, no stop below 10η₀ | Sozański Fig. 2b, small crowders: v/v₀ = 0.46–0.61 at η_eff = 2.8–3.9 η₀ (BSA, Dextran 10k, PEG 6k, sucrose), 0.12 at 3.2 (TetraEG); **no motion** for PEG 6k from 3.5 and PEG 18k at 4.5. Large crowders ≈ 1.0 up to 3.45 η₀. | **Fails** for small crowders (model 0.70–0.96 at those η). Every member underpredicts the slowdown 1 − v/v₀: the MaxEnt members by 3–13×, the floppy edge by 1.3–2.6×. The floppy edge comes closest, as B.4 anticipated, but it too misses. |
| P7 [v3] | 1:1 load ∝ T (+0.34%/K) | Taniguchi Table 1 (1:1 load = kT·ln(k_f0/k_b0)/(d_f − d_b)): 9.47, 10.51, 9.22, 8.92 pN at 280–308 K, slope **−0.36 ± 0.28%/K**. Hong's force 5.2 ± 0.2 vs 5.3 ± 0.2 pN at 280.5 and 295 K (ratio 0.98 ± 0.05; model 0.951). | **Fails** at 2.5 s.e.m. (Taniguchi, bovine-brain kinesin), but fragile: without the 287-K point the slope is −0.21 ± 0.30%/K (1.8 s.e.m.), and the errors ignore the covariance of k_f0 with d_f. Hong is consistent but weak. The velocity part of P7 was mis-scaled (C.6). |
| P8 [bet] | Commitment ≤ 0.38 kT, within the budget | none direct | Untested. Post hoc, the budget statement used a lower bound (B.2): the MaxEnt members stay affordable with the exact cost, but most of the valley does not. Context: Taniguchi's 6 k_BT directional bias is 3–6× the docking energy (p.345). |
| P9 [bet] | Imprint 0.013–0.165 kT per step at 2 pN (≤ 1% of hidden) | Ariga's ~16 kT hidden (second-hand) | consistent, not a test |
| P10 [bet] | Imprint rises toward stall (I-projection members) | no per-step hidden dissipation against load | Untested. Post hoc, the committed JSON shows the imprint also peaks under assisting load for some members (Part B, error 2). |
| P11 [bet] | Differences in which chemistry cancels ≤ 0.5 kT | no data | untested |

**What the scorecard says.**
* **v3 fails near most of the limits the data reach:**
  - past stall;
  - in the backward/forward dwell ratio at intermediate loads;
  - in small crowders (truncated GFP-kinesin, TIRF);
  - against temperature (bovine-brain kinesin; 1.8–2.5 s.e.m.);
  - near stall for Visscher's squid kinesin (a different force–velocity
    curve);
  - under [ATP] changes. Only suggestive, against the authors' reading.
* **The failures do not discriminate between members.** At every failing
  feature the seven members agree with each other to within the data's
  resolution; in viscosity they differ, but none approaches the data.
* **What is falsified is v3's specific choices:**
  - a load-independent ATP-binding rate;
  - a backstep gate that fires only during the ATP-bound race and slows under
    load (one Bell gate, δ_b fitted at 3–9 pN);
  - the search as the only viscosity-dependent step;
  - an entropic temperature scaling.

  A diffusive search raced against a gate is not falsified in general. And
  the bet's size is itself estimated inside v3, so these results do not
  exonerate the bet either: a revised scaffold may need a different bet.
* **The bet-specific predictions (P8–P11) cannot be tested** with any dataset
  in the folder.
* **"The signal is strongest near the limits, in the character of the
  curves" cannot be supported here.** In the model, the near-limit character
  is fixed by the constraints and the scaffold, not by the bet. The measured
  character differs from the model's for reasons every member shares.

### C.2 Per dataset (every Part 0 dataset)

**Load.**
* **C&C Fig. 2c** (1 mM and 10 µM, −15 to +15 pN).
  - Inside the 3–9 pN window the velocity bins are not independent of the fit.
    The model is 7–49% faster than the bins at 1.6–4.6 pN; its 1:1 load is
    7.0–7.2 pN, against a velocity zero crossing at ≈ 6.5 pN.
  - **Past stall** (independent): data −16 to −31 nm/s at 1 mM and −4.5 to
    −22 nm/s at 10 µM (8.4–14.5 pN). Model −4.5 to −12.1 and −0.7 to −2.0. The model's
    backward speed past stall collapses because the fitted gate slows under
    load (δ_b = 0.4–0.8 nm) and, in v3, fires only after ATP binds.
  - **Assisting loads, 10 µM:** stepping is 2–4× faster than v3 allows. v3's
    ATP wait (≈ 94 ms at 10 µM) is load-independent; C&C's dwells under
    assisting load are 22–50 ms. This fits **load-dependent ATP binding** (or
    ATP-independent forward steps under assisting load). The folder already
    reported load-dependent binding:
    - Visscher's abstract: loads "raise the apparent Michaelis–Menten
      constant"; K_M = 88 ± 7, 140 ± 6 and 312 ± 49 µM at 1.05, 3.59 and
      5.63 pN (Fig. 2 inset, p.185).
    - Block 2003, Fig. 4B and p.2354: K_M rises under backward load.

    v3's load-independent k_on, inherited from v2, conflicted with both
    before any comparison.
* **Visscher 1999 Fig. 3a** (squid kinesin, 2 mM and 5 µM).
  - At 2 mM the model, fitted to C&C, is 120–420 nm/s too slow at every load
    (z 8–32 for the MaxEnt members). Visscher's motor holds 478 nm/s at 5 pN, where C&C's gives
    ≈ 87 (interpolated between its 4.6 and 6.5-pN bins). The two force–velocity curves have different shapes, so no single
    parameter set fits both.
  - At 5 µM the model matches from 3.6 to 5.6 pN (z −3.0 to +2.3) but is
    16–42% low at ≤ 1 pN.
* **Block 2003 Fig. 4A** (1.6 mM and 4.2 µM, −8 to +5 pN).
  - At 1.6 mM: 100–270 nm/s too slow at low and assisting loads for the MaxEnt
    members (the known unloaded miss; only the floppy edge reaches the
    assisting-load data). Close near stall (310 vs 267–273 nm/s at 3.7 pN; 128 vs 138–149 at
    4.75 pN).
  - At 4.2 µM: data fall 8× from −1 to +4.25 pN, the model only 2–2.7×.
    Load-dependent ATP binding again.
* **Nishiyama 2002** (Table 1).
  - v(3.8)/v(0) = 0.25, against the model's 0.43–0.61.
  - Zero-load odds R₀ = 221, inside the model's 135–338.
  - Total load distance 2.9 nm, against the model's 3.3–3.7 nm (the slope of
    ln odds over 2–8 pN; a secant over 0–7.6 pN gives 2.8–3.4 nm).
  - 1:1 at 7.6 pN, against 7.0–7.2.
* **Kondo 2023.** See Part A: the forward-rate shape is shared and the
  backstep-rate structure differs. The 1:1 load is 9.3 pN (KIF5A; a different
  motor).
* **Khalil 2008.** The K401 stall force is 4.96 pN, a different construct
  from C&C's full-length KHC. Not a test.

**Randomness.**
* **Visscher Fig. 4b.**
  - r is 0.41–0.48 to 4.2 pN, where the model agrees (z within ±2.5) up to
    3.6 pN.
  - At 4.2–5.8 pN the model is too high (z −4 to −7) on an absolute-load
    axis.
  - Compared at equal v/v(0), however, the model is too **low**. At
    v ≈ 0.45 v(0), Visscher has r = 1.13 while the model has 0.57–0.63.
    Visscher's randomness rises while the motor still runs at half speed.
    No member does that.
* **Visscher Fig. 4a.**
  - At 5.69 pN: 0.74–0.83 against 1.44–1.53 (same construct caveat).
  - At low [ATP], r = 1.09–1.35 (six points) against 0.77–1.07. The data sit
    above 1, which a single-state race without backsteps cannot produce.
* **Block Fig. 4C.**
  - At 1.6 mM the near-stall character is right (3.7 and 4.8 pN). The model is
    0.02–0.08 too high at low and assisting loads.
  - At 4.2 µM the data rise from 0.77 (assisting) to 1.20 (2 pN); the model
    goes from 0.88 to 0.93.
  - My CSV misses 2 of the 10 points at 4.2 µM. The reviewer's by-eye
    readings are r ≈ 1.07 at 1 pN hindering (model ≈ 0.90) and ≈ 0.89 at
    4.4 pN assisting (model 0.88). Neither changes the picture.
* **Yildiz 2008.** r = 0.57 for one labelled head (16-nm steps, low ATP): a
  different observable, not compared.

**Odds.**
* **C&C.** In the fit window the 1 mM residuals are within ±0.7 in ln.
  - At 10 µM the points above 6 pN lie above the 1 mM points, by e^(0.95 ± 0.29)
    in a post-hoc selection. Over all bins the difference is +0.25 ± 0.28.
  - v3 cannot produce any ATP dependence, because its race is decided after
    ATP binds.
  - But the authors read their data as ATP-independent (p.309), and
    Nishiyama fit both ATPs with one relation (Fig. 4b). Without counting
    errors, P2 is undecided.
* **Nishiyama.** See above.
* **Taniguchi.** P(forward) = 0.5 at "about 8 pN" at all temperatures (p.343),
  counting detachments. See C.6.
* **Kondo.** 44:1 unloaded and 1:1 at 9.3 pN. The model gives 135–338 at
  zero load: the zero-load odds are an extrapolation and are not
  constrained.

**Dwells.**
* **C&C.**
  - Backward dwells are 1.3–4.3× the forward ones at intermediate loads;
    about equal near stall.
  - Past stall they are flat at ≈ 300 ms (1 mM), against a model that rises to
    0.7–1.8 s.
  - C&C report exponential dwells at every load (p.309). The model's CV is
    0.92–0.98 near and past stall: consistent.
* **Nishiyama (p.793).** Dwells before back steps follow the forward-step
  curve: consistent with P5.
* **Taniguchi (p.343).** Two rate-limiting transitions even at 35 °C and
  5.5–7.5 pN. The model's dwell has the rest of the cycle (two stages) plus
  the race. Its CV at 35 °C and 5 pN is 0.66–0.78 (1/CV² = 1.6–2.3):
  consistent.
* **Kondo.** Slow-event dwells are single-exponential, with a
  backward/forward ratio of 1.22 ± 0.23: consistent.

**Viscosity (Sozański).**
* **The control.** v₀ = 0.80 µm/s, from "about 800 nm/s" (p.218102-2).
  - Sozański's own Michaelis–Menten fit (V_max 103.4 s⁻¹ with 8-nm steps,
    K_M 311 µM) gives 631 nm/s at 1 mM instead.
  - A cluster of points at η_eff/η₀ ≈ 1.0 (0.72–0.84 µm/s, read by the
    reviewer; not in my CSV) supports ≈ 0.8.
  - With v₀ = 0.72 the small-crowder ratios rise by 11% (0.51–0.67). They are
    still below every member except, marginally, the floppy edge.
* **Small crowders** slow the motor far more than the model at the same
  η_eff: v/v₀ 0.46–0.61 at 2.8–3.9 η₀, against 0.70–0.96. The MaxEnt members
  underpredict the slowdown 3–13×, the floppy edge 1.3–2.6×.
* **Stops.** PEG 6k stops at η_eff 3.5 and PEG 18k at 4.5. PEG 18k runs at
  full speed at 3.45 η₀ and then stops, so the stop is not a viscous
  slowdown; Sozański note loss of runs (p.218102-3).
* **Large crowders** (PEG 1000k, Dextran 500k, PEG 18k ≤ 3.45) show no
  slowdown, and neither does the model at their low η_eff. Consistent, but
  weak.
* **Summary.** In v3 the diffusive search is too fast (600–3200 s⁻¹) to limit
  the unloaded cycle. The small-crowder data need a diffusive or
  viscosity-coupled step that limits the cycle at zero load. That is
  incompatible with the C&C-constrained search rate unless another
  η-sensitive transition exists; the η-sensitive gate variant makes it worse.

**Temperature.** See C.6. Taniguchi's 1:1 load does not rise with T. Hong's
force data are too coarse to decide.

**Docking and neck-linker mutants.** No prediction of mutant kinetics was
committed, so these are comparisons of scale only:
* **Budaitis 2019.** CNB and latch mutants detach at 0.8–0.9 pN (WT
  4.6 ± 0.8), and are faster and more processive unloaded.
* **Khalil 2008.** The 2G cover-strand mutant stalls at 61% of WT.
* **The contrast with the bet.** Weakening docking costs much of the force.
  The bet the C&C arrival times need is ≤ 0.4 kT and can be zero (κ₀ 0.21,
  data weights). So in the state-bet reading, docking's measured mechanical
  importance is not carried by the bet.
  - Either the undocked tether is floppy, which makes the no-bet landscape
    fail C&C by Δχ² 21–84;
  - or docking does something the arrival-time constraint does not see,
    through the gate, detachment or force transmission.
* **Yildiz 2008.** Neck-linker-extended mutants hydrolyse futilely (k_cat
  similar, motion 5–10× slower). Their chemistry per step is not WT's (C.5).

**Inputs, not tests.** Mickolajczyk 2015 (p: mean and substep spread),
Guydosh & Block 2009 and Kutys 2010 (the κ₀ bracket), Rice via Block 2007 and
Xu 2021 (the budget), and Taniguchi's enthalpies.

### C.3 The TUR

**Form.** For a Markov process in a non-equilibrium steady state (jump or
overdamped diffusion, time-independent rates, all dissipative degrees of
freedom included), any time-integrated current J_t obeys
Var(J_t)/⟨J_t⟩² ≥ 2k_B/Σ_t, where Σ_t is the total entropy produced in time t
(Barato & Seifert 2015; proved by Gingrich et al. 2016, and for finite t by
Horowitz & Gingrich 2017).

Take J = the motor's displacement, with ⟨X⟩ = vt and Var X = 2D_eff t at long
times. Then σ̇ ≥ k_B v²/D_eff. Dividing by the net stepping rate v/d gives,
per net forward step,

  Σ_step ≥ k_B·v·d/D_eff = 2k_B/r, with r = 2D_eff/(v·d).

The brief's form is right, provided "per step" means per **net forward**
step. Per ATP hydrolysed the bound is (2/r)·(P_f − P_b)/(P_f + P_b), which is
smaller.

**Assumptions, and how the data meet them.**
1. **Steady state at constant force.** A force clamp holds the force by
   moving the trap with feedback on the bead position. A feedback-controlled
   system can in principle beat the TUR through information flow. With clamp
   bandwidths far below the bead's relaxation this should be negligible, but
   it is not proved.
2. **Markovian, overdamped dynamics.** Hidden states only add to Σ, so the
   bound holds for the total, hidden dissipation included. Inertia is
   negligible for the bead.
3. **The measured r is the long-time randomness of the same current.**
   Instrument noise and drift inflate r and so weaken the bound, which keeps
   it conservative.
4. **Σ is total entropy production:** Δμ per ATP minus the work Fd per step
   (plus the assisting work when F < 0), with tight coupling (one ATP per
   step in either direction, C&C and Coy 1999).

**Application** (`results/tables_partC.md`, C.7):

| load range | bound 2/r (k_BT per net step) | fraction of the total |
|---|---|---|
| Visscher, 2 mM, 0.5–4.2 pN | 4.1–4.9 | 0.24–0.30 |
| Visscher, 2 mM, 5.0 / 5.4 / 5.76 pN | 2.95 / 2.60 / 1.78 | 0.17 / 0.13 / 0.07 |
| Block, 1.6 mM, −4.7 to +2.7 pN | 4.75–5.55 | 0.16–0.35 |
| Block, 1.6 mM, 3.7 / 4.8 pN | 3.62 / 2.32 | 0.24 / 0.14 |

The totals use Δμ = 20.5 kT and **C&C's odds, applied to Visscher's and
Block's motors (an assumption)**: 15–20 kT per net step at low hindering
loads, rising to 26 kT at 5.8 pN. They are larger under assisting load, where
the assisting work is dissipated too. Δμ = 20.5 kT is probably low for
buffers with no added ADP or Pi; a larger Δμ lowers every fraction.

**Informative, or slack?**
* **At low load it is informative about the total.** Any motor as precise as
  kinesin (r ≈ 0.4) must dissipate at least ≈ 5 kT per net step, a quarter to
  a third of Δμ.
* **Near stall it goes slack.** v → 0 with D finite, so the bound falls
  (to 1.8 kT at 5.76 pN) while the true dissipation per net step diverges
  (backsteps waste ATP).
  - At Visscher's 5.76 pN it captures 7% with C&C's odds (R = 3.4); 13% if
    that motor's odds were ≈ 10; ≈ 5% with Δμ = 25 kT.
  - Visscher's motor still runs at 45% of its unloaded speed there, so this
    point is not yet near its stall.
* **At low ATP** (r ≈ 1–1.3) it gives ≈ 1.5–2 kT per net step, against a Δμ
  itself lowered by kT·ln([ATP] ratio). Slack again.

**Against the model, the imprint and Ariga, per net step at 2 pN:**

| quantity | k_BT per net step |
|---|---|
| predicted imprint, MaxEnt members | 0.014–0.17 |
| TUR bound, measured (Visscher; Block) | 4.8 (4.2–5.5); 5.55 (4.9–6.4) |
| TUR bound, model's own randomness | 4.8–5.4 |
| Ariga's hidden dissipation | ≈ 16 |
| total (model; data) | 17.0–17.3; 16.8 |

The TUR bounds the total dissipation, not the hidden part, so comparing it
with Ariga's value is not a test. The imprint sits more than an order of
magnitude below the bound. **The TUR cannot say anything about the imprint.**

### C.4 Differences in which chemistry cancels

**Predicted** (PREDICTIONS §3.7), per committed step:
* hidden(6 pN) − hidden(2 pN) from the imprint alone: −0.08 to +0.48 kT for
  the MaxEnt members (−0.51 to +0.76 for all members);
* WT − docking-abolished mutant: +0.009 to +0.161 kT at 2 pN, +0.002 to
  +0.23 at 0 pN (MaxEnt members, κ₀ 0.21).

**Assumption.** The chemical part of the hidden dissipation is the same in the
two conditions, so the difference isolates mechanical changes such as the
imprint.

*Evidence for:*
* Δμ is a property of the solution, independent of load and mutation.
* Coupling is tight at the loads in question: one 8-nm step per ATP (Coy
  1999), and C&C's backsteps are ATP-dependent. So ATP per step does not
  change with load.
* Taniguchi's fit has a load-independent biochemical transition (k_c,
  p.344).
* For the 2G mutant, Khalil infer an unaffected ATPase from unloaded
  velocity and run length (p.19249). That is an inference, not a
  measurement.

*Evidence against:*
* **Load-dependent ATP binding** is reported in the folder:
  - Visscher's abstract and Fig. 2 inset: K_M 88 → 312 µM from 1.05 to
    5.63 pN;
  - Block 2003, Fig. 4B.

  Two datasets here are consistent with it:
  - C&C at 10 µM under assisting load steps with 22–50-ms dwells, against a
    ≈ 94-ms ATP wait in any load-independent-binding fit;
  - Block's 4.2 µM velocity falls 8× from −1 to +4.25 pN, where
    load-independent binding gives 2–2.7×.
* **Load-dependent chemical dissipation at fixed chemical free-energy
  drops.** Takaki's fitted cycle is one example: my reconstruction from their
  deposited code, not a number printed in their paper.
  - Its chemical free-energy drop per cycle is load-independent (7.34 kT),
    but its rate constants are not: the ADP-release pair k₁₂, k₂₁ carries a
    common load factor 2/(1 + e^(0.4·Fd/kT)).
  - Its chemical entropy production per cycle falls from 11.30 to 6.23 kT
    between 2 and 6 pN (1 mM), because the load shifts the occupancies.
  - That change of **5.1 kT** is 11–320× the predicted imprint difference
    over the same span (|Δ| = 0.016–0.48 kT for the MaxEnt members).
* **Mutants.**
  - Budaitis's mutants have no ATPase measurement, and their MD suggests
    changed catalytic-site closure (p.8). They are also 25% faster unloaded,
    so something in the cycle is not WT's.
  - Yildiz's neck-linker-extended mutants hydrolyse 5–10 ATP per productive
    step, changing the per-step chemistry by tens of kT.

**Verdict.** Between loads, the chemistry does not cancel at the needed
precision. Between WT and a docking mutant there is no evidence either way,
and some against.

### C.5 Circularity: where the separation of constraint and imprint is imperfect

* **p is used twice.** It is the search's start in the fit and the reference
  in D(p‖q). The fitted q adapts to p: the κ₀ 0.21 MaxEnt members have
  D(p‖q) = 0.006–0.007 partly because q ≈ p. The ±2-nm label offset moves
  the imprint to 0.37–0.56 kT (Gaussian MaxEnt, κ₀ 0.21), q held fixed. A
  refit with the shifted start would move q too; that was not done.
* **The out-of-window bins are not independent experiments.** C&C's bins past
  stall, under assisting load and at 10 µM come from the same traces and
  construct as the fitted relations. They are separate bins, not separate
  experiments.
* **Taniguchi is used on both sides.** Their enthalpies are inputs, and their
  Table 1 then tests the 1:1-load scaling that those enthalpies imply. The
  failure (C.6) is an internal tension in Taniguchi's own parametrisation
  (d_f rises with T) as much as a model failure.
* **Ariga's value is independent** of the fit (different construct and
  measurement), but it is second-hand and a single condition.

### C.6 Temperature, and an error in the committed temperature model

* **The 1:1 load.** From Taniguchi's Table 1, the 1:1 load
  kT·ln(k_f0/k_b0)/(d_f − d_b) is 9.47 ± 0.63, 10.51 ± 0.71, 9.22 ± 0.56 and
  8.92 ± 0.47 pN at 280, 287, 298 and 308 K. The weighted relative slope is
  −0.36 ± 0.28%/K.
  - The committed default (entropic landscape, equal enthalpies) predicts
    +0.34%/K. The difference is 2.5 s.e.m.
  - The result is fragile. Without the 287-K point the slope is
    −0.21 ± 0.30%/K (1.8 s.e.m.), and the errors ignore the covariance
    between the fitted k_f0 and d_f.
  - ln(k_f0/k_b0) is constant (5.83–5.88), as the entropic picture requires.
    But Taniguchi's fitted d_f rises from 2.4 to 2.8 nm, which cancels the kT
    growth.
  - None of the committed variants produces a flat 1:1 load; the enthalpic
    tether makes it rise even faster.
* **Error found after the commit.**
  - Taniguchi's ΔH‡ = 18.3 and 18.2 k_BT₀ are fitted after removing the
    diffusion-limited prefactor A_T ∝ T/η(T) (Methods, p.347). The raw
    Arrhenius enthalpies of k_f0 and k_b0 in their Table 1 are 26.7 ± 1.1 and
    26.4 ± 1.4 k_BT₀.
  - The committed temperature model instead gave k_f(0) a **total** enthalpy
    of 18.3 k_BT₀: 7.9 from D ∝ T/η, plus a 10.4 prefactor.
  - What this does not affect: anything at 23 °C, the 1:1-load prediction
    (which depends only on the forward–gate difference, 0.1 k_BT₀ in both
    readings) and the imprint.
  - What it does affect: every rate's temperature dependence is too weak by
    ≈ 8 k_BT₀. The committed default's unloaded-velocity activation energy is
    18.2 k_BT; correctly read, the stepping rates would carry ≈ 26 k_BT.
  - That matches Hong's velocity Ea ≈ 26 k_BT (65 kJ/mol). The committed
    "rest_H26" variant (24.9–25.5 k_BT for the unloaded velocity) comes
    closest. The velocity part of P7 is void; PREDICTIONS.md is left as
    committed.
  - Taniguchi's biochemical k_c has 13.8 ± 0.7 k_BT₀, below the 18.3 assumed
    for the rest of the cycle.

---

## Part D. The decisive Harada–Sasa experiment

**What would decide it.** The hypothesis's only signature that the scaffold's
failures cannot hide is the imprint, a dissipation. So the test is a
difference of hidden dissipation between two conditions of the **same**
molecules, in which the imprint changes and nothing else should.

**Design.**
1. **Measurement.** A force clamp on full-length KHC at 1 mM ATP and 23 °C
   (C&C's construct and conditions, so that the fitted valley applies). Each
   molecule alternates between 2 pN and 6 pN in interleaved segments. The
   bead's velocity spectrum C̃(f) and response R̃′(f), the latter by a small
   sinusoidal force modulation, give the probe dissipation by Harada–Sasa.
   Hidden dissipation per committed step = Δμ − work − probe dissipation,
   divided by the counted steps.
2. **Controls.**
   - Shift Δμ by known amounts, through [ADP] and [Pi] at fixed load, to
     check that the hidden-dissipation difference tracks ΔΔμ (an accuracy
     calibration of the whole chain).
   - Repeat the protocol on a docking-impaired mutant at 0.5 pN; mutants
     detach near 0.8–0.9 pN (Budaitis).

**Predicted differences** (per committed step):

| contrast | MaxEnt members | all members |
|---|---|---|
| hidden(6 pN) − hidden(2 pN), imprint only | −0.08 to +0.48 kT | −0.51 to +0.76 kT |
| WT − docking-abolished mutant at 0–2 pN | +0.002 to +0.23 kT | up to +1.0 kT |

For the MaxEnt members with p as tracked (κ₀ 0.21), the expected contrasts
are ≈ 0.01–0.15 kT.

**Precision needed.** Resolving a 0.1 kT difference at 2σ needs σ ≈ 0.035 kT
per step in each condition. For the largest MaxEnt contrast (0.48 kT),
σ ≈ 0.17 kT per condition. Both are against a hidden dissipation of ≈ 16 kT
per step: a relative precision of 0.2–1% on each.

**Data volume** (statistical floor; my estimate, corrected after the review):
* **Noise model.** The Harada–Sasa integral's noise is dominated by the
  bead's thermal velocity noise. Its two-sided spectrum is C̃ ≈ 2kT/γ up to
  the band limit f_max.
  - A periodogram estimate of γ∫C̃ df over ±f_max, from a record of length
    T_rec, has a standard error of 4kT·√(f_max/T_rec) (the ± frequencies are
    not independent).
  - The independently estimated response adds a comparable term, so I take
    σ_J ≈ 4√2·kT·√(f_max/T_rec) for the heat rate, and σ_J/k_step per step.
* **Assumed values.** f_max = 5 kHz (a few trap corner frequencies; a guess),
  and the model's step rates (55 s⁻¹ at 2 pN, 9 s⁻¹ at 6 pN).
* **Resolving 0.1 kT at 2σ** needs about:
  - 42,000 s of record at 2 pN (≈ 2×10⁶ steps);
  - 1.6×10⁶ s at 6 pN (≈ 1.4×10⁷ steps; 18 days of continuous clamp).
* **Resolving the largest MaxEnt contrast (0.48 kT)** needs ≈ 1,800 s
  (≈ 10⁵ steps) and ≈ 68,000 s (≈ 6×10⁵ steps).
* **In runs.** With Andreasson 2015's WT unbinding rate (1.11 s⁻¹, 0.60 nm;
  session 1), a run lasts ≈ 0.7 s at 2 pN and ≈ 0.4 s at 6 pN. So these
  targets need about 10³–10⁵ runs at 2 pN and 10⁵–10⁶ at 6 pN, each bead
  with its own calibration.

**Biggest confounds.**
1. **Chemistry that does not cancel.**
   - In Takaki's fitted cycle (my reconstruction from their deposited code,
     not a printed number), the chemical entropy production changes by
     5.1 kT between 2 and 6 pN at fixed chemical free-energy drops.
   - Load-dependent ATP binding is reported by Visscher and by Block (C.4).
   - That is 11–320× the predicted imprint difference, and a Harada–Sasa
     measurement alone cannot separate it from the imprint.
2. **Calibration.**
   - The work difference between 2 and 6 pN is 8.0 kT per step. A 3–5%
     force-calibration error (a guess) adds 0.24–0.40 kT of systematic error,
     comparable to the largest predicted signal.
   - Harada–Sasa also relies on the FDT holding at high frequency to
     calibrate the temperature scale, which S.-W. Wang 2018 shows is biased
     when the perturbation is asymmetric.
3. **The mutant contrast.** The mutant's ATPase is unmeasured (Budaitis), so
   the chemistry confound is unbounded there.

**Conclusion.** The predicted signal is 0.01–0.5 kT per step.
* A Harada–Sasa measurement of the Ariga type, which infers the hidden part
  as Δμ − work − probe dissipation, **cannot isolate the imprint from the
  chemistry and calibration confounds** as they stand.
* Statistically, resolving 0.1-kT contrasts needs 10⁶–10⁷ steps, and the
  largest MaxEnt contrast needs 10⁵–10⁶, under the assumed 5-kHz band.
* These confounds are argued, not proved insurmountable: a chemistry model
  validated to ≈ 0.1 kT per cycle and a force calibration to ≈ 1% would
  change the verdict.
* With what the folder supports, the imprint is **below realistic
  precision**. That is the result of Part D.

**A better, non-thermodynamic test** (outside the brief's Harada–Sasa frame)
is direct high-speed tracking of the free head before and after ATP binding
(p against q) under a 3–6 pN load. The two families make different
predictions:
* **I-projection MaxEnt members:** docking removes the rear tail
  (x < −4 to −5 nm). The mean moves by 0.06 nm (κ₀ 0.21) or 1.2 nm
  (κ₀ 0.03).
* **Gaussian MaxEnt members:** docking shifts the whole distribution forward,
  by 0.4 nm (κ₀ 0.21) to 1.6 nm (κ₀ 0.03). The other members shift it
  0.8–3.9 nm.

A before/after comparison with the same label cancels the ±2-nm label offset.
The required precision is ≈ 0.2 nm on the mean and ≈ 10% on the tail
fraction beyond −4 nm.

---

## Open issues

1. **The scaffold's failures past stall and at low ATP** need, at least:
   - a backstep route whose rate does not fall with load;
   - a backstep route that is not confined to the ATP-bound race;
   - load-dependent ATP binding.

   None was tried: the brief asked for the predictions to be fixed first.
   Any such model must be refitted and the valley redone.
2. **Construct dependence.** Visscher's squid kinesin and C&C's Drosophila
   KHC have different force–velocity shapes, so cross-construct
   comparisons test universality as much as the model. Block's near-stall
   agreement may be partly fortuitous.
3. **The unloaded speed** is still too low. The model gives 428–612 nm/s
   across members and ATP levels, against the trap-off or unloaded speeds of
   C&C (824), Block (668 ± 8.5, p.2353) and Visscher (~820 nm/s).
4. **The MaxEnt search is incomplete.** The I-projection family is MaxEnt only
   if the data fix exactly ⟨e^(−fx)⟩_q at the basis loads. The fitted
   transition state moves with the barrier (x_TS 6.1–7.7 nm, B.0), so they do
   not. A general minimum-D search over landscapes (a non-parametric
   I-projection against the exact capture-rate constraints) was not done.
   The reported D values are upper bounds.
5. **Ariga 2018, Takaki 2022 and Rice 2003 are not in the folder.** Their key
   numbers are second-hand.
6. **The temperature model's enthalpy reading** needs correcting (C.6) before
   any further temperature prediction.
7. **Sozański's stops** (PEG 6k, PEG 18k) look like loss of runs, which v3
   does not model (no detachment exit).
8. **The error model of the fit is nominal.** The 21 constraint points are
   samples of three fitted curves. A bootstrap over C&C's traces (not
   available) would be needed to calibrate the valley's confidence levels.
9. **Visscher's two randomness panels disagree** near 5.7 pN: 0.79 ± 0.12
   (Fig. 4a) against 1.13 ± 0.10 (Fig. 4b).
10. **Digitisation gaps.** Block 2003 Fig. 4C at 4.2 µM misses 2 of 10
    points, and Sozański Fig. 2b misses the η ≈ 1 cluster. Both were found by
    the reviewer; neither changes a conclusion.
11. **The imprint under assisting load** assumes that the head fully relaxes
    in q before capture. At −5 pN the search takes microseconds, so a
    non-quasi-static treatment (the dissipation of a quench interrupted by
    capture) is needed before those values are used.

## Guesses and figure-read values

| item | value | type | where used |
|---|---|---|---|
| head radius | 2.5 nm | guess (v1) | D |
| rest-of-cycle shape | 2 equal stages (CV² 0.5) | choice after reading Mickolajczyk and Taniguchi | B, C |
| enthalpy of clock, k_on, rest | 18.3 k_BT₀ (variants 10, 26.2) | guess | B temperature |
| Taniguchi ΔH‡ read as total enthalpy | 18.3 | **misreading**, see C.6 | B temperature |
| gate viscosity | η-independent (variant ∝ 1/η) | guess | B viscosity |
| p | N(0.2 nm, 1/0.21 nm²); mean ±2 nm label offset | mean from Mickolajczyk (p.E7187); width a WLC guess | A, B |
| κ₀ bracket | 0.01–0.30 kT/nm² | 0.21 from Kutys (p.4, WLC); 0.01 and 0.03 are **guesses** loosely based on Guydosh (p.126: ~23 nm at ±1.7 pN; "similar transitions" at ±0.4 pN, a non-Hookean pair); 0.30 a guess | B |
| docking budget | 1.2 kT (1–2) | second-hand (Block 2007 p.2991; Xu 2021 p.2627) | B |
| Δμ | 20.5 kT | Takaki's code reference line | C, D |
| Ariga hidden dissipation | ≈ 16 kT per step at 2 pN | 0.8 × 20.5, second-hand | C, D |
| C&C σ(ln) per bin | 0.20 / 0.33 (RMS 0.31 / 0.33 measured) | read off the vector figure | A, B, C |
| Sozański v₀ | 0.80 µm/s (range ±5% mine) | text, p.218102-2 | C |
| measurement precisions in B.4 | v 15 nm/s, 1:1 load 0.3 pN, r 0.06–0.10, ln odds 0.3, CV 0.05, v ratio 0.05 | guesses from the folder's error bars | B.4 |
| Harada–Sasa band f_max | 5 kHz | guess | D |
| force-calibration error | 3–5% | guess | D |
| C&C Fig. 2a, 2b, 2c | ratio, dwells, velocity bins | figure-read (vector extraction; s.e.m. not extracted) | A, C |
| Visscher Fig. 3a, 4a, 4b | velocity, randomness (s.e.m. from bar extents; Fig. 4a markers for 1.05 and 3.59 pN merged) | figure-read (raster) | C |
| Block Fig. 4A, 4C | velocity, randomness | figure-read (raster) | C |
| Sozański Fig. 2b | velocity against η_eff (colour blobs; overlapping markers possible) | figure-read (raster) | C |
| Taniguchi Tables 1–2 | rates, distances, enthalpies | table text (exact) | B, C |
| Takaki's chemical entropy production (14.35 / 11.30 / 8.69 / 6.23 kT at 0, 2, 4, 6 pN) | my evaluation of their deposited code | reconstruction, not a printed value | 0, C, D |
| Harada–Sasa noise | σ_J ≈ 4√2·kT·√(f_max/T_rec) | estimate (periodogram noise plus a response term of the same size) | D |
| run durations | from Andreasson 2015's WT unbinding (1.11 s⁻¹, 0.60 nm) | session-1 value, applied to C&C's construct | D |
| Sozański control cluster at η ≈ 1 | 0.72–0.84 µm/s | the reviewer's reading of the rendered figure; not digitised | C |
| Block Fig. 4C, the two missing 4.2-µM points | r ≈ 1.07 (1 pN hindering), ≈ 0.89 (4.4 pN assisting) | the reviewer's by-eye reading | C |

## Review fixes

An independent, fresh-context review re-derived the two key results:
* the I-projection minimises D(q‖p₀) under linear constraints (the
  Pythagorean identity), and D is monotone along rays in λ;
* the TUR per net step is ≥ 2k_B/r.

It reran the tests (all passed), checked the commit order, and checked about
85 numbers against the PDFs, tables and JSON. Its findings, and what was done:

| # | finding (severity) | fix |
|---|---|---|
| 1 | The budget compared D(q‖p₀), a lower bound on the cost of a restriction-only docking potential, with the 1.2 kT budget; the exact least cost is ln max(q/p₀) (must-fix) | `bet.restriction_cost` added and tested; C.8 of `results/tables_partC.md` recomputes affordability; B.2 and the summary corrected; PREDICTIONS error 1. The MaxEnt members stay affordable; most of the valley does not for κ₀ ≥ 0.21. |
| 2 | PREDICTIONS §4 says the imprint is < 1 kT at ≤ 2 pN; the committed JSON has 1.02 kT (floppy edge, F = 0) and 0.96–1.52 kT under assisting load (must-fix) | PREDICTIONS error 2; B.3–B.4 and P10 corrected; caveat that capture may outrun relaxation (open issue 11) |
| 3 | P2 verdict ignored C&C's own reading (p.309) and Nishiyama's single fit; the > 6 pN selection was post hoc (must-fix) | P2 is now "not decided", with the quotations, Visscher's opposite finding and the post-hoc label; figure 7c retitled |
| 4 | Pre-registration disclosure incomplete: untracked skeleton, README line, "before comparing anything" (must-fix) | Disclosure paragraph in Part B; summary and deliverables reworded |
| 5 | "No Harada–Sasa experiment can decide" overclaimed (must-fix) | Part D and the summary scoped to Ariga-type measurements and to the confounds as they stand, with what would change the verdict |
| 6 | The small-bet headline and "refutes the scaffold" were not conditional on v3 (must-fix) | "Within v3" throughout; the conclusion now lists the specific v3 choices that fail and says the results neither convict nor clear the bet |
| 7 | The MaxEnt claim needed a derivation and conditions; x_TS is member-dependent (6.1–7.7 nm) (should-fix) | B.0 added, with the derivation, the conditions and the x_TS check reproduced; "two-parameter I-projection family"; `bet.py` docstring |
| 8 | TUR "7%" depends on C&C's odds and Δμ; "consistent with Ariga" is not a test (should-fix) | Range 5–13% with the assumptions stated; sentence replaced |
| 9 | Sozański's control: their MM fit gives 631 nm/s; an η ≈ 1 cluster supports 0.8 (should-fix) | Both noted; conclusions checked at v₀ = 0.72; "several times" reworded to the computed factors |
| 10 | Cross-construct failures in the summary without caveats; P7 fragile (should-fix) | Constructs named; P7 without the 287-K point is −0.21 ± 0.30%/K (1.8 s.e.m.) |
| 11 | Load-dependent ATP binding was already in the folder (Visscher abstract and Fig. 2 inset; Block Fig. 4B) (should-fix) | Cited; "direct evidence" removed |
| 12 | Visscher Figs. 4a and 4b disagree near 5.7 pN (should-fix) | Noted in P3 and P4 and in open issue 9 |
| 13 | Takaki: the rates are load-dependent (only the chemical drops are not); the load span did not match; not labelled as a reconstruction (should-fix) | Corrected to 5.1 kT over 2→6 pN (11–320× the contrast), labelled as my reconstruction |
| 14 | The Harada–Sasa noise factor was too small (should-fix) | 4√2·kT·√(f_max/T_rec); record lengths and run counts recomputed |
| 15 | Δχ² levels are nominal; MaxEnt members sit at the edge; the no-bet fit is marginal (should-fix) | B.1 and B.2 and open issue 8; 68% values given |
| 16 | κ₀ = 0.01 and 0.03 "from Guydosh" are guesses (should-fix) | Relabelled (PREDICTIONS error 7; guesses table) |
| 17 | The per-committed-step accounting is an assumption; v3's x-independent clock makes bind/unbind loops driven (should-fix) | PREDICTIONS error 6 |
| 18 | "≈16 kT that Ariga report" is my 0.8 × 20.5 (should-fix) | Summary reworded |
| 19–20 | Page citations (Taniguchi A_T p.347; Guydosh quote p.127); Block's sign convention in Part 0; Block's unloaded speed 668 ± 8.5 (nits) | Corrected |
| 21 | Misquoted ranges: 1.4–5.7× past stall; model −4.5 to −12.1 nm/s; 10 µM dwell ratio 1.3–4.3 with an anomalous bin and C&C's 5–8 pN distortion note; C&C v(5 pN) ≈ 87; model r 0.57–0.63 at 0.45 v₀; Nishiyama's distance method; R4's wall height (nits) | Corrected |
| 22 | `bet.py` documented reverse = D(q‖p); the wall potential dips to −0.014 kT; the ±15 pN race check had no committed test (nits) | Docstring fixed; dip noted; test added (`tests/test_bet.py`) |
| 23 | Digitisation omissions (Block 4C, 2 points; Sozański η ≈ 1 cluster) (nit) | Noted with the reviewer's readings (C.2, guesses table, open issue 10); CSVs not edited |

