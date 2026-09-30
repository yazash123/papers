# Competing exits, session 2: is neck-linker docking a state bet?

*(Plain summary: to be written when Parts A–D are complete.)*

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
| Block et al. 2003, *PNAS* 100:2351, Fig. 4A (p.2354) | 1.6 mM, 4.2 µM ATP; 2-D force clamp | 1.6 mM: ~650–700 → 128 nm/s at 4.7 pN hindering *(fig)*; "near zero between −4 and −6 pN" (p.2354) |
| Nishiyama et al. 2002, *NCB* 4:790, Table 1 (p.790) | bovine brain kinesin, 298 K | v = 930, 230, 0 nm/s at 0, 3.8, 7.6 pN (model rates from dwells) |
| Kondo et al. 2023, *Traffic* 24:463 | mouse KIF5A, 25 °C | zero velocity ~9.1 pN; loads to >40 pN (session 1) |
| Khalil et al. 2008, *PNAS* 105:19247, Table 1 (p.19250) | *Drosophila* K401 WT and cover-strand mutants | stall forces 4.96 ± 0.05 (WT), 3.02 ± 0.03 (2G), 1.37 ± 0.04 pN (DEL, biased high) |

**(b) Randomness against load and ATP**

| dataset | conditions | how close to the limit |
|---|---|---|
| Visscher 1999, Fig. 4b (p.187) | 2 mM ATP | r = 0.41–0.48 from 0.5 to 4.2 pN; **0.68 ± 0.06 (5.0 pN), 0.77 ± 0.07 (5.4 pN), 1.13 ± 0.10 (5.76 pN)** *(fig)*. Last point ~1.2 pN below stall. |
| Visscher 1999, Fig. 4a (p.187) | r vs [ATP] at 1.05, 3.59, 5.69 pN | 5.69 pN: r = 0.74, 0.79, 0.83 at 0.4, 2, 5 mM *(fig)*; "unable to measure r for 5.69 pN at more than a few ATP levels" (p.187) |
| Block 2003, Fig. 4C (p.2354) | 1.6 mM, 4.2 µM | 1.6 mM: r = 0.36–0.42 from −4.7 (assisting) to 2.7 pN hindering; **0.55 ± 0.04 at 3.7 pN, 0.86 ± 0.06 at 4.8 pN** *(fig)*; text: r = 0.38 ± 0.01 over "−3 to 5 pN". 4.2 µM: 0.8–1.2 |
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
| Taniguchi 2005, Tables 1–2 (p.345) | 280–308 K, bovine brain, 1 mM ATP | k_f0 = 100 → 1353 s⁻¹, k_b0 = 0.28 → 3.8 s⁻¹; ΔH‡ = 18.3 ± 1.1 (forward), 18.2 ± 1.4 k_BT₀ (back); ΔS‡ = 9.9 ± 1.1, 3.9 ± 1.4 k_B; slow dwell phase doubles per 10 °C drop, fast phase less (p.344) |
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
| Guydosh & Block 2009, *Nature* 461:125, p.126 | with a bead on one head, ±1.7 pN alternating loads move the unbound head through ~23 nm; "similar transitions" at ±0.4 pN; "the unbound head is mobile, and can be readily pulled about its bound partner and even rotated" |
| Kutys et al. 2010, *PLoS Comput. Biol.*, p.4 | WLC neck linker: L_p = 0.7 nm fits MD; 0.364 nm per residue (14 residues ≈ 5.1 nm) |
| Xu et al. 2021, p.2629 | WLC with l_c = 10 nm, l_p = 0.5 nm |

---

*(Parts A–D follow in later commits.)*
