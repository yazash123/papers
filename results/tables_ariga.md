### Ariga 2018 re-run (post hoc; conditions: 2 pN, 25 °C, Δμ = 84.5 pN·nm)

Ariga's hidden dissipation per step, from Table I: Δμ − F·d − J_x·d/⟨v⟩ (d = 8 nm, as they use)

| condition | ⟨v⟩ (nm/s) | work per step (pN·nm) | J_x per step (pN·nm) | hidden per step (pN·nm) | fraction of Δμ |
|---|---|---|---|---|---|
| 1 mM ATP | 575 | 16.0 | 0.89 | 67.6 ± 2.6 | 0.800 |
| 10 µM ATP | 205 | 16.0 | 0.16 | 68.3 ± 2.5 | 0.809 |

#### 1 mM ATP, native model, 23 °C

Ariga (experiment): power 1150 ± 120, Δμ/τ 6160 ± 560 pN·nm/s; authors' head-motion estimate ≈ 500 pN·nm/s

| member | v (nm/s) | power F·v (pN·nm/s) | Δμ × step rate (pN·nm/s) | backstep fraction | total dissipated per step (kT) | imprint per step (kT) | imprint rate (pN·nm/s) | imprint / total | TUR 2/r (kT per net step) |
|---|---|---|---|---|---|---|---|---|---|
| iproj_k0=0.03 | 392 | 784 | 4132 | 0.0113 | 16.75 | 0.165 | 32.9 | 0.0098 | 4.78 |
| iproj_k0=0.21 | 400 | 799 | 4219 | 0.0118 | 16.76 | 0.017 | 3.4 | 0.0010 | 4.86 |
| maxent_k0=0.03 | 418 | 836 | 4436 | 0.0144 | 16.78 | 0.101 | 21.7 | 0.0060 | 5.12 |
| maxent_k0=0.21 | 407 | 815 | 4308 | 0.0127 | 16.76 | 0.013 | 2.7 | 0.0008 | 4.95 |
| best | 449 | 899 | 4746 | 0.0122 | 16.76 | 0.472 | 108.4 | 0.0282 | 5.19 |
| middle | 442 | 884 | 4684 | 0.0139 | 16.77 | 0.060 | 13.6 | 0.0036 | 5.25 |
| floppy_edge | 463 | 927 | 4952 | 0.0178 | 16.81 | 0.649 | 155.4 | 0.0386 | 5.41 |

#### 10 µM ATP, native model, 23 °C

Ariga (experiment): power 410 ± 60, Δμ/τ 2190 ± 310 pN·nm/s; authors' head-motion estimate ≈ 5 pN·nm/s

| member | v (nm/s) | power F·v (pN·nm/s) | Δμ × step rate (pN·nm/s) | backstep fraction | total dissipated per step (kT) | imprint per step (kT) | imprint rate (pN·nm/s) | imprint / total | TUR 2/r (kT per net step) |
|---|---|---|---|---|---|---|---|---|---|
| iproj_k0=0.03 | 57 | 114 | 600 | 0.0113 | 16.75 | 0.165 | 4.8 | 0.0098 | 2.55 |
| iproj_k0=0.21 | 58 | 115 | 609 | 0.0118 | 16.76 | 0.017 | 0.5 | 0.0010 | 2.54 |
| maxent_k0=0.03 | 60 | 121 | 640 | 0.0144 | 16.78 | 0.101 | 3.1 | 0.0060 | 2.51 |
| maxent_k0=0.21 | 59 | 118 | 621 | 0.0127 | 16.76 | 0.013 | 0.4 | 0.0008 | 2.53 |
| best | 65 | 129 | 683 | 0.0122 | 16.76 | 0.472 | 15.6 | 0.0282 | 2.53 |
| middle | 64 | 127 | 673 | 0.0139 | 16.77 | 0.060 | 2.0 | 0.0036 | 2.51 |
| floppy_edge | 66 | 133 | 708 | 0.0178 | 16.81 | 0.649 | 22.2 | 0.0386 | 2.44 |

#### 1 mM ATP, 25 °C, temperature model as committed

Ariga (experiment): power 1150 ± 120, Δμ/τ 6160 ± 560 pN·nm/s; authors' head-motion estimate ≈ 500 pN·nm/s

| member | v (nm/s) | power F·v (pN·nm/s) | Δμ × step rate (pN·nm/s) | backstep fraction | total dissipated per step (kT) | imprint per step (kT) | imprint rate (pN·nm/s) | imprint / total | TUR 2/r (kT per net step) |
|---|---|---|---|---|---|---|---|---|---|
| iproj_k0=0.03 | 444 | 887 | 4677 | 0.0112 | 16.63 | 0.164 | 37.5 | 0.0099 | 4.78 |
| iproj_k0=0.21 | 453 | 905 | 4776 | 0.0117 | 16.64 | 0.017 | 3.9 | 0.0010 | 4.85 |
| maxent_k0=0.03 | 474 | 947 | 5023 | 0.0143 | 16.66 | 0.102 | 24.8 | 0.0061 | 5.11 |
| maxent_k0=0.21 | 461 | 923 | 4877 | 0.0126 | 16.64 | 0.013 | 3.1 | 0.0008 | 4.95 |
| best | 509 | 1018 | 5376 | 0.0121 | 16.64 | 0.470 | 123.1 | 0.0283 | 5.19 |
| middle | 501 | 1002 | 5306 | 0.0137 | 16.65 | 0.060 | 15.4 | 0.0036 | 5.24 |
| floppy_edge | 525 | 1051 | 5613 | 0.0177 | 16.68 | 0.651 | 178.1 | 0.0390 | 5.42 |

#### 10 µM ATP, 25 °C, temperature model as committed

Ariga (experiment): power 410 ± 60, Δμ/τ 2190 ± 310 pN·nm/s; authors' head-motion estimate ≈ 5 pN·nm/s

| member | v (nm/s) | power F·v (pN·nm/s) | Δμ × step rate (pN·nm/s) | backstep fraction | total dissipated per step (kT) | imprint per step (kT) | imprint rate (pN·nm/s) | imprint / total | TUR 2/r (kT per net step) |
|---|---|---|---|---|---|---|---|---|---|
| iproj_k0=0.03 | 64 | 129 | 680 | 0.0112 | 16.63 | 0.164 | 5.4 | 0.0099 | 2.55 |
| iproj_k0=0.21 | 65 | 131 | 689 | 0.0117 | 16.64 | 0.017 | 0.6 | 0.0010 | 2.54 |
| maxent_k0=0.03 | 68 | 137 | 725 | 0.0143 | 16.66 | 0.102 | 3.6 | 0.0061 | 2.52 |
| maxent_k0=0.21 | 67 | 133 | 703 | 0.0126 | 16.64 | 0.013 | 0.4 | 0.0008 | 2.53 |
| best | 73 | 146 | 773 | 0.0121 | 16.64 | 0.470 | 17.7 | 0.0283 | 2.53 |
| middle | 72 | 144 | 763 | 0.0137 | 16.65 | 0.060 | 2.2 | 0.0036 | 2.51 |
| floppy_edge | 75 | 150 | 802 | 0.0177 | 16.68 | 0.651 | 25.5 | 0.0390 | 2.44 |

#### 1 mM ATP, 25 °C, Taniguchi's enthalpies read correctly

Ariga (experiment): power 1150 ± 120, Δμ/τ 6160 ± 560 pN·nm/s; authors' head-motion estimate ≈ 500 pN·nm/s

| member | v (nm/s) | power F·v (pN·nm/s) | Δμ × step rate (pN·nm/s) | backstep fraction | total dissipated per step (kT) | imprint per step (kT) | imprint rate (pN·nm/s) | imprint / total | TUR 2/r (kT per net step) |
|---|---|---|---|---|---|---|---|---|---|
| iproj_k0=0.03 | 434 | 868 | 4576 | 0.0112 | 16.63 | 0.164 | 36.7 | 0.0099 | 4.73 |
| iproj_k0=0.21 | 443 | 887 | 4677 | 0.0117 | 16.64 | 0.017 | 3.8 | 0.0010 | 4.79 |
| maxent_k0=0.03 | 466 | 932 | 4942 | 0.0143 | 16.66 | 0.102 | 24.4 | 0.0061 | 5.04 |
| maxent_k0=0.21 | 453 | 905 | 4784 | 0.0126 | 16.64 | 0.013 | 3.1 | 0.0008 | 4.88 |
| best | 501 | 1001 | 5287 | 0.0121 | 16.64 | 0.470 | 121.1 | 0.0283 | 5.11 |
| middle | 493 | 987 | 5228 | 0.0137 | 16.65 | 0.060 | 15.2 | 0.0036 | 5.17 |
| floppy_edge | 523 | 1046 | 5589 | 0.0177 | 16.68 | 0.651 | 177.3 | 0.0390 | 5.40 |

#### 10 µM ATP, 25 °C, Taniguchi's enthalpies read correctly

Ariga (experiment): power 410 ± 60, Δμ/τ 2190 ± 310 pN·nm/s; authors' head-motion estimate ≈ 5 pN·nm/s

| member | v (nm/s) | power F·v (pN·nm/s) | Δμ × step rate (pN·nm/s) | backstep fraction | total dissipated per step (kT) | imprint per step (kT) | imprint rate (pN·nm/s) | imprint / total | TUR 2/r (kT per net step) |
|---|---|---|---|---|---|---|---|---|---|
| iproj_k0=0.03 | 64 | 129 | 680 | 0.0112 | 16.63 | 0.164 | 5.4 | 0.0099 | 2.57 |
| iproj_k0=0.21 | 65 | 131 | 690 | 0.0117 | 16.64 | 0.017 | 0.6 | 0.0010 | 2.56 |
| maxent_k0=0.03 | 69 | 137 | 727 | 0.0143 | 16.66 | 0.102 | 3.6 | 0.0061 | 2.53 |
| maxent_k0=0.21 | 67 | 133 | 704 | 0.0126 | 16.64 | 0.013 | 0.5 | 0.0008 | 2.55 |
| best | 73 | 147 | 776 | 0.0121 | 16.64 | 0.470 | 17.8 | 0.0283 | 2.55 |
| middle | 72 | 145 | 766 | 0.0137 | 16.65 | 0.060 | 2.2 | 0.0036 | 2.53 |
| floppy_edge | 76 | 152 | 811 | 0.0177 | 16.68 | 0.651 | 25.7 | 0.0390 | 2.46 |

#### Ariga's two-state fits (Fig. 3 legend), implied at 25 °C

| condition | odds f:b at 2 pN | 1:1 load (pN) |
|---|---|---|
| 1 mM | 6.9 | 4.11 |
| 10 µM | 312.4 | 9.46 |
