# Research Roadmap: opto-sim — Physical-Layer Fiber-Optic Simulator

---

## Core Principle

Each component is validated against the analytic theory, published literature,
or datasheet from which it is derived. The paper demonstrates that these
independently validated models compose correctly into an end-to-end
optical-system simulation.

---

## Publication Plan

| # | Paper | Target | When |
|---|-------|--------|------|
| A | **Software validation**: component-wise validation + modular composition | **Journal** (see §Journal Recommendations below) | Now |
| B | QKD environmental sensitivity | OE / JLT / EPJ QT | After A |
| C | Physical-layer security analysis | IEEE T-IFS / NPJ QI / PRX Quantum | After B |
| D | Full QKD system | Nature Comms / PRX Quantum | Optional |

### Journal Recommendations for Paper A

Based on the paper's scope — open-source software, physically validated
optical channel model for QKD, component-wise literature verification —
the following journals are strong fits, ranked by suitability:

| Rank | Journal | IF | Why it fits | Similar published work |
|------|---------|-----|-------------|----------------------|
| **1** | **SoftwareX** (Elsevier) | ~3.8 | Purpose-built for open-source research software. Accepts papers that introduce software tools validated against benchmarks. The paper's emphasis on open-source code + reproducibility + validation aligns exactly with SoftwareX scope. | OptiEnvelope (2026), FBG_SiMul, PyLops, OptiGUI DataCollector, SPIROS (2026) |
| **2** | **Optics Express** (Optica) | ~3.8 | Broad optics journal; accepts simulation tools with experimental validation. QKD fiber-channel modeling is within scope. Open-access, high visibility. | QOSST (2024), CATNIP (2025) |
| **3** | **Journal of Lightwave Technology** (IEEE/Optica) | ~4.7 | Premier journal for fiber-optic systems and components. Strong fit for the fiber channel model (CD, PMD, birefringence, attenuation). Higher prestige, tighter scope. | Frequent QKD + fiber channel papers |
| **4** | **Optics Letters** (Optica) | ~3.6 | For concise, high-impact results. Viable if the paper is condensed to focus on the validation framework + key results. | CV-QKD fiber channel papers |
| **5** | **Quantum** (epjQT) | ~4.6 | Open-access quantum journal. Good fit for QKD simulation tools. More quantum-information focused than optics focused. | QOSST (2024) |

**Recommendation:** **SoftwareX** is the strongest fit. The paper is
primarily a software contribution (open-source framework + validation
pipeline), not a new physics result. SoftwareX reviewers will evaluate
the code quality, documentation, reproducibility, and validation — which
is exactly what we deliver. The 6-panel validation figures map cleanly to
the SoftwareX format.

**Alternative:** If the optics/physics contribution is emphasised over
the software contribution, **Optics Express** or **JLT** are viable.
JLT has the highest IF but expects deeper optical-systems insight;
Optics Express is broader and accepts simulation tools readily.

---

## Tier 0 — Infrastructure (DONE)

pytest suite, seeded reproducibility, 48 tests. No further work.

---

## Tier 1 — Software Validation Journal Paper

**Paper claim, in one sentence:** Each component matches its source literature,
and they compose into a working end-to-end simulation.

**Format:** Full journal paper (not conference). 6-panel validation
figures per impairment (CD, PMD, Attenuation, Birefringence, APD,
CWLaser, MZM, System). Literature cross-check table as supplementary.

### 1.1 Fiber Channel (DONE — validation scripts exist)

| Validation | Reference | Error | Script |
|------------|-----------|-------|--------|
| CD: Gaussian pulse broadening | Agrawal Fig 2.6 | 0.0000 % | `analysis/validation/validate_cd.py` |
| PMD: DGD histogram vs Maxwellian | Razavi Fig 2.11 | KS p = 0.31 | `analysis/validation/validate_pmd.py` |
| Attenuation: α×L at 1550 nm | SMF-28 spec (0.182 dB/km) | 0.0000 % | `analysis/validation/validate_attenuation.py` |
| Birefringence: Δn ∝ (r_f/R)² | Yuan Fig 1, Ulrich Eq 1 | 0.0000 % | `analysis/validation/validate_birefringence.py` |

### 1.2 Laser (NEW — validation scripts needed)

| Validation | Reference | What to check | Effort |
|------------|-----------|---------------|--------|
| RIN spectrum | Coldren Eq 5.3.38 | PSD peak frequency ±10 %, roll-off 1/f² | ~1 session |
| Phase noise → Lorentzian linewidth | Henry Eq 18 | Fitted FWHM = input Δν ±10 % | ~1 session |

No new physics. The RIN and phase noise models are already implemented
(`src/lasers/cwlaser.py`). These are scripts that confirm the output matches
the published formula.

### 1.3 MZM (NEW — validation script needed)

| Validation | Reference | What to check | Effort |
|------------|-----------|---------------|--------|
| Vπ from crystal geometry | Weis & Gaylord 1985 | Vπ within ±5 % of published value for given d, L, n_o, r | Minimal |
| Extinction ratio | MZM model | On/off ratio matches user-specified ER | Minimal |

Already tested in `tests/test_mzm.py`. One script to produce a publication figure.

### 1.4 APD (NEW — validation script needed)

| Validation | Reference | What to check | Effort |
|------------|-----------|---------------|--------|
| Responsivity R = η·e·λ/(h·c) | Kasap Eq 4.19 | R values match formula | Minimal |
| Excess noise factor F on shot noise only | Kasap Eq 4.45 | I_noise² = 2eIB·F, thermal not multiplied | Minimal |

Already tested in `tests/test_apd.py`. One script to produce a publication figure.

### 1.5 System Demonstration — Modular Composition (NEW)

| Demo | What it proves | Effort |
|------|----------------|--------|
| BB84 QBER vs distance with CD+PMD | All 4 components (laser, MZM, fiber, APD) working together | Already exists |
| Link budget: laser power → fiber → detected power | Power flows correctly through pipeline | ~1 session |

### Effort Summary

| What | Sessions | Status |
|------|----------|--------|
| Fiber validation scripts (4) | — | DONE |
| Laser validation scripts (2) | ~2 | TODO |
| MZM/APD validation scripts (2) | ~1 | TODO |
| System demo script | ~1 | TODO |
| Paper writing | — | TODO |
| **Total** | **~4 sessions** | |

No new physics. No spectral models. No Sellmeier. No D(λ) or α(λ) curves.
Just validation scripts that confirm each component matches its published source,
and a demonstration that they work together.

---

## Tier 2 — QKD Environmental Sensitivity (DEFERRED)

Requires Tier 1. Same scope as before.

---

## Tier 3 — Physical-Layer Security Analysis (DEFERRED)

Requires Tiers 1 and 2.

---

## Tier 4 — Full QKD System Optimisation (OPTIONAL)

---

## Summary

| Tier | Effort | Paper | Dependency |
|------|--------|-------|------------|
| 0 | DONE | — | — |
| **1** | **~4 sessions** | **Journal (SoftwareX / OE / JLT)** | **—** |
| 2 | 6 w | OE / JLT / EPJ QT | Tier 1 |
| 3 | 5 m | IEEE T-IFS / NPJ QI / PRX Q. | Tiers 1–2 |
| 4 | 10 m | Nature Comms / PRX Q. | Tiers 1–3 |

### Competitive Position

| Advantage | opto-sim | NetSquid | qkdX | Comsis (1998) |
|-----------|----------|----------|------|---------------|
| Physical-layer field propagation | ✅ | ❌ | ❌ | ✅ (Fortran) |
| Each component validated against its source literature | **in progress** | ❌ | ❌ | ❌ |
| Modular composition demonstrated end-to-end | **in progress** | ❌ | ❌ | ❌ |
| Open-source (Python, NumPy) | ✅ | ✅ | ❌ | ❌ |
| Laser RIN + phase noise | ✅ | ❌ | ❌ | ❌ |
| Temperature- and bend-dependent birefringence | ✅ | ❌ | ❌ | ❌ |
| FFT-based CD + PMD | ✅ | ❌ | ❌ | ❌ |
| APD with excess noise (Kasap) | ✅ | ❌ | ❌ | ❌ |
| Seeded reproducibility | ✅ | ❌ | ❌ | ❌ |
