# Research Roadmap: opto-sim — Physical-Layer Fiber-Optic Simulator

---

## Core Principle

> This paper models the physical components of an optical system sufficiently,
> and it aligns with published literature or actual physical components.
> Additionally these components are modular and can be attached with one another
> to produce relatively proper real-like graphs and parameters.

Every component is validated against the published formula or datasheet it
implements — no more, no less. The paper demonstrates that the components
compose correctly in a modular pipeline.

---

## Publication Plan

| # | Paper | Target | When |
|---|-------|--------|------|
| A | **Software validation**: component-wise validation + modular composition | OE / JLT / SoftwareX | Now |
| B | QKD environmental sensitivity | OE / JLT / EPJ QT | After A |
| C | Physical-layer security analysis | IEEE T-IFS / NPJ QI / PRX Quantum | After B |
| D | Full QKD system | Nature Comms / PRX Quantum | Optional |

---

## Tier 0 — Infrastructure (DONE)

pytest suite, seeded reproducibility, 48 tests. No further work.

---

## Tier 1 — Software Validation Paper

**Paper claim, in one sentence:** Each component matches its source literature,
and they compose into a working end-to-end simulation.

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
| **1** | **~4 sessions** | **OE / JLT / SoftwareX** | **—** |
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
