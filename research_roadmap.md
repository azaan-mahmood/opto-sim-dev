# Research Roadmap: opto-sim — Physical-Layer Fiber-Optic Simulator

A tiered publication plan with effort estimates, experimental-validation requirements, and scope boundaries.

---

## Perspective

opto-sim is a **general-purpose, open-source, validated physical-layer fiber-optic simulator**.
The complex-envelope electric field (`mean(|E|²) = Watts`) is the single source of truth.
BB84 QKD is one downstream application protocol, not the defining purpose.

This roadmap is organised around **four papers**. The first (Tier 1) validates the simulator
itself against experimental data. Tiers 2–4 are QKD-focused and depend on Tier 1 being published.

---

## Publication Plan

| # | Paper | Target | When | Depends on |
|---|-------|--------|------|------------|
| A | **Software validation**: validated physical-layer simulator | OE / JLT / SoftwareX | Now | — |
| B | **QKD environmental sensitivity**: temperature, bend, time-varying channel | OE / JLT / EPJ QT | After A | A |
| C | **Physical-layer security**: side-channel analysis | IEEE T-IFS / NPJ QI / PRX Quantum | After B | A+B |
| D | **Full QKD system**: decoy-state, SKR, open-source release | Nature Comms / PRX Quantum | Optional | A+B+C |

---

## Tier 0 — Infrastructure (COMPLETE)

| Effort | Deliverable | Paper |
|--------|-------------|-------|
| 1–2 weeks | pytest suite + seeded reproducibility + 48 tests | None (scaffolding) |

**Status: DONE.** No further work planned.

---

## Tier 1 — Software Validation Paper

**Goal:** Publish a validated physical-layer fiber-optic simulator that reproduces published
experimental data across all components.

**Novelty:** No open-source simulator validates CD, PMD, birefringence, attenuation, laser RIN,
phase noise, MZM transfer, and APD response against *both* analytic theory and
published experimental/datasheet curves.

### 1.1 Fiber Channel (DONE)

| Validation | Reference | Error | Script |
|------------|-----------|-------|--------|
| CD: Gaussian pulse broadening | Agrawal Fig 2.6 | 0.0000 % | `validate_cd.py` |
| PMD: DGD histogram vs Maxwellian | Razavi Fig 2.11 | KS p = 0.31 | `validate_pmd.py` |
| Attenuation: α = 0.182 dB/km at 1550 nm | SMF-28 spec | 0.0000 % | `validate_attenuation.py` |
| Birefringence: Δn ∝ (r_f/R)² | Yuan Fig 1 | 0.0000 % | `validate_birefringence.py` |

**What these prove:** The implementation matches the underlying analytic formula to
machine precision. These are *necessary* but not *sufficient* for a software paper —
they verify coding correctness, not physical fidelity against real-world data.

### 1.2 Fiber — Experimental Overlay (NEW)

**D(λ) dispersion curve against Corning SMF-28 datasheet**

| Aspect | Detail |
|--------|--------|
| What to add | Sellmeier coefficients for fused silica (Malitson 1965) → n(λ) → d²n/dλ² → D_material(λ) + D_waveguide |
| Effort | ~2 sessions |
| Good enough | |D_sim − D_datasheet| < 1 ps/(nm·km) across O–L bands |
| Too much | Full waveguide mode solver for arbitrary index profiles; temperature-dependent Sellmeier |
| Priority | **High** — D(λ) is the standard fiber datasheet curve; absence would be noticed |

**α(λ) attenuation spectrum against Corning SMF-28 datasheet**

| Aspect | Detail |
|--------|--------|
| What to add | Rayleigh (C_R/λ⁴) + IR tail (C_IR·exp(−α_IR/λ)) + UV tail + OH⁻ peak at 1383 nm |
| Effort | ~1 session |
| Good enough | |α_sim − α_datasheet| < 0.02 dB/km outside OH peak, < 0.1 dB/km within |
| Too much | Temperature-dependent α; radiation-induced loss; bend-loss spectrum |
| Priority | **High** — attenuation spectrum is the second standard fiber curve |

**Total fiber effort: ~3 sessions**

### 1.3 Laser Validation (NEW)

**RIN spectrum — S_RIN(f) against relaxation-oscillation model**

| Aspect | Detail |
|--------|--------|
| What to add | Validation script: generate RIN time series → Welch PSD → plot against Coldren Eq 5.3.38 |
| Effort | ~1 session |
| Good enough | Peak frequency within ±10 %; roll-off slope matches 1/f²; DC level matches RIN_0 |
| Too much | High-resolution RIN measurement with shot-noise calibration; relative intensity noise at 1/f |
| Priority | **High** — the relaxation-oscillation model is a unique feature of this simulator |

**Phase noise — Lorentzian linewidth verification**

| Aspect | Detail |
|--------|--------|
| What to add | Validation script: simulate field → complex-envelope PSD via Welch → fit Lorentzian → compare FWHM to input Δν |
| Effort | ~1 session |
| Good enough | Fitted FWHM within ±10 % of input linewidth |
| Too much | Gaussian line shape decomposition; 1/f noise divergence at low offset frequencies |
| Priority | **High** — linewidth validation is essential for any laser model |

**Total laser effort: ~2 sessions**

### 1.4 MZM Validation (NEW)

| Validation | What to add | Effort | Good enough |
|------------|-------------|--------|-------------|
| Vπ from crystal parameters vs literature | Validation script reporting Vπ for default geometry | None (test exists) | Within ±5 % of published values for given d, L, n_o, r |
| Extinction ratio | Plot E_out vs V, measure on/off ratio | None (test exists) | Matches user-specified ER parameter |
| Push-pull vs single-drive chirp | Compare phase of output for both modes | ~1 session | Single-drive chirp = 0.5× Vπ induced phase (Koyama & Iga) |

**Priority: Medium** — MZM model is well-tested. One validation script to publish alongside.

### 1.5 APD Validation (NEW)

| Validation | What to add | Effort | Good enough |
|------------|-------------|--------|-------------|
| Responsivity R(λ) | Spectral QE model (eta(λ) for InGaAs), plot R vs λ against datasheet | ~1 session | Peak responsivity ±10 %, roll-off shape matches |
| Shot-noise scaling | Verify I_noise² = 2eIB with F factor | None (test exists) | Within numerical tolerance |

**Priority: Medium** — simpler than fiber/laser, adds breadth to the paper.

### 1.6 System Demonstration

| Demonstration | Status | Effort |
|---------------|--------|--------|
| BB84 QBER under CD+PMD vs distance | Existing (`qber_vs_distance_dispersion.png`) | 0 |
| Full link budget: laser → fiber → detector | ~1 session | Low |

**Priority: Medium** — vendors like to see a complete system story.

### Tier 1 — Scope Boundaries

**IN SCOPE:**
- Fiber: D(λ) dispersion curve, α(λ) attenuation spectrum (Sellmeier + Rayleigh + OH⁻)
- Laser: RIN spectrum plot, Lorentzian linewidth verification
- MZM: Vπ, extinction, chirp-mode validation scripts
- APD: Spectral responsivity, shot-noise validation
- System: BB84 QBER demonstration, link budget

**OUT OF SCOPE (rejected for this paper — may belong in later tiers or separate work):**
- Nonlinear effects (Kerr, SPM, XPM, FWM) — 3+ months, separate paper
- EDFA / optical amplifier — not needed for passive-fiber validation
- PDL (polarization-dependent loss) — second-order effect, no experimental data to validate against
- Splice / connector loss — negligible publication value
- Temperature-dependent fiber (Sellmeier T-coefficients) — belongs in Tier 2
- Waveguide mode solvers — defeats the purpose of a macroscopic simulator
- Adding more protocols (E91, MDI-QKD, TF-QKD) — strong diminishing returns
- Quantum security analysis — belongs in Tier 3

**Effort summary for Tier 1:**

| Sub-tier | Sessions | New code | Status |
|----------|----------|----------|--------|
| 1.1 Fiber (analytic) | — | — | DONE |
| 1.2 Fiber (experimental) | ~3 | Sellmeier, α(λ) model | TODO |
| 1.3 Laser | ~2 | RIN/linewidth validation scripts | TODO |
| 1.4 MZM | ~1 | Validation script | TODO |
| 1.5 APD | ~1 | Spectral QE model | TODO |
| 1.6 System | ~1 | Link budget script | TODO |
| **Total** | **~8 sessions** | | |

**Estimated calendar time: 3–4 weeks combined coding + paper writing.**

---

## Tier 2 — QKD Environmental Sensitivity (DEFERRED)

**Requires:** Tier 1 published first.

**Goal:** Publish a journal paper on how environmental parameters (temperature, bend radius)
affect QKD system performance — a study only possible with a physical-layer simulator.

**Tasks** (unchanged from previous roadmap):
1. Temperature sweep (0–50 °C, 25 km spans) → Δn(T), beat length, QBER
2. Bend sensitivity (R = 10–100 mm) → Δn(R), coupled T×R regimes
3. Time-varying channel (1 °C/min ramp) → QBER time series, protocol recovery
4. Aerial vs buried comparison (ASHRAE weather data) → QBER variance

**Scope boundaries for Tier 2:**

| IN | OUT |
|----|-----|
| Temperature-dependent birefringence and phase drift | Temperature-dependent Sellmeier (negligible at 1 °C resolution) |
| Static and slow-ramp temperature profiles | Thermo-mechanical stress modeling |
| Single-span fiber (25 km) | Multi-span with amplifiers |
| BB84 phase-encoding | Other protocols |
| Bend-induced birefringence | Bend-induced attenuation (separate model) |

**Effort: ~6 weeks. Stop if QBER delta < 0.5 %.**

---

## Tier 3 — Physical-Layer Security Analysis (DEFERRED)

**Requires:** Tiers 1 and 2 published.

**Goal:** Publish a high-impact paper on a QKD side channel that *only a physical-layer
simulator can reveal*.

**Tasks** (unchanged from previous roadmap):
1. Laser-noise side-channel (RIN resonance → amplitude leakage, bits/symbol)
2. Finite-key effects with correlated noise (non-i.i.d. RIN/phase noise)
3. Detector side-channel (dead time, afterpulsing, time-shift attacks)
4. Phase modulator flaw analysis (drive ripple → deterministic phase errors)

**Effort: ~5 months. A single strong result beats a catalogue of re-runs.**

---

## Tier 4 — Full QKD System Optimisation (OPTIONAL)

**Requires:** Tiers 1–3.

**Goal:** Decoy-state protocol, SKR optimisation, open-source release.

**Effort: ~10 months. Strongly diminishing returns after Tier 3.**

---

## Summary

| Tier | What | Effort | Paper target | Dependency |
|------|------|--------|--------------|------------|
| 0 | Infrastructure | DONE | — | — |
| 1 | **Software validation** | **3–4 w** | **OE / JLT / SoftwareX** | **—** |
| 2 | QKD environmental | 6 w | OE / JLT / EPJ QT | Tier 1 |
| 3 | QKD security | 5 m | IEEE T-IFS / NPJ QI / PRX Q. | Tiers 1–2 |
| 4 | Full QKD system | 10 m | Nature Comms / PRX Q. | Tiers 1–3 |

### Recommended route

1. **Tier 1 now** (3–4 weeks). Publish a validated physical-layer simulator paper.
2. **If Tier 2's QBER signal is strong** (≥0.5 % delta), publish an environmental journal paper (6 weeks).
3. **Tier 3** (5 months). The laser-noise side-channel is the highest-impact result this platform enables.
4. **Stop** after Tier 3 unless you want a flagship capstone. Diminishing returns past this point are steep.

| Advantage | opto-sim | NetSquid | qkdX | Comsis (1998) |
|-----------|----------|----------|------|---------------|
| Physical-layer field propagation | ✅ | ❌ | ❌ | ✅ (Fortran) |
| Validated against multiple experimental datasets | *in progress* | ❌ | ❌ | ❌ |
| Open-source | ✅ | ✅ | ❌ | ❌ |
| Modern language (Python, NumPy) | ✅ | ✅ | ❌ | ❌ |
| Laser RIN + phase noise | ✅ | ❌ | ❌ | ❌ |
| Temperature- and bend-dependent birefringence | ✅ | ❌ | ❌ | ❌ |
| FFT-based CD + PMD | ✅ | ❌ | ❌ | ❌ |
| APD with excess noise factor (Kasap) | ✅ | ❌ | ❌ | ❌ |
| Seeded reproducibility | ✅ | ❌ | ❌ | ❌ |
| Modular component architecture | ✅ | ✅ | ✅ | ❌ |

This roadmap makes opto-sim the first **open-source, multi-experiment-validated,
physical-layer fiber-optic simulator** — a unique position in the landscape.
