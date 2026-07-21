# Manuscript Revision Report (Updated)

## Status Summary

Most recommendations from the original report have been **implemented** in the current `manuscript.tex` (1408 lines). This document tracks what was done, what changed in the process, and the **remaining items** that still need attention.

---

## Part A: Manuscript (.tex) — Implementation Status

| # | Recommendation | Status | Notes |
|---|---|---|---|
| A1 | Abstract rewrite with novelty framing | ✅ Done | Lines 47–83 |
| A2 | Motivation paragraph (NLSE, not sinusoid) | ✅ Done | Lines 108–123 |
| A3 | Commercial simulators critique fixed | ✅ Done | Lines 180–207 |
| A4 | Gap analysis item (g) added | ✅ Done | Lines 229–237 |
| A5 | Contributions restructured | ✅ Done | Lines 286–326 |
| A6 | Beat length note | ✅ Partially | `Δn₀` changed to `5.0×10⁻⁸` for multi-section model (`L_B ≈ 31 m`). But phenomenological model still uses old `0.87×10⁻⁵`. See A15 below. |
| A7 | σ convention in Panel B | ✅ Done | Line 1051 |
| A8 | Temperature panel discussion | ⚠️ Revised differently | See A8a below |
| A9 | Subtitle | ✅ Done | Line 37 |
| A10 | Agrawal citation consistency | ✅ Done | All now use `\cite{agrawal5}` |
| A11 | RIN implementation note | ✅ Done | Lines 396–400 |

---

### A8a. Temperature Panel — What Was Implemented vs. Suggested

**My original suggestion** included a strong caveat that the null is an artefact:
> "...would not occur in deployed fibre... should not be interpreted as a quantitative prediction at the null point."

**What was actually written** (lines 1058–1069):
```
The QBER at 75~km exhibits a broad minimum around 40°C, where the
temperature-induced birefringence crosses zero
(Δn(T ≈ 42°C) ≈ 0). At this point the phenomenological rotation
angle θ ≈ 0, the polarisation is minimally disturbed, and the QBER
drops to a few percent.
```

**Assessment**: Improved from the original (which said QBER → 0%), but still presents the cancellation as a physical feature without caveat. The temperature null is a property of the phenomenological model's deterministic `Δn(T)`, not a physical prediction. In real fibre, thermal and residual birefringence have different physical origins and do not cancel. **Recommend adding a caveat sentence.**

---

## Part A: Remaining Items

### A12. Conclusion Contradicts the Revised Framing

**Lines 1205–1213**:
```latex
The primary contribution is the validation methodology: every equation
cites a specific numbered reference, and every component is verified
independently, providing a level of traceability absent from existing
quantum network and optical system simulators.
```

**Problem**: This says "the primary contribution is the validation methodology." The revised abstract and contributions section correctly frame the contribution as the *combination* of four capabilities (field + validation + composition + reproducibility). The conclusion needs to match.

**Suggested rewrite**:
```latex
We have presented an open-source optical channel simulation framework for
QKD that simultaneously provides complex-envelope field propagation,
per-component literature validation, verified system-level composition,
and seeded reproducibility --- a combination absent from existing
protocol-level and commercial simulators.  Five components --- CW laser,
Mach--Zehnder modulator, fibre channel (birefringence, CD, PMD,
attenuation), and avalanche photodiode --- are independently validated
against published literature (34~sources).  The system-level BB84
demonstration confirms that these validated components compose correctly,
producing physically sensible QBER under combined impairments.
```

---

### A13. Table 1 Column Count vs. Criteria Count

**Lines 229–237** list 7 criteria (a–g), but **Table 1** (lines 249–264) has only 6 columns and does not include a "System composition" column. A reviewer may notice criteria (g) is listed but not tabulated.

**Options**:
1. Add a 7th column "Sys. comp." to Table 1 (would show `No` for all rows except `Yes` for "This work")
2. Note in the caption that (g) is demonstrated in Section 5 and not tabulated here

**Recommend option 1** for completeness, but it adds table width. If space is tight, add a sentence after Table 1:
```latex
Criterion~(g) (system-level composition) is demonstrated in Section~5
and is not included in the table because no existing simulator provides
it in conjunction with the other six criteria.
```

---

### A14. σ Convention Consistency

**Line 1051** (Panel B) now correctly states:
```
Narrow pulses ($\sigma < 10$~ps RMS, where $\sigma = T_0/\sqrt{2}$
for a Gaussian pulse)
```

However, **line 826** (CD validation) uses the same convention but does not define it:
```
\sigma(z) = \sigma_0 \sqrt{1 + (z/L_D)^2},
\qquad \sigma_0 = \frac{T_0}{\sqrt{2}}
```

**Fix**: Add a brief footnote at line 809 where `σ(z)` is introduced:
```latex
(Agrawal Eq.~3.2.6, where $\sigma$ is the RMS width
and $\sigma_0 = T_0/\sqrt{2}$).
```

---

### A15. Two `Δn₀` Values Must Be Explicitly Distinguished

The manuscript now uses **two different `Δn₀` values** for the two birefringence models:

| Model | `Δn₀` | `L_B` at 1550 nm | Purpose |
|---|---|---|---|
| Multi-section (L < 2 km) | `5.0×10⁻⁸` (line 520) | ~31 m | Physically consistent with SMF-28 (Agrawal §4.1) |
| Phenomenological (L ≥ 2 km) | `0.87×10⁻⁵` (line 566) | ~0.18 m | Fitted to produce `L_char = 75 km` for long-haul QBER |

These differ by a **factor of 174×**. The paper does not explicitly state that the phenomenological `Δn₀` is a *fit parameter*, not a material measurement. A reviewer familiar with fibre optics will flag this.

**Fix**: Insert after line 566:
```latex
Note that this phenomenological $\Delta n_0 = 0.87\times10^{-5}$ differs
from the material beat-length value used in the multi-section model
($5.0\times10^{-8}$).  The phenomenological value is a fit parameter
calibrated so that $L_0 = 75$~km produces diffusive polarisation
scrambling over the 10--200~km range relevant to long-haul BB84,
rather than a measurement of the fibre's intrinsic birefringence.
```

---

### A16. Boundary Continuity at 2 km

The dual-model birefringence dispatches automatically at `L = 2000 m` (line 574–579). The paper does not discuss whether:
- The QBER is continuous across this boundary
- A discontinuity would be observable in a sweep
- The user can override the dispatch

**Fix**: Add after line 578:
```latex
The two models produce consistent polarisation statistics in the
vicinity of the 2~km boundary: at $L = 2$~km, both the ordered product
of $\sim$20 multi-section matrices and the phenomenological single
rotation produce a net SO(3) rotation of approximately $0.1$--$0.3$~rad,
so the transition is smooth and no discontinuity appears in distance
sweeps.
```

(Note: verify this claim numerically before adding.)

---

## Part B: Birefringence Model — What Was Actually Implemented

The original report proposed replacing the entire `apply_birefringence()` with a pure multi-section model. The actual implementation chose a **dual-model dispatch** approach instead. This section describes what was actually built and its implications.

### B1. Architecture

The current `manuscript.tex` describes (lines 483–579):

1. **Multi-section model** (`L < 2000 m`):
   - `Δz = 100 m` sections with `Δn₀ = 5.0×10⁻⁸` → `L_B ≈ 31 m` → `Δz/L_B ≈ 3.2`
   - Each section: `Δβ·Δz = 2π·|Δn|·Δz/λ ≈ 2.8 rad` (SO(3) rotation)
   - Ordered product: `J_total = J_N · ... · J_1`
   - Converges to uniform SU(2) (Haar) distribution within ~1 km

2. **Phenomenological model** (`L ≥ 2000 m`):
   - Single SU(2) rotation with `θ = min(π, √(L/L_char)·π/2)`
   - `L_char = 75 km · (Δn₀/|Δn|)²` using `Δn₀ = 0.87×10⁻⁵`
   - Same temperature/bend sensitivity as before
   - Preserves distance-dependent structure for long-haul BB84

3. **Dispatch**: Automatic at 2000 m, overridable via `model` parameter

### B2. Why Not a Pure Multi-Section Model

The pure multi-section model converges to uniform SU(2) within ~1 km (confirmed by Fig. 4, the Poincaré sphere convergence panel). This means:
- For `L > 1 km`, the output polarisation is uniformly distributed on the Poincaré sphere
- QBER saturates at ~50% for all distances beyond ~1 km
- No distance-dependent structure is visible in QBER sweeps

This makes the pure multi-section model **unsuitable for long-haul BB84 simulation**. The dual-model approach is a pragmatic compromise: the multi-section model is correct for short fibres but the phenomenological model is needed to produce a meaningful QBER-vs-distance curve.

### B3. Key Weakness: Section-to-Beat-Length Ratio

```
Δz / L_B = 100 m / 31 m ≈ 3.2
```

Each section produces an SO(3) rotation of:
```
Δβ·Δz = 2π · Δn · Δz / λ = 2π · 5.0×10⁻⁸ · 100 / 1.55×10⁻⁶ ≈ 2.8 rad
```

This is a **large-angle rotation**, not a small-angle diffusive step. The "multi-section random walk" is coarse-grained: with only ~20 sections in 2 km, each producing ~2.8 rad, the convergence to uniform SU(2) happens after a few sections, not after many small-angle steps.

The manuscript acknowledges this in the limitations (lines 1168–1177) as a "practical compromise." This is honest but understates the limitation: with `Δz > L_B`, the model is not resolving the physical random walk.

**Potential improvement**: Reduce `Δz` to match `L_B` (e.g., `Δz = 10 m`, giving `Δz/L_B ≈ 0.32` and `Δβ·Δz ≈ 0.28 rad` per section), and increase the dispatch threshold from 2 km to allow more sections. This would require ~200 sections for 2 km, which is computationally feasible.

### B4. Impact on Validation Script

The validation script (`validate_birefringence.py`) has been rewritten to validate **both** models:

| Panel | Content | Validates |
|---|---|---|
| A | Mean \|Ex\|² vs distance (0–1.5 km) | Multi-section model (diffusive walk) |
| B | Mean \|Ex\|² vs distance (0–200 km) | Phenomenological model (√L law) |
| C | Mean \|Ex\|² vs temperature at 50 km | Phenomenological model |
| D | Mean \|Ex\|² vs bend radius at 50 km | Phenomenological model + Ulrich |
| E | Beat length vs wavelength | `L_B = λ/Δn₀` (material parameter check) |
| F | Δn components vs T for various R | Ulrich bend model |

Additionally, a **new Poincaré sphere figure** (`val_birefringence_poincare--seed42.png`) shows convergence to uniform SU(2) at 0.01, 0.1, and 1 km.

### B5. Impact on System Demo

The system demo (`val_system.py`) uses the **phenomenological model** because all sweeps are at distances ≥ 75 km. The multi-section model is only active for `L < 2 km`, which none of the demo panels probe.

**Implication for Panel C (QBER vs temperature)**: The phenomenological model does **not** include the stochastic `Δn_stoch` term or the `|Δn|` clamping — those were only added to the multi-section model description. The temperature null is therefore still deterministic in the system demo. If you want the softened minimum in Panel C, the stochastic term must be added to the phenomenological model's `Δn` calculation as well.

---

### B6. Impact on Test Suite

The `tests/test_fiber.py` suite was updated to test the `propagate()` entry point with the `model` parameter:

| Test | Status |
|---|---|
| `test_birefringence_power_conservation` | ✅ Passes (both models unitary) |
| `test_zero_length_edge` | ✅ Passes |
| `test_birefringence_dispatch_sectional` | ✅ Tests `model='sectional'` at 1 km |
| `test_birefringence_dispatch_phenomenological` | ✅ Tests `model='phenomenological'` at 100 km |
| `test_birefringence_dispatch_auto_short` | ✅ Auto selects sectional for L < 2 km |
| `test_birefringence_dispatch_auto_long` | ✅ Auto selects phenomenological for L ≥ 2 km |
| `test_temperature_sensitivity` | ✅ Both models |
| `test_bend_sensitivity` | ✅ Both models |
| `test_seeded_reproducibility` | ⚠️ May need separate seeds for each model |
| `test_impairment_toggles` | ✅ Tests `birefringence=False` bypass |

---

## Part C: Updated Implementation Order

Since many changes are already in the manuscript, the remaining items in priority order:

1. **`paperwork/manuscript.tex`** — Fix conclusion (A12)
2. **`paperwork/manuscript.tex`** — Add note explaining dual `Δn₀` values (A15)
3. **`paperwork/manuscript.tex`** — Add boundary continuity note (A16)
4. **`paperwork/manuscript.tex`** — Add Table 1 caveat for criterion (g) (A13)
5. **`paperwork/manuscript.tex`** — Add σ convention footnote in CD section (A14)
6. **`paperwork/manuscript.tex`** — Add caveat sentence to Panel C temperature discussion (A8a)
7. **`src/channel/fiber.py`** — Optionally add stochastic `Δn_stoch` to the phenomenological model so Panel C shows a softened minimum
8. **Regenerate figures** — If code changes are made
9. **Recompile**: `xelatex manuscript.tex`

---

## Appendix: Δn₀ Value Cross-Reference

| Location in manuscript | Value | Beat length `L_B` at 1550 nm | Role |
|---|---|---|---|
| Line 520, multi-section model | `5.0 × 10⁻⁸` | ~31 m | Physical — consistent with SMF-28 (Agrawal §4.1) |
| Line 566, phenomenological model | `0.87 × 10⁻⁵` | ~0.18 m | Fit parameter — calibrated to give `L_char = 75 km` |
| Old manuscript (v1) | `0.87 × 10⁻⁵` (only value) | ~0.18 m | Same as phenomenological, no distinction made |

The multi-section `Δn₀` is a **material property** (measured from the fibre's intrinsic birefringence). The phenomenological `Δn₀` is a **numerical fit parameter** (chosen so that `L_char = 75 km` produces scrambling over 10–200 km). These serve different purposes and should not be confused. The manuscript currently does not state this distinction explicitly — see A15.
