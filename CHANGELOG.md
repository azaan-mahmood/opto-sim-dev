# Changelog

All timestamps are local time (UTC+5).

---

## 2026-07-20 — Birefringence recalibration (SMF-28), fiber.py → fiber_sectional.py

### Session: Recalibrate Δn₀ to SMF-28 literature, rename file

| Change | Files | Rationale |
|---|---|---|
| Recalibrated birefringence to SMF-28 | `src/channel/fiber.py` → `src/channel/fiber_sectional.py` | Δn₀ = 8.7e-6 → 5.0e-8 (L_B = 31 m, Agrawal §4.1 SMF-28 range 10–100 m). T_coeff = -5e-7 → -3e-9/°C (same ~6%/°C ratio). Clamping floor 1e-8 → 5e-10. Stochastic residual 10% of Δn₀. |
| Renamed fiber.py → fiber_sectional.py | `src/channel/fiber_sectional.py` | Multi-section model now has explicit filename. All imports updated. |
| Updated imports | `src/channel/__init__.py`, 5 validation scripts, `test_fiber.py` | s/src.channel.fiber/src.channel.fiber_sectional/g |
| Updated manuscript parameter values | `paperwork/manuscript.tex` | Δn₀ 8.7e-6 → 5e-8, L_B 0.18 m → 31 m. Removed "compressed correlation length" language — now SMF-28 realistic. |
| Updated validation Panel D units | `analysis/validation/validate_birefringence.py` | Beat length in metres (not mm) for new L_B ≈ 31 m. |
| Updated documentation | `AGENTS.md`, `journal_paper_outline.tex`, `CHANGELOG.md` | Fiber references point to fiber_sectional.py. |

**Tests:** 48/48 pass.

## 2026-07-20 — Rename fiber_sectional → fiber, cable → propagate, individual impairment flags, rework birefringence validation

### Session: renames, independent impairment toggles, split-model validation

**File rename:** `src/channel/fiber_sectional.py` → `src/channel/fiber.py`
- Updated all imports: `__init__.py`, `test_fiber.py`, 4 validation scripts
- Updated references in `journal_paper_outline.tex` and `AGENTS.md`

**Function rename:** `cable()` → `propagate()` in `src/channel/fiber.py`
- Updated all call sites: `analysis/val_system.py`, `tests/test_fiber.py`, all 3 BB84 protocol scripts, `main.py`
- Updated `__init__.py` export
- Updated AGENTS.md and CHANGELOG.md

**Independent impairment flags in `propagate()`:**
- `birefringence` (default True), `cd` (default None → uses `dispersion`), `pmd` (default None → uses `dispersion`), `attenuation` (default True)
- Backward compatible: `dispersion=True` still enables cd+pmd; `dispersion=False` (default) disables them
- New patterns: `propagate(..., birefringence=False, attenuation=False)` → no impairments;
  `propagate(..., birefringence=False)` → attenuation only;
  `propagate(..., attenuation=False)` → birefringence only
- `apply_birefringence()` also gets `enabled=True` parameter
- `dt` required only when `cd=True` or `pmd=True` (not just `dispersion=True`)

**Birefringence validation reworked (`validate_birefringence.py`):**
- Now validates both models explicitly: sectional (L < 2 km) and phenomenological (L ≥ 2 km) in separate test functions
- Added auto-dispatch test: verifies model selection at boundary
- Self-consistency checks now labelled per-model

**Tests:** 48/48 pass, validation clean.

## 2026-07-20 — Hybrid birefringence dispatch (sectional + phenomenological)

### Session: dual-model dispatch, bug fixes, performance tuning

**Hybrid birefringence model:**
- `apply_birefringence()` dispatches automatically via `model='auto'`:
  - **Short fibres** (L < `SECTIONAL_LIMIT` = 2 km, `model='sectional'`): multi-section ordered product of random-axis SU(2) matrices, L_B ≈ 31 m (Agrawal §4.1). For DV-QKD and DPS QKD.
  - **Long fibres** (L >= 2 km, `model='phenomenological'`): single SU(2) rotation with `θ = min(π, √(L/L_char)·π/2)` (Menyuk & Wai 1994), L₀ = 75 km. For long-haul BB84 with distance-dependent QBER.
  - Auto-dispatch required because multi-section model converges to uniform SU(2) within ~1 km regardless of parameters, producing flat ~50% QBER.
- Fixed: duplicate return statement, `section_length` default to 1.0, `propagate()` indentation, unused `L0` in validation

**System demo results (hybrid mode):**
- 0–70 km: 0% QBER → 80 km: 23% → 200 km: 68% (peak) → 500+ km: ~45–53% (dark count floor)
- Temperature: V-shaped null at ~35–50°C; pulse width: 21.7% (5 ps) → 1.7% (30 ps); bend radius: 69.2% (2 mm) → 19.6% (5 cm)

**Performance:** `section_length=1.0` → 67 ms/call (2 hrs sweep), `section_length=100.0` → 1 ms/call (2 min sweep). Sectional model defaults to 100 m sections.

## 2026-07-19 — Random-axis birefringence model, system demo 0–1000 km, manuscript 21 pp

### Session: random birefringence + extended demo + manuscript

**Random-axis birefringence model:**
- `src/channel/fiber.py`: new `_random_su2_rotation()` and updated `apply_birefringence()`:
  - SU(2) rotation around uniformly random Poincaré-sphere axis (per-bit axis varies)
  - Rotation angle follows diffusive walk: `θ = min(π, √(L/L_char)·π/2)`
  - Characteristic length: `L_char = L₀·(Δn₀/|Δn|)²`, `L₀ = 75 km` at base Δn = 1.2e-7
  - Temperature (`T_coeff = -5e-7/°C`) and bend radius (`0.135·(r_clad/R)²`) modulate Δn → change L_char → scrambling rate
  - Literature: Menyuk & Wai (1994), Wai & Menyuk (1996) for diffusive polarization model
- `analysis/validation/validate_birefringence.py`: rewritten for random model — 6 self-consistency checks (power conservation, zero-length identity, temp/wavelength dependence, seeded reproducibility, polarization variance at long distance)

**System-level demo extended:**
- `analysis/val_system.py`: distance sweep 0–200 km (10 km steps) + 250–1000 km (50 km steps) = 37 points
- Three regimes: 0% QBER (0–80 km) → peak ~68% at 200 km (birefringence-dominated) → decay to ~48% dark-count floor (500–1000 km)
- Effective bit rate: 250 MHz (4 ns window per bit, 4000 samples × 1 ps); cited to Takesue (2007, 1.6 GHz BB84 over 200 km) and Dixon (2008, GHz-clock QKD)

**Manuscript updated:**
- Section 3.3 (Birefringence): describes SU(2) random rotation, diffusive angle, per-bit axis variation
- Section 4 (Birefringence Validation): rewritten for new validation figure (6 panels)
- Section 5 (System-Level Demo): three-regime QBER description, 0–1000 km, literature citations
- Model Limitations #1: phenomenological birefringence noted
- 3 new bibliography entries: Gobby (2004, APL), Takesue (2007, Nature Photonics), Dixon (2008, APL); total 37 references
- Manuscript now 21 pages (up from 18), compiles cleanly

**Output files:**
- `val_system/val_system--seed42.png` (265 KB, 200 DPI) — full 1000 km sweep with bit-rate title
- `val_system/val_system--seed42.csv` — QBER data at all 37 distances
- `analysis/val_birefringence/val_birefringence--seed42.png` — updated 6-panel validation figure
- `analysis/val_birefringence/val_birefringence--seed42.csv` — validation data

---

## 2026-07-12 — LaTeX compilation fix: mdframed → colorbox + scope clarification

### Session: LaTeX fix + project scope broadened

**Background:** `journal_paper_outline.tex` failed with "Not in outer par mode" on every `\begin{table}[!ht]`.

**Root cause:** [`mdframed`](https://ctan.org/pkg/mdframed) patches `\@xfloat` using `\color@vbox`, which was removed from the LaTeX kernel in 2024–2025.

**Fix:**
- Removed `mdframed` + `caption` packages
- Replaced abstract box: `\mdframed{...}` → `\colorbox{highlightblue!5}{...}`
- `journal_paper_outline.pdf` now compiles (13 pages, xelatex)

**Scope clarification:** The project is a general-purpose physical-layer fiber-optic simulator, not QKD-specific. BB84 is one protocol on top. AGENTS.md updated accordingly.

---

## 2026-07-08 — Literature verification fixes + Ulrich bend model

### Session: literature verification (commit 25cf9b9)

| Fix | File | Detail |
|-----|------|--------|
| 1 — Stokes S3 clip | `src/visualization/stokes.py` | Added `np.clip(S3, -1.0, 1.0)` before `np.arcsin()` — prevented NaN at near-circular polarization |
| 2 — Symmetric Jones + standard Δβ | `src/channel/fiber.py` | Switched to `diag(exp(±j·Δβ·L/2))` with Δβ = 2π·Δn/λ (Agrawal [6] Eq 4.1.2). Relative phase = Δβ·L = 2π·L/L_B as before. |
| 3 — PMD sign randomized | `src/channel/fiber.py` | `if np.random.rand() < 0.5: Hx, Hy = Hy, Hx` — 50:50 fast/slow axis per realization |
| 4 — Vπ doc + refs | `src/channel/phase_modulator.py` | Added numbered refs [1] Weis & Gaylord 1985, [2] Alferness 1988; clarified Vπ is MZM push-pull effective value |
| 5 — Fields title | `src/visualization/fields.py` | Changed to `plt.suptitle()` spanning all subplots |
| — Yuan ref corrected | `src/channel/fiber.py` | Vol 27, 2019 → vol 24, no. 2, pp. 1062-1071, 2016 + DOI |

### Session: Ulrich bend model (this session)

**Bend model rebuilt — `num_bends` → `bend_radius`:**

| File | Change |
|------|--------|
| `src/channel/fiber.py` | Replaced `num_bends=0` (Yuan stress rod factor 2.4e-4) with `bend_radius=None` using Δn_bend = 0.135·(r_fiber/R)². References [7] Ulrich 1980, [8] Smith 1980, [9] Shibata 1986. |
| `analysis/validation/validate_birefringence.py` | Sweeps `bend_radius` (2 mm–2 cm), fits Δn vs (r_f/R)². 0.0000% error on slope coefficient 0.135. |
| `analysis/validation/validate_cd.py`, `validate_pmd.py`, `validate_attenuation.py` | `num_bends=0` → `bend_radius=None` |
| `tests/test_fiber.py` | `num_bends=0` → `bend_radius=None` |
| `src/protocols/bb84_*.py` | `num_bends=10` → `bend_radius=None` |
| `literature_verification_report.md` | Updated table, issue summary, recommendations (all fixed). |
| `validation_report.md` | Updated Fix 2 description. |
| `AGENTS.md` | Updated parameters, birefringence description, removed known issue. |

**Cumulative status:** All 5 literature verification issues resolved. All 4 Tier-1 validations pass with 0.0000% error. Bend model blocker eliminated.

---

## 2026-07-02 — Tier 1 Channel Validation

### Session: ~11:00–11:30 UTC+5

**Restructured validation into `analysis/validation/` with per-task output dirs:**
- `analysis/validation/validate_cd.py` → outputs to `analysis/val_cd/`
- `analysis/validation/validate_pmd.py` → outputs to `analysis/val_pmd/`
- `analysis/validation/validate_attenuation.py` → outputs to `analysis/val_att/`
- `analysis/validation/validate_birefringence.py` → outputs to `analysis/val_biref/`

**Task 1 — CD validation (Agrawal [6] Fig 2.6, Eq 3.2.6):**
- 30 ps Gaussian pulse at z/L_D = 0.0, 0.5, 1.0, 2.0
- RMS width ratio vs analytic: **0.0000% error** at all points
- Output: `analysis/val_cd/val_cd--seed42.{png,csv}`

**Task 2 — PMD validation (Razavi [5] Fig 2.11):**
- Fixed `fiber.py:145`: `np.random.rayleigh` → `scipy.stats.maxwell.rvs`
  Old: Rayleigh (2D) gave RMS = 1.128 × target
  New: Maxwell (3D) gives RMS = target exactly
- 20k realizations: RMS DGD = 22.36 ps (target 22.36 ps), mean = 20.62 ps
- KS test p = 0.82 — data consistent with Maxwellian
- All 12 fiber tests continue to pass

**Task 3 — Attenuation validation (Keiser [1] Eq 3.6, SMF-28):**
- Sweep 0–200 km, 0.182 dB/km @ 1550 nm
- **0.0000% error** at all distances

**Task 4 — Birefringence validation (Yuan [4] Fig 1):**
- Physics-based bend birefringence: Δn_bend = 8.762e-10 / R²
- L_B vs R for Δn₀ = {0.5, 1.0, 2.0} × 10⁻⁵
- Confirmed L_B ∝ R / Δn₀ scaling in bend-dominated regime
- **Finding**: fiber.py uses fixed `bend_effect_factor = 2.4e-4` (not R-dependent);
  correct Yuan model requires 1/R² scaling — flagged for future update

**Task 5 — Reproducibility pass:**
- `run_all.py` at project root runs all 4 validations with `--seed N`
- All outputs tagged with seed in filename: `{name}--seed{N}.{ext}`
- Total runtime: ~17.5s for all 4 scripts

## 2026-06-24 — QBER vs distance dispersion graph

### Session: ~16:00–16:10 UTC+5

**New analysis script — `analysis/qber_vs_distance_dispersion.py`:**
- Sweeps fiber length 10–200 km with 5 ps MZM-carved pulses (300 bits/point).
- Two curves: dispersion ON (CD+PMD active) and dispersion OFF (baseline).
- Dispersion OFF: flat 0% QBER at all distances.
- Dispersion ON: QBER climbs from ~0% at 10 km to ~42% at 200 km (CD + PMD accumulate with distance).
- Output: `analysis/qber_vs_distance_dispersion.png`

## 2026-06-24 — BB84 migration to `sample_field` + dispersion test

### Session: ~15:30–16:00 UTC+5

**BB84 scripts now use `sample_field()` (not `instantaneous_field`):**
- `bb84_ideal.py`, `bb84_high_bitrate.py`: `alice_laser.instantaneous_field(normalize=False, over_period=True)` → `alice_laser.sample_field(dt=1e-12, n_samples=1000)`
- This returns the complex envelope over 1 ns (one bit at 1 Gbaud), unblocking chromatic dispersion and PMD.
- `--dispersion` CLI flag added to both scripts (default False for backward compatibility).
- `dt=1e-12` is now passed to `cable()` so that the FFT frequency grid is valid.

**New file — `src/protocols/bb84_test_dispersion.py`:**
- MZM-carved Gaussian pulses for broadband field generation (5–30 ps σ).
- Laser initialized Ey-only so X-cut MZM modulates the entire field.
- MZM biased at V_pi (V=0 → null, V=V_pi → full transmission).
- Higher laser power (+10 dBm) compensates for low pulsed duty cycle.
- `dispersion=True` by default (this is the test file's purpose).
- CLI flags: `--pulse-sigma`, `--short-pulse`, `--no-dispersion`.

**CD/PMD now produces measurable QBER (100 km, seed=42):**
| Pulse σ | dispersion | QBER |
|---|---|---|
| 30 ps | False | 0.00 % |
| 30 ps | True  | 0.00 % (PMD < pulse width) |
| 5 ps  | True  | 15.00 % (z/LD ≈ 87, PMD >> pulse) |

## 2026-06-24 — `get_electric_field` → `instantaneous_field`

### Session: ~15:00–15:10 UTC+5

**Rename:**
- `CWLaser.get_electric_field()` renamed to `CWLaser.instantaneous_field()` to make its purpose unambiguous — returns the full optical field over one ~5 fs period for fast single-bit polarisation/phase validation.
- `instantaneous_field` docstring explicitly warns: NOT for CD, PMD, or baud-rate physics. Use `sample_field()` instead.
- `sample_field` docstring updated to note that `instantaneous_field` is available for quick validation.

**All callers updated:**
- `src/protocols/bb84_ideal.py`, `src/protocols/bb84_high_bitrate.py`
- `analysis/laser_characterization.py`
- `tests/test_cwlaser.py`
- `src/channel/fiber.py` (docstring reference)
- `AGENTS.md`, `README.md`

**Unchanged (legacy, still uses `get_electric_field`):**
- `src/deprecated/sslaser.py` — SolidStateLaser's method kept as-is.
- `main.py`, `main.ipynb`, `scripts/` — all use SolidStateLaser, not CWLaser.

---

## 2026-06-24 — Repository restructure & Tier 0 (testing/reproducibility)

### Session: ~10:30–11:45 UTC+5

**Structural changes:**
- `src/__init__.py` removed — `src` is now a namespace root, not a package (PEP 420).
- `src/opto_eq/` → `src/channel/` — clearer name for optical channel components.
- `src/viewers/` → `src/visualization/` — descriptive name for plotting utilities.
- `src/protocols/examples/` flattened → `src/protocols/` — BB84 scripts are the main protocols, not examples.
- `src/lasers/sslaser.py`, `src/lasers/ndyag.py` → `src/deprecated/` — broken/unused lasers, out of sight but kept.
- `src/lasers/__init__.py` now only exports `CWLaser`.
- 7 root-level loose scripts moved to `scripts/` (except `main.py`).
- `opto-sim.rar` deleted.

**Tier 0 — Testing & reproducibility:**
- `tests/` directory created with `conftest.py` (auto-seeds `random` + `np.random`, `--seed` CLI arg).
- 48 unit tests across 4 files:
  - `test_cwlaser.py` (11): power convention, phase noise, RIN scaling, seeded reproducibility.
  - `test_mzm.py` (13): Vpi, null/peak, quadrature bias, push-pull vs single-drive, insertion loss.
  - `test_fiber.py` (10): attenuation, birefringence, temperature, CD power conservation.
  - `test_apd.py` (11): responsivity, noise scaling, detect_photons, thermal floor.
- `--seed` CLI argument added to `bb84_ideal.py` and `bb84_high_bitrate.py`.
- `pytest>=8.0` added to `requirements.txt`.

---

## 2026-06-04 — Physics-informed detector overhaul & bug fixes

### Session: 11:30–11:45 UTC+5

---

### 1. `src/channel/fiber.py`

| Item | Detail |
|---|---|
| **File** | `src/channel/fiber.py` |
| **Total lines** | 119 (was 101) |
| **Change type** | Edit |

**Change A — Literature sources (lines 1–16)**
- **Old:** Lines 1–4: three bare comment lines (`# Gerd Keiser Book Chapter 3`, `# Behzad Razavi`, `# Thorlabs`)
- **New:** Lines 1–16: formatted literature block citing Keiser [1], Hui [2], Keck [3], Yuan [4], Razavi [5], Agrawal [6]

**Change B — Attenuation formula (lines 105–117)**
- **Old** (was line ~94): `pout = pin / (10**(-attenuation_factor * fiber_length / 10))`
  - Bug: `pin / (10^(-αL/10))` = `pin * 10^(+αL/10)`, so power *increased* with loss.
- **New** (lines 105–117):
  ```python
  att_lin = 10 ** (-attenuation_factor * fiber_length / 10)
  pout = pin * att_lin
  E = E * np.sqrt(att_lin)
  ```
  - Fix: `pout = pin * 10^(-αL/10)` for correct physical loss (Keiser [1] Eq 3.6).
  - Addition: electric field scaled by `sqrt(att_lin)` so `|E_out|² / |E_in|² = att_lin`, keeping the field consistent with the power budget.

---

### 2. `src/detectors/apd.py`

| Item | Detail |
|---|---|
| **File** | `src/detectors/apd.py` |
| **Total lines** | 129 (was 96) |
| **Change type** | Rewrite |

**Complete rewrite of the APD detector class:**

| Aspect | Old | New | Reason |
|---|---|---|---|
| **`__init__` signature** | `(self, wavelength, excess_noise_factor, load_resistance, temperature, gain=12, frequency=3e8/1550e-09, quantum_efficiency=0.9, dark_current=10e-6)` | Removed `frequency` parameter | Frequency is derived from `c/λ` — it is not a user-tunable knob. The old code's `frequency=40` in BB84 scripts was a unit-compensation hack, not physical. |
| **Duplicate constants** | `self.charge = 1.6e-19` and `self.q = 1.602e-19` (both e) | Single `self.charge = 1.602e-19`. `self.q` removed. | Prevent numerical inconsistency. |
| **`detect_photons()`** | `mpn = field_energy/h * frequency` → units of `1/m³`; then `mpn * area * exposure_time` gives `s/m` (not dimensionless). Missing factor of `c`. | `power / (h·ν) * exposure_time * η` where `power` is the actual optical power in Watts. Units: `(W / J) · s = 1` (dimensionless). | Old formula was dimensionally wrong. New formula follows Agrawal [2] Eq 4.1.2 and Saleh & Teich [3] Eq 17.1-10. |
| **`detect_photons()` args** | `(self, field_energy, area, exposure_time=1e-9)` | `(self, power, exposure_time)` | Power passed directly in Watts; no need to derive from field energy density. |
| **`calculate_output_current()`** | `I_signal = self.gain * self.R * power` (correct math), accepted unused `frequency` arg | Removed `frequency` arg. Same formula, now `(self, power)`. | `frequency` was dead code. |
| **`calculate_noise()`** | `I_noise = sqrt(i_d² + i_q² + i_th²) * F` — excess noise factor F applied to *all* noise terms including thermal. | `I_noise = sqrt(F·(i_d² + i_q²) + i_th²)` — F applies only to shot-noise terms (Kasap [1] Eq 4.45). | Physical correction: thermal noise is not multiplied by the avalanche excess noise factor. |
| **`output()`** | `(self, E, area, bandwidth, details=False)` — took electric field array, returned `detected_photons` (int) | `(self, power, bandwidth, area=1, exposure_time=None, details=False)` — takes optical power in Watts, returns noisy `I_total` (float, Amperes) | Current-based output is physically meaningful for receiver design. Includes Gaussian noise realization `I_total = normal(I_signal, I_noise)`. |
| **Photon noise sampling** | Clipped negative counts but didn't add noise to output current | Always returns `I_total` with additive Gaussian noise; photon count still available via `details=True` | Realistic receiver output includes both signal and noise. |
| **Literature citations** | None | Lines 3–9: Kasap [1], Agrawal [2], Saleh & Teich [3] | Traceability for every formula. |

---

### 3. `src/detectors/apd_v2.py`

| Item | Detail |
|---|---|
| **File** | `src/detectors/apd_v2.py` |
| **Change type** | Deleted |

- `apd_v2` was generated by an external coding agent, not by the repository owner.
- It was never exported in `src/detectors/__init__.py` (which only lists `apd`).
- The `bb84_high_bitrate.py` script that imported it was updated to use the regular `apd` instead (see below).

---

### 4. `src/visualization/stokes.py`

| Item | Detail |
|---|---|
| **File** | `src/visualization/stokes.py` |
| **Total lines** | 77 (was 66) |
| **Change type** | Edit |

| Lines | Change |
|---|---|
| 4–10 | Added literature block: Collett [1], Hecht [2], Born & Wolf [3] |
| 27–37 | Reordered computations: `S0` computed first, guard `if S0 == 0` moved before division. Added inline citations to Collett [1] Eq 2.12–2.15 |
| 39–43 | **Removed** lines `psi = 0.5 * arctan2(S2, S1)` and `chi = 0.5 * arcsin(S3/S0)`. The `chi` line was dividing the already-normalized `S3` by `S0` again, causing `|S3/S0| >> 1` → `arcsin(NaN)`. This code was dead (computed but never returned). Replaced with a commented block noting the correct formula and the clipping needed to prevent floating-point NaNs. |

---

### 5. `src/visualization/fields.py`

| Item | Detail |
|---|---|
| **File** | `src/visualization/fields.py` |
| **Total lines** | 37 (was 29) |
| **Change type** | Edit |

| Lines | Change |
|---|---|
| 1–6 | Added `import numpy as np` and `import matplotlib.pyplot as plt` (were missing, causing `NameError` at runtime). Added literature citation Hecht [1]. |

---

### 6. `src/lasers/cwlaser.py` (NEW FILE)

| Item | Detail |
|---|---|
| **File** | `src/lasers/cwlaser.py` |
| **Total lines** | 199 |
| **Change type** | Created |

New physics-informed continuous-wave laser model for QKD simulation:

| Feature | Implementation | Source |
|---|---|---|
| Phase noise | Wiener process with diffusion coefficient `D_φ = 2π·Δν` | Henry [1] Eq 18 |
| RIN | White noise low-pass filtered at `rin_bandwidth` | Kikuchi [5], Coldren [2] Ch. 5 |
| Polarization | Jones vector from azimuth `ψ` and ellipticity `χ` | Yariv [3] Ch. 6 |
| Power | User-specified in dBm, converted to Watts | Steady-state CW model |
| Field output | `get_electric_field(t, over_period, normalize)` — same interface as `SolidStateLaser` | Backward compatibility |

**Literature sources in file headers:**
- [1] Henry, C. H., "Theory of the Linewidth of Semiconductor Lasers", IEEE JQE 1982
- [2] Coldren, Corzine & Mashanovitch, "Diode Lasers and Photonic Integrated Circuits", 2nd ed., Wiley 2012
- [3] Yariv, A., "Optical Electronics", 4th ed., Saunders 1991, Ch. 6
- [4] Schawlow & Townes, "Infrared and Optical Masers", Phys. Rev. 1958
- [5] Kikuchi, K., "Characterization of semiconductor-laser phase noise", Opt. Express 2012

---

### 7. `src/lasers/__init__.py`

| Item | Detail |
|---|---|
| **File** | `src/lasers/__init__.py` |
| **Total lines** | 5 (was 4) |
| **Change type** | Edit |

| Lines | Change |
|---|---|
| 3, 5 | Added `from .cwlaser import CWLaser` and `'CWLaser'` to `__all__` |

---

### 8. `src/protocols/examples/bb84_ideal.py`

| Item | Detail |
|---|---|
| **File** | `src/protocols/examples/bb84_ideal.py` |
| **Total lines** | 146 (was 135) |
| **Change type** | Edit |

| Lines | Change |
|---|---|
| 30–33 | Removed `frequency=40` from `apd()` constructor call (APD no longer accepts this parameter) |
| 71–72 | Added `pout = alice_laser.power_out * 1e-3` — convert SolidStateLaser's mW output to Watts for the detector chain |
| 81–85 | `E, _ = cable(...)` → `E, pout = cable(...)` — now captures attenuated power from cable (in Watts) |
| 96–120 | **Rewrote detection section:** |
| 97–106 | New: derive PBS arm power from field ratio × `pout` (calibration-independent) |
| 108–110 | `detector.output(power=power_x, bandwidth=1e6)` — new APD signature |
| 112–120 | Old: ratio-based comparison of photon counts (`if abs(px-py)/(px+py) > 0.001`). New: differential detection with 3σ noise floor threshold (`if I_x > threshold or I_y > threshold: bob_bit = 0 if I_x > I_y else 1`) |

---

### 9. `src/protocols/examples/bb84_high_bitrate.py`

| Item | Detail |
|---|---|
| **File** | `src/protocols/examples/bb84_high_bitrate.py` |
| **Total lines** | 132 (was 135 initially, then rewrites) |
| **Change type** | Edit |

| Lines | Change |
|---|---|
| 6 | `from src.detectors.apd_v2 import APD_v2` → `from src.detectors import apd` (apd_v2 was deleted) |
| 20–23 | Removed `frequency=40` from `apd()` constructor call |
| 53–55 | Added `pout = alice_laser.power_out * 1e-3` — mW → W conversion |
| 66–69 | `E, _ = cable(...)` → `E, pout = cable(...)` — captures attenuated power |
| 82–105 | **Rewrote detection section** matching `bb84_ideal.py`: power from field ratio; current-based differential detection with 3σ threshold |

---

### 10. Generated graphs

| File | Description | Command |
|---|---|---|
| `analysis/qber_vs_distance.png` | QBER vs fiber length (10–200 km), dispersion off, corrected APD | `py run_experiments.py` (with dispersion patch) |
| `analysis/qber_vs_bitrate.png` | QBER vs detector bandwidth (10 MHz–10 GHz), 100 km fiber, corrected APD | `py run_bitrate_experiments.py` (with dispersion patch) |

---

### Summary of bugs fixed

| Bug | File | Impact | Root cause |
|---|---|---|---|
| Output power increased with fiber loss | `fiber.py:96` | All QKD results were wildly wrong | `pin / 10^(-αL/10)` instead of `pin * 10^(-αL/10)` |
| Dark count attribute `self.darkcount` undefined | `apd.py:80` | `details=True` would raise `AttributeError` | Typo: should be `self.dcr` |
| Missing imports | `fields.py:1` | `plot_field()` would raise `NameError` | `import numpy` / `pyplot` were absent |
| `arcsin(S3/S0)` returned NaN | `stokes.py:29` | RuntimeWarning whenever Stokes were computed | `S3` already normalized; `S3/S0` exceeded 1 due to FP noise |
| `frequency` used as unit-compensation hack | `apd.py:20` | Physical meaning of parameter was ambiguous | Missing `c` in `detect_photons()` forced a workaround |
| Excess noise factor applied to thermal noise | `apd.py:68` | Noise current overestimated | F should only scale shot-noise terms (Kasap Eq 4.45) |
| Photon detection units broken | `apd.py:31–33` | Photon count had units of `s/m` | Missing factor of `c` in `detect_photons()` |
| Electric field not scaled by attenuation | `fiber.py:94` | Detector received pre-attenuation field amplitude | Field was never multiplied by `sqrt(att_lin)` |

---

## 2026-06-04 — Field-as-power convention (removed separate pin/pout tracking)

### Session: 13:30–13:45 UTC+5

**Design rationale:** The E-field is the single source of truth for both polarization and optical power. Previously, `fiber.py` and `apd.py` tracked power as a separate float (`pin`/`pout`), creating a dichotomy with `channel/` which already derives power from the field. Now, power is always derived from `mean(|E|²)`.

**Field convention:** `mean(|E|²)` = optical power in **Watts**. The field is calibrated once at laser output in each script (e.g., BB84 scripts scale `E` so that `mean(|E|²) = laser.power_out * 1e-3`).

---

### 1. `src/channel/fiber.py` — cable() returns field only

| Lines | Change |
|---|---|
| 17–18 | Removed `pin` from `cable()` signature. Function is now `cable(fiber_length, E, dispersion, attenuation_factor, temperature, num_bends, pm_dispersion)` |
| 105–117 | Removed `pout = pin * att_lin`. Attenuation applied directly to field: `E = E * sqrt(att_lin)`. Convention: `mean(|E_out|²) = att_lin · mean(|E_in|²)`. Return value is just `E` (was `E, pout`). |

### 2. `src/detectors/apd.py` — output() takes E-field

| Lines | Change |
|---|---|
| 83–84 | `output(self, power, ...)` → `output(self, E, ...)`. Power derived internally: `power = mean(|E|²)`. Convention assumes mean(|E|²) = power in Watts (caller must pre-calibrate). |

### 3. `src/protocols/examples/bb84_ideal.py`

| Lines | Change |
|---|---|
| 70–72 | Added field calibration: `power_W = laser.power_out * 1e-3; E *= sqrt(power_W / mean(|E|²))` |
| 82–85 | `E, pout = cable(... pin=pout ...)` → `E = cable(...)` |
| 99–106 | Removed `field_power`, `power_x`, `power_y` computation (was deriving power from field × pout ratio). PBS arms `Ex, Ey` passed directly to detector. |
| 109–110 | `detector.output(power=power_x, ...)` → `detector.output(E=Ex, ...)` |

### 4. `src/protocols/examples/bb84_high_bitrate.py`

| Lines | Change |
|---|---|
| 53–55 | Same field calibration as bb84_ideal |
| 66–69 | `E, pout = cable(... pin=pout ...)` → `E = cable(...)` |
| 83–86 | Removed field_power/power_x/power_y computation |
| 89–90 | `detector.output(power=power_x, ...)` → `detector.output(E=Ex, ...)` |

### 5. `main.py`

| Lines | Change |
|---|---|
| 28 | Removed `pout = source.power_out` (was only used as cable input) |
| 60 | `E, _ = cable(100, E, pout, dispersion=True)` → `E = cable(100, E, dispersion=True)` |

### Files affected (this session)

| File | Lines changed | Change type |
|---|---|---|
| `src/channel/fiber.py` | 17–18, 105–117 | Edit |
| `src/detectors/apd.py` | 83–84 | Edit |
| `src/protocols/examples/bb84_ideal.py` | 70–72, 82–85, 99–110 | Edit |
| `src/protocols/examples/bb84_high_bitrate.py` | 53–55, 66–69, 83–90 | Edit |
| `main.py` | 28, 60 | Edit |

---

## 2026-06-04 — CWLaser RIN: relaxation-oscillation model (Coldren Eq 5.3.38)

### Session: 14:00–14:30 UTC+5

**`src/lasers/cwlaser.py`** — Completely rewrote `_sample_rin()` and `_generate_rin()`:

| Aspect | Old | New | Reason |
|---|---|---|---|
| **RIN model** | White noise → 1st-order LP filter at `rin_bandwidth` | White noise → shape by sqrt(S_RIN(f)) from linearized rate equations (Coldren [2] Eq 5.3.38) | The old model had no physical basis — missed the relaxation oscillation resonance entirely. |
| **Parameters** | `rin_bandwidth` (removed) | `relaxation_frequency`, `damping_rate` | f_RO and γ are physically meaningful; rin_bandwidth was arbitrary. |
| **Implementation** | Time-domain IIR (first-order) | Frequency-domain: rFFT → shape → irFFT | Avoids bilinear-transform instability at optical sampling rates. |
| **RIN PSD** | Flat LP | `S_RIN(f) ∝ (γ²+ω²) / | ω_R²-ω²+jγω|²` | Correct resonance: flat below f_RO, peak at f_RO, 1/f² roll-off above. |
| **Citations** | Kikuchi [5] only | Coldren [2] §5.3.3 (Eq 5.3.38), Petermann [5] Ch. 7 | Traceable physics. |

**Bug fixed:** `H_sq` was the raw `num/den` without normalization by the DC value, causing H_sq values ~10⁻²¹ and making the RIN output have ~0 variance. Fixed by `H_sq = H_sq_raw / H_sq_raw[0]`.

**Other improvements:**
- Added `_rin_dt_min = 1/(10·f_RO)` to decouple RIN time resolution from optical-time sampling. When `dt` is < `_rin_dt_min`, RIN is generated at the coarser rate and interpolated, preventing FFT artifacts at extreme sample rates (~10¹⁷ Hz for optical periods).
- Single-sample `get_electric_field(t, over_period=False)` now uses `_sample_rin` for the RIN value instead of `np.random.normal(0, sqrt(RIN_lin * 1e9))` which had a hardcoded 1 GHz bandwidth.

---

## 2026-06-04 — Laser characterisation script with eye diagrams

### Session: 15:00–15:15 UTC+5

**Created `analysis/laser_characterization.py`** — comprehensive verification suite for the CWLaser:

| Feature | Method | Verification |
|---|---|---|
| **Power convention** | `mean(|E|²)` vs `_power_w` | Bar chart showing error < 1% |
| **Optical linewidth** | Complex-envelope PSD via Welch, Lorentzian fit (curve_fit) | FWHM from fit vs specified Δν |
| **Phase noise** | Structure function `D_φ(τ) = ⟨[φ(t+τ)-φ(t)]²⟩` | Slope `2π·Δν` for Wiener process (Henry [1]) |
| **RIN spectrum** | Welch PSD of `δP(t)` vs Coldren Eq 5.3.38 | Resonance peak at `f_RO`, DC level at `RIN_0` |
| **Polarisation** | Stokes parameters via `compute_stokes_parameters()`; Poincaré sphere via existing `poincare()` | ψ, χ match laser settings |
| **Eye diagram (NRZ-OOK)** | PRBS → intensity modulation → RIN + phase noise → direct detection | Eye closure at increasing bitrates |

**Output files:**
| File | Content |
|---|---|
| `analysis/laser_characterization.png` | 8-panel dashboard |
| `analysis/poincare_sphere.png` | Poincaré sphere (via `stokes.py:poincare()`) |
| `analysis/eye_diagrams.png` | 5/10/25 Gbaud eyes side-by-side |

**Note on OOK modulation:** A phase modulator alone cannot produce OOK (amplitude modulation). The script uses an idealised intensity-modulator model (`E_mod = E_cw · wfm`). Real OOK requires direct laser current modulation or a Mach-Zehnder modulator.

---

## 2026-06-05 — CWLaser API redesign: `sample_field()`, MZM device model, characterisation rewrite

### Session: 08:30–10:00 UTC+5

**Three interconnected changes:**

| Change | Files | Rationale |
|---|---|---|
| `CWLaser.sample_field(dt, n_samples)` | `src/lasers/cwlaser.py` | New public method returning complex-envelope field (n_samples, 2) with all physical effects (power, phase noise, RIN, polarisation). The laser owns the physics — no more ad-hoc noise generation in characterisation scripts. |
| `get_electric_field` API cleanup | `src/lasers/cwlaser.py` | `t=0` → `dt=1e-12` (descriptive). Hardcoded `1000` → parameter `n_samples=1000`. Backward-compatible for all existing callers that use `over_period=True`. |
| MZM physical model | `src/channel/mzm.py` (new), `src/channel/__init__.py` | Push-pull MZM: `E_out = E_in·cos(π·V/V_pi)·exp(j·π·V_bias/V_pi)`. Replaces idealised `E·wfm` intensity modulator with a physically correct interferometric model (Agrawal [1] §4.2). |
| Characterisation rewritten | `analysis/laser_characterization.py` | Removed `_field_complex_envelope` / `_field_series` helpers. All plots now use `laser.sample_field()`. Eye diagrams use `MZM`. Phase noise and RIN measured from full field output (end-to-end verification). Added `matplotlib.use('Agg')` for headless operation. |

**Impact on call sites:**
- `get_electric_field(dt=..., over_period=True, n_samples=...)` — all existing callers use keyword args `over_period=True, normalize=False` which are unchanged.
- `get_electric_field(t=...)` for single-sample → now `get_electric_field(dt=...)`. No existing callers use single-sample mode, so no breakage.
- `sample_field()` is the recommended API going forward for both characterisation and BB84 scripts.
- MZM is importable as `from src.channel.mzm import MZM` or `from src.channel import MZM`.

### BB84 scripts updated to CWLaser — `src/protocols/examples/bb84_ideal.py` and `bb84_high_bitrate.py`

| Change | Rationale |
|---|---|
| `SolidStateLaser` → `CWLaser` | CWLaser provides physics-informed RIN, phase noise, and correct power scaling. |
| Removed `E *= sqrt(power_W / mean(|E|²))` calibration | CWLaser's `get_electric_field()` already has `mean(|E|²) = P_W` — no rescaling needed. |
| Removed unused `eve_laser` (bb84_ideal) | Dead code — variable was defined but never referenced in the simulation loop. |
| `dispersion=True` → `dispersion=False` (bb84_ideal) | The FFT-based dispersion function (`apply_dispersion` in fiber.py:53) computes `f = (1/100)*(1/(D·L·0.2e-9))` which produces ~10¹⁷ Hz — unphysical and corrupts the 193 THz CWLaser field. This is a pre-existing blocked issue documented in AGENTS.md. The high-bitrate script already used `dispersion=False`. |

**Verified:** 0% QBER at 10–200 km (bb84_ideal), 0% QBER at 1 MHz–5 GHz bandwidth / 100 km (bb84_high_bitrate).

---

## 2026-06-08 — Physics-based MZM rewrite; Vpi calibration fix in BB84

### Session

| Change | Files | Rationale |
|---|---|---|
| MZM rewritten as MZI + PhaseModulator | `src/channel/mzm.py` | MZM now internally uses `PhaseModulator` instances per arm. Y-branch splitter/combiner model with configurable insertion loss and extinction ratio. Supports push-pull (zero chirp) and single-drive (residual chirp) modes. `V_pi` derived from crystal parameters — no more empirical `cos(π·V/V_pi)` with a user-specified V_pi. |
| Hardcoded `Vpi = 3.757` removed from BB84 scripts | `bb84_ideal.py`, `bb84_high_bitrate.py` | Both scripts now derive V_pi from `pm_alice.Vpi` at runtime. The stale hardcoded value (3.757 V) differed by 3.2 % from the PhaseModulator's crystal-computed Vpi (3.8826 V), causing mismatched-basis QBER to drift to ~38 % instead of 50 %. With the fix: sifted QBER = 0 %, total QBER = 25 %. |
| Characterisation script updated | `analysis/laser_characterization.py` | `MZM(V_pi=5.0)` → `MZM()` (uses default PhaseModulator). Docstring updated for new switching voltage convention. |

**Verified:** `analysis/laser_characterization.py` runs clean, all three output files generated.

---

## 2026-06-09 — FFT-based chromatic dispersion, birefringence/PMD fixes in fiber.py

### Session

| Change | Files | Rationale |
|---|---|---|
| Chromatic dispersion rewritten with FFT | `src/channel/fiber.py` | Old `apply_dispersion()` used `f = (1/100)*(1/(D·L·0.2e-9))` producing ~10¹⁷ Hz. Replaced with `H(Ω) = exp(-j·β₂·Ω²·L/2)` via `np.fft.fftfreq` (Agrawal [6] §2.4). Applied to both Ex/Ey (CD is isotropic). Verified against Gaussian pulse: broadening ratio error < 0.06 % at 0.5–2.0× LD. |
| Unit fix in GVD calculation | `src/channel/fiber.py` | D stored as `17e-12` (ps/(nm·km)) but β₂ formula needs s/m². Added conversion `D_SI = D × 1e-6`. Previously β₂ was 6 orders too small. |
| Birefringence Jones matrix fixed | `src/channel/fiber.py` | Removed spurious `del_T = pmd_sd²` factor that zeroed out the beat-length phase (~10⁻¹⁶ rad instead of the correct ~10⁶ rad). Now uses `Δβ = 4π·Δn/λ`, Jones = `diag(exp(j·Δβ·L/2), 1)`, preserving the L/2 beat-length convention (SM vs PM fibre discrimination). |
| PMD rewritten (frequency-domain) | `src/channel/fiber.py` | Old `apply_pmd()` added random per-sample phase (phase noise, not PMD). Replaced with frequency-domain DGD: Maxwellian-distributed Δτ, Jones matrix `diag(exp(∓j·ω·Δτ/2), exp(±j·ω·Δτ/2))`. |
| New parameters | `src/channel/fiber.py` | Added `dt` (required for dispersion), `wavelength` (no longer hardcoded). Default `attenuation_factor` changed from 0.25 → 0.182 dB/km (SMF-28 at 1550 nm). |
| Updated `main.py` | `main.py` | Added `dt=1e-12` to `cable()` call. |

**Verified:** `analysis/laser_characterization.py`, `bb84_ideal.py` (0 % sifted QBER), `bb84_high_bitrate.py` (0 % sifted QBER). Gaussian pulse broadening matches Agrawal theory.
