# Project Context: opto-sim — Physical-Layer Fiber-Optic Simulator
# CRITICAL RULES - MUST FOLLOW

## Responses
- Keep responses concise and to the point - unless the user asks for more elaboration

## PLANNING MODE

- Always ask clarifying questions
- Never assume design, tech stack or features
- Use deep-dive sub-agents to assist with research where possible
- All physics based implementations MUST be literature proof and cross-checked. Never assume and hallucinate to prove something
- If an implementation does not have concrete literature, provide the user with this info concisely but do not remove crucial information 
- All literature sources must be mentioned
- Use deep-dive sub-agents to review the different aspects of your plan before presenting to the user where possible

## CHANGE / EDIT MODE

- Use the best model for the task
- After completing a feature, always run commands like lint, type check and next build to check for code quality
- After completing a physics-based feature, always double check with the literature
- When using sub-agents to implement features, act as a co-ordinator only
- When using sub-agents to implement physics-based features, act as a supervisor


## Summary
- **332 tests all pass** (per-file counts in the Test Suite table)
- **Gobby replication COMPLETE (GOBBY-1..7e)** — zero fitted parameters. Link budget taken from the paper: α = 0.2 dB/km, μ = 0.1 photons/clock, η_Bob = 0.045 (incl. Bob's 5 dB), gate 3.5 ns, clock 2 MHz, pulse 80 ps, P_e = 8.5e-7/clock/detector (dark 3.2e-7 + stray 5.3e-7, deliberately lumped — their sum is what was measured). Nine-point sweep (`--seed 42 --target-sifted 3000`): floor **3.088 ± 0.247 % vs the paper's 3.3 %** (afterpulse = 0, §19.5; modulation error `e_mod` on the PhaseModulator, GOBBY-6), 122 km **9.50 ± 0.43 % vs 8.9 %** (two-clock drift defect fixed, GOBBY-7d), χ²/dof = 1.95 with p = 0.099 — **never significant** (§28.3). Visibility is an *output* (`V = S/(S + 2·P_e)`), not an input. OPEN-1/2/3/5 all closed
- **Duplinskiy chain brought to standard (DUPL-1)** — O(1) sifting (five O(N) lists removed), 8-outcome precompute (29,938 → 386,245 pulses/s, 12.9×, bit-identical), detector parameters (η/DCR/dead time/afterpulse) no longer literals, `cd`/`pmd` exposed. 50 km: **2.07 ± 0.40 %** (4e6 pulses, 1,258 sifted) vs the paper's ~2 % — but this is **occurrence #7**: one mechanism (afterpulsing) produces a total the paper decomposes into three (afterpulse + extinction + drift). A1 (p_ap = 0.05 net of dead time) is load-bearing for 100 % of the error budget; DUPL-2 (extinction term, spec at §29) is the discriminating test. Drift/recalibration deliberately **not** implemented until A1 settles
- **A4 SETTLED TRUE by algebra (§28.6a, 2026-08-11, zero runs)**: `circular_analyser` (QWP+PBS) ≡ the paper's PC3 (≈HWP 22.5°)+PBS readout — `J_circ = diag(1,−i)·H·diag(1,−i)`, port statistics agree state-for-state. `cd`/`pmd` on this chain are **inert by construction** (single time sample, §27.3)
- **VALID-1 validation architecture (§28)**: three levels of claim — Level 1 replication (Gobby ✅; Duplinskiy total ✅ but mechanism misattributed), Level 2 coverage, Level 3 two-sided constraint (null in one protocol + swing in another, same unchanged code). **Level 3 holds for zero mechanisms today**: the birefringence pair is one properly-powered sweep away (§28.8 #1 — the Duplinskiy swing side is a 2–50 sifted-bit probe, unquotable); the linewidth pair is blocked on DPS-QKD; the CD/PMD pair on the time-resolved field (measured free to run, 0.004 % of a sweep row)
- **Next work (project priority §28.8)**: ① birefringence sweep on Duplinskiy at ≥3,000 sifted/row (first Level-3 pair), ② A4 ✅ done, ③ DUPL-2 extinction (spec §29), ④ time-resolved field, ⑤ DPS-QKD, ⑥ polarisation-paper survey (§28.5)
- **Key files**: `analysis/val_gobby/validate_gobby.py` (GOBBY sweep: ALPHA_dB=0.2, ETA_BOB=0.045, P_E=8.5e-7, AFTERPULSE_PROB=0.0, `--target-sifted 3000`); `src/analytic/gobby_model.py` (analytic reference: `QBER(L) = e_mod + (1−V)/2`, `T_INT`); `src/protocols/bb84_time_bin.py` (PERF-2 8-outcome table, ~10 µs/pulse); `src/protocols/bb84_duplinskiy.py` (DUPL-1 standard); `opto-sim-issues-and-fixes.md` (authoritative issue/status log — §§18–29 are the current passes); `paper/` (paper PDFs, untracked)
- **Performance**: Gobby chain ~10 µs/pulse (PERF-2 memoisation), 10M-pulse point ≈ 1.7 min; Duplinskiy chain 386,245 pulses/s

## Goal
Open-source, validated physical-layer fiber-optic simulator where the complex-envelope electric field is the single source of truth. QKD (BB84) is one application protocol; the platform is general-purpose: CW lasers, MZMs, CD, PMD, birefringence, attenuation, APDs — any classical or quantum optical link that can be modeled as a linear pipeline of physically parametrised impairments.

## Field Convention
- `mean(|E|²)` = optical power in **Watts**
- Field is calibrated once at the laser output in each script
- All downstream components (optics, fiber, detector) derive power from the field — no separate `pin`/`pout` tracking
- `Ez` component is excluded by design (standard fiber-optic simulators do the same)

## Active Components

### CW Laser (`src/lasers/cwlaser.py`)
- Steady-state model with Wiener phase noise (Henry [1]) and relaxation-oscillation RIN (Coldren Eq 5.3.38)
- Parameters: `wavelength`, `power_dbm`, `linewidth`, `rin_density`, `polarization_azimuth`, `polarization_ellipticity`, `relaxation_frequency` (default 5 GHz), `damping_rate` (default 1.88e10 rad/s)
- RIN PSD: `S_RIN(f) = RIN_0 · (γ² + 4π²f²) / |(2πf_RO)² - (2πf)² + jγ(2πf)|²`
- Frequency-domain implementation (rFFT → shape → irFFT) with coarser internal time resolution when optical sampling is extreme
- **`sample_field(dt, n_samples)`** — public method returning complex-envelope field (n_samples, 2) with power, phase noise, RIN, and polarisation. Primary API for downstream components.
- `instantaneous_field(dt=1e-12, over_period=False, n_samples=1000, normalize=True)` — single-period (~5 fs) field for fast polarisation/phase validation; NOT for CD/PMD/baud-rate physics.
- BB84 scripts use `sample_field(dt=1e-12, n_samples=1000)` — complex envelope unlocks CD/PMD.
- `instantaneous_field` is reserved for quick per-bit polarisation/phase checks.

### APD (`src/detectors/apd.py`)
- Current-based output: `output(E, bandwidth)` derives `power = mean(|E|²)` from the field
- Kasap Eq 4.45: excess noise factor `F` applied **only** to shot noise, not thermal
- `detect_photons()` uses `power/(h·ν) · t · η` (Agrawal Eq 4.1.2, Saleh & Teich Eq 17.1-10)
- No `frequency` parameter — optical frequency derived from `c/λ`

### SPAD (`src/detectors/spad.py`)
- Geiger-mode single-photon detector (inherits APD): dead time, DCR, afterpulsing, gated detection; ID230 specs (η = 0.10, dead time 13 µs, DCR 15 Hz, p_ap = 0.05, gate 20 ns in the Duplinskiy chain)
- **GOBBY-7c**: click probability corrected to the Poisson form `1 − exp(−η·μ)` (was `η·(1 − exp(−μ))`) — affects every protocol chain

### Fiber (`src/channel/fiber.py`)
- `propagate()` no longer takes `pin` or returns `pout`
- Signal impairments applied in order: birefringence → chromatic dispersion → PMD → attenuation
- **Attenuation**: `E *= sqrt(10^(-α·L/10))` (Keiser [1] Eq 3.6), default 0.182 dB/km (SMF-28 at 1550 nm)
- **Birefringence**: single multi-section model at all lengths (the former phenomenological model was deleted in the 5th pass, PHYS-5): ordered product of random-axis SU(2) matrices. `N = round(L / section_length)` sections, each with independent random axis and phase `Δβ·Δz = 2π·|Δn|·Δz/λ`. Physically correct beat length (L_B ≈ 31 m, Δn₀ = 5×10⁻⁸) with temperature (`T_coeff = -3×10⁻⁹/°C`) and bend (Ulrich [7]) modulation, stochastic residual, and clamping. Quasi-static: `FiberRealization` draws the matrix once per fibre. Performance scales as O(N) via pairwise tree reduction — ~0.16 ms/apply at 122 km (N = 2440).
- **Chromatic dispersion** (disabled by default): FFT-based, `H(Ω) = exp(-j·β₂·Ω²·L/2)` applied to both Ex and Ey (Agrawal [6] Eq 2.4.11). Requires `dt` (sampling interval) and assumes the field is the complex envelope.
- **PMD** (disabled by default): Frequency-domain DGD with Maxwellian-distributed differential group delay (Razavi [5]). Applied alongside CD.
- Parameters: `fiber_length (km)`, `E`, `dt` (required for cd/pmd), `wavelength`. Impairments independently toggled: `birefringence`, `cd`, `pmd`, `attenuation` (all bool, defaults True/None/None/True). Legacy `dispersion` flag sets both `cd` and `pmd` when not explicitly provided. `attenuation_factor`, `temperature`, `bend_radius`, `pm_dispersion`, `section_length`, `model`.
- **`FiberRealization.birefringence_matrix()`** — accessor returning the quasi-static Jones matrix (copy; None if disabled). Basis of the Duplinskiy compensation (`U_comp = J⁻¹`, ARCH-3).

### Mach-Zehnder Modulator (`src/channel/mzm.py`)
- Physics-based MZI built from `PhaseModulator` instances per arm
- Y-branch splitter/combiner model with insertion loss and extinction ratio
- Two electrode configurations:
  - **push-pull** (default): both arms modulated with opposite voltages, zero chirp
  - **single-drive**: one arm modulated, other is passive reference; residual frequency chirp (Koyama & Iga [2])
- `V_pi` derived from the internal PhaseModulator's crystal parameters (`pm.Vpi`)
- `switching_voltage = V_pi` (voltage for ON→OFF); `bias_voltage` shifts operating point (Vb=Vpi/2 → quadrature, 50 % transmission)
- Transfer function (per polarisation component, X-cut modulates Ey):
  - push-pull: `E_out ∝ cos(π·(V+V_bias)/(2·V_pi))`
  - single-drive: `E_out ∝ exp(j·π·(V+V_bias)/(2·V_pi)) · cos(π·(V+V_bias)/(2·V_pi))`
- References: Agrawal [1] §4.2, Koyama & Iga [2], Weis & Gaylord [3]

### PhaseModulator (`src/channel/phase_modulator.py`)
- X-cut LiNbO3 crystal; `V_pi` from material parameters (Zelmon, Small & Jundt 1997 — indices verified to 0.5 %, V_π to 1.5 %; Δn off 18 % and inert)
- GOBBY-6/7: carries the modulation error (`phase_error_rad`), arm-length drift (`phase_drift_rad_s`) and static bias (`bias_offset_v`) that reproduce Gobby's 3.3 % floor and 122 km point — all from the paper, nothing fitted

### Interferometers (`src/channel/interferometer.py`)
- **PolarizationMultiplexedAMZI** (Gobby topology, GOBBY-3/4): encoder 50:50 split fixed; `coupler_split` real convention; `T_INT` is the polarisation-multiplexed transmission and must NOT be applied to the balanced topology
- **AsymmetricMZI** generalised (delay, phase) — the basis for the recommended next protocol, DPS-QKD (§27.9)

### Optics (`src/channel/optics.py`)
- `circular_analyser(E)` — QWP(45°)+PBS circular-basis analyser (`J = (1/√2)[[1,−i],[−i,1]]`); **proven equivalent (§28.6a) to the Duplinskiy paper's PC3(≈HWP 22.5°)+PBS readout**: `J_circ = diag(1,−i)·H·diag(1,−i)`, with the λ/4 basis-selection phase at PM2 in the paper and inside the fixed analyser in the chain
- `pbs(E)` — true H/V projector (PHYS-6; the phase-sensitive element is `circular_analyser`)
- `coupler_combine` — real convention (verified unreachable before the flip, 944c24e); interference is **cos**, not sin
- `halfwave` / `quarterwave` / `polarization_rotator` / `polarization_controller` / `voa` / `polarizer` / `hadamard` — no absolute-phase prefactors (RETARDER PHASE CONVENTION in module header)

### Stokes viewer (`src/visualization/stokes.py`)
- `compute_stokes_parameters(E)` returns `[S0,S1,S2,S3], [psi, chi]`
- S1,S2,S3 normalized by S0; S0 set to 1.0
- `chi = 0.5·arcsin(S3)` (S3 already normalized, S0=1)

### Chromatic Dispersion (`src/channel/fiber.py`)
- **FIXED**: FFT-based model `H(Ω) = exp(-j·β₂·Ω²·L/2)` via `np.fft.fftfreq` (Agrawal [6] Eq 2.4.11)
- Verified against Gaussian pulse broadening: ratio error < 0.06 % at z = 0.5–2.0× LD
- Requires `dt` (sampling interval) — callers must pass this for `dispersion=True`
- PMD uses frequency-domain DGD (Maxwellian) applied alongside CD

## Known Issues / Blockers

### SolidStateLaser (`src/lasers/sslaser.py`)
- ODE rate equations never achieve population inversion because `Rp = 1/τ₂` (transparency threshold). Needs `Rp > 1/τ₂` by at least 2-10×.
- `power_dbm` only sets initial photon density `I_0`, doesn't affect pump — changing it changes none of the dynamics.
- Electric field units are mW-based (not Watts). BB84 scripts calibrate with `E *= sqrt(power_W / mean(|E|²))` at laser output.
- Retained in `src/deprecated/` for reference; CWLaser is the active replacement.

### LaTeX compilation: `mdframed` + `\begin{table}` (FIXED)
- `mdframed` patches `\@xfloat` using `\color@vbox`, which is no longer defined in the current LaTeX kernel (2024+). Causes "Not in outer par mode" error on every `\begin{table}`.
- **Fix:** Removed `mdframed` package; replaced abstract box with `\colorbox{highlightblue!5}{\begin{minipage}{...}}`.
- Also removed `caption` package (was only needed for `\captionof`, no longer used).
- `journal_paper_outline.pdf` now compiles at 13 pages with xelatex.

### Fiber Birefringence — single multi-section model (PHYS-5)
- One model at all lengths: multi-section ordered product of random-axis SU(2) matrices, L_B ≈ 31 m (Agrawal §4.1), L_c = 50 m default (Menyuk & Wai [10])
- `model='auto'` ≡ `model='sectional'`; `model='phenomenological'` raises ValueError (removed 5th pass — fitted `θ = min(π, √(L/L_char)·π/2)`, L₀ = 75 km, Δn₀ = 0.87e-5 had no literature backing and no speed advantage; 0.16 ms/apply at 122 km)
- Ensemble mean polarization saturates to uniform SU(2) within ~1 correlation cell (per-section retardance δ ≈ 10 rad wraps mod 2π)

### Gobby Validation (Time-bin BB84) — GOBBY-1..7e COMPLETE
- `analysis/val_gobby/validate_gobby.py` runs QBER vs distance (0–122 km, 9 points), `--seed 42 --target-sifted 3000`
- Link budget from the paper, nothing fitted: α = 0.2 dB/km, μ = 0.1, η_Bob = 0.045 (incl. Bob's 5 dB), gate 3.5 ns, clock 2 MHz, pulse 80 ps, P_e = 8.5e-7/clock/detector (dark 3.2e-7 + stray 5.3e-7, lumped)
- `AFTERPULSE_PROB = 0.0` (GOBBY-7b; §19.5 — afterpulsing is not in Gobby's budget; at 2 MHz gating with 13 µs dead time the detector is off while trapped carriers release)
- `VISIBILITY` is an **output** (V = S/(S + 2·P_e)); the 3.3 % floor comes from the PhaseModulator modulation error (`e_mod`), the 122 km point from the link budget + drift (`phase_drift_rad_s`/`bias_offset_v`, GOBBY-7/7d)
- Results: floor 3.088 ± 0.247 % vs 3.3 %; 122 km 9.50 ± 0.43 % vs 8.9 %; χ²/dof = 1.95, p = 0.099 (never significant, §28.3)
- Analytic reference: `src/analytic/gobby_model.py` — `QBER(L) = e_mod + (1−V)/2`, `V = S/(S + 2·P_e)`, `S = μ_eff·10^(−αL/10)·η_Bob`; `T_INT` applied (polmux transmission)
- `cd`/`pmd`/birefringence are an **exact null** on this chain (§26.6 — time-bin immunity, demonstrated bit-identically at field-perturbation level)

### Duplinskiy (Polarisation BB84) — DUPL-1 done; A1 and DUPL-2 open
- `src/protocols/bb84_duplinskiy.py` — SPAD chain, 45° input → PM1 (encode) → fibre → U_comp (inverse Jones, default ON, `--no-compensation` control) → PM2 (basis) → `circular_analyser` → 2× SPAD
- 50 km: 2.07 ± 0.40 % (4e6 pulses, 1,258 sifted) vs paper ~2 % — but **occurrence #7**: afterpulsing alone produces a total the paper decomposes into three (afterpulse + extinction + drift). The "0.98 % vs 2 %" gap was underpowered noise (140 sifted @400k)
- **A1 load-bearing**: `afterpulse_prob = 0.05` assumed net of dead-time suppression — 100 % of the error budget. DUPL-2 (extinction term, spec §29) is the discriminating test
- Drift/recalibration cycle (§27.5) **not implemented** — the paper's ~2 % average over an 80 % duty cycle is citable but would double-count until A1 settles
- `cd`/`pmd` exposed (DUPL-1) but **inert by construction** — single time sample (§27.3); the time-resolved field is the only route to a CD/PMD pair
- A4 settled TRUE by algebra (§28.6a) — the readout (`circular_analyser` vs the paper's PC3+PBS) is not an open question

### Open work (project priority, §28.8)
1. **Birefringence sweep** on Duplinskiy at ≥3,000 sifted/row — completes the first Level-3 pair (one run; current probe is 2–50 sifted bits, unquotable)
2. A4 ✅ done (2026-08-11)
3. **DUPL-2** extinction term (spec §29: 3×3 matrix, decision rules 1–4, 8 acceptance criteria)
4. **Time-resolved field** — free to run (0.004 % of a sweep row); only route to the CD/PMD Level-3 pair
5. **DPS-QKD** — the linewidth Level-3 pair; reuses `AsymmetricMZI`
6. **Polarisation-paper survey** (§28.5) — pre-empts the cherry-picking objection

### BB84 Scripts
- `examples/bb84_ideal.py` / `examples/bb84_high_bitrate.py` (legacy demos, moved to `src/protocols/examples/`): use CWLaser → `sample_field()` → polarizer → phase modulator → propagate (dispersion flag) → PBS → APD
- Both accept `--dispersion` CLI flag (default False for backward compatibility in ideal/bitrate scripts)
- `bb84_test_dispersion.py`: MZM-carved Gaussian pulses (5–30 ps σ) for broadband field generation; dispersion=True by default
- `bb84_duplinskiy.py`: SPAD-based BB84 replication (Duplinskiy et al. 2017). Polarization compensation (inverse J, default ON; `--no-compensation` control). 50 km: **2.07 ± 0.40 %** vs paper ~2 % (1,258 sifted @4e6; the older 0.98 % figure was underpowered). Uncompensated control: 25.9 % @ 58 sifted (probe only — full sweep pending, §28.8 #1). Detector constants parameterised (DUPL-1); `cd`/`pmd` exposed but inert (single time sample).
- QBER with CW laser (1 MHz linewidth) is 0% regardless of dispersion flag — near-monochromatic field is CD/PMD-agnostic
- QBER with 5 ps MZM-carved pulses + dispersion at 100 km: **15.0 %** (z/LD ≈ 87, PMD >> pulse width)

## Generated Graphs
- `analysis/val_cd/val_cd--seed42.png`: CD validation, Agrawal Fig 2.6 (0.0000% error)
- `analysis/val_pmd/val_pmd_dgd--seed42.png`: PMD DGD histogram vs Maxwellian (Razavi Fig 2.11)
- `analysis/val_att/val_att--seed42.png`: Attenuation vs distance SMF-28 (0.0000% error)
- `analysis/val_biref/val_biref--seed42.png`: Birefringence L_B vs R (Yuan Fig 1) [old model]
- `analysis/val_birefringence/val_birefringence--seed42.png`: Random-axis birefringence validation (6 panels)
- `analysis/val_gobby/val_gobby--seed42.png` + `val_gobby_table.tex`: Gobby nine-point sweep, `--seed 42 --target-sifted 3000` (GOBBY-1..7e)
- `val_system/val_system--seed42.png`: System-level time-bin BB84 demo (ARCH-1 rebuild): QBER vs distance/pulse width/visibility/μ/DCR, SPAD path
- `val_system/val_system_scenarios--seed42.{csv,tex}`: BLOCK-2 impairment table, `--seed 42 --target-sifted 3000` per scenario (OPEN-3 closed; time-bin immunity demonstrated, §26.6)

- `analysis/qber_vs_distance.png`: 0% QBER 10-190 km
- `analysis/qber_vs_bitrate.png`: 0% at 215 MHz → 35% at 10 GHz
- `analysis/qber_vs_distance_dispersion.png`: QBER climb 0→42% at 10→200 km with dispersion (5 ps pulse)
- `analysis/laser_characterization.png`: 8-panel CWLaser dashboard
- `analysis/poincare_sphere.png`: Poincaré sphere from Stokes parameters
- `analysis/eye_diagrams.png`: NRZ-OOK eye diagrams at 5/10/25 Gbaud

## Test Suite (Tier 0)

### Running tests
```bash
python -m pytest tests/ -v                    # all tests
python -m pytest tests/ -v --seed=123         # custom RNG seed
python -m pytest tests/test_cwlaser.py -v     # single file
```

### Current coverage (332 tests, all passing)
| File | Tests | What they check |
|---|---|---|
| `tests/test_analytic_gobby.py` | 30 | GOBBY-2 analytic reference: residuals within 0.05 pp of §19.3, stated visibilities (0.9906/0.8840), sensitivity rows, 165.8 km out-of-sample prediction, `T_INT` handling |
| `tests/test_apd.py` | 11 | Responsivity formula, signal current scaling, zero-power, bandwidth scaling of shot noise, thermal noise floor, `output()` return type and details dict, `detect_photons` edge cases, DCR formula, seeded reproducibility |
| `tests/test_cwlaser.py` | 20 | Power convention, `sample_field` shape, phase noise scaling, polarisation vector, seeded reproducibility, `instantaneous_field` shapes, `power_out`, zero-linewidth edge case, RIN scaling |
| `tests/test_fiber.py` | 59 | Attenuation formula and distance scaling, birefringence (unitarity, phase shift, Ulrich bend law, sectional model, matrix accessor), temperature sensitivity, dispersion `dt` requirement and power conservation, output shape, zero-length edge case, wavelength dependence, seeded reproducibility |
| `tests/test_interferometer.py` | 37 | AMZI delay/phase response, visibility, polarisation behaviour, power conservation; AsymmetricMZI generalisation (GOBBY-3) |
| `tests/test_mzm.py` | 13 | `V_pi` derivation, null/peak transmission, quadrature bias, transfer symmetry, push-pull vs single-drive, mode validation, X-cut vs Y-cut `V_pi`, array modulation, insertion loss, crystal-cut component selection |
| `tests/test_optics.py` | 68 | VOA attenuation, polariser/rotator/WP behaviour, beam-splitter power balance, `circular_analyser` unitarity, PBS phase-blindness (PHYS-6 regression), `coupler_combine` real convention (GOBBY-4/944c24e) |
| `tests/test_phase_modulator.py` | 28 | Phase response, Vπ convention, `sample_field` shapes; modulation error / drift / bias knobs (GOBBY-6/7) |
| `tests/test_polmux_interferometer.py` | 31 | PolarizationMultiplexedAMZI (GOBBY-3): 50:50 encoder split, PBC, transmission `T_INT`, interference as cos |
| `tests/test_protocols.py` | 10 | Duplinskiy compensation-flag regression (uncompensated QBER 49.6% vs compensated 1.5% @10 km μ=2); DUPL-1 chain standard (bit-identity, negative control) |
| `tests/test_spad.py` | 25 | Dead time, DCR, afterpulsing, gated detection, ID230 specs; Poisson click form `1−exp(−η·μ)` (GOBBY-7c) |
| `tests/conftest.py` | — | Auto-seeds `random` and `np.random` at configured start; `--seed` CLI option (default 42); `rng_seed` fixture |

### Seed convention
- Each test file seeds `np.random` before RNG-dependent tests.
- `conftest.py` pins both `random` and `np.random` at session start (`pytest_configure`).
- `--seed` on the CLI overrides the default (42).
- BB84 scripts accept `--seed` for per-run reproducibility.

## Directory Structure (post-refactor)

```
opto-sim-dev/
├── main.py                    # Quick demo (uses deprecated SolidStateLaser)
├── run_all.py                 # Runs all Tier-1 validations with --seed N
├── scripts/                   # Experiment & diagnostic scripts
│   ├── compare_dc_rf.py
│   ├── laser_characterization.py   # SolidStateLaser characterisation (deprecated)
│   ├── laser_diagnostic.py
│   ├── ndyag_characterization.py
│   ├── run_experiments.py
│   └── run_bitrate_experiments.py
├── analysis/                  # Generated figures & validation scripts
│   ├── validation/            # Tier-1 channel validation scripts
│   │   ├── validate_cd.py     # CD: Agrawal Fig 2.6, 0.0000% error
│   │   ├── validate_pmd.py    # PMD: Razavi Fig 2.11, Maxwellian DGD
│   │   ├── validate_attenuation.py  # Attenuation: SMF-28, 0.0000% error
│   │   └── validate_birefringence.py # Birefringence: sectional model (6-panel)
│   ├── val_cd/                # CD outputs (seed-tagged PNG+CSV)
│   ├── val_pmd/               # PMD outputs
│   ├── val_att/               # Attenuation outputs
│   ├── val_biref/             # Birefringence L_B vs R (old model)
│   ├── val_birefringence/     # Random-axis birefringence validation (6-panel)
│   ├── val_mzm/               # MZM validation outputs
│   ├── val_apd/  val_attenuation/  val_cwlaser/   # component outputs
│   ├── val_gobby/             # Gobby sweep (GOBBY-1..7e): validate_gobby.py + CSV/PNG/TeX, --seed 42 --target-sifted 3000
│   ├── val_system.py          # System-level time-bin BB84 demo (SPAD path)
│   ├── val_system_scenarios.py  # Impairment-table generator (BLOCK-2, CSV+TeX)
│   ├── laser_characterization.py   # Active: CWLaser dashboard (Agg, headless)
│   ├── qber_vs_distance_dispersion.py  # Dispersion QBER sweep
│   ├── *.png
│   └── *.tex / *.pdf
├── tests/                     # pytest suite (332 tests)
│   ├── conftest.py            # --seed CLI, auto-seed at session start
│   ├── test_apd.py
│   ├── test_cwlaser.py
│   ├── test_fiber.py
│   ├── test_interferometer.py
│   ├── test_mzm.py
│   ├── test_optics.py
│   ├── test_phase_modulator.py
│   ├── test_polmux_interferometer.py  # PolarizationMultiplexedAMZI (GOBBY-3)
│   ├── test_protocols.py       # Duplinskiy compensation + DUPL-1 regressions
│   ├── test_spad.py
│   └── test_analytic_gobby.py  # GOBBY-2 analytic reference
├── src/
│   ├── lasers/
│   │   ├── __init__.py        # Exports: CWLaser
│   │   └── cwlaser.py         # Active laser model
│   ├── analytic/
│   │   ├── __init__.py
│   │   └── gobby_model.py     # GOBBY-2 analytic reference (QBER(L), T_INT)
│   ├── deprecated/            # Broken/unused, kept for reference
│   │   ├── __init__.py        # Exports: SolidStateLaser, NdYAGLaser
│   │   ├── sslaser.py
│   │   └── ndyag.py
│   ├── detectors/
│   │   ├── __init__.py
│   │   ├── apd.py
│   │   └── spad.py              # Geiger-mode SPAD (ID230 specs)
│   ├── channel/               # Optical channel components (renamed from opto_eq)
│   │   ├── __init__.py
│   │   ├── fiber.py
│   │   ├── interferometer.py   # AMZI + PolarizationMultiplexedAMZI (GOBBY-3/4)
│   │   ├── mzm.py
│   │   ├── optics.py           # circular_analyser, pbs, WPs, voa, couplers (real convention)
│   │   └── phase_modulator.py  # phase_error_rad / phase_drift_rad_s / bias_offset_v
│   ├── visualization/         # Plotting & Poincaré (renamed from viewers)
│   │   ├── __init__.py
│   │   ├── fields.py
│   │   ├── polarimeter.py
│   │   └── stokes.py
│   ├── protocols/
│   │   ├── __init__.py
│   │   ├── bb84_time_bin.py        # Time-bin phase-encoding BB84 (active, PERF-2)
│   │   ├── bb84_test_dispersion.py # MZM-carved pulses, dispersion=on by default
│   │   ├── bb84_duplinskiy.py      # Duplinskiy et al. replication (SPAD, DUPL-1 standard)
│   │   └── examples/
│   │       ├── bb84_ideal.py       # Legacy CW-based BB84 demo
│   │       └── bb84_high_bitrate.py# Legacy bitrate-sweep variant
│   └── common/                # README images only
├── paper/                     # Source-paper PDFs (untracked): Duplinskiy 2017, etc.
├── val_system/                # System demo + scenario outputs (PNG+CSV+TeX, seed-tagged)
├── research_roadmap.md
├── opto-sim-issues-and-fixes.md   # authoritative issue/status log (§§18-29 current passes; NOT committed)
├── AGENTS.md
├── CHANGELOG.md
├── pyproject.toml
├── requirements.txt
└── README.md
```

## Files Changed (recent sessions — most recent first)
| File | Change |
|---|---|
| `opto-sim-issues-and-fixes.md` | **A4 (2026-08-11)**: §28.6a SETTLED TRUE by algebra, zero runs — `J_circ = diag(1,−i)·H·diag(1,−i)` proves `circular_analyser` (QWP+PBS) ≡ the Duplinskiy paper's PC3 (≈HWP 22.5°)+PBS readout; register row → SETTLED; §27.5/§28.8/§29.2/§29.8 downstream edits |
| `CHANGELOG.md` | 2026-08-11: A4 settlement entry (docs-only) |
| `src/protocols/bb84_duplinskiy.py` | DUPL-1: O(1) sifting (five O(N) lists removed), 8-outcome precompute (29,938 → 386,245 pulses/s, 12.9×, bit-identical), detector parameters (η/DCR/dead time/afterpulse) no longer literals |
| `src/detectors/spad.py` | GOBBY-7c: click probability corrected to Poisson form `1 − exp(−η·μ)` (was `η·(1 − exp(−μ))`) — affects every protocol chain |
| `src/channel/phase_modulator.py` | GOBBY-6/7: `phase_error_rad` / `phase_drift_rad_s` / `bias_offset_v` knobs — reproduce Gobby's 3.3 % floor and 122 km point from the paper, nothing fitted |
| `src/channel/interferometer.py` | GOBBY-3: `PolarizationMultiplexedAMZI` (50:50 encoder split fixed, `T_INT` polmux transmission) + `AsymmetricMZI` generalisation (delay, phase) |
| `src/channel/optics.py` | GOBBY-4: `coupler_combine` real convention (interference is cos; unreachability of the pre-flip form verified, 944c24e) |
| `src/analytic/gobby_model.py` | GOBBY-2 Step 1 (NEW): analytic reference implementation — `QBER(L) = e_mod + (1 − V_fringe)/2`, `V_fringe = S/(S + 2·P_e)`, `S = μ_eff·10^(−αL/10)·η_Bob`; all parameters from the paper, nothing fitted (μ_eff = 0.0793 via visibility inversion, §19.2); `p_e`/`mu_eff` overridable for §19.11 and sensitivity rows |
| `src/analytic/__init__.py` | GOBBY-2 Step 1 (NEW): package exports |
| `tests/test_analytic_gobby.py` | GOBBY-2 Step 1 (NEW, 19 tests): §19.3 four residuals within 0.05 pp, stated visibilities (0.9906/0.8840), sensitivity rows (0.26/0.52/2.01 pp), acceptance bar ≤0.30 pp, 165.8 km out-of-sample prediction (V = 0.729) |
| `analysis/val_gobby/validate_gobby.py` | GOBBY-1 (10th pass): link budget re-parameterised to paper values — `ALPHA_dB` 0.182→0.2 (and now actually passed to the MC), `ETA`→`ETA_BOB` 0.10→0.045, `GATE_WIDTH` 1→3.5 ns, `REP_RATE` 2.5→2 MHz, `PULSE_WIDTH` 100→80 ps, `DCR`→`P_E` 8.5e-7/clock, `VISIBILITY` 0.934→1.0 (output), `--p-e`/`--afterpulse` flags, `PILOT_BITS` 2M, `CEILING` 1e9; provenance headers in artifacts; sweep: slope +4.99 pp/100 km, mean |res| 1.22 pp, 0 fitted params |
| `src/protocols/bb84_time_bin.py` | PERF-2 (7th pass): 8-outcome gate-power table built once per point (12 `modulate` calls), per-pulse loop = lookup + SPAD MC; bitwise-equal to the old chain, ~40× faster (~10 µs/pulse) |
| `analysis/val_gobby/validate_gobby.py` | OPEN-1 (7th pass): `VISIBILITY = 0.934` (cited to Gobby's 3.3 % floor) actually passed to `simulate_bb84_time_bin`; `--visibility` CLI flag; V=0.934+no-noise → 3.36 % at 0 km |
| `analysis/val_system_scenarios.py` | BLOCK-2 (6th pass): NEW impairment-table generator — 8 explicit configs @ 100 km, recorded seed, CSV + `.tex` with script/seed/commit in caption; `simulate_point()` gained independent birefringence/attenuation/cd/pmd toggles |
| `src/channel/fiber.py` | PHYS-5 (5th pass): phenomenological model deleted (`_build_jones_phenomenological`, `_apply_birefringence_phenomenological`, dead `_random_su2_rotation*`); `SECTIONAL_LIMIT` removed; `model='auto'` ≡ `'sectional'` at all lengths; `'phenomenological'` raises ValueError |
| `analysis/validation/validate_birefringence.py` | Panels B–D rewritten on the single multi-section model: long-distance ensemble plateau (uniform SU(2) within ~200 m), temperature & bend sensitivity at 2 m (single correlation cell, visible δ(T)/δ(R)) |
| `tests/test_fiber.py` | 5 phenomenological tests → long-distance sectional (power at 100 km, T/λ/seed dependence, output variation); dispatch test → `auto ≡ sectional at all lengths` + `phenomenological raises` |
| `src/protocols/bb84_duplinskiy.py` | PHYS-5: `--birefringence-model` choices reduced to auto/sectional. Bug fix: `compensate` flag now actually disables compensation (was a silent no-op) |
| `tests/test_protocols.py` | NEW: compensation-flag regression (uncompensated QBER 49.6% vs compensated 1.5% @10 km μ=2) |
| `analysis/val_system.py` | ARCH-1: full rebuild on time-bin SPAD chain (CWLaser→MZM carve→encoder AMZI→fiber→decoder AMZI→2×SPAD); linearity shortcut (~5 s/point); panels A distance, B pulse σ, C visibility, D μ, E DCR; floor ≈ 6.4% @75 km (V 3.3% + jitter 0.8% + afterpulse 1.5% + double-click 0.5%) |
| `src/channel/fiber.py` | ARCH-3: `FiberRealization.birefringence_matrix()` accessor (copy of quasi-static J, None if disabled) |
| `src/protocols/bb84_duplinskiy.py` | ARCH-3: polarization compensation via inverse J (default ON, `--no-compensation` control); 50 km QBER 0.98% (1M) vs paper ~2% |
| `tests/test_fiber.py` | +5 `TestBirefringenceMatrixAccessor` tests (unitarity, match-apply, roundtrip, disabled→None, quasi-static) |
| `val_system/val_system--seed42.png` | ARCH-1 figure regenerated: 2×3 grid (panels A–E + notes), 1M pulses/point |
| `opto-sim-issues-and-fixes.md` | Section 13: ARCH-1 (SPAD rebuild) + ARCH-3 (compensation) resolved |
| `CHANGELOG.md` | Fourth pass: ARCH-1/ARCH-3 entry |
| `README.md` | Third pass: fixed Gobby path, counts, exact invocation table |
| `analysis/val_gobby/val_gobby--seed42.png` | Regenerated with 10M pulses (was 200k). Clean curve beyond 80 km. |
| `analysis/val_gobby/validate_gobby.py` | Run with 10M pulses — numbers updated in manuscript. |
| `paperwork/tables/val_gobby_table.tex` | Updated with 10M-pulse results (23546 → 154 sifted bits, 2.43% → 4.55% QBER). |
| `paperwork/manuscript.tex` | Gobby section: removed "statistical noise" language, added systematic offset discussion. |
| `src/detectors/spad.py` | NEW: Geiger-mode SPAD, inherits from APD. Dead time, DCR, afterpulsing, gated detection. ID230 specs. |
| `src/detectors/__init__.py` | Added `spad` export |
| `src/channel/optics.py` | Added `voa(E, attenuation_dB)` for variable optical attenuation |
| `src/protocols/bb84_duplinskiy.py` | SPAD-based BB84 example. 0 km back-to-back validated (2.6% vs paper ~2%). |
| `paperwork/manuscript.tex` | Updated system-level scenarios, removed old model references |
| `paperwork/tables/val_system_table.tex` | 9-row system impairment table |
| `src/channel/fiber.py` | Hybrid dispatch: `apply_birefringence()` auto-routes to `_sectional` (L < 2 km) or `_phenomenological` (L >= 2 km). Fixed duplicate return, section_length=1.0 default, propagate() indentation. Added `model` param to propagate(). |
| `analysis/validation/validate_birefringence.py` | Removed unused `L0 = 75e3`. Validates both models via auto-dispatch. |
| `analysis/val_system.py` | Updated outputs (QBER: 0%→68% at 200 km, dark-count floor ~45% at 1000 km) |
| `AGENTS.md` | Hybrid birefringence description, dispatch threshold, updated files changed |
| `CHANGELOG.md` | Added hybrid dispatch entry |
| `src/channel/fiber.py` | New `_random_su2_rotation()` + updated `apply_birefringence()` — random-axis SU(2) rotation, diffusive angle, per-bit varying axis |
| `analysis/validation/validate_birefringence.py` | Dual-model validation: 13 self-consistency checks (sectional + phenomenological + auto-dispatch) |
| `analysis/val_system.py` | NEW: system-level demo with 0–1000 km sweep, 250 MHz bit rate |
| `paperwork/manuscript.tex` | Updated Sections 3.3/4/5; 3 new refs; 21 pages |
| `val_system/val_system--seed42.png` | Demo figure with full 1000 km sweep |
| `analysis/val_birefringence/val_birefringence--seed42.png` | Updated 6-panel validation figure |
| `journal_paper_outline.tex` | Removed `mdframed` + `caption` packages; replaced abstract box with `\colorbox{highlightblue!5}{...}`. Compiles to 13 pages (was failing with "Not in outer par mode"). |
| `AGENTS.md` | Updated files changed, known issues. |
| `CHANGELOG.md` | Added mdframed fix + journal outline compilation entry. |
| `src/channel/fiber.py` | `num_bends` → `bend_radius` (Ulrich 1980 bend model); symmetric Jones + standard Δβ; PMD sign randomized; refs [7][8][9] added |
| `analysis/validation/validate_birefringence.py` | Sweeps `bend_radius` (2 mm–2 cm); fits Δn vs (r_f/R)²; 0.0000% error |
| `analysis/validation/validate_cd.py` | `num_bends=0` → `bend_radius=None` |
| `analysis/validation/validate_pmd.py` | `num_bends=0` → `bend_radius=None` |
| `analysis/validation/validate_attenuation.py` | `num_bends=0` → `bend_radius=None` |
| `tests/test_fiber.py` | `num_bends=0` → `bend_radius=None` |
| `src/protocols/bb84_ideal.py` | `num_bends=10` → `bend_radius=None` |
| `src/protocols/bb84_high_bitrate.py` | `num_bends=10` → `bend_radius=None` |
| `src/protocols/bb84_test_dispersion.py` | `num_bends=10` → `bend_radius=None` |
| `src/visualization/stokes.py` | Added `np.clip(S3, -1.0, 1.0)` before arcsin (Fix 1) |
| `src/visualization/fields.py` | Changed to `plt.suptitle()` (Fix 5) |
| `src/channel/phase_modulator.py` | Added refs [1][2]; clarified Vπ convention (Fixes 4/6/7) |
| `literature_verification_report.md` | All issues marked FIXED; bend validation updated |
| `validation_report.md` | Fix 2 updated for Ulrich bend model |
| `AGENTS.md` | Updated params, known issues, files changed |
| `CHANGELOG.md` | Added literature verification + bend model entries |
