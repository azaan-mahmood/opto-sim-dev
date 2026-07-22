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
- **Time-bin phase-encoding BB84 implemented**: pulsed laser, AsymmetricMZI, time-bin protocol, Gobby replication validation
- **77 tests all pass** (48 original + 8 pulsed laser + 21 AMZI)
- **Gobby 2004 replication**: 0 km QBER = 3.19% (paper: ~3.3%) — excellent baseline match
- **Manuscript updated**: 28 pages, 7 validated components, 38 references — new AsymmetricMZI, SPAD, time-bin BB84 validation sections
- **New files**: `src/channel/interferometer.py`, `src/protocols/bb84_time_bin.py`, `tests/test_interferometer.py`, `analysis/val_gobby/validate_gobby.py`
- **Performance**: ~2600 pulses/s on 3rd Gen i5 — 2.5M pulses in ~16 min, 25M pulses in ~2.7h
- **Next**: Long-distance statistics need more pulses; consider increasing mu or using batched simulation for speed

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

### Fiber (`src/channel/fiber.py`)
- `propagate()` no longer takes `pin` or returns `pout`
- Signal impairments applied in order: birefringence → chromatic dispersion → PMD → attenuation
- **Attenuation**: `E *= sqrt(10^(-α·L/10))` (Keiser [1] Eq 3.6), default 0.182 dB/km (SMF-28 at 1550 nm)
- **Birefringence**: Hybrid dual-model dispatch (automatic based on fibre length):
  - **Short fibres** (L < `SECTIONAL_LIMIT` = 2 km, `model='sectional'`): Multi-section ordered product of random-axis SU(2) matrices. `N = round(L / section_length)` sections, each with independent random axis and phase `Δβ·Δz = 2π·|Δn|·Δz/λ`. Physically correct beat length (L_B ≈ 31 m, Δn₀ = 5×10⁻⁸) with temperature (`T_coeff = -3×10⁻⁹/°C`) and bend (Ulrich [7]) modulation, stochastic residual, and clamping. Suitable for DV-QKD and DPS QKD. Performance scales as O(N) — fast at default section_length=100 m.
  - **Long fibres** (L >= 2 km, `model='phenomenological'`): Single SU(2) rotation with angle `θ = min(π, √(L/L_char)·π/2)` (Menyuk & Wai 1994). `L_char = L₀·(Δn₀/|Δn|)²` with L₀ = 75 km, Δn₀ = 0.87×10⁻⁵, `T_coeff = -5×10⁻⁷/°C`. Produces gradual distance-dependent QBER for long-haul BB84. Switched automatically because the multi-section model converges to uniform SU(2) within ~1 km regardless of parameters.
  - Dispatch: `apply_birefringence(model='auto')` or `propagate(model='auto')`; override with `'sectional'` or `'phenomenological'`. `SECTIONAL_LIMIT = 2000 m`.
- **Chromatic dispersion** (disabled by default): FFT-based, `H(Ω) = exp(-j·β₂·Ω²·L/2)` applied to both Ex and Ey (Agrawal [6] Eq 2.4.11). Requires `dt` (sampling interval) and assumes the field is the complex envelope.
- **PMD** (disabled by default): Frequency-domain DGD with Maxwellian-distributed differential group delay (Razavi [5]). Applied alongside CD.
- Parameters: `fiber_length (km)`, `E`, `dt` (required for cd/pmd), `wavelength`. Impairments independently toggled: `birefringence`, `cd`, `pmd`, `attenuation` (all bool, defaults True/None/None/True). Legacy `dispersion` flag sets both `cd` and `pmd` when not explicitly provided. `attenuation_factor`, `temperature`, `bend_radius`, `pm_dispersion`, `section_length`, `model`.

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

### Fiber Birefringence — Hybrid dispatch
- Two models:
  - **Sectional** (L < 2 km): first-principles multi-section, L_B ≈ 31 m (Agrawal §4.1)
  - **Phenomenological** (L >= 2 km): `θ = min(π, √(L/L_char)·π/2)`, L₀ = 75 km
- Auto-dispatch via `model='auto'`; override with `model='sectional'` or `'phenomenological'`
- Old symmetric Jones matrix `diag(exp(±j·Δβ·L/2))` superseded (retained in git history)
- Bend model `Δn_bend = 0.135·(r_fiber/R)²` (Ulrich [7]) shared by both models

### BB84 Scripts
- `bb84_ideal.py` / `bb84_high_bitrate.py`: use CWLaser → `sample_field()` → polarizer → phase modulator → propagate (dispersion flag) → PBS → APD
- Both accept `--dispersion` CLI flag (default False for backward compatibility in ideal/bitrate scripts)
- `bb84_test_dispersion.py`: MZM-carved Gaussian pulses (5–30 ps σ) for broadband field generation; dispersion=True by default
- `bb84_duplinskiy.py`: SPAD-based BB84 example (Duplinskiy et al. 2017). 0 km back-to-back QBER validates SPAD + PM + PBS chain (2.6% vs paper ~2%).
- QBER with CW laser (1 MHz linewidth) is 0% regardless of dispersion flag — near-monochromatic field is CD/PMD-agnostic
- QBER with 5 ps MZM-carved pulses + dispersion at 100 km: **15.0 %** (z/LD ≈ 87, PMD >> pulse width)

## Generated Graphs
- `analysis/val_cd/val_cd--seed42.png`: CD validation, Agrawal Fig 2.6 (0.0000% error)
- `analysis/val_pmd/val_pmd_dgd--seed42.png`: PMD DGD histogram vs Maxwellian (Razavi Fig 2.11)
- `analysis/val_att/val_att--seed42.png`: Attenuation vs distance SMF-28 (0.0000% error)
- `analysis/val_biref/val_biref--seed42.png`: Birefringence L_B vs R (Yuan Fig 1) [old model]
- `analysis/val_birefringence/val_birefringence--seed42.png`: Random-axis birefringence validation (6 panels)
- `val_system/val_system--seed42.png`: System-level demo: QBER vs distance 0–1000 km, 250 MHz bit rate

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

### Current coverage (48 tests, all passing)
| File | Tests | What they check |
|---|---|---|
| `tests/test_cwlaser.py` | 11 | Power convention, `sample_field` shape, phase noise scaling, polarisation vector, seeded reproducibility, `instantaneous_field` shapes, `power_out`, zero-linewidth edge case, RIN scaling |
| `tests/test_mzm.py` | 13 | `V_pi` derivation, null/peak transmission, quadrature bias, transfer symmetry, push-pull vs single-drive, mode validation, X-cut vs Y-cut `V_pi`, array modulation, insertion loss, crystal-cut component selection |
| `tests/test_fiber.py` | 10 | Attenuation formula and distance scaling, birefringence unitarity and phase shift, temperature sensitivity, dispersion `dt` requirement and power conservation, output shape, zero-length edge case, wavelength dependence, seeded reproducibility |
| `tests/test_apd.py` | 11 | Responsivity formula, signal current scaling, zero-power, bandwidth scaling of shot noise, thermal noise floor, `output()` return type and details dict, `detect_photons` edge cases, DCR formula, seeded reproducibility |
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
├── scripts/                   # Experiment & diagnostic scripts
│   ├── compare_dc_rf.py
│   ├── laser_characterization.py   # SolidStateLaser characterisation (deprecated)
│   ├── laser_diagnostic.py
│   ├── ndyag_characterization.py
│   ├── run_experiments.py
│   └── run_bitrate_experiments.py
├── run_all.py                 # Runs all Tier-1 validations with --seed N
├── analysis/                  # Generated figures & TeX papers
│   ├── validation/            # Tier-1 channel validation scripts
│   │   ├── validate_cd.py     # CD: Agrawal Fig 2.6, 0.0000% error
│   │   ├── validate_pmd.py    # PMD: Razavi Fig 2.11, Maxwellian DGD
│   │   ├── validate_attenuation.py  # Attenuation: SMF-28, 0.0000% error
│   │   └── validate_birefringence.py # Birefringence: Yuan Fig 1, L_B vs R
│   ├── val_cd/                # CD outputs (seed-tagged PNG+CSV)
│   ├── val_pmd/               # PMD outputs
│   ├── val_att/               # Attenuation outputs
│   ├── val_biref/             # Birefringence L_B vs R (old model)
│   ├── val_birefringence/     # Random-axis birefringence validation (6-panel)
│   ├── val_mzm/               # MZM validation outputs

│   ├── val_system.py          # System-level demo
│   ├── laser_characterization.py   # Active: CWLaser dashboard (Agg, headless)
│   ├── qber_vs_distance_dispersion.py  # Dispersion QBER sweep
│   ├── *.png
│   └── *.tex / *.pdf
├── tests/                     # pytest suite (48 tests)
│   ├── conftest.py            # --seed CLI, auto-seed at session start
│   ├── test_apd.py
│   ├── test_cwlaser.py
│   ├── test_fiber.py
│   └── test_mzm.py
├── src/
│   ├── lasers/
│   │   ├── __init__.py        # Exports: CWLaser
│   │   └── cwlaser.py         # Active laser model
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
│   │   ├── mzm.py
│   │   ├── optics.py
│   │   └── phase_modulator.py
│   ├── visualization/         # Plotting & Poincaré (renamed from viewers)
│   │   ├── __init__.py
│   │   ├── fields.py
│   │   ├── polarimeter.py
│   │   └── stokes.py
│   ├── protocols/
│   │   ├── __init__.py
│   │   ├── bb84_ideal.py           # CW laser, optional dispersion
│   │   ├── bb84_high_bitrate.py    # Bitrate-sweep variant
│   │   ├── bb84_test_dispersion.py # MZM-carved pulses, dispersion=on by default
│   │   └── bb84_duplinskiy.py      # Duplinskiy et al. replication (SPAD, VOA)
│   └── common/                # README images only
├── research_roadmap.md
├── AGENTS.md
├── CHANGELOG.md
├── pyproject.toml
├── requirements.txt
├── val_system/               # System demo outputs (PNG+CSV)
└── README.md
```

## Files Changed (recent sessions — most recent first)
| File | Change |
|---|---|
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
