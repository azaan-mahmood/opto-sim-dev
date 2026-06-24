# Project Context: opto-sim QKD Simulation
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


## Goal
Refactor the QKD simulation codebase with physics-informed models for APD detection, fiber attenuation, and laser sources.

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
- `cable()` no longer takes `pin` or returns `pout`
- Signal impairments applied in order: birefringence → chromatic dispersion → PMD → attenuation
- **Attenuation**: `E *= sqrt(10^(-α·L/10))` (Keiser [1] Eq 3.6), default 0.182 dB/km (SMF-28 at 1550 nm)
- **Birefringence**: Beat-length Jones matrix `[[exp(j·Δβ·L/2), 0], [0, 1]]` where `Δβ = 4π·Δn/λ` and Δn depends on temperature and bends (Yuan [4]). Beat length `L_B = λ/Δn` distinguishes SM (long L_B) vs PM (short L_B) fibre.
- **Chromatic dispersion** (disabled by default): FFT-based, `H(Ω) = exp(-j·β₂·Ω²·L/2)` applied to both Ex and Ey (Agrawal [6] Eq 2.4.11). Requires `dt` (sampling interval) and assumes the field is the complex envelope.
- **PMD** (disabled by default): Frequency-domain DGD with Maxwellian-distributed differential group delay (Razavi [5]). Applied alongside CD.
- Parameters: `fiber_length (km)`, `E`, `dt` (required for dispersion), `wavelength` (default 1550e-9), `dispersion`, `attenuation_factor`, `temperature`, `num_bends`, `pm_dispersion`.

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

### BB84 Scripts
- `bb84_ideal.py` / `bb84_high_bitrate.py`: use CWLaser → `sample_field()` → polarizer → phase modulator → cable (dispersion flag) → PBS → APD
- Both accept `--dispersion` CLI flag (default False for backward compatibility in ideal/bitrate scripts)
- `bb84_test_dispersion.py`: MZM-carved Gaussian pulses (5–30 ps σ) for broadband field generation; dispersion=True by default
- QBER with CW laser (1 MHz linewidth) is 0% regardless of dispersion flag — near-monochromatic field is CD/PMD-agnostic
- QBER with 5 ps MZM-carved pulses + dispersion at 100 km: **15.0 %** (z/LD ≈ 87, PMD >> pulse width)

## Generated Graphs
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
├── analysis/                  # Generated figures & TeX papers
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
│   │   └── apd.py
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
│   │   └── bb84_test_dispersion.py # MZM-carved pulses, dispersion=on by default
│   └── common/                # README images only
├── research_roadmap.md
├── AGENTS.md
├── CHANGELOG.md
├── pyproject.toml
├── requirements.txt
└── README.md
```

## Files Changed (recent sessions — most recent first)
| File | Change |
|---|---|
| `analysis/qber_vs_distance_dispersion.py` | **New**: QBER vs distance sweep with MZM-carved 5 ps pulses; dispersion ON vs OFF |
| `src/protocols/bb84_test_dispersion.py` | **New**: MZM-carved pulses for CD/PMD testing; dispersion=True by default |
| `src/protocols/bb84_ideal.py` | Migrated to `sample_field()`, added `--dispersion` CLI flag |
| `src/protocols/bb84_high_bitrate.py` | Migrated to `sample_field()`, added `--dispersion` CLI flag |
| `src/lasers/cwlaser.py` | `get_electric_field()` → `instantaneous_field()`; updated docstrings for both methods |
| `tests/test_cwlaser.py` | Updated test names for `instantaneous_field` |
| `src/channel/fiber.py` | Updated docstring reference |
| `README.md` | Updated to CWLaser + `instantaneous_field` example |
| `CHANGELOG.md` | Added BB84 migration + dispersion test entries |
| `AGENTS.md` | Updated for all sessions |
