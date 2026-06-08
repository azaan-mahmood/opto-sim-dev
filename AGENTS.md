# Project Context: opto-sim QKD Simulation

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
- `get_electric_field(dt=1e-12, over_period=False, n_samples=1000, normalize=True)` — backward-compatible, `t`→`dt`, hardcoded 1000→`n_samples` param.
- Now used in BB84 example scripts (replaced SolidStateLaser)
- BB84 scripts call `get_electric_field(normalize=False, over_period=True)` — power already calibrated, no extra rescaling needed
- `dispersion=False` in cable(): the FFT dispersion model is blocked (unphysical frequencies), disabled to avoid field corruption

### APD (`src/detectors/apd.py`)
- Current-based output: `output(E, bandwidth)` derives `power = mean(|E|²)` from the field
- Kasap Eq 4.45: excess noise factor `F` applied **only** to shot noise, not thermal
- `detect_photons()` uses `power/(h·ν) · t · η` (Agrawal Eq 4.1.2, Saleh & Teich Eq 17.1-10)
- No `frequency` parameter — optical frequency derived from `c/λ`

### Fiber (`src/opto_eq/fiber.py`)
- `cable()` no longer takes `pin` or returns `pout`
- Attenuation applied directly to field: `E *= sqrt(10^(-α·L/10))` (Keiser Eq 3.6)
- Birefringence Jones matrix `[[exp(j·δβ·L/2), 0], [0, 1]]` shifts Ex phase

### Mach-Zehnder Modulator (`src/opto_eq/mzm.py`)
- Push-pull configuration: `E_out = E_in · cos(π·V/V_pi) · exp(j·π·V_bias/V_pi)`
- `switching_voltage = V_pi / 2` (voltage for ON→OFF transition)
- Parameters: `V_pi` (single-arm π voltage, default 5 V), `bias_voltage` (common-mode phase, default 0)
- Agrawal [1] §4.2: External Modulation and Mach-Zehnder Modulators

### Stokes viewer (`src/viewers/stokes.py`)
- `compute_stokes_parameters(E)` returns `[S0,S1,S2,S3], [psi, chi]`
- S1,S2,S3 normalized by S0; S0 set to 1.0
- `chi = 0.5·arcsin(S3)` (S3 already normalized, S0=1)

### Chromatic Dispersion (`src/opto_eq/fiber.py:apply_dispersion()`)
- **BLOCKED**: `f = (1/100)*(1/(D*L*0.2e-9))` produces unphysical frequencies ~10¹⁷ Hz
- Needs rewrite with proper FFT-based model (deferred to later session)

## Known Issues / Blockers

### SolidStateLaser (`src/lasers/sslaser.py`)
- ODE rate equations never produce population inversion because `Rp = 1/τ₂` (transparency threshold). Needs `Rp > 1/τ₂` by at least 2-10×.
- `power_dbm` only sets initial photon density `I_0`, doesn't affect pump — changing it changes none of the dynamics.
- Electric field units are mW-based (not Watts). BB84 scripts calibrate with `E *= sqrt(power_W / mean(|E|²))` at laser output.
- Awaiting replacement with CWLaser once dispersion is fixed.

### BB84 Scripts
- `bb84_ideal.py` / `bb84_high_bitrate.py`: use CWLaser (no power calibration needed) → polarizer → phase modulator → cable (dispersion=False) → PBS → APD
- QBER is 0% at 1-200 km (correct: differential detection, signal dominates noise)
- QBER is 0% at all tested bandwidths (1 MHz–5 GHz) at 100 km: CWLaser's -5 dBm with M=10 provides ~3.5 mA signal current, far exceeding the noise floor

## Generated Graphs
- `analysis/qber_vs_distance.png`: 0% QBER 10-190 km
- `analysis/qber_vs_bitrate.png`: 0% at 215 MHz → 35% at 10 GHz
- `analysis/laser_characterization.png`: 8-panel CWLaser dashboard
- `analysis/poincare_sphere.png`: Poincaré sphere from Stokes parameters
- `analysis/eye_diagrams.png`: NRZ-OOK eye diagrams at 5/10/25 Gbaud

## Files Changed (recent sessions)
| File | Change |
|---|---|
| `src/lasers/cwlaser.py` | Relaxation-oscillation RIN model (Coldren Eq 5.3.38), rewritten `_sample_rin()` |
| `src/opto_eq/fiber.py` | Removed `pin`/`pout`, attenuation on field only |
| `src/detectors/apd.py` | `output(E, ...)` derives power from field |
| `src/protocols/examples/bb84_ideal.py` | Switched to CWLaser, removed power calibration, dispersion=False |
| `src/protocols/examples/bb84_high_bitrate.py` | Switched to CWLaser, removed power calibration |
| `main.py` | Updated cable/stokes calls |
| `src/viewers/stokes.py` | Normalized S0=1.0, returns `[S0,S1,S2,S3], [psi,chi]` |
| `CHANGELOG.md` | Full audit trail with line ranges and timestamps |
| `analysis/laser_characterization.py` | Rewritten: uses `laser.sample_field()`, MZM model, no inline noise generation |
| `src/lasers/cwlaser.py` | Added `sample_field()`, renamed `t`→`dt`, `n_samples` param |
| `src/opto_eq/mzm.py` | New: physical MZM device model (push-pull) |
| `src/opto_eq/__init__.py` | Exports MZM |
