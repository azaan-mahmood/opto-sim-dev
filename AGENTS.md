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
- Not yet used in BB84 scripts (still uses SolidStateLaser)

### APD (`src/detectors/apd.py`)
- Current-based output: `output(E, bandwidth)` derives `power = mean(|E|²)` from the field
- Kasap Eq 4.45: excess noise factor `F` applied **only** to shot noise, not thermal
- `detect_photons()` uses `power/(h·ν) · t · η` (Agrawal Eq 4.1.2, Saleh & Teich Eq 17.1-10)
- No `frequency` parameter — optical frequency derived from `c/λ`

### Fiber (`src/opto_eq/fiber.py`)
- `cable()` no longer takes `pin` or returns `pout`
- Attenuation applied directly to field: `E *= sqrt(10^(-α·L/10))` (Keiser Eq 3.6)
- Birefringence Jones matrix `[[exp(j·δβ·L/2), 0], [0, 1]]` shifts Ex phase

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
- `bb84_ideal.py` / `bb84_high_bitrate.py`: use SolidStateLaser with field calibration → cable → PBS → APD
- QBER is 0% at 1-150 km with dispersion off (correct: differential detection, signal dominates noise)
- QBER rises from 0% at ~200 MHz to 35% at 10 GHz (100 km): noise scales with √bandwidth, 11% abort threshold at ~4.6 GHz

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
| `src/protocols/examples/bb84_ideal.py` | Field calibration to Watts, new cable/apd signatures |
| `src/protocols/examples/bb84_high_bitrate.py` | Same |
| `main.py` | Updated cable/stokes calls |
| `src/viewers/stokes.py` | Normalized S0=1.0, returns `[S0,S1,S2,S3], [psi,chi]` |
| `CHANGELOG.md` | Full audit trail with line ranges and timestamps |
| `analysis/laser_characterization.py` | CWLaser verification suite: power, linewidth, phase noise, RIN spectrum, Stokes, eye diagrams |
