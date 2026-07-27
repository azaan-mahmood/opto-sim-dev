# opto-sim: Physical-Layer Fiber-Optic Simulation Framework

**opto-sim** is an open-source Python framework for physical-layer simulation of optical fibre channels, with a focus on quantum key distribution (QKD) systems. The complex-envelope electric field $\mathbf{E}(t)=[E_x(t),E_y(t)]^\mathsf{T}$ is the single source of truth — optical power is derived directly from the field as $\langle|E|^2\rangle$, calibrated once at the laser output and tracked through every downstream component.

Seven components are independently validated against published literature (38 references), then composed into a system-level time-bin BB84 QKD demonstration that reproduces Gobby et al. (2004) with 0 km QBER = 2.43% versus the paper's 3.3%.

## Key Features

- **Complex-envelope field propagation** — full polarisation state (Stokes parameters), phase noise, RIN, and power tracked through every component
- **7 validated components**, each against its own literature source:
  - **CW Laser** with switchable pulsed mode, Wiener phase noise (Henry 1982), relaxation-oscillation RIN (Coldren 2012)
  - **Mach–Zehnder Modulator** with push-pull and single-drive configurations (Agrawal, Koyama & Iga)
  - **Asymmetric Mach–Zehnder Interferometer** for time-bin phase encoding/decoding
  - **Fibre Channel** with birefringence (sectional + phenomenological models), chromatic dispersion (FFT-based, Agrawal), PMD (Maxwellian DGD, Razavi), and attenuation (Keiser)
  - **Avalanche Photodiode** with shot noise, thermal noise, excess noise factor (Kasap)
  - **Geiger-mode SPAD** with dead time, dark count rate, afterpulsing, gated detection (ID230 specs)
  - **Time-bin BB84 Protocol** with pulsed laser, AMZI encoder/decoder, dual SPAD detection
- **Seeded reproducibility** — all RNG seeded at session start (`--seed` default 42)
- **77 unit tests** — component-level validation of power convention, noise scaling, polarisation, interference fringes, edge cases
- **No GUI, no proprietary dependencies** — pure Python + NumPy + SciPy + Matplotlib

## Architecture

```
src/
├── lasers/cwlaser.py          # CW laser with pulsed mode option
├── channel/
│   ├── fiber.py               # birefringence, CD, PMD, attenuation
│   ├── mzm.py                 # Mach-Zehnder modulator
│   ├── interferometer.py      # asymmetric MZI for time-bin encoding
│   ├── optics.py              # VOA, beam splitters, etc.
│   └── phase_modulator.py     # LiNbO3 phase modulator
├── detectors/
│   ├── apd.py                 # linear-mode APD
│   └── spad.py                # Geiger-mode SPAD
├── protocols/
│   ├── bb84_ideal.py          # CW-based BB84
│   ├── bb84_high_bitrate.py   # bitrate-sweep BB84
│   ├── bb84_test_dispersion.py# MZM-carved pulses with CD/PMD
│   ├── bb84_duplinskiy.py     # SPAD-based BB84 replication
│   └── bb84_time_bin.py       # time-bin phase-encoding BB84
└── visualization/
    ├── stokes.py               # Stokes parameter computation
    ├── polarimeter.py          # Poincaré sphere plotting
    └── fields.py               # field visualization
```

The field convention is consistent: `mean(|E|²)` = optical power in Watts. Every component reads and writes the same field array `E` with shape `(n_samples, 2)`.

## Validation

Each component is validated against a published source — see `analysis/validation/` for the full scripts and `paperwork/main.tex` for the manuscript.

| Component | Method | Error |
|---|---|---|
| CW Laser | Power convention, phase noise (Henry), RIN (Coldren) | < 1 % |
| MZM | Transfer function, null/peak, quadrature (Agrawal) | < 0.1 % |
| AMZI | Fringe visibility cos²(Δφ/2) | < 0.1 % |
| CD | Gaussian pulse broadening (Agrawal Fig 2.6) | < 0.06 % |
| PMD | Maxwellian DGD histogram (Razavi Fig 2.11) | KS test pass |
| Attenuation | SMF-28 distance sweep (Keiser Eq 3.6) | < 1e-12 |
| Birefringence | Beat length vs bend radius (Ulrich 1980) | < 0.1 % |
| APD | Responsivity, noise floor (Kasap Eq 4.45) | < 1 % |
| SPAD | Gobby replication (Appl. Phys. Lett. 84, 2004) | 0 km QBER 2.43 % |

## Installation

```bash
pip install -r requirements.txt
python -m pytest tests/ -v      # run all tests
```

## Quick Start

```python
from src.lasers import CWLaser
from src.channel import fiber
from src.detectors import APD
import numpy as np

laser = CWLaser(wavelength=1550e-9, power_dbm=0)
E = laser.sample_field(dt=1e-12, n_samples=1000)

f = fiber.Fiber(fiber_length=50)
E = f.propagate(E, dt=1e-12, wavelength=1550e-9)

apd = APD()
result = apd.output(E, bandwidth=1e9)
print(f"Signal current: {result['I_signal']:.2e} A")
```

## Running BB84

```bash
python -m src.protocols.bb84_ideal --fiber-length 50
python -m src.protocols.bb84_time_bin --distance 10
```

All scripts accept `--seed` for reproducibility.

## Reproducibility

Every test and validation script pins both `random` and `np.random` at session start via `conftest.py`. The `--seed` CLI flag (default 42) controls all RNG. Validation outputs are tagged with the seed in the filename (e.g., `val_gobby--seed42.png`).

## Citation

If you use this framework in your research, please cite the accompanying manuscript:

```
A. Mahmood, "An Open-Source Physically Validated Optical Channel Simulation
Framework for Quantum Key Distribution," arXiv preprint, 2026.
```

## License

MIT License — see LICENSE file.

## Contact

A. Mahmood — azaan.mahmood@dsu.edu.pk
