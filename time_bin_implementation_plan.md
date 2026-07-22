# Time-Bin Phase-Encoding QKD — Implementation Plan

## Goal
Replicate **Gobby, Yuan & Shields (2004)** — "Quantum key distribution over 122 km of standard telecom fiber" (Appl. Phys. Lett. 84, 3762) — the first QKD demonstration beyond 100 km.

## Why time-bin phase encoding?
| Aspect | Polarization encoding (current) | Time-bin phase encoding (target) |
|--------|-------------------------------|----------------------------------|
| Fiber birefringence | Scrambles signal → QBER ~50% at >100 m | Same rotation on both pulses → relative phase preserved |
| Polarization compensation | Required, not implemented | Not needed |
| Gobby et al. match | No | Yes |
| Practical relevance | Limited | Industry standard for fiber QKD |

---

## Part 1 — Pulsed Laser Operation (`src/lasers/cwlaser.py`)

### What changes
Add a switchable pulsed mode to `CWLaser.sample_field()`. When `pulsed=True`, the output is a train of Gaussian pulses instead of CW.

### Parameters (new, all with defaults so existing code is unaffected)
| Parameter | Default | Description |
|-----------|---------|-------------|
| `pulsed` | `False` | Enable pulsed operation |
| `pulse_width` | `100e-12` | FWHM of Gaussian pulse envelope (s) |
| `repetition_rate` | `2.5e6` | Pulse repetition frequency (Hz) |
| `timing_jitter_rms` | `0.0` | RMS timing jitter (s) |

### Algorithm
```
sample_field(dt, n_samples):
    1. Generate CW field as currently (phase noise, RIN, amplitude)
    2. If not pulsed: return CW field (unchanged path)
    3. If pulsed:
       a. Build pulse train envelope:
          - sigma = pulse_width / (2*sqrt(2*ln(2)))  (Gaussian sigma from FWHM)
          - T_period = 1/repetition_rate
          - n_period = int(T_period / dt)  (samples per period)
          - For each sample index i:
              t = i * dt
              nearest_pulse_center = round(t / T_period) * T_period
              if timing_jitter_rms > 0:
                  perturb center by Gaussian(0, timing_jitter_rms)
              envelope[i] = exp(-(t - t_center)^2 / (2*sigma^2))
       b. Apply envelope: E *= envelope[:, np.newaxis]
       c. Scale power so that energy per pulse = P_mean * T_period
          (i.e., average power matches power_dbm setting)
    4. Return field
```

### Backward compatibility
- `pulsed=False` by default → all existing tests and scripts unchanged
- Existing parameters (`power_dbm`, `linewidth`, `RIN`, polarization) all apply in pulsed mode
- Power calibration: `mean(|E|²)` over a full period = `_power_w`

### Files changed
- `src/lasers/cwlaser.py` — add `pulsed`, `pulse_width`, `repetition_rate`, `timing_jitter_rms` parameters; modify `sample_field()`

---

## Part 2 — Asymmetric Mach-Zehnder Interferometer (`src/channel/interferometer.py`)

### What's needed
A new component that creates and recombines time-bin pairs. Two modes:

**Encoder mode** (Alice): Single input pulse → two time-separated pulses (early + late), optionally with a phase shift applied to the late bin.

**Decoder mode** (Bob): Two time-bin pulses → recombine to interfere, outputting to constructive and destructive ports.

### Design
```python
class AsymmetricMZI:
    """Unbalanced Mach-Zehnder interferometer for time-bin encoding/decoding.

    Encoder: splits input, delays one copy by `delay`, recombines.
    Decoder: splits input, delays one copy by `delay`, applies phase via
             an internal PhaseModulator, recombines via a 50:50 coupler
             into two output ports (constructive / destructive).

    Parameters
    ----------
    delay : float — differential delay between arms (seconds).
    mode : {'encoder', 'decoder'} — operating mode.
    pm : PhaseModulator or None — phase modulator in the delayed arm.
    insertion_loss_db : float or None.
    """
    def __init__(self, delay, mode='encoder', pm=None,
                 insertion_loss_db=None):
        ...

    def modulate(self, E, dt, phase=None):
        """Apply MZI transformation.

        Parameters
        ----------
        E : ndarray (N, 2) — input field.
        dt : float — sampling interval (required for delay → samples).
        phase : float or None — phase shift applied to delayed arm.

        Returns (encoder mode)
        -------
        ndarray (N, 2) — field with two time bins.

        Returns (decoder mode)
        -------
        tuple (E_constructive, E_destructive) — two output ports.
        """
```

### Implementation details
```
modulate(E, dt, phase=None):
    delay_samples = int(self.delay / dt)
    N = len(E)

    # 50:50 split (Hadamard)
    E_short = E / sqrt(2)
    E_long = E / sqrt(2)

    # Delay long arm (shift right, zero-fill left)
    E_long = np.roll(E_long, delay_samples, axis=0)
    E_long[:delay_samples] = 0

    if mode == 'encoder':
        # Apply phase to long arm
        if phase is not None:
            E_long *= exp(1j * phase)
        # Recombine (just add — single output)
        return E_short + E_long

    else:  # decoder
        # Apply phase to delayed arm
        if phase is not None:
            E_long *= exp(1j * phase)
        # 50:50 combiner → two outputs
        E_constructive = (E_short + E_long) / sqrt(2)
        E_destructive = (E_short - E_long) / sqrt(2)
        return E_constructive, E_destructive
```

### Backward compatibility
- New file, no existing code changes

### Files changed
- `src/channel/interferometer.py` (NEW)
- `src/channel/__init__.py` — export `AsymmetricMZI`

---

## Part 3 — SPAD Time-Gated Detection (`src/detectors/spad.py`)

### What changes
The current SPAD `detect(power, t)` method already supports gated detection. For time-bin QKD, we use it the same way — but we need to call `detect()` for each detector at the interference time (t_pulse + Δt).

**No changes needed to the SPAD class itself.** The protocol code will handle:
1. Two SPAD instances (one per interferometer output port)
2. One detection per pulse (at the time-bin overlap time)
3. The SPAD's dead time, DCR, and afterpulsing already work correctly for this

### What we might want
A convenience method `detect_time_bin(powers_constructive, powers_destructive, times)` that:
- Takes arrays of constructive/destructive port powers
- Calls `self.detect(power_c, t)` and `self.detect(power_d, t)` for each pulse
- Returns click pairs

**Low priority** — can be added to the protocol directly.

### Files changed
- `src/detectors/spad.py` — optional convenience method (or skip and handle in protocol)

---

## Part 4 — Phase-Encoding BB84 Protocol (`src/protocols/bb84_time_bin.py`)

### New protocol script

```
signal chain:
  CWLaser(pulsed=True)
    → AsymmetricMZI(mode='encoder')    # create time bins
    → PhaseModulator (encodes φ_A on late bin)
    → propagate()                       # fiber channel
    → AsymmetricMZI(mode='decoder')    # recombine
    → 2× SPAD                          # detect constructive/destructive
```

### Encoding scheme (Gobby et al.)
| Basis | Bit 0 | Bit 1 | Bob phase |
|-------|-------|-------|-----------|
| X     | φ=0   | φ=π   | φ_B = 0   |
| Y     | φ=π/2 | φ=3π/2| φ_B = π/2 |

### Detection logic
```
I_constructive ∝ (1 + cos(φ_A - φ_B)) / 2
I_destructive  ∝ (1 - cos(φ_A - φ_B)) / 2

Same basis (Δφ = 0 or π): one detector fires → deterministic bit
Diff basis (Δφ = ±π/2): both fire equally → random (sifted out)
```

### Key parameters (Gobby et al.)
| Parameter | Value |
|-----------|-------|
| Wavelength | 1550 nm |
| Pulse width | 100 ps (Gaussian FWHM) |
| Repetition rate | 2.5 MHz |
| Mean photons/pulse (μ) | ~0.1 |
| AMZI delay | 5.8 ns (~1.2 m fiber) |
| SPAD detection efficiency | 10% |
| SPAD dark count prob/gate | 8.5×10⁻⁷ |
| Fiber attenuation | 0.182 dB/km (SMF-28) |

### Files changed
- `src/protocols/bb84_time_bin.py` (NEW)

---

## Part 5 — Gobby et al. Replication (`analysis/val_gobby/`)

### Validation script
```
analysis/val_gobby/validate_gobby.py:
  1. Set up pulsed laser with Gobby parameters
  2. Sweep distance: 0, 10, 20, 40, 65, 101, 122 km
  3. For each distance: run BB84 protocol, measure QBER
  4. Plot QBER vs distance overlaid on Gobby Fig 3 data
  5. Plot sifted key rate vs distance overlaid on Gobby Fig 4
  6. Output: val_gobby--seed{N}.png, val_gobby_table.tex
```

### Expected results (Gobby paper)
| Distance (km) | QBER (%) |
|---------------|----------|
| 4.4 | ~3.3 |
| 65.0 | ~3.3 |
| 101.0 | ~6.0 |
| 122.0 | 8.9 |

### Analytical check
QBER formula from Gobby Eq 1-2:
```
QBER = (P_dark * ν_gate) / (μ * η * 10^(-αL/10)) + QBER_opt

where:
  P_dark = 8.5e-7 (dark count prob per gate)
  ν_gate = 2.5 MHz (gating rate)
  μ = 0.1 (mean photons/pulse)
  η = 0.10 (detection efficiency)
  α = 0.2 dB/km (fiber attenuation)
  QBER_opt ≈ 3.3% (baseline optical error)
```

### Files changed
- `analysis/val_gobby/validate_gobby.py` (NEW)
- `analysis/val_gobby/` output directory

---

## Part 6 — Testing Plan

### New tests

| File | Test | What it checks |
|------|------|----------------|
| `tests/test_cwlaser.py` | `test_pulsed_mode_shape` | Pulsed `sample_field()` returns correct shape (N,2) |
| `tests/test_cwlaser.py` | `test_pulsed_power_conservation` | Mean power in pulsed mode = `_power_w` |
| `tests/test_cwlaser.py` | `test_pulsed_pulse_width` | Gaussian FWHM matches `pulse_width` |
| `tests/test_cwlaser.py` | `test_pulsed_repetition_rate` | Pulse spacing matches `repetition_rate` |
| `tests/test_cwlaser.py` | `test_pulsed_cw_identical` | `pulsed=False` gives same output as before |
| `tests/test_interferometer.py` | `test_encoder_creates_two_bins` | Encoder output has two distinct time bins |
| `tests/test_interferometer.py` | `test_decoder_interference` | Decoder output interference depends correctly on phase |
| `tests/test_interferometer.py` | `test_encoder_decoder_roundtrip` | Encoder + decoder with same delay reconstructs original |
| `tests/test_interferometer.py` | `test_phase_response` | Output varies as cos²(Δφ/2) |
| `tests/test_interferometer.py` | `test_power_conservation` | Total power conserved through AMZI |
| `tests/test_spad.py` | (existing) | No changes needed |

### Existing tests that must continue to pass
- All 48 existing tests
- Pulsed mode must not affect any `pulsed=False` behavior

### Validation criteria (what constitutes "success")
1. **Pulsed laser**: Pulse envelope FWHM within 5% of `pulse_width` parameter
2. **AMZI encoder**: Two time bins separated by `delay` ± 1 sample
3. **AMZI decoder**: Constructive/destructive port power proportional to `(1 ± cos(Δφ))` with R² > 0.99
4. **Gobby replication**: QBER vs distance curve within 20% relative error of paper data points
5. **Zero-fiber case**: QBER ≈ QBER_opt (dominated by detector dark counts at μ = 0.1)

---

## Implementation Order

| Step | What | Depends on | Effort |
|------|------|-----------|--------|
| 1 | Pulsed laser mode | Nothing | ~1 session |
| 2 | AsymmetricMZI component | Nothing | ~1 session |
| 3 | Tests for 1 + 2 | 1, 2 | ~1 session |
| 4 | Time-bin BB84 protocol | 1, 2 | ~1 session |
| 5 | Gobby replication script | 4 | ~1 session |
| 6 | Validate vs paper + iterate | 5 | ~1 session |
| **Total** | | | **~5 sessions** |

## Risk Assessment

| Risk | Impact | Mitigation |
|------|--------|------------|
| Timing jitter model inaccurate | Low — Gobby uses gain-switched DFB with ~few ps jitter, negligible at 100 ps pulse width | Start with `timing_jitter_rms=0` |
| AMZI delay resolution limited by dt | Medium — if dt is too coarse, delay_samples may be off by 1 | Use `dt` matching the laser's internal resolution |
| QBER at 122 km may be dominated by different physics (afterpulsing) than paper | Medium — paper uses id210 APD, we model SPAD differently | Check afterpulsing contribution; add afterpulsing suppression if needed |
| Fiber CD may spread 100 ps pulse significantly at 122 km | Low — D = 17 ps/(nm·km), with typical DFB linewidth ~1 MHz, CD broadening is negligible | Validate with CD on/off |
| Power calibration in pulsed mode | Medium — need `mean(|E|²)` over full period = `_power_w` | Verify with pulse train power measurement test |
