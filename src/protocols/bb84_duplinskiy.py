"""Replication of Duplinskiy et al., Opt. Express 25(23), 28886, 2017.

"Low loss QKD optical scheme for fast polarization encoding"

Proof-of-concept BB84 QKD over 50 km fiber spool (10 dB loss)
using LiNbO3 phase modulators and InGaAs SPADs.

Signal chain (mirrors bb84_ideal.py but with VOA + SPAD):
  Alice: polarizer('45') → PM1 (encode)
  Channel: propagate (birefringence + attenuation)
  Bob: PM2 (basis select) → circular_analyser → 2x SPAD

Detection uses the same 50:50 BS with π/2 phase shift
(optics.circular_analyser — not a true PBS, see PHYS-6 in
opto-sim-issues-and-fixes.md):
  P_x − P_y = sin(Δφ_alice + Δφ_bob)
  Same basis: |sin| = 1 → deterministic bit
  Diff basis: sin ≈ 0 → random (sifted out)

Encoding (same as bb84_ideal.py):
  X basis (linear):  bit 0 → V=0,       bit 1 → V=Vpi
  C basis (circular): bit 0 → V=Vpi/2,  bit 1 → V=3·Vpi/2
  Bob X basis: V_bob = Vpi/2
  Bob C basis: V_bob = 0

Validation note
---------------
The 0 km (back-to-back) QBER validates the SPAD + phase modulator +
circular-analyser detection chain against the paper's baseline
measurement.  At 0 km there
is no fibre, so birefringence does not enter — only detector dark counts,
afterpulsing, and the intrinsic sifting loss contribute.

Distance-dependent QBER values are comparable to the paper *with
polarization compensation enabled* (the default): the channel is
quasi-static, so Bob applies the inverse of the fibre's single SU(2)
Jones matrix (obtained from `FiberRealization.birefringence_matrix()`).
This mirrors the paper's calibration procedure, which maintains
alignment with polarization controllers and automatic recalibration —
for a quasi-static channel, that loop converges to exactly this
inverse.  Pass ``--no-compensation`` to see the uncompensated channel:
because the fibre applies a fixed random rotation, the received SOP no
longer matches the encoding axes and QBER collapses toward ~50 %
(25.9 % at 50 km in our seed-42 run), which is what a polarization-
encoded system *without* any active stabilization sees.

References
----------
[1] Duplinskiy et al., Opt. Express 25(23), 28886-28897, 2017.
[2] ID Quantique, ID230 InGaAs SPAD datasheet.
"""
import argparse
import numpy as np
import random
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from src.channel import optics, FiberRealization
from src.channel.phase_modulator import PhaseModulator
from src.detectors.spad import spad


def simulate_bb84_duplinskiy(num_bits, fiber_length=50, alpha_dB=0.2,
                              mu=0.1, bob_loss_dB=2.0,
                              gate_width=20e-9, rep_rate=10e6,
                              compensate=True, model='auto',
                              seed=None, verbose=False):
    """BB84 simulation matching the Duplinskiy et al. experimental setup.

    Parameters
    ----------
    num_bits : int — number of BB84 pulses to simulate.
    fiber_length : float — fiber length in km (default 50).
    alpha_dB : float — fiber attenuation in dB/km (default 0.2, SMF-28).
    mu : float — mean photon number per pulse (default 0.1).
    bob_loss_dB : float — Bob's internal loss in dB (default 2.0).
    gate_width : float — SPAD gate window in seconds (default 20 ns).
    rep_rate : float — laser repetition rate in Hz (default 10 MHz).
    compensate : bool — apply Bob's polarization compensation (the
        inverse of the fibre's quasi-static Jones matrix) before
        decoding, mirroring the paper's calibration loop (default True).
    seed : int or None — RNG seed.
    model : str — birefringence model for the fibre: 'auto' (default) or
        'sectional'. 'phenomenological' was removed in the fifth pass
        (PHYS-5) and raises ValueError.
    verbose : bool — print progress.

    Returns
    -------
    dict with keys: qber, n_total, n_sifted, n_errors, sifted_key_rate, etc.
    """
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)

    h = 6.626e-34
    c = 3e8
    wavelength = 1550e-9
    photon_energy = h * c / wavelength

    # Calibrate field power so that mean(|E|^2) * gate_width gives mu photons
    # power_per_pulse = mu * h * nu / gate_width
    power_per_pulse = mu * photon_energy / gate_width

    # Phase modulators (LiNbO3, X-cut)
    pm_alice = PhaseModulator(crystal_cut='X', modulation='DC')
    pm_bob = PhaseModulator(crystal_cut='X', modulation='DC')
    Vpi = pm_alice.Vpi

    # SPAD detectors (ID230 specs: 10% QE, 13 us dead time, 15 Hz DCR, 5% afterpulse)
    spd1 = spad(wavelength=wavelength, quantum_efficiency=0.10,
                dead_time=13e-6, dark_count_rate=15.0,
                afterpulse_prob=0.05, gate_width=gate_width)
    spd2 = spad(wavelength=wavelength, quantum_efficiency=0.10,
                dead_time=13e-6, dark_count_rate=15.0,
                afterpulse_prob=0.05, gate_width=gate_width)

    dt_pulse = 1.0 / rep_rate

    # One physical fibre for the whole run: birefringence is quasi-static,
    # so every pulse must see the same Jones matrix (see ROOT-1 in
    # opto-sim-issues-and-fixes.md). Built once here, reused per pulse
    # below. CD/PMD stay off, matching this script's original explicit
    # propagate(cd=False, pmd=False) call.
    fibre = FiberRealization(L_m=fiber_length * 1000, temperature=25,
                             bend_radius=None, attenuation_factor=alpha_dB,
                             cd=False, pmd=False, model=model, seed=seed)

    # Bob's polarization compensation (default on): the channel's Jones
    # matrix is unitary, so the inverse is its conjugate transpose.
    # Applied per pulse before decoding, mirroring the paper's
    # calibration loop (see module docstring).
    J_channel = fibre.birefringence_matrix()
    U_comp = None if J_channel is None else J_channel.conj().T

    alice_bits, alice_bases = [], []
    bob_bits, bob_bases = [], []
    has_click = []

    t_start = time.time()

    for pulse_idx in range(num_bits):
        # --- Alice's encoding ---
        alice_basis = random.choice(['C', 'X'])
        alice_bit = random.randint(0, 1)

        if alice_basis == 'C':
            if alice_bit == 0:
                phase_alice = Vpi / 2
            else:
                phase_alice = 3 * Vpi / 2
        else:  # 'X'
            if alice_bit == 0:
                phase_alice = 0
            else:
                phase_alice = Vpi

        # --- Field: 45° polarized pulse with correct power ---
        E = np.sqrt(power_per_pulse / 2.0) * np.ones((1, 2), dtype=complex)
        E = pm_alice.modulate(E_field=E, V=phase_alice)

        # --- Channel: birefringence + attenuation ---
        E = fibre.apply(E, dt=dt_pulse)

        # --- Bob's polarization compensation (calibration inverse) ---
        if compensate and U_comp is not None:
            E = np.transpose(U_comp @ np.transpose(E))

        # --- Bob's internal loss (optics, connectors, coupler) ---
        E = optics.voa(E, bob_loss_dB)

        # --- Bob's decoding ---
        bob_basis = random.choice(['C', 'X'])
        if bob_basis == 'C':
            phase_bob = 0
        else:
            phase_bob = Vpi / 2

        E = pm_bob.modulate(E_field=E, V=phase_bob)
        # circular_analyser, not pbs: detection depends on the relative
        # phase between Ex/Ey (PHYS-6 in opto-sim-issues-and-fixes.md).
        Ex, Ey = optics.circular_analyser(E)

        # --- Detection ---
        power_x = np.mean(np.abs(Ex) ** 2)
        power_y = np.mean(np.abs(Ey) ** 2)

        t_pulse = pulse_idx * dt_pulse
        click_x = spd1.detect(power_x, t_pulse)
        click_y = spd2.detect(power_y, t_pulse)

        if click_x and not click_y:
            bob_bit = 0
        elif click_y and not click_x:
            bob_bit = 1
        elif click_x and click_y:
            bob_bit = random.randint(0, 1)
        else:
            bob_bit = -1  # no click — discard this pulse

        alice_bits.append(alice_bit)
        alice_bases.append(alice_basis)
        bob_bits.append(bob_bit)
        bob_bases.append(bob_basis)
        has_click.append(click_x or click_y)

        if verbose and (pulse_idx + 1) % 10000 == 0:
            elapsed = time.time() - t_start
            rate = (pulse_idx + 1) / elapsed
            print(f"  Pulse {pulse_idx+1}/{num_bits} ({rate:.0f} pulses/s)", flush=True)

    # --- Sifting: same basis AND at least one click ---
    sifted_indices = [i for i in range(num_bits)
                      if alice_bases[i] == bob_bases[i] and has_click[i]]
    sifted_alice = [alice_bits[i] for i in sifted_indices]
    sifted_bob = [bob_bits[i] for i in sifted_indices]

    n_sifted = len(sifted_alice)
    n_errors = sum(a != b for a, b in zip(sifted_alice, sifted_bob))
    qber = n_errors / n_sifted if n_sifted > 0 else 0.0

    fiber_loss_dB = alpha_dB * fiber_length
    total_loss_dB = fiber_loss_dB + bob_loss_dB
    sifted_key_rate = n_sifted / (num_bits * dt_pulse) if num_bits > 0 else 0.0
    elapsed = time.time() - t_start

    results = {
        'qber': qber,
        'n_total': num_bits,
        'n_sifted': n_sifted,
        'n_errors': n_errors,
        'sifted_key_rate': sifted_key_rate,
        'fiber_length_km': fiber_length,
        'total_loss_dB': total_loss_dB,
        'mu': mu,
        'elapsed_s': elapsed,
    }

    if verbose:
        n_clicks = sum(has_click)
        print(f"\nDuplinskiy et al. BB84 — {fiber_length} km")
        print(f"  Polarization compensation: {'ON' if U_comp is not None else 'OFF'}")
        print(f"  Total loss: {total_loss_dB:.1f} dB (fiber {fiber_loss_dB:.1f} + Bob {bob_loss_dB:.1f})")
        print(f"  Mu: {mu} photons/pulse")
        print(f"  Total pulses: {num_bits}")
        print(f"  Clicks: {n_clicks} ({n_clicks/num_bits*100:.2f}%)")
        print(f"  Sifted (same basis + click): {n_sifted} ({n_sifted/num_bits*100:.2f}%)")
        print(f"  Errors: {n_errors}")
        print(f"  QBER: {qber*100:.2f}%")
        print(f"  Time: {elapsed:.1f}s")

    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Replicate Duplinskiy et al. BB84 experiment")
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--bits', type=int, default=100000,
                        help='Number of pulses (default 100k)')
    parser.add_argument('--fiber-length', type=float, default=50,
                        help='Fiber length in km (default 50)')
    parser.add_argument('--mu', type=float, default=0.1,
                        help='Mean photons per pulse (default 0.1)')
    parser.add_argument('--no-compensation', action='store_true',
                        help='Disable Bob polarization compensation (the '
                             'fibre\'s fixed random SU(2) then scrambles '
                             'the encoding)')
    parser.add_argument('--birefringence-model', choices=['auto', 'sectional'],
                        default='auto',
                        help='Birefringence model (default: auto)')
    parser.add_argument('--verbose', action='store_true')
    args = parser.parse_args()

    results = simulate_bb84_duplinskiy(
        num_bits=args.bits, fiber_length=args.fiber_length,
        mu=args.mu, seed=args.seed, verbose=True,
        compensate=not args.no_compensation,
        model=args.birefringence_model)
