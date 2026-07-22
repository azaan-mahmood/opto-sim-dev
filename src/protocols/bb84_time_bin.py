"""Time-bin phase-encoding BB84 QKD protocol.

Signal chain:
  Pulsed laser → AsymmetricMZI(encoder) → PhaseModulator (φ_A)
    → propagate() (fiber) → AsymmetricMZI(decoder) → 2× SPAD

Based on Gobby, Yuan & Shields (2004), Appl. Phys. Lett. 84, 3762.

References
----------
[1] Gobby, C., Yuan, Z. L., & Shields, A. J. (2004). Quantum key
    distribution over 122 km of standard telecom fiber. Appl. Phys.
    Lett. 84(19), 3762-3764.
[2] Bennett, C. H. & Brassard, G. (1984). Quantum cryptography:
    Public key distribution and coin tossing. Proc. IEEE Int. Conf.
    on Computers, Systems and Signal Processing, 175-179.
[3] Townsend, P. D. et al. (1993). Single photon interference in
    10 km long optical fibre interferometer. Electron. Lett. 29(7),
    634-635.
"""
import argparse
import numpy as np
import random
import sys

sys.path.insert(0, 'src')
from src.channel.interferometer import AsymmetricMZI
from src.detectors.spad import spad


def gaussian_pulse(t, sigma, A=1.0):
    """Gaussian pulse envelope, unit peak amplitude."""
    return A * np.exp(-t ** 2 / (2 * sigma ** 2))


def simulate_bb84_time_bin(num_bits, fiber_length=0, alpha_dB=0.182,
                            mu=0.1, wavelength=1550e-9,
                            pulse_width=100e-12, repetition_rate=2.5e6,
                            delay=5.8e-9, gate_width=1e-9,
                            spad_eta=0.10, dark_count_rate=15.0,
                            afterpulse_prob=0.05, dead_time=13e-6,
                            seed=None, verbose=False):
    """BB84 time-bin phase-encoding simulation.

    Parameters
    ----------
    num_bits : int — number of BB84 pulses to simulate.
    fiber_length : float — fiber length in km (default 0).
    alpha_dB : float — fiber attenuation in dB/km (default 0.182).
    mu : float — mean photon number per pulse (default 0.1).
    wavelength : float — centre wavelength (m, default 1550e-9).
    pulse_width : float — FWHM of Gaussian pulse (s, default 100e-12).
    repetition_rate : float — pulse repetition rate (Hz, default 2.5e6).
    delay : float — AMZI differential delay (s, default 5.8e-9).
    gate_width : float — SPAD gate window (s, default 1e-9).
    spad_eta : float — SPAD quantum efficiency (default 0.10).
    dark_count_rate : float — SPAD dark count rate (Hz, default 15.0).
    afterpulse_prob : float — afterpulse probability (default 0.05).
    dead_time : float — SPAD dead time (s, default 13e-6).
    seed : int or None — RNG seed.
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
    photon_energy = h * c / wavelength

    # Gaussian sigma from FWHM
    sigma = pulse_width / (2.0 * np.sqrt(2.0 * np.log(2.0)))

    # Sampling: must resolve pulse_width (dt <= pulse_width / 10)
    dt = min(pulse_width / 10.0, delay / 20.0)

    # Samples needed: cover original pulse + encoder + decoder + margin
    n_samples = int(np.ceil((2.0 * delay + 5.0 * pulse_width) / dt))
    t = np.arange(n_samples, dtype=float) * dt

    # Build a single Gaussian pulse centred at t = delay/2 so that the
    # interference (at t = delay/2 + delay) falls within the window.
    pulse_center = delay / 2.0
    energy_per_pulse = mu * photon_energy
    # Peak field amplitude: energy = A² ∫exp(-t²/σ²) dt = A² · σ·√π (continuous)
    # In discrete: Σ|A·g|²·dt = A²·σ·√π  →  A² = energy / (σ·√π)
    pulse_amplitude = np.sqrt(energy_per_pulse / (sigma * np.sqrt(np.pi)))
    pulse = gaussian_pulse(t - pulse_center, sigma, A=pulse_amplitude)
    E_pulse = pulse[:, np.newaxis] * np.array([1.0, 0.0], dtype=complex)  # X-polarised

    # Build AMZI components
    enc = AsymmetricMZI(delay=delay, mode='encoder')
    dec = AsymmetricMZI(delay=delay, mode='decoder')

    # SPAD detectors (constructive and destructive ports)
    spad_c = spad(wavelength=wavelength, quantum_efficiency=spad_eta,
                  dead_time=dead_time, dark_count_rate=dark_count_rate,
                  afterpulse_prob=afterpulse_prob, gate_width=gate_width)
    spad_d = spad(wavelength=wavelength, quantum_efficiency=spad_eta,
                  dead_time=dead_time, dark_count_rate=dark_count_rate,
                  afterpulse_prob=afterpulse_prob, gate_width=gate_width)

    # Fiber loss factor
    fiber_loss_lin = 10.0 ** (-alpha_dB * fiber_length / 20.0)  # field factor (sqrt of power)

    # SPAD gate in samples (half-window for power extraction)
    gate_half_samples = max(1, int(gate_width / dt / 2.0))
    delay_samples = int(delay / dt)
    interference_idx = pulse_center_idx = int(pulse_center / dt) + delay_samples

    alice_bits, alice_bases = [], []
    bob_bits, bob_bases = [], []
    has_click = []

    for pulse_idx in range(num_bits):
        # --- Alice's encoding ---
        alice_basis = random.choice(['X', 'Y'])
        alice_bit = random.randint(0, 1)

        if alice_basis == 'X':
            phi_A = 0.0 if alice_bit == 0 else np.pi
        else:  # Y basis
            phi_A = np.pi / 2.0 if alice_bit == 0 else 3.0 * np.pi / 2.0

        # --- Encode time bins with phase φ_A ---
        E_field = enc.modulate(E_pulse, dt, phase=phi_A)

        # --- Fiber channel (attenuation only; time-bin immune to birefringence) ---
        E_field *= fiber_loss_lin

        # --- Bob's decoding ---
        bob_basis = random.choice(['X', 'Y'])
        phi_B = 0.0 if bob_basis == 'X' else np.pi / 2.0

        E_c, E_d = dec.modulate(E_field, dt, phase=phi_B)

        # --- Extract interference power ---
        start_i = max(0, interference_idx - gate_half_samples)
        end_i = min(len(E_c), interference_idx + gate_half_samples + 1)
        P_c = float(np.mean(np.sum(np.abs(E_c[start_i:end_i]) ** 2, axis=1)))
        P_d = float(np.mean(np.sum(np.abs(E_d[start_i:end_i]) ** 2, axis=1)))

        # --- Detection ---
        click_c = spad_c.detect(P_c, float(pulse_idx) / repetition_rate)
        click_d = spad_d.detect(P_d, float(pulse_idx) / repetition_rate)

        # --- Determine Bob's bit ---
        if click_c and not click_d:
            bob_bit = 0
        elif click_d and not click_c:
            bob_bit = 1
        elif click_c and click_d:
            bob_bit = random.randint(0, 1)
        else:
            bob_bit = -1  # no click — discard

        alice_bits.append(alice_bit)
        alice_bases.append(alice_basis)
        bob_bits.append(bob_bit)
        bob_bases.append(bob_basis)
        has_click.append(click_c or click_d)

        if verbose and (pulse_idx + 1) % max(1, num_bits // 10) == 0:
            print(f"  Pulse {pulse_idx+1}/{num_bits}", flush=True)

    # --- Sifting ---
    sifted_indices = [i for i in range(num_bits)
                      if alice_bases[i] == bob_bases[i] and has_click[i]]
    sifted_alice = [alice_bits[i] for i in sifted_indices]
    sifted_bob = [bob_bits[i] for i in sifted_indices]

    n_sifted = len(sifted_alice)
    n_errors = sum(a != b for a, b in zip(sifted_alice, sifted_bob))
    qber = n_errors / n_sifted if n_sifted > 0 else 0.0

    fiber_loss_dB = alpha_dB * fiber_length
    sifted_key_rate = n_sifted / (num_bits / repetition_rate) if num_bits > 0 else 0.0

    results = {
        'qber': qber,
        'n_total': num_bits,
        'n_sifted': n_sifted,
        'n_errors': n_errors,
        'sifted_key_rate': sifted_key_rate,
        'fiber_length_km': fiber_length,
        'total_loss_dB': fiber_loss_dB,
        'mu': mu,
    }

    if verbose:
        n_clicks = sum(has_click)
        print(f"\nTime-bin BB84 — {fiber_length} km")
        print(f"  Fibre loss: {fiber_loss_dB:.1f} dB")
        print(f"  Mu: {mu} photons/pulse")
        print(f"  Total pulses: {num_bits}")
        print(f"  Clicks: {n_clicks} ({n_clicks/num_bits*100:.2f}%)")
        print(f"  Sifted: {n_sifted} ({n_sifted/num_bits*100:.2f}%)")
        print(f"  Errors: {n_errors}")
        print(f"  QBER: {qber*100:.2f}%")

    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Time-bin phase-encoding BB84 QKD (Gobby et al. 2004)")
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--bits', type=int, default=10000,
                        help='Number of pulses (default 10k)')
    parser.add_argument('--fiber-length', type=float, default=0,
                        help='Fiber length in km (default 0)')
    parser.add_argument('--mu', type=float, default=0.1,
                        help='Mean photons per pulse (default 0.1)')
    parser.add_argument('--verbose', action='store_true')
    args = parser.parse_args()

    results = simulate_bb84_time_bin(
        num_bits=args.bits, fiber_length=args.fiber_length,
        mu=args.mu, seed=args.seed, verbose=True)
