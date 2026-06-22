import numpy as np
from src.lasers.cwlaser import CWLaser
from src.opto_eq import cable, optics
from src.opto_eq.phase_modulator import PhaseModulator
from src.viewers import fields, polarimeter
from src.detectors import apd
import random


def simulate_bb84(num_bits, fiber_length=100, show_pol = False):
    """
    Simulation for an ideal BB84, that does not contain Eve.
    # Example usage
    # num_bits = 1000
    # params = simulate_bb84(num_bits, show_pol = False)
    # print(f"Quantum Bit Error Rate (QBER): {params[7] * 100:.2f}%")
    :param num_bits: Number of bits for the simulation
    :param show_pol: Show polarization for each bit.
    :return: list containing [alice_bits, alice_bases, bob_bits, bob_bases, sifted_indices, sifted_alice, sifted_bob, qber]


    """
    # Initialize lists to store Alice's and Bob's data
    alice_bits, alice_bases = [], []
    bob_bits, bob_bases = [], []

    detector = apd(
        wavelength=1550e-9, quantum_efficiency=0.9, gain=10, excess_noise_factor=10,
        load_resistance=50, temperature=25
    )
    alice_laser = CWLaser(
        wavelength=1550e-9,
        polarization_azimuth=np.pi,
        polarization_ellipticity=np.pi / 2,
        power_dbm=-5,
        linewidth=1e6,
        rin_density=-140,
    )

    pm_alice = PhaseModulator(crystal_cut='X', modulation="DC")
    pm_bob = PhaseModulator(crystal_cut='X', modulation="DC")
    # Use the modulator's own Vpi so that the phase voltages match its
    # crystal-calculated half-wave voltage.
    Vpi = pm_alice.Vpi

    for _ in range(num_bits):
        # Alice's random choices
        alice_basis = random.choice(['C', 'X'])
        alice_bit = random.randint(0, 1)
        # phase_alice = alice_phase_shift(alice_basis, alice_bit)
        if alice_basis == 'C':
            if alice_bit == 0:
                phase_alice = Vpi / 2
            else:
                phase_alice = 3 * Vpi / 2
        else:  # 'X' basis
            if alice_bit == 0:
                phase_alice = 0
            else:
                phase_alice = Vpi

        # Generate photon with Alice's encoding
        E = alice_laser.get_electric_field(normalize=False, over_period=True)
        # mean(|E|²) is already calibrated to optical power in Watts (CWLaser convention)

        # Apply Alice Phase Modulation
        E = optics.polarizer(E, '45')  # Initial polarization
        E = pm_alice.modulate(E_field=E, V=phase_alice)  # Alice's phase shift

        if show_pol:
            polarimeter(E, title=f"Bit Number {num_bits},Alice Bit/Basis = {alice_bit} / {alice_basis}")

        # Channel transmission (QC). Dispersion disabled — the FFT-based
        # dispersion model produces unphysical frequencies ~10^17 Hz and
        # corrupts the field at optical carrier rates (known blocked issue).
        E = cable(
            fiber_length=fiber_length, E=E, dispersion=False,
            attenuation_factor=0.182, temperature=25, num_bends=10
        )

        # Bob's random basis choice
        bob_basis = random.choice(['C', 'X'])
        if bob_basis == 'C':
            phase_bob = 0
        else:
            phase_bob = Vpi / 2

        # Apply Bob Phase Modulation
        E = pm_bob.modulate(E_field=E, V=phase_bob)  # Bob's phase shift
        # Measurement: PBS splits into two spatial modes
        Ex, Ey = optics.pbs(E)

        # Noisy photocurrent from each detector (power derived from field)
        I_x = detector.output(E=Ex, bandwidth=1e6)
        I_y = detector.output(E=Ey, bandwidth=1e6)

        # Differential detection: ensure signal exceeds noise floor,
        # then compare photocurrents to determine the bit.
        noise_floor = detector.calculate_noise(0, 1e6)
        threshold = 3 * noise_floor

        if I_x > threshold or I_y > threshold:
            bob_bit = 0 if I_x > I_y else 1
        else:
            bob_bit = random.randint(0, 1)

        if show_pol:
            polarimeter(E, title=f"Bit Number {num_bits},Bob Bit/Basis = {bob_bit} / {bob_basis}")
        # Store results
        alice_bits.append(alice_bit)
        alice_bases.append(alice_basis)
        bob_bits.append(bob_bit)
        bob_bases.append(bob_basis)

    # Sift keys to retain matching bases
    sifted_indices = [i for i in range(num_bits) if alice_bases[i] == bob_bases[i]]
    sifted_alice = [alice_bits[i] for i in sifted_indices]
    sifted_bob = [bob_bits[i] for i in sifted_indices]

    # Calculate QBER using a portion of the sifted key
    qber = 0.0
    if len(sifted_alice) > 0:
        sample_size = min(len(sifted_alice) // 2, 100)  # Check up to 100 bits
        if sample_size > 0:
            error_count = 0
            for i in random.sample(range(len(sifted_alice)), sample_size):
                if sifted_alice[i] != sifted_bob[i]:
                    error_count += 1
            qber = error_count / sample_size

    return [alice_bits, alice_bases, bob_bits, bob_bases, sifted_indices, sifted_alice, sifted_bob, qber]


if __name__ == "__main__":
    alice_bits, alice_bases, bob_bits, bob_bases, sifted_indices, sifted_alice, sifted_bob, qber=simulate_bb84(500, 100)
    print(f"Alice Bits: {alice_bits}")
    print(f"Alice Bases: {alice_bases}")
    print(f"Bob Bits: {bob_bits}")
    print(f"Bob Bases: {bob_bases}")
    print(f"Sifted Indices: {sifted_indices}")
    print(f"Sifted Alice: {sifted_alice}")
    print(f"Sifted Bob: {sifted_bob}")
    print(f"QBER: {qber}")