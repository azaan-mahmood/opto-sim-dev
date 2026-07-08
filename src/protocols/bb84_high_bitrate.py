import argparse
import numpy as np
from src.lasers.cwlaser import CWLaser
from src.channel import cable, optics
from src.channel.phase_modulator import PhaseModulator
from src.visualization import fields, polarimeter
from src.detectors import apd
import random


def simulate_bb84_high_bitrate(num_bits, fiber_length=100, bandwidth=1e9,
                               dispersion=False, show_pol=False, seed=None):
    """
    Simulation for a BB84 protocol using realistic detector noise.

    Parameters
    ----------
    seed : int or None
        RNG seed for reproducibility. If None, uses default seeding.
    """
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)

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
    Vpi = pm_alice.Vpi

    dt = 1e-12     # time step (1 ps) — 1 ns total per bit at 1000 samples
    n_samples = 1000

    for _ in range(num_bits):
        # Alice's random choices
        alice_basis = random.choice(['C', 'X'])
        alice_bit = random.randint(0, 1)

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

        # Generate field for one bit (complex envelope, 1 ns at 1 Gbaud)
        E = alice_laser.sample_field(dt=dt, n_samples=n_samples)

        # Apply Alice Phase Modulation
        E = optics.polarizer(E, '45')  # Initial polarization
        E = pm_alice.modulate(E_field=E, V=phase_alice)  # Alice's phase shift

        if show_pol:
            polarimeter(E, title=f"Bit Number {num_bits},Alice Bit/Basis = {alice_bit} / {alice_basis}")

        # Channel transmission (QC). Dispersion now works because
        # sample_field returns the complex envelope (not the optical carrier).
        E = cable(
            fiber_length=fiber_length, E=E, dt=dt, dispersion=dispersion,
            attenuation_factor=0.182, temperature=25, bend_radius=None
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
        I_x = detector.output(E=Ex, bandwidth=bandwidth)
        I_y = detector.output(E=Ey, bandwidth=bandwidth)

        # Differential detection: ensure signal exceeds noise floor,
        # then compare photocurrents to determine the bit.
        noise_floor = detector.calculate_noise(0, bandwidth)
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
        sample_size = min(len(sifted_alice), 100)
        if sample_size > 0:
            error_count = 0
            for i in random.sample(range(len(sifted_alice)), sample_size):
                if sifted_alice[i] != sifted_bob[i]:
                    error_count += 1
            qber = error_count / sample_size

    return qber


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='BB84 high-bitrate simulation')
    parser.add_argument('--seed', type=int, default=None,
                        help='RNG seed for reproducibility')
    parser.add_argument('--num-bits', type=int, default=500,
                        help='Number of bits to simulate')
    parser.add_argument('--fiber-length', type=float, default=100,
                        help='Fiber length in km')
    parser.add_argument('--bandwidth', type=float, default=1e9,
                        help='Detector bandwidth in Hz')
    parser.add_argument('--dispersion', action='store_true',
                        help='Enable chromatic dispersion and PMD')
    args = parser.parse_args()

    qber = simulate_bb84_high_bitrate(
        args.num_bits, args.fiber_length, bandwidth=args.bandwidth,
        dispersion=args.dispersion, seed=args.seed
    )
    print(f"QBER: {qber}")
