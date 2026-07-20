import argparse
import numpy as np
from src.lasers.cwlaser import CWLaser
from src.channel import propagate, optics
from src.channel.phase_modulator import PhaseModulator
from src.channel.mzm import MZM
from src.visualization import fields, polarimeter
from src.detectors import apd
import random


def simulate_bb84_dispersion(num_bits, fiber_length=100, pulse_sigma=30e-12,
                              dispersion=True, show_pol=False, seed=None):
    """
    BB84 with MZM-carved pulses to make CD and PMD measurable.

    Unlike the ideal BB84 (which uses a near-monochromatic CW field), this
    version carves short Gaussian-like pulses via an MZM, giving the field
    enough bandwidth for chromatic dispersion and PMD to induce QBER.

    Parameters
    ----------
    pulse_sigma : float
        Standard deviation of the Gaussian voltage pulse driving the MZM (s).
        Smaller values → broader bandwidth → stronger CD/PMD effects.
        Default 30e-12 (30 ps, FWHM ~71 ps).
    dispersion : bool
        Enable chromatic dispersion and PMD in the fiber.
        Default True (this is the whole point of the test).
    """
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)

    alice_bits, alice_bases = [], []
    bob_bits, bob_bases = [], []

    detector = apd(
        wavelength=1550e-9, quantum_efficiency=0.9, gain=10,
        excess_noise_factor=10, load_resistance=50, temperature=25
    )

    # Higher power to compensate for the low duty cycle of pulsed operation
    alice_laser = CWLaser(
        wavelength=1550e-9,
        polarization_azimuth=np.pi,
        polarization_ellipticity=np.pi / 2,  # Ey-only output
        power_dbm=10,
        linewidth=1e6,
        rin_density=-140,
    )

    pm_alice = PhaseModulator(crystal_cut='X', modulation="DC")
    pm_bob = PhaseModulator(crystal_cut='X', modulation="DC")
    Vpi = pm_alice.Vpi

    dt = 1e-12
    n_samples = 1000
    t = np.arange(n_samples) * dt
    t_center = n_samples * dt / 2

    # MZM carves a pulse: bias at V_pi so V=0 → null, V=V_pi → full tx
    mzm = MZM(mode='push-pull', bias_voltage=Vpi)

    # Gaussian voltage pulse — carves a short optical pulse from the CW field
    V_pulse = Vpi * np.exp(-0.5 * ((t - t_center) / pulse_sigma) ** 2)

    for _ in range(num_bits):
        alice_basis = random.choice(['C', 'X'])
        alice_bit = random.randint(0, 1)

        if alice_basis == 'C':
            phase_alice = Vpi / 2 if alice_bit == 0 else 3 * Vpi / 2
        else:
            phase_alice = 0 if alice_bit == 0 else Vpi

        # CW field (laser output is Ey-only for this polarization)
        E = alice_laser.sample_field(dt=dt, n_samples=n_samples)

        # Carve a short pulse via MZM (modulates Ey, the only component)
        E = mzm.modulate(E_in=E, V=V_pulse)

        # 45° polarizer spreads pulsed Ey equally into Ex and Ey
        E = optics.polarizer(E, '45')

        # Alice's phase encoding
        E = pm_alice.modulate(E_field=E, V=phase_alice)

        if show_pol:
            polarimeter(E, title=f"Alice: bit {alice_bit}, basis {alice_basis}")

        # Fiber with dispersion (now physically meaningful for broadband pulses)
        E = propagate(
            fiber_length=fiber_length, E=E, dt=dt, dispersion=dispersion,
            attenuation_factor=0.182, temperature=25, bend_radius=None
        )

        bob_basis = random.choice(['C', 'X'])
        phase_bob = 0 if bob_basis == 'C' else Vpi / 2
        E = pm_bob.modulate(E_field=E, V=phase_bob)

        Ex, Ey = optics.pbs(E)
        I_x = detector.output(E=Ex, bandwidth=1e6)
        I_y = detector.output(E=Ey, bandwidth=1e6)

        noise_floor = detector.calculate_noise(0, 1e6)
        threshold = 3 * noise_floor

        if I_x > threshold or I_y > threshold:
            bob_bit = 0 if I_x > I_y else 1
        else:
            bob_bit = random.randint(0, 1)

        alice_bits.append(alice_bit)
        alice_bases.append(alice_basis)
        bob_bits.append(bob_bit)
        bob_bases.append(bob_basis)

    sifted_indices = [i for i in range(num_bits) if alice_bases[i] == bob_bases[i]]
    sifted_alice = [alice_bits[i] for i in sifted_indices]
    sifted_bob = [bob_bits[i] for i in sifted_indices]

    qber = 0.0
    if len(sifted_alice) > 0:
        sample_size = min(len(sifted_alice) // 2, 100)
        if sample_size > 0:
            errors = sum(
                sifted_alice[i] != sifted_bob[i]
                for i in random.sample(range(len(sifted_alice)), sample_size)
            )
            qber = errors / sample_size

    return qber


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='BB84 with MZM-carved pulses for CD/PMD testing'
    )
    parser.add_argument('--seed', type=int, default=None,
                        help='RNG seed')
    parser.add_argument('--num-bits', type=int, default=500,
                        help='Number of bits')
    parser.add_argument('--fiber-length', type=float, default=100,
                        help='Fiber length (km)')
    parser.add_argument('--pulse-sigma', type=float, default=30e-12,
                        help='Gaussian pulse sigma in seconds (default 30 ps)')
    parser.add_argument('--no-dispersion', action='store_true',
                        help='Disable CD and PMD')
    parser.add_argument('--short-pulse', action='store_true',
                        help='Use 5 ps pulse for stronger CD/PMD effects')
    args = parser.parse_args()

    pulse_sigma = 5e-12 if args.short_pulse else args.pulse_sigma
    dispersion = not args.no_dispersion

    qber = simulate_bb84_dispersion(
        args.num_bits, args.fiber_length, pulse_sigma=pulse_sigma,
        dispersion=dispersion, seed=args.seed
    )
    print(f"QBER: {qber * 100:.2f}%")
