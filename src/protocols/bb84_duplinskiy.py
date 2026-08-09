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
                              bias_offset_v=0.0,
                              spad_eta=0.10, dead_time=13e-6,
                              dark_count_rate=15.0, afterpulse_prob=0.05,
                              cd=False, pmd=False,
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
    spad_eta : float — SPAD quantum efficiency (default 0.10, ID230 [2]).
    dead_time : float — SPAD dead time in seconds (default 13e-6, ID230).
    dark_count_rate : float — SPAD dark count rate in Hz (default 15.0,
        ID230).  Unlike the Gobby replication, ID230 is *contemporary*
        with this 2017 experiment, so the datasheet figures are era-
        appropriate rather than anachronistic (contrast §17.2b).
    afterpulse_prob : float — afterpulse probability per click (default
        0.05, ID230).  §19.5 showed this value was standing in for the
        modulation error in the Gobby chain and belonged at 0 there; that
        argument is specific to Gobby's stated error budget and does not
        transfer.  Exposed here so the question can be asked rather than
        assumed — see DUPL-1.
    cd : bool — apply chromatic dispersion (default False).
    pmd : bool — apply polarisation-mode dispersion (default False).
        Both were hardcoded off before DUPL-1.  This chain is the one
        worth sweeping them on: a polarisation-encoded observable responds
        to them, where the time-bin chain is invariant by construction
        (§26.6).
    bias_offset_v : float — static bias error on Alice's phase modulator,
        as a drive voltage (default 0 = perfectly biased).  Converted by
        `PhaseModulator` through its crystal-derived V_pi.  Setting a
        modulator's bias imperfectly is universal to phase-modulated QKD,
        not specific to any one experiment.
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
    pm_alice = PhaseModulator(crystal_cut='X', modulation='DC',
                              bias_offset_v=bias_offset_v)
    pm_bob = PhaseModulator(crystal_cut='X', modulation='DC')
    Vpi = pm_alice.Vpi

    # SPAD detectors.  These were four hardcoded literals until DUPL-1;
    # they are parameters now so the detector budget can be interrogated
    # the way Gobby's was.  Defaults are the ID230 datasheet figures [2]
    # and reproduce the previous behaviour exactly.
    spd1 = spad(wavelength=wavelength, quantum_efficiency=spad_eta,
                dead_time=dead_time, dark_count_rate=dark_count_rate,
                afterpulse_prob=afterpulse_prob, gate_width=gate_width)
    spd2 = spad(wavelength=wavelength, quantum_efficiency=spad_eta,
                dead_time=dead_time, dark_count_rate=dark_count_rate,
                afterpulse_prob=afterpulse_prob, gate_width=gate_width)

    dt_pulse = 1.0 / rep_rate

    # One physical fibre for the whole run: birefringence is quasi-static,
    # so every pulse must see the same Jones matrix (see ROOT-1 in
    # opto-sim-issues-and-fixes.md). Built once here, reused per pulse
    # below.  `cd` and `pmd` default off, matching this script's original
    # explicit propagate(cd=False, pmd=False) call; they are exposed as of
    # DUPL-1 so this chain can carry an impairment sweep.  It is the right
    # chain for one: polarisation encoding is *not* blind to birefringence
    # the way the time-bin chain is (§26.6).
    fibre = FiberRealization(L_m=fiber_length * 1000, temperature=25,
                             bend_radius=None, attenuation_factor=alpha_dB,
                             cd=cd, pmd=pmd, model=model, seed=seed)

    # Bob's polarization compensation (default on): the channel's Jones
    # matrix is unitary, so the inverse is its conjugate transpose.
    # Applied per pulse before decoding, mirroring the paper's
    # calibration loop (see module docstring).
    J_channel = fibre.birefringence_matrix()
    U_comp = None if J_channel is None else J_channel.conj().T

    # --- Precomputed response table (PERF-2 analogue, DUPL-1) -------------
    #
    # The field chain is DETERMINISTIC given (alice_basis, alice_bit,
    # bob_basis) -- 8 combinations -- because the fibre Jones matrix is
    # sampled once for the whole run (quasi-static, ROOT-1) and no stage
    # between the source and the detectors consumes randomness.  Verified:
    # `pm.modulate`, `fibre.apply`, `optics.voa` and `circular_analyser`
    # all leave both RNG streams untouched, and `fibre.apply` is repeatable
    # on identical input.
    #
    # So the whole per-pulse chain -- modulate, propagate, compensate, VOA,
    # modulate, analyse -- collapses to a table lookup.  This is exactly
    # PERF-2's argument (see `bb84_time_bin.py`), and it is what makes a
    # target-sifted polarisation sweep affordable: the chain measured
    # ~30,000 pulses/s walking the fields per pulse.
    #
    # The table is built by running the REAL chain, never by re-deriving
    # the physics by hand, so there is no second expression of it to drift
    # out of step.
    #
    # NOTE: this is exact only while the response stays deterministic.  A
    # per-pulse random phase (`phase_noise_rad` on either modulator) would
    # break it, exactly as it broke PERF-2's 8-outcome form in GOBBY-6.
    # Neither modulator is given one here.
    def _response(a_basis, a_bit, b_basis):
        """Run the full field chain once; return the gated (P_x, P_y)."""
        if a_basis == 'C':
            v_a = Vpi / 2 if a_bit == 0 else 3 * Vpi / 2
        else:
            v_a = 0 if a_bit == 0 else Vpi

        E = np.sqrt(power_per_pulse / 2.0) * np.ones((1, 2), dtype=complex)
        E = pm_alice.modulate(E_field=E, V=v_a)
        E = fibre.apply(E, dt=dt_pulse)
        if compensate and U_comp is not None:
            E = np.transpose(U_comp @ np.transpose(E))
        E = optics.voa(E, bob_loss_dB)

        v_b = 0 if b_basis == 'C' else Vpi / 2
        E = pm_bob.modulate(E_field=E, V=v_b)
        # circular_analyser, not pbs: detection depends on the relative
        # phase between Ex/Ey (PHYS-6 in opto-sim-issues-and-fixes.md).
        Ex, Ey = optics.circular_analyser(E)
        return (float(np.mean(np.abs(Ex) ** 2)),
                float(np.mean(np.abs(Ey) ** 2)))

    RESPONSE = {(ab, bit, bb): _response(ab, bit, bb)
                for ab in ('C', 'X')
                for bit in (0, 1)
                for bb in ('C', 'X')}

    # Sifting is accumulated inline rather than into per-pulse lists.  This
    # previously kept FIVE of them -- alice_bits, alice_bases, bob_bits,
    # bob_bases, has_click -- each O(num_bits), then sifted by list
    # comprehension.  Fine at the 100k runs this was written for, several GB
    # at the ~30e6 pulses per row a target-sifted polarisation sweep needs.
    # Same defect the Gobby chain carried until §17.5.  Counting inline is
    # O(1) and touches no RNG draw, so results are bit-identical.
    n_sifted = 0
    n_errors = 0
    n_clicks = 0

    t_start = time.time()

    for pulse_idx in range(num_bits):
        # --- Alice's encoding ---
        alice_basis = random.choice(['C', 'X'])
        alice_bit = random.randint(0, 1)

        # --- Bob's decoding basis ---
        bob_basis = random.choice(['C', 'X'])

        # --- Response from the precomputed table (see above) ---
        power_x, power_y = RESPONSE[(alice_basis, alice_bit, bob_basis)]

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

        if click_x or click_y:
            n_clicks += 1
            if alice_basis == bob_basis:
                n_sifted += 1
                if alice_bit != bob_bit:
                    n_errors += 1

        if verbose and (pulse_idx + 1) % 10000 == 0:
            elapsed = time.time() - t_start
            rate = (pulse_idx + 1) / elapsed
            print(f"  Pulse {pulse_idx+1}/{num_bits} ({rate:.0f} pulses/s)", flush=True)

    # Sifting (same basis AND at least one click) was accumulated in the
    # loop above.
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
