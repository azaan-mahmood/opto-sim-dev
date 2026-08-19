"""Tests for the protocol-level scripts (src/protocols/).

Currently covers the ARCH-3 polarization-compensation regression: the
`compensate` flag of `bb84_duplinskiy` must actually disable the inverse-
Jones correction, and the correction must exist at all (a fibre with a
scrambling quasi-static Jones matrix must raise the QBER when it is
disabled).
"""
import math
import random

import numpy as np
import pytest
from src.protocols.bb84_duplinskiy import simulate_bb84_duplinskiy


class TestDuplinskiyCompensation:

    def test_no_compensation_scrambles_encoding(self):
        """At 10 km the sectional model draws a near-uniform SU(2)
        rotation; with compensation off, that rotation must scramble the
        phase encoding (QBER jumps well above the ~1.5% floor), and with
        compensation on it must be undone. Fixed seed → deterministic.
        (Regression for the ARCH-3 control arm: the flag was silently
        ignored, so '--no-compensation' runs were bit-identical to
        compensated ones.)"""
        comp = simulate_bb84_duplinskiy(
            num_bits=20000, fiber_length=10, mu=2.0, bob_loss_dB=0.0,
            model='sectional', compensate=True, seed=42)
        uncomp = simulate_bb84_duplinskiy(
            num_bits=20000, fiber_length=10, mu=2.0, bob_loss_dB=0.0,
            model='sectional', compensate=False, seed=42)

        assert comp["n_sifted"] > 50
        assert uncomp["n_sifted"] > 50
        assert comp["qber"] < 0.05
        assert uncomp["qber"] > comp["qber"] + 0.10


class TestDuplinskiyResponseTable:
    """The 8-outcome precompute must be exact, not merely fast.

    The field chain is deterministic given `(alice_basis, alice_bit,
    bob_basis)` because the fibre Jones matrix is sampled once per run
    (quasi-static, ROOT-1) and no stage between source and detectors
    consumes randomness.  So the whole per-pulse chain collapses to a table
    lookup — PERF-2's argument, applied to the polarisation chain (DUPL-1).

    These pin the two things that make that legitimate: the stages really
    are RNG-free, and the table really does reproduce the walked chain.
    """

    @staticmethod
    def _rng_state():
        return np.random.get_state()[2], random.getstate()[1][0]

    def test_field_stages_consume_no_randomness(self):
        """If any stage drew from the RNG, precomputing would shift the
        stream and silently change every downstream detector draw."""
        from src.channel import optics as _optics
        from src.channel.fiber import FiberRealization
        from src.channel.phase_modulator import PhaseModulator
        pm = PhaseModulator(crystal_cut='X', modulation='DC')
        f = FiberRealization(L_m=50e3, temperature=25, bend_radius=None,
                             attenuation_factor=0.2, cd=False, pmd=False,
                             model='auto', seed=42)
        E = np.sqrt(1e-9 / 2) * np.ones((1, 2), dtype=complex)
        for fn in (lambda: pm.modulate(E_field=E.copy(), V=1.0),
                   lambda: f.apply(E.copy(), dt=1e-7),
                   lambda: _optics.voa(E.copy(), 2.0),
                   lambda: _optics.circular_analyser(E.copy())):
            before = self._rng_state()
            fn()
            assert self._rng_state() == before

    def test_fibre_apply_is_repeatable(self):
        """The table is built once and reused, so the channel must give the
        same field every time it is asked."""
        from src.channel.fiber import FiberRealization
        f = FiberRealization(L_m=50e3, temperature=25, bend_radius=None,
                             attenuation_factor=0.2, cd=False, pmd=False,
                             model='auto', seed=42)
        E = np.sqrt(1e-9 / 2) * np.ones((1, 2), dtype=complex)
        assert np.array_equal(f.apply(E.copy(), dt=1e-7),
                              f.apply(E.copy(), dt=1e-7))

    @pytest.mark.parametrize("length,compensate,qber,sifted,errors", [
        (0, True, 0.0120481928, 83, 1),
        (0, False, 0.0120481928, 83, 1),
        (10, True, 0.0273972603, 73, 2),
        (10, False, 0.4833333333, 60, 29),
        (50, True, 0.0000000000, 14, 0),
        (50, False, 0.3125000000, 16, 5),
    ])
    def test_matches_the_walked_chain(self, length, compensate, qber,
                                      sifted, errors):
        """Frozen from the pre-optimisation implementation, which walked the
        full field chain per pulse.  The precompute is a performance change
        and must not move a single bit."""
        from src.protocols.bb84_duplinskiy import simulate_bb84_duplinskiy
        r = simulate_bb84_duplinskiy(40000, fiber_length=length,
                                     compensate=compensate, seed=42)
        assert r['n_sifted'] == sifted
        assert r['n_errors'] == errors
        assert r['qber'] == pytest.approx(qber, abs=1e-9)

    def test_the_equivalence_gate_can_fail(self):
        """Negative control.  A gate that cannot fail proves nothing — the
        lesson from GOBBY-6 §23.1, where flipping Bob's sign returned the
        passing value.  Swapping the analyser outputs is a real physics
        change and must break the frozen numbers above."""
        from src.channel import optics as _optics
        import src.protocols.bb84_duplinskiy as _dupl
        original = _optics.circular_analyser

        def swapped(E):
            a, b = original(E)
            return b, a

        _dupl.optics.circular_analyser = swapped
        try:
            r = _dupl.simulate_bb84_duplinskiy(40000, fiber_length=10,
                                               compensate=True, seed=42)
        finally:
            _dupl.optics.circular_analyser = original
        assert r['qber'] != pytest.approx(0.0273972603, abs=1e-9)


class TestDuplinskiyFibreDrift:
    """The same fibre clock on the polarisation chain.

    `calibration_temperature` / `calibration_bend_radius` already give a
    fixed two-state mismatch. Drift is the time dependence they cannot
    express: the residual grows during the run, which is what actually
    forces the paper's recalibrations.
    """

    BASE = dict(num_bits=20_000, fiber_length=50, mu=2.0, bob_loss_dB=0.0,
                model='sectional', seed=42)

    @staticmethod
    def _key(r):
        return r['qber'], r['n_sifted'], r['n_errors']

    def test_drift_blocks_do_not_perturb_a_static_run(self):
        """The drift partition must be inert when nothing drifts, and must
        stay independent of `block_size`, which slices the run for
        reporting rather than for physics."""
        ref = simulate_bb84_duplinskiy(compensate=True, **self.BASE)
        for blocks in (1, 250):
            got = simulate_bb84_duplinskiy(compensate=True,
                                           drift_blocks=blocks, **self.BASE)
            assert self._key(got) == self._key(ref)

    def test_drift_degrades_qber(self):
        """Bob aligns at t=0 and the fibre walks away from it."""
        ref = simulate_bb84_duplinskiy(compensate=True, **self.BASE)
        drift = simulate_bb84_duplinskiy(
            compensate=True, run_duration=120.0,
            drift_temperature_rate_C_s=1e-3, drift_blocks=50, **self.BASE)
        assert drift['n_sifted'] > 50
        assert drift['qber'] > ref['qber'] + 0.15

    def test_run_duration_decouples_drift_from_the_pulse_budget(self):
        """Two budgets sampling the same 120 s experiment must agree on the
        QBER, or asking for tighter error bars would silently be asking for
        a longer experiment -- the bug `run_duration` exists to prevent."""
        common = dict(compensate=True, run_duration=120.0,
                      drift_temperature_rate_C_s=1e-3, drift_blocks=50,
                      fiber_length=50, mu=2.0, bob_loss_dB=0.0,
                      model='sectional', seed=42)
        short = simulate_bb84_duplinskiy(num_bits=20_000, **common)
        long = simulate_bb84_duplinskiy(num_bits=60_000, **common)
        # Both estimate the same expectation; 3x the pulses must not move it
        # beyond the binomial spread of the smaller run.
        spread = 3.0 * math.sqrt(
            short['qber'] * (1 - short['qber']) / short['n_sifted'])
        assert abs(short['qber'] - long['qber']) < spread


class TestTimeBinFibreDrift:
    """`bb84_time_bin`'s fibre and its clock.

    The time-bin chain models the channel with a closed form that fixes
    the interference coefficients for one fibre state, so a drifting fibre
    has to be blocked. These pin the properties that makes safe.
    """

    # Gobby's link budget at a short distance, sized for speed rather than
    # for statistical power -- these test invariants, not physics values.
    BASE = dict(num_bits=60_000, fiber_length=10, alpha_dB=0.2, mu=0.1,
                repetition_rate=2e6, pulse_width=80e-12, spad_eta=0.045,
                dark_count_rate=8.5e-7 / 3.5e-9, afterpulse_prob=0.0,
                dead_time=13e-6, gate_width=3.5e-9, split_ratio=1.6,
                interferometer='polarisation_multiplexed',
                run_duration=120.0, seed=42)

    @staticmethod
    def _key(r):
        return r['qber'], r['n_sifted'], r['n_errors']

    def test_impairments_off_is_bit_identical(self):
        """The whole fibre block must be inert when every flag is off, or
        the existing Gobby sweep stops being comparable to what has run."""
        from src.protocols.bb84_time_bin import simulate_bb84_time_bin
        plain = simulate_bb84_time_bin(**self.BASE)
        explicit = simulate_bb84_time_bin(birefringence=False, cd=False,
                                          pmd=False, **self.BASE)
        assert self._key(plain) == self._key(explicit)

    def test_static_fibre_with_alignment_is_an_exact_null(self):
        """U_comp @ J = I for any unitary, so this null is arithmetic and
        must be reported as such -- it is not evidence about physics."""
        from src.protocols.bb84_time_bin import simulate_bb84_time_bin
        ref = simulate_bb84_time_bin(**self.BASE)
        aligned = simulate_bb84_time_bin(birefringence=True, compensate=True,
                                         **self.BASE)
        assert self._key(aligned) == self._key(ref)

    def test_uncompensated_static_fibre_cuts_the_rate(self):
        """The Jones matrix is SU(2), so both interfering arms are scaled by
        the same |U00| and the sifted rate falls as |U00|^2 -- Gobby's
        "polarisation drift reduces the bit rate". A run that shows no rate
        change has not applied the impairment."""
        from src.channel.fiber import _build_jones_matrix
        from src.protocols.bb84_time_bin import simulate_bb84_time_bin
        J = _build_jones_matrix(np.random.default_rng(42), 10_000.0)
        ref = simulate_bb84_time_bin(**self.BASE)
        unc = simulate_bb84_time_bin(birefringence=True, compensate=False,
                                     **self.BASE)
        assert unc['n_sifted'] < ref['n_sifted']
        # Predicted rate factor, generous tolerance for binomial noise at
        # this deliberately small pulse budget.
        assert unc['n_sifted'] / ref['n_sifted'] == pytest.approx(
            abs(J[0, 0]) ** 2, rel=0.35)

    def test_drift_blocks_do_not_perturb_a_static_run(self):
        """`drift_blocks` is ignored when nothing drifts. Without this the
        knob would silently change results on runs that have no drift."""
        from src.protocols.bb84_time_bin import simulate_bb84_time_bin
        a = simulate_bb84_time_bin(birefringence=True, drift_blocks=1,
                                   **self.BASE)
        b = simulate_bb84_time_bin(birefringence=True, drift_blocks=500,
                                   **self.BASE)
        assert self._key(a) == self._key(b)

    def test_residual_after_alignment_is_mostly_phase(self):
        """The deterministic half of the claim, asserted on the operator
        rather than through Monte Carlo noise.

        Bob aligns at t=0, so what light meets is the residual
        R = U_comp @ J(t). R is SU(2) like its factors, so it splits the
        same way: |R00|^2 scales both arms (rate) and 2*arg(R11) separates
        them (QBER). At these rates the phase runs far ahead of the
        amplitude, which is why drift shows up in QBER first."""
        from src.channel.fiber import FiberRealization
        fibre = FiberRealization(L_m=10_000, seed=42, attenuation=False,
                                 drift_temperature_rate_C_s=1e-3)
        U = fibre.birefringence_matrix().conj().T
        R = U @ fibre.at(120.0).birefringence_matrix()
        implied_qber = (1.0 - np.cos(2.0 * np.angle(R[1, 1]))) / 2.0
        assert abs(R[0, 0]) ** 2 > 0.8          # rate largely spared
        assert implied_qber > 0.05              # QBER clearly moved

    def test_drift_degrades_qber(self):
        """The Monte Carlo half. Only the QBER is asserted here: the effect
        size is large, whereas the rate claim above needs more sifted bits
        than this deliberately small budget produces and is pinned
        deterministically instead."""
        from src.protocols.bb84_time_bin import simulate_bb84_time_bin
        ref = simulate_bb84_time_bin(**self.BASE)
        drift = simulate_bb84_time_bin(
            birefringence=True, compensate=True,
            drift_temperature_rate_C_s=1e-3, drift_blocks=50, **self.BASE)
        assert drift['qber'] > ref['qber'] + 0.05
        # Distinguishes "rate largely spared" from the uncompensated case,
        # which lands near 0.19 at this distance.
        assert drift['n_sifted'] / ref['n_sifted'] > 0.7

    def test_drift_without_birefringence_raises(self):
        """Drift acts only through the Jones matrix, so this combination
        would silently do nothing."""
        from src.protocols.bb84_time_bin import simulate_bb84_time_bin
        with pytest.raises(ValueError, match="birefringence"):
            simulate_bb84_time_bin(birefringence=False, cd=True,
                                   drift_temperature_rate_C_s=1e-3,
                                   **self.BASE)


class TestTimeBinPhaseServo:
    """Bob's phase servo: the piezo stretcher that holds the operating point.

    Without it the model contradicts the paper it replicates. A residual
    rotation is SU(2), so it costs O(eps^2) in amplitude against O(eps) in
    phase, and QBER therefore moves before the rate does -- the reverse of
    "polarisation drift reduces the bit rate, but does not degrade the
    QBER". The servo is what makes the ordering come out right.
    """

    BASE = dict(num_bits=200_000, fiber_length=10, alpha_dB=0.2, mu=0.1,
                repetition_rate=2e6, pulse_width=80e-12, spad_eta=0.045,
                dark_count_rate=8.5e-7 / 3.5e-9, afterpulse_prob=0.0,
                dead_time=13e-6, gate_width=3.5e-9, split_ratio=1.6,
                interferometer='polarisation_multiplexed',
                run_duration=120.0, seed=42)
    # 3e-3 C/s, not the 1e-3 used elsewhere: at the lower rate the amplitude
    # loss is inside the counting noise, so a test there could not tell
    # "corrects phase only" from "corrects everything".
    DRIFT = dict(birefringence=True, compensate=True,
                 drift_temperature_rate_C_s=3e-3, drift_blocks=100)

    @staticmethod
    def _key(r):
        return r['qber'], r['n_sifted'], r['n_errors']

    def _run(self, **over):
        from src.protocols.bb84_time_bin import simulate_bb84_time_bin
        return simulate_bb84_time_bin(**{**self.BASE, **over})

    def test_servo_off_is_inert(self):
        """None must mean exactly what leaving it out means."""
        assert self._key(self._run(**self.DRIFT)) == self._key(
            self._run(phase_servo_interval_s=None, **self.DRIFT))

    def test_servo_without_drift_does_nothing(self):
        """Nothing to track, so nothing to correct."""
        a = self._run(birefringence=True, compensate=True)
        b = self._run(birefringence=True, compensate=True,
                      phase_servo_interval_s=1.0)
        assert self._key(a) == self._key(b)

    def test_the_servo_signal_is_the_residual_phase(self):
        """The claim the whole thing rests on.

        `arg(S)` shifts by exactly the phase the fibre puts between the
        interfering arms, so the servo reads its error signal off an
        extraction that already happens. With one block the fibre is frozen
        at the run midpoint, so the correction is a single known number and
        the servo must be indistinguishable from cancelling it by hand.
        """
        from src.channel.fiber import FiberRealization
        fibre = FiberRealization(L_m=10_000.0, seed=42, attenuation=False,
                                 drift_temperature_rate_C_s=3e-3)
        R = (fibre.birefringence_matrix().conj().T
             @ fibre.at(60.0).birefringence_matrix())
        theta = 2.0 * np.angle(R[1, 1])
        frozen = dict(self.DRIFT, drift_blocks=1)
        servo = self._run(phase_servo_interval_s=0.0, **frozen)
        by_hand = self._run(phase_error_rad=-theta, **frozen)
        assert self._key(servo) == self._key(by_hand)

    def test_restores_qber_without_restoring_rate(self):
        """Gobby's sentence, as a test.

        The rate loss must survive the servo. If it does not, the servo is
        correcting amplitude as well and is the wrong model of a device
        that only turns a phase.
        """
        ref = self._run()
        unserved = self._run(**self.DRIFT)
        served = self._run(phase_servo_interval_s=1.2, **self.DRIFT)

        # QBER: wrecked without, back to baseline with
        assert unserved['qber'] > 0.20
        assert served['qber'] < ref['qber'] + 0.02

        # Rate: lost in both, and lost by the same amount
        r_un = unserved['n_sifted'] / ref['n_sifted']
        r_sv = served['n_sifted'] / ref['n_sifted']
        assert r_un < 0.8 and r_sv < 0.8
        assert abs(r_sv - r_un) < 0.10

    def test_a_slower_servo_tracks_worse(self):
        """Finite bandwidth is the point of the interval. A servo that
        re-locks rarely lets the residual grow between locks, so the QBER
        climbs back towards the unserved value."""
        fast = self._run(phase_servo_interval_s=1.2, **self.DRIFT)
        slow = self._run(phase_servo_interval_s=60.0, **self.DRIFT)
        unserved = self._run(**self.DRIFT)
        assert fast['qber'] < slow['qber'] < unserved['qber']
