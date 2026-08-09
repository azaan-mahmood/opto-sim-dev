"""Tests for the protocol-level scripts (src/protocols/).

Currently covers the ARCH-3 polarization-compensation regression: the
`compensate` flag of `bb84_duplinskiy` must actually disable the inverse-
Jones correction, and the correction must exist at all (a fibre with a
scrambling quasi-static Jones matrix must raise the QBER when it is
disabled).
"""
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
