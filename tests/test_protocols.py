"""Tests for the protocol-level scripts (src/protocols/).

Currently covers the ARCH-3 polarization-compensation regression: the
`compensate` flag of `bb84_duplinskiy` must actually disable the inverse-
Jones correction, and the correction must exist at all (a fibre with a
scrambling quasi-static Jones matrix must raise the QBER when it is
disabled).
"""
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
