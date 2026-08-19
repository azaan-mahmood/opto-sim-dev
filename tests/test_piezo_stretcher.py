"""Tests for `src/channel/piezo_stretcher.py`.

The device is a Thorlabs FVP155P by default, so most of these are checks
that the model still agrees with the datasheet it claims to implement.
"""
import math

import numpy as np
import pytest

from src.channel import PiezoFibreStretcher

TWO_PI = 2.0 * math.pi


class TestDatasheetAgreement:
    """The defaults must stay the part they say they are (DOC-103641 B)."""

    def test_defaults_match_the_datasheet(self):
        st = PiezoFibreStretcher()
        assert st.v_pi == 20.0            # "Half-Wave Voltage  < 20 V"
        assert st.v_max == 150.0          # "Drive Voltage Range  0 - 150 V"
        assert st.stroke_rad == pytest.approx(7.0 * math.pi)
        assert st.resonant_hz == 80e3
        assert st.insertion_loss_db == 0.1
        assert st.residual_am == 0.0015

    def test_the_stroke_is_reachable_within_the_drive_range(self):
        """7*pi at 150 V is the headline spec. With v_pi at the 20 V bound
        the linear law reaches 7.5*pi, so the modelled part covers the
        stroke rather than falling short of it."""
        st = PiezoFibreStretcher()
        assert st.phase_for(st.v_max) >= st.stroke_rad

    def test_a_full_fringe_costs_two_v_pi(self):
        """2*pi is the most any phase correction can need after wrapping,
        so this is the worst-case demand on the part."""
        st = PiezoFibreStretcher()
        assert st.voltage_for(TWO_PI - 1e-12) == pytest.approx(
            2 * st.v_pi, rel=1e-9)
        assert 2 * st.v_pi < st.v_max


class TestVoltagePhaseRoundTrip:

    def test_round_trip_is_exact_modulo_a_fringe(self):
        st = PiezoFibreStretcher()
        for phi in (0.1, 1.0, math.pi, 5.0, 5 * math.pi, -1.0, -math.pi):
            back = st.phase_for(st.voltage_for(phi))
            assert (back - phi) % TWO_PI == pytest.approx(0.0, abs=1e-12) or \
                   (back - phi) % TWO_PI == pytest.approx(TWO_PI, abs=1e-12)

    def test_wrapping_is_what_makes_a_finite_stroke_enough(self):
        """A correction tracking a steadily drifting fibre would grow
        without bound. Wrapping means the device flies back to the nearer
        fringe, which is what a real one does when it runs out of travel."""
        st = PiezoFibreStretcher()
        assert st.voltage_for(5 * math.pi) == pytest.approx(st.v_pi)
        assert st.voltage_for(101 * math.pi) == pytest.approx(st.v_pi)

    def test_over_range_drive_raises(self):
        """The datasheet warns that exceeding 150 V shortens the device's
        life and that reverse bias can destroy it, so this refuses rather
        than clipping to the limit."""
        st = PiezoFibreStretcher()
        with pytest.raises(ValueError, match="outside the device range"):
            st.phase_for(200.0)
        with pytest.raises(ValueError, match="outside the device range"):
            st.phase_for(-1.0)

    def test_a_part_too_slow_for_a_fringe_raises(self):
        """v_pi above half the drive range cannot reach 2*pi at all, and
        that is a device limit worth hearing about rather than clipping."""
        st = PiezoFibreStretcher(v_pi=100.0, v_max=150.0)
        assert not st.delivers(TWO_PI - 0.01)
        with pytest.raises(ValueError, match="cannot reach a full fringe"):
            st.voltage_for(TWO_PI - 0.01)


class TestInsertionLoss:

    def test_applied_by_default(self):
        """A real device always has its loss. One that silently lost it
        would flatter every budget it appeared in."""
        st = PiezoFibreStretcher()
        assert st.apply_insertion_loss is True
        E = np.ones((4, 2), dtype=complex)
        out = st.apply(E, 10.0)
        assert np.sum(np.abs(out) ** 2) == pytest.approx(
            np.sum(np.abs(E) ** 2) * 10 ** (-0.1 / 10.0), rel=1e-12)

    def test_suppressible(self):
        """The suppression is a link-budget decision made at the call
        site, not a property of the part."""
        st = PiezoFibreStretcher(apply_insertion_loss=False)
        E = np.ones((4, 2), dtype=complex)
        assert np.sum(np.abs(st.apply(E, 10.0)) ** 2) == pytest.approx(
            np.sum(np.abs(E) ** 2), rel=1e-12)

    def test_the_phase_is_applied_either_way(self):
        for suppressed in (True, False):
            st = PiezoFibreStretcher(apply_insertion_loss=not suppressed)
            E = np.ones((1, 2), dtype=complex)
            got = np.angle(st.apply(E, st.v_pi)[0, 0])
            assert got == pytest.approx(math.pi) or \
                   got == pytest.approx(-math.pi)
