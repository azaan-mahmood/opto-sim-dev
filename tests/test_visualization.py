"""Tests for `src/visualization/`.

Covers the eye diagram, which is the only module there carrying physics
rather than drawing: it runs a real `MZM` and a real source, so it can be
wrong in ways a plotting helper cannot.
"""
import matplotlib
matplotlib.use('Agg')

import numpy as np
import pytest

from src.lasers.cwlaser import CWLaser
from src.lasers.dfblaser import DFBLaser
from src.lasers.laser_driver import DriveParams, LaserDriver
from src.visualization.eye import _align, eye_diagram


def _cw(azimuth=np.pi / 4, seed=42):
    """A CW source with its noise pinned.

    `CWLaser` draws from the global `numpy.random` state and takes no seed
    of its own, so pinning it is the caller's job — the same thing
    `validate_cwlaser.py` does at import.
    """
    np.random.seed(seed)
    return CWLaser(1550e-9, power_dbm=0, linewidth=1e6, rin_density=-130,
                   polarization_azimuth=azimuth)


def _dfb(seed=11):
    return LaserDriver(DFBLaser(n_sections=15, seed=seed),
                       DriveParams(mode='cw', i_bias=0.120), seed=seed)


class TestAlign:
    """The polarisation controller in front of the modulator."""

    def test_conserves_power(self):
        """It is a unitary, not a projection. A projection would throw away
        the component it cannot use and silently halve the power."""
        E = _cw().sample_field(1e-12, 200)
        for axis in (0, 1):
            assert np.sum(np.abs(_align(E, axis)) ** 2) == pytest.approx(
                np.sum(np.abs(E) ** 2), rel=1e-12)

    def test_puts_all_light_on_the_requested_axis(self):
        E = _cw().sample_field(1e-12, 200)
        total = np.sum(np.abs(E) ** 2)
        for axis in (0, 1):
            other = 1 - axis
            assert np.sum(np.abs(_align(E, axis)[:, other]) ** 2) / total < 1e-20

    def test_zero_field_is_left_alone(self):
        """No polarisation to rotate, and normalising would divide by zero."""
        E = np.zeros((10, 2), dtype=complex)
        np.testing.assert_array_equal(_align(E, 1), E)


class TestEyeDiagram:

    def test_repeats_with_the_source_pinned(self):
        """A figure that differs every run cannot be compared against a
        previous version of itself, which is the point of the seed."""
        a = eye_diagram(_cw(), 10e9, n_bits=64, spb=32, seed=7)
        b = eye_diagram(_cw(), 10e9, n_bits=64, spb=32, seed=7)
        np.testing.assert_array_equal(a, b)

    def test_a_different_pattern_gives_a_different_eye(self):
        a = eye_diagram(_cw(), 10e9, n_bits=64, spb=32, seed=7)
        c = eye_diagram(_cw(), 10e9, n_bits=64, spb=32, seed=8)
        assert not np.array_equal(a, c)

    def test_serves_both_sources(self):
        """The contract is `sample_field(dt, n)` and nothing more, so the
        CW source and the gain-switched DFB both qualify."""
        for laser, rate in ((_cw(), 10e9), (_dfb(), 2e9)):
            eye = eye_diagram(laser, rate, n_bits=32, spb=32, seed=7)
            assert eye.shape[1] == 64
            assert eye.min() >= 0.0

    @pytest.mark.parametrize('laser_fn,rate', [(_cw, 10e9), (_dfb, 2e9)])
    def test_alignment_is_what_makes_the_modulator_bite(self, laser_fn, rate):
        """Regression for the reason `align` exists.

        The MZM is X-cut and modulates Ey. `LaserDriver` emits entirely on
        Ex, so unaligned its eye is a flat line with the modulator doing
        nothing; `CWLaser` at 45 degrees leaves half the power unmodulated
        as a floor. Neither is a fact about the modulator, which is what
        the figure would appear to be reporting.
        """
        def er(eye):
            return eye.max() / max(eye.min(), 1e-30)

        unaligned = eye_diagram(laser_fn(), rate, n_bits=32, spb=32, seed=7,
                                align=False)
        aligned = eye_diagram(laser_fn(), rate, n_bits=32, spb=32, seed=7,
                              align=True)
        assert er(unaligned) < 10.0           # under 10 dB, and often ~1
        assert er(aligned) > 1e6              # ideal MZM, so effectively total

    def test_extinction_ratio_sets_the_off_rail(self):
        """The floor is the modulator's, and it has to be reachable —
        otherwise the off rail is exactly zero and a reader cannot tell a
        device limit from a source-only picture."""
        ideal = eye_diagram(_cw(), 10e9, n_bits=64, spb=32, seed=7)
        real = eye_diagram(_cw(), 10e9, n_bits=64, spb=32, seed=7,
                           extinction_ratio_db=20.0)
        assert ideal.min() < 1e-12
        assert real.min() > ideal.min()
        # 20 dB nominal, measured against the peak rather than the mean rail
        assert 10 * np.log10(real.max() / real.min()) == pytest.approx(
            20.0, abs=3.0)

    def test_band_limiting_is_what_opens_the_eye(self):
        """With an unfiltered drive the edges are one sample wide, so there
        is no opening to look at — which is what the ported version did.

        Measured as the share of samples caught in transit near half
        height. A step edge has none there; a band-limited one spends
        several samples crossing, and a narrower drive spends more.
        """
        def in_transit(eye):
            lo, hi = np.percentile(eye, 5), np.percentile(eye, 95)
            return np.mean(np.abs(eye - 0.5 * (lo + hi)) < 0.10 * (hi - lo))

        step = eye_diagram(_cw(), 10e9, n_bits=64, spb=32, seed=7,
                           drive_bandwidth=0)
        band = eye_diagram(_cw(), 10e9, n_bits=64, spb=32, seed=7)
        narrow = eye_diagram(_cw(), 10e9, n_bits=64, spb=32, seed=7,
                             drive_bandwidth=4e9)
        assert in_transit(step) < 0.005          # measured 0.0000
        assert in_transit(band) > 0.015          # measured 0.0292
        assert in_transit(narrow) > in_transit(band)   # 0.0544 > 0.0292

    def test_too_short_a_pattern_raises(self):
        """Returning quietly would leave empty axes looking like a result."""
        with pytest.raises(ValueError, match="at least 2"):
            eye_diagram(_cw(), 10e9, n_bits=1, spb=32, seed=7)
