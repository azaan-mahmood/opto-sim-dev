import numpy as np
import pytest
from src.lasers.cwlaser import CWLaser


@pytest.fixture
def laser():
    return CWLaser(
        wavelength=1550e-9,
        power_dbm=0.0,
        linewidth=1e6,
        rin_density=-150.0,
        polarization_azimuth=np.pi / 4,
        polarization_ellipticity=0.0,
    )


class TestCWLaser:

    def test_power_convention(self, laser):
        """mean(|E|^2) from sample_field should equal _power_w."""
        E = laser.sample_field(1e-12, 5000)
        P_meas = np.mean(np.sum(np.abs(E)**2, axis=1))
        assert np.isclose(P_meas, laser._power_w, rtol=0.05)

    def test_sample_field_shape(self, laser):
        """sample_field returns (n_samples, 2) complex array."""
        E = laser.sample_field(1e-12, 100)
        assert E.shape == (100, 2)
        assert np.iscomplexobj(E)

    def test_sample_field_power_scales_with_power_dbm(self):
        """Higher power_dbm should give higher mean|E|^2."""
        lo = CWLaser(1550e-9, power_dbm=-10, linewidth=1e6)
        hi = CWLaser(1550e-9, power_dbm=10, linewidth=1e6)
        E_lo = lo.sample_field(1e-12, 1000)
        E_hi = hi.sample_field(1e-12, 1000)
        P_lo = np.mean(np.sum(np.abs(E_lo)**2, axis=1))
        P_hi = np.mean(np.sum(np.abs(E_hi)**2, axis=1))
        assert P_hi > P_lo

    def test_phase_noise_increases_with_linewidth(self):
        """Phase variance should increase with linewidth."""
        nw = CWLaser(1550e-9, power_dbm=0, linewidth=1e6)
        ww = CWLaser(1550e-9, power_dbm=0, linewidth=100e6)
        dt, N = 1e-12, 20000
        phi_nw = np.unwrap(np.angle(nw.sample_field(dt, N)[:, 0]))
        phi_ww = np.unwrap(np.angle(ww.sample_field(dt, N)[:, 0]))
        var_nw = np.var(phi_nw[-5000:])
        var_ww = np.var(phi_ww[-5000:])
        assert var_ww > var_nw

    def test_polarization_vector_normalized(self, laser):
        """Jones vector should have unit norm."""
        pol = laser._polarization_vector()
        assert np.isclose(np.linalg.norm(pol), 1.0)

    def test_polarization_linear(self):
        """Linear polarization (chi=0) should give real-only Jones."""
        laser = CWLaser(1550e-9, power_dbm=0, linewidth=1e6,
                        polarization_azimuth=0, polarization_ellipticity=0)
        pol = laser._polarization_vector()
        assert np.allclose(pol.imag, 0, atol=1e-15)

    def test_regression_seeded_reproducibility(self):
        """Same seed should produce identical sample_field."""
        np.random.seed(12345)
        a = CWLaser(1550e-9, power_dbm=0, linewidth=1e6).sample_field(1e-12, 50)
        np.random.seed(12345)
        b = CWLaser(1550e-9, power_dbm=0, linewidth=1e6).sample_field(1e-12, 50)
        assert np.allclose(a, b)

    def test_instantaneous_field_over_period_shape(self, laser):
        """instantaneous_field(over_period=True) returns (n_samples, 2)."""
        E = laser.instantaneous_field(over_period=True, n_samples=100)
        assert E.shape == (100, 2)

    def test_instantaneous_field_single_sample(self, laser):
        """instantaneous_field(over_period=False) returns (2,)."""
        E = laser.instantaneous_field(over_period=False)
        assert E.shape == (2,)

    def test_power_out_matches_mw(self, laser):
        """power_out returns power in mW."""
        assert np.isclose(laser.power_out, 1.0)

    def test_narrow_linewidth_phase_coeff(self):
        """Zero linewidth should give zero phase diffusion."""
        laser = CWLaser(1550e-9, power_dbm=0, linewidth=0)
        assert laser._phase_diff_coeff == 0.0
        phi = laser._sample_phase_noise(1e-12, 1000)
        assert np.all(phi == 0)

    def test_rin_zero_for_highly_negative_density(self):
        """RIN should be much smaller at very negative rin_density."""
        low = CWLaser(1550e-9, power_dbm=0, linewidth=1e6, rin_density=-200)
        high = CWLaser(1550e-9, power_dbm=0, linewidth=1e6, rin_density=-130)
        rin_low = low._sample_rin(1e-12, 5000)
        rin_high = high._sample_rin(1e-12, 5000)
        assert np.std(rin_low) < np.std(rin_high)
