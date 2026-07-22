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


class TestCWLaserPulsed:

    PULSE_LASER_KW = dict(
        wavelength=1550e-9, power_dbm=0, linewidth=1e6,
        pulsed=True, pulse_width=100e-12, repetition_rate=2.5e6,
    )

    def test_pulsed_field_shape(self):
        """Pulsed sample_field returns (n_samples, 2)."""
        las = CWLaser(**self.PULSE_LASER_KW)
        E = las.sample_field(1e-12, 5000)
        assert E.shape == (5000, 2)
        assert np.iscomplexobj(E)

    def test_pulsed_power_conservation(self):
        """mean(|E|^2) in pulsed mode should equal _power_w."""
        las = CWLaser(**self.PULSE_LASER_KW)
        E = las.sample_field(1e-12, 5000)
        P = np.mean(np.sum(np.abs(E) ** 2, axis=1))
        assert np.isclose(P, las._power_w, rtol=0.05)

    def test_pulsed_pulse_width(self):
        """Gaussian FWHM should match pulse_width parameter."""
        las = CWLaser(**self.PULSE_LASER_KW)
        dt = 2e-12
        # Generate envelope with many periods to find an interior pulse
        n = int(3.0 / las.repetition_rate / dt)  # 3 periods
        env = las._pulse_envelope(dt, n)
        # Find the highest peak not at the array edges
        margin = int(las._pulse_sigma * 5 / dt)
        interior = env[margin:-margin]
        peak_idx = np.argmax(interior) + margin
        half = 0.5 * env[peak_idx]
        # Search left and right from peak
        left = np.where(env[:peak_idx] <= half)[0]
        right = np.where(env[peak_idx:] <= half)[0]
        if len(left) == 0 or len(right) == 0:
            pytest.skip("Cannot resolve half-max points at this resolution")
        fwhm_meas = ((right[0] + peak_idx) - left[-1]) * dt
        assert np.isclose(fwhm_meas, las.pulse_width, rtol=0.15)

    def test_pulsed_repetition_rate(self):
        """Pulse spacing should match repetition_rate."""
        las = CWLaser(**self.PULSE_LASER_KW)
        dt = 5e-12
        env = las._pulse_envelope(dt, 200000)
        # Find peaks
        threshold = 0.5 * np.max(env)
        peaks = []
        for i in range(1, len(env) - 1):
            if env[i] > threshold and env[i] > env[i - 1] and env[i] > env[i + 1]:
                peaks.append(i)
        if len(peaks) < 2:
            pytest.skip("Fewer than 2 peaks found")
        spacings = np.diff(peaks) * dt
        mean_spacing = np.mean(spacings)
        expected = 1.0 / las.repetition_rate
        assert np.isclose(mean_spacing, expected, rtol=0.02)

    def test_pulsed_zero_power_when_pulse_disabled(self):
        """pulsed=False should give same output as before (CW)."""
        cw = CWLaser(1550e-9, power_dbm=0, linewidth=1e6, pulsed=False)
        E_cw = cw.sample_field(1e-12, 2000)
        # All samples should have non-zero power in CW mode
        P_per_sample = np.sum(np.abs(E_cw) ** 2, axis=1)
        assert np.all(P_per_sample > 0)

    def test_pulsed_inter_pulse_zeros(self):
        """Samples far from any pulse centre should be near zero."""
        las = CWLaser(**self.PULSE_LASER_KW)
        E = las.sample_field(1e-12, 5000)
        P = np.sum(np.abs(E) ** 2, axis=1)
        # Mid-point within first period — within array bounds
        T = 1.0 / las.repetition_rate
        half_period = int(0.5 * T / 1e-12)
        if half_period >= len(P):
            half_period = len(P) // 2
        # Power at mid-point should be << peak power
        extinction = P[half_period] / np.max(P)
        assert extinction < 0.05  # at least 13 dB extinction

    def test_pulsed_energy_per_pulse(self):
        """Energy per pulse should equal P_avg / repetition_rate."""
        las = CWLaser(**self.PULSE_LASER_KW)
        dt = 1e-12
        T = 1.0 / las.repetition_rate
        n_period = int(T / dt)
        E = las.sample_field(dt, n_period * 20)
        P = np.sum(np.abs(E) ** 2, axis=1)
        energy_per_period = np.sum(P.reshape(-1, n_period), axis=1) * dt
        expected = las._power_w / las.repetition_rate
        assert np.isclose(np.mean(energy_per_period), expected, rtol=0.05)

    def test_pulsed_with_jitter(self):
        """Timing jitter should produce different pulse positions."""
        las = CWLaser(**self.PULSE_LASER_KW, timing_jitter_rms=20e-12)
        E1 = las.sample_field(1e-12, 20000)
        P1 = np.sum(np.abs(E1) ** 2, axis=1)
        E2 = las.sample_field(1e-12, 20000)
        P2 = np.sum(np.abs(E2) ** 2, axis=1)
        # Jitter should cause different peak positions
        assert not np.array_equal(np.argmax(P1), np.argmax(P2))
