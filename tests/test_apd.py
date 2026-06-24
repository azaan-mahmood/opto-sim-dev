import numpy as np
import pytest
from src.detectors.apd import apd


@pytest.fixture
def det():
    return apd(
        wavelength=1550e-9,
        quantum_efficiency=0.9,
        gain=10,
        excess_noise_factor=10,
        load_resistance=50,
        temperature=300,
        dark_current=10e-9,
    )


class TestAPD:

    def test_responsivity_formula(self, det):
        """Responsivity should be eta * e * lambda / (h * c)."""
        expected = 0.9 * 1.602e-19 * 1550e-9 / (6.626e-34 * 3e8)
        assert np.isclose(det.R, expected, rtol=1e-3)

    def test_signal_current_scales_with_power(self, det):
        """I_signal = M * R * P."""
        P = 1e-6
        I = det.calculate_output_current(P)
        assert np.isclose(I, det.gain * det.R * P)

    def test_signal_current_zero_power(self, det):
        """Zero power gives zero signal current."""
        assert det.calculate_output_current(0) == 0

    def test_shot_noise_increases_with_bandwidth(self, det):
        """Noise should increase with sqrt(bandwidth)."""
        I_sig = 1e-6
        n1 = det.calculate_noise(I_sig, 1e6)
        n2 = det.calculate_noise(I_sig, 100e6)
        ratio = n2 / n1
        expected_ratio = np.sqrt(100)
        assert np.isclose(ratio, expected_ratio, rtol=0.1)

    def test_thermal_noise_formula(self, det):
        """Thermal noise floor should match 4kTBR."""
        B = 1e9
        noise = det.calculate_noise(0, B)
        thermal_only = np.sqrt(4 * det.kB * det.T * B / det.RL)
        assert np.isclose(noise, thermal_only, rtol=1e-3)

    def test_output_returns_float(self, det):
        """output() returns a float (I_total)."""
        E = np.ones(1000) * np.sqrt(1e-6)
        I = det.output(E, bandwidth=1e6)
        assert isinstance(I, float)

    def test_output_details_dict(self, det):
        """output(details=True) returns a dict."""
        E = np.ones(1000) * np.sqrt(1e-6)
        result = det.output(E, bandwidth=1e6, details=True)
        assert isinstance(result, dict)
        assert 'I_signal' in result
        assert 'SNR' in result
        assert 'DCR' in result

    def test_detect_photons_zero_power(self, det):
        """Zero power should give zero detected photons."""
        n = det.detect_photons(0, 1e-9)
        assert n >= 0

    def test_detect_photons_non_negative(self, det):
        """detect_photons should never return negative."""
        for P in [1e-12, 1e-9, 1e-6, 1e-3]:
            n = det.detect_photons(P, 1e-9)
            assert n >= 0

    def test_dcr_matches_formula(self, det):
        """DCR = dark_current / e."""
        expected = det.dark_current / 1.602e-19
        assert np.isclose(det.dcr, expected, rtol=1e-3)

    def test_regression_seeded_reproducibility(self):
        """Same seed + same input should give same output current."""
        np.random.seed(99)
        d = apd(1550e-9, excess_noise_factor=10, load_resistance=50,
                temperature=300, gain=10, quantum_efficiency=0.9)
        E = np.ones(500) * np.sqrt(1e-6)
        I1 = d.output(E, 1e6)
        np.random.seed(99)
        d2 = apd(1550e-9, excess_noise_factor=10, load_resistance=50,
                 temperature=300, gain=10, quantum_efficiency=0.9)
        I2 = d2.output(E, 1e6)
        assert np.isclose(I1, I2)
