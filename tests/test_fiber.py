import numpy as np
import pytest
from src.channel.fiber_sectional import cable


@pytest.fixture
def field():
    """10 samples of Ex-only field at ~1 mW total power."""
    np.random.seed(42)
    E = np.zeros((1000, 2), dtype=complex)
    E[:, 0] = np.sqrt(0.5) * np.random.randn(1000)
    return E


class TestFiber:

    def test_attenuation_reduces_power(self, field):
        """Cable should reduce optical power per dB/km."""
        P_in = np.mean(np.abs(field)**2)
        E_out = cable(50, np.copy(field))
        P_out = np.mean(np.abs(E_out)**2)
        assert P_out < P_in

    def test_attenuation_scales_with_distance(self, field):
        """Longer fibre should attenuate more."""
        E_10 = cable(10, np.copy(field))
        E_50 = cable(50, np.copy(field))
        P_10 = np.mean(np.abs(E_10)**2)
        P_50 = np.mean(np.abs(E_50)**2)
        assert P_50 < P_10

    def test_attenuation_formula_accuracy(self, field):
        """Attenuation should match exp(-alpha*L) within tolerance."""
        alpha = 0.182  # dB/km
        L = 50
        att_lin = 10 ** (-alpha * L / 10)
        P_in = np.mean(np.abs(field)**2)
        E_out = cable(L, np.copy(field), attenuation_factor=alpha)
        P_out = np.mean(np.abs(E_out)**2)
        assert np.isclose(P_out / P_in, att_lin, rtol=0.01)

    def test_birefringence_preserves_power(self, field):
        """Birefringence should be unitary (power-conserving)."""
        P_in = np.mean(np.abs(field)**2)
        E_out = cable(10, np.copy(field))
        P_out = np.mean(np.abs(E_out)**2)
        # power is slightly reduced by attenuation, but compare to
        # a pure-attenuation run
        E_att = cable(10, np.copy(field), temperature=25, bend_radius=None)
        P_att = np.mean(np.abs(E_att)**2)
        assert np.isclose(P_out, P_att, rtol=0.01)

    def test_birefringence_changes_phase(self, field):
        """Birefringence should introduce a phase shift in Ex."""
        phase_in = np.angle(field[:, 0])
        E_out = cable(10, np.copy(field), temperature=25, bend_radius=None)
        phase_out = np.angle(E_out[:, 0])
        avg_shift = np.mean(np.unwrap(phase_out - phase_in))
        assert np.abs(avg_shift) > 1e-6

    def test_birefringence_vs_temperature(self, field):
        """Different temperatures should give different phase shifts."""
        E_T1 = cable(10, np.copy(field), temperature=0, bend_radius=None)
        E_T2 = cable(10, np.copy(field), temperature=50, bend_radius=None)
        assert not np.allclose(E_T1, E_T2)

    def test_dispersion_requires_dt(self, field):
        """dispersion=True without dt should raise."""
        with pytest.raises(ValueError):
            cable(10, np.copy(field), dispersion=True, dt=None)

    def test_dispersion_preserves_power(self, field):
        """CD should be unitary (no power loss)."""
        P_in = np.mean(np.abs(field)**2)
        E_out = cable(10, np.copy(field), dt=1e-12, dispersion=True)
        P_out = np.mean(np.abs(E_out)**2)
        # power may change slightly from attenuation, but CD alone is lossless
        E_ref = cable(10, np.copy(field))  # no dispersion
        P_ref = np.mean(np.abs(E_ref)**2)
        assert np.isclose(P_out, P_ref, rtol=0.02)

    def test_cable_output_shape(self, field):
        """cable should preserve input shape."""
        E_out = cable(10, np.copy(field))
        assert E_out.shape == field.shape

    def test_cable_zero_length(self, field):
        """Zero-length cable should return field unchanged."""
        E_out = cable(0, np.copy(field))
        assert np.allclose(E_out, field)

    def test_wavelength_parameter(self, field):
        """Different wavelength should affect birefringence phase."""
        E_1550 = cable(10, np.copy(field), wavelength=1550e-9)
        E_1310 = cable(10, np.copy(field), wavelength=1310e-9)
        assert not np.allclose(E_1550, E_1310)

    def test_regression_seeded_reproducibility(self):
        """cable should be deterministic when np.random is seeded."""
        E = np.ones((100, 2), dtype=complex)
        np.random.seed(42)
        A = cable(10, np.copy(E))
        np.random.seed(42)
        B = cable(10, np.copy(E))
        assert np.allclose(A, B)
