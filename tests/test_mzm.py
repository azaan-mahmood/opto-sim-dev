import numpy as np
import pytest
from src.channel.mzm import MZM
from src.channel.phase_modulator import PhaseModulator


@pytest.fixture
def mzm():
    return MZM()


@pytest.fixture
def field():
    """Single field sample at 45° polarization, 1 mW."""
    E = np.array([1.0, 1.0], dtype=complex) * np.sqrt(1e-3)
    return E


@pytest.fixture
def field_ey():
    """Field with only Ey component (for X-cut modulation)."""
    return np.array([0.0, 1.0], dtype=complex) * np.sqrt(1e-3)


@pytest.fixture
def field_ex():
    """Field with only Ex component (for Y-cut modulation)."""
    return np.array([1.0, 0.0], dtype=complex) * np.sqrt(1e-3)


class TestMZM:

    def test_vpi_derived_from_phase_modulator(self, mzm):
        """V_pi should be > 0 and match a default PhaseModulator."""
        pm = PhaseModulator()
        assert np.isclose(mzm.V_pi, pm.Vpi)
        assert mzm.V_pi > 0

    def test_switching_voltage_alias(self, mzm):
        """switching_voltage should equal V_pi."""
        assert mzm.switching_voltage == mzm.V_pi

    def test_transmission_at_null(self, mzm, field_ey):
        """V = V_pi should fully extinguish the modulated component."""
        E_out = mzm.modulate(field_ey, mzm.V_pi)
        P_out = np.mean(np.abs(E_out)**2)
        assert P_out < 1e-10

    def test_transmission_at_peak(self, mzm, field):
        """V = 0 should give maximum output."""
        E_out = mzm.modulate(field, 0.0)
        P_out = np.mean(np.abs(E_out)**2)
        E_in = np.mean(np.abs(field)**2)
        assert P_out > 0.9 * E_in

    def test_quadrature_bias(self, mzm, field_ey):
        """V_bias = V_pi/2 should give 50% transmission on modulated axis."""
        mzm_bq = MZM(bias_voltage=mzm.V_pi / 2)
        E_out = mzm_bq.modulate(field_ey, 0.0)
        P_out = np.mean(np.abs(E_out)**2)
        E_in = np.mean(np.abs(field_ey)**2)
        assert np.isclose(P_out / E_in, 0.5, rtol=0.02)

    def test_transfer_function_symmetry(self, mzm, field):
        """Transfer should be symmetric: V and -V give same intensity."""
        E_pos = mzm.modulate(field, 1.0)
        E_neg = mzm.modulate(field, -1.0)
        P_pos = np.mean(np.abs(E_pos)**2)
        P_neg = np.mean(np.abs(E_neg)**2)
        assert np.isclose(P_pos, P_neg, rtol=1e-6)

    def test_push_pull_vs_single_drive(self, mzm, field_ey):
        """Single-drive should give different phase than push-pull."""
        pp = MZM(mode='push-pull')
        sd = MZM(mode='single-drive')
        V_test = pp.V_pi / 4
        E_pp = pp.modulate(field_ey, V_test)
        E_sd = sd.modulate(field_ey, V_test)
        assert not np.allclose(E_pp, E_sd)

    def test_mode_setter_valid(self, mzm):
        """mode should accept valid values."""
        mzm.mode = 'single-drive'
        assert mzm.mode == 'single-drive'
        mzm.mode = 'push-pull'
        assert mzm.mode == 'push-pull'

    def test_mode_setter_invalid(self, mzm):
        """mode should reject invalid values."""
        with pytest.raises(ValueError):
            mzm.mode = 'invalid'

    def test_vpi_xcut_vs_ycut(self):
        """Different cuts should give different V_pi."""
        x = MZM(pm=PhaseModulator(crystal_cut='X'))
        y = MZM(pm=PhaseModulator(crystal_cut='Y'))
        assert not np.isclose(x.V_pi, y.V_pi)

    def test_modulate_ndarray_input(self, mzm):
        """modulate should handle (N, 2) field with V array."""
        E_in = np.ones((100, 2), dtype=complex) * np.sqrt(1e-3)
        V = np.linspace(0, mzm.V_pi, 100)
        E_out = mzm.modulate(E_in, V)
        assert E_out.shape == (100, 2)
        P = np.mean(np.abs(E_out)**2, axis=1)
        assert P[0] > P[-1]

    def test_insertion_loss(self, field):
        """Insertion loss should reduce output power."""
        ideal = MZM(insertion_loss_db=None)
        lossy = MZM(insertion_loss_db=3)
        E_ideal = ideal.modulate(np.copy(field), 0.0)
        E_lossy = lossy.modulate(np.copy(field), 0.0)
        P_ideal = np.mean(np.abs(E_ideal)**2)
        P_lossy = np.mean(np.abs(E_lossy)**2)
        assert P_lossy < P_ideal

    def test_crystal_cut_applies_phase_to_correct_component(self, mzm, field):
        """X-cut should only modulate Ey."""
        E_out = mzm.modulate(field, mzm.V_pi)
        assert np.abs(E_out[1]) < 1e-10
