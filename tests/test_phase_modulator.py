"""Tests for src/channel/phase_modulator.py.

Verifies the V_π formula (Alferness 1988, Weis & Gaylord 1985), the
X-cut/Y-cut modulation-axis behaviour, DC/RF phase application, and
parameter validation.
"""
import numpy as np
import pytest
from src.channel.phase_modulator import PhaseModulator

LAMBDA = 1550e-9
D = 24.133e-6
N_O = 2.2
GAMMA = 0.8
L = 5.588e-2


def vpi_x():
    return LAMBDA * D / (2 * N_O ** 3 * 10.12e-12 * GAMMA * L)


def vpi_y():
    return LAMBDA * D / (2 * N_O ** 3 * 3.4e-12 * GAMMA * L)


@pytest.fixture
def pm_x():
    return PhaseModulator(crystal_cut='X', modulation='DC')


@pytest.fixture
def pm_y():
    return PhaseModulator(crystal_cut='Y', modulation='DC')


@pytest.fixture
def field():
    N = 10
    E = np.zeros((N, 2), dtype=complex)
    E[:, 0] = np.exp(1j * np.linspace(0, 0.4, N))
    E[:, 1] = np.exp(1j * np.linspace(0.2, 1.2, N))
    return E


class TestVpi:

    def test_x_cut_formula(self, pm_x):
        """X-cut V_π must match lambda*d / (2*n_o^3*r13*Gamma*L)."""
        assert np.isclose(pm_x.Vpi, vpi_x(), rtol=1e-9)

    def test_y_cut_formula(self, pm_y, pm_x):
        """Y-cut V_π must use r22 instead of r13."""
        assert np.isclose(pm_y.Vpi, vpi_y(), rtol=1e-9)
        assert pm_y.Vpi > pm_x.Vpi

    def test_vpi_scales_with_wavelength(self):
        pm = PhaseModulator(crystal_cut='X', params={'wavelength': 2 * LAMBDA})
        assert np.isclose(pm.Vpi, 2 * vpi_x(), rtol=1e-9)

    def test_vpi_frozen_after_construction(self):
        """V_π is computed once at construction; later param edits must
        not change it (cached in __Vpi)."""
        pm = PhaseModulator(crystal_cut='X')
        Vpi = pm.Vpi
        pm.L = 10 * L
        assert pm.Vpi == Vpi

    def test_zero_denominator_raises(self):
        with pytest.raises(ZeroDivisionError):
            PhaseModulator(crystal_cut='X', params={'Gamma': 0.0})

    def test_unknown_param_key_raises(self):
        with pytest.raises(RuntimeError):
            PhaseModulator(params={'Vpi': 5.0})


class TestGetPhi:

    def test_vpi_gives_pi(self, pm_x):
        assert np.isclose(pm_x.get_phi(pm_x.Vpi), np.pi, rtol=1e-12)

    def test_half_vpi_gives_half_pi(self, pm_x):
        assert np.isclose(pm_x.get_phi(pm_x.Vpi / 2), np.pi / 2, rtol=1e-12)

    def test_zero_voltage_gives_zero_phase(self, pm_x):
        assert pm_x.get_phi(0.0) == 0.0


class TestModulateDC:

    def test_x_cut_modulates_ey_only(self, pm_x, field):
        """X-cut must apply phase to Ey, leaving Ex untouched."""
        E_out = pm_x.modulate(field, pm_x.Vpi)  # phi = pi
        assert np.allclose(E_out[:, 0], field[:, 0])
        assert np.allclose(E_out[:, 1], -field[:, 1])

    def test_y_cut_modulates_ex_only(self, pm_y, field):
        """Y-cut must apply phase to Ex, leaving Ey untouched."""
        E_out = pm_y.modulate(field, pm_y.Vpi)  # phi = pi
        assert np.allclose(E_out[:, 0], -field[:, 0])
        assert np.allclose(E_out[:, 1], field[:, 1])

    def test_half_vpi_gives_j_phase(self, pm_x, field):
        """phi = pi/2 must multiply the modulated component by j."""
        E_out = pm_x.modulate(field, pm_x.Vpi / 2)
        assert np.allclose(E_out[:, 1], 1j * field[:, 1], atol=1e-12)

    def test_power_conserved(self, pm_x, field):
        """Phase modulation is unitary: power must be conserved."""
        E_out = pm_x.modulate(field, 2.0)
        assert np.isclose(np.sum(np.abs(E_out) ** 2),
                          np.sum(np.abs(field) ** 2), rtol=1e-12)

    def test_zero_voltage_is_identity(self, pm_x, field):
        E_out = pm_x.modulate(field, 0.0)
        assert np.allclose(E_out, field)

    def test_bad_field_shape_raises(self, pm_x):
        with pytest.raises(ValueError):
            pm_x.modulate(np.ones(4, dtype=complex), 1.0)


class TestModulateRF:

    @pytest.fixture
    def pm_rf(self):
        return PhaseModulator(crystal_cut='X', modulation='RF')

    def test_per_sample_phase(self, pm_rf, field):
        """RF mode must apply phase per time step, e.g. V = Vpi/4 -> e^{j pi/4}."""
        N = field.shape[0]
        V = np.full(N, pm_rf.Vpi / 4)
        E_out = pm_rf.modulate(field, V)
        assert np.allclose(E_out[:, 1],
                           field[:, 1] * np.exp(1j * np.pi / 4), atol=1e-12)
        assert np.allclose(E_out[:, 0], field[:, 0])

    def test_length_mismatch_raises(self, pm_rf, field):
        with pytest.raises(ValueError):
            pm_rf.modulate(field, np.array([1.0, 2.0]))


class TestValidation:

    def test_unknown_crystal_cut_raises(self):
        with pytest.raises(RuntimeError):
            PhaseModulator(crystal_cut='Z')

    def test_unknown_modulation_raises(self):
        with pytest.raises(RuntimeError):
            PhaseModulator(modulation='AM')

    def test_custom_params_override_defaults(self):
        pm = PhaseModulator(crystal_cut='X', params={'L': 2 * L})
        assert np.isclose(pm.Vpi, vpi_x() / 2, rtol=1e-9)


class TestBiasOffsetVoltage:
    """`bias_offset_v` -- the same static bias error as `phase_error_rad`,
    stated in the units a modulator is actually set in.

    Gobby et al. (2004) attribute their QBER floor first to "slight
    inaccuracies of the phase modulator biases".  Expressing that as a
    voltage puts it on the scale the hardware lives on rather than leaving
    it a bare angle: the conversion runs through this modulator's own
    crystal-derived V_pi, so nothing external is assumed.
    """

    def test_converts_through_vpi(self, pm_x):
        pm = PhaseModulator(crystal_cut='X', bias_offset_v=pm_x.Vpi / 4)
        assert np.isclose(pm.phase_error_rad, np.pi / 4, rtol=1e-12)

    def test_full_vpi_gives_pi(self, pm_x):
        pm = PhaseModulator(crystal_cut='X', bias_offset_v=pm_x.Vpi)
        assert np.isclose(pm.phase_error_rad, np.pi, rtol=1e-12)

    def test_round_trip_against_gobby_floor(self, pm_x):
        """20.93 deg <-> 451.5 mV, the offset reproducing a 3.3% floor."""
        deg = np.degrees(np.arccos(1.0 - 2.0 * 0.033))
        volts = np.radians(deg) * pm_x.Vpi / np.pi
        assert np.isclose(volts, 0.4515, atol=5e-4)
        pm = PhaseModulator(crystal_cut='X', bias_offset_v=volts)
        assert np.isclose(np.degrees(pm.phase_error_rad), deg, rtol=1e-9)

    def test_y_cut_uses_its_own_vpi(self, pm_y):
        """Conversion is per-device, not a shared constant."""
        pm = PhaseModulator(crystal_cut='Y', bias_offset_v=pm_y.Vpi / 2)
        assert np.isclose(pm.phase_error_rad, np.pi / 2, rtol=1e-12)

    def test_both_units_at_once_raises(self):
        """One mechanism, two unit systems: accepting both would sum them
        silently, which is how a bias error gets double-counted."""
        with pytest.raises(RuntimeError, match="same static bias error"):
            PhaseModulator(crystal_cut='X', bias_offset_v=0.45,
                           phase_error_rad=0.36)

    def test_default_is_inert(self, pm_x, field):
        """Zero bias must leave the component exactly as it was."""
        assert pm_x.bias_offset_v == 0.0
        assert pm_x.phase_error_rad == 0.0
        base = PhaseModulator(crystal_cut='X')
        assert np.array_equal(pm_x.modulate(field, 1.0),
                              base.modulate(field, 1.0))

    def test_bias_shifts_the_modulated_axis(self, pm_x, field):
        """The offset must actually reach the field, on the X-cut axis."""
        pm = PhaseModulator(crystal_cut='X', bias_offset_v=pm_x.Vpi / 2)
        out = pm.modulate(field, 0.0)
        assert np.allclose(out[:, 1], field[:, 1] * np.exp(1j * np.pi / 2),
                           atol=1e-12)
        assert np.allclose(out[:, 0], field[:, 0])

    def test_equivalent_to_radian_form(self, pm_x, field):
        """Both spellings are the same mechanism, so they must agree."""
        a = PhaseModulator(crystal_cut='X', bias_offset_v=pm_x.Vpi / 3)
        b = PhaseModulator(crystal_cut='X', phase_error_rad=np.pi / 3)
        assert np.allclose(a.modulate(field, 0.7), b.modulate(field, 0.7),
                           atol=1e-12)
