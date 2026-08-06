"""Tests for src/channel/optics.py.

Covers the PHYS-6 fix surface (circular_analyser vs pbs), unitarity of
the wave-plate / rotation transforms, projector idempotence of the
polarisers, VOA attenuation, and coupler power bookkeeping.
"""
import numpy as np
import pytest
from src.channel.optics import (
    circular_analyser, pbs, hadamard, polarizer, voa,
    halfwave, quarterwave, polarization_rotator,
    coupler_split, coupler_combine,
)


@pytest.fixture
def field():
    """Unit-power mixed-polarisation field (N, 2)."""
    N = 1000
    E = np.zeros((N, 2), dtype=complex)
    E[:, 0] = np.exp(1j * np.linspace(0, 0.3, N)) / np.sqrt(2)
    E[:, 1] = np.exp(1j * np.linspace(0.1, 1.0, N)) / np.sqrt(2)
    return E


def power(E):
    return np.sum(np.abs(E) ** 2, axis=1) if E.ndim == 2 else np.abs(E) ** 2


class TestHadamard:

    def test_unitary(self, field):
        """Hadamard must conserve power (unitary)."""
        E_out = hadamard(field)
        assert np.allclose(np.sum(power(E_out)), np.sum(power(field)),
                           rtol=1e-12)

    def test_splits_equally(self, field):
        """Hadamard of H input must split 50/50."""
        E = np.zeros((10, 2), dtype=complex)
        E[:, 0] = 1.0
        E_out = hadamard(E)
        assert np.allclose(np.sum(np.abs(E_out[:, 0]) ** 2), 5.0, rtol=1e-12)
        assert np.allclose(np.sum(np.abs(E_out[:, 1]) ** 2), 5.0, rtol=1e-12)


class TestPBS:

    def test_h_input_to_h_port_only(self):
        """Pure-H input must go entirely to the H port."""
        E = np.zeros((5, 2), dtype=complex)
        E[:, 0] = 1.0
        Ex, Ey = pbs(E)
        assert np.allclose(Ex, 1.0)
        assert np.allclose(Ey, 0.0)

    def test_v_input_to_v_port_only(self):
        """Pure-V input must go entirely to the V port."""
        E = np.zeros((5, 2), dtype=complex)
        E[:, 1] = 1.0
        Ex, Ey = pbs(E)
        assert np.allclose(Ex, 0.0)
        assert np.allclose(Ey, 1.0)

    def test_power_conserved(self, field):
        """PBS must conserve total power across the two ports."""
        Ex, Ey = pbs(field)
        P_out = np.sum(np.abs(Ex) ** 2) + np.sum(np.abs(Ey) ** 2)
        assert np.isclose(P_out, np.sum(power(field)), rtol=1e-12)

    def test_blind_to_relative_phase(self):
        """A true PBS is blind to the phase between Ex and Ey (PHYS-6)."""
        N = 4
        out_ports = []
        for phi in [0.0, np.pi / 2, np.pi]:
            E = np.zeros((N, 2), dtype=complex)
            E[:, 0] = 1.0 / np.sqrt(2)
            E[:, 1] = np.exp(1j * phi) / np.sqrt(2)
            Ex, Ey = pbs(E)
            out_ports.append((np.sum(np.abs(Ex) ** 2),
                              np.sum(np.abs(Ey) ** 2)))
        for i in range(1, len(out_ports)):
            assert np.allclose(out_ports[i], out_ports[0], atol=1e-12), \
                "PBS port powers must not depend on relative phase"


class TestCircularAnalyser:

    def test_unitary(self, field):
        """circular_analyser must conserve total power (unitary)."""
        Ex, Ey = circular_analyser(field)
        P_out = np.sum(np.abs(Ex) ** 2) + np.sum(np.abs(Ey) ** 2)
        assert np.isclose(P_out, np.sum(power(field)), rtol=1e-12)

    def test_h_input_splits_5050(self):
        """Pure-H input must split 50/50 (QWP+PBS behaviour)."""
        E = np.zeros((10, 2), dtype=complex)
        E[:, 0] = 1.0
        Ex, Ey = circular_analyser(E)
        assert np.isclose(np.sum(np.abs(Ex) ** 2), 5.0, rtol=1e-12)
        assert np.isclose(np.sum(np.abs(Ey) ** 2), 5.0, rtol=1e-12)

    def test_phase_discrimination(self):
        """circular_analyser must be sensitive to the Ex-Ey phase — the
        exact property that makes it an analyser and not a PBS (PHYS-6)."""
        N = 10
        # RCP (Ex = 1, Ey = -j): all power must land in the Ey port
        E_rcp = np.zeros((N, 2), dtype=complex)
        E_rcp[:, 0] = 1.0 / np.sqrt(2)
        E_rcp[:, 1] = -1j / np.sqrt(2)
        Ex, Ey = circular_analyser(E_rcp)
        assert np.sum(np.abs(Ex) ** 2) < 1e-12
        assert np.isclose(np.sum(np.abs(Ey) ** 2), N, rtol=1e-12)

    def test_matches_documented_ports(self, field):
        """Ex = (Ex_in - i*Ey_in)/sqrt(2), Ey = (-i*Ex_in + Ey_in)/sqrt(2)."""
        Ex, Ey = circular_analyser(field)
        assert np.allclose(Ex, (field[:, 0] - 1j * field[:, 1]) / np.sqrt(2),
                           rtol=1e-12)
        assert np.allclose(Ey, (-1j * field[:, 0] + field[:, 1]) / np.sqrt(2),
                           rtol=1e-12)

    def test_differs_from_pbs(self):
        """PHYS-6 regression: analyser output must differ from pbs for a
        phase-modulated 45-degree field."""
        N = 10
        phi = 0.7
        E = np.zeros((N, 2), dtype=complex)
        E[:, 0] = 1.0 / np.sqrt(2)
        E[:, 1] = np.exp(1j * phi) / np.sqrt(2)
        Ex_a, Ey_a = circular_analyser(E)
        Ex_p, Ey_p = pbs(E)
        assert not np.allclose(np.abs(Ex_a) ** 2, np.abs(Ex_p) ** 2)


class TestPolarizer:

    @pytest.mark.parametrize("pol", ['45', '-45'])
    def test_projector_idempotence(self, pol):
        """Applying the same projector twice must equal applying it once."""
        E = np.ones((10, 2), dtype=complex)
        E1 = polarizer(E, pol)
        E2 = polarizer(E1, pol)
        assert np.allclose(E1, E2, atol=1e-12)

    @pytest.mark.parametrize("pol", ['H', 'V'])
    def test_cascaded_polarizers_double_extinction(self, pol):
        """H/V polarisers pass 1 % of the orthogonal component; two
        cascaded passes must suppress it quadratically (0.01^2)."""
        E = np.zeros((10, 2), dtype=complex)
        E[:, 1 if pol == 'H' else 0] = 1.0   # pure orthogonal input
        E1 = polarizer(E, pol)
        E2 = polarizer(E1, pol)
        leak1 = power(E1).sum()
        leak2 = power(E2).sum()
        assert np.isclose(leak2, leak1 * 1e-4, rtol=1e-12)

    def test_cascaded_polarizers_pass_loss(self):
        """Dominant component scales by the 1 % pass-through loss per pass."""
        E = np.zeros((10, 2), dtype=complex)
        E[:, 0] = 1.0
        E1 = polarizer(E, 'H')
        E2 = polarizer(E1, 'H')
        assert np.allclose(E2[:, 0], 0.99 * E1[:, 0], atol=1e-12)

    def test_h_extinguishes_v(self):
        """H polariser must extinguish pure-V input."""
        E = np.zeros((10, 2), dtype=complex)
        E[:, 1] = 1.0
        E_out = polarizer(E, 'H')
        assert np.sum(power(E_out)) < 0.05 * 10  # ~1% leakage per component

    def test_45_passes_45(self):
        """45-degree projector passes a 45-degree input without loss."""
        E = np.ones((10, 2), dtype=complex)
        E_out = polarizer(E, '45')
        P_out = np.sum(power(E_out))
        P_in = np.sum(power(E))
        assert np.isclose(P_out, P_in, rtol=1e-12)

    def test_invalid_polarization_raises(self):
        with pytest.raises(Exception):
            polarizer(np.ones((2, 2), dtype=complex), 'X')


class TestVOA:

    def test_attenuation_scales_power(self, field):
        """10 dB VOA must divide power by 10."""
        E_out = voa(field, 10.0)
        assert np.isclose(np.sum(power(E_out)), np.sum(power(field)) / 10.0,
                          rtol=1e-12)

    def test_zero_dB_is_identity(self, field):
        E_out = voa(field, 0.0)
        assert np.allclose(E_out, field)

    def test_negative_dB_amplifies(self, field):
        E_out = voa(field, -3.0)
        assert np.sum(power(E_out)) > np.sum(power(field))


class TestWaveplates:

    @pytest.mark.parametrize("func", [halfwave, quarterwave,
                                      lambda E, **kw: polarization_rotator(E, 30)])
    def test_unitary(self, field, func):
        E_out = func(field)
        assert np.allclose(np.sum(power(E_out)), np.sum(power(field)),
                           rtol=1e-12)

    def test_halfwave_rotates_h_to_v_at_45(self):
        """HWP at 45 degrees rotates H polarisation to V."""
        E = np.zeros((10, 2), dtype=complex)
        E[:, 0] = 1.0
        E_out = halfwave(E, theta=45, rotation=True)
        P_x = np.sum(np.abs(E_out[:, 0]) ** 2)
        P_y = np.sum(np.abs(E_out[:, 1]) ** 2)
        assert P_x < 1e-12
        assert np.isclose(P_y, 10.0, rtol=1e-12)

    def test_quarterwave_h_to_circular(self):
        """QWP at 45 degrees converts H polarisation to circular (RCP)."""
        E = np.zeros((10, 2), dtype=complex)
        E[:, 0] = 1.0
        E_out = quarterwave(E, theta=45, rotation=True)
        P_x = np.sum(np.abs(E_out[:, 0]) ** 2)
        P_y = np.sum(np.abs(E_out[:, 1]) ** 2)
        assert np.isclose(P_x, 5.0, rtol=1e-12)
        assert np.isclose(P_y, 5.0, rtol=1e-12)
        assert np.allclose(E_out[:, 1], -1j * E_out[:, 0], atol=1e-12)

    def test_quarterwave_zero_angle_leaves_h(self):
        """QWP aligned with H leaves H unchanged up to a global phase."""
        E = np.zeros((10, 2), dtype=complex)
        E[:, 0] = 1.0
        E_out = quarterwave(E, rotation=False)
        assert np.allclose(E_out[:, 0], np.exp(-1j * np.pi / 4), atol=1e-12)
        assert np.allclose(E_out[:, 1], 0.0, atol=1e-12)

    def test_rotator_rotates_linear(self):
        """90-degree rotation must swap H and V."""
        E = np.zeros((10, 2), dtype=complex)
        E[:, 0] = 1.0
        E_out = polarization_rotator(E, 90)
        P_x = np.sum(np.abs(E_out[:, 0]) ** 2)
        P_y = np.sum(np.abs(E_out[:, 1]) ** 2)
        assert P_x < 1e-12
        assert np.isclose(P_y, 10.0, rtol=1e-12)


class TestCoupler:

    def test_split_power_conserved(self, field):
        """port_power + tap_power must equal input power."""
        P_in = np.sum(power(field))
        port_p, _, tap_p, _ = coupler_split(P_in, field, ratio=0.3)
        assert np.isclose(port_p + tap_p, P_in, rtol=1e-12)

    def test_split_ratio_applied(self, field):
        P_in = np.sum(power(field))
        port_p, _, tap_p, _ = coupler_split(P_in, field, ratio=0.3)
        assert np.isclose(port_p, P_in * 0.3, rtol=1e-12)
        assert np.isclose(tap_p, P_in * 0.7, rtol=1e-12)

    def test_split_invalid_ratio_raises(self, field):
        with pytest.raises(Exception):
            coupler_split(1.0, field, ratio=1.5)

    def test_combine_single_port_matrix(self, field):
        """Single-port combine applies one arm of the ideal 3 dB coupler:
        E_out = (E1 + j*E2)/sqrt(2), with power derived from the field."""
        E1 = field / np.sqrt(2)
        E2 = field / np.sqrt(2)
        P1 = np.sum(power(E1))
        P2 = np.sum(power(E2))
        pout, E_out = coupler_combine(P1, E1, P2, E2, out_ports=1)
        assert np.allclose(E_out, (E1 + 1j * E2) / np.sqrt(2), atol=1e-12)
        assert np.isclose(pout, np.sum(np.abs(E_out) ** 2), rtol=1e-12)

    def test_combine_two_ports_matrix(self, field):
        """Two-port combine applies the ideal 3 dB coupler scattering
        matrix: E_out1 = (E1 + j*E2)/sqrt(2), E_out2 = (j*E1 + E2)/sqrt(2)."""
        E1 = field / np.sqrt(2)
        E2 = field / np.sqrt(2)
        _, E_out1, _, E_out2 = coupler_combine(1.0, E1, 1.0, E2, out_ports=2)
        assert np.allclose(E_out1, (E1 + 1j * E2) / np.sqrt(2), atol=1e-12)
        assert np.allclose(E_out2, (1j * E1 + E2) / np.sqrt(2), atol=1e-12)

    def test_combine_two_ports_unitary(self, field):
        """Two-port combine must conserve power: the scattering matrix is
        unitary, so |E_out1|^2 + |E_out2|^2 = |E1|^2 + |E2|^2 sample-wise."""
        E1 = field / np.sqrt(2)
        E2 = field / np.sqrt(2)
        P1 = np.sum(power(E1))
        P2 = np.sum(power(E2))
        pout1, E_out1, pout2, E_out2 = coupler_combine(P1, E1, P2, E2,
                                                       out_ports=2)
        P_out = np.sum(np.abs(E_out1) ** 2) + np.sum(np.abs(E_out2) ** 2)
        assert np.isclose(P_out, P1 + P2, rtol=1e-12)
        assert np.isclose(pout1 + pout2, P1 + P2, rtol=1e-12)

    def test_combine_invalid_ports_raises(self, field):
        with pytest.raises(Exception):
            coupler_combine(1.0, field, 1.0, field, out_ports=3)
