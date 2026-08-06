import numpy as np
import pytest
from src.channel.fiber import propagate, apply_pmd, apply_birefringence


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
        E_out = propagate(50, np.copy(field))
        P_out = np.mean(np.abs(E_out)**2)
        assert P_out < P_in

    def test_attenuation_scales_with_distance(self, field):
        """Longer fibre should attenuate more."""
        E_10 = propagate(10, np.copy(field))
        E_50 = propagate(50, np.copy(field))
        P_10 = np.mean(np.abs(E_10)**2)
        P_50 = np.mean(np.abs(E_50)**2)
        assert P_50 < P_10

    def test_attenuation_formula_accuracy(self, field):
        """Attenuation should match exp(-alpha*L) within tolerance."""
        alpha = 0.182  # dB/km
        L = 50
        att_lin = 10 ** (-alpha * L / 10)
        P_in = np.mean(np.abs(field)**2)
        E_out = propagate(L, np.copy(field), attenuation_factor=alpha)
        P_out = np.mean(np.abs(E_out)**2)
        assert np.isclose(P_out / P_in, att_lin, rtol=0.01)

    def test_birefringence_preserves_power(self, field):
        """Birefringence should be unitary (power-conserving)."""
        P_in = np.mean(np.abs(field)**2)
        E_out = propagate(10, np.copy(field))
        P_out = np.mean(np.abs(E_out)**2)
        # power is slightly reduced by attenuation, but compare to
        # a pure-attenuation run
        E_att = propagate(10, np.copy(field), temperature=25, bend_radius=None)
        P_att = np.mean(np.abs(E_att)**2)
        assert np.isclose(P_out, P_att, rtol=0.01)

    def test_birefringence_changes_phase(self, field):
        """Birefringence should introduce a phase shift in Ex."""
        phase_in = np.angle(field[:, 0])
        E_out = propagate(10, np.copy(field), temperature=25, bend_radius=None)
        phase_out = np.angle(E_out[:, 0])
        avg_shift = np.mean(np.unwrap(phase_out - phase_in))
        assert np.abs(avg_shift) > 1e-6

    def test_birefringence_vs_temperature(self, field):
        """Different temperatures should give different phase shifts."""
        E_T1 = propagate(10, np.copy(field), temperature=0, bend_radius=None)
        E_T2 = propagate(10, np.copy(field), temperature=50, bend_radius=None)
        assert not np.allclose(E_T1, E_T2)

    def test_dispersion_requires_dt(self, field):
        """dispersion=True without dt should raise."""
        with pytest.raises(ValueError):
            propagate(10, np.copy(field), dispersion=True, dt=None)

    def test_dispersion_preserves_power(self, field):
        """CD should be unitary (no power loss)."""
        P_in = np.mean(np.abs(field)**2)
        E_out = propagate(10, np.copy(field), dt=1e-12, dispersion=True)
        P_out = np.mean(np.abs(E_out)**2)
        # power may change slightly from attenuation, but CD alone is lossless
        E_ref = propagate(10, np.copy(field))  # no dispersion
        P_ref = np.mean(np.abs(E_ref)**2)
        assert np.isclose(P_out, P_ref, rtol=0.02)

    def test_propagate_output_shape(self, field):
        """propagate should preserve input shape."""
        E_out = propagate(10, np.copy(field))
        assert E_out.shape == field.shape

    def test_propagate_zero_length(self, field):
        """Zero-length propagate should return field unchanged."""
        E_out = propagate(0, np.copy(field))
        assert np.allclose(E_out, field)

    def test_wavelength_parameter(self, field):
        """Different wavelength should affect birefringence phase."""
        E_1550 = propagate(10, np.copy(field), wavelength=1550e-9)
        E_1310 = propagate(10, np.copy(field), wavelength=1310e-9)
        assert not np.allclose(E_1550, E_1310)

    def test_regression_seeded_reproducibility(self):
        """propagate should be deterministic when np.random is seeded."""
        E = np.ones((100, 2), dtype=complex)
        np.random.seed(42)
        A = propagate(10, np.copy(E))
        np.random.seed(42)
        B = propagate(10, np.copy(E))
        assert np.allclose(A, B)


class TestPMD:
    """PHYS-1: pmd_coeff_ps_sqrt_km is in ps/sqrt(km) (Corning SMF-28
    Ultra spec <= 0.1 ps/sqrt(km) [12] in fiber.py's literature list).
    The DGD is Maxwell-distributed with scale a = sigma/sqrt(3), where
    sigma = pmd_coeff_ps_sqrt_km * sqrt(L_km) (Razavi [5] Fig 2.11).
    A Maxwell distribution's mean is 2*a*sqrt(2/pi).
    """

    def test_mean_dgd_matches_maxwellian_analytic_formula(self):
        """0.1 ps/sqrt(km) at 100 km: sigma = 1.0 ps, mean DGD = 0.921 ps."""
        E = np.ones((10, 2), dtype=complex)
        dt = 1e-12
        L_m = 100_000  # 100 km
        pmd_coeff = 0.1  # ps/sqrt(km)

        sigma_ps = pmd_coeff * np.sqrt(L_m / 1000.0)
        a_ps = sigma_ps / np.sqrt(3)
        expected_mean_ps = 2 * a_ps * np.sqrt(2 / np.pi)
        assert np.isclose(expected_mean_ps, 0.921, atol=0.001)

        np.random.seed(7)
        n = 20_000
        dgds_ps = np.array([
            apply_pmd(E.copy(), dt, L_m, pmd_coeff_ps_sqrt_km=pmd_coeff)[1] * 1e12
            for _ in range(n)
        ])
        assert np.isclose(dgds_ps.mean(), expected_mean_ps, rtol=0.03)


class TestBirefringenceDepolarization:
    """PHYS-4: for a product of N i.i.d. random-axis SU(2) rotations of
    fixed angle alpha, the ensemble-mean Stokes vector decays exactly as

        <S_out> = p^N * S_in,   p = (1 + 2*cos(alpha)) / 3

    (Menyuk & Wai [10], JOSA B 11(7), 1994 — random birefringence axis
    model; see PHYS-4 in opto-sim-issues-and-fixes.md, which verifies this
    against 400 realizations per point at Delta_n=5e-8, Delta_z=1m,
    alpha=0.2027 rad).
    """

    def test_ensemble_mean_stokes_matches_pN_law(self):
        wavelength = 1550e-9
        delta_n = 5.0e-8   # Agrawal [6] Sec 4.1 nominal birefringence @ 25C
        dz = 1.0           # metres — correlation_length, small-angle regime
        alpha = 2 * np.pi * delta_n * dz / wavelength
        p = (1 + 2 * np.cos(alpha)) / 3

        E_in = np.array([[1.0, 0.0]], dtype=complex)  # S1_in = 1
        n_realizations = 400
        lengths_m = [10, 50, 100, 200]

        np.random.seed(0)
        for L_m in lengths_m:
            N = round(L_m / dz)
            S1_vals = []
            for _ in range(n_realizations):
                out = apply_birefringence(E_in.copy(), L_m, wavelength=wavelength,
                                          correlation_length=dz, model='sectional')
                S1_vals.append(abs(out[0, 0]) ** 2 - abs(out[0, 1]) ** 2)
            mean_S1 = np.mean(S1_vals)
            predicted = p ** N
            # Finite-sample floor from 400 realizations, per the doc's own
            # measurement (~1/sqrt(400) ~= 0.05).
            tol = 0.08
            assert abs(mean_S1 - predicted) < tol, (
                f"L={L_m} m, N={N}: measured mean(S1)={mean_S1:.4f}, "
                f"predicted p^N={predicted:.4f}"
            )


class TestOrderedProductTreeReduction:
    """PERF-1: _ordered_product replaces an O(N) Python loop with an
    O(log N) vectorised pairwise tree reduction. Matrix multiplication is
    associative, so this must be exact (not an approximation) against the
    naive left-fold loop, for both even and odd N.
    """

    @staticmethod
    def _naive_ordered_product(J):
        """J[N-1] @ ... @ J[1] @ J[0], computed the slow obvious way."""
        J_total = np.eye(2, dtype=complex)
        for k in range(J.shape[0]):
            J_total = J[k] @ J_total
        return J_total

    @pytest.mark.parametrize("N", [1, 2, 3, 4, 7, 8, 15, 100, 101, 1001])
    def test_matches_naive_loop_exactly(self, N):
        from src.channel.fiber import _ordered_product

        rng = np.random.default_rng(N)  # vary the seed with N for coverage
        # Random unitary-ish 2x2 complex matrices (need not be physical
        # Jones matrices — associativity is a property of matrix
        # multiplication in general, not specific to SU(2)).
        J = rng.normal(size=(N, 2, 2)) + 1j * rng.normal(size=(N, 2, 2))

        expected = self._naive_ordered_product(J)
        actual = _ordered_product(J)
        assert np.allclose(actual, expected, atol=1e-10)
