import numpy as np
import pytest
from src.channel.fiber import (propagate, apply_pmd, apply_birefringence,
                               bend_birefringence, FiberRealization)


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


class TestUlrichBendLaw:
    """REPRO-4 fix #1: the bend-induced birefringence must match Ulrich's
    published law Delta_n_bend = 0.135 * (r_fiber / R)^2 (Ulrich,
    Rashleigh & Eickhoff, Opt. Lett. 5(6), 1980, Eq. 1 — ref [7] in
    fiber.py). Tested both as a pure formula unit and end-to-end through
    the sectional Jones matrix, where the bend term must appear in the
    section retardance.
    """

    R_FIBER = 62.5e-6   # SMF-28 cladding radius used across the codebase

    def test_bend_birefringence_matches_published_law(self):
        """Exact match against 0.135*(r/R)^2 across a bend-radius sweep
        (stated tolerance: rtol = 1e-12)."""
        radii = np.array([0.002, 0.003, 0.005, 0.01, 0.02, 0.05, 0.1])
        expected = 0.135 * (self.R_FIBER / radii) ** 2
        measured = np.array([bend_birefringence(R) for R in radii])
        assert np.allclose(measured, expected, rtol=1e-12)

    def test_bend_birefringence_reference_case(self):
        """2 mm bend on SMF-28: Delta_n_bend ~ 1.32e-4 (the canonical
        tight-bend value used in the validation figure)."""
        assert np.isclose(bend_birefringence(2e-3), 0.135 * (62.5e-6 / 2e-3) ** 2,
                          rtol=1e-12)

    def test_bend_birefringence_custom_fiber_radius(self):
        """r_fiber is a parameter, not a hardcoded constant."""
        assert np.isclose(bend_birefringence(5e-3, r_fiber=125e-6),
                          0.135 * (125e-6 / 5e-3) ** 2, rtol=1e-12)

    @staticmethod
    def _section_jones(bend_radius, length_m=50.0, seed=1234):
        """Full 2x2 Jones matrix of a single-section fibre (L < L_c, so
        N = 1), via the public API only (FiberRealization.apply with the
        two orthogonal basis inputs)."""
        fr = FiberRealization(length_m, wavelength=1550e-9, temperature=25,
                              bend_radius=bend_radius,
                              model='sectional', attenuation=False, seed=seed)
        E_h = np.array([[1.0, 0.0]], dtype=complex)
        E_v = np.array([[0.0, 1.0]], dtype=complex)
        U = np.stack([fr.apply(E_h)[0], fr.apply(E_v)[0]], axis=1)
        return U

    @pytest.mark.parametrize("R_mm, rtol", [
        (6, 0.05), (8, 0.05), (10, 0.05), (15, 0.08), (30, 0.12), (60, 0.18),
    ])
    def test_bend_law_flows_through_sectional_jones(self, R_mm, rtol):
        """Recover Delta_n_bend from the single-section retardance
        delta = pi * |Delta_n| * L / lambda (eigenphase of the SU(2)
        matrix, axis-independent) and compare against Ulrich's law.
        The fibre length is scaled to each radius so the retardance is
        exactly pi/2 — an SU(2) matrix only carries retardance mod 2*pi,
        so an unwrapped comparison needs delta < pi. Tolerance absorbs
        the model's 10 % stochastic residual on the base birefringence."""
        wavelength = 1550e-9
        base_dn = 5.0e-8     # Agrawal [6] Sec 4.1 nominal, 25 C (no temp term)

        dn_bend_expected = bend_birefringence(R_mm * 1e-3)
        L_m = 0.5 * wavelength / (base_dn + dn_bend_expected)  # half = pi/2
        assert L_m < 50.0, "length must stay within one correlation cell"

        U = self._section_jones(R_mm * 1e-3, length_m=L_m)
        delta = np.abs(np.angle(np.linalg.eigvals(U)[0]))
        assert 0 < delta <= np.pi / 2 + 1e-9, "retardance must be unwrapped"

        dn_measured = delta * wavelength / (np.pi * L_m)
        dn_bend_measured = dn_measured - base_dn
        assert np.isclose(dn_bend_measured, dn_bend_expected, rtol=rtol), (
            f"R = {R_mm} mm: measured Delta_n_bend = {dn_bend_measured:.3e}, "
            f"Ulrich law = {dn_bend_expected:.3e}"
        )

    def test_bend_radius_changes_jones_matrix(self):
        """No-bend vs 6 mm-bend single-section fibres must differ."""
        U_flat = self._section_jones(None)
        U_bent = self._section_jones(6e-3)
        assert not np.allclose(U_flat, U_bent, atol=1e-6)


class TestBirefringenceSelfConsistency:
    """REPRO-4 fix #2: the 13 self-consistency checks formerly living in
    analysis/validation/validate_birefringence.py moved into the test
    suite, where they belong. They are internal invariants of the two
    birefringence models (not literature comparisons — those live in
    TestBirefringenceDepolarization and TestUlrichBendLaw)."""

    WAVELENGTH = 1550e-9

    @staticmethod
    def _field(n=1000):
        return np.random.randn(n, 2) + 1j * np.random.randn(n, 2)

    def test_power_conservation_sectional(self):
        E = self._field()
        P_in = np.mean(np.abs(E) ** 2)
        for L_m in [1, 10, 100, 1000]:
            E_out = apply_birefringence(E.copy(), L_m, wavelength=self.WAVELENGTH,
                                        model='sectional')
            P_out = np.mean(np.abs(E_out) ** 2)
            assert abs(P_out - P_in) / P_in < 1e-12

    def test_temperature_dependence_sectional(self):
        E = np.ones((100, 2), dtype=complex)
        np.random.seed(43)
        Js = [apply_birefringence(E.copy(), 1000, wavelength=self.WAVELENGTH,
                                  temperature=T, model='sectional')[0, 0]
              for T in [0, 25, 50]]
        assert not np.allclose(Js[0], Js[1])

    def test_wavelength_dependence_sectional(self):
        E = np.ones((100, 2), dtype=complex)
        np.random.seed(44)
        Js = [apply_birefringence(E.copy(), 1000, wavelength=lam,
                                  model='sectional')[0, 0]
              for lam in [1310e-9, 1550e-9]]
        assert not np.allclose(Js[0], Js[1])

    def test_seed_dependence_sectional(self):
        E = np.ones((100, 2), dtype=complex)
        np.random.seed(42)
        J1 = apply_birefringence(E.copy(), 100, wavelength=self.WAVELENGTH,
                                 model='sectional')[0, 0]
        np.random.seed(142)
        J2 = apply_birefringence(E.copy(), 100, wavelength=self.WAVELENGTH,
                                 model='sectional')[0, 0]
        assert not np.allclose(J1, J2)

    def test_output_variation_sectional(self):
        E = np.array([[1.0, 0.0]], dtype=complex)
        np.random.seed(42)
        outputs = [apply_birefringence(E.copy(), 1500, wavelength=self.WAVELENGTH,
                                       model='sectional')[0] for _ in range(50)]
        ex_powers = [np.abs(o[0]) ** 2 for o in outputs]
        assert np.allclose([np.abs(o[0]) ** 2 + np.abs(o[1]) ** 2
                            for o in outputs], 1.0, atol=1e-12)
        assert np.std(ex_powers) > 0.05

    def test_power_conservation_long_distance_sectional(self):
        """The multi-section model must conserve power at long distances
        (former phenomenological regime; 5th-pass PHYS-5 — the model was
        removed, sectional serves all lengths)."""
        E = self._field()
        P_in = np.mean(np.abs(E) ** 2)
        for L_m in [5000, 50000, 100000]:
            E_out = apply_birefringence(E.copy(), L_m, wavelength=self.WAVELENGTH,
                                        model='sectional')
            P_out = np.mean(np.abs(E_out) ** 2)
            assert abs(P_out - P_in) / P_in < 1e-12

    def test_zero_length_identity(self):
        E = self._field(100)
        E_out = apply_birefringence(E.copy(), 0, wavelength=self.WAVELENGTH)
        assert np.allclose(E_out, E)

    def test_temperature_dependence_long_distance_sectional(self):
        E = np.ones((100, 2), dtype=complex)
        np.random.seed(45)
        Js = [apply_birefringence(E.copy(), 50000, wavelength=self.WAVELENGTH,
                                  temperature=T, model='sectional')[0, 0]
              for T in [0, 25, 50]]
        assert not np.allclose(Js[0], Js[1])

    def test_wavelength_dependence_long_distance_sectional(self):
        E = np.ones((100, 2), dtype=complex)
        np.random.seed(46)
        Js = [apply_birefringence(E.copy(), 50000, wavelength=lam,
                                  model='sectional')[0, 0]
              for lam in [1310e-9, 1550e-9]]
        assert not np.allclose(Js[0], Js[1])

    def test_seed_dependence_long_distance_sectional(self):
        E = np.ones((100, 2), dtype=complex)
        np.random.seed(42)
        J1 = apply_birefringence(E.copy(), 50000, wavelength=self.WAVELENGTH,
                                 model='sectional')[0, 0]
        np.random.seed(142)
        J2 = apply_birefringence(E.copy(), 50000, wavelength=self.WAVELENGTH,
                                 model='sectional')[0, 0]
        assert not np.allclose(J1, J2)

    def test_output_variation_long_distance_sectional(self):
        """At 100 km the multi-section model is fully scrambled (uniform
        SU(2) per realization), so a fresh draw must vary the output
        strongly."""
        E = np.array([[1.0, 0.0]], dtype=complex)
        np.random.seed(42)
        outputs = [apply_birefringence(E.copy(), 100e3, wavelength=self.WAVELENGTH,
                                       model='sectional')[0]
                   for _ in range(50)]
        ex_powers = [np.abs(o[0]) ** 2 for o in outputs]
        assert np.std(ex_powers) > 0.05

    def test_auto_equals_sectional_at_all_lengths(self):
        """'auto' and 'sectional' are the same model at every length
        (5th-pass PHYS-5: the phenomenological model was removed)."""
        E = np.ones((10, 2), dtype=complex)
        for L_m in [100, 100000, 122000]:
            np.random.seed(52)
            out_auto = apply_birefringence(E.copy(), L_m,
                                           wavelength=self.WAVELENGTH,
                                           model='auto')
            np.random.seed(52)
            out_sec = apply_birefringence(E.copy(), L_m,
                                          wavelength=self.WAVELENGTH,
                                          model='sectional')
            assert np.allclose(out_auto, out_sec), \
                f"auto should equal sectional at L={L_m} m"

    def test_phenomenological_model_removed_raises(self):
        """The phenomenological model was deleted (PHYS-5, 5th pass); its
        name must fail loudly rather than silently dispatch."""
        E = np.ones((10, 2), dtype=complex)
        for L_m in [100, 50000]:
            with pytest.raises(ValueError):
                apply_birefringence(E.copy(), L_m, wavelength=self.WAVELENGTH,
                                    model='phenomenological')
        with pytest.raises(ValueError):
            FiberRealization(L_m=50000, model='phenomenological')

    def test_enabled_false_returns_unchanged(self):
        E = self._field(100)
        E_out = apply_birefringence(E.copy(), 10000, wavelength=self.WAVELENGTH,
                                    enabled=False)
        assert np.allclose(E_out, E)


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
        # Jones matrices �?" associativity is a property of matrix
        # multiplication in general, not specific to SU(2)).
        J = rng.normal(size=(N, 2, 2)) + 1j * rng.normal(size=(N, 2, 2))

        expected = self._naive_ordered_product(J)
        actual = _ordered_product(J)
        assert np.allclose(actual, expected, atol=1e-10)


class TestBirefringenceMatrixAccessor:
    """FiberRealization.birefringence_matrix() — the quasi-static Jones
    matrix, exposed for receiver-side polarization compensation."""

    def test_returns_unitary_matrix(self):
        fibre = FiberRealization(L_m=50_000, cd=False, pmd=False,
                                 attenuation=False, seed=7)
        J = fibre.birefringence_matrix()
        assert J is not None
        assert J.shape == (2, 2)
        # SU(2): J^dagger @ J = I to float precision
        np.testing.assert_allclose(J.conj().T @ J, np.eye(2), atol=1e-10)

    def test_matches_apply_with_other_impairments_off(self):
        """With CD/PMD/attenuation disabled, apply() must equal J @ E."""
        fibre = FiberRealization(L_m=50_000, cd=False, pmd=False,
                                 attenuation=False, seed=7)
        J = fibre.birefringence_matrix()
        rng = np.random.default_rng(3)
        E = rng.normal(size=(5, 2)) + 1j * rng.normal(size=(5, 2))
        np.testing.assert_allclose(
            fibre.apply(E),
            np.transpose(J @ np.transpose(E)),
            rtol=1e-12, atol=1e-14,
        )

    def test_compensation_roundtrip_recovers_field(self):
        """Applying J^dagger after J must restore the input exactly —
        the basis of active polarization compensation (Duplinskiy et
        al. 2017 calibration loop)."""
        fibre = FiberRealization(L_m=50_000, cd=False, pmd=False,
                                 attenuation=False, seed=7)
        J = fibre.birefringence_matrix()
        Jinv = J.conj().T
        rng = np.random.default_rng(3)
        E = rng.normal(size=(5, 2)) + 1j * rng.normal(size=(5, 2))
        E_fibre = np.transpose(J @ np.transpose(E))
        E_back = np.transpose(Jinv @ np.transpose(E_fibre))
        np.testing.assert_allclose(E_back, E, rtol=1e-10, atol=1e-12)

    def test_returns_none_when_birefringence_disabled(self):
        fibre = FiberRealization(L_m=50_000, birefringence=False, seed=7)
        assert fibre.birefringence_matrix() is None

    def test_quasi_static_same_matrix_every_call(self):
        """The matrix must not change between apply() calls."""
        fibre = FiberRealization(L_m=50_000, cd=False, pmd=False,
                                 attenuation=False, seed=7)
        J1 = fibre.birefringence_matrix()
        rng = np.random.default_rng(3)
        E = rng.normal(size=(3, 2)) + 1j * rng.normal(size=(3, 2))
        for _ in range(5):
            fibre.apply(E)
        J2 = fibre.birefringence_matrix()
        np.testing.assert_array_equal(J1, J2)
