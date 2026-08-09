"""Regression tests for the GOBBY-2 analytic reference model.

These pin two things:

**Traceability.** No constant may be fitted to the data the model exists to
reproduce.  `MU_EFF` in particular must remain *derived* from stated
quantities (`MU * T_INT`) rather than assigned a literal — an earlier
version set it to 0.0793 by inverting Gobby's measured fringe visibilities,
which is fitting to the target, and that produced a phantom discrepancy
that cost real effort to chase.  `test_mu_eff_is_derived_not_fitted` is the
guard against that returning.

**The decomposition.** QBER(L) = e_mod + (1 - V_fringe)/2 against Gobby's
four *measured* points.  If a later refactor breaks the field chain, the
MC-vs-analytic comparison in validate_gobby.py catches it; if this module
itself changes, these catch it.
"""
import numpy as np
import pytest
from src.analytic.gobby_model import (qber, fringe_visibility, erroneous_counts,
                                      signal, MU, SPLIT_RATIO, T_INT, MU_EFF,
                                      P_E, P_E_DARK, E_MOD)

# Gobby's four measured points, and what the derived model gives at each.
# z, model total (%), Gobby measured (%), residual (pp)
MEASURED = [
    (4.4, 3.33, 3.3, +0.03),
    (65.0, 3.79, 3.3, +0.49),
    (101.0, 5.75, 6.0, -0.25),
    (122.0, 9.26, 8.9, +0.36),
]

# Stated fringe visibilities the derivation is checked against.  The model
# is *not* fitted to these -- they are an independent cross-check, so the
# tolerance is looser than a fit would give.
PAPER_VISIBILITIES = [(65.0, 0.9906), (122.0, 0.8840)]

# Sensitivity of the mean |residual| to the signal reading.  The derived
# value is the middle row; the other two are the readings section 19.2
# excluded.
SENSITIVITY_MEAN_ABS_RESIDUAL = {0.1000: 0.52, MU_EFF: 0.28, 0.0385: 2.00}


class TestTraceability:
    """Every constant stated or derived; nothing fitted."""

    def test_mu_eff_is_derived_not_fitted(self):
        """MU_EFF must be computed from MU and T_INT, not assigned.

        This is the guard the module docstring's no-fitting rule needs to
        be more than a comment.  Exact equality, not approx: if someone
        writes `MU_EFF = 0.0793` it fails even though the numbers are
        within 3% of each other.
        """
        assert MU_EFF == MU * T_INT

    def test_t_int_is_derived_from_split_ratio(self):
        """T_INT = 2/(1+r): equalising the arms discards the excess
        reference power.  Derived from the stated 1.6:1, not measured."""
        assert T_INT == 2.0 / (1.0 + SPLIT_RATIO)
        assert T_INT == pytest.approx(0.769231, abs=1e-6)

    def test_stated_constants_match_paper(self):
        assert MU == pytest.approx(0.1)
        assert SPLIT_RATIO == pytest.approx(1.6)
        assert P_E == pytest.approx(8.5e-7)
        assert P_E_DARK == pytest.approx(3.2e-7)
        assert E_MOD == pytest.approx(0.033)

    def test_geometry_agrees_with_the_visibility_inversion(self):
        """Cross-check, NOT an input.

        Inverting Gobby's stated visibilities gives mu_eff/mu = 0.793; the
        geometry predicts 0.769.  Two independent routes agreeing to 3% is
        what makes the derived value credible.  If this ever drifts, one of
        them has been tuned.
        """
        assert abs(T_INT - 0.793) / 0.793 < 0.05


class TestDecomposition:

    @pytest.mark.parametrize("z,total,gobby,residual", MEASURED)
    def test_total_qber(self, z, total, gobby, residual):
        assert qber(z) == pytest.approx(total, abs=0.05)

    @pytest.mark.parametrize("z,total,gobby,residual", MEASURED)
    def test_residual_vs_gobby(self, z, total, gobby, residual):
        assert (qber(z) - gobby) == pytest.approx(residual, abs=0.05)

    @pytest.mark.parametrize("z,v", PAPER_VISIBILITIES)
    def test_fringe_visibility_near_stated_values(self, z, v):
        """Independent cross-check against the paper's stated V.

        Tolerance is 0.005, not the 0.0005 a fitted model would hit: the
        signal is now derived from the geometry, so agreement here is a
        prediction rather than a construction.
        """
        assert fringe_visibility(z) == pytest.approx(v, abs=0.005)

    def test_erroneous_counts_is_fringe_complement(self):
        z = 122.0
        assert erroneous_counts(z) == pytest.approx(
            (1.0 - fringe_visibility(z)) / 2.0 * 100.0)

    def test_mean_abs_residual_within_acceptance(self):
        """Mean |residual| <= 0.30 pp over the four published points.

        The derived signal scores 0.282 pp against the old fitted value's
        0.256 pp -- 0.026 pp is the price of removing the last free
        parameter, and it is worth paying.
        """
        resid = np.abs([qber(z) - gobby for z, _, gobby, _ in MEASURED])
        assert resid.mean() <= 0.30


class TestSignalSensitivity:

    @pytest.mark.parametrize("mu_eff,expected",
                             sorted(SENSITIVITY_MEAN_ABS_RESIDUAL.items()))
    def test_mean_abs_residual_matches_sensitivity_table(self, mu_eff, expected):
        resid = np.abs([qber(z, mu_eff=mu_eff) - gobby
                        for z, _, gobby, _ in MEASURED])
        assert resid.mean() == pytest.approx(expected, abs=0.02)

    def test_naive_and_encoded_only_readings_excluded(self):
        """0.10 (naive total, ignoring the interferometer) and 0.0385
        (the 1.6:1 encoded-only reading) both exceed the 0.30 pp bar; the
        derived value does not."""
        for mu_eff in (0.1000, 0.0385):
            resid = np.abs([qber(z, mu_eff=mu_eff) - gobby
                            for z, _, gobby, _ in MEASURED])
            assert resid.mean() > 0.30
        resid = np.abs([qber(z) - gobby for z, _, gobby, _ in MEASURED])
        assert resid.mean() <= 0.30


class TestOutOfSamplePrediction:
    """Section 19.11: 165.8 km with stray light removed.

    A documented predicted *disagreement* -- Gobby states V = 0.86 and the
    model gives less.  It is recorded, not tuned away; the likely
    explanation (a narrower gate when the detector was electronically
    synchronised) is in section 19.11.
    """

    def test_signal_at_165_8_km(self):
        expected = MU_EFF * 10.0 ** (-0.2 * 165.8 / 10.0) * 0.045
        assert signal(165.8) == pytest.approx(expected, rel=1e-9)

    def test_predicted_visibility_disagrees_with_stated(self):
        v = fringe_visibility(165.8, p_e=P_E_DARK)
        assert v < 0.86, "the documented disagreement must persist"
        assert v == pytest.approx(0.72, abs=0.02)


def _validate_gobby():
    """Import the Gobby validation script, which is an analysis entry point
    rather than a package module and so is not importable by default."""
    import os
    import sys
    d = os.path.join(os.path.dirname(__file__), '..', 'analysis', 'val_gobby')
    d = os.path.abspath(d)
    if d not in sys.path:
        sys.path.insert(0, d)
    import validate_gobby
    return validate_gobby


class TestErroneousCountSaturation:
    """The erroneous-count share cannot explode, whatever the signal does.

    `P_e / (S + 2*P_e)` saturates at 1/2 as `S -> 0`: a background-dominated
    link is a coin flip, not an unbounded error rate.  This is the reason a
    signal-probability error shifts the QBER within a bounded range rather
    than blowing it up, which is what made the 41% `signal_click_prob()`
    omission (GOBBY-7 §24.5) tolerable while it went unnoticed.

    Asserted rather than argued in prose, so a future rewrite of the error
    weighting cannot quietly lose the property.
    """

    @staticmethod
    def share(s, p_e=P_E):
        return p_e / (s + 2.0 * p_e)

    def test_bounded_by_one_half_across_many_decades(self):
        for s in np.logspace(-12, -1, 200):
            assert self.share(s) <= 0.5

    def test_tends_to_one_half_as_signal_vanishes(self):
        assert self.share(0.0) == pytest.approx(0.5)
        assert self.share(1e-15) == pytest.approx(0.5, abs=1e-6)

    def test_monotonically_decreasing_in_signal(self):
        vals = [self.share(s) for s in np.logspace(-12, -1, 100)]
        assert all(a >= b for a, b in zip(vals, vals[1:]))

    def test_a_large_signal_error_cannot_explode_the_share(self):
        """Even a 41% signal error moves it by a bounded amount."""
        for s in np.logspace(-6, -2, 50):
            assert abs(self.share(s * 0.71) - self.share(s)) < 0.15


class TestSignalClickProbCarriesTInt:
    """`validate_gobby.signal_click_prob` must apply the interferometer
    transmission.

    It did not, for the whole of GOBBY-1 through GOBBY-6, which put it 41%
    above what the chain delivers and made `model_qber` and
    `predicted_visibility` optimistic with it.  This is the regression
    guard: `T_INT` is derived from the stated 1.6:1 split and is the same
    value the analytic model carries as `MU_EFF`, so the two must agree.
    """

    def test_matches_mu_eff_at_zero_distance(self):
        G = _validate_gobby()
        assert G.signal_click_prob(0) == pytest.approx(MU_EFF * G.ETA_BOB,
                                                       rel=1e-12)

    def test_is_not_the_bare_mu_form(self):
        """The specific bug: MU where MU_EFF belongs."""
        G = _validate_gobby()
        assert G.signal_click_prob(0) < 0.99 * MU * G.ETA_BOB

    def test_reproduces_the_stated_visibilities(self):
        """Gobby state >0.99 at 65 km and 0.884 at 122 km.  The 122 km
        figure is the sharper test, being a value rather than a bound, and
        is what confirmed the omission."""
        G = _validate_gobby()
        assert G.predicted_visibility(65) > 0.99
        assert G.predicted_visibility(122) == pytest.approx(0.884, abs=0.01)

    def test_agrees_with_the_monte_carlo_within_the_stated_bound(self):
        """Guards against the 41% class of regression, not a measurement.

        The reference figure is 1.019 +/- 0.012, from a 4e6-pulse run with
        dead time and afterpulsing off (7,056 sifted).  A unit-test-sized
        run cannot reach that precision -- at 1e6 pulses the sifted count
        is ~1,700, so the statistical error alone is ~2.4% -- hence the
        tolerance is looser than the agreement.

        Two defects have driven this ratio away from unity and both are
        guarded here: the missing `T_INT` put it at 0.71, and the SPAD's
        `eta*(1-exp(-mu))` detection probability put it at 0.95.
        """
        G = _validate_gobby()
        from src.protocols.bb84_time_bin import simulate_bb84_time_bin as sim
        N = 1_000_000
        r = sim(num_bits=N, fiber_length=0, alpha_dB=G.ALPHA_dB,
                repetition_rate=G.REP_RATE, mu=G.MU, spad_eta=G.ETA_BOB,
                gate_width=G.GATE_WIDTH, pulse_width=G.PULSE_WIDTH,
                interferometer='polarisation_multiplexed',
                split_ratio=G.SPLIT_RATIO, dark_count_rate=0.0,
                afterpulse_prob=0.0, dead_time=0.0, seed=7)
        s_mc = r['n_sifted'] / (N * 0.5)
        ratio = s_mc / G.signal_click_prob(0)
        assert ratio == pytest.approx(1.0, rel=0.10)
        # And explicitly above both known-defect values, so a revert of
        # either is caught even if the tolerance above is ever widened.
        assert ratio > 0.97
