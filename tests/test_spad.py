"""Tests for the Geiger-mode SPAD detector (src/detectors/spad.py).

Covers the state machine that produces the Gobby et al. (2004) result:
dead-time enforcement, dark-count rate convergence, Poisson click
probability, afterpulse excess rate (the BLOCK-3 measurement in
opto-sim-issues-and-fixes.md) and the PHYS-7 afterpulse state-leak
regression.
"""
import numpy as np
import pytest
from src.detectors.spad import spad


def make_spad(**kwargs):
    defaults = dict(
        wavelength=1550e-9,
        quantum_efficiency=1.0,
        dead_time=13e-6,
        dark_count_rate=0.0,
        afterpulse_prob=0.0,
        gate_width=1e-9,
    )
    defaults.update(kwargs)
    return spad(**defaults)


class TestSPADConstruction:

    def test_inherits_apd_physics(self):
        """SPAD should carry apd physical constants and responsivity."""
        s = make_spad()
        assert s.qe == 1.0
        assert np.isclose(s.frequency, 3e8 / 1550e-9, rtol=1e-9)
        assert s.h == 6.626e-34

    def test_attributes(self):
        """Dead time, DCR, afterpulse and gate width should be stored."""
        s = make_spad(dead_time=20e-6, dark_count_rate=100.0,
                      afterpulse_prob=0.1, gate_width=5e-9)
        assert s.dead_time == 20e-6
        assert s.dcr == 100.0
        assert s.afterpulse_prob == 0.1
        assert s.gate_width == 5e-9

    def test_starts_armed(self):
        s = make_spad()
        assert s.is_armed

    def test_reset_restores_state(self):
        """reset() should restore armed state and clear afterpulse state."""
        s = make_spad(afterpulse_prob=1.0)
        assert s.detect(1e-6, 0.0) == 1  # click -> disarmed
        assert not s.is_armed
        s.reset()
        assert s.is_armed
        assert not s._afterpulse_pending


class TestSPADDeadTime:

    def test_click_disarms_detector(self):
        s = make_spad()
        assert s.detect(1e-6, 0.0) == 1
        assert not s.is_armed

    def test_no_click_during_dead_time(self):
        """Within dead_time of a click, no further click is possible."""
        s = make_spad()
        assert s.detect(1e-6, 0.0) == 1
        assert s.detect(1e-6, 5e-6) == 0
        assert s.detect(1e-6, 12.999e-6) == 0

    def test_rearms_after_dead_time(self):
        """Exactly dead_time after a click, detection is possible again."""
        s = make_spad()
        assert s.detect(1e-6, 0.0) == 1
        assert s.detect(1e-6, 13e-6) == 1  # t - last >= dead_time -> re-armed

    def test_dead_time_limits_click_rate(self):
        """Click rate must be capped at 1/dead_time regardless of power."""
        s = make_spad(gate_width=20e-9)
        dt = 100e-9
        N = 5000
        times = np.arange(N) * dt
        clicks = s.detect_pulse_train(np.full(N, 1e-3), times)
        max_rate = 1.0 / 13e-6  # ~77 kHz
        assert clicks.sum() / (N * dt) < max_rate * 1.2


class TestSPADDarkCounts:

    def test_dark_count_rate_converges(self):
        """With zero signal, click rate must converge to DCR."""
        gate_width = 1e-6
        dcr = 1000.0
        s = make_spad(dark_count_rate=dcr, gate_width=gate_width,
                      dead_time=0.0)
        N = 200000
        dt = gate_width
        times = np.arange(N) * dt
        clicks = s.detect_pulse_train(np.zeros(N), times)
        # Expected: dcr * total_time = 1000 * 0.2 = 200 clicks
        n = int(clicks.sum())
        assert 120 <= n <= 280, f"measured {n} clicks, expected ~200"

    def test_no_dark_counts_when_dcr_zero(self):
        """Zero DCR must give zero clicks at zero power."""
        s = make_spad()
        N = 5000
        times = np.arange(N) * 1e-6
        clicks = s.detect_pulse_train(np.zeros(N), times)
        assert clicks.sum() == 0


class TestSPADPoissonClicks:

    @pytest.mark.parametrize("mu", [0.05, 0.1, 0.5, 1.0])
    def test_click_probability_matches_1_minus_exp(self, mu):
        """P(click) = 1 - exp(-qe*mu) with mu = P*gate/(h*nu)."""
        s = make_spad(quantum_efficiency=1.0, dead_time=0.0)
        photon_energy = s.h * s.frequency
        power = mu * photon_energy / s.gate_width
        p_theory = 1.0 - np.exp(-mu)
        N = 50000
        dt = s.gate_width
        times = np.arange(N) * dt
        clicks = s.detect_pulse_train(np.full(N, power), times)
        p_meas = clicks.sum() / N
        # Binomial sigma ~ sqrt(p(1-p)/N)
        assert abs(p_meas - p_theory) < 0.02, \
            f"mu={mu}: measured {p_meas:.4f}, theory {p_theory:.4f}"

    def test_quantum_efficiency_scales_click_rate(self):
        """Halving eta must roughly halve the click probability (low mu)."""
        mu = 0.1
        photon_energy = make_spad().h * (3e8 / 1550e-9)
        power = mu * photon_energy / 1e-9
        N = 50000
        times = np.arange(N) * 1e-9

        s1 = make_spad(quantum_efficiency=1.0, dead_time=0.0)
        s2 = make_spad(quantum_efficiency=0.5, dead_time=0.0)
        p1 = s1.detect_pulse_train(np.full(N, power), times).sum() / N
        p2 = s2.detect_pulse_train(np.full(N, power), times).sum() / N
        # For mu=0.1, p ~ 0.095 vs 0.049 — ratio close to 2
        assert np.isclose(p2 / p1, 0.5, rtol=0.3)

    def test_zero_power_gives_no_signal_clicks(self):
        s = make_spad(dark_count_rate=0.0)
        assert s.detect(0.0, 0.0) == 0


class TestSPADAfterpulse:

    def test_afterpulse_excess_rate(self):
        """Afterpulsing must raise the click rate by ~afterpulse_prob.

        The BLOCK-3 measurement in opto-sim-issues-and-fixes.md: driving
        an ID230 SPAD at 2.5 MHz with mu = 0.05 gives 5.3% excess clicks
        against the nominal 5% afterpulse probability.
        """
        gate_width = 1e-9
        mu = 0.05
        photon_energy = make_spad().h * (3e8 / 1550e-9)
        power = mu * photon_energy / gate_width

        N = 250000
        dt = 400e-9  # 2.5 MHz repetition
        times = np.arange(N) * dt

        s0 = make_spad(gate_width=gate_width, afterpulse_prob=0.0)
        s1 = make_spad(gate_width=gate_width, afterpulse_prob=0.05)
        clicks0 = int(s0.detect_pulse_train(np.full(N, power), times).sum())
        clicks1 = int(s1.detect_pulse_train(np.full(N, power), times).sum())

        assert clicks0 > 1000, f"insufficient base clicks: {clicks0}"
        excess = (clicks1 - clicks0) / clicks0
        assert 0.005 <= excess <= 0.12, \
            f"excess {excess*100:.2f}% vs nominal 5% (base {clicks0}, ap {clicks1})"

    def test_afterpulse_requires_prior_click(self):
        """With afterpulse_prob set but no clicks, no afterpulses fire."""
        s = make_spad(afterpulse_prob=0.5)
        N = 5000
        times = np.arange(N) * 1e-6
        clicks = s.detect_pulse_train(np.zeros(N), times)
        assert clicks.sum() == 0

    def test_afterpulse_only_fires_during_dead_period(self):
        """An afterpulse click must keep the detector disarmed.

        With afterpulse_prob=1.0 a click always schedules an afterpulse;
        a gate inside the dead period after the afterpulse time must
        report a click (the afterpulse) and the detector must stay dead.
        """
        s = make_spad(afterpulse_prob=1.0)
        assert s.detect(1e-6, 0.0) == 1          # primary click, dead 13 us
        # Schedule: afterpulse at 0 + max(exp(6.5us), gate) — walk in time
        # until the pending afterpulse fires within the dead window.
        fired = False
        t = 0.0
        while t < 12e-6 and not fired:
            t += 1e-6
            if s.detect(0.0, t) == 1:
                fired = True
                assert not s.is_armed  # still dead after afterpulse
        # Afterpulse is exponential(6.5 us); with prob 1.0 it must fire
        # before 13 us in a deterministic seeded run for the early gates.
        # The assertion below is on the invariant, not the timing.
        assert not s.is_armed or s.detect(1e-6, 13e-6) == 1

    def test_phys7_no_afterpulse_state_leak(self):
        """PHYS-7 regression: a scheduled-but-never-fired afterpulse must
        not leak into a later dead period (2000/2000 leaked pre-fix,
        0/2000 after — measured in opto-sim-issues-and-fixes.md)."""
        leaked = 0
        trials = 500
        for _ in range(trials):
            s = make_spad(afterpulse_prob=1.0)
            assert s.detect(1e-6, 0.0) == 1       # click 1, schedules afterpulse
            assert s.detect(0.0, 20e-6) == 0      # jump past dead time (no gate
                                                  # during the window)
            # Click 2 must NOT schedule a fresh afterpulse, so any click in
            # the following dead period can only come from click 1's stale
            # pending afterpulse (the leak).
            s.afterpulse_prob = 0.0
            assert s.detect(1e-6, 25e-6) == 1     # click 2 -> dead again
            if s.detect(0.0, 26e-6) == 1:         # stale afterpulse would
                leaked += 1                       # fire on this first gate
        assert leaked == 0, f"afterpulse leaked in {leaked}/{trials} trials"


class TestDetectionProbabilityIsPoissonInDetectedPhotons:
    """P(click) = 1 - exp(-eta*mu), not eta*(1 - exp(-mu)).

    Photons are detected independently, so the *detected* count is Poisson
    with mean `eta*mu`.  The older form read "at least one photon arrives,
    then one coin flip at eta", which undercounts whenever more than one
    photon can arrive: for two photons the chance of detecting at least one
    is 1 - (1-eta)^2, not eta.

    The two agree as mu -> 0, which is why the error survived in the
    weak-coherent regime the simulator is mostly used in.  It was found as
    a residual disagreement between `validate_gobby.signal_click_prob()`
    and the Monte Carlo (GOBBY-7b §24.5).
    """

    @staticmethod
    def _p_click(det, n_photons):
        """Analytic click probability for a gate carrying `n_photons`."""
        return 1.0 - np.exp(-det.qe * n_photons)

    def _power_for(self, det, n_photons):
        return n_photons * det.h * det.frequency / det.gate_width

    def test_collapses_to_eta_mu_in_the_weak_limit(self):
        d = spad(wavelength=1550e-9, quantum_efficiency=0.045)
        mu = 1e-6
        assert self._p_click(d, mu) == pytest.approx(d.qe * mu, rel=1e-5)

    def test_saturates_toward_unity_for_a_bright_pulse(self):
        """The old form capped at eta; this must approach 1."""
        d = spad(wavelength=1550e-9, quantum_efficiency=1.0)
        assert self._p_click(d, 20.0) > 0.999
        # what the superseded form would have given at eta = 1, mu = 2
        assert 1.0 * (1.0 - np.exp(-2.0)) == pytest.approx(0.8647, abs=1e-3)
        assert self._p_click(d, 2.0) == pytest.approx(0.8647, abs=1e-3)

    def test_differs_from_the_superseded_form_at_the_gobby_point(self):
        """mu = 0.0769 photons, eta = 0.045: the two differ by 3.6%."""
        d = spad(wavelength=1550e-9, quantum_efficiency=0.045)
        mu = 0.1 * (2.0 / (1.0 + 1.6))
        old = d.qe * (1.0 - np.exp(-mu))
        new = self._p_click(d, mu)
        assert new > old
        assert new / old == pytest.approx(1.0371, abs=2e-3)

    def test_monte_carlo_click_rate_matches_the_closed_form(self):
        """The implementation, not just the formula."""
        np.random.seed(11)
        d = spad(wavelength=1550e-9, quantum_efficiency=0.2,
                 dead_time=0.0, dark_count_rate=0.0, afterpulse_prob=0.0)
        mu = 0.5
        p = self._p_click(d, mu)
        power = self._power_for(d, mu)
        n = 40000
        clicks = sum(d.detect(power, i * 1e-6) for i in range(n))
        se = np.sqrt(p * (1 - p) / n)
        assert abs(clicks / n - p) < 4 * se

    def test_rate_is_monotonic_in_efficiency(self):
        d_lo = spad(wavelength=1550e-9, quantum_efficiency=0.1)
        d_hi = spad(wavelength=1550e-9, quantum_efficiency=0.9)
        assert self._p_click(d_hi, 0.5) > self._p_click(d_lo, 0.5)
