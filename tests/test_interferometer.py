import numpy as np
import pytest
from src.channel.interferometer import AsymmetricMZI


@pytest.fixture
def gaussian_pulse():
    """Short Gaussian pulse field (100 ps FWHM), shape (N, 2)."""
    dt = 2e-12
    sigma = 50e-12 / (2 * np.sqrt(2 * np.log(2)))  # 50 ps FWHM
    n = 2000  # 4 ns total
    t = np.arange(n, dtype=float) * dt
    envelope = np.exp(-0.5 * (t - 1e-9) ** 2 / sigma ** 2)
    E = np.zeros((n, 2), dtype=complex)
    E[:, 0] = envelope * np.sqrt(0.5)
    E[:, 1] = envelope * np.sqrt(0.5)
    return E, dt


class TestAsymmetricMZIConstruction:

    def test_positive_delay_required(self):
        """delay <= 0 should raise ValueError."""
        with pytest.raises(ValueError):
            AsymmetricMZI(delay=0)
        with pytest.raises(ValueError):
            AsymmetricMZI(delay=-1e-9)

    def test_default_mode_is_encoder(self):
        mzi = AsymmetricMZI(delay=1e-9)
        assert mzi.mode == 'encoder'

    def test_decoder_mode(self):
        mzi = AsymmetricMZI(delay=1e-9, mode='decoder')
        assert mzi.mode == 'decoder'

    def test_insertion_loss_default(self):
        mzi = AsymmetricMZI(delay=1e-9)
        assert mzi._il_lin == 1.0

    def test_insertion_loss_conversion(self):
        mzi = AsymmetricMZI(delay=1e-9, insertion_loss_db=3.0)
        expected = 10.0 ** (-3.0 / 10.0)
        assert np.isclose(mzi._il_lin, expected, rtol=1e-6)


class TestAsymmetricMZIEncoder:

    def test_output_shape(self, gaussian_pulse):
        """Encoder output should match input shape."""
        E, dt = gaussian_pulse
        mzi = AsymmetricMZI(delay=1e-9)
        E_out = mzi.modulate(E, dt)
        assert E_out.shape == E.shape
        assert np.iscomplexobj(E_out)

    def test_creates_two_time_bins(self, gaussian_pulse):
        """Encoder output should have two distinct time-separated pulses."""
        E, dt = gaussian_pulse
        delay = 1e-9
        mzi = AsymmetricMZI(delay=delay)
        E_out = mzi.modulate(E, dt)
        P = np.sum(np.abs(E_out) ** 2, axis=1)
        # Find peaks (> 10% of max)
        threshold = 0.1 * np.max(P)
        peaks = []
        for i in range(1, len(P) - 1):
            if P[i] > threshold and P[i] > P[i - 1] and P[i] > P[i + 1]:
                peaks.append(i)
        assert len(peaks) >= 2, f"Expected ≥2 peaks, found {len(peaks)}"

    def test_delay_between_bins(self, gaussian_pulse):
        """Time between output peaks should match delay parameter."""
        E, dt = gaussian_pulse
        delay = 1e-9
        mzi = AsymmetricMZI(delay=delay)
        E_out = mzi.modulate(E, dt)
        P = np.sum(np.abs(E_out) ** 2, axis=1)
        # Find first two peaks
        threshold = 0.1 * np.max(P)
        peaks = []
        for i in range(1, len(P) - 1):
            if P[i] > threshold and P[i] > P[i - 1] and P[i] > P[i + 1]:
                peaks.append(i)
                if len(peaks) >= 2:
                    break
        assert len(peaks) >= 2
        measured_delay = (peaks[1] - peaks[0]) * dt
        assert np.isclose(measured_delay, delay, rtol=0.05)

    def test_power_conservation(self, gaussian_pulse):
        """Encoder should conserve total power (within numerical precision)."""
        E, dt = gaussian_pulse
        P_in = np.mean(np.sum(np.abs(E) ** 2, axis=1))
        mzi = AsymmetricMZI(delay=1e-9)
        E_out = mzi.modulate(E, dt)
        P_out = np.mean(np.sum(np.abs(E_out) ** 2, axis=1))
        # 3 dB intrinsic loss (each half of light goes to different path;
        # at peak times, only one arm's light is present in the output)
        assert np.isclose(P_out, P_in, rtol=0.01)

    def test_phase_applies_to_delayed_bin(self, gaussian_pulse):
        """Phase shift should only affect the delayed (second) bin."""
        E, dt = gaussian_pulse
        mzi = AsymmetricMZI(delay=1e-9)
        E_ref = mzi.modulate(E, dt, phase=0)
        E_ph = mzi.modulate(E, dt, phase=np.pi)
        P_ref = np.sum(np.abs(E_ref) ** 2, axis=1)
        P_ph = np.sum(np.abs(E_ph) ** 2, axis=1)
        # First bin should be identical (same power)
        delay_samples = int(1e-9 / dt)
        first_bin_end = delay_samples // 2
        assert np.allclose(P_ref[:first_bin_end], P_ph[:first_bin_end])
        # Second bin should differ (phase changes interference at combiner?)
        # Actually the phase is applied to the long arm only, which becomes
        # the second time bin.  Power should be the same (phase doesn't
        # change |exp(jφ)|^2 = 1) but the field values differ.
        E_long_ref = np.roll(E / np.sqrt(2), delay_samples, axis=0)
        E_long_ref[:delay_samples] = 0
        E_long_ph = np.roll(E / np.sqrt(2), delay_samples, axis=0)
        E_long_ph[:delay_samples] = 0
        E_long_ph *= np.exp(1j * np.pi)
        expected_ref = E / np.sqrt(2) + E_long_ref
        expected_ph = E / np.sqrt(2) + E_long_ph
        assert np.allclose(E_ref, expected_ref)
        assert np.allclose(E_ph, expected_ph)

    def test_varying_delay(self, gaussian_pulse):
        """Different delays should produce different bin spacings."""
        E, dt = gaussian_pulse
        delays = [0.5e-9, 1.0e-9, 2.0e-9]
        spacings = []
        for d in delays:
            mzi = AsymmetricMZI(delay=d)
            E_out = mzi.modulate(E, dt)
            P = np.sum(np.abs(E_out) ** 2, axis=1)
            threshold = 0.1 * np.max(P)
            peaks = []
            for i in range(1, len(P) - 1):
                if P[i] > threshold and P[i] > P[i - 1] and P[i] > P[i + 1]:
                    peaks.append(i)
                    if len(peaks) >= 2:
                        break
            assert len(peaks) >= 2
            spacings.append((peaks[1] - peaks[0]) * dt)
        # Verify monotonic: longer delay → larger spacing
        assert spacings[0] < spacings[1] < spacings[2]

    def test_no_phase_output(self, gaussian_pulse):
        """Output should be real when phase=None and input is real."""
        E, dt = gaussian_pulse
        mzi = AsymmetricMZI(delay=1e-9)
        E_out = mzi.modulate(E, dt, phase=None)
        E_zero = mzi.modulate(E, dt, phase=0.0)
        assert np.allclose(E_out, E_zero)


class TestAsymmetricMZIDecoder:

    def test_output_is_tuple(self, gaussian_pulse):
        """Decoder should return (E_c, E_d) tuple."""
        E, dt = gaussian_pulse
        mzi = AsymmetricMZI(delay=1e-9, mode='decoder')
        result = mzi.modulate(E, dt)
        assert isinstance(result, tuple)
        assert len(result) == 2

    def test_output_shapes(self, gaussian_pulse):
        """Both decoder outputs should match input shape."""
        E, dt = gaussian_pulse
        mzi = AsymmetricMZI(delay=1e-9, mode='decoder')
        E_c, E_d = mzi.modulate(E, dt)
        assert E_c.shape == E.shape
        assert E_d.shape == E.shape

    def test_interference_fringes(self, gaussian_pulse):
        """Decoder interference region follows P_c ∝ 1+cos(Δφ)."""
        E, dt = gaussian_pulse
        delay = 0.5e-9
        enc = AsymmetricMZI(delay=delay)
        E_in = enc.modulate(E, dt, phase=0)

        dec = AsymmetricMZI(delay=delay, mode='decoder')
        phases = np.linspace(0, 2 * np.pi, 13)
        P_c_vals = []
        P_d_vals = []
        for phi_b in phases:
            E_c, E_d = dec.modulate(E_in, dt, phase=phi_b)
            pc, pd = self._interference_power(E_c, E_d, E, dt, delay)
            P_c_vals.append(pc)
            P_d_vals.append(pd)

        P_c = np.array(P_c_vals)
        P_d = np.array(P_d_vals)

        # Normalise to [0, 1]
        P_c_norm = (P_c - P_c.min()) / (P_c.max() - P_c.min())
        P_d_norm = (P_d - P_d.min()) / (P_d.max() - P_d.min())

        expected_c = (1.0 + np.cos(phases)) / 2.0
        expected_d = (1.0 - np.cos(phases)) / 2.0

        assert np.allclose(P_c_norm, expected_c, atol=0.15)
        assert np.allclose(P_d_norm, expected_d, atol=0.15)

    def test_interference_sum_conserved(self, gaussian_pulse):
        """Total power (P_c + P_d) should be conserved."""
        E, dt = gaussian_pulse
        delay = 0.5e-9
        enc = AsymmetricMZI(delay=delay)
        E_in = enc.modulate(E, dt, phase=0)

        P_in = np.mean(np.sum(np.abs(E_in) ** 2, axis=1))
        dec = AsymmetricMZI(delay=delay, mode='decoder')
        for phi_b in [0, np.pi / 4, np.pi / 2, np.pi]:
            E_c, E_d = dec.modulate(E_in, dt, phase=phi_b)
            P_out = np.mean(
                np.sum(np.abs(E_c) ** 2, axis=1) +
                np.sum(np.abs(E_d) ** 2, axis=1)
            )
            assert np.isclose(P_out, P_in, rtol=0.01)

    @staticmethod
    def _interference_power(E_c, E_d, E_original, dt, delay):
        """Compute power at the interference peak (overlap region)."""
        delay_samples = int(delay / dt)
        P_c = np.sum(np.abs(E_c) ** 2, axis=1)
        P_d = np.sum(np.abs(E_d) ** 2, axis=1)
        # Find the max of P_c + P_d (the interference region)
        total = P_c + P_d
        peak = np.argmax(total)
        # Window around the peak (3× the pulse sigma)
        half = max(1, delay_samples // 4)
        start = max(0, peak - half)
        end = min(len(total), peak + half + 1)
        return np.mean(P_c[start:end]), np.mean(P_d[start:end])

    def test_constructive_all_in_one_port(self, gaussian_pulse):
        """When Δφ = 0, constructive port should dominate at overlap."""
        E, dt = gaussian_pulse
        delay = 0.5e-9
        enc = AsymmetricMZI(delay=delay)
        E_in = enc.modulate(E, dt, phase=0)
        dec = AsymmetricMZI(delay=delay, mode='decoder')
        E_c, E_d = dec.modulate(E_in, dt, phase=0)
        P_c, P_d = self._interference_power(E_c, E_d, E, dt, delay)
        assert P_c > 3 * P_d

    def test_destructive_all_in_other_port(self, gaussian_pulse):
        """When Δφ = π, destructive port should dominate at overlap."""
        E, dt = gaussian_pulse
        delay = 0.5e-9
        enc = AsymmetricMZI(delay=delay)
        E_in = enc.modulate(E, dt, phase=0)
        dec = AsymmetricMZI(delay=delay, mode='decoder')
        E_c, E_d = dec.modulate(E_in, dt, phase=np.pi)
        P_c, P_d = self._interference_power(E_c, E_d, E, dt, delay)
        assert P_d > 3 * P_c


class TestAsymmetricMZIRoundtrip:

    def test_encoder_decoder_identity(self, gaussian_pulse):
        """Encoder + matching decoder fringe should follow 1+cos(Δφ)."""
        E, dt = gaussian_pulse
        delay = 0.5e-9
        enc = AsymmetricMZI(delay=delay)
        dec = AsymmetricMZI(delay=delay, mode='decoder')

        for phi_a in [0, np.pi / 2, np.pi, 3 * np.pi / 2]:
            E_in = enc.modulate(E, dt, phase=phi_a)
            for phi_b in [0, np.pi / 4, np.pi / 2]:
                E_c, E_d = dec.modulate(E_in, dt, phase=phi_b)
                P_c, P_d = TestAsymmetricMZIDecoder._interference_power(
                    E_c, E_d, E, dt, delay
                )
                total = P_c + P_d
                if total > 0:
                    ratio_c = P_c / total
                    expected = (1 + np.cos(phi_a - phi_b)) / 2
                    assert np.isclose(ratio_c, expected, atol=0.15)


class TestAsymmetricMZIVisibility:

    def test_visibility_default_is_ideal(self):
        """Default visibility should be 1.0 (ideal)."""
        mzi = AsymmetricMZI(delay=1e-9, mode='decoder')
        assert mzi.visibility == 1.0

    def test_visibility_validation(self):
        """Visibility outside (0, 1] should raise ValueError."""
        with pytest.raises(ValueError):
            AsymmetricMZI(delay=1e-9, mode='decoder', visibility=0.0)
        with pytest.raises(ValueError):
            AsymmetricMZI(delay=1e-9, mode='decoder', visibility=-0.5)
        with pytest.raises(ValueError):
            AsymmetricMZI(delay=1e-9, mode='decoder', visibility=1.5)

    def test_visibility_matches_fringe_contrast(self, gaussian_pulse):
        """Measured fringe contrast (max-min)/(max+min) must equal V."""
        E, dt = gaussian_pulse
        delay = 0.5e-9
        enc = AsymmetricMZI(delay=delay)

        for V in [1.0, 0.934, 0.8, 0.5]:
            dec = AsymmetricMZI(delay=delay, mode='decoder', visibility=V)
            E_in = enc.modulate(E, dt, phase=0)
            P_c_vals = []
            for phi_b in np.linspace(0, 2 * np.pi, 25):
                E_c, E_d = dec.modulate(E_in, dt, phase=phi_b)
                pc, _ = TestAsymmetricMZIDecoder._interference_power(
                    E_c, E_d, E, dt, delay)
                P_c_vals.append(pc)
            P_c = np.array(P_c_vals)
            contrast = (P_c.max() - P_c.min()) / (P_c.max() + P_c.min())
            assert np.isclose(contrast, V, atol=0.02), \
                f"V={V}: measured contrast {contrast:.4f}"

    def test_visibility_optical_misalignment_error(self, gaussian_pulse):
        """At Delta_phi = 0, e_opt = P_d/(P_c+P_d) must equal (1-V)/2."""
        E, dt = gaussian_pulse
        delay = 0.5e-9
        enc = AsymmetricMZI(delay=delay)
        E_in = enc.modulate(E, dt, phase=0)

        for V in [1.0, 0.934, 0.8, 0.5]:
            dec = AsymmetricMZI(delay=delay, mode='decoder', visibility=V)
            E_c, E_d = dec.modulate(E_in, dt, phase=0)
            P_c, P_d = TestAsymmetricMZIDecoder._interference_power(
                E_c, E_d, E, dt, delay)
            e_opt = P_d / (P_c + P_d)
            assert np.isclose(e_opt, (1.0 - V) / 2.0, atol=0.02), \
                f"V={V}: e_opt {e_opt:.4f} vs {(1.0-V)/2.0:.4f}"

    def test_visibility_power_conserved(self, gaussian_pulse):
        """Finite visibility must not change total power at the fringe."""
        E, dt = gaussian_pulse
        delay = 0.5e-9
        enc = AsymmetricMZI(delay=delay)
        E_in = enc.modulate(E, dt, phase=0)

        for V in [1.0, 0.934, 0.5]:
            dec = AsymmetricMZI(delay=delay, mode='decoder', visibility=V)
            E_c, E_d = dec.modulate(E_in, dt, phase=0)
            P_c, P_d = TestAsymmetricMZIDecoder._interference_power(
                E_c, E_d, E, dt, delay)
            total = P_c + P_d
            # Ideal decoder at same delay gives the reference total
            dec_ideal = AsymmetricMZI(delay=delay, mode='decoder')
            E_ci, E_di = dec_ideal.modulate(E_in, dt, phase=0)
            P_ci, P_di = TestAsymmetricMZIDecoder._interference_power(
                E_ci, E_di, E, dt, delay)
            assert np.isclose(total, P_ci + P_di, rtol=0.01)

    def test_visibility_one_is_backward_compatible(self, gaussian_pulse):
        """visibility=1.0 must reproduce the ideal 50:50 combiner exactly."""
        E, dt = gaussian_pulse
        delay = 0.5e-9
        enc = AsymmetricMZI(delay=delay)
        E_in = enc.modulate(E, dt, phase=np.pi / 3)

        dec_ideal = AsymmetricMZI(delay=delay, mode='decoder')
        dec_v1 = AsymmetricMZI(delay=delay, mode='decoder', visibility=1.0)
        for phi_b in [0, np.pi / 4, np.pi / 2]:
            E_c1, E_d1 = dec_ideal.modulate(E_in, dt, phase=phi_b)
            E_c2, E_d2 = dec_v1.modulate(E_in, dt, phase=phi_b)
            assert np.allclose(E_c1, E_c2, atol=1e-12)
            assert np.allclose(E_d1, E_d2, atol=1e-12)

    def test_phase_error_shifts_fringe(self, gaussian_pulse):
        """A static phase error must shift the fringe by -phi_err."""
        E, dt = gaussian_pulse
        delay = 0.5e-9
        enc = AsymmetricMZI(delay=delay)
        E_in = enc.modulate(E, dt, phase=0)

        for phi_err in [0.0, 0.3, 1.0]:
            dec = AsymmetricMZI(delay=delay, mode='decoder',
                                phase_error=phi_err)
            P_c_vals = []
            for phi_b in np.linspace(0, 2 * np.pi, 49):
                E_c, E_d = dec.modulate(E_in, dt, phase=phi_b)
                pc, _ = TestAsymmetricMZIDecoder._interference_power(
                    E_c, E_d, E, dt, delay)
                P_c_vals.append(pc)
            P_c = np.array(P_c_vals)
            # Fringe peak moves to phi_b = -phi_err (mod 2pi)
            phi_grid = np.linspace(0, 2 * np.pi, 49)
            peak_phi = phi_grid[np.argmax(P_c)]
            expected = (-phi_err) % (2 * np.pi)
            assert np.isclose(peak_phi, expected, atol=0.15), \
                f"phi_err={phi_err}: peak at {peak_phi:.3f}, expected {expected:.3f}"

    def test_phase_error_default_zero(self):
        """Default phase_error should be 0.0."""
        mzi = AsymmetricMZI(delay=1e-9, mode='decoder')
        assert mzi.phase_error == 0.0


class TestAsymmetricMZILoss:

    def test_insertion_loss_reduces_power(self, gaussian_pulse):
        """Insertion loss should scale output power correctly."""
        E, dt = gaussian_pulse
        P_in = np.mean(np.sum(np.abs(E) ** 2, axis=1))
        mzi = AsymmetricMZI(delay=1e-9, mode='encoder', insertion_loss_db=3.0)
        E_out = mzi.modulate(E, dt)
        P_out = np.mean(np.sum(np.abs(E_out) ** 2, axis=1))
        # With 3 dB loss, output should be ~half
        assert np.isclose(P_out, P_in * 0.5, rtol=0.02)

    def test_insertion_loss_decoder(self, gaussian_pulse):
        """Insertion loss should affect both decoder ports."""
        E, dt = gaussian_pulse
        delay = 0.5e-9
        enc = AsymmetricMZI(delay=delay)
        E_in = enc.modulate(E, dt)
        dec = AsymmetricMZI(
            delay=delay, mode='decoder', insertion_loss_db=3.0
        )
        E_c, E_d = dec.modulate(E_in, dt)
        P_out = np.mean(
            np.sum(np.abs(E_c) ** 2, axis=1) +
            np.sum(np.abs(E_d) ** 2, axis=1)
        )
        P_in = np.mean(np.sum(np.abs(E_in) ** 2, axis=1))
        assert np.isclose(P_out, P_in * 0.5, rtol=0.02)
