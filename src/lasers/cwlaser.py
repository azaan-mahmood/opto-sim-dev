import numpy as np
from numpy.typing import NDArray

# --- LITERATURE SOURCES ---
# [1] Henry, C. H., "Theory of the Linewidth of Semiconductor Lasers",
#     IEEE J. Quantum Electron., vol. QE-18, no. 2, pp. 259-264, 1982.
# [2] Coldren, L. A., Corzine, S. W., & Mashanovitch, M. L.,
#     "Diode Lasers and Photonic Integrated Circuits", 2nd ed., Wiley, 2012.
#     §5.3: Small-signal modulation response and RIN (Eq. 5.3.30-5.3.38).
# [3] Yariv, A., "Optical Electronics", 4th ed., Saunders, 1991, Ch. 6.
# [4] Schawlow, A. L. & Townes, C. H., "Infrared and Optical Masers",
#     Phys. Rev., vol. 112, no. 6, pp. 1940-1949, 1958.
# [5] Petermann, K., "Laser Diode Modulation and Noise", Kluwer, 1988.
#     Ch. 7: Relaxation oscillations and RIN spectral density.
# [6] Gobby, C., Yuan, Z. L. & Shields, A. J., "Quantum key distribution
#     over 122 km of standard telecom fiber", Appl. Phys. Lett. 84,
#     3762-3764, 2004.  — gain-switched DFB laser: 100 ps pulses, 2.5 MHz.
# [7] Agrawal, G. P., "Fiber-Optic Communication Systems", 5th ed., Wiley,
#     2021, §3.4: Gain-switched semiconductor lasers and pulse generation.


class CWLaser:
    """
    Physics-informed continuous-wave (CW) laser model for QKD simulation.

    Models steady-state output with:
      - Arbitrary polarization state (azimuth + ellipticity)
      - Phase noise via Wiener process (finite linewidth, Henry [1])
      - RIN from the linearized rate equations, including the relaxation
        oscillation resonance (Coldren [2] §5.3)
      - **Optional pulsed mode** (switchable): Gaussian pulse train at
        configurable repetition rate and pulse width (Gobby [6], Agrawal [7])

    The electric field at time t is:
        E(t) = sqrt(P(t)) * exp(j * (omega * t + phi_noise(t))) * E_pol

    In pulsed mode, a normalised Gaussian envelope g(t) is applied so that
    ``mean(|g|^2) = 1`` over the sample window, preserving average power.

    The RIN power spectral density follows the linearized rate-equation
    result (Coldren [2] Eq 5.3.38):

        S_RIN(f) = RIN_0 * (gamma^2 + (2pi f)^2)
                   ------------------------------------------------
                   |(2pi f_RO)^2 - (2pi f)^2 + j*gamma*(2pi f)|^2

    where f_RO is the relaxation oscillation frequency, gamma is the
    damping rate (rad/s), and RIN_0 is the DC RIN density (1/Hz).
    Below f_RO the spectrum is flat; it peaks at f_RO and rolls off
    as 1/f^2 above, in agreement with the classical resonance in the
    laser photon-number response to pump fluctuations (Petermann [5]).

    Parameters
    ----------
    wavelength : float
        Center wavelength in metres.
    power_dbm : float
        Output optical power in dBm.
    linewidth : float
        FWHM optical linewidth in Hz.  Sets the phase diffusion rate:
        D_phi = 2*pi*linewidth  (rad^2/s)  [1, eq. 18].
    rin_density : float
        Relative Intensity Noise spectral density at DC, in dB/Hz.
        Typical values: -140 to -160 dB/Hz for DFB lasers.
        Converts to linear: RIN_lin = 10^(rin_density/10).
    polarization_azimuth : float
        Polarization azimuth angle in radians (psi).
    polarization_ellipticity : float
        Polarization ellipticity angle in radians (chi).
        chi=0  -> linear;  |chi|=pi/4 -> circular.
    relaxation_frequency : float
        Relaxation oscillation frequency f_RO in Hz.  Typical 1-10 GHz
        for DFB lasers at moderate power.  Scales as sqrt(P - P_th)
        (Coldren [2] §5.3.1).  Default 5 GHz.
    damping_rate : float
        Damping rate gamma in rad/s.  Typical (1-5) x 10^10 rad/s for
        DFB lasers.  Related to f_RO via the K factor:
        gamma = K * f_RO^2 + gamma_0  (Coldren [2] §5.3.4).
        Default 1.88e10 rad/s (3 GHz in Hz units).
    pulsed : bool
        If True, output a train of Gaussian pulses instead of CW.
        Default False (backward compatible).
    pulse_width : float
        FWHM of the Gaussian pulse envelope in seconds (default 100e-12).
        Only used when pulsed=True.
    repetition_rate : float
        Pulse repetition frequency in Hz (default 2.5e6).
        Only used when pulsed=True.
    timing_jitter_rms : float
        RMS timing jitter in seconds (default 0.0).
        Each pulse centre is perturbed by N(0, jitter_rms^2).
        Only used when pulsed=True.
    """
    def __init__(
        self,
        wavelength: float,
        power_dbm: float = 0.0,
        linewidth: float = 1e6,
        rin_density: float = -140.0,
        polarization_azimuth: float = 0.0,
        polarization_ellipticity: float = 0.0,
        relaxation_frequency: float = 5e9,
        damping_rate: float = 1.88e10,
        pulsed: bool = False,
        pulse_width: float = 100e-12,
        repetition_rate: float = 2.5e6,
        timing_jitter_rms: float = 0.0,
    ) -> None:
        self.wavelength = wavelength
        self.c = 3e8
        self.h = 6.626e-34
        self.frequency = self.c / self.wavelength
        self.omega = 2 * np.pi * self.frequency

        # Power
        self.power_dbm = power_dbm
        self.power_mw = 10 ** (power_dbm / 10)
        self._power_w = self.power_mw * 1e-3

        # Polarization
        self.polarization_azimuth = polarization_azimuth
        self.polarization_ellipticity = polarization_ellipticity

        # Phase noise (linewidth) — Henry [1], Coldren [2] Ch. 5
        self.linewidth = linewidth
        self._phase_diff_coeff = 2.0 * np.pi * linewidth  # D_phi (rad^2/s)

        # RIN — Coldren [2] §5.3.3, Petermann [5] Ch. 7
        self.rin_density = rin_density
        self._rin_linear = 10.0 ** (rin_density / 10.0)
        self._relaxation_freq = relaxation_frequency
        self._damping_rate = damping_rate

        # RIN correlation time ~ 1/(2*pi*f_RO).  The internal generation
        # uses a coarse timestep when the requested dt is far finer.
        self._rin_dt_min = 1.0 / (10.0 * self._relaxation_freq)

        # Pulsed-mode parameters — Gobby [6], Agrawal [7] §3.4
        self.pulsed = pulsed
        self.pulse_width = pulse_width
        self.repetition_rate = repetition_rate
        self.timing_jitter_rms = timing_jitter_rms

        # Precompute Gaussian sigma from FWHM
        self._pulse_sigma = self.pulse_width / (2.0 * np.sqrt(2.0 * np.log(2.0)))

    @property
    def power_out(self) -> float:
        """Output power in mW (follows SolidStateLaser convention)."""
        return self._power_w * 1e3

    def _polarization_vector(self) -> NDArray[np.complex128]:
        """
        Jones vector for the polarization state.

        From chi (ellipticity) and psi (azimuth):
            E_pol = [cos(chi)*cos(psi) - j*sin(chi)*sin(psi),
                     cos(chi)*sin(psi) + j*sin(chi)*cos(psi)]

        See Yariv [3] Ch. 6.
        """
        chi = self.polarization_ellipticity
        psi = self.polarization_azimuth
        Ex = (np.cos(chi) * np.cos(psi) - 1j * np.sin(chi) * np.sin(psi))
        Ey = (np.cos(chi) * np.sin(psi) + 1j * np.sin(chi) * np.cos(psi))
        return np.array([Ex, Ey])

    def _sample_phase_noise(self, dt: float, n_samples: int) -> NDArray[np.float64]:
        """
        Wiener process phase noise.

        Increments are Gaussian with variance D_phi * dt = 2*pi*linewidth*dt.
        This is the standard model for laser phase diffusion (Henry [1]).
        """
        std_step = np.sqrt(self._phase_diff_coeff * dt)
        increments = np.random.normal(0, std_step, n_samples)
        return np.cumsum(increments)

    def _sample_rin(self, dt: float, n_samples: int) -> NDArray[np.float64]:
        """
        RIN noise with relaxation-oscillation resonance (Coldren [2] Eq 5.3.38).

        The RIN power spectral density from the linearized rate equations is:

            S_RIN(f) = RIN_0 * (gamma^2 + (2pi f)^2)
                       --------------------------------------------------
                       |(2pi*f_RO)^2 - (2pi*f)^2 + j*gamma*(2pi*f)|^2

        where RIN_0 is the DC density (linear), f_RO is the relaxation
        oscillation frequency, and gamma is the damping rate (rad/s).

        Implementation: white noise -> rFFT -> shape by sqrt(S_RIN(f)) ->
        correct PSD scaling -> irFFT.  When the requested dt is smaller
        than the RIN correlation time (~1/f_RO), the noise is generated
        on a coarser internal grid and linearly interpolated to avoid
        numerical issues at extreme sample rates.
        """
        if n_samples < 2:
            return np.zeros(n_samples)

        # RIN varies on timescales ~ 1/f_RO; avoid sub-picosecond sampling.
        internal_dt = dt if dt >= self._rin_dt_min else self._rin_dt_min
        total_time = dt * n_samples
        n_internal = max(2, int(np.ceil(total_time / internal_dt)))

        rin = self._generate_rin(internal_dt, n_internal)

        if n_internal == n_samples and internal_dt == dt:
            return rin
        t_internal = np.arange(n_internal) * internal_dt
        t_requested = np.arange(n_samples) * dt
        return np.interp(t_requested, t_internal, rin)

    def _generate_rin(self, dt: float, n_samples: int) -> NDArray[np.float64]:
        """
        Generate RIN with the relaxation-oscillation spectrum.

        Frequency-domain method:
          1. Unit-variance white noise in time
          2. rFFT to frequency domain
          3. Shape amplitude by sqrt(S_RIN(f)) with correct PSD scaling
          4. irFFT back to time domain
        """
        if n_samples < 2:
            return np.zeros(n_samples)

        fs = 1.0 / dt
        omega_R = 2.0 * np.pi * self._relaxation_freq
        gamma = self._damping_rate

        # RIN PSD shape (Coldren Eq 5.3.38), normalised to 1 at DC
        freqs = np.fft.rfftfreq(n_samples, dt)
        omega = 2.0 * np.pi * freqs
        num = gamma**2 + omega**2
        den = (omega_R**2 - omega**2)**2 + (gamma * omega)**2
        H_sq_raw = np.divide(num, den, out=np.zeros_like(num), where=den > 0)
        H_sq = H_sq_raw / H_sq_raw[0]  # normalise so PSD shape = 1 at DC
        H = np.sqrt(H_sq)

        # White noise in frequency domain
        white = np.random.randn(n_samples)
        W = np.fft.rfft(white)

        # Scale to obtain the correct one-sided PSD.
        #   Y[k] = H[k] · scale · W[k]
        #
        # For k > 0:  G_Y(f_k) = 2·|Y[k]|^2 / (N·fs)  (one-sided, AC)
        # For k = 0:  G_Y(0)   =   |Y[0]|^2 / (N·fs)  (one-sided, DC)
        #
        # E[|W[0]|^2] = N,  E[|W[k]|^2] = N  for k>0.
        #
        # We want G_Y(f) = |H(f)|^2 · RIN_lin.
        #   DC:  E[G_Y(0)] = |H(0)|^2 · scale_dc^2 · N / (N·fs)
        #                    = scale_dc^2 / fs  =  RIN_lin
        #        -> scale_dc = sqrt(RIN_lin · fs)
        #   AC:  E[G_Y(f)]  = |H(f)|^2 · scale_ac^2 · 2·N / (N·fs)
        #                    = |H(f)|^2 · scale_ac^2 · 2/fs  =  |H(f)|^2 · RIN_lin
        #        -> scale_ac = sqrt(RIN_lin · fs / 2)
        scale_ac = np.sqrt(self._rin_linear * fs / 2.0)
        Y = W * H * scale_ac
        Y = Y.astype(np.complex128, copy=False)
        Y[0] *= np.sqrt(2.0)
        if n_samples % 2 == 0:
            Y[-1] *= np.sqrt(2.0)

        rin = np.fft.irfft(Y, n=n_samples)
        return rin

    def _pulse_envelope(self, dt: float, n_samples: int) -> NDArray[np.float64]:
        """Gaussian pulse train envelope normalised so ``mean(|g|^2) = 1``.

        Constructs a train of Gaussian pulses at ``self.repetition_rate``,
        each with sigma = ``self._pulse_sigma``.  Normalisation preserves
        the average optical power set by ``power_dbm``.

        When ``timing_jitter_rms > 0``, each pulse centre is randomly
        perturbed by ``N(0, timing_jitter_rms^2)``, modelling the timing
        jitter of a gain-switched semiconductor laser (Agrawal [7] §3.4).

        Parameters
        ----------
        dt : float
            Sampling interval in seconds.
        n_samples : int
            Number of samples.

        Returns
        -------
        numpy.ndarray, shape (n_samples,)
            Normalized pulse envelope.
        """
        T = 1.0 / self.repetition_rate
        sigma = self._pulse_sigma
        t = np.arange(n_samples, dtype=float) * dt

        half = max(1, int(3.0 * sigma / dt))
        n_pulses = int(np.ceil(n_samples * dt / T))

        envelope = np.zeros(n_samples)
        for k in range(n_pulses):
            t_center = k * T
            if self.timing_jitter_rms > 0:
                t_center += np.random.normal(0, self.timing_jitter_rms)
            center_idx = int(t_center / dt)
            start = max(0, center_idx - half)
            end = min(n_samples, center_idx + half + 1)
            if start >= end:
                continue
            t_local = t[start:end] - t_center
            envelope[start:end] += np.exp(-0.5 * (t_local / sigma) ** 2)

        # Normalise to preserve average power
        mean_sq = np.mean(envelope ** 2)
        if mean_sq > 0:
            envelope /= np.sqrt(mean_sq)
        return envelope

    def sample_field(self, dt: float, n_samples: int) -> NDArray[np.complex128]:
        """
        Primary API — generate field samples with all physical effects.

        Returns the complex envelope (no optical carrier), shape (n_samples, 2)
        for [Ex, Ey], with power, phase noise, RIN, and polarisation.

        When ``pulsed=True`` a normalised Gaussian pulse-train envelope is
        applied (Gobby [6], Agrawal [7]).

        This is the method to use for physics-based simulations involving:
          - Multi-bit sequences at arbitrary baud rate
          - Chromatic dispersion (requires complex envelope, not carrier)
          - PMD
          - Bit-rate-dependent detector response

        For quick single-bit polarisation/phase validation (5 fs window with
        optical carrier), use ``instantaneous_field(over_period=True)`` instead.

        Parameters
        ----------
        dt : float
            Time step between samples in seconds.
        n_samples : int
            Number of samples to generate.

        Returns
        -------
        numpy.ndarray, shape (n_samples, 2)
            Complex envelope field [Ex, Ey] at each sample.
        """
        phi = self._sample_phase_noise(dt, n_samples)
        rin = self._sample_rin(dt, n_samples)
        amp = np.sqrt(np.maximum(self._power_w * (1.0 + rin), 0.0))
        E_pol = self._polarization_vector()
        E = np.outer(amp * np.exp(1j * phi), E_pol)

        if self.pulsed:
            env = self._pulse_envelope(dt, n_samples)
            E *= env[:, np.newaxis]

        return E

    def instantaneous_field(
        self, dt: float = 1e-12, over_period: bool = False,
        n_samples: int = 1000, normalize: bool = True
    ) -> NDArray[np.complex128]:
        """
        Quick validation field — one optical period (~5 fs at 1550 nm).

        Returns the full optical field including the carrier (exp(j·ω·t)) at
        ultra-high temporal resolution.  Designed for single-bit polarisation
        and phase encoding checks only.

        NOT suitable for:
          - High-baud-rate or multi-bit sequences (no modulation bandwidth)
          - Chromatic dispersion (5 fs FFT grid is unphysical)
          - PMD
          - Bit-rate-dependent effects

        For physics-based simulations use ``sample_field(dt, n_samples)``.

        Parameters
        ----------
        dt : float
            Ignored when over_period=True (step derived from period/n_samples).
            Used as phase-noise diffusion step when over_period=False.
        over_period : bool
            If True, return field over one optical period.
        n_samples : int
            Number of samples when over_period=True (default 1000).
        normalize : bool
            If True, return unit-amplitude field (direction only).

        Returns
        -------
        numpy.ndarray
            Shape (2,) if over_period=False, or (n_samples, 2) if over_period=True.
        """
        E_pol = self._polarization_vector()

        if over_period:
            T = 1.0 / self.frequency
            dt = T / n_samples
            t_arr = np.linspace(0, T, n_samples, endpoint=False)

            phi_noise = self._sample_phase_noise(dt, n_samples)
            delta_p = self._sample_rin(dt, n_samples)
            amp = np.sqrt(np.maximum(self._power_w * (1.0 + delta_p), 0))

            carrier = np.exp(1j * (self.omega * t_arr + phi_noise))
            E = np.outer(amp * carrier, E_pol)

            if normalize:
                nrm = np.linalg.norm(E)
                if nrm > 0:
                    E = E / nrm
            return E
        else:
            phi_std = np.sqrt(self._phase_diff_coeff * abs(dt) + 1e-30)
            phi_noise = np.random.normal(0, phi_std)

            delta_p = self._sample_rin(max(dt, 1e-12), 2)[0]
            amp = np.sqrt(max(self._power_w * (1.0 + delta_p), 0))

            E = amp * np.exp(1j * (self.omega * dt + phi_noise)) * E_pol
            if normalize:
                nrm = np.linalg.norm(E)
                if nrm > 0:
                    E = E / nrm
            return E

    def __str__(self) -> str:
        s = (
            f"CWLaser: lambda={self.wavelength:.2e}m, "
            f"P={self.power_mw:.2f}mW, "
            f"dnu={self.linewidth/1e6:.2f}MHz, "
            f"RIN={self.rin_density:.0f}dB/Hz, "
            f"f_RO={self._relaxation_freq/1e9:.2f}GHz, "
            f"psi={self.polarization_azimuth:.4f}rad"
        )
        if self.pulsed:
            s += (
                f", pulsed={self.pulsed}"
                f", width={self.pulse_width:.1e}s"
                f", rate={self.repetition_rate:.1e}Hz"
            )
            if self.timing_jitter_rms > 0:
                s += f", jitter={self.timing_jitter_rms:.1e}s"
        return s
