"""Single-Photon Avalanche Diode (SPAD) — Geiger-mode detector model.

Models gated single-photon detection with dead time, dark counts, and
afterpulsing.  Inherits physical constants and responsivity from the
linear APD class; overrides the detection method to return binary
click/no-click output.

References
----------
[1] ID Quantique, "ID230 InGaAs SPAD datasheet", 2020.
    Dead time 13 us, efficiency 10%, DCR 15 Hz, afterpulsing 5%.
[2] Saleh, B. E. A. & Teich, M. C., "Fundamentals of Photonics",
    3rd ed., Wiley, 2019.  Ch. 17: Photon detection statistics.
[3] Yuan, Z. L. et al., "High speed single photon detection in the
    near infrared", Appl. Phys. Lett. 91(4), 041114, 2007.
"""
import numpy as np
from src.detectors.apd import apd


class spad(apd):
    """Geiger-mode single-photon avalanche diode (SPAD).

    Inherits wavelength, quantum efficiency, and physical constants from
    the linear APD class.  Adds dead time, dark count rate, afterpulsing,
    and gated detection.

    Parameters
    ----------
    wavelength : float — centre wavelength (m).
    quantum_efficiency : float — detection efficiency eta (default 0.10).
    dead_time : float — dead time after each click (s, default 13e-6).
    dark_count_rate : float — dark counts per second (Hz, default 15).
    afterpulse_prob : float — probability of afterpulse after each click
                         (default 0.05).
    gate_width : float — detection gate window (s, default 20e-9).
    temperature : float — detector temperature (K, default 298).
    """

    def __init__(self, wavelength, quantum_efficiency=0.10,
                 dead_time=13e-6, dark_count_rate=15.0,
                 afterpulse_prob=0.05, gate_width=20e-9,
                 temperature=298):
        # Linear APD parameters not used in Geiger mode — pass defaults
        super().__init__(
            wavelength=wavelength,
            excess_noise_factor=1.0,
            load_resistance=50.0,
            temperature=temperature,
            gain=1,
            quantum_efficiency=quantum_efficiency,
            dark_current=0.0,
        )

        self.dead_time = dead_time
        self.dcr = dark_count_rate          # Hz
        self.afterpulse_prob = afterpulse_prob
        self.gate_width = gate_width

        # State
        self._armed = True
        self._last_click_time = -np.inf
        self._afterpulse_pending = False
        self._afterpulse_time = -np.inf

    # ------------------------------------------------------------------
    # State management
    # ------------------------------------------------------------------
    def reset(self):
        """Reset detector state for a new pulse train."""
        self._armed = True
        self._last_click_time = -np.inf
        self._afterpulse_pending = False
        self._afterpulse_time = -np.inf

    @property
    def is_armed(self):
        return self._armed

    # ------------------------------------------------------------------
    # Core detection
    # ------------------------------------------------------------------
    def detect(self, power, t):
        """Process a single detection gate.

        Parameters
        ----------
        power : float — mean optical power during the gate window (W).
        t : float — centre time of this gate (s).

        Returns
        -------
        int — 1 (click) or 0 (no click).
        """
        # --- Enforce dead time ------------------------------------------------
        if not self._armed:
            if t - self._last_click_time >= self.dead_time:
                self._armed = True
                # Clear any afterpulse scheduled but never fired during
                # this dead period.  Otherwise it stays pending with a
                # stale timestamp and can fire on the first gate of a
                # *later* dead period.
                self._afterpulse_pending = False
                self._afterpulse_time = -np.inf
            else:
                # Still dead — check for scheduled afterpulse
                if self._afterpulse_pending and t >= self._afterpulse_time:
                    self._afterpulse_pending = False
                    self._armed = False
                    self._last_click_time = t
                    self._schedule_afterpulse(t)
                    return 1
                return 0

        click = 0

        # --- Dark count (Poisson within gate) --------------------------------
        p_dark = self.dcr * self.gate_width
        if np.random.random() < p_dark:
            click = 1

        # --- Signal photon detection -----------------------------------------
        if not click and power > 0:
            photon_energy = self.h * self.frequency
            # Expected photons in gate: mu = P * t_gate / (h*nu)
            mu = power * self.gate_width / photon_energy
            # Each photon is detected INDEPENDENTLY with probability eta, so
            # the detected count is Poisson(eta*mu) and
            #
            #     P(click) = 1 - exp(-eta*mu)
            #
            # Not `qe * (1 - exp(-mu))`, which reads "at least one photon
            # arrives, THEN one coin flip at eta".  That undercounts: if two
            # photons arrive the chance of detecting at least one is
            # 1 - (1-eta)^2, not eta.  The two forms agree as mu -> 0 and
            # diverge with intensity -- -0.05 % at mu = 0.001, -20.4 % at
            # mu = 0.5, -54.8 % at mu = 2.0 -- so the difference hides at
            # single-photon levels and matters as soon as mu is raised.
            #
            # It was found as the residual disagreement between
            # `validate_gobby.signal_click_prob()` and the Monte Carlo:
            # predicted 0.9625 against a measured 0.9533 +/- 0.0123, 0.75
            # sigma.
            #
            # One RNG draw either way, so pulse-for-pulse stream alignment
            # is unchanged and any difference is physical.
            p_click = 1.0 - np.exp(-self.qe * mu)
            if np.random.random() < p_click:
                click = 1

        # --- Update state on click -------------------------------------------
        if click:
            self._armed = False
            self._last_click_time = t
            self._schedule_afterpulse(t)

        return click

    def _schedule_afterpulse(self, t):
        """With probability afterpulse_prob, schedule a false click.

        Always assigns both fields so a failed roll clears any
        previously pending afterpulse rather than leaving its stale state
        in place.
        """
        if np.random.random() < self.afterpulse_prob:
            self._afterpulse_pending = True
            delay = np.random.exponential(self.dead_time * 0.5)
            self._afterpulse_time = t + max(delay, self.gate_width)
        else:
            self._afterpulse_pending = False
            self._afterpulse_time = -np.inf

    # ------------------------------------------------------------------
    # Convenience: detect a full pulse train
    # ------------------------------------------------------------------
    def detect_pulse_train(self, powers, times):
        """Process a sequence of gates.

        Parameters
        ----------
        powers : array-like — power in each gate (W).
        times : array-like — centre time of each gate (s).

        Returns
        -------
        ndarray of int — click (1) / no-click (0) for each gate.
        """
        powers = np.asarray(powers, dtype=float)
        times = np.asarray(times, dtype=float)
        clicks = np.zeros(len(powers), dtype=int)
        for i in range(len(powers)):
            clicks[i] = self.detect(powers[i], times[i])
        return clicks
