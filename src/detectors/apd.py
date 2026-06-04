import numpy as np

# --- LITERATURE SOURCES ---
# [1] Kasap, S. O., "Optoelectronics and Photonics", 2nd ed., Pearson, 2013.
#     Ch. 4: APD gain, responsivity, shot/thermal noise (Eq. 4.19, 4.23, 4.42-4.46).
# [2] Agrawal, G. P., "Fiber-Optic Communication Systems", 5th ed., Wiley, 2021.
#     Ch. 4: Optical receiver noise, dark current, SNR (Eq. 4.1.2-4.4.3).
# [3] Saleh, B. E. A. & Teich, M. C., "Fundamentals of Photonics", 3rd ed., Wiley, 2019.
#     Ch. 17: Photon detection statistics, Poisson model (Eq. 17.1-10).


class apd:
    def __init__(self, wavelength, excess_noise_factor, load_resistance,
                 temperature, gain=12, quantum_efficiency=0.9, dark_current=10e-6):
        # Optical frequency is derived from wavelength — not a user-tunable knob
        self.c = 3e8           # Speed of light in vacuum (m/s)
        self.h = 6.626e-34     # Planck's constant (J·s)
        self.kB = 1.38e-23     # Boltzmann's constant (J/K)
        self.charge = 1.602e-19  # Elementary charge (C)
        self.epsilon = 8.854e-12  # Permittivity of free space (F/m)

        self.wavelength = wavelength
        self.frequency = self.c / self.wavelength   # Optical frequency (Hz)
        self.qe = quantum_efficiency    # Quantum efficiency (η)
        self.gain = gain                # APD gain (M)
        self.enf = excess_noise_factor  # Excess noise factor (F)
        self.RL = load_resistance       # Load resistance (Ω)
        self.T = temperature            # Temperature (K)

        # Responsivity: R = η·e·λ/(h·c)  [A/W]  (Kasap [1] Eq. 4.19)
        self.R = self.qe * self.charge * self.wavelength / (self.h * self.c)

        self.dark_current = dark_current  # Dark current (A)
        self.dcr = self.dark_current / self.charge  # Dark count rate (Hz)

    def detect_photons(self, power, exposure_time):
        """
        Poisson-distributed photon detection from incident optical power.

        photon_rate = P / (h·ν)          (Agrawal [2] Eq. 4.1.2)
        expected    = photon_rate · t · η  (Saleh & Teich [3] Eq. 17.1-10)

        Gaussian approximation is used when expected > 1e6.
        """
        photon_energy = self.h * self.frequency
        if photon_energy <= 0:
            return 0

        expected = (power / photon_energy) * exposure_time * self.qe

        if expected > 1e6:
            detected = int(np.random.normal(expected, np.sqrt(expected)))
        else:
            detected = np.random.poisson(expected)

        return max(detected, 0)

    def calculate_output_current(self, power):
        """
        APD signal current: I_signal = M · R · P  (Kasap [1] Eq. 4.23)
        """
        return self.gain * self.R * power

    def calculate_noise(self, I_signal, bandwidth):
        """
        RMS noise current of the APD receiver.

        Noise sources (Kasap [1] Eq. 4.42–4.46, Agrawal [2] Eq. 4.4.3):

          i_d²  = 2·e·I_dark·B       — dark current shot noise
          i_q²  = 2·e·I_signal·B     — quantum (signal) shot noise
          i_th² = 4·k·T·B / R_L      — Johnson-Nyquist thermal noise

        The excess noise factor F applies only to the shot-noise terms:
          i_total² = F · (i_d² + i_q²) + i_th²   (Kasap [1] Eq. 4.45)
        """
        shot_dark_sq = 2 * self.charge * self.dark_current * bandwidth
        shot_signal_sq = 2 * self.charge * I_signal * bandwidth
        thermal_sq = 4 * self.kB * self.T * bandwidth / self.RL

        return np.sqrt(self.enf * (shot_dark_sq + shot_signal_sq) + thermal_sq)

    def output(self, E, bandwidth, area=1, exposure_time=None, details=False):
        """
        Compute the APD output with realistic noise.

        The electric field carries the optical power (field convention:
        mean(|E|²) = optical power in Watts). The power flows through:
          1. Poisson photon detection (discrete statistics)
          2. Responsivity + gain → signal current
          3. Additive Gaussian noise (shot + thermal)

        Parameters
        ----------
        E : ndarray
            Electric field time samples (1D array for one PBS arm).
            Must be pre-calibrated so that mean(|E|²) = power in Watts.
        bandwidth : float
            Receiver electrical bandwidth in Hz.
        area : float
            Detector area in m² (reserved for future use).
        exposure_time : float or None
            Integration time in seconds. If None, uses 1/(2·B) (Nyquist).
        details : bool
            If True, return a dict with all intermediate values.

        Returns
        -------
        float or dict
            details=False: noisy output current I_total (A).
            details=True:  dict with photon counts, currents, SNR.
        """
        if exposure_time is None:
            exposure_time = 1.0 / (2 * bandwidth)

        # Optical power derived from the field (Watts convention)
        power = np.mean(np.abs(E)**2)

        detected_photons = self.detect_photons(power, exposure_time)
        I_signal = self.calculate_output_current(power)
        I_noise = self.calculate_noise(I_signal, bandwidth)
        I_total = np.random.normal(I_signal, max(I_noise, 1e-30))
        SNR = I_signal / I_noise if I_noise > 0 else float('inf')

        if details:
            return {
                'Detected Photons': detected_photons,
                'DCR': self.dcr,
                'I_signal': I_signal,
                'noise_current': I_noise,
                'I_total': I_total,
                'SNR': SNR,
            }
        return I_total
