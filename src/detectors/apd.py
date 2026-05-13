import numpy as np


class apd:
    def __init__(self,wavelength,
                 excess_noise_factor, load_resistance,
                 temperature, gain=12, frequency=3e8/1550e-09, quantum_efficiency=0.9, dark_current=10e-6):
        self.qe = quantum_efficiency  # Quantum efficiency (η)
        self.gain = gain  # APD gain (M)
        self.charge = 1.6e-19   # Electron Charge (e)
        self.enf = excess_noise_factor  # Excess noise factor (F)
        self.RL = load_resistance  # Load resistance (ohms)
        self.T = temperature  # Temperature (K)
        self.q = 1.602e-19  # Elementary charge (C)
        self.kB = 1.38e-23  # Boltzmann's constant (J/K)
        self.wavelength = wavelength    # Wavelength
        self.epsilon = 8.854e-12    # Permittivity of Free Space
        self.c = 3e8
        self.h = 6.626e-34  # Planck's constant (J·s)
        self.frequency = frequency
        self.R = self.qe*(self.charge*self.wavelength)/(self.h*self.c) #AW^-1
        self.dark_current = dark_current
        self.dcr = self.dark_current/self.charge

    def detect_photons(self, field_energy, area, exposure_time=1e-9):
        # Calculate the energy of a single photon
        # energy_per_photon = self.h * self.c / self.wavelength

        # Calculate the number of photons based on the input power and energy per photon
        # photon_flux = power / energy_per_photon  # Photon flux (photons per second) Fundamental of Photonics Chapter 13
        mpn = field_energy/self.h * self.frequency
        # Calculate the expected number of photons detected in the given area over the exposure time
        expected_photons = mpn * area * exposure_time * self.qe  # Photons detected over exposure time

        # Check if the expected value is too large for Poisson; use Normal approximation if necessary
        if expected_photons > 1e6:  # Thresholds
            detected_photons = int(np.random.normal(expected_photons, np.sqrt(expected_photons)))
        else:
            detected_photons = np.random.poisson(expected_photons)

        # Remove negative distribution
        detected_photons = max(detected_photons, 0)

        return detected_photons

    def calculate_output_current(self, power, frequency):
        # Calculate the primary photo-current
        # Optoelectronics and Photonics 2nd ed Kasap
        I_ph = self.R * power

        # Amplified signal current
        I_signal = self.gain * I_ph

        return I_signal

    def calculate_noise(self, I_signal, bandwidth):
        #Page 426 of Optoelectronics and Photonics 2nd ed Kasap
        # Shot noise
        shot_noise = np.sqrt(2 * self.charge * self.dark_current * bandwidth)

        # Quantum noise
        quantum_noise = np.sqrt(2 * self.charge * I_signal * bandwidth)

        # Thermal noise (Johnson-Nyquist noise)
        thermal_noise = np.sqrt(4 * self.kB * self.T * bandwidth / self.RL)

        # Total noise current
        I_noise = np.sqrt(shot_noise ** 2 + thermal_noise ** 2 + quantum_noise**2)*self.enf

        return I_noise

    def output(self, E, area, bandwidth, details=False):
        field_energy = 0.5 * self.epsilon * np.real(E.max())      # Energy of an electric field
        detected_photons = self.detect_photons(field_energy, area)   # Detect Photons over an Area
        energy_per_photon = self.h * self.c / self.wavelength   # Calculate typical energy of a photon
        power = detected_photons * energy_per_photon    # Detected Photons times the energy of a single photon
        I_signal = self.calculate_output_current(power, self.frequency)     # Output signal current value
        I_noise = self.calculate_noise(I_signal, bandwidth)     # Output signal current value along with noise

        SNR = I_signal / I_noise  # Signal-to-noise ratio
        if details:
            return {
                'Detected Photons': detected_photons,
                'DCR': self.darkcount,
                'I_signal': I_signal,
                'noise_current': I_noise,
                'SNR': SNR
            }
        else:
            return detected_photons
