import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

#https://www.thorlabs.com/newgrouppage9.cfm?objectgroup_id=14199
# E is a column matrix, use column stack E = np.column_stack((Ex, Ey))
# Modeling and simulation of pulsed Er:YAG laser system
# Spectroscopic properties of Er:Yb by Lincong Rao


class SolidStateLaser:
    def __init__(self, wavelength, polarization_azimuth, polarization_ellipticity=None, frequency=None, power_dbm=1.0, noise_std=0.0):
        self.tau1 = 1.071e-3 # Lifetime of lower energy level (Handbook Series on Semiconductor Parameters), also called radiative life time
        self.tau2 = 6.381e-3  # Lifetime of upper energy level
        self.tau_c = 14.4e-6 # Photon Decay or also called as Fluorescence Lifetime or Photon Lifetime
        self.N0 = 8.59e20 # Total number of atoms of Er+3:Yb+3
        self.sigma12 = 0.98e-20  # Absorption cross-section
        self.sigma21 = 1.01e-20  # Emission cross-section
        self.g2byg1 = 1 # Degeneracy Rate.
        self.Rp = self.g2byg1*(1/self.tau2)  # Pump rate (Minimum required for population inversion)
        
        self.alpha = 0.1  # Gain/Loss coefficient (Lowered to allow for visible dynamics)
        self.beta = 1e-4  # Spontaneous emission coupling factor (Physically required for laser startup)
        self.Gamma = 1.0  # Confinement factor, confinement of photons in the laser
        self.c = 3e8  # Speed of light in vacuum (m/s)
        self.h = 6.626e-34  # Planck's constant (J·s)
        self.noise_std = noise_std # Langevin noise scaling factor

        self.wavelength = wavelength
        self.polarization_azimuth = polarization_azimuth
        if frequency is None:
            self.frequency = self.c / self.wavelength
        else:
            self.frequency = frequency
        if polarization_ellipticity is None:
            self.polarization_ellipticity = 0
        else:
            self.polarization_ellipticity = polarization_ellipticity
        self.power_dbm = power_dbm
        self.power_mw = 10**(power_dbm/10)
        # --- MODIFIED: 2026-05-07 14:49:41 ---
        # Bug Fix: I_0 must use the optical photon frequency (c/lambda), not the
        # user-supplied RF/modulation frequency. Using RF frequency inflated I_0 by ~1e8x.
        self.optical_frequency = self.c / self.wavelength
        self.I_0 = self.power_mw / (self.h * self.optical_frequency)
        # --- END MODIFIED ---
        self.power_out, self.final_photons = self.out_pow()

    def out_pow(self):
        # --- MODIFIED: 2026-05-07 14:49:41 ---
        # Bug Fix: N1_0 was hardcoded to 2e23, which is 232x larger than N0 (total atoms).
        # Physically, at t=0 all atoms are in the ground state, so N1_0 = N0.
        N1_0 = self.N0  # All atoms start in lower energy level (ground state)
        # --- END MODIFIED ---
        N2_0 = 0  # Initial population of upper energy level
        y0 = [N1_0, N2_0, self.I_0]
        # --- MODIFIED: 2026-05-07 14:49:41 ---
        # Bug Fix: Simulation window extended from 1us to 20ms.
        # tau1=1.07ms, tau2=6.38ms — the 1us window was far shorter than any
        # relevant atomic timescale, causing no observable laser dynamics.
        t_span = [0, 20e-3]
        t_eval = np.linspace(0, 20e-3, 2000)
        # --- END MODIFIED ---
        # sol = solve_ivp(self.rate, t_span, y0, t_eval=t_eval, method='RK45')
        sol = solve_ivp(
            self.rate,
            t_span,
            y0,
            t_eval=t_eval,
            method='BDF',
            rtol=1e-6,
            atol=1e-9
        )
        # Integrate photon density over time
        integrated_I = abs(np.trapezoid(sol.y[2], sol.t))
        average_I = integrated_I / (t_span[1] - t_span[0])
        # --- MODIFIED: 2026-05-12 08:45:31 ---
        # Source: Saleh & Teich, Fundamentals of Photonics, Eq. 14.0-1
        # Bug Fix: Use optical frequency (c/lambda) for photon energy, not user RF frequency.
        power_out = average_I * self.h * self.optical_frequency
        # --- END MODIFIED ---
        return power_out, abs(sol.y[2][-1])

    def rate(self, t, y):
        # --- MODIFIED: 2026-05-12 08:45:31 ---
        # Literature-Grounded Quasi-3-Level Model
        # Sources:
        # [1] Bjarklev, A., "Optical Fiber Amplifiers", Sec 3.2 (Atom Conservation)
        # [2] Desurvire, E., "Erbium-Doped Fiber Amplifiers", Eq 1.4.1-1.4.3
        # [3] Saleh & Teich, "Fundamentals of Photonics", Eq 15.1-7 & 15.1-8
        
        N1_dynamic, N2, I = y
        n_fiber = 1.45  # Group index: Desurvire Eq 1.3.6
        vg = self.c / n_fiber # Group velocity
        
        # Enforce Atom Conservation: N1 + N2 = N0 (Bjarklev Sec 3.2)
        # We derive N1 from N0 to ensure physical consistency
        N1 = self.N0 - N2
        
        # Stimulated Transition Rate (Emission - Absorption)
        # Saleh & Teich Eq 15.1-7, Bjarklev Sec 3.4
        # Note: we use sigma21 for emission and sigma12 for absorption
        W_st = (self.sigma21 * N2 - self.sigma12 * N1) * vg * I
        
        # Population Dynamics: Desurvire Eq 1.4.2
        # N2 grows by pump, decays by spontaneous emission and stimulated transitions
        dN2_dt = self.Rp * N1 - N2/self.tau2 - W_st
        
        # Photon Dynamics: Saleh & Teich Eq 15.1-8
        # I grows by stimulated transitions and decays by cavity loss.
        # Added a spontaneous emission term (beta * N2 / tau2) to allow laser startup.
        
        # Langevin-style noise term (Stochastic fluctuation)
        # To ensure numerical stability with adaptive-step solvers like BDF,
        # we use a pseudo-random seed based on 't' so the derivative is consistent.
        noise_term = 0
        if self.noise_std > 0:
            # Consistent pseudo-random value for time t
            seed = int(t * 1e9) % 2**32
            np.random.seed(seed)
            noise_term = np.random.normal(0, self.noise_std * np.sqrt(max(0, I)))
            
        dI_dt = self.Gamma * W_st - self.alpha * I / self.tau_c + self.beta * (N2 / self.tau2) + noise_term
        
        # dN1_dt is strictly the opposite of dN2_dt due to conservation
        return [-dN2_dt, dN2_dt, dI_dt]
        # --- END MODIFIED ---

    def get_electric_field(self, t=0, over_period=False, normalize = True):
        """
        Get the electric field vector at time t or over one period of the oscillation.

        Parameters:
        t (float): Time in seconds. Ignored if over_period  is True.
        over_period (bool): If True, return the electric field vectors over one period.
        normalize (bool): If True, returns electric field as a norm, otherwise as product of
                          E = E*sqrt(power)

        Returns:
        numpy array: If over_period is False, return the electric field vector [Ex, Ey, Ez] at time t.
                     If over_period is True, return a tuple (t, Ex, Ey) where t is an array of time points,
                     Ex and Ey are arrays of the electric field components over one period.
                     Returns a normalized E Matrix which only contains direction and not amplitude
        """
        E0 = np.sqrt(2 * self.power_out)  # amplitude of the electric field
        omega = 2 * np.pi * self.frequency
        phi = self.polarization_azimuth
        chi = self.polarization_ellipticity

        if over_period:
            t = np.linspace(0, 2 * np.pi / self.frequency, 1000)
            E = np.array([self._calculate_electric_field(ti, E0, omega, phi, chi) for ti in t])
            if normalize:
                E_normalized = E / np.linalg.norm(E)
                return E_normalized
            elif not normalize:
                return E
        else:
            return self._calculate_electric_field(t, E0, omega, phi, chi)

    def _calculate_electric_field(self, t, E0, omega, phi, chi):
        """
        Calculate the electric field vector at time t.

        Parameters:
        t (float): Time in seconds.
        E0 (float): Amplitude of the electric field.
        omega (float): Angular frequency of the electric field.
        phi (float): Polarization azimuth angle.
        chi (float): Polarization ellipticity angle.

        Returns:
        numpy array: Electric field vector [Ex, Ey].
        """
        # From Gerd Keiser Book (Optical Fiber Communications)
        # Electric Field Ex (Generally called slow axis or ordinary ray).
        # Electric Field Ey (Generally called fast axis or extraordinary ray).
        # Assume Ez = 0 because of its orthogonal propagation
        # Not using exp(), was causing logical errors!
        # Might figure out later #ADD TO_DO!
        Ex = E0 * (np.cos(omega * t + phi) + 1j*np.sin(omega*t + phi))
        Ey = E0 * (np.cos(omega * t + chi) + 1j*np.sin(omega*t + chi))
        return np.array([Ex, Ey])

    def plot_photon_density(self, t, I):
        """Plot the photon density over time."""
        plt.figure(figsize=(8, 6))
        plt.plot(t, I, label='Photon Density (I)')
        plt.xlabel('Time (s)')
        plt.ylabel('Photon Density')
        plt.title('Photon Density vs Time')
        plt.grid(True)
        plt.legend()
        plt.show()

    def __str__(self):
        return (f"Polarized Light Source: λ={self.wavelength:.2e} m,"
                f" φ={self.polarization_azimuth:.4f} rad, f={self.frequency:.2e} Hz, "
                f"Pdbm={10 * np.log10(self.power_out):.2f} dBm, "
                f"Pout={self.power_out:.3f} mW, "
                f"Photon Density={self.I_0}")