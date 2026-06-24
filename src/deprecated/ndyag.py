import numpy as np
from scipy.integrate import solve_ivp

# --- LITERATURE SOURCES ---
# 1. Saleh, B. E. A., & Teich, M. C. "Fundamentals of Photonics", 4th-Level Systems, Chapter 15.
# 2. Koechner, W. "Solid-State Laser Engineering", Nd:YAG Spectroscopic Properties.
# 3. Siegman, A. E. "Lasers", Rate Equations for 4-Level Systems.

class NdYAGLaser:
    """
    Nd:YAG (1064nm) 4-Level Laser Model.
    
    In a 4-level system, the lower laser level (N1) is significantly above the ground state (N0)
    and decays almost instantaneously. Thus, population inversion is achieved easily (N2 > N1 ~ 0).
    """
    def __init__(self, wavelength=1064e-9, power_dbm=0.0, noise_std=1e3):
        # --- Spectroscopic Constants (Koechner / Saleh & Teich) ---
        self.tau2 = 230e-6    # Upper state lifetime (s)
        self.tau1 = 30e-9     # Lower state lifetime (s) - Very fast, often treated as 0
        self.tau_c = 10e-9    # Cavity photon lifetime (s)
        
        # Ion Concentration (1% doping in YAG crystal)
        # N0 ~ 1.38e20 cm^-3 = 1.38e26 m^-3
        self.N0 = 1.38e26 
        
        # Stimulated Emission Cross Section (m^2)
        self.sigma = 2.8e-23 
        
        self.wavelength = wavelength
        self.c = 3e8
        self.n = 1.82         # Refractive index of YAG crystal
        self.vg = self.c / self.n
        
        # Simulation parameters
        self.alpha = 0.05     # Cavity loss coefficient
        self.beta = 1e-6      # Spontaneous emission coupling factor
        self.noise_std = noise_std
        
        # Initial Conditions
        self.I_0 = 10**(power_dbm/10) * 1e-3 / (6.626e-34 * (self.c/self.wavelength) * self.vg)
        
        # Pump Rate (R_p)
        # For Nd:YAG, the threshold pump rate R_th = 1 / (sigma * vg * tau2 * tau_c)
        # We set it to 5x threshold for clear dynamics.
        self.R_th = 1.0 / (self.sigma * self.vg * self.tau2 * self.tau_c)
        self.Rp = 5.0 * self.R_th # Normalized Pump Rate

    def rate(self, t, y):
        """
        4-Level Rate Equations.
        y[0] = N2 (Upper level population)
        y[1] = I  (Photon density)
        """
        N2, I = y
        
        # Stimulated Transition Rate
        # W_st = sigma * N2 * vg * I
        W_st = self.sigma * N2 * self.vg * I
        
        # 1. Upper Level Population (N2)
        # Grows by pump, decays by spontaneous emission and stimulated emission.
        # Saleh & Teich Eq. 15.1-4
        dN2_dt = self.Rp - (N2 / self.tau2) - W_st
        
        # 2. Photon Density (I)
        # Grows by stimulated emission and spontaneous seed, decays by cavity loss.
        # Added Langevin-style noise for research-grade realism.
        
        # Consistent pseudo-random noise for solver stability
        seed = int(t * 1e11) % 2**32
        np.random.seed(seed)
        noise_term = np.random.normal(0, self.noise_std * np.sqrt(max(0, I)))
        
        dI_dt = W_st - (I / self.tau_c) + self.beta * (N2 / self.tau2) + noise_term
        
        return [dN2_dt, dI_dt]

    def out_pow(self, t_span=[0, 1e-3]):
        """Solve and return output power."""
        y0 = [0, self.I_0]
        t_eval = np.linspace(t_span[0], t_span[1], 1000)
        
        sol = solve_ivp(
            self.rate, t_span, y0, 
            t_eval=t_eval, method='BDF', rtol=1e-6, atol=1e-9
        )
        
        # Power = I * h * nu * vg * CoreArea (Assuming 1mm^2 for Nd:YAG rod)
        nu = self.c / self.wavelength
        area = 1e-6 # 1mm^2
        power = sol.y[1] * (6.626e-34 * nu) * self.vg * area
        
        return power, sol
