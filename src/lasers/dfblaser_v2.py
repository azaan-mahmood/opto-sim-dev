"""Split-step time-domain (SS-TDM) dynamic model of a DFB laser diode.

Implements Kim, Chung & Lee, "An Efficient Split-Step Time-Domain Dynamic
Modeling of DFB/DBR Laser Diodes", IEEE J. Quantum Electron. 36(7), 787-794
(2000).  Equation numbers below refer to that paper.

Units are SI converted. The paper's Table I uses CGS system.

Physics summary
---------------


Deviations from the paper, and why
----------------------------------
N/A
"""
import numpy as np
from scipy.linalg import expm

class Laser:
    def __repr__(self):
        pass

    def __init__(self,
                 laser_order:int = 1,                       #mth order laser
                 grating_length:float = 600e-6,             #Length of DFB laser grating
                 n_sections:int = 15,                       #Number of sections DFB is divided into
                 wavelength:float = 1.55e-6,                #Optical Lasing Wavelength
                 i_bias:float = 100e-3,                     #Bias operating current of the Laser
                 run_time = 5e-12                           #Laser observation time or running time
                 ):
        self.grating_length = grating_length
        self.n = n_sections
        self.wavelength = wavelength
        self.i_bias = i_bias
        self.run_time = run_time
        # Parameter List:
        self.w_waveguide = 2e-6
        self.confinement = 0.3
        self.B = 1e-10*1e-6                                 #Spontaneous Recombination, cm^-3/s into m^3/s
        self.C = 0.75e-28*1e-12                             #Auger Carrier cm^6/s into m^6/s
        self.g_N = 2.5e-16*1e-4                                  #Differential Gain
        self.tau = 10e-9                                    #Carrier Lifetime, in seconds
        self.N_0 = 1.8e18*1e6                                   #Transparency Carrier cm^-3
        self.n_eff_0 = 3.283                                #Effective Phase refractive index with injection
        self.n_g = 3.7                                      #Effective Group Refractive Index
        self.d_act = 0.2e-6                                 #Thickness of active layer
        self.A_act = self.w_waveguide * self.d_act          #Area of the Active Layer
        self.beta = 0.5e-4                                  #Spontaneous Coupling Factor
        self.alpha = 5                                      #Linewidth Enhancement Factor
        self.alpha_m = 40*100                                   #Waveguide Loss cm^-1 into m^-1
        self.epsilon = 2e-17*1e-6                                #Nonlinear Saturation Coefficient cm^3 into m^3
        self.c = 3e8                                        #Speed of Light
        self.v_g = self.c/self.n_g                          #Typical group velocity equation based on group refractive index
        self.dz = self.grating_length / self.n
        self.dt = self.dz / self.v_g                        #dz = v_g*dt
        self.total_time_step = max(1, int(self.run_time / self.dt))
        self.Fx = np.zeros(self.n + 1, dtype=complex)
        self.Rx = np.zeros(self.n + 1, dtype=complex)
        self.Fy = np.zeros(self.n + 1, dtype=complex)
        self.Ry = np.zeros(self.n + 1, dtype=complex)
        self.N = np.ones(self.n) * 1e24                     #Carrier density defined AT THE CENTERS of the 15 sections, m^-3
        self.kappa = 50 * 100                               #Coupling Coefficient cm^-1 into m^-1

        #Constants
        self.h = 6.626e-34
        self.nu = self.c/self.wavelength
        self.E_photon = self.h * self.nu
        self.m_order = laser_order
        self.eff_mode = 3.4
        self.ar = 0.01

        #Bragg Condition
        self.bragg_condition = self.m_order*self.wavelength/2*self.eff_mode

        def _simulate_(self):
            self.Fx[:] = 1e-3
            self.Rx[:] = 1e-3

            for _ in range(self.total_time_step):
                Rx_next = np.zeros(self.n + 1, dtype=complex)
                Fx_next = np.zeros(self.n + 1, dtype=complex)
                Fy_next = np.zeros(self.n + 1, dtype=complex)
                Ry_next = np.zeros(self.n + 1, dtype=complex)

                #Compute Center Power by Averaging at Bounderies at each time step
                Px_center = 0.5*((np.abs(self.Fx[:-1])**2 + np.abs(self.Fx[1:])**2) +
                                (np.abs(self.Rx[:-1])**2 + np.abs(self.Rx[1:])**2))
                Py_center = 0.5*((np.abs(self.Fy[:-1])**2 + np.abs(self.Fy[1:])**2) +
                                (np.abs(self.Ry[:-1])**2 + np.abs(self.Ry[1:])**2))
                
                #Computer optical photon densities from optical power
                Sx = Px_center / (self.v_g*self.A_act*self.E_photon)
                Sy = Py_center / (self.v_g*self.A_act*self.E_photon)

                #Get field Gains
                gx = (self.confinement*self.g_N*(self.N - self.N_0))/(2*(1+self.epsilon*Sx) - self.alpha/2)
                gy = (self.confinement*self.g_N*(self.N - self.N_0))/(2*(1+self.epsilon*Sy) - self.alpha/2)

                #Get Change in refractive index
                del_n = -(self.wavelength/4*np.pi)*self.confinement*self.alpha_m*self.g_N*self.N

                #Get detuning
                detuning = 0

                for i in range(self.n):
                    # Fx_in = self.Fx[i]
                    # Rx_in = self.Rx[i+1]
                    # Fy_in = self.Fy[i]
                    # Ry_in = self.Ry[i+1]

                    M_Matrix_x = np.array([
                        [gx[i] - 1j*detuning, -1j*self.kappa],
                        [1j*self.kappa, -(gx[i] - 1j*detuning)]
                    ])

                    M_Matrix_y = np.array([
                        [gy[i] - 1j*detuning, -1j*self.kappa],
                        [1j*self.kappa, -(gy[i] - 1j*detuning)]
                    ])

                    T_mat_x = expm(M_Matrix_x*self.dz)
                    T_mat_y = expm(M_Matrix_y*self.dz)

                    boundary_in_x = np.array([self.Fx[i], self.Rx[i+1]])
                    boundary_out_x = T_mat_x @ boundary_in_x
                    boundary_in_y = np.array([self.Fx[i], self.Ry[i+1]])
                    boundary_out_y = T_mat_y @ boundary_in_y

                    Fx_next[i+1] = boundary_out_x[0]
                    Rx_next[i] =   boundary_out_x[1]

                    Fy_next[i+1] = boundary_out_y[0]
                    Ry_next[i] =   boundary_out_y[1]

                # Apply Facet
                Fx_next[0] = self.ar*Rx_next[0]
                Rx_next[self.n] = self.ar*Fx_next[self.n]

                Fy_next[0] = self.ar*Ry_next[0]
                Ry_next[self.n] = self.ar*Fy_next[self.n]
                
                # Update Carrier densities
                q = 1.602e-19
                V_act = self.A_act*self.dz 
                J = self.i_bias/self.n  #Injection Current Density

                dN_dt_x = J/q*V_act - self.B*self.N**2 - self.C*self.N**3 - (self.v_g*self.g_N*(self.N-self.N_0)*Sx)/(1+self.epsilon*Sx)
                dN_dt_y = J/q*V_act - self.B*self.N**2 - self.C*self.N**3 - (self.v_g*self.g_N*(self.N-self.N_0)*Sy)/(1+self.epsilon*Sy)
                dN_dt = dN_dt_x + dN_dt_y

                self.N +=dN_dt * self.dt
                self.Fx = np.copy(Fx_next)
                self.Rx = np.copy(Rx_next)
                self.Fy = np.copy(Fy_next)
                self.Ry = np.copy(Ry_next)

        return self.Fx, self.Rx, self.Fy, self.Ry
            


