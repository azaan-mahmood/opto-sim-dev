import numpy as np

# --- LITERATURE SOURCES ---
# [1] Keiser, G., "Optical Fiber Communications", 5th ed., McGraw-Hill, 2015.
#     Ch. 3: Fiber attenuation (dB/km), chromatic dispersion, PMD.
# [2] Hui, R. & O'Sullivan, M., "Fiber-Optic Measurement Techniques",
#     Academic Press, 2009. Material dispersion coefficients for silica.
# [3] Keck, D. B. et al., "Waveguide dispersion in single-mode fibers",
#     IEEE J. Quantum Electron., vol. QE-21, no. 6, 1985.
# [4] Yuan, L. et al., "Stress-induced birefringence and fabrication of
#     in-fiber polarization structures", Opt. Express, vol. 27, 2019.
# [5] Razavi, B., "Design of Integrated Circuits for Optical Communications",
#     2nd ed., Wiley, 2012. PMD and birefringence effects in fiber links.
# [6] Agrawal, G. P., "Fiber-Optic Communication Systems", 5th ed., Wiley, 2021.
#     Ch. 2: Signal degradation in optical fibers.

def cable(fiber_length, E,
                dispersion=False, attenuation_factor=0.25, temperature=0, num_bends=0, pm_dispersion=0.1e-12):
    """
    Transmission through cable, applies dispersion effects on x polarization only

    Parameters:
    cable_type (string): PM, SM
    fiber_length (float): Length of the Fiber in kilometers
    attenuation_factor (float): Known as alpha, expressed in decibels or dB per km. Default 0.182 dB/km for 1550nm
    E (numpy array): X and Y component of the electric field. Accepts as a 2D numpy column array or matrix
    temperature (float): Temperature in Celsius
    num_bends (int): Number of bends in the cable
    pm_dispersion (float): Defaults to 0.1e-12. Is the polarization mode dispersion of the fiber

    Returns:
    numpy array: Electric field (carrying power information)
    """
    # Constants and initial parameters
    L = fiber_length * 1000 # converting kilometers to meters
    T0 = 25  # reference temperature in Celsius

    birefringence_T0 = 0.87e-5  # base birefringence at T0 silica glass (source: the effect of temperature and pressure on the refractive index of some oxide glasses)
    temperature_coefficient = -5e-7  # change in birefringence due temperature coefficient gamma, ThorLabs, pure silica glass
    bend_effect_factor = 2.4e-4  # change in birefringence per bend based on stress (Stress-induced birefringence and fabrication of in-fiber polarization, by Lei Yuan)

    D_material = 17e-12  # ps/(nm*km) for material dispersion at 1550 nm (Chromatic) Hui2009
    D_waveguide = -3e-12 # ps/(nm*km) for waveguide dispersion at 1550 nm. Keck 1985 IEEE
                         # Negative due to self phase modulation dispersion cancellation

    wavelength = 1550e-9  # central wavelength in meters
    pmd_sd = pm_dispersion * np.sqrt(L)
    del_T = pmd_sd**2

    # Calculate total chromatic dispersion
    D_total = D_material + D_waveguide  # ps/(nm*km)

    def apply_dispersion(E, D_total_per_meter, L):
        beta2 = D_total_per_meter * (wavelength ** 2) / (2 * np.pi * 3e8)  # s^2/m, GVD or group velocity dispersion
        f = (1/100)*(1/(D_total_per_meter*L*0.2e-9))
        omega = 2 * np.pi * f
        # omega = 2 * np.pi * 3e8 / wavelength  # angular frequency

        H = np.array([[np.exp(1j * beta2 * L * omega ** 2/2), 0],  # dispersion transfer function L/2 because of two waves
                          [0, 1]])

        E_dispersed = np.transpose(H @ np.transpose(E))

        return E_dispersed

    # Function to calculate birefringence change due to temperature
    def birefringence_temp_change(temp):
        return birefringence_T0 + (temp - T0) * temperature_coefficient

    # Function to calculate birefringence change due to bends
    def birefringence_bend_change(num_bends):
        return num_bends * bend_effect_factor


    # # Function to apply attenuation
    # def apply_attenuation(Ex, Ey, alpha, length):
    #     attenuation = np.exp((alpha * length / 10))
    #     return Ex * attenuation, Ey * attenuation

    def apply_pmd(E, pmd_sd):
        # PMD-induced broadening: model delay as a Gaussian distribution
        # Apply PMD effect to the signal, assuming the distribution is Gaussian at 0 mean and length of Electric field
        Ex = E[:, 0]
        Ey = E[:, 1]
        delay_const = np.random.normal(0, pmd_sd, Ex.shape)
        pmd_matrix = np.exp(1j*delay_const)

        E_pmd_x = pmd_matrix * Ex
        E_pmd = np.column_stack((E_pmd_x, Ey))
        return E_pmd

    birefringence = birefringence_temp_change(temperature) + birefringence_bend_change(num_bends)
    delta_beta = del_T * birefringence / 1.55e-6  # Beat Length for 1550 nm wavelength
    jones = np.array([[np.exp(1j * delta_beta * L / 2), 0],  # Beat Length Formulation
                      [0, 1]])
    # Calculate the Jones matrix for the fiber
    output_signals = np.transpose(jones @ np.transpose(E))
    E = output_signals
    if dispersion is True:
        E = apply_dispersion(output_signals, D_total, L)
        E = apply_pmd(E, pmd_sd)
    else:
        pass
    # Attenuation (Keiser [1] Eq 3.6, Agrawal [6] Eq 2.1.4):
    # Applied directly to the field so that mean(|E_out|²) = att_lin · mean(|E_in|²),
    # keeping power carried within the field (no separate pin/pout).
    att_lin = 10 ** (-attenuation_factor * fiber_length / 10)
    E = E * np.sqrt(att_lin)

    return E


