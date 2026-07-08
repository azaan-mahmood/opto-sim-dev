import numpy as np
from scipy import stats

# --- LITERATURE SOURCES ---
# [1] Keiser, G., "Optical Fiber Communications", 5th ed., McGraw-Hill, 2015.
#     Ch. 3: Fiber attenuation (dB/km), chromatic dispersion, PMD.
# [2] Hui, R. & O'Sullivan, M., "Fiber-Optic Measurement Techniques",
#     Academic Press, 2009. Material dispersion coefficients for silica.
# [3] Keck, D. B. et al., "Waveguide dispersion in single-mode fibers",
#     IEEE J. Quantum Electron., vol. QE-21, no. 6, 1985.
# [4] Yuan, L. et al., "Stress-induced birefringence and fabrication of
#     in-fiber polarization devices by controlled femtosecond laser
#     irradiations", Opt. Express, vol. 24, no. 2, pp. 1062-1071, 2016.
#     DOI: 10.1364/oe.24.001062.
# [5] Razavi, B., "Design of Integrated Circuits for Optical Communications",
#     2nd ed., Wiley, 2012. PMD and birefringence effects in fiber links.
# [6] Agrawal, G. P., "Fiber-Optic Communication Systems", 5th ed., Wiley, 2021.
#     Ch. 2: Signal degradation in optical fibers (CD, PMD).
#     Ch. 4: Fiber birefringence and beat length.

# Module-level dispersion parameters (Hui [2], Keck [3]) —
# exposed for validation scripts to import dynamically.
D_MATERIAL = 17.0      # ps/(nm·km) material @ 1550 nm
D_WAVEGUIDE = -3.0     # ps/(nm·km) waveguide @ 1550 nm
D_TOTAL = D_MATERIAL + D_WAVEGUIDE   # ps/(nm·km)

# Validation: stores the DGD sampled in each cable() call when
# dispersion=True.  Cleared by the validation script before its loop.
_dgd_sampled = []


def cable(fiber_length, E, dt=None, wavelength=1550e-9,
          dispersion=False, attenuation_factor=0.182,
          temperature=25, num_bends=0, pm_dispersion=0.1e-12):
    """
    Transmission through an optical fiber.

    Applies, in order:
      1. Birefringence (beat-length Jones matrix)  — polarisation rotation
      2. Chromatic dispersion (FFT-based, Agrawal [6] §2.4)
      3. PMD (frequency-domain DGD)
      4. Attenuation (Keiser [1] Eq 3.6)

    Parameters
    ----------
    fiber_length : float
        Fiber length in kilometres.
    E : ndarray, shape (N, 2)
        Complex-envelope optical field [Ex, Ey] at N time samples.
        **Must be the complex envelope** (no optical carrier) when
        ``dispersion=True`` — use ``laser.sample_field()`` or
        ``laser.sample_field()`` or ``instantaneous_field(over_period=False)``.
    dt : float or None
        Sampling interval in seconds.  **Required when dispersion=True**
        so that the FFT frequency grid can be constructed.
    wavelength : float
        Centre wavelength in metres.  Default 1550e-9 (1550 nm).
    dispersion : bool
        Apply chromatic dispersion and PMD?  Default False — only
        birefringence and attenuation are applied.
    attenuation_factor : float
        Power attenuation in dB/km (Keiser [1] Eq 3.6).
        Default 0.182 dB/km (standard SMF-28 at 1550 nm).
    temperature : float
        Ambient temperature in °C.  Affects birefringence.
        Default 25 °C.
    num_bends : int
        Number of discrete bends.  Each bend contributes to the
        stress-induced birefringence (Yuan [4]).
        Default 0.
    pm_dispersion : float
        PMD coefficient in s/sqrt(m).  Default 0.1e-12.
        The RMS DGD after length L is ``pm_dispersion * sqrt(L)``.

    Returns
    -------
    ndarray, shape (N, 2)
        Transmitted field.

    References
    ----------
    ...as listed in the module header.
    """
    L = fiber_length * 1000                # km → m
    T0 = 25.0                              # reference temperature °C

    # --- Physical constants ---
    c0 = 2.99792458e8                      # speed of light m/s
    omega0 = 2 * np.pi * c0 / wavelength   # carrier angular frequency

    # --- Birefringence (beat-length model) ---
    # Δn = n_slow - n_fast, dimensionless.
    # Base birefringence at T0 for silica glass (source: effect of
    # temperature and pressure on the refractive index of some oxide
    # glasses).  Temperature coefficient from ThorLabs, pure silica.
    # Bend contribution from Yuan [4].
    birefringence_T0 = 0.87e-5
    temperature_coefficient = -5e-7
    bend_effect_factor = 2.4e-4

    birefringence = (
        birefringence_T0
        + temperature_coefficient * (temperature - T0)
        + bend_effect_factor * num_bends
    )

    # Propagation-constant difference (Agrawal [6] Eq 4.1.2):
    #   Δβ = β_slow - β_fast = 2π · Δn / λ
    dbeta = 2.0 * np.pi * birefringence / wavelength   # rad/m
    # Symmetric Jones matrix: Ex advances by Δβ·L/2, Ey retards by Δβ·L/2.
    # Relative phase = Δβ·L = 2π·L·Δn/λ = 2π·L/L_B.
    jones = np.array([[np.exp(1j * dbeta * L / 2), 0],
                      [0, np.exp(-1j * dbeta * L / 2)]])
    E = np.transpose(jones @ np.transpose(E))

    # --- Chromatic dispersion (Agrawal [6] §2.4.1) ---
    if dispersion:
        if dt is None:
            raise ValueError(
                "dt (sampling interval) is required when dispersion=True."
            )

        # Total dispersion parameter (module-level constants)
        D_SI = D_TOTAL * 1e-6   # ps/(nm·km) → s/m²

        # GVD parameter (Agrawal [6] Eq 2.4.11)
        #   β₂ = -D · λ² / (2πc)   [s²/m]
        beta2 = -D_SI * wavelength**2 / (2 * np.pi * c0)

        # Frequency grid (baseband; E is the complex envelope)
        N = E.shape[0]
        f = np.fft.fftfreq(N, d=dt)        # Hz, FFT-native order
        omega = 2.0 * np.pi * f            # rad/s, baseband angular freq

        # Dispersion transfer function (Agrawal [6] Eq 2.4.11):
        #   H(ω) = exp(-j · β₂ · ω² · L / 2)
        # ω here is the baseband frequency (deviation from ω₀).
        H = np.exp(-1j * beta2 * omega**2 * L / 2)

        # Apply in frequency domain to both polarisations (CD is isotropic)
        E_f = np.fft.fft(E, axis=0)
        E_f = E_f * H[:, np.newaxis]
        E = np.fft.ifft(E_f, axis=0)

        # --- PMD (Razavi [5] Fig 2.11, Agrawal [6] §4.5) ---
        # DGD follows a Maxwellian distribution (3D PMD vector magnitude).
        # scale a = pmd_sd / √3 gives RMS(DGD) = pmd_sd and
        # mean(DGD) = 2·a·√(2/π) ≈ 0.921·pmd_sd.
        pmd_sd = pm_dispersion * np.sqrt(L)          # seconds, RMS DGD
        maxwell_scale = pmd_sd / np.sqrt(3)
        dgd = stats.maxwell.rvs(scale=maxwell_scale)  # Maxwellian DGD
        _dgd_sampled.append(dgd)                     # record for validation

        # Frequency-dependent Jones matrix for PMD:
        #   J_pmd = diag(exp(-j·ω·Δτ/2), exp(+j·ω·Δτ/2))
        # The fast/slow axis orientation is random per realization.
        phase_pmd = omega * dgd / 2.0
        Hx = np.exp(-1j * phase_pmd)
        Hy = np.exp(+1j * phase_pmd)
        if np.random.rand() < 0.5:
            Hx, Hy = Hy, Hx

        E_f_after_cd = np.fft.fft(E, axis=0)    # re-FFT after CD + IFFT
        E_f_after_cd[:, 0] *= Hx
        E_f_after_cd[:, 1] *= Hy
        E = np.fft.ifft(E_f_after_cd, axis=0)

    # --- Attenuation (Keiser [1] Eq 3.6, Agrawal [6] Eq 2.1.4) ---
    # Applied directly to the field so that mean(|E_out|²) keeps the
    # power in the field — no separate pin/pout tracking.
    att_lin = 10.0 ** (-attenuation_factor * fiber_length / 10.0)
    E = E * np.sqrt(att_lin)

    return E
