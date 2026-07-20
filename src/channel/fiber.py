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
# [7] Ulrich, R. et al., "Bending-induced birefringence in single-mode
#     fibers", Opt. Lett., vol. 5, no. 6, pp. 273-275, 1980.
#     DOI: 10.1364/ol.5.000273.
# [8] Smith, A. M., "Bend-induced birefringence in single-mode optical
#     fibers", Appl. Opt., vol. 19, no. 15, pp. 2606-2611, 1980.
#     DOI: 10.1364/ao.19.002606.
# [9] Shibata, N. et al., "Bend-induced birefringence of a single-mode
#     fiber evaluated by a heterodyne method", J. Opt. Soc. Am. A,
#     vol. 3, no. 11, pp. 1935-1939, 1986.  DOI: 10.1364/josaa.3.001935.
#     Confirms Ulrich's model down to R = 2 mm.
# [10] Menyuk, C. R. & Wai, P. K. A., "Random birefringence and polarization
#      mode dispersion", JOSA B, vol. 11, no. 7, pp. 1288-1299, 1994.
# [11] Wai, P. K. A. & Menyuk, C. R., "Anisotropic nonlinear pulse
#      propagation in randomly birefringent fibers", JLT, vol. 14, no. 2,
#      pp. 148-157, 1996.

# Module-level dispersion parameters (Hui [2], Keck [3]) —
# exposed for validation scripts to import dynamically.
D_MATERIAL = 17.0      # ps/(nm·km) material @ 1550 nm
D_WAVEGUIDE = -3.0     # ps/(nm·km) waveguide @ 1550 nm
D_TOTAL = D_MATERIAL + D_WAVEGUIDE   # ps/(nm·km)

# Birefringence model threshold: fibres shorter than this (in metres) use
# the multi-section model; longer fibres use the phenomenological model.
SECTIONAL_LIMIT = 2000  # 2 km


def _random_su2_rotation(angle):
    """Random SU(2) rotation matrix for a given rotation angle (radians).

    Generates a uniformly random axis on the Poincaré sphere and rotates
    by `angle` around it.  angle = 0 -> identity, angle ~ pi -> fully mixed.

    References
    ----------
    Menyuk & Wai, JOSA B 11(7), 1994.
    Wai & Menyuk, JLT 14(2), 1996.
    """
    theta = np.random.uniform(0, 2 * np.pi)
    phi = np.arccos(2 * np.random.uniform() - 1)
    n = np.array([np.sin(phi) * np.cos(theta),
                  np.sin(phi) * np.sin(theta),
                  np.cos(phi)])
    sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
    sigma_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)
    n_sigma = n[0] * sigma_x + n[1] * sigma_y + n[2] * sigma_z
    half = angle / 2.0
    U = np.cos(half) * np.eye(2, dtype=complex) - 1j * np.sin(half) * n_sigma
    det = np.linalg.det(U)
    return U / np.sqrt(det)


def _apply_birefringence_sectional(E, L, wavelength=1550e-9, temperature=25,
                                   bend_radius=None, section_length=1.0):
    """Multi-section birefringence via ordered product of random-axis
    Jones matrices (first-principles model).

    The fibre is divided into N = round(L / section_length) sections.
    Each section has an independent random birefringence axis and
    a deterministic phase retardation delta_beta * dz = 2*pi*|delta_n|*dz/lambda.
    The total Jones matrix is the ordered product of all section matrices.

    Suitable for short fibres (< SECTIONAL_LIMIT), DV-QKD, and DPS QKD
    where the detailed polarisation evolution is resolved.

    Parameters
    ----------
    E : ndarray (N, 2) — complex envelope [Ex, Ey].
    L : float — fibre length in metres.
    wavelength : float — centre wavelength (default 1550 nm).
    temperature : float — ambient temperature in C (default 25).
    bend_radius : float or None — bend radius in metres.
    section_length : float — section length in metres (default 1.0).

    Returns
    -------
    ndarray (N, 2) — field after birefringence rotation.

    References
    ----------
    Menyuk & Wai, JOSA B 11(7), 1994.
    Wai & Menyuk, JLT 14(2), 1996.
    Agrawal, Fiber-Optic Comm. Systems, 5th ed., §4.1–4.2.
    """
    if L <= 0:
        return E.copy()

    T0 = 25.0
    r_fiber = 62.5e-6
    birefringence_T0 = 5.0e-8      # Agrawal §4.1: L_B ~ 31 m for SMF-28
    temperature_coefficient = -3.0e-9
    bend_effect_factor = 0.135

    delta_n = (
        birefringence_T0
        + temperature_coefficient * (temperature - T0)
    )
    if bend_radius is not None:
        delta_n += bend_effect_factor * (r_fiber / bend_radius) ** 2

    delta_n += np.random.normal(0, 0.1 * birefringence_T0)
    delta_n = np.sign(delta_n) * max(abs(delta_n), 5e-10)

    N = max(1, int(np.round(L / section_length)))
    dz = L / N

    delta_beta_dz = 2 * np.pi * abs(delta_n) * dz / wavelength
    half = delta_beta_dz / 2.0
    cos_half = np.cos(half)
    sin_half = np.sin(half)

    theta = np.random.uniform(0, 2 * np.pi, size=N)
    phi = np.arccos(2 * np.random.uniform(size=N) - 1)
    n0 = np.sin(phi) * np.cos(theta)
    n1 = np.sin(phi) * np.sin(theta)
    n2 = np.cos(phi)

    J = np.empty((N, 2, 2), dtype=complex)
    J[:, 0, 0] = cos_half - 1j * sin_half * n2
    J[:, 0, 1] = -1j * sin_half * (n0 - 1j * n1)
    J[:, 1, 0] = -1j * sin_half * (n0 + 1j * n1)
    J[:, 1, 1] = cos_half + 1j * sin_half * n2

    det = J[:, 0, 0] * J[:, 1, 1] - J[:, 0, 1] * J[:, 1, 0]
    J /= np.sqrt(det)[:, np.newaxis, np.newaxis]

    J_total = np.eye(2, dtype=complex)
    for k in range(N):
        J_total = J[k] @ J_total

    return np.transpose(J_total @ np.transpose(E))


def _apply_birefringence_phenomenological(E, L, wavelength=1550e-9,
                                           temperature=25, bend_radius=None):
    """Phenomenological birefringence via single SU(2) rotation with
    sqrt(L)-dependent angle (Menyuk & Wai 1994).

    The rotation angle is theta = min(pi, sqrt(L / L_char) * pi / 2)
    where L_char = L_0 * (delta_n_0 / |delta_n|)^2, L_0 = 75 km.

    Suitable for long fibres where the multi-section model's rapid
    scrambling converges to uniform SU(2).

    Parameters
    ----------
    E : ndarray (N, 2) — complex envelope [Ex, Ey].
    L : float — fibre length in metres.
    wavelength : float — centre wavelength (default 1550 nm).
    temperature : float — ambient temperature in C (default 25).
    bend_radius : float or None — bend radius in metres.

    Returns
    -------
    ndarray (N, 2) — field after birefringence rotation.

    References
    ----------
    Menyuk & Wai, JOSA B 11(7), 1994.
    Wai & Menyuk, JLT 14(2), 1996.
    Ulrich, Opt. Lett. 5(6), 1980.
    """
    if L <= 0:
        return E.copy()

    T0 = 25.0
    r_fiber = 62.5e-6
    delta_n_0 = 0.87e-5            # base birefringence (phenomenological)
    temperature_coefficient = -5e-7
    bend_effect_factor = 0.135
    L0 = 75e3                       # characteristic length at delta_n_0 (75 km)

    delta_n = (
        delta_n_0
        + temperature_coefficient * (temperature - T0)
    )
    if bend_radius is not None:
        delta_n += bend_effect_factor * (r_fiber / bend_radius) ** 2

    delta_n = np.sign(delta_n) * max(abs(delta_n), 1e-10)

    L_char = L0 * (delta_n_0 / abs(delta_n)) ** 2
    theta = min(np.pi, np.sqrt(L / L_char) * np.pi / 2)

    return np.transpose(_random_su2_rotation(theta) @ np.transpose(E))


def apply_birefringence(E, L, wavelength=1550e-9, temperature=25,
                        bend_radius=None, section_length=1.0, model='auto',
                        enabled=True):
    """Apply birefringence to optical field.

    Dispatches to the multi-section or phenomenological model based on
    fibre length and the `model` parameter.

    Parameters
    ----------
    E : ndarray (N, 2) — complex envelope [Ex, Ey].
    L : float — fibre length in metres.
    wavelength : float — centre wavelength (default 1550 nm).
    temperature : float — ambient temperature in C (default 25).
    bend_radius : float or None — bend radius in metres.
    section_length : float — section length for multi-section model (m).
    model : str — 'auto' (default, dispatch by L), 'sectional', or
           'phenomenological'.
    enabled : bool — if False, return E.copy() unchanged (default True).

    Returns
    -------
    ndarray (N, 2) — field after birefringence rotation.
    """
    if not enabled:
        return E.copy()
    if model == 'sectional':
        return _apply_birefringence_sectional(
            E, L, wavelength, temperature, bend_radius, section_length)
    elif model == 'phenomenological':
        return _apply_birefringence_phenomenological(
            E, L, wavelength, temperature, bend_radius)
    else:  # 'auto'
        if L < SECTIONAL_LIMIT:
            return _apply_birefringence_sectional(
                E, L, wavelength, temperature, bend_radius, section_length)
        else:
            return _apply_birefringence_phenomenological(
                E, L, wavelength, temperature, bend_radius)


def apply_cd(E, dt, L, wavelength=1550e-9):
    """Apply chromatic dispersion via FFT (Agrawal [6] Eq 2.4.11).

    Parameters
    ----------
    E : ndarray (N, 2) — complex envelope [Ex, Ey].
    dt : float — sampling interval in seconds.
    L : float — fibre length in metres.
    wavelength : float — centre wavelength (default 1550 nm).

    Returns
    -------
    ndarray (N, 2) — CD-dispersed field.
    """
    c0 = 2.99792458e8
    D_SI = D_TOTAL * 1e-6
    beta2 = -D_SI * wavelength**2 / (2 * np.pi * c0)

    N = E.shape[0]
    f = np.fft.fftfreq(N, d=dt)
    omega = 2.0 * np.pi * f
    H = np.exp(-1j * beta2 * omega**2 * L / 2)

    E_f = np.fft.fft(E, axis=0)
    E_f = E_f * H[:, np.newaxis]
    return np.fft.ifft(E_f, axis=0)


def apply_pmd(E, dt, L, pm_dispersion=0.1e-12):
    """Apply PMD via frequency-domain DGD (Razavi [5] Fig 2.11).

    Parameters
    ----------
    E : ndarray (N, 2) — complex envelope [Ex, Ey].
    dt : float — sampling interval in seconds.
    L : float — fibre length in metres.
    pm_dispersion : float — PMD coefficient in s/sqrt(m) (default 0.1e-12).

    Returns
    -------
    E_out : ndarray (N, 2) — field after PMD.
    dgd : float — the sampled differential group delay (seconds).
    """
    if L <= 0:
        return E.copy(), 0.0
    pmd_sd = pm_dispersion * np.sqrt(L)
    maxwell_scale = pmd_sd / np.sqrt(3)
    dgd = stats.maxwell.rvs(scale=maxwell_scale)

    N = E.shape[0]
    f = np.fft.fftfreq(N, d=dt)
    omega = 2.0 * np.pi * f
    phase_pmd = omega * dgd / 2.0
    Hx = np.exp(-1j * phase_pmd)
    Hy = np.exp(+1j * phase_pmd)
    if np.random.rand() < 0.5:
        Hx, Hy = Hy, Hx

    E_f = np.fft.fft(E, axis=0)
    E_f[:, 0] *= Hx
    E_f[:, 1] *= Hy
    return np.fft.ifft(E_f, axis=0), dgd


def apply_attenuation(E, L_km, attenuation_factor=0.182):
    """Apply power attenuation (Keiser [1] Eq 3.6).

    Parameters
    ----------
    E : ndarray (N, 2) — complex envelope [Ex, Ey].
    L_km : float — fibre length in kilometres.
    attenuation_factor : float — dB/km (default 0.182, SMF-28 @ 1550 nm).

    Returns
    -------
    ndarray (N, 2) — attenuated field.
    """
    att_lin = 10.0 ** (-attenuation_factor * L_km / 10.0)
    return E * np.sqrt(att_lin)


def propagate(fiber_length, E, dt=None, wavelength=1550e-9,
              dispersion=False, attenuation_factor=0.182,
              temperature=25, bend_radius=None, pm_dispersion=0.1e-12,
              section_length=1.0, model='auto',
              birefringence=True, cd=None, pmd=None, attenuation=True):
    """Transmission through an optical fibre.

    Each impairment can be independently enabled/disabled via individual
    flags.  The legacy `dispersion` flag is an alias for setting both
    ``cd`` and ``pmd`` when they are not explicitly provided.

    Impairments applied in order:
      1. Birefringence (sectional or phenomenological)
      2. Chromatic dispersion (FFT-based, Agrawal [6] §2.4)
      3. PMD (frequency-domain DGD)
      4. Attenuation (Keiser [1] Eq 3.6)

    Parameters
    ----------
    fiber_length : float — fibre length in kilometres.
    E : ndarray (N, 2) — complex-envelope field.
    dt : float or None — sampling interval (required when cd=True or pmd=True).
    wavelength : float — centre wavelength (default 1550 nm).
    dispersion : bool — legacy alias for cd + pmd (default False).
    attenuation_factor : float — dB/km (default 0.182).
    temperature : float — C (default 25).
    bend_radius : float or None — m (default None).
    pm_dispersion : float — s/sqrt(m) (default 0.1e-12).
    section_length : float — birefringence section length in m (default 1.0).
    model : str — birefringence model: 'auto', 'sectional', 'phenomenological'.
    birefringence : bool — apply birefringence? (default True).
    cd : bool or None — apply CD?  None uses ``dispersion`` value.
    pmd : bool or None — apply PMD?  None uses ``dispersion`` value.
    attenuation : bool — apply attenuation? (default True).

    Returns
    -------
    ndarray (N, 2) — transmitted field.
    """
    L = fiber_length * 1000

    # Resolve individual CD/PMD flags from the legacy dispersion flag
    enable_cd = cd if cd is not None else dispersion
    enable_pmd = pmd if pmd is not None else dispersion

    if birefringence:
        E = apply_birefringence(E, L, wavelength, temperature, bend_radius,
                                section_length, model)

    if enable_cd or enable_pmd:
        if dt is None:
            raise ValueError(
                "dt (sampling interval) is required when cd=True or pmd=True."
            )
        if enable_cd:
            E = apply_cd(E, dt, L, wavelength)
        if enable_pmd:
            E, _ = apply_pmd(E, dt, L, pm_dispersion)

    if attenuation:
        E = apply_attenuation(E, fiber_length, attenuation_factor)
    return E
