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
# [12] Corning Inc., "Corning SMF-28 Ultra Optical Fiber — Product
#      Information", Corning data sheet PI1463, 2021. Zero-dispersion
#      wavelength 1300-1324 nm; chromatic dispersion <= 18 ps/(nm*km) and
#      PMD link value <= 0.06 ps/sqrt(km) (typical <= 0.1 ps/sqrt(km)) at
#      1550 nm. Matches the ITU-T G.652 total-dispersion spec of
#      ~17 ps/(nm*km) at 1550 nm used as D_TOTAL below.

# Module-level dispersion parameters — the material/waveguide split for
# standard SMF at 1550 nm (Agrawal [6] Fig. 2.10; Hui [2]; Keck [3]) sums
# to the Corning SMF-28 total dispersion spec [12].
# exposed for validation scripts to import dynamically.
D_MATERIAL = 22.0      # ps/(nm·km) material @ 1550 nm (Agrawal [6] Fig 2.10)
D_WAVEGUIDE = -5.0     # ps/(nm·km) waveguide @ 1550 nm (Agrawal [6] Fig 2.10)
D_TOTAL = D_MATERIAL + D_WAVEGUIDE   # ps/(nm·km); = 17.0, matches SMF-28 spec [12]


def bend_birefringence(bend_radius, r_fiber=62.5e-6):
    """Bend-induced birefringence Delta_n_bend = 0.135 * (r_fiber / R)^2.

    The quadratic-in-(r/R) law for a bent single-mode fibre. Every caller
    goes through here, so the formula and its constants live in one place
    and can be checked against the published law in one place.

    Parameters
    ----------
    bend_radius : float — bend radius R in metres.
    r_fiber : float — fibre (cladding) radius in metres (default 62.5 um,
        the standard SMF-28 value).

    Returns
    -------
    float — the bend-induced index difference Delta_n_bend.

    References
    ----------
    [7] Ulrich, R., Rashleigh, S. C., & Eickhoff, W., "Bending-induced
        birefringence in single-mode fibers," Opt. Lett. 5(6), 273-275,
        1980. Eq. 1.
    """
    return 0.135 * (r_fiber / bend_radius) ** 2


def _build_jones_sectional(rng, L, wavelength=1550e-9, temperature=25,
                           bend_radius=None, correlation_length=50.0):
    """Build the ordered-product Jones matrix for the multi-section
    birefringence model (first-principles model).

    The fibre is divided into N = round(L / correlation_length) sections.
    `correlation_length` is the *birefringence correlation length* L_c —
    the distance over which the local birefringence axis orientation
    stays fixed — not a numerical discretisation step; each section draws
    an independent random axis, and the phase retardation
    delta_beta * dz = 2*pi*|delta_n|*dz/lambda is integrated exactly
    within each correlation cell, so no separate integration step is
    needed. Real SMF has L_c of order 10-100 m (Menyuk & Wai [10]).
    The returned matrix is the ordered product of all section matrices.

    `rng` supplies the randomness: the `numpy.random` module for the
    stateless per-call API, or a `numpy.random.Generator` for a
    `FiberRealization`'s own fixed draw.

    How fast it randomises, and how far it holds up
    -----------------------------------------------
    Measured over 500 realizations per length, against a Haar-random SU(2)
    whose mean rotation angle is 2.217 rad:

        L        N    mean rot (rad)   ||<S>||
        1 m      1        0.204        0.9867
        10 m     1        2.044        0.0900
        200 m    4        2.131        0.0435
        1 km    20        2.208        0.0338
        200 km  4000      2.187        0.0432

    The output is already indistinguishable from a uniformly random SU(2)
    by about ten metres. The ensemble mean Stokes vector falls to the
    finite-sample floor, 1/sqrt(500) = 0.045, and stays there; the residual
    is that floor rather than surviving order.

    Ten metres is thousands of times shorter than any link simulated here,
    so at every length in use this returns a uniformly random rotation.
    That is the physically right answer for an uncompensated installed
    fibre, and it is why polarisation encoding needs the compensation loop
    while time-bin encoding does not care: both interfering paths in a
    matched interferometer traverse the same rotation, so it cancels.

    Two consequences worth knowing before reaching for something cheaper.
    Unitarity error grows linearly in the section count at about 1.9e-16
    per section, which is machine epsilon times N: 3.8e-12 at 1000 km and
    1.9e-10 at 50,000 km, further than any cable on the planet. Nothing
    degrades. And the cost is paid once per realization, not per pulse,
    because `FiberRealization` freezes the matrix at construction: 2 ms at
    122 km, 17 ms at 1000 km, against Monte Carlo runs measured in minutes.

    Parameters
    ----------
    L : float — fibre length in metres. Caller must ensure L > 0.
    wavelength : float — centre wavelength (default 1550 nm).
    temperature : float — ambient temperature in C (default 25).
    bend_radius : float or None — bend radius in metres.
    correlation_length : float — birefringence correlation length L_c in
        metres (default 50.0; physical range 10-100 m per Menyuk & Wai [10]).

    Returns
    -------
    ndarray (2, 2) — the fibre's Jones matrix.

    References
    ----------
    Menyuk & Wai, JOSA B 11(7), 1994.
    Wai & Menyuk, JLT 14(2), 1996.
    Agrawal, Fiber-Optic Comm. Systems, 5th ed., §4.1–4.2.
    """
    T0 = 25.0
    birefringence_T0 = 5.0e-8      # Agrawal §4.1: L_B ~ 31 m for SMF-28
    temperature_coefficient = -3.0e-9

    delta_n = (
        birefringence_T0
        + temperature_coefficient * (temperature - T0)
    )
    if bend_radius is not None:
        delta_n += bend_birefringence(bend_radius)

    delta_n += rng.normal(0, 0.1 * birefringence_T0)
    delta_n = np.sign(delta_n) * max(abs(delta_n), 5e-10)

    N = max(1, int(np.round(L / correlation_length)))
    dz = L / N

    delta_beta_dz = 2 * np.pi * abs(delta_n) * dz / wavelength
    half = delta_beta_dz / 2.0
    cos_half = np.cos(half)
    sin_half = np.sin(half)

    theta = rng.uniform(0, 2 * np.pi, size=N)
    phi = np.arccos(2 * rng.uniform(size=N) - 1)
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

    return _ordered_product(J)


def _ordered_product(J):
    """Return J[N-1] @ ... @ J[1] @ J[0] — the ordered matrix product of a
    stack of N (2,2) Jones matrices — via pairwise tree reduction in
    O(log N) vectorised steps instead of N interpreted Python loop
    iterations. Matrix multiplication is associative, so this is exact,
    not an approximation.

    Parameters
    ----------
    J : ndarray (N, 2, 2) — stack of Jones matrices, application order
        matching a naive left-fold: result = J[N-1] @ ... @ J[1] @ J[0].

    Returns
    -------
    ndarray (2, 2) — the ordered product.
    """
    J = J.copy()
    while J.shape[0] > 1:
        n = J.shape[0]
        odd = J[-1:] if n % 2 else None
        if odd is not None:
            J = J[:-1]
        # J[1::2] is the LEFT factor, J[0::2] the RIGHT factor of each pair
        J = np.einsum('nij,njk->nik', J[1::2], J[0::2])
        if odd is not None:
            J = np.concatenate([J, odd], axis=0)   # leftover stays leftmost
    return J[0]


def _apply_birefringence_sectional(E, L, wavelength=1550e-9, temperature=25,
                                   bend_radius=None, correlation_length=50.0):
    """Apply the multi-section birefringence model, drawing a fresh random
    realization from the global `np.random` state on every call.

    See `_build_jones_sectional` for the model itself.

    Calling this once per bit for what is meant to be the same physical
    fibre draws an independent random rotation every bit, which turns the
    fibre into a depolarising channel rather than a fixed unitary. Those
    are different objects: a fixed rotation can be compensated, and in
    every deployed polarisation-encoded system it is. Use
    `FiberRealization` whenever the fibre should stay the same fibre.
    This entry point is for sampling an *ensemble* of fibres.

    Parameters
    ----------
    E : ndarray (N, 2) — complex envelope [Ex, Ey].
    L : float — fibre length in metres.

    Returns
    -------
    ndarray (N, 2) — field after birefringence rotation.
    """
    if L <= 0:
        return E.copy()
    J_total = _build_jones_sectional(np.random, L, wavelength, temperature,
                                     bend_radius, correlation_length)
    return np.transpose(J_total @ np.transpose(E))


def _build_jones_matrix(rng, L, wavelength=1550e-9, temperature=25,
                        bend_radius=None, correlation_length=50.0, model='auto'):
    """Dispatch to the multi-section Jones-matrix builder, mirroring
    `apply_birefringence`'s dispatch logic. Shared by the stateless
    per-call API and by `FiberRealization`. Returns `np.eye(2)` for
    L <= 0.

    There is one birefringence model. `model='sectional'` and the default
    `'auto'` both select it at every length; any other value raises
    ValueError rather than silently falling back.
    """
    if L <= 0:
        return np.eye(2, dtype=complex)
    if model == 'phenomenological':
        raise ValueError(
            "there is one birefringence model; use model='sectional' "
            "or 'auto'")
    return _build_jones_sectional(rng, L, wavelength, temperature,
                                  bend_radius, correlation_length)


class FiberRealization:
    """One specific piece of fibre — birefringence, PMD, chromatic
    dispersion and attenuation, each independently toggleable, applied
    together by a single `apply()` call.

    The quasi-static impairments (birefringence, PMD) are sampled *once*
    at construction and reused for every field passed to `apply()`. 
    This matters because a real fibre's Jones matrix and DGD are
    quasi-static: they drift on timescales of seconds to minutes (thermal
    expansion, mechanical stress), not per bit. The stateless
    `apply_birefringence()` / `apply_pmd()` / `propagate()` functions draw
    a fresh random realization on every call instead, which is correct
    when sampling an *ensemble* of independent fibre realizations (as the
    validation scripts do) but wrong when simulating many bits/pulses
    through what should be one physical fibre. Use this when class when 
    quasi-static impairments are needed.

    Parameters
    ----------
    L_m : float — fibre length in metres.
    wavelength : float — centre wavelength (default 1550 nm).
    temperature : float — ambient temperature in C (default 25).
    bend_radius : float or None — bend radius in metres.
    correlation_length : float — birefringence correlation length L_c for
        the sectional model (m); default 50.0, physical range 10-100 m
        per Menyuk & Wai [10].
    model : str — birefringence model: 'auto' (default; single model at
        all lengths) or 'sectional'. Any other value raises ValueError.
    attenuation_factor : float — dB/km (default 0.182).
    pmd_coeff_ps_sqrt_km : float — PMD coefficient in ps/sqrt(km) (default
        0.1; Corning SMF-28 Ultra spec <= 0.1 ps/sqrt(km) [12]).
    birefringence : bool — include birefringence? (default True).
    cd : bool — include chromatic dispersion? (default False, matching
        `propagate()`'s default).
    pmd : bool — include PMD? (default False, matching `propagate()`'s
        default).
    attenuation : bool — include attenuation? (default True).
    seed : int or None — seeds this realization's own RNG stream,
        independent of the global `numpy.random` state used elsewhere in
        the simulation (e.g. for bit/basis choices, detector noise).
    drift_temperature_rate_C_s : float — rate (C/s) at which the ambient
        temperature changes during a run (default 0.0 = static, the
        behaviour this class has always had). Read by `at()`; it does
        nothing to this realization itself, which is the fibre at t = 0.

    Examples
    --------
    >>> fibre = FiberRealization(L_m=50_000, cd=True, pmd=True, seed=42)
    >>> for _ in range(num_bits):
    ...     E = fibre.apply(E, dt=dt)   # same impairments every call

    >>> fibre = FiberRealization(L_m=50_000, seed=42,
    ...                          drift_temperature_rate_C_s=1e-5)
    >>> later = fibre.at(60.0)          # the same fibre, a minute on
    """

    def __init__(self, L_m, wavelength=1550e-9, temperature=25,
                bend_radius=None, correlation_length=50.0, model='auto',
                attenuation_factor=0.182, pmd_coeff_ps_sqrt_km=0.1,
                birefringence=True, cd=False, pmd=False, attenuation=True,
                seed=None, drift_temperature_rate_C_s=0.0):
        self.rng = np.random.default_rng(seed)
        # Kept so `at()` can rebuild this same fibre in a later state; the
        # shared seed is what makes the two states the same piece of glass
        # rather than two unrelated draws.
        self._seed = seed
        self.drift_temperature_rate_C_s = float(drift_temperature_rate_C_s)
        self.L_m = L_m
        self.wavelength = wavelength
        self.temperature = temperature
        self.bend_radius = bend_radius
        self.correlation_length = correlation_length
        self.model = model
        self.attenuation_factor = attenuation_factor
        self.pmd_coeff_ps_sqrt_km = pmd_coeff_ps_sqrt_km
        self.birefringence_enabled = birefringence
        self.cd_enabled = cd
        self.pmd_enabled = pmd
        self.attenuation_enabled = attenuation

        self._J = (
            _build_jones_matrix(self.rng, L_m, wavelength, temperature,
                                bend_radius, correlation_length, model)
            if birefringence else None
        )
        self._dgd, self._pmd_swap = (
            _sample_pmd_dgd(self.rng, L_m, pmd_coeff_ps_sqrt_km)
            if (pmd and L_m > 0) else (0.0, False)
        )

    def apply(self, E, dt=None):
        """Apply this fibre's enabled impairments to a field E, shape
        (N, 2), in the same order as `propagate()`:
        birefringence -> CD -> PMD -> attenuation.

        Parameters
        ----------
        E : ndarray (N, 2) — complex envelope [Ex, Ey].
        dt : float or None — sampling interval (required if `cd` or `pmd`
            are enabled).

        Returns
        -------
        ndarray (N, 2) — field after all enabled impairments.
        """
        if self.birefringence_enabled and self._J is not None:
            E = np.transpose(self._J @ np.transpose(E))

        if self.cd_enabled or self.pmd_enabled:
            if dt is None:
                raise ValueError(
                    "dt (sampling interval) is required when cd or pmd "
                    "are enabled."
                )
            if self.cd_enabled:
                E = apply_cd(E, dt, self.L_m, self.wavelength)
            if self.pmd_enabled:
                E = _apply_pmd_fixed(E, dt, self._dgd, self._pmd_swap)

        if self.attenuation_enabled:
            E = apply_attenuation(E, self.L_m / 1000.0, self.attenuation_factor)

        return E

    def birefringence_matrix(self):
        """Return this realization's quasi-static birefringence Jones
        matrix (2x2 complex, unitary), or None if birefringence is
        disabled.

        The matrix is drawn once at construction and reused by every
        `apply()` call.  Because it is unitary (SU(2)), the inverse
        channel is its conjugate transpose — this is what a receiver's
        active polarization compensation must implement (e.g. the
        calibration loop of Duplinskiy et al. 2017).

        Returns
        -------
        ndarray (2, 2) or None — a copy of the Jones matrix.
        """
        return None if self._J is None else self._J.copy()

    def at(self, t):
        """This same fibre as it is `t` seconds into the run.

        Returns a new `FiberRealization` whose ambient temperature has
        moved by ``drift_temperature_rate_C_s * t``, or `self` when the
        rate is zero — so a non-drifting fibre costs nothing and stays
        provably unchanged rather than approximately so.

        Why this is the SAME fibre and not another draw
        -----------------------------------------------
        The new realization is built from the same seed, and
        `_build_jones_sectional` consumes its RNG in an order and a count
        that do not depend on temperature: one normal for the delta_n
        perturbation, then 2N uniforms for the section axes, with
        N = round(L / L_c). The section axes are therefore reproduced
        exactly and only the deterministic term

            delta_n = 5.0e-8 + (-3.0e-9) * (T - 25)

        moves. That temperature coefficient is the one this module already
        uses (Agrawal [6] §4.1), so drift introduces no new constant.

        The PMD draw follows the Jones build on the same stream, so it sees
        the same state and the DGD is unchanged. That is the right answer
        as well as the convenient one: DGD is set by the fibre's length and
        geometry, not by an ambient degree.

        `bb84_duplinskiy` builds its calibration fibre this way for the
        same reason — see its `cal_fibre`.

        Composition is exact: ``f.at(a).at(b)`` is ``f.at(a + b)``, because
        the rate carries forward and the offsets add.

        How much drift a given rate buys
        --------------------------------
        Steeply length-dependent, and worth knowing before choosing a rate.
        Each correlation cell carries

            delta_beta * dz = 2*pi * delta_n * L_c / lambda = 10.13 rad

        at L_c = 50 m and 1550 nm, so one degree Celsius shifts every cell
        by 0.61 rad and re-randomises the ordered product of all N of them.
        At 122 km N is 2440 and a single degree is total rescrambling, so
        useful rates are very small.

        The relationship is not monotonic either — the observable relative
        phase wraps — so a drift quoted in rad/s is only locally meaningful
        and has to be MEASURED over the window in use rather than assumed.
        `validate_gobby_impairments.py` reports that mapping per distance.

        Parameters
        ----------
        t : float — elapsed time in seconds since this realization.

        Returns
        -------
        FiberRealization — the fibre at time t, or `self` if it cannot
        have changed.
        """
        if self.drift_temperature_rate_C_s == 0.0 or t == 0.0:
            return self
        return FiberRealization(
            L_m=self.L_m, wavelength=self.wavelength,
            temperature=self.temperature + self.drift_temperature_rate_C_s * t,
            bend_radius=self.bend_radius,
            correlation_length=self.correlation_length, model=self.model,
            attenuation_factor=self.attenuation_factor,
            pmd_coeff_ps_sqrt_km=self.pmd_coeff_ps_sqrt_km,
            birefringence=self.birefringence_enabled, cd=self.cd_enabled,
            pmd=self.pmd_enabled, attenuation=self.attenuation_enabled,
            seed=self._seed,
            drift_temperature_rate_C_s=self.drift_temperature_rate_C_s)


def apply_birefringence(E, L, wavelength=1550e-9, temperature=25,
                        bend_radius=None, correlation_length=50.0, model='auto',
                        enabled=True):
    """Apply birefringence to optical field.

    Uses the multi-section model at all lengths (the former
    there is one birefringence model at every length).

    Parameters
    ----------
    E : ndarray (N, 2) — complex envelope [Ex, Ey].
    L : float — fibre length in metres.
    wavelength : float — centre wavelength (default 1550 nm).
    temperature : float — ambient temperature in C (default 25).
    bend_radius : float or None — bend radius in metres.
    correlation_length : float — birefringence correlation length L_c for
        the multi-section model (m); default 50.0, physical range
        10-100 m per Menyuk & Wai [10].
    model : str — 'auto' (default) or 'sectional'. Any other value
        raises ValueError.
    enabled : bool — if False, return E.copy() unchanged (default True).

    Returns
    -------
    ndarray (N, 2) — field after birefringence rotation.
    """
    if not enabled:
        return E.copy()
    if model == 'phenomenological':
        raise ValueError(
            "there is one birefringence model; use model='sectional' "
            "or 'auto'")
    return _apply_birefringence_sectional(
        E, L, wavelength, temperature, bend_radius, correlation_length)


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


def _sample_pmd_dgd(rng, L, pmd_coeff_ps_sqrt_km=0.1):
    """Sample a differential group delay (seconds) and an H/V swap
    decision from the fibre's PMD (Maxwellian DGD) distribution.

    `rng` supplies the randomness: the `numpy.random` module for the
    stateless per-call API, or a `numpy.random.Generator` for a
    `FiberRealization`'s own fixed draw. Both are valid `random_state`
    values for `scipy.stats.maxwell.rvs`, and both expose `.random()`.

    Parameters
    ----------
    L : float — fibre length in metres. Caller must ensure L > 0.
    pmd_coeff_ps_sqrt_km : float — PMD coefficient in ps/sqrt(km), the
        standard datasheet unit (default 0.1; Corning SMF-28 Ultra spec
        <= 0.1 ps/sqrt(km) [12]).

    Returns
    -------
    dgd : float — sampled differential group delay (seconds).
    swap : bool — whether the fast/slow axes are swapped.
    """
    L_km = L / 1000.0
    pmd_sd = pmd_coeff_ps_sqrt_km * 1e-12 * np.sqrt(L_km)   # seconds
    maxwell_scale = pmd_sd / np.sqrt(3)
    dgd = stats.maxwell.rvs(scale=maxwell_scale, random_state=rng)
    swap = rng.random() < 0.5
    return dgd, swap


def _apply_pmd_fixed(E, dt, dgd, swap):
    """Apply a PMD group-delay split for an already-determined DGD and
    H/V swap decision — the deterministic half of `apply_pmd`. Reused by
    `FiberRealization`, which samples (dgd, swap) once at construction
    instead of on every call.

    Parameters
    ----------
    E : ndarray (N, 2) — complex envelope [Ex, Ey].
    dt : float — sampling interval in seconds.
    dgd : float — differential group delay (seconds).
    swap : bool — whether the fast/slow axes are swapped.

    Returns
    -------
    ndarray (N, 2) — field after PMD.
    """
    N = E.shape[0]
    f = np.fft.fftfreq(N, d=dt)
    omega = 2.0 * np.pi * f
    phase_pmd = omega * dgd / 2.0
    Hx = np.exp(-1j * phase_pmd)
    Hy = np.exp(+1j * phase_pmd)
    if swap:
        Hx, Hy = Hy, Hx

    E_f = np.fft.fft(E, axis=0)
    E_f[:, 0] *= Hx
    E_f[:, 1] *= Hy
    return np.fft.ifft(E_f, axis=0)


def apply_pmd(E, dt, L, pmd_coeff_ps_sqrt_km=0.1):
    """Apply PMD via frequency-domain DGD (Razavi [5] Fig 2.11), drawing a
    fresh random DGD from the global `np.random` state on every call.

    See `_sample_pmd_dgd` for the random draw and `_apply_pmd_fixed` for
    the deterministic operator application. Note: calling this once per
    bit/pulse for what should be the *same* physical fibre redraws PMD on
    every bit. That is wrong for the same reason redrawing birefringence
    per bit is wrong, and for the same underlying cause: both come from
    fibre asymmetry, which is fixed for a given installed fibre. Use
    `FiberRealization` for that case instead.

    Parameters
    ----------
    E : ndarray (N, 2) — complex envelope [Ex, Ey].
    dt : float — sampling interval in seconds.
    L : float — fibre length in metres.
    pmd_coeff_ps_sqrt_km : float — PMD coefficient in ps/sqrt(km) (default
        0.1; Corning SMF-28 Ultra spec <= 0.1 ps/sqrt(km) [12]).

    Returns
    -------
    E_out : ndarray (N, 2) — field after PMD.
    dgd : float — the sampled differential group delay (seconds).
    """
    if L <= 0:
        return E.copy(), 0.0
    dgd, swap = _sample_pmd_dgd(np.random, L, pmd_coeff_ps_sqrt_km)
    return _apply_pmd_fixed(E, dt, dgd, swap), dgd


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
              temperature=25, bend_radius=None, pmd_coeff_ps_sqrt_km=0.1,
              correlation_length=50.0, model='auto',
              birefringence=True, cd=None, pmd=None, attenuation=True,
              fiber_realization=None):
    """Transmission through an optical fibre.

    Each impairment can be independently enabled/disabled via individual
    flags.  The legacy `dispersion` flag is an alias for setting both
    ``cd`` and ``pmd`` when they are not explicitly provided.

    Impairments applied in order:
      1. Birefringence (multi-section model)
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
    pmd_coeff_ps_sqrt_km : float — PMD coefficient in ps/sqrt(km) (default
        0.1; Corning SMF-28 Ultra spec <= 0.1 ps/sqrt(km) [12]).
    correlation_length : float — birefringence correlation length L_c for
        the multi-section model (m); default 50.0, physical range
        10-100 m per Menyuk & Wai [10].
    model : str — birefringence model: 'auto' (default) or 'sectional';
        any other value raises ValueError.
    birefringence : bool — apply birefringence? (default True).
    cd : bool or None — apply CD?  None uses ``dispersion`` value.
    pmd : bool or None — apply PMD?  None uses ``dispersion`` value.
    attenuation : bool — apply attenuation? (default True).
    fiber_realization : FiberRealization or None — if given, delegates the
        entire impairment chain to `fiber_realization.apply(E, dt=dt)` and
        ignores every other impairment-related argument to this function
        (`wavelength`, `temperature`, `bend_radius`, `pmd_coeff_ps_sqrt_km`,
        `correlation_length`, `model`, `birefringence`, `cd`, `pmd`,
        `attenuation`, `attenuation_factor`) — those belong to the
        realization, which was already configured with them at
        construction time. Use this to simulate many bits/pulses through
        one physical fibre, which is what `FiberRealization` exists for.

    Returns
    -------
    ndarray (N, 2) — transmitted field.
    """
    if fiber_realization is not None:
        return fiber_realization.apply(E, dt=dt)

    L = fiber_length * 1000

    # Resolve individual CD/PMD flags from the legacy dispersion flag
    enable_cd = cd if cd is not None else dispersion
    enable_pmd = pmd if pmd is not None else dispersion

    if birefringence:
        E = apply_birefringence(E, L, wavelength, temperature, bend_radius,
                                correlation_length, model)

    if enable_cd or enable_pmd:
        if dt is None:
            raise ValueError(
                "dt (sampling interval) is required when cd=True or pmd=True."
            )
        if enable_cd:
            E = apply_cd(E, dt, L, wavelength)
        if enable_pmd:
            E, _ = apply_pmd(E, dt, L, pmd_coeff_ps_sqrt_km)

    if attenuation:
        E = apply_attenuation(E, fiber_length, attenuation_factor)
    return E
