import numpy as np

# --- LITERATURE SOURCES ---
# [1] Zeilinger, A., "General properties of lossless beam splitters in
#     interferometry", Am. J. Phys. 49(9), pp. 882-883, 1981.
#     The general 2x2 unitary for a lossless four-port, expressed via the
#     Pauli matrices.  Source of the coupler scattering matrices below.
# [2] Collett, E., "Field Guide to Polarization", SPIE Press, 2005.
#     Jones calculus; ideal polarisers and polarising beam splitters as
#     the diagonal operators diag(1,0) and diag(0,1).
# [3] Saleh, B. E. A. & Teich, M. C., "Fundamentals of Photonics",
#     3rd ed., Wiley, 2019.  Beam splitters, couplers, wave retarders.
#
# COUPLER PHASE CONVENTION -- read before adding a coupler.
# Zeilinger [1] shows the lossless 2x2 unitary is
#     U = exp(i*phi) * [[r, i*t], [i*conj(t), conj(r)]],  |r|^2 + |t|^2 = 1,
# and that the "imaginary" form [[t, i*r], [i*r, t]] and the "real" form
# [[t, r], [r, -t]] are both unitary, differing only in where the phase
# reference sits.  They are physically equivalent *only if one convention
# is used consistently*.
#
# THIS MODULE USES THE REAL FORM, matching `AsymmetricMZI`.  Mixing the
# two is not a cosmetic difference: with the imaginary form the two-beam
# interference term goes as sin(delta_phi) rather than cos(delta_phi), so
# a protocol that encodes in {0, pi} -- as BB84 phase encoding does --
# lands exactly on the zeros and the fringe vanishes entirely.  That bug
# produced a flat ~50% QBER and is documented in GOBBY-2 (section 19) of
# opto-sim-issues-and-fixes.md.
#
# RETARDER PHASE CONVENTION -- read before adding a wave plate.
# `halfwave` and `quarterwave` use the standard Jones forms [2],
#     J_HWP(t) = [[cos 2t, sin 2t], [sin 2t, -cos 2t]],
# with NO absolute-phase prefactor.  Some sources write these with a
# leading `j` (or exp(-i*pi/2), exp(-i*pi/4)); that is the absolute-phase
# variant and it multiplies the whole matrix, carrying no physical content.
#
# Both functions used to carry such a prefactor.  It made `halfwave` return
# a *purely imaginary* field for a real input (measured |Im|/|Re| = 1.6e16)
# -- invisible in |E|^2 for a single path, but a real relative phase as
# soon as two interfering paths pass through different numbers of
# retarders.  Removed; see GOBBY-4 (section 21) in
# opto-sim-issues-and-fixes.md.
#
# What is NOT a global phase: the relative `i` between the fast and slow
# axes of `quarterwave`.  That is the retardance -- the entire point of the
# component -- and it stays.


def coupler_split(power, E, ratio=0.5):
    """Split one input across the two output ports of a 2x2 coupler.

    Amplitude splitting for a lossless coupler [1]: with a single
    populated input, the outputs are sqrt(ratio) and sqrt(1 - ratio)
    times the input field, so the *powers* divide as ratio : (1 - ratio)
    and the total is conserved.  Real convention (see the module note).

    Parameters
    ----------
    power : float — input power (bookkeeping; the returned powers are
        `power * ratio` and `power * (1 - ratio)`).
    E : ndarray (N, 2) — complex-envelope field.
    ratio : float in [0, 1] — fraction of the power sent to the port arm
        (default 0.5, an ideal 3 dB coupler).

    Returns
    -------
    (port_power, port_E, tap_power, tap_E)

    Notes
    -----
    An earlier version returned the *unscaled* input as `port_E` and a
    Hadamard-mixed copy as `tap_E`, so `ratio` affected only the returned
    powers while the fields ignored it entirely — the two were mutually
    inconsistent.  Any caller relying on the old field behaviour was
    relying on a bug.
    """
    if not 0 <= ratio <= 1:
        raise Exception("Incorrect Ratio")
    port_power = power * ratio
    tap_power = power * (1 - ratio)
    port_E = E * np.sqrt(ratio)
    tap_E = E * np.sqrt(1.0 - ratio)
    return port_power, port_E, tap_power, tap_E


def coupler_combine(power_port, port_E, power_tap, tap_E, out_ports=1):
    """Ideal 3 dB coupler combine, REAL convention.

        E_out1 = (E1 + E2)/sqrt(2)
        E_out2 = (E1 - E2)/sqrt(2)

    Unitary, so total power is conserved, and the interference term goes as
    cos(delta_phi).

    CONVENTION (GOBBY-7d).  This previously returned the imaginary form,
    `(E1 + 1j*E2)/sqrt(2)` and `(1j*E1 + E2)/sqrt(2)`.  Both forms are
    legitimate and unitary -- see the COUPLER PHASE CONVENTION note in the
    module header -- but this module and `AsymmetricMZI` use the real one,
    and mixing them is not cosmetic: with the imaginary form the
    interference goes as sin(delta_phi), so a protocol encoding in {0, pi}
    lands on the zeros and the fringe vanishes.  That produced a flat ~50%
    QBER once already (GOBBY-2 §19).

    It was left inconsistent through GOBBY-3 to GOBBY-7b as a known,
    recorded deferral (§21.3) on the grounds that a legitimate-but-
    inconsistent convention is a different category from a prefactor with
    no physical content.  It is aligned here because the module now claims
    the real form in its own header, and a function contradicting its
    module's stated convention is a trap for the next reader.

    **No production call site exists** -- this function is reachable only
    through `src.channel.__init__` and its own tests, which is why the
    change is safe to make in one step.  `AsymmetricMZI` implements its own
    combiner (`E_c = r*E_s + s*E_l`, `E_d = r*E_s - s*E_l`) and is
    untouched by this; the two now agree.

    Parameters
    ----------
    power_port, power_tap : float — input powers (bookkeeping only; the
        returned powers are always derived from the output fields).
    port_E, tap_E : ndarray (N, 2) — input complex-envelope fields.
    out_ports : int — 1 returns the constructive arm only, 2 returns both.

    Returns
    -------
    out_ports=1 : (pout, E_out)
    out_ports=2 : (pout1, E_out1, pout2, E_out2)
    """
    if out_ports == 1:
        E_out = (port_E + tap_E) / np.sqrt(2.0)
        pout = np.sum(np.abs(E_out) ** 2)
        return pout, E_out
    elif out_ports == 2:
        E_out1 = (port_E + tap_E) / np.sqrt(2.0)
        E_out2 = (port_E - tap_E) / np.sqrt(2.0)
        pout1 = np.sum(np.abs(E_out1) ** 2)
        pout2 = np.sum(np.abs(E_out2) ** 2)
        return pout1, E_out1, pout2, E_out2
    else:
        raise Exception("Incorrect Number of Ports. Argument accepted is 1 or 2.")

def halfwave(E, theta=0, rotation=True):
    """Half-wave plate, standard Jones form (no absolute-phase prefactor).

        J_HWP(t) = [[cos 2t,  sin 2t],
                    [sin 2t, -cos 2t]]

    which is what the matrix below is, via cos^2 - sin^2 = cos 2t and
    2 sin t cos t = sin 2t.  Real, so a real input gives a real output --
    a HWP is a reflection on the Poincare sphere and introduces no phase.

    See the RETARDER PHASE CONVENTION note in the module header for why
    the `exp(-i*pi/2)` prefactor this used to carry was removed.
    """
    if not rotation:
        half_matrix = np.array([
            [1, 0],
            [0, -1]
        ], dtype=complex)
        E_half = half_matrix @ np.transpose(E)
        return np.transpose(E_half)
    elif rotation:
        theta = np.radians(theta)
        half_matrix = np.array([
            [np.cos(theta)**2 - np.sin(theta)**2, 2*np.cos(theta)*np.sin(theta)],
            [2*np.cos(theta)*np.sin(theta), -np.cos(theta)**2 + np.sin(theta)**2]
        ], dtype=complex)

        E_half = half_matrix @ np.transpose(E)
        return np.transpose(E_half)
    else:
        raise Exception("Incorrect Rotation of waveplate")


def quarterwave(E, theta=0, rotation=True):
    """Quarter-wave plate, standard Jones form (no absolute-phase prefactor).

    The relative `i` between the fast and slow axes is the **retardance**
    and is the entire physical content of the component -- it stays.  Only
    the overall `exp(-i*pi/4)` prefactor was removed; see the RETARDER
    PHASE CONVENTION note in the module header.

    Consequently an H-polarised input now emerges real rather than rotated
    45 degrees into the complex plane, while H at 45 degrees still becomes
    circular with E_y = -i*E_x.
    """
    if not rotation:
        quarter_matrix = np.array([
            [1, 0],
            [0, -1j]
        ], dtype=complex)
        E_quarter = quarter_matrix @ np.transpose(E)
        return np.transpose(E_quarter)
    elif rotation:
        theta = np.radians(theta)
        quarter_matrix = np.array([
            [np.cos(theta) ** 2 + 1j * np.sin(theta) ** 2, (1 - 1j) * np.sin(theta) * np.cos(theta)],
            [(1 - 1j) * np.sin(theta) * np.cos(theta), np.sin(theta) ** 2 + 1j * np.cos(theta) ** 2]
        ], dtype=complex)
        E_quarter = quarter_matrix @ np.transpose(E)
        return np.transpose(E_quarter)
    else:
        raise Exception("Incorrect Rotation of waveplate")


def polarization_rotator(E, rot_angle):
    theta = np.radians(rot_angle)
    rot_matrix = np.array([
        [np.cos(theta), np.sin(theta)],
        [-np.sin(theta), np.cos(theta)]
    ])

    E_rot = rot_matrix @ np.transpose(E)
    return np.transpose(E_rot)

def polarization_controller(E, qwp1=0, hwp=0,qwp2=0):
    E = quarterwave(E, qwp1)
    E = halfwave(E, hwp)
    E = quarterwave(E, qwp2)
    return E


def voa(E, attenuation_dB):
    """Variable optical attenuator.

    Parameters
    ----------
    E : ndarray (N, 2) — complex-envelope field.
    attenuation_dB : float — attenuation in dB (positive = loss).

    Returns
    -------
    ndarray (N, 2) — attenuated field.
    """
    atten_lin = 10.0 ** (-attenuation_dB / 10.0)
    return E * np.sqrt(atten_lin)


def polarizer(E, polarization):
    if polarization == 'H':
        J_H = np.array([
            [0.99, 0],
            [0, 0.01]
        ])
        return np.transpose(J_H @ np.transpose(E))
    elif polarization == 'V':
        J_V = np.array([
            [0.001, 0],
            [0, 0.99]
        ])
        return np.transpose(J_V @ np.transpose(E))
    elif polarization == '45':
        J_45 = np.array([
            [1 / 2, 1 / 2],
            [1 / 2, 1 / 2]
        ])
        # J_45 = 0.5*np.array([
        #     [1, 1],
        #     [1, 1]
        # ])
        return np.transpose(J_45 @ np.transpose(E))
    elif polarization == '-45':
        J_i45 = np.array([
            [1 / 2, -1 / 2],
            [-1 / 2, 1 / 2]
        ])
        # J_i45 = np.array([
        #     [1, -1],
        #     [-1, 1]
        # ])
        return np.transpose(J_i45 @ np.transpose(E))
    else:
        raise Exception("Wrong Value")

#Version 1
# def beam_splitter(E):
#     # E = E * np.sqrt(power)
#     J_x =(1/np.sqrt(2)) * np.array([[0.95, 0],
#                                     [0, 0.05]])
#
#     J_y = (1/np.sqrt(2)) * np.array([[0.05, 0],
#                                     [0, 0.95]])
#
#     # Calculate the transmitted and reflected components
#     E_x = np.transpose(J_x @ np.transpose(E))
#     E_y = np.transpose(J_y @ np.transpose(E))
#     E_y = halfwave(E_y, rotation=False)
#
#     # Calculate the power in each component
#     P_x = np.linalg.norm(E_x) ** 2
#     P_y = np.linalg.norm(E_y) ** 2
#     print(P_x, P_y)
#     # return P_x, E_x, P_y, E_y
#     return E_x, E_y


def hadamard(E):
    # E is assumed to be an (N,2) array.
    J = (1 / np.sqrt(2)) * np.array([[1, 1],
                                     [1, -1]])
    E_out = E @ J.T  # Apply transformation along the last axis
    return E_out

def circular_analyser(E):
    """Circular-basis analyser: quarter-wave plate followed by a
    polarising beam splitter (QWP+PBS), projecting onto the L/R circular
    basis rather than H/V.

    This is *not* a polarising beam splitter (PHYS-6 in
    opto-sim-issues-and-fixes.md — this function used to be misnamed
    `pbs`). A true PBS projects onto H and V (`diag(1,0)` / `diag(0,1)`,
    see `pbs()` below) and is blind to the relative phase between Ex and
    Ey. This function instead converts that relative phase into an
    intensity imbalance between its two output ports:

        E_x = (Ex - i*Ey) / sqrt(2),   E_y = (-i*Ex + Ey) / sqrt(2)

    which is what phase-encoded schemes (e.g. Section 5's 45°-polarised,
    phase-modulated BB84) need for detection — a true H/V PBS would give
    a fixed 50/50 split regardless of phase and discriminate nothing.
    Any caller relying on this behaviour must call `circular_analyser()`
    explicitly, not `pbs()`.

    Parameters
    ----------
    E : ndarray (N, 2) — complex envelope [Ex, Ey].

    Returns
    -------
    E_x, E_y : ndarray (N,) — the two output-port fields.
    """
    J = (1 / np.sqrt(2)) * np.array([[1, -1j],
                                    [-1j, 1]])
    E_out = E @ J.T  # Apply transformation along the last axis
    E_x = E_out[:, 0]
    E_y = E_out[:, 1]
    return E_x, E_y


def apply_extinction(P_a, P_b, epsilon):
    """Finite analyser extinction, as symmetric power cross-talk.

    A real polarisation analyser does not send every photon to the right
    port.  Duplinskiy et al. (Opt. Express 25(23), 28886, 2017) state
    calibration goal 3 as "Bob's measurements differentiate BB84
    orthogonal states with extinction higher than 98 %", and this is that
    imperfection.

    Parameters
    ----------
    P_a, P_b : float or ndarray — the two analyser port powers.
    epsilon : float — fraction of each port's power appearing at the
        other.  ``0`` is a perfect analyser.

    Returns
    -------
    (P_a', P_b') — power-conserving: ``P_a' + P_b' == P_a + P_b``.

    For a perfectly discriminated state (``P_b = 0``) this leaks exactly
    ``epsilon`` of the power to the wrong port, so the parameter maps
    directly onto the stated quantity with no intermediate calibration.

    Reading the paper's 98 %
    ------------------------
    Two defensible readings differ by exactly a factor of two, and the
    ambiguity falls on the quantity a test of the afterpulse budget turns
    on, so both are carried rather than one chosen (register entry A7):

    - power-fraction: 98 % of the power reaches the correct port, so
      ``epsilon = 1 - 0.98 = 0.0200``;
    - visibility-like: ``(I_max - I_min)/(I_max + I_min) = 0.98``, so
      ``epsilon = I_min/I_max = (1-V)/(1+V) = 0.0101``.

    This is a LUMPED term, and must be described as one
    ------------------------------------------------------
    The paper states extinction end-to-end for Bob's whole measurement
    (PC2 -> PM2 -> PC3 -> PBS), not as a component spec, so one parameter
    stands for at least three mechanisms: PC3 tuning residual, PM2 voltage
    error and PBS leakage.  That is faithful to what is published -- the
    same choice `P_E` represents for Gobby -- but it is not "the PBS
    extinction" and calling it that would be wrong.

    If it is ever decomposed, the components are *not* equivalent: PBS
    leakage is basis-symmetric while a polarisation-controller tuning
    residual is an angular misalignment and therefore basis-dependent, so
    the two are distinguishable by comparing QBER across the X and C
    bases.  That is the route; it is not taken here, because the source
    gives one number.

    Where to apply it
    -----------------
    At the analyser output, **before** detection.  Dark counts and
    afterpulses are generated inside the SPAD and are already
    uncorrelated with Alice's bit; mixing them into each other would be
    meaningless and would inflate the term.
    """
    if epsilon == 0.0:
        return P_a, P_b
    return ((1.0 - epsilon) * P_a + epsilon * P_b,
            (1.0 - epsilon) * P_b + epsilon * P_a)


def pbs(E):
    """Polarising beam splitter: projects onto H and V.

    A true PBS transmits the H component to one port and reflects the V
    component to the other, and is blind to the relative phase between Ex
    and Ey — unlike `circular_analyser()` (formerly misnamed `pbs`, see
    PHYS-6 in opto-sim-issues-and-fixes.md), which projects onto the
    circular basis and *does* depend on that phase.

    Parameters
    ----------
    E : ndarray (N, 2) — complex envelope [Ex, Ey].

    Returns
    -------
    E_x, E_y : ndarray (N,) — the H-port and V-port fields.
    """
    E_x = E[:, 0]
    E_y = E[:, 1]
    return E_x, E_y


def pbc(E_h, E_v):
    """Polarising beam combiner — the inverse of `pbs()`.

    Places two independent fields onto orthogonal polarisation axes of a
    single spatial mode: `E_h` onto H (column 0), `E_v` onto V (column 1).
    In Jones terms the two inputs are acted on by the complementary
    projectors diag(1,0) and diag(0,1) and summed [2], so the operation is
    lossless and `pbs(pbc(a, b)) == (a, b)` exactly.

    This is what makes *deterministic* polarisation routing expressible.
    A non-polarising combiner (`beam_combiner`) adds two fields into the
    same mode, which interferes them and costs half the light at the
    complementary port; a PBC keeps both in full because they occupy
    orthogonal states.

    Parameters
    ----------
    E_h, E_v : ndarray (N,) or (N, 2)
        Fields for the H and V axes.  (N,) is the natural pairing with
        `pbs()`'s return.  An (N, 2) input must be single-polarisation —
        its populated column is taken — since a PBC cannot combine two
        already-orthogonal states without loss.

    Returns
    -------
    ndarray (N, 2) — complex envelope carrying both inputs.

    Raises
    ------
    ValueError — if an (N, 2) input has both columns populated, or if the
        two inputs differ in length.
    """
    def _scalar(E, name):
        E = np.asarray(E)
        if E.ndim == 1:
            return E
        if E.ndim == 2 and E.shape[1] == 2:
            occupied = [c for c in (0, 1) if np.any(np.abs(E[:, c]) > 0)]
            if len(occupied) > 1:
                raise ValueError(
                    f"pbc: {name} carries both polarisations; a polarising "
                    f"beam combiner needs single-polarisation inputs, one "
                    f"per axis."
                )
            return E[:, occupied[0]] if occupied else E[:, 0]
        raise ValueError(f"pbc: {name} must be (N,) or (N, 2), got {E.shape}")

    h = _scalar(E_h, "E_h")
    v = _scalar(E_v, "E_v")
    if h.shape != v.shape:
        raise ValueError(
            f"pbc: inputs must be the same length, got {h.shape} and {v.shape}")

    out = np.zeros((h.shape[0], 2), dtype=complex)
    out[:, 0] = h
    out[:, 1] = v
    return out
#
# Version 2
# def beam_splitter(E, power):
#     # Scale electric field by the square root of power
#     # E = E * np.sqrt(power)
#
#     # Jones matrices for transmitted and reflected components
#     J_t = (1 / np.sqrt(2)) * np.array([[1, 0],
#                                        [0, 1]])
#     J_r = (1 / np.sqrt(2)) * np.array([[1, 0],
#                                        [0, -1]])
#
#     # Calculate transmitted (E_x) and reflected (E_y) fields
#     E_x = np.transpose(J_t @ np.transpose(E))  # Transmitted
#     E_y = np.transpose(J_r @ np.transpose(E))  # Reflected
#
#     # Optional: apply a half-wave plate or other transformations to E_y
#     # E_y = halfwave(E_y, rotation=False) if needed
#
#     # Calculate the power in each component
#     P_x = np.linalg.norm(E_x) ** 2
#     P_y = np.linalg.norm(E_y) ** 2
#
#     return P_x, E_x, P_y, E_y

# Version 3
# def beam_splitter(E):
#     E_x_normalized = E[:, 0]
#     E_y_normalized = E[:, 1]
#     E_x_normalized = np.asarray(E_x_normalized, dtype=complex)
#     E_y_normalized = np.asarray(E_y_normalized, dtype=complex)
#     print(f"E_x_normalized[0]: {E_x_normalized[0]} (Type: {type(E_x_normalized[0])})")
#     print(f"E_y_normalized[0]: {E_y_normalized[0]} (Type: {type(E_y_normalized[0])})")
#     theta = np.angle(E_y_normalized[0]) - np.angle(E_x_normalized[0])
#     print(f'Theta = {theta}')
#
#     T = np.cos(theta)**2  # Transmission coefficient
#     R = np.sin(theta)**2  # Reflection coefficient
#     print(f'T = {T}, R = {R}')
#     # T = 0.5
#     # R = 0.5
#     T_matrix = np.array([[T, 0],
#                          [0, R]])
#
#     E_out = np.transpose(T_matrix @ np.transpose(E))
#     E_x = E_out[:, 0]
#     E_y = E_out[:, 1]
#     P_H = np.linalg.norm(E_x) ** 2
#     P_V = np.linalg.norm(E_y) ** 2
#     return P_H, E_out[:, 0], P_V, E_out[:, 1]

def beam_combiner(P_x, E_x, P_y, E_y, normalized=True):
    if normalized and np.linalg.norm(E_x) ** 2 < 1 and np.linalg.norm(E_y) ** 2 < 1:
        E = E_x + E_y
        pout = P_x+P_y
        return E, pout
    elif not normalized:
        E = E_x + E_y
        P_out = np.linalg.norm(E) ** 2
        return E, P_out
