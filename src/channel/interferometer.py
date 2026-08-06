import numpy as np

# --- LITERATURE SOURCES ---
# [1] Gobby, C., Yuan, Z. L. & Shields, A. J., "Quantum key distribution
#     over 122 km of standard telecom fiber", Appl. Phys. Lett. 84,
#     3762-3764, 2004.  — time-division Mach-Zehnder interferometer
#     with 5.8 ns delay for phase-encoded BB84.
# [2] Bennett, C. H. & Brassard, G., "Quantum cryptography: Public key
#     distribution and coin tossing", Proc. IEEE Int. Conf. Comput. Syst.
#     Signal Process., pp. 175-179, Bangalore, 1984.
#     — BB84 protocol implemented via the AMZI.
# [3] Townsend, P. D., "Quantum key distribution over distances up to
#     30 km", Electron. Lett. 29(14), 1291-1293, 1993.
#     — first demonstration of fibre-based phase-encoded QKD.
# [4] Gobby, C., Yuan, Z. L. & Shields, A. J., "Quantum key distribution
#     over 122 km of standard telecom fiber", Appl. Phys. Lett. 84(19),
#     3762-3764, 2004. — Gobby's 3.3% QBER floor at short range is
#     interferometer visibility: e_opt = (1 - V)/2, so V = 0.934.


class AsymmetricMZI:
    """Unbalanced Mach-Zehnder interferometer for time-bin QKD.

    Two operating modes:

    **Encoder** (Alice): splits an input pulse into two time-separated
    copies (early / late) via a 50:50 splitter, delays one arm by
    ``delay`` seconds, optionally applies a phase shift to the delayed
    arm, and recombines.  The output contains two non-overlapping time
    bins separated by ``delay``.

    **Decoder** (Bob): receives a two-time-bin field, splits it, delays
    one copy so the early pulse aligns with the late pulse, applies a
    basis-dependent phase shift, and outputs to constructive /
    destructive ports via a second 50:50 combiner.  The power at each
    port follows the standard interference pattern:

        P_c  ∝  1 + V·cos(phi_alice - phi_bob + phi_err)
        P_d  ∝  1 - V·cos(phi_alice - phi_bob + phi_err)

    where ``visibility`` V ∈ (0, 1] is the interferometric contrast
    (implemented as a combiner amplitude imbalance, see [4]) and
    ``phase_error`` is a static decoder phase offset representing
    imperfect AMZI path-length matching.  Both apply to the decoder
    only; the encoder is an ideal 50:50 splitter.

    Parameters
    ----------
    delay : float
        Differential delay between arms (seconds).  Must be positive.
    mode : {'encoder', 'decoder'}
        Operating mode (default 'encoder').
    insertion_loss_db : float or None
        Total device insertion loss in dB.  None (default) means ideal
        (0 dB, no insertion loss).
    visibility : float
        Interferometric visibility in (0, 1], default 1.0 (ideal).
        Decoder only.  Optical misalignment error e_opt = (1 - V)/2
        (Gobby et al. 2004 [4]: 3.3% floor <-> V ≈ 0.934).
    phase_error : float
        Static phase offset (radians) applied to the delayed arm in
        decoder mode, default 0.0.  Represents imperfect path-length
        matching (delta = 2*pi*dL/lambda).

    References
    ----------
    Gobby et al. (2004) [1] — time-division MZI for 122 km QKD.
    Bennett & Brassard (1984) [2] — BB84 protocol.
    Townsend (1993) [3] — first fibre phase-encoded QKD.
    Gobby et al. (2004) [4] — visibility floor e_opt = (1 - V)/2.
    """

    def __init__(self, delay, mode='encoder', insertion_loss_db=None,
                 visibility=1.0, phase_error=0.0):
        if delay <= 0:
            raise ValueError(f"delay must be positive, got {delay}")
        if not 0.0 < visibility <= 1.0:
            raise ValueError(f"visibility must be in (0, 1], got {visibility}")
        self.delay = delay
        self.mode = mode
        self.visibility = float(visibility)
        self.phase_error = float(phase_error)
        self._il_lin = 1.0 if insertion_loss_db is None else \
            10.0 ** (-insertion_loss_db / 10.0)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def modulate(self, E, dt, phase=None):
        """Apply the AMZI transformation.

        Parameters
        ----------
        E : ndarray, shape (N, 2)
            Complex-envelope field [Ex, Ey].
        dt : float
            Sampling interval in seconds.  Used to convert ``delay``
            to an integer number of samples.
        phase : float or None
            Phase shift (radians) applied to the **delayed** arm.
            None means no phase shift.

        Returns
        -------
        encoder mode : ndarray, shape (N, 2)
            Field with two time bins (early + delayed).
        decoder mode : tuple (ndarray, ndarray), each shape (N, 2)
            (E_constructive, E_destructive) — two output ports.
        """
        delay_samples = int(self.delay / dt)
        if delay_samples == 0:
            delay_samples = 1  # minimum one-sample delay for AMZI operation

        # 50:50 split (Hadamard-like)
        E_short = E / np.sqrt(2)
        E_long = E / np.sqrt(2)

        # Delay the long arm (shift right, zero-fill)
        E_long = np.roll(E_long, delay_samples, axis=0)
        E_long[:delay_samples] = 0.0

        # Phase shift on the delayed arm
        if phase is not None:
            E_long *= np.exp(1j * phase)
        if self.mode == 'decoder' and self.phase_error != 0.0:
            E_long *= np.exp(1j * self.phase_error)

        if self.mode == 'encoder':
            E_out = E_short + E_long
            if self._il_lin != 1.0:
                E_out *= np.sqrt(self._il_lin)
            return E_out
        else:  # decoder
            # Finite-visibility combiner: unequal arm amplitudes r, s with
            # r^2 + s^2 = 1 and 2rs = V.  At interference this yields
            # P_c ∝ (1 + V·cos Δφ), P_d ∝ (1 - V·cos Δφ), so the optical
            # misalignment error is e_opt = (1 - V)/2.  For V = 1 this
            # reduces exactly to the ideal 50:50 combiner.
            v = self.visibility
            r = np.sqrt(0.5 + 0.5 * np.sqrt(max(0.0, 1.0 - v * v)))
            s = np.sqrt(0.5 - 0.5 * np.sqrt(max(0.0, 1.0 - v * v)))
            E_c = r * E_short + s * E_long
            E_d = r * E_short - s * E_long
            if self._il_lin != 1.0:
                E_c *= np.sqrt(self._il_lin)
                E_d *= np.sqrt(self._il_lin)
            return E_c, E_d

    # ------------------------------------------------------------------
    # Utility
    # ------------------------------------------------------------------

    def __repr__(self):
        return (f"AsymmetricMZI(delay={self.delay:.2e}s, "
                f"mode={self.mode!r}, visibility={self.visibility:.3f}, "
                f"phase_error={self.phase_error:.3f})")
