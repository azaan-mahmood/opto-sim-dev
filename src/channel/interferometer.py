import numpy as np

from src.channel.optics import coupler_split

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
#     3762-3764, 2004. — source of the 5.8 ns delay and the visibility
#     relation e_opt = (1 - V)/2.
#
#     `visibility` is a decoder imperfection: finite fringe contrast in
#     the device itself.  It is not the place to reproduce a link's
#     measured error rate.  Where a link budget carries background counts,
#     the observed fringe visibility follows from them as
#     V = S/(S + 2*P_e), with S the signal click probability and P_e the
#     background per clock — so injecting a visibility derived from a
#     measured QBER on top of that budget applies the same physics twice.


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

    Unbalanced splitting and arm access
    -----------------------------------
    ``split_ratio`` sets the input coupler: the *short* arm receives that
    fraction of the power and the long arm the rest, via `coupler_split`
    (Zeilinger [5], real convention).  The default 0.5 is an ideal 3 dB
    coupler and reproduces the balanced behaviour exactly.

    The arms can also be taken out and put back:

    * encoder with ``recombine=False`` returns ``(E_short, E_long)``
      instead of their sum, so the two arms can be routed independently —
      for example onto orthogonal polarisations with `optics.pbc`;
    * decoder accepts either a single field (which it splits itself, as
      usual) or an already-separated ``(E_a, E_b)`` pair, e.g. straight
      out of `optics.pbs`.

    Together these express interferometers whose two arms travel by
    distinguishable routes rather than sharing one fibre — polarisation
    multiplexing being the common case.

    Parameters
    ----------
    delay : float
        Differential delay between arms (seconds).  Must be positive.
    mode : {'encoder', 'decoder'}
        Operating mode (default 'encoder').
    split_ratio : float in (0, 1)
        Power fraction into the short arm (default 0.5 = 3 dB).
    phase_arm : {'long', 'short'}
        Which arm the phase modulator sits in (default 'long').  A real
        interferometer has its modulator in one specific arm, and which
        one it is fixes the *sign* of the relative phase seen at the
        output.  For a symmetric device the choice is unobservable —
        intensity goes as cos(delta), which is even — but it stops being
        unobservable the moment the arms differ (unequal loss, a second
        modulator, chirp), so it is a property of the device rather than a
        convention.  `bb84_time_bin.py` wires a chain where the encoded
        pulse travels the short arm.
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
    phase_drift_rad_s : float
        Rate (rad/s) at which the arm imbalance drifts, default 0.0 (no
        drift).  Applied on the delayed arm alongside ``phase_error``, as
        ``phase_error + phase_drift_rad_s * t``, with ``t`` passed to
        `modulate`.  Decoder only, by the same convention as
        ``phase_error``: the two devices in a link both drift, and the
        observable is their *net* relative phase, so a chain carries the
        pair's combined rate on the decoder rather than double-counting it
        at both ends.

        This is **arm-length** drift -- thermal and convective change in
        the relative length of the two arms -- which is a property of the
        interferometer.  It is not modulator bias drift; see
        `PhaseModulator` for that mechanism and why it is not modelled.

        The model is a linear ramp, which is the conservative reading of a
        drift *rate bound*: a source that states "less than r per second"
        constrains the accumulated phase to r*t, and a ramp saturates that
        bound where a random walk would sit below it.  Gobby et al. (2004)
        [1] measure < 0.05 deg/s (8.727e-4 rad/s) with both setups cased in
        enclosures to prevent air convection, over key transfers of about
        two minutes -- 6.0 deg accumulated, contributing 0.091% QBER.

    References
    ----------
    Gobby et al. (2004) [1] — time-division MZI for 122 km QKD.
    Bennett & Brassard (1984) [2] — BB84 protocol.
    Townsend (1993) [3] — first fibre phase-encoded QKD.
    Gobby et al. (2004) [4] — visibility floor e_opt = (1 - V)/2.
    Zeilinger (1981) [5] — lossless 2x2 coupler; see the phase-convention
        note in `optics.py`.  This class uses the real convention, so the
        interference term goes as cos(delta_phi).
    """

    def __init__(self, delay, mode='encoder', split_ratio=0.5,
                 phase_arm='long', insertion_loss_db=None, visibility=1.0,
                 phase_error=0.0, phase_drift_rad_s=0.0):
        if delay <= 0:
            raise ValueError(f"delay must be positive, got {delay}")
        if not 0.0 < visibility <= 1.0:
            raise ValueError(f"visibility must be in (0, 1], got {visibility}")
        if not 0.0 < split_ratio < 1.0:
            raise ValueError(
                f"split_ratio must be in (0, 1), got {split_ratio}")
        if phase_arm not in ('long', 'short'):
            raise ValueError(
                f"phase_arm must be 'long' or 'short', got {phase_arm!r}")
        self.delay = delay
        self.mode = mode
        self.split_ratio = float(split_ratio)
        self.phase_arm = phase_arm
        self.visibility = float(visibility)
        self.phase_error = float(phase_error)
        self.phase_drift_rad_s = float(phase_drift_rad_s)
        self._il_lin = 1.0 if insertion_loss_db is None else \
            10.0 ** (-insertion_loss_db / 10.0)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def arm_phase_offset(self, t=0.0):
        """Delayed-arm phase imbalance at elapsed time ``t`` (seconds).

        The single source of truth for this device's static offset plus
        accumulated drift.  Callers that cannot afford to push fields
        through `modulate` per shot -- the time-bin chain evaluates a closed
        form instead, see `bb84_time_bin` -- must obtain the offset here
        rather than re-expressing ``phase_error + rate * t`` themselves, so
        there is only ever one copy of the law to keep correct.

        Returns 0.0 in encoder mode, matching `modulate`.
        """
        if self.mode != 'decoder':
            return 0.0
        return self.phase_error + self.phase_drift_rad_s * t

    def modulate(self, E, dt, phase=None, recombine=True, t=0.0):
        """Apply the AMZI transformation.

        Parameters
        ----------
        E : ndarray (N, 2), or tuple of two such arrays
            Complex-envelope field [Ex, Ey].  A ``(E_a, E_b)`` tuple is
            accepted in **decoder** mode only and is taken as the two arms
            already separated upstream (e.g. by `optics.pbs`); the input
            coupler is then skipped, `E_a` becoming the short arm and
            `E_b` the long one.
        dt : float
            Sampling interval in seconds.  Used to convert ``delay``
            to an integer number of samples.
        phase : float or None
            Phase shift (radians) applied to the **delayed** arm.
            None means no phase shift.
        recombine : bool
            Encoder mode only.  True (default) returns the summed field.
            False returns ``(E_short, E_long)`` unrecombined, so the arms
            can be routed separately — e.g. onto orthogonal polarisations
            via `optics.pbc`.
        t : float
            Elapsed time in seconds since the start of the run, used only
            to accumulate ``phase_drift_rad_s``.  Default 0.0 leaves the
            device exactly as it was before drift existed.

        Returns
        -------
        encoder, recombine=True : ndarray (N, 2)
            Field with two time bins (early + delayed).
        encoder, recombine=False : (E_short, E_long), each (N, 2)
        decoder : (E_constructive, E_destructive), each (N, 2)
        """
        pre_split = isinstance(E, tuple)
        if pre_split and self.mode != 'decoder':
            raise ValueError(
                "a pre-split (E_a, E_b) input is only meaningful in decoder "
                "mode; the encoder is what performs the split.")

        ref = E[0] if pre_split else E
        delay_samples = int(self.delay / dt)
        if delay_samples == 0:
            delay_samples = 1  # minimum one-sample delay for AMZI operation

        if pre_split:
            # Arms arrived separately (e.g. routed by a PBS); the input
            # coupler has already been applied wherever they were split.
            E_short, E_long = E[0], E[1].copy()
        else:
            # Input coupler: short arm takes `split_ratio` of the power.
            _, E_short, _, E_long = coupler_split(1.0, E,
                                                  ratio=self.split_ratio)

        # Delay the long arm (shift right, zero-fill)
        E_long = np.roll(E_long, delay_samples, axis=0)
        E_long[:delay_samples] = 0.0

        # Phase shift, on whichever arm carries the modulator
        if phase is not None:
            if self.phase_arm == 'long':
                E_long = E_long * np.exp(1j * phase)
            else:
                E_short = E_short * np.exp(1j * phase)
        # Path-length mismatch, static and drifting, is a property of the
        # delayed arm.  `arm_phase_offset` owns the law; see its docstring.
        arm_offset = self.arm_phase_offset(t)
        if arm_offset != 0.0:
            E_long = E_long * np.exp(1j * arm_offset)

        if self.mode == 'encoder':
            if not recombine:
                if self._il_lin != 1.0:
                    root = np.sqrt(self._il_lin)
                    return E_short * root, E_long * root
                return E_short, E_long
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
                f"mode={self.mode!r}, split_ratio={self.split_ratio:.4f}, "
                f"phase_arm={self.phase_arm!r}, "
                f"visibility={self.visibility:.3f}, "
                f"phase_error={self.phase_error:.3f}, "
                f"phase_drift={self.phase_drift_rad_s:.3e} rad/s)")
