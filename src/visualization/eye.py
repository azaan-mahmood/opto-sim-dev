"""NRZ on-off-keying eye diagram, generated through a real modulator.

The eye is the standard way to look at what a source does to a bit
stream: every two-unit-interval slice of the intensity is overlaid, so
jitter shows as horizontal blur at the crossings, amplitude noise as
vertical thickness of the rails, and finite bandwidth as a closing of the
opening in the middle.

Why the modulator is a real one
-------------------------------
The signal is carved by an `MZM` driven between 0 V and its own
`switching_voltage`, not by multiplying the field with a 0/1 waveform.
Those are different pictures. A multiply gives perfectly flat rails and
instantaneous edges by construction, so the eye can only show what the
*source* contributes. The device's transfer is
`cos^2(pi*V / (2*V_pi))`, so the rail levels, the extinction between
them, and the shape of the transition all come from the modulator, and
the eye shows the pair.

Source-agnostic by contract
---------------------------
The only thing required of `laser` is `sample_field(dt, n)` returning a
complex `(n, 2)` field in Watts under `sum(|E|**2, axis=1)`. Both
`CWLaser` and `LaserDriver` provide it -- `LaserDriver.sample_field` says
so in as many words, "so a protocol can take either source without
knowing which it has" -- which is why this serves the CW source and the
gain-switched DFB unchanged.

References
----------
[1] Agrawal, G. P., "Fiber-Optic Communication Systems", 5th ed., Wiley,
    2021, §4.3.  Eye-diagram interpretation for digital lightwave
    systems.
[2] ITU-T G.957 / IEEE 802.3 optical eye-mask test conditions: reference
    receiver is a fourth-order Bessel-Thomson response at 0.75 times the
    line rate.
"""
import numpy as np

from src.channel.mzm import MZM


def _align(E, axis):
    """Rotate the field's polarisation onto `axis` (0 = Ex, 1 = Ey).

    The polarisation controller a real bench puts in front of its
    modulator, and the same thing Duplinskiy's PC1 does before PM1.

    It is needed rather than decorative.  The MZM is X-cut, so it phase-
    modulates **Ey** and leaves Ex untouched, while the sources here do not
    agree on where they put their light: `LaserDriver` emits entirely on
    Ex, so an unaligned DFB eye is a flat line with the modulator doing
    nothing at all, and `CWLaser` at 45 degrees splits the power evenly, so
    the unmodulated half sits underneath as a CW floor and caps the
    extinction near 3 dB.  Neither is a fact about the modulator, which is
    what the figure would appear to be showing.

    Implemented as a unitary, so it is lossless and reversible rather than
    a projection that quietly throws half the light away: with the source's
    normalised Jones vector `v`, the matrix `[[v0*, v1*], [-v1, v0]]` has
    unit determinant and maps `v` to the requested axis.  Exact for a
    source whose polarisation is fixed across the window, which both of
    these are.
    """
    v = E.mean(axis=0)
    n = np.linalg.norm(v)
    if n == 0.0:
        return E
    v = v / n
    U = np.array([[np.conj(v[0]), np.conj(v[1])],
                  [-v[1], v[0]]], dtype=complex)
    if axis == 1:
        U = U[::-1]
    return np.transpose(U @ np.transpose(E))


def eye_diagram(laser, bitrate, n_bits=128, spb=64, ax=None, title=None,
                seed=None, align=True, extinction_ratio_db=None,
                drive_bandwidth=None):
    """Overlay every 2-UI slice of an NRZ-OOK intensity waveform.

    Parameters
    ----------
    laser : object exposing `sample_field(dt, n_samples) -> (n, 2)` complex.
        `CWLaser` and `LaserDriver` both qualify.
    bitrate : float — baud, e.g. 10e9.
    n_bits : int — bits in the pattern (default 128).
    spb : int — samples per bit (default 64).  Sets the horizontal
        resolution of the eye; the sampling interval is `1/(bitrate*spb)`.
    ax : matplotlib Axes or None — drawn on a new figure if None.
    title : str or None — overrides the default title.
    seed : int or None — seeds the BIT PATTERN, and nothing else.

        **Pass one.**  The pattern is random, so without a seed the figure
        differs on every run and cannot be compared against a previous
        version of itself.  The draw uses its own `default_rng`, never the
        global `numpy.random` state, so seeding it here cannot shift the
        bit and basis choices a protocol makes from that global stream.

        It does **not** make the figure reproducible on its own, because
        the source's own noise -- phase noise, RIN, spontaneous emission --
        is the source's business and is drawn before this function sees
        anything.  Pinning that is the caller's job, and the two callers
        here already do it: `LaserDriver` takes a `seed` argument, while
        `CWLaser` draws from the global `numpy.random` state, which
        `validate_cwlaser.py` seeds at import.  With the source pinned, a
        fixed `seed` reproduces the eye exactly.
    align : bool — rotate the source's polarisation onto the modulated axis
        first (default True).  This is the polarisation controller a real
        bench has, and without it the figure shows the source's
        polarisation rather than the modulator's extinction -- see
        `_align`.  Pass False to see that unaligned case deliberately.
    extinction_ratio_db : float or None — passed to `MZM`, whose own
        default this matches.  None means the ideal device: infinite
        extinction, so the off rail sits at exactly zero and **anything
        visible above it came from the source**, which is usually what the
        figure is for.  A real LiNbO3 modulator reaches 20-30 dB; set one
        to see a device floor instead of a source-only picture.
    drive_bandwidth : float or None — 3 dB bandwidth of the drive
        electronics, in Hz.  None (default) uses 0.75 * bitrate, the
        reference-receiver bandwidth the eye-mask standards specify [2].
        0 disables the filter, which gives vertical edges and no eye
        opening -- useful only to see what the limit is doing.

    Returns
    -------
    ndarray (n_eyes, 2*spb) — the overlaid intensity traces, in Watts.

    Raises
    ------
    ValueError — if the pattern is too short to form at least two traces.
        Returning quietly would leave an empty pair of axes looking like a
        measurement.
    """
    dt = 1.0 / (bitrate * spb)
    Tbit = 1.0 / bitrate

    mzm = MZM(extinction_ratio_db=extinction_ratio_db)
    V_off = mzm.switching_voltage      # drives the modulator to its null

    rng = np.random.default_rng(seed)
    bits = rng.integers(0, 2, n_bits)
    wfm = np.repeat(bits, spb).astype(float)

    E_cw = laser.sample_field(dt, len(wfm))
    if align:
        # The X-cut MZM modulates Ey; put the source's light there.  See
        # `_align` -- without it the DFB eye is flat and the CW eye carries
        # a floor that reads as the modulator's extinction limit.
        E_cw = _align(E_cw, 1 if mzm._pm.crystal_cut == 'X' else 0)

    V_signal = np.where(wfm > 0.5, 0.0, V_off)

    # Band-limit the DRIVE, which is where the limit physically sits: the
    # RF amplifier and the modulator's electrode have finite bandwidth, the
    # optical field does not care directly.
    #
    # Without this the eye has nothing to show.  An ideal rectangular drive
    # steps between rails inside one sample, so the edges are vertical, the
    # crossings are a single point and there is no opening -- the figure
    # becomes two flat rails and a line, which is what the version this was
    # ported from produced.  Rise time is the whole content of an eye.
    #
    # Fourth-order Bessel-Thomson at 0.75 * baud is the reference-receiver
    # response the eye-mask standards specify [2], so it is a convention
    # with a citation rather than a filter chosen to look right.  Causal
    # (`lfilter`, not `filtfilt`): a zero-phase filter would let the rails
    # start moving before the bit they belong to.
    lag = 0
    if drive_bandwidth is None:
        drive_bandwidth = 0.75 * bitrate
    if drive_bandwidth > 0:
        from scipy.signal import bessel, lfilter, group_delay
        wn = min(drive_bandwidth / (0.5 / dt), 0.99)
        b, a = bessel(4, wn, btype='low', norm='mag')
        V_signal = lfilter(b, a, V_signal)
        # A causal filter delays everything, so without correcting for it
        # the crossings drift off the window edges and the eye opening
        # straddles the boundary.  Undoing it is what a scope's clock
        # trigger does.
        #
        # The delay is read from the filter itself at DC rather than
        # measured off the waveform.  Bessel is the maximally-flat-delay
        # family, which is why the standards pick it, so one number is
        # accurate across the band.
        #
        # Discarded rather than rolled.  `np.roll` would wrap the tail of
        # the record onto the front and draw one trace joining samples that
        # are a whole pattern apart, which shows up as a stray vertical
        # line through the eye.  One extra bit goes with it, because the
        # filter starts from rest and the first edge rises from a zero
        # initial condition rather than from the previous bit.
        lag = int(round(float(group_delay((b, a), w=[0.0])[1][0]))) + spb

    E_mod = mzm.modulate(E_cw, V_signal)
    I = np.sum(np.abs(E_mod) ** 2, axis=1)
    if lag:
        I = I[lag:]

    eye_len = 2 * spb
    n_eyes = len(I) // eye_len
    if n_eyes < 2:
        raise ValueError(
            f"n_bits={n_bits} at spb={spb} gives {n_eyes} trace(s); an eye "
            f"needs at least 2 to overlay. Raise n_bits.")

    eye = I[:n_eyes * eye_len].reshape(n_eyes, eye_len)
    t_eye = np.arange(eye_len) * dt / Tbit

    import matplotlib.pyplot as plt
    if ax is None:
        _, ax = plt.subplots()
    for row in eye:
        ax.plot(t_eye, row * 1e3, 'steelblue', lw=0.4, alpha=0.5)

    ax.set_xlabel('Time (UI)')
    ax.set_ylabel('Intensity (mW)')
    ax.set_title(title or f'Eye diagram @ {bitrate / 1e9:.1f} Gbaud (MZM)')
    ax.set_xlim(0, 2)
    return eye
