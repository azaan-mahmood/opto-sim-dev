"""Laser driver: houses the SS-TDM DFB device and drives its injection current.

The driver holds a :class:`~src.lasers.dfblaser.DFBLaser` instance and a
:class:`DriveParams` waveform specification.  The only difference between
CW and gain-switched operation is the injection-current waveform applied
to the same device.

Modes
-----
- ``cw``: constant bias current ``i_bias``; the device settles to its
  steady-state output.
- ``gain_switched``: a train of current pulses (``i_bias`` baseline +
  ``i_peak`` pulses at ``period``) short compared with the carrier
  lifetime, producing relaxation-oscillation optical pulses (Agrawal,
  "Fiber-Optic Communication Systems", 5th ed., 2021, Sec. 3.4; Gobby et
  al., APL 84, 3762 (2004): 100 ps gain-switched DFB pulses at 2.5 MHz).

That is the whole distinction: the device is the same, and only the
injection-current waveform changes.

What is NOT set here
--------------------
Relative intensity noise, chirp and phase statistics are not parameters of
this class, or of the device.  They come out of the spontaneous-emission
seeding and the carrier rate equation in ``dfblaser`` and are *measured*

The drive mode decides which of them is worth looking at:

- under ``cw`` the device settles and shows relative intensity noise with
  a relaxation-oscillation resonance, plus a slowly wandering optical
  phase
- under ``gain_switched`` each pulse grows from spontaneous emission, so
  it carries a fresh random phase and a large chirp across its envelope.

``analysis/validation/validate_dfb_drive.py`` measures both.

Waveforms (``gain_switched`` mode)
----------------------------------
- ``gaussian``: current pulses ``i_peak * exp(-0.5*(tau/sigma)**2)`` with
  ``sigma`` from the FWHM ``width``, centred on each multiple of
  ``period``.  ``tau`` is the distance to the *nearest* period boundary,
  so the leading half of each pulse is carried by the tail end of the
  preceding period and the pulse is symmetric.  Measuring ``tau`` as
  ``t % period`` instead would put the peak at the start of the period
  with nothing before it, giving only the falling half and therefore
  half the requested FWHM.
- ``trapezoidal``: pulses of length ``width`` with linear edges of
  duration ``t_rise`` (square when ``t_rise = 0``), starting at each
  period boundary.  This branch is one-sided by construction, so it
  takes ``tau`` unwrapped.

The drive current is clamped at zero, not rectified.

What gain switching does to the pulse, measured
-----------------------------------------------
Driving 60 -> 200 mA with a 100 ps current pulse, at the default device:

- the **optical pulse comes out narrower than the drive** -- 42 ps FWHM
  from the gaussian, 51 ps from the trapezoidal, against 100 ps of
  current.  That compression is the point of gain switching;
- it also comes out **late**: 76 ps after the gaussian's peak and 143 ps
  after the trapezoidal's.  For the trapezoidal drive that is longer than
  the current pulse itself, so the light forms after the drive has already
  returned to bias;
- each pulse grows from spontaneous emission, so its optical phase is
  uncorrelated with the last one's.

Because of the turn-on delay the device **stops producing separated
pulses below a period of about 210 ps** (4.8 GHz).  Measured extinction
between pulses: 29 dB at 250 ps, 21.7 dB at 220 ps, 17.7 dB at 200 ps,
9.6 dB at 150 ps, and by 50 ps the output is a rippled CW level at the
average current.  The pulse occupies roughly 200 ps of whatever period it
is given, and it cannot be given less.

No warning is raised for this.  The limit moves with drive amplitude and
width, so a fixed threshold in code would be a number nobody sourced.  It
is also the expected ceiling for the device being modelled: Kim's Table I
gives a 0.2 um active layer, a bulk InGaAsP region rather than a quantum
well, and a few GHz is what such a device does.  A faster source is a
different parameter set, which this model accepts.

``analysis/validation/validate_dfb_drive.py`` measures all of the above,
including ``--period-sweep`` which reproduces the extinction table.
"""

from dataclasses import dataclass
from typing import Callable

import numpy as np

from .cwlaser import jones_vector
from .dfblaser import DFBLaser, SimResult

_MODE_CHOICES = ("cw", "gain_switched")
_WAVEFORM_CHOICES = ("gaussian", "trapezoidal")

# Measured: the device settles about 30 ns after turn-on from zero field.
# Runs shorter than that are measuring the transient, not steady operation.
SETTLE_TIME = 40e-9


@dataclass
class DriveParams:
    """Injection-current waveform specification for :class:`LaserDriver`.

    Parameters
    ----------
    mode : str
        ``"cw"`` (constant ``i_bias``) or ``"gain_switched"`` (pulse
        train).
    waveform : str
        ``"gaussian"`` or ``"trapezoidal"``; used in ``gain_switched``
        mode only.  Both are gain switching -- trapezoidal is the squarer
        current drive, gaussian the smooth one.
    i_bias : float
        Baseline current (A).  The CW operating point in ``cw`` mode.
    i_peak : float
        Pulse amplitude above ``i_bias`` (A) in ``gain_switched`` mode.
    period : float
        Pulse repetition period (s) in ``gain_switched`` mode.
    width : float
        Pulse FWHM (``gaussian``) or full base width (``trapezoidal``), in
        s.  The two differ: a trapezoid of base ``width`` with linear edges
        of ``t_rise`` has FWHM ``width - t_rise``.
    t_rise : float
        Rise/fall time of the trapezoidal edges (s); ``0`` gives a square
        pulse.  Used in ``trapezoidal`` mode only.
    """

    mode: str = "cw"
    waveform: str = "gaussian"
    i_bias: float = 100e-3
    i_peak: float = 200e-3
    period: float = 400e-9
    width: float = 100e-12
    t_rise: float = 20e-12

    def __post_init__(self) -> None:
        if self.mode not in _MODE_CHOICES:
            raise ValueError(f"mode must be one of {_MODE_CHOICES}, got {self.mode!r}")
        if self.waveform not in _WAVEFORM_CHOICES:
            raise ValueError(f"waveform must be one of {_WAVEFORM_CHOICES}, got {self.waveform!r}")
        if self.i_bias < 0 or self.i_peak < 0:
            raise ValueError("i_bias and i_peak must be non-negative")
        if self.mode == "gain_switched":
            if self.period <= 0:
                raise ValueError(f"period must be positive in gain_switched mode, got {self.period}")
            if self.width <= 0:
                raise ValueError(f"width must be positive in gain_switched mode, got {self.width}")
            if self.waveform == "trapezoidal" and not 0 <= self.t_rise <= self.width / 2:
                raise ValueError(f"t_rise must be in [0, width/2], got {self.t_rise}")
        self._pulse_sigma = self.width / (2.0 * np.sqrt(2.0 * np.log(2.0)))

    def i(self, t: float | np.ndarray) -> np.ndarray:
        """Injection current at time ``t`` (s), scalar or array; >= 0."""
        t = np.asarray(t, dtype=float)
        if self.mode == "cw":
            return np.full_like(t, self.i_bias)
        tau = t % self.period
        if self.waveform == "gaussian":
            # Distance to the nearest period boundary, not to the start of
            # the period: the pulse is centred on the boundary and its
            # leading half comes from the tail of the preceding period.
            # Without the wrap the peak sits at tau = 0 with nothing before
            # it, so only the falling half survives and the effective FWHM
            # is half of `width`.
            tau = np.minimum(tau, self.period - tau)
            current = self.i_bias + self.i_peak * np.exp(-0.5 * (tau / self._pulse_sigma) ** 2)
        else:  # trapezoidal
            if self.t_rise <= 0:
                amp = (tau < self.width).astype(float)
            else:
                amp = np.where(tau < self.t_rise, tau / self.t_rise,
                               np.where(tau < self.width - self.t_rise, 1.0,
                                        np.where(tau < self.width,
                                                 (self.width - tau) / self.t_rise, 0.0)))
            current = self.i_bias + self.i_peak * np.clip(amp, 0.0, 1.0)
        return np.clip(current, 0.0, None)

    def __call__(self, t: float | np.ndarray) -> np.ndarray:
        """Alias for :meth:`i` -- usable directly as the device current callable."""
        return self.i(t)


class LaserDriver:
    """Drives a :class:`DFBLaser` device with a :class:`DriveParams` waveform.

    ``run()`` resets the device to its initial state (seeded with
    ``seed`` when provided) and drives the injection current through the
    waveform for ``t_end``, returning the recorded facet powers and
    complex envelopes on the *device* timebase.

    ``sample_field()`` is the interface the rest of the project uses: a
    complex ``[Ex, Ey]`` field on a requested grid, matching
    :meth:`CWLaser.sample_field` so the two sources are interchangeable.

    Polarisation lives here rather than on the device.  The device settles
    a *state* -- one pure linear polarisation, the TE mode.  How that axis
    sits relative to lab x and y is mounting and pigtail alignment, and
    that is what ``polarization_azimuth`` and ``polarization_ellipticity``
    express, with the same names, defaults and convention as ``CWLaser``.
    The default (0, 0) puts all the power on Ex.
    """

    def __init__(self, laser: DFBLaser, drive: DriveParams,
                 seed: int | None = None,
                 polarization_azimuth: float = 0.0,
                 polarization_ellipticity: float = 0.0) -> None:
        self.laser = laser
        self.drive = drive
        self.seed = seed
        self.polarization_azimuth = polarization_azimuth
        self.polarization_ellipticity = polarization_ellipticity
        # Set by the most recent sample_field() call: the lasing mode's
        # offset from the Bragg reference the device integrates against.
        self.last_mode_offset_hz: float | None = None
        self._field_cache: dict = {}

    def run(self, t_end: float | None = None, record_every: int = 1) -> SimResult:
        """Drive the device and return the recorded output (SimResult)."""
        if self.seed is not None:
            self.laser.seed = self.seed
        return self.laser.simulate(current=self.drive.i, t_end=t_end, record_every=record_every)

    def sample_field(self, dt: float, n_samples: int,
                     settle: float = SETTLE_TIME,
                     facet: str = 'right') -> np.ndarray:
        """Complex ``[Ex, Ey]`` field on the requested grid, shape (n, 2).

        Mirrors :meth:`CWLaser.sample_field`, so a protocol can take
        either source without knowing which it has.  ``sum(|E|**2,
        axis=1)`` is the optical power in Watts.

        Four things happen between the device and the returned array.

        **The turn-on is discarded.**  The device starts from zero field
        and builds up out of spontaneous emission; measured, it takes
        about 30 ns to settle, so ``settle`` (40 ns by default) is run and
        thrown away before the requested window begins.  That is a fixed
        cost of a few seconds per distinct request, which is why results
        are cached.

        **The lasing mode's carrier offset is removed.**  The device
        integrates against the Bragg reference, and the mode sits well
        off it -- measured 549 GHz at 100 mA, 586 at 120, 591 at 150,
        moving with current because the carrier density shifts the index.
        The device grid carries that fine (Nyquist ~1013 GHz at N=15) but
        a 2 ps chain grid has a Nyquist of 250 GHz and would alias it into
        nonsense.  The offset is estimated from the single-lag
        autocorrelation, divided out, and left on
        ``last_mode_offset_hz``.  This is also the physically right frame:
        a complex envelope is defined against the optical carrier, and the
        carrier here is the lasing mode, not the Bragg reference.

        **The result is decimated by averaging, not slicing.**  The device
        step is always finer, so slicing would alias.  Averaging
        band-limits, which is correct but has a consequence: a coarse grid
        cannot carry fast structure.  The chirp across a gain-switched
        pulse is hundreds of GHz, so asking for it on a 2 ps grid will
        lose most of it -- ask for a finer ``dt`` if the chirp matters.

        **The TE axis is placed on lab axes** with the Jones vector from
        ``polarization_azimuth`` / ``polarization_ellipticity``.  This is a
        rotation of one pure state, not a split into two modes: the device
        emits a single amplitude and the two components keep a fixed phase
        relation.

        Parameters
        ----------
        dt : float
            Output sample interval (s).  Must be at least the device step
            ``laser.dt``.
        n_samples : int
            Number of output samples.
        settle : float
            Turn-on time to run and discard (s).
        facet : str
            ``'right'`` or ``'left'``.

        Returns
        -------
        numpy.ndarray, shape (n_samples, 2), complex
        """
        if facet not in ('right', 'left'):
            raise ValueError(f"facet must be 'right' or 'left', got {facet!r}")
        if n_samples < 1:
            raise ValueError(f"n_samples must be positive, got {n_samples}")
        dt_dev = self.laser.dt
        if dt < dt_dev:
            raise ValueError(
                f"requested dt = {dt:.4e} s is finer than the device step "
                f"{dt_dev:.4e} s; sample_field only decimates. Raise dt, or "
                f"raise n_sections to shorten the device step.")

        key = (dt, n_samples, settle, facet, self.polarization_azimuth,
               self.polarization_ellipticity)
        if key in self._field_cache:
            amp, offset = self._field_cache[key]
            self.last_mode_offset_hz = offset
        else:
            amp, offset = self._scalar_envelope(dt, n_samples, settle, facet)
            self._field_cache[key] = (amp, offset)
            self.last_mode_offset_hz = offset

        e_pol = jones_vector(self.polarization_azimuth,
                             self.polarization_ellipticity)
        return np.outer(amp, e_pol)

    def _scalar_envelope(self, dt: float, n_samples: int, settle: float,
                         facet: str) -> tuple[np.ndarray, float]:
        """Device field on the output grid, carrier removed. Returns (amp, offset_hz)."""
        dt_dev = self.laser.dt
        span = dt * n_samples
        res = self.run(t_end=settle + span, record_every=1)

        keep = res.t >= settle
        e = res.E_right[keep, 0] if facet == 'right' else res.E_left[keep, 0]
        t = res.t[keep] - settle

        # Single-lag autocorrelation frequency estimate.  Preferred over
        # unwrapping the phase because the product conj(E[k])*E[k+1] is
        # weighted by |E[k]|*|E[k+1]|, so the near-zero stretches between
        # gain-switched pulses -- where the phase is pure noise -- carry
        # almost no weight.  Valid while |omega*dt_dev| < pi, i.e. below
        # the device Nyquist, which the measured 1.7-1.83 rad/step is.
        z = np.sum(np.conj(e[:-1]) * e[1:])
        offset_hz = 0.0
        if np.abs(z) > 0.0:
            offset_hz = float(np.angle(z) / dt_dev / (2.0 * np.pi))
            e = e * np.exp(-2j * np.pi * offset_hz * t)

        # Bin-average onto the output grid.  Bins are half-open
        # [k*dt, (k+1)*dt); dt is not generally an integer multiple of the
        # device step, so counts per bin can differ by one and the mean has
        # to be taken per bin rather than by reshaping.
        idx = np.floor(t / dt).astype(int)
        np.clip(idx, 0, n_samples - 1, out=idx)
        counts = np.bincount(idx, minlength=n_samples)
        acc = (np.bincount(idx, weights=e.real, minlength=n_samples) +
               1j * np.bincount(idx, weights=e.imag, minlength=n_samples))
        amp = np.zeros(n_samples, dtype=complex)
        nz = counts > 0
        amp[nz] = acc[nz] / counts[nz]
        return amp, offset_hz
