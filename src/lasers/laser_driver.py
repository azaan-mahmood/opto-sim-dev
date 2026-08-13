"""Laser driver: houses the SS-TDM DFB device and drives its injection current.

The driver holds a :class:`~src.lasers.dfblaser.Laser` instance and a
:class:`DriveParams` waveform specification.  The only difference between
CW and gain-switched operation is the injection-current waveform applied
to the same device.

Modes
-----
- ``cw``: constant bias current ``i_bias``; the device settles to its
  steady-state output.
- ``pulsed``: gain-switched operation -- a train of current pulses
  (``i_bias`` baseline + ``i_peak`` pulses at ``period``) short compared
  with the carrier lifetime, producing relaxation-oscillation optical
  pulses (Agrawal, "Fiber-Optic Communication Systems", 5th ed., 2021,
  Sec. 3.4; Gobby et al., APL 84, 3762 (2004): 100 ps gain-switched DFB
  pulses at 2.5 MHz).

Waveforms (``pulsed`` mode)
---------------------------
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

The drive current is clamped to be non-negative everywhere.
"""

from dataclasses import dataclass
from typing import Callable

import numpy as np

from .dfblaser import Laser, SimResult

_MODE_CHOICES = ("cw", "pulsed")
_WAVEFORM_CHOICES = ("gaussian", "trapezoidal")


@dataclass
class DriveParams:
    """Injection-current waveform specification for :class:`LaserDriver`.

    Parameters
    ----------
    mode : str
        ``"cw"`` (constant ``i_bias``) or ``"pulsed"`` (gain-switched
        pulse train).
    waveform : str
        ``"gaussian"`` or ``"trapezoidal"``; used in ``pulsed`` mode only.
    i_bias : float
        Baseline current (A).  The CW operating point in ``cw`` mode.
    i_peak : float
        Pulse amplitude above ``i_bias`` (A) in ``pulsed`` mode.
    period : float
        Pulse repetition period (s) in ``pulsed`` mode.
    width : float
        Pulse FWHM (``gaussian``) or pulse width (``trapezoidal``), in s.
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
        if self.mode == "pulsed":
            if self.period <= 0:
                raise ValueError(f"period must be positive in pulsed mode, got {self.period}")
            if self.width <= 0:
                raise ValueError(f"width must be positive in pulsed mode, got {self.width}")
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
    """Drives a :class:`Laser` device with a :class:`DriveParams` waveform.

    ``run()`` resets the device to its initial state (seeded with
    ``seed`` when provided) and drives the injection current through the
    waveform for ``t_end``, returning the recorded facet powers and
    complex envelopes.
    """

    def __init__(self, laser: Laser, drive: DriveParams, seed: int | None = None) -> None:
        self.laser = laser
        self.drive = drive
        self.seed = seed

    def run(self, t_end: float | None = None, record_every: int = 1) -> SimResult:
        """Drive the device and return the recorded output (SimResult)."""
        if self.seed is not None:
            self.laser.seed = self.seed
        return self.laser.simulate(current=self.drive.i, t_end=t_end, record_every=record_every)
