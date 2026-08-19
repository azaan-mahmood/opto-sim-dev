"""Piezo-driven fibre stretcher: an optical phase shifter you can buy.

A length of fibre is bonded to a piezoelectric ceramic.  Driving the
ceramic elongates it, stretching the fibre, and the optical path changes
by more than the mechanical extension alone -- the strain also changes the
refractive index through the photoelastic effect.  Butter & Hocker [1]
measured the combined coefficient at about 1.2e7 radians per unit strain
per metre at 633 nm, matching their theory; it scales as 1/lambda, so
1550 nm light sees roughly 4.9e6 m^-1.

This is the device Gobby et al. [2] use to hold Bob's operating point:
*"The phase was varied by applying a DC bias to the piezo-driven fiber
stretcher in the long arm of Bob's interferometer."*  It is why their
drift costs bit rate and not QBER -- a stretcher turns a phase and cannot
return photons, so it corrects the half of a fibre rotation that is a
phase and leaves the half that is a loss.

Defaults: Thorlabs FVP155P
--------------------------
Every default here is one line of DOC-103641 Rev B, so there is a single
place to correct if the part changes:

    half-wave voltage      < 20 V
    phase stroke           7*pi +/- 1*pi at 150 V, 1 kHz
    drive range            0 - 150 V
    resonant frequency     80 kHz +/- 15 %
    insertion loss         < 0.1 dB (without connectors)
    residual AM            < 0.15 %
    operating wavelength   1290 - 1625 nm
    fibre                  PM, equivalent to PM1300-XP

Two of those are bounds rather than values, and are treated as the worst
device the specification permits.  See `v_pi` and `residual_am`.

What is NOT modelled
--------------------
The datasheet's stroke-versus-voltage plot is visibly nonlinear: roughly
straight to about 75 V, then a plateau near 90-125 V, then rising again to
7.2*pi at 150 V.  This models the law as linear through `v_pi`.  Digitising
the published curve would be fitting a figure, which is how a number ends
up wearing a device's name without being a measurement of it.

Residual amplitude modulation is carried as a specification and not
applied, because the datasheet bounds its size without giving its
dependence on drive.  Applying it would mean inventing that shape.  At
0.15 % it is also well under what a Monte Carlo run here can resolve.

References
----------
[1] Butter, C. D. & Hocker, G. B., "Fiber optics strain gauge",
    Appl. Opt. 17(18), 2867-2869, 1978.  Strain-optic phase coefficient.
[2] Gobby, C., Yuan, Z. L., & Shields, A. J., "Quantum key distribution
    over 122 km of standard telecom fiber", Appl. Phys. Lett. 84(19),
    3762-3764, 2004.
[3] Thorlabs, "Piezo-Based PM Fiber Phase Shifter, 150 V, 7pi Phase
    Stroke, FC/PC (FVP155P)", DOC-103641 Rev B, 22 January 2026.
"""
import math

import numpy as np

from .optics import voa

TWO_PI = 2.0 * math.pi


class PiezoFibreStretcher:
    """A fibre stretcher driven in volts, delivering optical phase.

    Interface mirrors `PhaseModulator`, which is the component it most
    resembles: a drive voltage in, a phase out, with `v_pi` setting the
    scale.

    Parameters
    ----------
    v_pi : float — drive voltage for pi radians (default 20.0 V).

        The datasheet says "< 20 V", a BOUND.  Taking 20.0 models the
        worst device the specification permits, which is the conservative
        reading and is a documented assumption rather than a measurement.
        A faster part needs less voltage for the same phase, so nothing
        here becomes unreachable by assuming the slow end.
    v_max : float — maximum drive voltage (default 150.0 V).  The
        datasheet warns that exceeding it shortens the device's life and
        that reverse bias can destroy it, so the range is [0, v_max] and
        going outside raises rather than clipping.
    stroke_rad : float — phase reachable at `v_max` (default 7*pi).
    resonant_hz : float — piezo resonance (default 80e3).  Sets the
        bandwidth: a servo cannot usefully re-lock faster than this, and
        the datasheet derates the usable stroke well below it -- 8*pi at
        2 kHz down to 2*pi at 20 kHz -- to keep the thermal cutout from
        tripping.  Carried so a caller can check its update rate against
        the part rather than against nothing.
    insertion_loss_db : float — fixed optical loss (default 0.1 dB).
    apply_insertion_loss : bool — whether `apply()` imposes it
        (default True).

        True is the honest default: a real device always has the loss, and
        one that silently lost it would flatter every budget it appeared
        in.  Suppress it at the CALL SITE, where the reason lives.  The
        Gobby chain does exactly that, because its `ETA_BOB = 0.045`
        already folds in "5 dB of loss in Bob's apparatus" and the
        stretcher sits inside that apparatus -- applying it again would
        double-count against a number taken from the paper.
    residual_am : float — residual amplitude modulation (default 0.0015).
        Carried, reported, and deliberately not applied; see the module
        docstring.
    """

    def __init__(self, v_pi=20.0, v_max=150.0, stroke_rad=7.0 * math.pi,
                 resonant_hz=80e3, insertion_loss_db=0.1,
                 apply_insertion_loss=True, residual_am=0.0015):
        if v_pi <= 0:
            raise ValueError(f"v_pi must be positive, got {v_pi}")
        if v_max <= 0:
            raise ValueError(f"v_max must be positive, got {v_max}")
        self.v_pi = float(v_pi)
        self.v_max = float(v_max)
        self.stroke_rad = float(stroke_rad)
        self.resonant_hz = float(resonant_hz)
        self.insertion_loss_db = float(insertion_loss_db)
        self.apply_insertion_loss = bool(apply_insertion_loss)
        self.residual_am = float(residual_am)

    # ------------------------------------------------------------------

    def phase_for(self, v):
        """Optical phase delivered by drive voltage `v`, in radians."""
        if not 0.0 <= v <= self.v_max:
            raise ValueError(
                f"drive voltage {v:g} V is outside the device range "
                f"[0, {self.v_max:g}] V. The datasheet warns that "
                f"exceeding it shortens the device's life and that reverse "
                f"bias can destroy it, so this raises rather than clipping.")
        return math.pi * v / self.v_pi

    def voltage_for(self, phase_rad):
        """Drive voltage that delivers `phase_rad`, wrapped into range.

        Wrapping is what makes a finite stroke sufficient.  An
        interferometer's operating point is periodic, so a correction of
        5*pi and one of pi are the same point and the device flies back to
        the nearer -- which is what a real servo does when it runs out of
        travel.  Without wrapping, a correction tracking a steadily
        drifting fibre would grow without bound and no stroke would be
        enough.

        Wrapped into [0, 2*pi), so the demand never exceeds `2 * v_pi`
        (40 V on the default part) however large the phase asked for.
        """
        wrapped = phase_rad % TWO_PI
        v = self.v_pi * wrapped / math.pi
        if v > self.v_max:
            raise ValueError(
                f"phase {phase_rad:g} rad wraps to {wrapped:g} rad, needing "
                f"{v:g} V against a {self.v_max:g} V limit. The part cannot "
                f"reach a full fringe: v_pi = {self.v_pi:g} V requires "
                f"{2 * self.v_pi:g} V for 2*pi.")
        return v

    def delivers(self, phase_rad):
        """Whether this device can reach `phase_rad` at all, after wrap."""
        return self.v_pi * (phase_rad % TWO_PI) / math.pi <= self.v_max

    def apply(self, E, v):
        """Apply the device to a field: its phase, and its loss if enabled.

        Parameters
        ----------
        E : ndarray (N, 2) — complex envelope [Ex, Ey].
        v : float — drive voltage.

        Returns
        -------
        ndarray (N, 2) — field after the stretcher.
        """
        out = np.asarray(E, dtype=complex) * np.exp(1j * self.phase_for(v))
        if self.apply_insertion_loss:
            out = voa(out, self.insertion_loss_db)
        return out

    def __repr__(self):
        return (f"PiezoFibreStretcher(v_pi={self.v_pi:g} V, "
                f"v_max={self.v_max:g} V, "
                f"stroke={self.stroke_rad / math.pi:g} pi, "
                f"f_res={self.resonant_hz / 1e3:g} kHz, "
                f"loss={self.insertion_loss_db:g} dB"
                f"{'' if self.apply_insertion_loss else ' [suppressed]'})")
