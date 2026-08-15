"""Split-step time-domain (SS-TDM) dynamic model of a DFB laser diode.

Implements Kim, Chung & Lee, "An Efficient Split-Step Time-Domain Dynamic
Modeling of DFB/DBR Laser Diodes", IEEE J. Quantum Electron. 36(7), 787-794
(2000).  Equation numbers below refer to that paper.

Units are SI converted. The paper's Table I uses CGS system.

One transverse mode, not two
----------------------------
The paper models a single transverse mode.  Its (1) is

    E(x, y, z, t) = phi(x,y) * [F(z,t) e^{-i beta_0 z}
                                + R(z,t) e^{+i beta_0 z}] * e^{i omega_0 t}

with one modal function ``phi(x,y)`` and one forward/backward amplitude
pair, and (7)'s carrier equation has one photon density.  That mode is TE:
in a 0.2 um active layer the TM confinement factor is far lower, so TM does
not reach threshold, and the paper states a single ``Gamma`` accordingly.

``SimResult`` still reports ``(n_rec, 2)`` fields for the project's
``[Ex, Ey]`` convention, with the TE amplitude in column 0 and column 1
identically zero.  Orienting that axis against lab x and y is a separate
question -- mounting and pigtail alignment, not device physics -- and lives
in ``LaserDriver.sample_field`` as a Jones vector, the same way
``CWLaser`` handles it.

Physics summary
---------------
- Counter-propagating field envelopes F/R for the TE mode, advanced
  section-by-section through the split-step operator (paper
  (10), (14), (16)): the constant coupling matrix
  ``K = [[sech(gamma*dz), j*tanh(gamma*dz)], [j*tanh(gamma*dz),
  sech(gamma*dz)]]`` with ``gamma = sqrt(kappa*conj(kappa)) = |kappa|``
  (paper (16); unconditionally stable, paper (17)), followed by the
  per-section gain/detuning scalar ``exp((G - j*delta)*dz)`` (paper
  (14))::

      F(z+dz) = exp((G - j*delta)*dz) * (sech(gamma*dz) * F(z)
                                         + j*tanh(gamma*dz) * R(z+dz))
- Field gain, paper (3)::

      G = (Gamma * g_N * (N - N_0)) / (2 * (1 + eps * S)) - alpha/2

  where ``alpha`` is the waveguide loss (paper: "alpha is the waveguide
  loss due to scattering and absorption").
- Detuning factor, paper (4) and (5): ``delta = (omega_0/c)*(n_eff_0 +
  del_n) - pi/bragg_condition`` with ``omega_0/c = 2*pi/lambda`` and the
  Bragg condition ``bragg_condition = Lambda`` the grating period
  ``Lambda = lambda/(2*n_eff_0)`` (first order).  The carrier-induced
  index change is derived from the linewidth enhancement factor
  ``alpha_m`` (paper (5))::

      del_n = -(lambda / (4*pi)) * Gamma * alpha_m * g_N * (N - N_0)
      delta = (2*pi/lambda) * (n_eff_0 + del_n) - pi/Lambda

  so that ``delta = 0`` exactly at transparency (N = N_0), matching the
  paper's statement that the reference index is the one "when the carrier
  density reaches the transparency".  ``n_eff_0 = 3.283`` is the paper's
  index (Fig. 8 caption: "the index without heating is assumed to be
  3.283").  Both ``del_n`` and ``delta`` are evaluated per section, so the
  real part of the index is spatially resolved (spatial hole burning in
  the index, the paper's stated advantage).
- Carrier rate equation, paper (7), solved per section: current-density
  injection ``J/(q*d_act)`` with ``J = i_k/(w*L)`` (uniform along the
  cavity) minus bimolecular and Auger recombination minus the
  stimulated recombination of the lasing mode.  The stimulated term
  uses the modal gain ``Gamma*g_N*(N-N_0)/(1+eps*S)`` -- the same gain
  that appears in the field equation (3) -- so energy is conserved
  between the fields and the carriers.  There is no linear (A N / tau)
  term, matching the paper's Eq. (7).
- Spontaneous emission, paper (2a)/(2b): F~ and R~ "represent the
  spontaneous emission noise, which operate as the driving sources for
  oscillation"; "the amplitude distribution ... is Gaussian and the
  phases ... are assumed to change randomly" (p. 788); "the spontaneous
  noise is generated from a Gaussian distributed random number generator"
  (p. 791).  Implemented per section per step as a seeded complex Gaussian
  (Rayleigh amplitude, uniform phase) whose mean square adds half of the
  modal spontaneous rate ``beta * B * N**2`` per direction, converted to
  the field-power normalisation (``|F|**2`` is optical power in W here):
  ::

      dP_sp  = 0.5 * beta * B * N**2 * dt * v_g * A_act * E_photon
      Fn     = Fn + sqrt(dP_sp) * (g1 + j*g2)/sqrt(2)

- Facets: amplitude reflectivity ``ar`` on both facets (AR coated,
  consistent with the paper's AR-coated DFB devices; the paper quotes
  facet reflectivities only for the DBR device, 81% / 0.1%).

Deviations from the paper, and why
----------------------------------
- The finite gain-bandwidth IIR filter (paper Sec. III, refs [11], [14])
  is not implemented yet.  With a ~40 nm (5 THz) gain bandwidth at
  1550 nm, the filter is wider than the simulation band and is a
  near-null on single-mode operation; deferred.
- ``tau`` (10 ns) is defined but unused: the paper's Eq. (7) has no
  linear recombination term.
- m-th order grating (``laser_order > 1``): the grating period scales
  with ``m`` (Lambda = m*lambda/(2*n_eff_0)) but the detuning uses the
  first-order reference ``pi/Lambda``; the paper's ``beta_B = m*pi/Lambda``
  equals ``pi/Lambda`` for the configured first order (m = 1).  The
  coupling amplitude ``kappa`` is kept at the first-order value.
"""

from dataclasses import dataclass
from typing import Callable
import warnings

import numpy as np


@dataclass
class SimResult:
    """Recorded output of `Laser.simulate`.

    ``t`` is the time at the end of each recorded step; ``i`` the drive
    current applied during that step; ``P_right``/``P_left`` the optical
    power leaving the right/left facet; ``E_right``/``E_left`` the complex
    envelope there as ``[Ex, Ey]``, so ``sum(|E|**2, axis=1)`` reproduces
    the matching power column exactly, in Watts per the project field
    convention.

    The device lases one TE mode, so column 0 carries it and column 1 is
    identically zero.  The ``(n_rec, 2)`` shape is kept because it is the
    convention every other component here uses (CWLaser, MZM, fibre,
    AMZI), so the output feeds them without a conversion step.  Column 1
    is zero because TM does not reach threshold, not because a second mode
    was left unmodelled.

    Putting the TE amplitude on lab x and y axes is an orientation
    question rather than a device one; ``LaserDriver.sample_field`` does
    it with a Jones vector.
    """

    t: np.ndarray            # (n_rec,) seconds
    i: np.ndarray            # (n_rec,) amperes, >= 0
    P_right: np.ndarray      # (n_rec,) watts
    P_left: np.ndarray       # (n_rec,) watts
    E_right: np.ndarray      # (n_rec, 2) complex [Ex, Ey]
    E_left: np.ndarray       # (n_rec, 2) complex [Ex, Ey]


class DFBLaser:
    """SS-TDM DFB laser diode (Kim, Chung & Lee 2000)."""

    def __init__(self,
                 laser_order: int = 1,                 # m-th order laser
                 grating_length: float = 600e-6,       # Length of DFB laser grating
                 n_sections: int = 15,                 # Number of sections DFB is divided into
                 wavelength: float = 1.55e-6,          # Optical Lasing Wavelength
                 i_bias: float = 100e-3,               # Bias operating current of the Laser
                 run_time: float = 5e-9,               # Laser observation time or running time
                 seed: int | None = None,              # RNG seed for spontaneous emission
                 ):
        self.grating_length = grating_length
        self.n = n_sections
        self.wavelength = wavelength
        self.i_bias = i_bias
        self.run_time = run_time
        self.seed = seed
        # Parameter List:
        self.w_waveguide = 2e-6
        self.confinement = 0.3                          # TE confinement factor; TM is far lower and does not lase
        self.B = 1e-10 * 1e-6                          # Spontaneous Recombination, cm^3/s into m^3/s
        self.C = 0.75e-28 * 1e-12                       # Auger Carrier cm^6/s into m^6/s
        self.g_N = 2.5e-16 * 1e-4                       # Differential Gain, cm^2 into m^2 (2.5e-20 m^2)
        self.tau = 10e-9                                # Carrier Lifetime, in seconds (unused; paper (7) has no linear term)
        self.N_0 = 1.8e18 * 1e6                         # Transparency Carrier cm^-3 into m^-3
        self.n_eff_0 = 3.283                            # Effective Phase refractive index at transparency (paper Fig. 8)
        self.n_g = 3.7                                  # Effective Group Refractive Index
        self.d_act = 0.2e-6                             # Thickness of active layer
        self.A_act = self.w_waveguide * self.d_act      # Area of the Active Layer
        self.beta = 0.5e-4                              # Spontaneous Coupling Factor
        self.alpha = 40 * 100                           # Waveguide Loss cm^-1 into m^-1 (paper (3): -alpha/2)
        self.alpha_m = 5                                # Linewidth Enhancement Factor (paper (5))
        self.epsilon = 2e-17 * 1e-6                     # Nonlinear Saturation Coefficient cm^3 into m^3
        self.c = 3e8                                    # Speed of Light
        self.v_g = self.c / self.n_g                    # Group velocity from the group refractive index
        self.q = 1.602e-19                              # Electron charge, C
        self.dz = self.grating_length / self.n
        self.dt = self.dz / self.v_g                    # dz = v_g*dt
        self.m_order = laser_order
        self.ar = 0.01                                  # Facet amplitude reflectivity (AR coated)

        # Constants
        self.h = 6.626e-34
        self.nu = self.c / self.wavelength
        self.E_photon = self.h * self.nu

        # Bragg condition (paper (1)): the grating period Lambda.  For
        # first order the reference propagation constant is pi/Lambda.
        self.bragg_condition = (self.m_order * self.wavelength) / (2 * self.n_eff_0)  # Lambda

        self.kappa = 50 * 100                           # Coupling Coefficient cm^-1 into m^-1

        self.kappa_dz = abs(self.kappa) * self.dz
        if self.kappa_dz > 0.2 * (1.0 + 1e-9):
            warnings.warn(
                f"kappa*dz = {self.kappa_dz:.3f} exceeds the 0.2 convergence "
                f"limit (paper Fig. 5); results will lose accuracy. Raise "
                f"n_sections above {int(np.ceil(abs(self.kappa) * self.grating_length / 0.2))} "
                f"or shorten grating_length.",
                UserWarning, stacklevel=2,
            )

        # Coupling entries are constant in kappa and dz, so they are built
        # once here rather than per step (see _coupling).
        self._c00, self._c01 = self._coupling()

        # Field and carrier state (carriers at the CENTERS of the sections).
        # F/R are the forward/backward amplitudes of the one TE mode.
        self.F = np.zeros(self.n + 1, dtype=complex)
        self.R = np.zeros(self.n + 1, dtype=complex)
        self.N = np.ones(self.n) * 1e24                 # Carrier density, m^-3
        self.reset()

    def __repr__(self) -> str:
        return (
            f"DFBLaser(L={self.grating_length * 1e6:.0f} um, n={self.n}, "
            f"lambda={self.wavelength * 1e9:.1f} nm, i_bias={self.i_bias * 1e3:.0f} mA, "
            f"run_time={self.run_time * 1e9:.1f} ns)"
        )

    def reset(self, seed: int | None = None) -> None:
        """Restore the initial state (fields zero, carriers at 1e24 m^-3).

        Re-seeds the spontaneous-emission RNG with ``seed`` if given,
        otherwise with the constructor seed.  With ``seed=None`` runs are
        non-deterministic.
        """
        if seed is not None:
            self.seed = seed
        self._rng = np.random.default_rng(self.seed)
        self.F[:] = 0.0
        self.R[:] = 0.0
        self.N[:] = 1e24

    def _coupling(self) -> tuple[float, complex]:
        """Constant coupling entries of the split-step operator (paper (16)).

        ``K = [[sech(gamma*dz), j*tanh(gamma*dz)],
               [j*tanh(gamma*dz), sech(gamma*dz)]]`` with
        ``gamma = sqrt(kappa*conj(kappa)) = |kappa|`` (5000 m^-1 for the
        index grating).  The matrix is the same for every section -- the
        per-section physics (gain, detuning) enters through the scalar
        ``exp((G - j*delta)*dz)`` multiplier in the sweep.  Successive
        multiplication with ``K`` is unconditionally stable (paper (17)):
        in the passive stopband it conserves power per section
        (``|sech|^2 + |j*tanh|^2 = 1``), where the bare full-flow
        ``expm`` entries (``cosh``/``sinh``) would amplify the evanescent
        mode and blow up a passive device.
        """
        gamma = np.sqrt(self.kappa * np.conj(self.kappa))
        x = gamma * self.dz
        return 1 / np.cosh(x), 1j * np.tanh(x)

    def simulate(self,
                 current: float | Callable[[float], float] | None = None,
                 t_end: float | None = None,
                 record_every: int = 10) -> SimResult:
        """Run the device for ``t_end`` at the injection ``current``.

        ``current`` is either a scalar (constant bias, defaults to
        ``i_bias``) or a callable ``i(t)`` evaluated at the end of every
        step (for the laser driver's current waveforms).  The current is
        clamped to be non-negative.

        Records every ``record_every``-th step (default 10; pass 1 for a
        dense record of every step): time, applied current, right/left
        facet power and complex envelope.

        The device takes roughly 30 ns to settle from its zero-field
        initial state, so a run shorter than that is measuring turn-on
        rather than steady operation.  ``LaserDriver.sample_field``
        discards a settle window for this reason.
        """
        if current is None:
            current = self.i_bias
        if t_end is None:
            t_end = self.run_time
        record_every = max(1, int(record_every))
        self.reset()

        n_steps = max(1, int(t_end / self.dt))
        n_rec = n_steps // record_every
        t_rec = np.empty(n_rec)
        i_rec = np.empty(n_rec)
        p_right = np.empty(n_rec)
        p_left = np.empty(n_rec)
        # Column 1 stays zero: the device lases TE only, and zeros() rather
        # than empty() is what makes that true without writing it each step.
        e_right = np.zeros((n_rec, 2), dtype=complex)
        e_left = np.zeros((n_rec, 2), dtype=complex)
        rec = 0

        # The coupling entries depend only on kappa and dz, both fixed at
        # construction, so they are read once here instead of being
        # recomputed (sqrt, cosh, tanh) on every one of n_steps iterations.
        c00, c01 = self._c00, self._c01

        # Double buffers for the sweep.  Every element of both arrays is
        # written each step -- F[0] from the left facet and F[1..n] from the
        # loop, R[0..n-1] from the loop and R[n] from the right facet -- so
        # the buffers can be reused rather than reallocated, as long as they
        # are SWAPPED with the state arrays at the end of the step rather
        # than assigned to them.  Assigning would alias the state and the
        # write target, and the next sweep would read values it had already
        # overwritten.
        F_next = np.empty(self.n + 1, dtype=complex)
        R_next = np.empty(self.n + 1, dtype=complex)

        for k in range(n_steps):
            t_now = (k + 1) * self.dt
            i_k = current(t_now) if callable(current) else float(current)
            i_k = max(i_k, 0.0)

            # Center power by averaging at the section boundaries
            P_center = 0.5 * ((np.abs(self.F[:-1]) ** 2 + np.abs(self.F[1:]) ** 2) +
                              (np.abs(self.R[:-1]) ** 2 + np.abs(self.R[1:]) ** 2))

            # Optical photon density from optical power
            S = P_center / (self.v_g * self.A_act * self.E_photon)

            # Field gain (paper (3)); waveguide loss outside the fraction
            g = (self.confinement * self.g_N * (self.N - self.N_0)) / (2 * (1 + self.epsilon * S)) - self.alpha / 2

            # Carrier-induced index change (paper (5)) and detuning (paper (4))
            del_n = -(self.wavelength / (4 * np.pi)) * self.confinement * self.alpha_m * self.g_N * (self.N - self.N_0)
            delta = (2 * np.pi / self.wavelength) * (self.n_eff_0 + del_n) - np.pi / self.bragg_condition

            # Spontaneous emission into the mode (Gaussian amplitude, random phase).
            # Half of beta*B*N^2 per direction, converted from photon density
            # (m^-3/s) to the field-power normalisation (|F|^2 in W).
            noise_amp = np.sqrt(0.5 * self.beta * self.B * self.N ** 2 * self.dt *
                                self.v_g * self.A_act * self.E_photon)
            xi = (self._rng.standard_normal((self.n, 2)) +
                  1j * self._rng.standard_normal((self.n, 2))) / np.sqrt(2)

            phase = np.exp((g - 1j * delta) * self.dz)

            # Split-step operator (paper (10), (14), (16)): the constant
            # coupling matrix K = [[sech, j*tanh], [j*tanh, sech]] times
            # the per-section gain/detuning scalar exp((G - j*delta)*dz),
            # transferred section by section.  The sweep is SEQUENTIAL
            # (paper: "transferring the fields section by section"): F is
            # advanced left-to-right using the just-updated F (facets on
            # the left), while R enters each section from the right with
            # its pre-sweep value.  A parallel (snapshot) application is
            # NOT equivalent and breaks energy conservation in the
            # stopband.  Boundary values are read into locals once per
            # section, which is worth doing because interpreted loop
            # overhead dominates this function.
            # Pre-sweep backward field; the swap below is what keeps this
            # distinct from the array being written.
            R_prev = self.R
            F_next[0] = self.ar * R_prev[0]
            for i in range(self.n):
                p = phase[i]
                f, r = F_next[i], R_prev[i + 1]
                F_next[i + 1] = p * (c00 * f + c01 * r)
                R_next[i] = p * (c01 * f + c00 * r)

            # Spontaneous emission driving sources (paper (2a)/(2b))
            F_next[1:] += noise_amp * xi[:, 0]
            R_next[:-1] += noise_amp * xi[:, 1]

            # Facet boundary conditions (AR coated).  Left: already applied
            # at the start of the sweep (F_next[0] = ar*R[0]).  Right:
            # the backward wave re-enters from the right facet.
            R_next[self.n] = self.ar * F_next[self.n]

            # Carrier rate equation (paper (7)): current-density injection
            # J/(q*d_act) with J = i_k/(w*L) (uniform along the cavity),
            # recombination terms, and stimulated recombination of the
            # lasing mode with the same modal gain as the field equation (3).
            J = i_k / (self.w_waveguide * self.grating_length)
            inj = J / (self.q * self.d_act)
            stim = (self.v_g * self.confinement * self.g_N * (self.N - self.N_0) * S) / (1 + self.epsilon * S)
            dN_dt = inj - self.B * self.N ** 2 - self.C * self.N ** 3 - stim

            self.N += dN_dt * self.dt
            self.F, F_next = F_next, self.F
            self.R, R_next = R_next, self.R

            if (k + 1) % record_every == 0:
                t_rec[rec] = t_now
                i_rec[rec] = i_k
                p_right[rec] = np.abs(self.F[-1]) ** 2
                p_left[rec] = np.abs(self.R[0]) ** 2
                e_right[rec, 0] = self.F[-1]
                e_left[rec, 0] = self.R[0]
                rec += 1

        return SimResult(t=t_rec, i=i_rec, P_right=p_right, P_left=p_left,
                         E_right=e_right, E_left=e_left)
