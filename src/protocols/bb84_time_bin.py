"""Time-bin phase-encoding BB84 QKD protocol.

Signal chain:
  Pulsed laser → AsymmetricMZI(encoder) → PhaseModulator (φ_A)
    → propagate() (fiber) → AsymmetricMZI(decoder) → 2× SPAD

Based on Gobby, Yuan & Shields (2004), Appl. Phys. Lett. 84, 3762.

Performance: the interference closed form
------------------------------------------------------------
**This is why a 1e8-pulse sweep is feasible at all.  Preserve it.**

The field chain is deterministic given the three
discrete choices (alice_basis, alice_bit, bob_basis), so 8 precomputed
outcomes replace one field propagation per pulse — a 2,000,000-pulse run
had been doing 4,000,000 propagations to obtain 8 distinct answers.

That form is exact only while the response is deterministic.  A per-pulse
random phase — laser linewidth against a residual path mismatch, or
modulator drive noise — makes it continuous, and 8 answers no longer
cover it.

**The same linearity argument extends, and the extension is exact rather
than an approximation.**  The chain is linear, so the two interfering
paths superpose and the gated power is exactly

    P(delta) = g0 + 2*Re(S * exp(i*delta)),    delta = phi_A - phi_B,

for constants (g0, S) fixed by the point's parameters.  Three evaluations
at delta = 0, pi/2, pi determine them:

    g0    = (P(0) + P(pi)) / 2
    Re(S) = (P(0) - P(pi)) / 4
    Im(S) = (g0 - P(pi/2)) / 2

So the cost is three field propagations per *point* — not per pulse — and
any phase thereafter costs two multiplies.  This generalises from
"8 fixed answers" to "one closed form, any phase": strictly more capable
at the same per-pulse cost.

Two properties worth not losing in a later refactor:

* The coefficients are extracted *from the field chain itself* rather than
  re-derived by hand.  There is therefore no second expression of the
  physics that can drift out of step with the components, and the form is
  correct for both interferometer topologies without special-casing.
* Because they are extracted by sweeping phi_A at phi_B = 0, they carry no
  information about how phi_B enters, and the per-pulse formula *assumes*
  delta = phi_A - phi_B.  A consistency assertion after the extraction
  checks that assumption against the chain at phi_B != 0, so the field
  chain stays authoritative — it fires on exactly the encoder
  phase_arm / Bob sign coupling that had to be verified by
  hand.

Verified bit-identical to the 8-entry table at 0/65/122 km (max absolute
error 1.8e-16 of full scale, i.e. the float64 floor), with negative
controls: a spurious 1 mrad phase offset shows up at 5e-4, and a flipped
Bob sign trips the consistency assertion.

References
----------
[1] Gobby, C., Yuan, Z. L., & Shields, A. J. (2004). Quantum key
    distribution over 122 km of standard telecom fiber. Appl. Phys.
    Lett. 84(19), 3762-3764.
[2] Bennett, C. H. & Brassard, G. (1984). Quantum cryptography:
    Public key distribution and coin tossing. Proc. IEEE Int. Conf.
    on Computers, Systems and Signal Processing, 175-179.
[3] Townsend, P. D. et al. (1993). Single photon interference in
    10 km long optical fibre interferometer. Electron. Lett. 29(7),
    634-635.
"""
import argparse
import cmath
import math
import numpy as np
import random
import sys

sys.path.insert(0, 'src')
from src.channel.interferometer import AsymmetricMZI
from src.channel.phase_modulator import PhaseModulator
from src.channel.optics import pbs, pbc, voa
from src.channel.fiber import FiberRealization
from src.channel.piezo_stretcher import PiezoFibreStretcher
from src.detectors.spad import spad


def gaussian_pulse(t, sigma, A=1.0):
    """Gaussian pulse envelope, unit peak amplitude."""
    return A * np.exp(-t ** 2 / (2 * sigma ** 2))


def field_grid(pulse_width=100e-12, delay=5.8e-9):
    """The time grid this chain builds internally, for a caller supplying
    its own source field.

    Returns ``(dt, n_samples, pulse_center)``.

    A source field must arrive already on this grid.  The chain does not
    resample one, deliberately: the DFB device step is 0.49 ps and this
    grid is typically 10 ps, so resampling here would alias unless it
    band-limited first -- and `LaserDriver.sample_field` already
    band-limits correctly by averaging.  Asking the driver for `dt`
    directly is both simpler and right, where a resample buried in the
    protocol would be neither.

    Note what that costs: at a 10 ps grid the Nyquist is 50 GHz and the
    gain-switched chirp is hundreds, so most of the chirp is averaged
    away.  For this chain that is not a loss, because a
    path-matched interferometer compares the pulse with a copy of itself
    at the same chirp phase.  It matters enormously in the polarisation
    chain, which is the contrast worth having.
    """
    dt = min(pulse_width / 10.0, delay / 20.0)
    n_samples = int(np.ceil((2.0 * delay + 5.0 * pulse_width) / dt))
    return dt, n_samples, delay / 2.0


def simulate_bb84_time_bin(num_bits, fiber_length=0, alpha_dB=0.182,
                            mu=0.1, wavelength=1550e-9,
                            pulse_width=100e-12, repetition_rate=2.5e6,
                            delay=5.8e-9, gate_width=1e-9,
                            spad_eta=0.10, dark_count_rate=15.0,
                            afterpulse_prob=0.05, dead_time=13e-6,
                            visibility=1.0, phase_error=0.0,
                            interferometer='balanced', split_ratio=1.6,
                            linewidth=0.0, path_mismatch=0.0,
                            phase_error_rad=0.0, phase_noise_rad=0.0,
                            bias_offset_v=0.0, phase_drift_rad_s=0.0,
                            run_duration=None, source_field=None,
                            birefringence=False, cd=False, pmd=False,
                            temperature=25.0, bend_radius=None,
                            pmd_coeff_ps_sqrt_km=0.1, compensate=True,
                            drift_temperature_rate_C_s=0.0, drift_blocks=100,
                            phase_servo_interval_s=None,
                            phase_servo_stretcher=None,
                            seed=None, verbose=False):
    """BB84 time-bin phase-encoding simulation.

    Parameters
    ----------
    num_bits : int — number of BB84 pulses to simulate.
    fiber_length : float — fiber length in km (default 0).
    alpha_dB : float — fiber attenuation in dB/km (default 0.182).
    mu : float — mean photon number per pulse (default 0.1).
    wavelength : float — centre wavelength (m, default 1550e-9).
    pulse_width : float — FWHM of Gaussian pulse (s, default 100e-12).
    repetition_rate : float — pulse repetition rate (Hz, default 2.5e6).
    delay : float — AMZI differential delay (s, default 5.8e-9).
    gate_width : float — SPAD gate window (s, default 1e-9).
    spad_eta : float — SPAD quantum efficiency (default 0.10).
    dark_count_rate : float — SPAD dark count rate (Hz, default 15.0).
    afterpulse_prob : float — afterpulse probability (default 0.05).
    dead_time : float — SPAD dead time (s, default 13e-6).
    visibility : float — decoder interferometric visibility in (0, 1],
        default 1.0 (ideal).  Optical misalignment error e_opt = (1 - V)/2.
    phase_error : float — static decoder phase offset (radians),
        default 0.0; represents imperfect AMZI path-length matching.
    interferometer : {'balanced', 'polarisation_multiplexed'}
        Which AMZI pair to use (default 'balanced', unchanged behaviour).

        'balanced' is the generic 50:50 `AsymmetricMZI`.  It produces four
        paths, two of which are non-interfering satellite bins at t = 0
        and t = 2*delay carrying half the launched energy, which the gate
        then discards -- so only mu/2 reaches the detectors.

        'polarisation_multiplexed' is `PolarizationMultiplexedAMZI`, the
        apparatus Gobby et al. actually describe: a polarising beam
        combiner/splitter routes Alice-short -> Bob-long and
        Alice-long -> Bob-short, so only two paths exist, there are no
        satellites, and the full mu reaches the single interference peak.
        Use this for the Gobby replication.
    split_ratio : float — Alice's reference:encoded intensity ratio,
        used only by 'polarisation_multiplexed' (default 1.6, per Gobby).
    linewidth : float — laser linewidth (Hz, default 0 = off).
    path_mismatch : float — RESIDUAL delay mismatch between the S-L and
        L-S routes after the delay line and fibre stretcher (s, default 0).

        Linewidth couples to the QBER *only* through this residual, as
        `sigma = sqrt(2*pi*linewidth*path_mismatch)`.  It does **not**
        couple through the full AMZI delay: the two interfering routes
        traverse the same total path, so the frequency-noise term cancels.
        That cancellation is why the delay line and stretcher exist, and
        why >99% fringe visibility is achievable from an 80 ps pulsed
        source whose coherence time is far shorter than the 5.8 ns delay.

        Both default to 0 so existing behaviour is unchanged.  1550 nm DFB
        diodes run from several hundred kHz to ~10 MHz (current parts cite
        2-3.2 MHz); Gobby state no value, so any choice here is a
        documented assumption, never a derived or fitted quantity.  At
        realistic trim the contribution is <0.02% — see the contribution
        budget.
    phase_error_rad : float — STATIC phase offset applied to every pulse
        (rad, default 0).  Models a modulator calibration offset.
    phase_noise_rad : float — per-pulse Gaussian phase jitter, sigma in
        radians (default 0).  Models modulator drive noise.  Adds in
        quadrature with the linewidth term above, both being zero-mean
        Gaussians entering the same relative phase.

        Drive electronics typically deliver 0.1-1% of V_pi, so this is a
        small term unless the modulator is unusually noisy.  A bias set
        once and left is the static mechanism instead; see
        `phase_error_rad` and `bias_offset_v`.
    bias_offset_v : float — the same static bias error as
        `phase_error_rad`, expressed as a modulator drive voltage (default
        0).  Converted by `PhaseModulator` through its crystal-derived
        V_pi; supplying both raises there.  This is the mechanism Gobby
        et al. (2004) name first, "slight inaccuracies of the phase
        modulator biases".
    phase_drift_rad_s : float — rate (rad/s) at which the decoder's arm
        imbalance drifts, default 0 (no drift).  This is *interferometer*
        arm-length drift, a distinct mechanism from modulator bias, and it
        is owned by `AsymmetricMZI.arm_phase_offset` — this function only
        supplies the clock.  Gobby measure < 0.05 deg/s = 8.727e-4 rad/s
        with both setups enclosed against air convection.

        It matters only for long runs.  Accumulated phase is rate * t, so
        a 2-minute key transfer (their figure) reaches 6.0 deg and 0.091%
        QBER, while a 3e6-pulse run at 2 MHz lasts 1.5 s and contributes
        essentially nothing.
    run_duration : float or None — physical duration of the experiment
        being simulated, in seconds (default None).  Sets the *drift*
        clock only; the detector clock always keeps true 1/repetition_rate
        spacing, because dead time and afterpulsing are defined against
        real elapsed time.

        None ties drift to the pulse budget (`num_bits / repetition_rate`),
        which is the historical behaviour and is bit-identical to it.
        **Prefer setting it.**  The pulse budget is chosen for statistical
        power, so leaving drift tied to it means asking for tighter error
        bars silently lengthens the simulated experiment.  That is what
        made the Gobby 122 km point read 13.52% against a stated 8.9%: it
        needed 1e9 pulses = 500 s at 2 MHz, accumulating 25 deg of drift
        where their two-minute transfer accumulates 6.

        With it set, the simulated pulses are a uniform sample across a run
        of that length.  The count may exceed what the apparatus actually
        sent in that time — that is the point, and what makes it a
        lower-variance estimate of the same expectation rather than a
        longer run.
    source_field : ndarray (n_samples, 2) complex, or None — the optical
        field entering the encoder, in place of the analytic Gaussian.
        Default None keeps the Gaussian this chain has always built.

        The array must already be on this chain's grid; call
        `field_grid(pulse_width, delay)` for `(dt, n_samples,
        pulse_center)`.  Nothing is resampled here on purpose -- see that
        function.  It is renormalised to the same `energy_per_pulse`, so mu
        stays calibrated once at the source and the swap changes the
        statistics of the light rather than its level.

        Measured: replacing the Gaussian with a gain-switched DFB
        leaves QBER alone across 0-122 km.  A path-matched AMZI interferes
        the pulse with a copy of itself at the same chirp phase, so a
        common chirp cancels, which is the linewidth argument reached
        from another direction.  Contrast the polarisation chain, where
        the same source's chirp
        turns PMD from nothing into +9.6 pp in the polarisation chain.

        Per-pulse ENERGY spread is not covered.  `bb84_duplinskiy` takes
        `pulse_energy_factors` because energy was the one
        thing that mattered there; there is no equivalent here yet.
    birefringence, cd, pmd : bool — fibre impairments, all default False.
        Names and meanings match `bb84_duplinskiy`.  Attenuation is NOT
        among them: it stays the scalar `alpha_dB` factor in every case, so
        with all three off the field chain is untouched and results are
        bit-identical to a run predating them.
    temperature : float — ambient temperature in C (default 25.0).
    bend_radius : float or None — bend radius in metres (default None =
        unbent).
    pmd_coeff_ps_sqrt_km : float — PMD coefficient in ps/sqrt(km) (default
        0.1, the Corning SMF-28 Ultra spec).
    compensate : bool — whether Bob aligns to the fibre (default True).
        Gobby et al. [1] describe "careful alignment of the polarisation
        maintaining optics" and measure 99.96 % classical fringe visibility
        over the full 122 km link, so the apparatus is an aligned one.
        Alignment is the channel matrix's conjugate transpose and removes a
        STATIC fibre exactly -- U_comp @ J = I for any unitary -- so a null
        measured with compensation on and nothing drifting is arithmetic
        rather than physics.  False gives the uncompensated control.
    drift_temperature_rate_C_s : float — rate (C/s) at which the fibre's
        ambient temperature changes during the run, default 0.0 (static).
        Requires `birefringence=True`, the only operator it acts through.

        This is the FIBRE's drift, distinct from the interferometer
        arm-length drift on `phase_drift_rad_s`.  In the
        polarisation-multiplexed topology the two are degenerate in their
        effect on QBER: the Jones matrix is SU(2), so the fibre contributes
        a common amplitude |U00|^2 and a relative phase 2*arg(U11), and
        that phase enters `delta` exactly as an arm-length offset would.
        They cannot be separated from QBER data alone.

        Off by default in the Gobby replication for that reason -- see
        PHASE_DRIFT_RAD_S in `validate_gobby.py`, where the paper's stated
        3.3 % floor is already assigned in full to modulator bias plus arm
        drift, and adding a third term would count it twice.
    drift_blocks : int — how many times the fibre is re-evaluated across the
        run while it drifts (default 100).  A COUNT rather than a pulse
        size, so raising `num_bits` for statistical power leaves the drift
        resolution alone.  Ignored when nothing drifts.
    phase_servo_interval_s : float or None — how often Bob re-locks his
        phase, in seconds.  None (default) means no servo.

        This is the piezo-driven fibre stretcher in Bob's long arm, which
        Gobby et al. [1] use to hold the operating point.  Without it the
        model contradicts the paper: a residual rotation is SU(2), so it
        costs O(eps^2) in amplitude against O(eps) in phase, and the QBER
        therefore moves before the rate does -- the reverse of
        "polarisation drift reduces the bit rate, but does not degrade the
        QBER".

        PHASE ONLY, which is the entire point.  Inverting the full Jones
        matrix would null the amplitude as well and the rate loss the
        paper reports would vanish with it.  Measured at 10 km, 3e-3 C/s
        over the paper's 120 s transfer: the sifted rate falls to 0.589 of
        reference with the servo on and 0.589 with it off, while the QBER
        goes from 41.6 % to 0.054 % against a 0.095 % baseline.

        Sample and hold, because a real servo has finite bandwidth.  The
        lock refreshes when it is older than this interval and the fibre
        keeps moving in between, so a short interval reproduces the paper
        and a long one degrades back to the unserved case -- measured
        0.00 %, 0.39 %, 4.09 % and 14.57 % at 1.2 s, 12 s, 60 s and never.

        Not offered on `bb84_duplinskiy`: there the residual is a full
        SU(2) acting on the encoding itself rather than a scalar phase
        between two arms, so a phase-only correction has nothing to
        correct.  That chain answers drift with the full recalibration
        `calibration_temperature` already models.
    phase_servo_stretcher : PiezoFibreStretcher or None — the device Bob
        actually turns.  None (default) builds a Thorlabs FVP155P with its
        insertion loss SUPPRESSED, which is correct for this chain and
        wrong in general -- see below.

        Routing the correction through a real part is what stops it being
        an unbounded number.  `voltage_for` wraps into one fringe first, so
        the demand never exceeds 2*v_pi (40 V on the default part) however
        far the fibre has drifted, comfortably inside the 150 V limit.  A
        demand the device cannot meet raises rather than clipping.

        **The suppressed loss is a link-budget decision, not a device
        one.**  `ETA_BOB = 0.045` in `validate_gobby.py` already folds in
        "5 dB of loss in Bob's apparatus", and the stretcher sits inside
        that apparatus, so charging its 0.1 dB again would double-count
        against a number taken from the paper.  A chain whose budget does
        not already contain Bob's optics should pass a stretcher with the
        loss left on, which is the component's own default.

        The loss lands on Bob's long arm, so it unbalances the interfering
        pair and moves visibility as well as rate.  Exaggerated to 3 dB it
        costs 23 % of the sifted rate and lifts the QBER to 4.7 %; at the
        datasheet's 0.1 dB neither is resolvable.

        NOTE ON REPRODUCIBILITY.  Turning the servo on is *statistically*
        identical with and without the device, but not bit-identical: the
        round trip through `v_pi` changes the delivered phase in its last
        bit, and this chain's RNG desynchronises on any such change,
        because a single flipped detection alters how many draws follow.
        Measured, a deliberate 1e-15 relative nudge moves the sifted count
        by about 13 % with no systematic trend -- a reshuffled sample of
        the same distribution, not a different answer.
    seed : int or None — RNG seed.
    verbose : bool — print progress.

    Returns
    -------
    dict with keys: qber, n_total, n_sifted, n_errors, sifted_key_rate, etc.
    """
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)

    h = 6.626e-34
    c = 3e8
    photon_energy = h * c / wavelength

    # Gaussian sigma from FWHM
    sigma = pulse_width / (2.0 * np.sqrt(2.0 * np.log(2.0)))

    # The grid, from the one place that defines it.  `field_grid` is what a
    # caller supplying its own source field is told to build against, and
    # the shape check below compares against this `n_samples` -- so the two
    # must be the same derivation rather than two copies of it.
    dt, n_samples, pulse_center = field_grid(pulse_width, delay)
    t = np.arange(n_samples, dtype=float) * dt

    # The pulse sits at t = delay/2 so the interference, at
    # t = delay/2 + delay, falls within the window.
    energy_per_pulse = mu * photon_energy
    # Peak field amplitude: energy = A² ∫exp(-t²/σ²) dt = A² · σ·√π (continuous)
    # In discrete: Σ|A·g|²·dt = A²·σ·√π  →  A² = energy / (σ·√π)
    pulse_amplitude = np.sqrt(energy_per_pulse / (sigma * np.sqrt(np.pi)))
    pulse = gaussian_pulse(t - pulse_center, sigma, A=pulse_amplitude)
    E_pulse = pulse[:, np.newaxis] * np.array([1.0, 0.0], dtype=complex)  # X-polarised

    # A real source, in place of the analytic Gaussian.  The field must
    # already be on this chain's grid -- `field_grid(pulse_width, delay)`
    # returns it -- because resampling here would alias: the DFB device
    # step is 0.49 ps against a 10 ps grid, and `LaserDriver.sample_field`
    # already band-limits by averaging when asked for the coarser dt.
    #
    # Renormalised to the same energy_per_pulse, so mu stays calibrated
    # once at the source and the swap changes the STATISTICS of the light
    # rather than its level.  Same convention as `bb84_duplinskiy`.
    if source_field is not None:
        E_src = np.asarray(source_field, dtype=complex)
        if E_src.shape != (n_samples, 2):
            raise ValueError(
                f"source_field must have shape ({n_samples}, 2) to match "
                f"this chain's grid, got {E_src.shape}. Build it with "
                f"field_grid(pulse_width={pulse_width!r}, delay={delay!r}), "
                f"which returns (dt, n_samples, pulse_center).")
        energy = float(np.sum(np.abs(E_src) ** 2)) * dt
        if energy <= 0:
            raise ValueError("source_field carries no energy")
        E_pulse = E_src * np.sqrt(energy_per_pulse / energy)

    # Build AMZI components.
    #
    # 'balanced'  -- one field carries both time bins, so Bob's coupler
    #   produces four paths: S-S and L-L land at t=0 and t=2*delay as
    #   non-interfering satellites carrying half the energy, which the gate
    #   discards.  Only mu/2 reaches the detectors.
    #
    # 'polarisation_multiplexed' -- Alice's arms leave on orthogonal
    #   polarisations (PBC) and Bob's PBS routes them into opposite arms,
    #   so only S-L and L-S exist.  No satellites; the full mu reaches the
    #   single interference peak.  Composed below from generic components.
    if interferometer not in ('balanced', 'polarisation_multiplexed'):
        raise ValueError(
            f"interferometer must be 'balanced' or "
            f"'polarisation_multiplexed', got {interferometer!r}")

    kappa = 1.0 / (1.0 + split_ratio)     # encoded (short) arm power share
    polmux = interferometer == 'polarisation_multiplexed'

    if polmux:
        if visibility != 1.0:
            raise ValueError(
                "visibility is not a free input for the "
                "polarisation-multiplexed topology: the arm amplitudes "
                "determine it. Injecting one would re-apply error physics "
                "the link budget already produces."
            )
        # Alice: unbalanced splitter, arms taken out separately so a PBC
        # can put them on orthogonal axes.  Her modulator sits in the
        # SHORT arm, because that is the arm carrying the encoded pulse
        # ("the encoded pulse (through Alice's short arm)" -- Gobby).  The
        # reference travels the long arm unmodulated.
        enc = AsymmetricMZI(delay=delay, mode='encoder', split_ratio=kappa,
                            phase_arm='short')
        # Bob: arms arrive pre-split from the PBS.
        dec = AsymmetricMZI(delay=delay, mode='decoder',
                            phase_error=phase_error,
                            phase_drift_rad_s=phase_drift_rad_s)
        # Equalising the arms is what buys the >0.99 visibility the source
        # reports; an unequal pair caps V at 2*sqrt(r)/(1+r) = 0.973 for
        # r = 1.6.  Attenuate the stronger (reference) arm with a real VOA.
        balance_dB = 10.0 * np.log10(split_ratio)
    else:
        enc = AsymmetricMZI(delay=delay, mode='encoder')
        dec = AsymmetricMZI(delay=delay, mode='decoder',
                            visibility=visibility, phase_error=phase_error,
                            phase_drift_rad_s=phase_drift_rad_s)
        balance_dB = 0.0

    # SPAD detectors (constructive and destructive ports)
    spad_c = spad(wavelength=wavelength, quantum_efficiency=spad_eta,
                  dead_time=dead_time, dark_count_rate=dark_count_rate,
                  afterpulse_prob=afterpulse_prob, gate_width=gate_width)
    spad_d = spad(wavelength=wavelength, quantum_efficiency=spad_eta,
                  dead_time=dead_time, dark_count_rate=dark_count_rate,
                  afterpulse_prob=afterpulse_prob, gate_width=gate_width)

    # Fiber loss factor
    fiber_loss_lin = 10.0 ** (-alpha_dB * fiber_length / 20.0)  # field factor (sqrt of power)

    # Fibre impairments, off by default.  Attenuation stays the scalar
    # factor above in every case, so with all three flags off the field
    # chain is untouched and results are bit-identical to a run without
    # this block.  `FiberRealization` therefore carries only the unitary
    # and dispersive operators here, built with attenuation=False.
    #
    # Where it goes matters.  In the polarisation-multiplexed topology the
    # fibre sits between Alice's beam combiner and Bob's splitter, which is
    # a real fibre's place in that apparatus and the only place a rotation
    # can act on the two arms at all.
    #
    # What a rotation does there is exactly two things, because the Jones
    # matrix is SU(2): |U00| = |U11| and arg(U00) = -arg(U11).  So both
    # interfering arms are scaled by the SAME |U00| -- there is no
    # imbalance available to collapse the fringe, and the sifted rate falls
    # as |U00|^2 -- while the pair acquires a relative phase 2*arg(U11),
    # which is degenerate with a modulator bias offset.
    #
    # The leaked light (reference into Bob's long arm, encoded into his
    # short) travels L-L and S-S, lands in the satellite bins two delays
    # late or not delayed at all, and is excluded by the gate.  Lost, not
    # misassigned: that is why the amplitude term is a pure rate loss.
    #
    # Both halves of Gobby et al. [1] follow -- "Polarisation drift reduces
    # the bit rate, but does not degrade the QBER provided that the signal
    # rate is significantly higher than the intrinsic error rate."  The
    # rate falls by |U00|^2; the QBER holds because the phase is calibrated
    # out, which their Bob does with the piezo-driven fibre stretcher in
    # his long arm; and the proviso is the background taking a larger share
    # of a reduced signal, the regime the 122 km point lives in.
    #
    # `validate_gobby_impairments.py` measures both closed forms.
    if drift_temperature_rate_C_s != 0.0 and not birefringence:
        raise ValueError(
            "drift_temperature_rate_C_s acts only through the birefringence "
            "Jones matrix, which is not built when birefringence=False, so "
            "this combination would drift nothing. Pass birefringence=True, "
            "or leave the rate at 0.")

    fibre = None
    if birefringence or cd or pmd:
        fibre = FiberRealization(
            L_m=fiber_length * 1000.0, wavelength=wavelength,
            temperature=temperature, bend_radius=bend_radius,
            attenuation_factor=alpha_dB,
            pmd_coeff_ps_sqrt_km=pmd_coeff_ps_sqrt_km,
            birefringence=birefringence, cd=cd, pmd=pmd,
            attenuation=False, seed=seed,
            drift_temperature_rate_C_s=drift_temperature_rate_C_s)

    # Polarisation alignment.  Gobby et al. [1] describe "careful alignment
    # of the polarisation maintaining optics", measure 99.96 % classical
    # fringe visibility over the full 122 km link, and report polarisation
    # stable for over 30 minutes there, so the apparatus is an ALIGNED one.
    # What it experiences is residual drift about that alignment, not the
    # free rotation an uncompensated fibre applies.
    #
    # That distinction is not cosmetic.  Beyond about ten metres the
    # sectional model delivers an essentially random SU(2), so 2*arg(U11)
    # is a uniformly random phase and the QBER lands wherever it lands --
    # which is neither what the paper reports nor what the apparatus does.
    #
    # Aligning is the inverse of the channel's own matrix, the same fixed
    # point `bb84_duplinskiy` uses for its calibration loop.  It removes
    # both terms exactly, U_comp @ J = I, so a null measured on this path
    # is arithmetic rather than physics.  Pass compensate=False for the
    # uncompensated control.
    U_comp = None
    if fibre is not None and compensate:
        J = fibre.birefringence_matrix()
        U_comp = None if J is None else J.conj().T

    # SPAD gate in samples (half-window for power extraction)
    gate_half_samples = max(1, int(gate_width / dt / 2.0))
    delay_samples = int(delay / dt)
    interference_idx = pulse_center_idx = int(pulse_center / dt) + delay_samples
    start_i = max(0, interference_idx - gate_half_samples)
    end_i = min(n_samples, interference_idx + gate_half_samples + 1)

    # Sifting is accumulated inline rather than into per-pulse lists.  The
    # lists were O(num_bits) in memory (five of them), which is fine at the
    # 1e6-pulse runs this was written for but reaches several GB at the
    # ~1e8 pulses a target-sifted sweep needs at 122 km.  Counting inline
    # is O(1) and touches no RNG draw, so results are bit-identical.
    n_sifted = 0
    n_errors = 0
    n_clicks = 0

    PHI_A = {'X': {0: 0.0, 1: np.pi},
             'Y': {0: np.pi / 2.0, 1: 3.0 * np.pi / 2.0}}
    PHI_B = {'X': 0.0, 'Y': np.pi / 2.0}

    def _gate_powers(phi_a, phi_b, op_fibre):
        """Run the field chain once and return the gated (P_c, P_d).

        `op_fibre` is the fibre as it is at the instant being evaluated --
        `fibre` itself when nothing drifts, `fibre.at(t)` when it does.
        Bob's `U_comp` is deliberately NOT re-derived from it: his
        compensator was aligned once against the fibre as it was, and the
        residual `U_comp @ J(t)` is the entire physics of drift.
        """
        if polmux:
            # Alice: unbalanced split, arms kept separate.  `enc` is
            # built with phase_arm='short', so phi_A lands on the
            # encoded pulse and the reference stays unmodulated.
            E_enc, E_ref = enc.modulate(E_pulse, dt, phase=phi_a,
                                        recombine=False)
            # Polarising beam combiner: reference -> H, encoded -> V.
            E_field = pbc(E_ref[:, 0], E_enc[:, 0])
        else:
            E_field = enc.modulate(E_pulse, dt, phase=phi_a)

        E_field = E_field * fiber_loss_lin
        if op_fibre is not None:
            E_field = op_fibre.apply(E_field, dt=dt)
            if U_comp is not None:
                E_field = np.transpose(U_comp @ np.transpose(E_field))

        if polmux:
            # Bob's PBS routes by polarisation: Alice's reference (H)
            # into Bob's short arm, Alice's encoded (V) into his long
            # arm -- the L-S and S-L paths.  Only these two exist, so
            # there are no satellite bins.
            h, v = pbs(E_field)
            arm_ref = np.column_stack([h, np.zeros_like(h)])
            arm_enc = np.column_stack([v, np.zeros_like(v)])
            # Equalise with a real attenuator before interfering.
            arm_ref = voa(arm_ref, balance_dB)
            # Bob's stretcher sits in his LONG arm, which is `arm_enc`, and
            # it is a real part with a real insertion loss.  Only the loss
            # goes here: it does not depend on drive voltage, so it is a
            # constant on one arm and needs no feedback.  The phase the
            # servo commands is voltage-dependent and per-block, and enters
            # `delta` further down with the other one-arm phases.
            #
            # Being on ONE arm, the loss unbalances the pair and so touches
            # visibility as well as total power -- which is why it belongs
            # in the field chain rather than as a scalar on the result.
            if _stretcher is not None and _stretcher.apply_insertion_loss:
                arm_enc = voa(arm_enc, _stretcher.insertion_loss_db)
            # Bob's modulator sits in his short arm, which carries
            # Alice's reference.  phi_A rides the encoded path and phi_B
            # the reference path, so the interfering pair differ by
            # phi_A - phi_B -- the convention PHI_A and PHI_B assume.
            #
            # The sign here is coupled to `enc`'s phase_arm: with phi_A
            # on the encoded arm the reference must carry +phi_B.
            # Flipping one without the other gives phi_A + phi_B and
            # inverts the Y basis while leaving X untouched -- verified
            # as the negative control for the gate-table equivalence
            # check.
            arm_ref = arm_ref * np.exp(1j * phi_b)
            E_c, E_d = dec.modulate((arm_ref, arm_enc), dt)
        else:
            E_c, E_d = dec.modulate(E_field, dt, phase=phi_b)

        P_c = float(np.mean(np.sum(np.abs(E_c[start_i:end_i]) ** 2, axis=1)))
        P_d = float(np.mean(np.sum(np.abs(E_d[start_i:end_i]) ** 2, axis=1)))
        return P_c, P_d

    # --- Interference coefficients (see module docstring) ---------------
    #
    # The 8 outcomes reachable from the discrete choices were precomputed
    # (a_basis, a_bit, b_basis) choices.  That is exact only while the
    # response is deterministic; a per-pulse random phase -- laser
    # linewidth against a residual path mismatch, or modulator noise --
    # makes it continuous and puts 8 answers out of reach.
    #
    # The same linearity argument extends.  The chain is linear, so the
    # two interfering paths superpose and the gated power is exactly
    #
    #     P(delta) = g0 + 2*Re(S * exp(i*delta)),   delta = phi_A - phi_B,
    #
    # for constants (g0, S) fixed by the point's parameters.  Three
    # evaluations at delta = 0, pi/2, pi determine them:
    #
    #     g0     = (P(0) + P(pi)) / 2
    #     Re(S)  = (P(0) - P(pi)) / 4
    #     Im(S)  = (g0 - P(pi/2)) / 2
    #
    # Extracting the coefficients from the field chain itself, rather than
    # re-deriving them by hand, keeps this exact by construction for both
    # topologies -- there is no second expression of the physics to drift
    # out of step.  Any delta then costs two multiplies.
    def _coeffs(index, op_fibre):
        P0 = _gate_powers(0.0, 0.0, op_fibre)[index]
        Ph = _gate_powers(np.pi / 2.0, 0.0, op_fibre)[index]
        Pp = _gate_powers(np.pi, 0.0, op_fibre)[index]
        g0 = 0.5 * (P0 + Pp)
        return g0, complex(0.25 * (P0 - Pp), 0.5 * (g0 - Ph))

    def _extract(op_fibre):
        """(g0_c, S_c, g0_d, S_d) for one fibre state, chain re-verified."""
        g0_c, S_c = _coeffs(0, op_fibre)
        g0_d, S_d = _coeffs(1, op_fibre)

        # The coefficients are extracted by sweeping phi_A at phi_B = 0, so
        # they carry no information about how phi_B enters.  The per-pulse
        # formula then *assumes* delta = phi_A - phi_B.  If Bob's sign in
        # `_gate_powers` were ever flipped, that assumption would silently
        # become wrong and the closed form would keep returning the old
        # answers -- the field chain would no longer be authoritative.
        #
        # So assert it: evaluate the chain at points with phi_B != 0 and
        # require the closed form to reproduce them.  Three extra field
        # evaluations per extraction, and it fails loudly on exactly the
        # phase-arm/sign coupling that had to be verified by hand.
        _scale = max(abs(g0_c) + 2.0 * abs(S_c), 1e-300)
        for _pa, _pb in ((0.0, np.pi / 2.0), (np.pi, np.pi / 2.0),
                         (np.pi / 2.0, np.pi / 2.0)):
            _ref_c, _ref_d = _gate_powers(_pa, _pb, op_fibre)
            _e = np.exp(1j * (_pa - _pb))
            _got_c = g0_c + 2.0 * (S_c * _e).real
            _got_d = g0_d + 2.0 * (S_d * _e).real
            if (abs(_got_c - _ref_c) > 1e-9 * _scale
                    or abs(_got_d - _ref_d) > 1e-9 * _scale):
                raise RuntimeError(
                    "interference closed form disagrees with the field chain "
                    f"at (phi_A={_pa:.4f}, phi_B={_pb:.4f}): the assumed "
                    "relative phase delta = phi_A - phi_B does not match the "
                    "chain. This usually means the encoder's phase_arm and "
                    "Bob's phase sign have been changed independently -- they "
                    "are coupled. See the phase-arm convention."
                )
        return g0_c, S_c, g0_d, S_d

    # Per-pulse random phase.  Both sources are zero-mean Gaussians and
    # add in the same delta, so one draw covers them:
    #   * laser linewidth against the RESIDUAL path mismatch left after
    #     the delay line and stretcher (the matched S-L / L-S routes
    #     cancel the full-delay term -- see the module docstring);
    #   * phase-modulator drive noise.
    theta_sigma = math.sqrt(2.0 * math.pi * linewidth * path_mismatch) \
        if (linewidth > 0.0 and path_mismatch > 0.0) else 0.0
    theta_sigma = math.hypot(theta_sigma, phase_noise_rad)

    # Alice's modulator.  Built even when every knob is zero so that the
    # bias law -- volts to radians through the crystal-derived V_pi -- and
    # the "not both units at once" check live in the component and are not
    # restated here.  Pushing pulses through `modulate` would cost a field
    # propagation each; reading the resolved offset costs nothing.
    pm_alice = PhaseModulator(crystal_cut='X', modulation='DC',
                              phase_error_rad=phase_error_rad,
                              bias_offset_v=bias_offset_v)
    static_phase_error = pm_alice.phase_error_rad

    # Bob's arm-length drift, likewise owned by the interferometer.  The
    # coefficients above were extracted with the chain at t = 0, so they
    # already contain `arm_phase_offset(0)`; only the increment since then
    # belongs in `delta`.  Subtracting the component's own t=0 value keeps
    # the ramp itself out of this file.
    arm_offset_0 = dec.arm_phase_offset(0.0)
    drifting = dec.phase_drift_rad_s != 0.0

    # TWO CLOCKS, and they are not the same clock.
    #
    #   detector clock  t_pulse = pulse_idx / repetition_rate
    #       Real elapsed time between gates.  Dead time and afterpulsing
    #       are defined against it and must keep true 1/f spacing.
    #
    #   drift clock     t_drift = (pulse_idx / num_bits) * run_duration
    #       Position within the *experiment's* drift profile.
    #
    # They coincide only when num_bits/repetition_rate == run_duration.
    # Keeping one clock for both is wrong, and wrong in a way that hides:
    # the pulse budget is chosen for statistical power, so tying drift to
    # it means a longer run silently becomes a longer *experiment*.  In the
    # Gobby sweep the 122 km point needs 1e9 pulses = 500 s at 2 MHz, which
    # accumulated 25 deg of drift against the 6 deg of their stated
    # two-minute transfer and inflated the effective modulation error from
    # 3.31% to 8.60% -- the QBER read 13.52% against their 8.9%, and the
    # whole excess was this.
    #
    # More pulses must mean a better estimate of the same experiment, not a
    # longer one.  With `run_duration` set, the simulated pulses are a
    # uniform sample across a run of that length; the count may exceed what
    # the apparatus actually sent in that time, which is exactly what makes
    # it a lower-variance estimate of the same expectation.
    #
    # Default None keeps the old behaviour (drift tied to the pulse budget)
    # so nothing changes for callers that do not declare a duration.
    if run_duration is not None and num_bits > 1:
        _drift_scale = float(run_duration) / (num_bits - 1)
    else:
        _drift_scale = 1.0 / repetition_rate

    # --- Fibre drift, and why it needs blocks -----------------------------
    #
    # The closed form fixes (g0, S) for ONE fibre state, so a fibre that
    # changes during the run puts it out of reach the same way a per-pulse
    # random phase does.  Re-extracting per pulse would cost the field
    # propagations the closed form exists to avoid, so the run is cut into
    # blocks: within a block the fibre is held still and the closed form is
    # exact, and between blocks it is re-extracted from `fibre.at(t)`.
    #
    # `drift_blocks` is a COUNT, not a pulse size, and that is deliberate.
    # A pulse size would tie the drift resolution to the statistical budget,
    # which is the mistake `run_duration` exists to prevent one level up --
    # asking for tighter error bars must not silently change the physics.
    # As a count, doubling num_bits leaves the resolution alone.
    #
    # Blocks advance on the DRIFT clock, never the detector clock.  See the
    # two-clocks note above; this must not introduce a third.
    #
    # Cost is three field propagations plus one Jones build (about 2 ms at
    # 122 km) per block, so the default 100 blocks is negligible against the
    # 1e8-pulse runs the long distances need.
    #
    # Not drifting means exactly one block and nothing rebuilt, so the
    # no-drift path is provably unchanged rather than approximately so.
    _fibre_drifts = (fibre is not None
                     and fibre.drift_temperature_rate_C_s != 0.0)
    _n_blocks = max(1, int(drift_blocks)) if _fibre_drifts else 1

    def _block_bounds(i):
        return (i * num_bits) // _n_blocks, ((i + 1) * num_bits) // _n_blocks

    def _extract_for(lo, hi):
        """Coefficients for the block [lo, hi), fibre taken at its midpoint."""
        t_mid = 0.5 * (lo + max(hi - 1, lo)) * _drift_scale
        return _extract(fibre if fibre is None else fibre.at(t_mid)), t_mid

    # --- Bob's phase servo ------------------------------------------------
    #
    # Gobby's Bob holds his operating point with the piezo-driven fibre
    # stretcher in his long arm.  Without it the model contradicts the
    # paper: a residual rotation is SU(2), so it costs O(eps^2) in
    # amplitude against O(eps) in phase, and the QBER therefore moves
    # before the rate does -- the reverse of "polarisation drift reduces
    # the bit rate, but does not degrade the QBER".
    #
    # The error signal costs nothing, because it is already being
    # computed.  The closed form is P = g0 + 2*Re(S*exp(i*delta)), so a
    # fibre phase theta on the interfering pair sends S -> S*exp(i*theta)
    # and shifts arg(S) by exactly theta.  `arg(S_c)` IS the fringe phase,
    # which is also the quantity a real lock-in servo measures, so reading
    # it here is the physical thing rather than a shortcut around it.
    #
    # Sample and hold, because a real servo has finite bandwidth: the lock
    # is refreshed when it is older than the interval, and between
    # refreshes the fibre keeps moving while the correction does not.  A
    # fast interval reproduces the paper; a slow one degrades back to the
    # unserved behaviour, and the crossover is measurable.
    #
    # PHASE ONLY, and that is the whole point.  Inverting the full Jones
    # matrix would null the amplitude too and the rate loss the paper
    # reports would vanish with it.
    _servo_on = (phase_servo_interval_s is not None
                 and phase_servo_interval_s >= 0.0 and _fibre_drifts)

    # The device behind the correction.  Default is the Thorlabs FVP155P,
    # whose insertion loss is SUPPRESSED here and nowhere else: ETA_BOB in
    # the Gobby replication already folds in "5 dB of loss in Bob's
    # apparatus", and the stretcher sits inside that apparatus, so applying
    # its 0.1 dB again would double-count against a number from the paper.
    # A chain whose budget does not already contain Bob's optics should
    # pass a stretcher with the loss left on.
    #
    # Its resonance, 80 kHz, is the floor on a useful re-lock interval.
    # Ours are of order a second, five orders of magnitude slower, so the
    # device's bandwidth never binds here.
    _stretcher = phase_servo_stretcher
    if _servo_on and _stretcher is None:
        _stretcher = PiezoFibreStretcher(apply_insertion_loss=False)
    _servo_held = 0.0        # the correction Bob is currently holding
    _servo_locked_at = None  # drift-clock time of the last lock

    def _servo_update(S, t_mid):
        """Re-lock if the held correction is older than the interval."""
        nonlocal _servo_held, _servo_locked_at
        if (_servo_locked_at is None
                or t_mid - _servo_locked_at >= phase_servo_interval_s):
            theta = cmath.phase(S) - _servo_ref
            # Through the device and back: what Bob can actually deliver,
            # not what he would like to.  `voltage_for` wraps into one
            # fringe first, which is why a 7*pi stroke is enough for a
            # correction that would otherwise grow without bound.
            _servo_held = _stretcher.phase_for(_stretcher.voltage_for(theta))
            _servo_locked_at = t_mid

    _blk = 0
    _blk_lo, _blk_end = _block_bounds(0)
    (g0_c, S_c, g0_d, S_d), _t_mid = _extract_for(_blk_lo, _blk_end)

    # Bob's calibration point: the fringe phase he aligned against, which
    # every later reading is measured relative to.
    #
    # Taken at t = 0 explicitly, not from the first block.  `U_comp` is
    # built from the fibre at t = 0, and the two references have to be the
    # same instant or the servo starts out holding a correction for drift
    # that alignment already removed.  The first block's midpoint is close
    # to zero at the default 100 blocks and not close at all at one, which
    # is exactly the sort of resolution-dependent physics `drift_blocks`
    # exists to avoid.
    #
    # Costs one extra extraction per run, and only when the servo is on.
    _servo_ref = cmath.phase(S_c)
    if _servo_on:
        _servo_ref = cmath.phase(_extract(fibre)[1])
        _servo_update(S_c, _t_mid)

    for pulse_idx in range(num_bits):
        if pulse_idx >= _blk_end:
            _blk += 1
            _blk_lo, _blk_end = _block_bounds(_blk)
            (g0_c, S_c, g0_d, S_d), _t_mid = _extract_for(_blk_lo, _blk_end)
            if _servo_on:
                _servo_update(S_c, _t_mid)

        # --- Alice's encoding ---
        alice_basis = random.choice(['X', 'Y'])
        alice_bit = random.randint(0, 1)

        # --- Bob's decoding basis ---
        bob_basis = random.choice(['X', 'Y'])

        # --- Interference power from the precomputed coefficients ---
        t_pulse = float(pulse_idx) / repetition_rate
        delta = (PHI_A[alice_basis][alice_bit] - PHI_B[bob_basis]
                 + static_phase_error)
        if drifting:
            # Drift clock, NOT the detector clock -- see the note above.
            delta += dec.arm_phase_offset(pulse_idx * _drift_scale) \
                - arm_offset_0
        if _servo_on:
            # Bob's stretcher, subtracting the correction he is currently
            # holding.  Same place and same convention as the arm-length
            # offset above: a phase applied to one arm enters `delta`.
            delta -= _servo_held
        if theta_sigma > 0.0:
            delta += random.gauss(0.0, theta_sigma)
        e = complex(math.cos(delta), math.sin(delta))
        P_c = g0_c + 2.0 * (S_c.real * e.real - S_c.imag * e.imag)
        P_d = g0_d + 2.0 * (S_d.real * e.real - S_d.imag * e.imag)

        # --- Detection ---
        click_c = spad_c.detect(P_c, t_pulse)
        click_d = spad_d.detect(P_d, t_pulse)

        # --- Determine Bob's bit ---
        if click_c and not click_d:
            bob_bit = 0
        elif click_d and not click_c:
            bob_bit = 1
        elif click_c and click_d:
            bob_bit = random.randint(0, 1)
        else:
            bob_bit = -1  # no click — discard

        # --- Sifting: same basis AND a click ---
        if click_c or click_d:
            n_clicks += 1
            if alice_basis == bob_basis:
                n_sifted += 1
                if alice_bit != bob_bit:
                    n_errors += 1

        if verbose and (pulse_idx + 1) % max(1, num_bits // 10) == 0:
            print(f"  Pulse {pulse_idx+1}/{num_bits}", flush=True)

    qber = n_errors / n_sifted if n_sifted > 0 else 0.0

    fiber_loss_dB = alpha_dB * fiber_length
    sifted_key_rate = n_sifted / (num_bits / repetition_rate) if num_bits > 0 else 0.0

    results = {
        'qber': qber,
        'n_total': num_bits,
        'n_sifted': n_sifted,
        'n_errors': n_errors,
        'sifted_key_rate': sifted_key_rate,
        'fiber_length_km': fiber_length,
        'total_loss_dB': fiber_loss_dB,
        'mu': mu,
    }

    if verbose:
        print(f"\nTime-bin BB84 — {fiber_length} km")
        print(f"  Fibre loss: {fiber_loss_dB:.1f} dB")
        print(f"  Mu: {mu} photons/pulse")
        print(f"  Total pulses: {num_bits}")
        print(f"  Clicks: {n_clicks} ({n_clicks/num_bits*100:.2f}%)")
        print(f"  Sifted: {n_sifted} ({n_sifted/num_bits*100:.2f}%)")
        print(f"  Errors: {n_errors}")
        print(f"  QBER: {qber*100:.2f}%")

    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Time-bin phase-encoding BB84 QKD (Gobby et al. 2004)")
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--bits', type=int, default=10000,
                        help='Number of pulses (default 10k)')
    parser.add_argument('--fiber-length', type=float, default=0,
                        help='Fiber length in km (default 0)')
    parser.add_argument('--mu', type=float, default=0.1,
                        help='Mean photons per pulse (default 0.1)')
    parser.add_argument('--visibility', type=float, default=1.0,
                        help='Decoder interferometric visibility, (0,1] (default 1.0)')
    parser.add_argument('--phase-error', type=float, default=0.0,
                        help='Static decoder phase offset in rad (default 0.0)')
    parser.add_argument('--verbose', action='store_true')
    args = parser.parse_args()

    results = simulate_bb84_time_bin(
        num_bits=args.bits, fiber_length=args.fiber_length,
        mu=args.mu, seed=args.seed, verbose=True,
        visibility=args.visibility, phase_error=args.phase_error)
