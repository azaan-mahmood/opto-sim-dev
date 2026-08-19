"""Replication of Duplinskiy et al., Opt. Express 25(23), 28886, 2017.

"Low loss QKD optical scheme for fast polarization encoding"

Proof-of-concept BB84 QKD over 50 km fiber spool (10 dB loss)
using LiNbO3 phase modulators and InGaAs SPADs.

Signal chain (mirrors bb84_ideal.py but with VOA + SPAD):
  Alice: polarizer('45') → PM1 (encode)
  Channel: propagate (birefringence + attenuation)
  Bob: PM2 (basis select) → circular_analyser → 2x SPAD

Detection uses the same 50:50 BS with π/2 phase shift
(optics.circular_analyser, not a true PBS):
  P_x − P_y = sin(Δφ_alice + Δφ_bob)
  Same basis: |sin| = 1 → deterministic bit
  Diff basis: sin ≈ 0 → random (sifted out)

Encoding (same as bb84_ideal.py):
  X basis (linear):  bit 0 → V=0,       bit 1 → V=Vpi
  C basis (circular): bit 0 → V=Vpi/2,  bit 1 → V=3·Vpi/2
  Bob X basis: V_bob = Vpi/2
  Bob C basis: V_bob = 0

Validation note
---------------
The 0 km (back-to-back) QBER validates the SPAD + phase modulator +
circular-analyser detection chain against the paper's baseline
measurement.  At 0 km there
is no fibre, so birefringence does not enter — only detector dark counts,
afterpulsing, and the intrinsic sifting loss contribute.

Distance-dependent QBER values are comparable to the paper *with
polarization compensation enabled* (the default): the channel is
quasi-static, so Bob applies the inverse of the fibre's single SU(2)
Jones matrix (obtained from `FiberRealization.birefringence_matrix()`).
This mirrors the paper's calibration procedure, which maintains
alignment with polarization controllers and automatic recalibration —
for a quasi-static channel, that loop converges to exactly this
inverse.  Pass ``--no-compensation`` to see the uncompensated channel:
because the fibre applies a fixed random rotation, the received SOP no
longer matches the encoding axes and QBER collapses toward ~50 %
(25.9 % at 50 km in our seed-42 run), which is what a polarization-
encoded system *without* any active stabilization sees.

References
----------
[1] Duplinskiy et al., Opt. Express 25(23), 28886-28897, 2017.
[2] ID Quantique, ID230 InGaAs SPAD datasheet.
"""
import argparse
import warnings

import numpy as np
import random
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from src.channel import optics, FiberRealization
from src.channel.phase_modulator import PhaseModulator
from src.detectors.spad import spad


class _Same:
    """Sentinel: 'this calibration setting matches the operating one'.

    `None` cannot carry that meaning for `calibration_bend_radius`,
    because `None` is itself a valid operating bend radius -- it means
    "unbent", which is a perfectly reasonable state to have calibrated in
    before the fibre was disturbed.  Overloading the two made
    `bend_radius=0.1, calibration_bend_radius=None` silently mean "no
    mismatch" when it reads like "calibrated straight, now bent".
    """

    def __repr__(self):
        return 'SAME_AS_OPERATING'


SAME_AS_OPERATING = _Same()


def simulate_bb84_duplinskiy(num_bits, fiber_length=50, alpha_dB=0.2,
                              mu=0.1, bob_loss_dB=2.0,
                              gate_width=20e-9, rep_rate=10e6,
                              compensate=True, model='auto',
                              bias_offset_v=0.0,
                              spad_eta=0.10, dead_time=13e-6,
                              dark_count_rate=15.0, afterpulse_prob=0.05,
                              cd=False, pmd=False,
                              extinction_epsilon=0.0,
                              temperature=25.0, bend_radius=None,
                              calibration_temperature=SAME_AS_OPERATING,
                              calibration_bend_radius=SAME_AS_OPERATING,
                              drift_temperature_rate_C_s=0.0, drift_blocks=100,
                              run_duration=None,
                              source_field=None, source_dt=None,
                              pulse_energy_factors=None,
                              block_size=None,
                              seed=None, verbose=False):
    """BB84 simulation matching the Duplinskiy et al. experimental setup.

    Parameters
    ----------
    num_bits : int — number of BB84 pulses to simulate.
    fiber_length : float — fiber length in km (default 50).
    alpha_dB : float — fiber attenuation in dB/km (default 0.2, SMF-28).
    mu : float — mean photon number per pulse (default 0.1).
    bob_loss_dB : float — Bob's internal loss in dB (default 2.0).
    gate_width : float — SPAD gate window in seconds (default 20 ns).
    rep_rate : float — laser repetition rate in Hz (default 10 MHz).
    compensate : bool — apply Bob's polarization compensation (the
        inverse of the fibre's quasi-static Jones matrix) before
        decoding, mirroring the paper's calibration loop (default True).
    spad_eta : float — SPAD quantum efficiency (default 0.10, ID230 [2]).
    dead_time : float — SPAD dead time in seconds (default 13e-6, ID230).
    dark_count_rate : float — SPAD dark count rate in Hz (default 15.0,
        ID230).  Unlike the Gobby replication, ID230 is *contemporary*
        with this 2017 experiment, so the datasheet figures are era-
        appropriate rather than anachronistic.
    afterpulse_prob : float — afterpulse probability per click (default
        0.05, ID230).  In the Gobby chain this value was found to be
        standing in for the modulation error and belonged at 0 there; that
        argument is specific to Gobby's stated error budget and does not
        transfer.  Exposed here so the question can be asked rather than
        assumed.
    cd : bool — apply chromatic dispersion (default False).
    pmd : bool — apply polarisation-mode dispersion (default False).
        Both were hardcoded off until this chain gained them.  It is the one
        worth sweeping them on: a polarisation-encoded observable responds
        to them, where the time-bin chain is invariant by construction.
    bias_offset_v : float — static bias error on Alice's phase modulator,
        as a drive voltage (default 0 = perfectly biased).  Converted by
        `PhaseModulator` through its crystal-derived V_pi.  Setting a
        modulator's bias imperfectly is universal to phase-modulated QKD,
        not specific to any one experiment.
    extinction_epsilon : float — finite analyser extinction, as the
        fraction of each port's power appearing at the other (default 0 =
        a perfect analyser, which is what this chain assumed before
        this chain assumed a perfect one).  Applied at the analyser
        output *before* detection, so it
        cannot mix dark counts and afterpulses into each other; those are
        generated inside the SPAD and are already uncorrelated with
        Alice's bit.

        The paper's calibration goal 3 is "Bob's measurements
        differentiate BB84 orthogonal states with extinction higher than
        98 %" (§5).  Two readings of that differ by exactly a factor of
        two — power-fraction gives 0.0200, visibility-like gives 0.0101 —
        and the ambiguity lands on the quantity the afterpulse test turns
        on, so both are run rather than one chosen (register A7).

        Note 98 % is a *threshold the tuning algorithm targets*, not an
        achieved operating value, so anything derived from it is an upper
        bound on the term (register A3).

        See `optics.apply_extinction` for why this is a lumped term
        standing for three mechanisms, and why calling it "the PBS
        extinction" would be wrong.
    temperature : float — ambient fibre temperature in C (default 25,
        the value formerly hardcoded here).  Shifts the sectional
        birefringence by -3e-9 per degree.
    bend_radius : float or None — fibre bend radius in metres (default
        None = unbent).  Adds the Ulrich bend birefringence
        `0.135*(r_fiber/R)^2`.

        **Both are exactly inert while Bob compensates the same fibre**,
        which is the default.  `U_comp` is the conjugate transpose of the
        channel's own Jones matrix, so `U_comp @ J = I` for any unitary
        `J`, and these two parameters only change `J`.  That is
        arithmetic, not a physical null, and must not be reported as one.
        They matter through `calibration_*` below, or with
        `compensate=False`.
    calibration_temperature : float — the temperature Bob's compensator
        was calibrated at.  Defaults to `SAME_AS_OPERATING`, i.e. perfect
        calibration, and the behaviour when the compensator is perfect.
    calibration_bend_radius : float or None — likewise for bend radius.
        **`None` here means "calibrated unbent", not "same as
        operating"** — use the `SAME_AS_OPERATING` sentinel for that,
        which is the default.  The two are different physical statements
        and `None` cannot carry both: an unbent calibration followed by a
        disturbed fibre is exactly the case worth simulating.

        Setting either makes Bob invert the fibre as it *was* while light
        travels through it as it *is*.  This is the paper's stated reason
        for spending 20 % of its duty cycle recalibrating: the channel is
        "not isolated from external influences, including mechanical and
        temperature ones" (§6).  Measured residual
        `||U_comp(25 C) @ J(25+dT) - I||_F` at 50 km: 0.0017 at dT = 1e-4,
        0.173 at 0.01, 1.771 at 0.1, saturating near 2.19 by 1 C.  The
        calibration fibre shares the operating fibre's seed, so the two
        differ only in the environmental term rather than being unrelated
        draws.
    drift_temperature_rate_C_s : float — rate (C/s) at which the fibre's
        ambient temperature changes DURING the run, default 0.0 (static).

        The `calibration_*` pair above gives a fixed two-state mismatch;
        this is the time dependence they cannot express, where the residual
        `U_comp @ J(t)` grows from the identity as the run proceeds. That
        growth is what forces recalibration, so it is the mechanism behind
        the paper's 20 % duty cycle rather than a proxy for it.
    drift_blocks : int — how many times the response table is rebuilt
        across the run while the fibre drifts (default 100). A count, not
        a pulse size, so raising `num_bits` for statistical power leaves
        the drift resolution alone. Independent of `block_size`, which
        slices the run for reporting rather than for physics. Ignored when
        nothing drifts.
    run_duration : float or None — physical duration of the experiment
        being simulated, in seconds (default None), used as the DRIFT clock
        only. The detector clock always keeps true 1/rep_rate spacing,
        because dead time and afterpulsing are defined against real elapsed
        time.

        None ties drift to the pulse budget, which is the historical
        behaviour and bit-identical to it. Prefer setting it: the budget is
        chosen for statistical power, so leaving drift tied to it means
        asking for tighter error bars silently lengthens the simulated
        experiment. Same argument, and the same fix, as in `bb84_time_bin`.
    source_field : ndarray (n, 2) complex, or None — the optical field
        entering Alice's phase modulator.  Default None keeps the flat
        analytic field this chain has always used,
        `sqrt(power_per_pulse/2) * ones((1, 2))`, which is a single time
        sample of diagonal (D) polarisation.  Pass a real source — e.g.
        `LaserDriver.sample_field(...)` from the DFB device model — to
        drive the chain from a laser instead.  It is renormalised to
        `power_per_pulse` on the way in, so mu is calibrated once at the
        source, matching the project convention.

        The field must arrive already polarised at 45 degrees to the
        modulator axes.  That is the paper's PC1: "the polarization
        controller (PC 1) is configured in such a way that the amplitudes
        of the field along the ordinary and extraordinary axes of the
        crystal inside the modulator (PM 1) are equal" (§2).  For the DFB
        that is `polarization_azimuth=pi/4` on the driver.

        Measured: swapping the flat field for the DFB changes the
        8-outcome response table by at most 2.6e-15 relative, and leaves
        the Stokes vectors of all four BB84 states exact with DOP = 1.
        That is structural rather than lucky -- `sample_field` returns one
        complex amplitude times a fixed Jones vector, so both components
        carry the *same* amplitude and the normalised Stokes parameters,
        which depend only on the ratio Ey/Ex, cannot see the source's RIN,
        phase noise or chirp at all.  A polarisation-encoding chain is
        therefore blind to everything a real source adds except pulse
        energy.  Same cancellation as the linewidth argument in the
        time-bin chain.
    source_dt : float or None — the sampling interval of `source_field`, in
        seconds.  Only CD and PMD need it: both are frequency-domain
        operators, and `FiberRealization.apply` builds its frequency grid
        from `np.fft.fftfreq(N, d=dt)`.

        Without this the chain passed the PULSE PERIOD (1/rep_rate, 100 ns
        at 10 MHz) as `dt`, which is the spacing between pulses rather than
        between samples within one.  On a 0.5 ps source grid that
        understates the bandwidth by five orders, and the PMD phase
        `omega*dgd/2` comes out around 1e-17 rad -- an exact-looking null
        produced by the wrong time base rather than physics.

        Required when `cd` or `pmd` is on and `source_field` carries more
        than one sample.  A single-sample field warns instead: CD and PMD
        cannot act on one time sample at all, so any null measured there is
        arithmetic.
    pulse_energy_factors : sequence of float, or None — per-pulse energy
        multipliers with mean 1, the one thing a real source does change
        here.  Applied to the response power before detection, which is
        exact because the whole chain is linear in the field.

        Consumed cyclically in index order, never sampled randomly.  Two
        reasons: a draw would consume RNG and break bit-identity with the
        frozen baseline even when every factor is 1.0, and the measured
        pulses may be correlated through carrier recovery, which cycling
        preserves and sampling would destroy.
    block_size : int or None — if set, also accumulate (n_sifted, n_errors)
        per block of this many pulses and return them under `blocks`.
        Default None keeps the single aggregate this chain has always
        returned, and consumes no RNG either way, so a frozen baseline
        is bit-identical with the parameter present.

        This exists for the paper's Fig. 7, a QBER-versus-time trace.  One
        call with blocking is ONE CONTINUOUS LINK: SPAD dead time and
        pending afterpulses carry across block boundaries, which is what a
        block-to-block variance measurement is actually about.  Running N
        separate calls instead would reset detector state at every
        boundary, and -- because `seed` also seeds `FiberRealization` --
        silently draw a different fibre per block.
    seed : int or None — RNG seed.
    model : str — birefringence model for the fibre: 'auto' (default) or
        'sectional'. Any other value raises ValueError.
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
    wavelength = 1550e-9
    photon_energy = h * c / wavelength

    # Calibrate field power so that mean(|E|^2) * gate_width gives mu photons
    # power_per_pulse = mu * h * nu / gate_width
    power_per_pulse = mu * photon_energy / gate_width

    # Phase modulators (LiNbO3, X-cut)
    pm_alice = PhaseModulator(crystal_cut='X', modulation='DC',
                              bias_offset_v=bias_offset_v)
    pm_bob = PhaseModulator(crystal_cut='X', modulation='DC')
    Vpi = pm_alice.Vpi

    # SPAD detectors.  These were four hardcoded literals;
    # they are parameters now so the detector budget can be interrogated
    # the way Gobby's was.  Defaults are the ID230 datasheet figures [2]
    # and reproduce the previous behaviour exactly.
    spd1 = spad(wavelength=wavelength, quantum_efficiency=spad_eta,
                dead_time=dead_time, dark_count_rate=dark_count_rate,
                afterpulse_prob=afterpulse_prob, gate_width=gate_width)
    spd2 = spad(wavelength=wavelength, quantum_efficiency=spad_eta,
                dead_time=dead_time, dark_count_rate=dark_count_rate,
                afterpulse_prob=afterpulse_prob, gate_width=gate_width)

    dt_pulse = 1.0 / rep_rate

    # One physical fibre for the whole run: birefringence is quasi-static,
    # so every pulse must see the same Jones matrix. Built once here,
    # reused per pulse
    # below.  `cd` and `pmd` default off, matching this script's original
    # explicit propagate(cd=False, pmd=False) call; they are exposed as of
    # Exposed so this chain can carry an impairment sweep.  It is the right
    # chain for one: polarisation encoding is *not* blind to birefringence
    # the way the time-bin chain is.
    # `temperature` and `bend_radius` were hardcoded at 25 C and None until
    # They reach the sectional birefringence model, which shifts
    # delta_n by -3e-9 per degree and adds the Ulrich bend term
    # 0.135*(r/R)^2 (see `src/channel/fiber.py`).
    fibre = FiberRealization(
        L_m=fiber_length * 1000, temperature=temperature,
        bend_radius=bend_radius, attenuation_factor=alpha_dB,
        cd=cd, pmd=pmd, model=model, seed=seed,
        drift_temperature_rate_C_s=drift_temperature_rate_C_s)

    # Bob's polarization compensation (default on): the channel's Jones
    # matrix is unitary, so the inverse is its conjugate transpose.
    # Applied per pulse before decoding, mirroring the paper's calibration
    # loop (see module docstring).
    #
    # By default Bob compensates the fibre light actually travels through,
    # which makes `temperature` and `bend_radius` EXACTLY inert here --
    # provably, not approximately: U_comp = J^dagger gives U_comp @ J = I
    # for any unitary J, and those two parameters only change J.  A null
    # measured on this path is therefore arithmetic, not physics, and must
    # not be reported as a finding.
    #
    # They become live under a calibration mismatch: Bob's compensator is
    # built from the fibre in one state while light propagates through it
    # in another.  That is the paper's own failure mode -- "the other 20%
    # has been required for recalibrations as the quantum channel has not
    # been isolated from external influences, including mechanical and
    # temperature ones" (paper §6) -- and it is what `calibration_temperature`
    # and `calibration_bend_radius` express.
    #
    # The calibration fibre uses the SAME seed, so the two realizations
    # share their random section axes and differ only in the environmental
    # term.  Drawing a fresh seed instead would compare two unrelated
    # fibres, which is a different and much less interesting experiment.
    cal_T = (temperature if calibration_temperature is SAME_AS_OPERATING
             else calibration_temperature)
    cal_R = (bend_radius if calibration_bend_radius is SAME_AS_OPERATING
             else calibration_bend_radius)

    if cal_T == temperature and cal_R == bend_radius:
        cal_fibre = fibre
    else:
        cal_fibre = FiberRealization(
            L_m=fiber_length * 1000, temperature=cal_T, bend_radius=cal_R,
            attenuation_factor=alpha_dB, cd=cd, pmd=pmd, model=model, seed=seed)

    J_channel = cal_fibre.birefringence_matrix()
    U_comp = None if J_channel is None else J_channel.conj().T

    # --- The source field -------------------------------------------------
    #
    # Default: the flat analytic field this chain has always used -- one
    # time sample, equal amplitude on x and y, i.e. diagonal (D)
    # polarisation, which is what the paper's PC1 delivers to PM1.
    #
    # With `source_field` given, a real laser drives the chain instead.  It
    # is renormalised to the same power_per_pulse, so mu stays calibrated
    # once at the source and the swap changes the *statistics* of the
    # light, not its mean level.  `bb84_time_bin` renormalises the same way.  Anything
    # passed here must already be polarised at 45 degrees to the modulator
    # axes; see the parameter docstring.
    if source_field is None:
        E_source = np.sqrt(power_per_pulse / 2.0) * np.ones((1, 2), dtype=complex)
    else:
        E_source = np.asarray(source_field, dtype=complex)
        if E_source.ndim != 2 or E_source.shape[1] != 2:
            raise ValueError(f"source_field must have shape (n, 2), got {E_source.shape}")
        p_mean = float(np.mean(np.sum(np.abs(E_source) ** 2, axis=1)))
        if p_mean <= 0:
            raise ValueError("source_field carries no power")
        E_source = E_source * np.sqrt(power_per_pulse / p_mean)

    # The time base CD and PMD are evaluated on.  See the `source_dt`
    # docstring: handing them the pulse period instead of the sample
    # interval produces a null that looks exact and means nothing.
    n_samples = E_source.shape[0]
    if cd or pmd:
        if n_samples > 1 and source_dt is None:
            raise ValueError(
                "cd or pmd is enabled and source_field carries "
                f"{n_samples} samples, but source_dt was not given. Falling "
                "back to the pulse period would evaluate both operators on "
                "the wrong time base and return a false null.")
        if n_samples == 1:
            warnings.warn(
                "cd or pmd is enabled but the source field is a single time "
                "sample. Both are frequency-domain operators, so neither can "
                "act: the result is a null by construction, not a "
                "measurement. Pass a time-resolved source_field with "
                "its source_dt.", UserWarning, stacklevel=2)
    dt_field = dt_pulse if source_dt is None else source_dt

    factors = None if pulse_energy_factors is None else np.asarray(
        pulse_energy_factors, dtype=float)
    if factors is not None and factors.size == 0:
        raise ValueError("pulse_energy_factors is empty")

    # --- Precomputed response table --------------------------------------
    #
    # The field chain is DETERMINISTIC given (alice_basis, alice_bit,
    # bob_basis) -- 8 combinations -- because the fibre Jones matrix is
    # sampled once for the whole run (quasi-static) and no stage
    # between the source and the detectors consumes randomness.  Verified:
    # `pm.modulate`, `fibre.apply`, `optics.voa` and `circular_analyser`
    # all leave both RNG streams untouched, and `fibre.apply` is repeatable
    # on identical input.
    #
    # So the whole per-pulse chain -- modulate, propagate, compensate, VOA,
    # modulate, analyse -- collapses to a table lookup.  This is exactly
    # the same argument as `bb84_time_bin.py`, and it is what makes a
    # target-sifted polarisation sweep affordable: the chain measured
    # ~30,000 pulses/s walking the fields per pulse.
    #
    # The table is built by running the REAL chain, never by re-deriving
    # the physics by hand, so there is no second expression of it to drift
    # out of step.
    #
    # NOTE: this is exact only while the response stays deterministic.  A
    # per-pulse random phase (`phase_noise_rad` on either modulator) would
    # break it, exactly as it broke the 8-outcome form in the time-bin chain.
    # Neither modulator is given one here.
    def _response(a_basis, a_bit, b_basis, op_fibre):
        """Run the full field chain once; return the gated (P_x, P_y).

        `op_fibre` is the fibre at the instant being evaluated. `U_comp`
        is deliberately not re-derived from it: Bob calibrated once, and
        the residual `U_comp @ J(t)` is the whole physics of drift.
        """
        if a_basis == 'C':
            v_a = Vpi / 2 if a_bit == 0 else 3 * Vpi / 2
        else:
            v_a = 0 if a_bit == 0 else Vpi

        E = E_source
        E = pm_alice.modulate(E_field=E, V=v_a)
        E = op_fibre.apply(E, dt=dt_field)
        if compensate and U_comp is not None:
            E = np.transpose(U_comp @ np.transpose(E))
        E = optics.voa(E, bob_loss_dB)

        v_b = 0 if b_basis == 'C' else Vpi / 2
        E = pm_bob.modulate(E_field=E, V=v_b)
        # circular_analyser, not pbs: detection depends on the relative
        # phase between Ex/Ey, which a true PBS would be blind to.
        Ex, Ey = optics.circular_analyser(E)
        P_x = float(np.mean(np.abs(Ex) ** 2))
        P_y = float(np.mean(np.abs(Ey) ** 2))
        # Finite analyser extinction, applied to the PORT POWERS and
        # before detection -- see the parameter docstring.  Power
        # conserving, so the sifted rate is untouched and only the error
        # rate moves; a shifted sifted count would mean it had been
        # applied in the wrong place.
        return optics.apply_extinction(P_x, P_y, extinction_epsilon)

    def _build_response(op_fibre):
        return {(ab, bit, bb): _response(ab, bit, bb, op_fibre)
                for ab in ('C', 'X')
                for bit in (0, 1)
                for bb in ('C', 'X')}

    RESPONSE = _build_response(fibre)

    # Sifting is accumulated inline rather than into per-pulse lists.
    # Keeping alice_bits, alice_bases, bob_bits, bob_bases and has_click
    # would be five arrays of O(num_bits): fine at 100k pulses, several GB
    # at the ~30e6 per row a target-sifted sweep needs.  Counting inline is
    # O(1) and touches no RNG draw, so results are unaffected.
    n_sifted = 0
    n_errors = 0
    n_clicks = 0

    # Block boundaries are recorded as CUMULATIVE snapshots and differenced
    # at the end, so the hot loop carries one comparison rather than a pair
    # of extra counters.  Nothing here touches either RNG stream, so a
    # blocked run is bit-identical to an unblocked one.
    blocked = block_size is not None
    if blocked and block_size < 1:
        raise ValueError(f"block_size must be >= 1, got {block_size}")
    snapshots = []

    # --- Fibre drift ------------------------------------------------------
    #
    # `calibration_temperature` / `calibration_bend_radius` give a fixed
    # two-state mismatch: Bob calibrated the fibre in one state, light
    # travels it in another, and neither state moves. Drift is the missing
    # time dependence -- the mismatch GROWS during the run, which is what
    # actually forces the paper's recalibrations.
    #
    # The response table is exact only for one fibre state, so the run is
    # cut into blocks: the fibre is held still within a block and the table
    # rebuilt between them. This partition is INDEPENDENT of `block_size`
    # above, which slices the run for statistical reporting; conflating the
    # two would tie the physics resolution to a reporting choice.
    #
    # `drift_blocks` is a count rather than a pulse size for the same
    # reason `run_duration` exists: raising `num_bits` for tighter error
    # bars must not silently change the simulated experiment.
    #
    # Not drifting means nothing here runs and the table stays exactly as
    # built above.
    _drift_on = fibre.drift_temperature_rate_C_s != 0.0
    _n_dblocks = max(1, int(drift_blocks)) if _drift_on else 1
    if run_duration is not None and num_bits > 1:
        _drift_scale = float(run_duration) / (num_bits - 1)
    else:
        _drift_scale = dt_pulse

    def _dblock_bounds(i):
        return (i * num_bits) // _n_dblocks, ((i + 1) * num_bits) // _n_dblocks

    _dblk = 0
    _dlo, _dblk_end = _dblock_bounds(0)
    if _drift_on:
        RESPONSE = _build_response(
            fibre.at(0.5 * (_dlo + max(_dblk_end - 1, _dlo)) * _drift_scale))

    t_start = time.time()

    for pulse_idx in range(num_bits):
        if _drift_on and pulse_idx >= _dblk_end:
            _dblk += 1
            _dlo, _dblk_end = _dblock_bounds(_dblk)
            RESPONSE = _build_response(
                fibre.at(0.5 * (_dlo + _dblk_end - 1) * _drift_scale))

        # --- Alice's encoding ---
        alice_basis = random.choice(['C', 'X'])
        alice_bit = random.randint(0, 1)

        # --- Bob's decoding basis ---
        bob_basis = random.choice(['C', 'X'])

        # --- Response from the precomputed table (see above) ---
        power_x, power_y = RESPONSE[(alice_basis, alice_bit, bob_basis)]

        # Per-pulse energy from a real source.  Exact rather than an
        # approximation: every stage between the source and the analyser
        # is linear in the field, so an energy factor scales straight
        # through the table.  Cycled, not sampled -- drawing would consume
        # RNG and move every result even at factor 1.0.
        if factors is not None:
            g = factors[pulse_idx % factors.size]
            power_x *= g
            power_y *= g

        t_pulse = pulse_idx * dt_pulse
        click_x = spd1.detect(power_x, t_pulse)
        click_y = spd2.detect(power_y, t_pulse)

        if click_x and not click_y:
            bob_bit = 0
        elif click_y and not click_x:
            bob_bit = 1
        elif click_x and click_y:
            bob_bit = random.randint(0, 1)
        else:
            bob_bit = -1  # no click — discard this pulse

        if click_x or click_y:
            n_clicks += 1
            if alice_basis == bob_basis:
                n_sifted += 1
                if alice_bit != bob_bit:
                    n_errors += 1

        if blocked and (pulse_idx + 1) % block_size == 0:
            snapshots.append((n_sifted, n_errors))

        if verbose and (pulse_idx + 1) % 10000 == 0:
            elapsed = time.time() - t_start
            rate = (pulse_idx + 1) / elapsed
            print(f"  Pulse {pulse_idx+1}/{num_bits} ({rate:.0f} pulses/s)", flush=True)

    # Sifting (same basis AND at least one click) was accumulated in the
    # loop above.
    qber = n_errors / n_sifted if n_sifted > 0 else 0.0

    fiber_loss_dB = alpha_dB * fiber_length
    total_loss_dB = fiber_loss_dB + bob_loss_dB
    sifted_key_rate = n_sifted / (num_bits * dt_pulse) if num_bits > 0 else 0.0
    elapsed = time.time() - t_start

    results = {
        'qber': qber,
        'n_total': num_bits,
        'n_sifted': n_sifted,
        'n_errors': n_errors,
        # n_clicks was computed for the verbose print and thrown away.  It
        # is the RAW key rate once divided by the run duration, and the
        # paper quotes one (265 bit/s), so there is a comparison here that
        # was unreachable purely because the number never left the function.
        'n_clicks': n_clicks,
        'raw_key_rate': n_clicks / (num_bits * dt_pulse) if num_bits > 0 else 0.0,
        'sifted_key_rate': sifted_key_rate,
        'fiber_length_km': fiber_length,
        'total_loss_dB': total_loss_dB,
        'mu': mu,
        'elapsed_s': elapsed,
    }

    if blocked:
        # Differences of the cumulative snapshots, so the blocks partition
        # the run exactly: sum(blocks) equals (n_sifted, n_errors) whenever
        # num_bits is a whole number of blocks.  A trailing partial block is
        # dropped rather than emitted short, because a block with fewer
        # pulses has a different variance and would contaminate any spread
        # measured across them.
        prev_s = prev_e = 0
        blocks = []
        for (cum_s, cum_e) in snapshots:
            blocks.append((cum_s - prev_s, cum_e - prev_e))
            prev_s, prev_e = cum_s, cum_e
        results['blocks'] = blocks
        results['block_size'] = block_size

    if verbose:
        print(f"\nDuplinskiy et al. BB84 — {fiber_length} km")
        print(f"  Polarization compensation: {'ON' if U_comp is not None else 'OFF'}")
        print(f"  Total loss: {total_loss_dB:.1f} dB (fiber {fiber_loss_dB:.1f} + Bob {bob_loss_dB:.1f})")
        print(f"  Mu: {mu} photons/pulse")
        print(f"  Total pulses: {num_bits}")
        print(f"  Clicks: {n_clicks} ({n_clicks/num_bits*100:.2f}%)")
        print(f"  Sifted (same basis + click): {n_sifted} ({n_sifted/num_bits*100:.2f}%)")
        print(f"  Errors: {n_errors}")
        print(f"  QBER: {qber*100:.2f}%")
        print(f"  Time: {elapsed:.1f}s")

    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Replicate Duplinskiy et al. BB84 experiment")
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--bits', type=int, default=100000,
                        help='Number of pulses (default 100k)')
    parser.add_argument('--fiber-length', type=float, default=50,
                        help='Fiber length in km (default 50)')
    parser.add_argument('--mu', type=float, default=0.1,
                        help='Mean photons per pulse (default 0.1)')
    parser.add_argument('--no-compensation', action='store_true',
                        help='Disable Bob polarization compensation (the '
                             'fibre\'s fixed random SU(2) then scrambles '
                             'the encoding)')
    parser.add_argument('--birefringence-model', choices=['auto', 'sectional'],
                        default='auto',
                        help='Birefringence model (default: auto)')
    parser.add_argument('--verbose', action='store_true')
    args = parser.parse_args()

    results = simulate_bb84_duplinskiy(
        num_bits=args.bits, fiber_length=args.fiber_length,
        mu=args.mu, seed=args.seed, verbose=True,
        compensate=not args.no_compensation,
        model=args.birefringence_model)
