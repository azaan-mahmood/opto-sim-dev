"""Fibre impairments in the Gobby chain: the null, the swing, and the drift.

The time-bin chain is often described as blind to birefringence, and for
the **balanced** topology it is: both time bins ride one polarisation, so
a rotation is common to the interfering pair and cancels. The Gobby
replication does not run that topology. It runs the
**polarisation-multiplexed** one the paper describes, where Alice's arms
leave on orthogonal polarisations through a beam combiner and Bob's
splitter routes them back, and there a rotation reaches the arms.

So the same impairment model has to produce three different required
outcomes in one chain, and this measures all three.

What a rotation can do here, and what it cannot
-----------------------------------------------
The fibre's Jones matrix is SU(2), which pins the answer exactly:

    |U00| = |U11|          both interfering arms scaled the SAME
    arg(U00) = -arg(U11)   so their relative phase shifts by 2*arg(U11)

There is therefore no arm imbalance available to collapse the fringe. The
whole effect splits into a common amplitude, which is a pure rate loss,
and a relative phase, which is degenerate with a modulator bias offset.
Both halves of Gobby et al. [1] follow:

    "Polarisation drift reduces the bit rate, but does not degrade the
     QBER provided that the signal rate is significantly higher than the
     intrinsic error rate."

The rate falls as |U00|^2.  The QBER holds because the phase is calibrated
out -- their Bob tunes it with the piezo-driven fibre stretcher in his
long arm -- and the proviso is the background claiming a larger share of a
reduced signal.

Why the field perturbation is reported beside every QBER
--------------------------------------------------------
A null with no field perturbation is a skipped code path wearing the
costume of physics.  `max|dE|/|E|` says the impairment reached the field;
only then does an unchanged QBER mean anything.  The same number appears
in `validate_duplinskiy_birefringence.py`, where the identical
perturbation becomes a 20+ pp swing, which is what makes the pair worth
more than either half.

Two nulls, and they are not the same kind of thing
--------------------------------------------------
* **Balanced, birefringence on** is physics: the rotation is applied and
  the topology absorbs it.
* **Polarisation-multiplexed, alignment on** is arithmetic:
  `U_comp = J^dagger` gives `U_comp @ J = I` for any unitary.  It is
  reported as a consistency check, never as evidence about impairments.

Only drift makes alignment interesting, because then Bob inverts the fibre
as it was while light travels through it as it is.

The drift sweep, and the two amplitudes it is easy to confuse
-------------------------------------------------------------
Two drift rows show the paper's sentence at two points; the sweep shows it
as a curve.  With Bob's servo holding the phase, what drift costs is bit
rate, and the cost grows with the rate while the QBER does not follow.

The prediction for a rate measured across the whole transfer is the
**time-averaged** `|R00|^2`, not its endpoint.  Bob aligns at t=0, so the
residual starts at the identity and walks away from it while pulses are
collected the entire time: at 10 km and 3e-3 C/s the endpoint amplitude is
0.095 and the mean over the run is 0.572, against a measured rate ratio of
0.603 +/- 0.027.  Both amplitudes are printed, because only one of them is
a prediction and they differ by a factor of six.

The swept range stops where the residual rotation completes about one
turn.  Past that `|R00|^2` wraps and rises again, and a curve drawn
through it would report where the rotation happened to land rather than
what drift costs -- the same trap `validate_duplinskiy_drift.py` records
for bend radius.  That the operator is monotone over exactly the swept
range is checked, not assumed.

References
----------
[1] Gobby, C., Yuan, Z. L., & Shields, A. J. (2004). Quantum key
    distribution over 122 km of standard telecom fiber. Appl. Phys.
    Lett. 84(19), 3762-3764.
"""
import argparse
import math
import os
import sys
from collections import namedtuple

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from src.channel.fiber import FiberRealization
from src.protocols.bb84_time_bin import simulate_bb84_time_bin
# The link budget has one home.  Restating it here would be a second
# expression of the same numbers, free to drift out of step with the
# replication this validator is about.
from analysis.val_gobby.validate_gobby import (
    ALPHA_dB, MU, LAM, SPLIT_RATIO, GATE_WIDTH, REP_RATE, PULSE_WIDTH,
    ETA_BOB, P_E, DEAD_TIME, AFTERPULSE_PROB, KEY_TRANSFER_S)

OUT_DIR = os.path.join(os.path.dirname(__file__), '..', 'val_gobby')

SEED = 42

# Two budgets, because two different kinds of claim are being made.
#
# TARGET_SIFTED sizes the cells whose claim is STATISTICAL: the
# uncompensated QBER is checked against a closed form at 4 sigma, and the
# two drift rates have to be told apart from each other and from the
# reference.  1200 sifted puts sigma near 1.3 pp at a 30 % QBER, which
# resolves all three comfortably.
#
# NULL_SIFTED sizes the cells whose claim is EXACT.  The balanced null and
# the aligned null are bit-identity, not agreement within error bars, so
# no budget makes them sharper -- the count only has to be large enough
# that the cell is not vacuously empty.  Sizing them like the statistical
# cells cost about six minutes per run and bought nothing.
TARGET_SIFTED = 1200
NULL_SIFTED = 300

# Distances for the Monte Carlo.  The deterministic section also runs at
# 122 km, because operator algebra costs nothing.
MC_DISTANCES = (10, 50)
OP_DISTANCES = (10, 50, 122)

# Measured sifted yield per pulse, used to budget each cell.  The
# polarisation-multiplexed topology delivers the full mu to one
# interference peak; the balanced one loses half to satellite bins that
# the gate discards, so it yields roughly half as much.
YIELD = {('polarisation_multiplexed', 10): 1.04e-3,
         ('polarisation_multiplexed', 50): 1.77e-4,
         ('balanced', 10): 5.3e-4,
         ('balanced', 50): 9.0e-5}

# Two drift rates, because the paper's claim carries a proviso and the
# point is to show both sides of it.  Neither is a claim about Gobby's
# apparatus, which reports polarisation stable for over 30 minutes at
# 122 km and is quieter than either.
#
#   QUIET  keeps the accumulated phase well under a radian across the
#          120 s transfer, the regime where "reduces the bit rate, but
#          does not degrade the QBER" holds.
#   LOUD   drives it past a radian, where the proviso stops holding and
#          the QBER goes with it.
#
# The crossover is the finding.  A single rate would land on one side of
# it and read as a verdict on the paper rather than on the regime.
#
# LOUD is 3e-3 rather than 1e-3 so that BOTH halves of the paper's sentence
# are visible.  At 1e-3 the amplitude loss is a couple of per cent, inside
# the counting noise, so the rate row could not distinguish a servo that
# corrects phase only from one that corrects everything.  At 3e-3 the mean
# rate factor is 0.57 at 10 km and 0.78 at 50 km, which is unmistakable.
DRIFT_QUIET_C_S = 1e-4
DRIFT_LOUD_C_S = 3e-3

# How often Bob re-locks his phase in the served row.  1.2 s is the default
# drift-block length across the paper's 120 s transfer, so the servo tracks
# at the resolution the model itself resolves.
SERVO_INTERVAL_S = 1.2

# The drift sweep: rate loss as a function of drift rate, servo on.
#
# The two rows above show the paper's sentence at two points.  This shows
# it as a curve -- with the phase held, what drift costs is BIT RATE, and
# the cost grows with the rate.
#
# The range stops at DRIFT_LOUD_C_S, and the reason is measurable rather
# than aesthetic.  |R00|^2 falls monotonically only while the residual
# rotation stays inside about one full turn.  Past that it wraps and comes
# back: at 10 km, 3e-3 C/s gives 0.095 and 1e-2 gives 0.963, with the
# accumulated phase going +4.41 rad to -5.92.  A sweep through that region
# would draw a curve that rises again, and the rise would be an accident of
# where the rotation happened to land -- the same trap
# `validate_duplinskiy_drift.py` records for bend radius, where 1 m is
# worse than 0.1 m.
#
# So the sweep covers the regime where the question is well posed, and the
# monotonicity of the operator over exactly this range is CHECKED below
# rather than assumed.  Both existing rates fall inside it, so the sweep
# extends the two rows rather than replacing them.
DRIFT_SWEEP_C_S = (0.0, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, DRIFT_LOUD_C_S)

# How many times the fibre is re-evaluated across the run.  Named because
# the sweep's prediction has to sample at the same points the model does.
DRIFT_BLOCKS = 100

# One distance for the Monte Carlo half of the sweep.  10 km yields about
# six times what 50 km does per pulse, so the same statistics cost a sixth
# of the runtime, and the claim -- rate falls with drift rate, QBER does
# not -- is not a statement about distance.  The operator half runs at
# every OP_DISTANCE anyway, because it is free.
SWEEP_KM = 10

BUDGET = dict(alpha_dB=ALPHA_dB, mu=MU, wavelength=LAM,
              repetition_rate=REP_RATE, pulse_width=PULSE_WIDTH,
              spad_eta=ETA_BOB, dark_count_rate=P_E / GATE_WIDTH,
              afterpulse_prob=AFTERPULSE_PROB, dead_time=DEAD_TIME,
              gate_width=GATE_WIDTH, split_ratio=SPLIT_RATIO,
              run_duration=KEY_TRANSFER_S, seed=SEED)

REFERENCE = 'no fibre (reference)'
UNCOMP = 'birefringence, uncompensated'
QUIET = 'aligned, drifting (quiet)'
LOUD = 'aligned, drifting (loud)'
SERVO = 'drifting (loud) + phase servo'
ORDER = (REFERENCE, 'birefringence, no align', 'birefringence, aligned',
         UNCOMP, QUIET, LOUD, SERVO)


def budget(km, topo, quick, target=None):
    n = int((target or TARGET_SIFTED) / YIELD[(topo, km)] * 1.10)
    return max(200_000, n // 12) if quick else n


def _run(km, topo, n, **over):
    """One Monte Carlo cell on the Gobby link budget."""
    return simulate_bb84_time_bin(num_bits=n, fiber_length=km,
                                  interferometer=topo, **BUDGET, **over)


def _sigma(q, n):
    return math.sqrt(max(q * (1 - q), 1e-12) / n) if n else float('nan')


def _same(a, b):
    return ((a['qber'], a['n_sifted'], a['n_errors'])
            == (b['qber'], b['n_sifted'], b['n_errors']))


def _fibre(km, rate=0.0):
    return FiberRealization(L_m=km * 1000.0, wavelength=LAM, temperature=25.0,
                            attenuation_factor=ALPHA_dB, attenuation=False,
                            seed=SEED, drift_temperature_rate_C_s=rate)


def _aligned(km, rate):
    """The drifting fibre, and Bob's inverse of it as it stood at t=0.

    One definition of the alignment convention.  Both callers below need
    it, and `mean_residual_amp` needs it hoisted out of a hundred-block
    loop, which is how it would otherwise become two definitions.
    """
    f = _fibre(km, rate=rate)
    return f, f.birefringence_matrix().conj().T


def residual(km, rate, t=KEY_TRANSFER_S):
    """`R = U_comp @ J(t)`: what Bob's t=0 alignment leaves after drift.

    Returns `(|R00|^2, 2*arg(R11))` -- the amplitude the rate follows and
    the phase the QBER follows, which is the whole SU(2) split.  One
    definition, used by both the operator table and the drift sweep, so
    the curve and the predictions it is checked against cannot drift out
    of step with each other.
    """
    f, U = _aligned(km, rate)
    R = U @ f.at(t).birefringence_matrix()
    return abs(R[0, 0]) ** 2, 2.0 * float(np.angle(R[1, 1]))


def mean_residual_amp(km, rate, num_bits, blocks=DRIFT_BLOCKS):
    """Time-averaged `|R00|^2` across the run, sampled as the model does.

    The ENDPOINT is not the prediction for a rate measured over the whole
    transfer, and getting that wrong is easy: Bob aligns at t=0, so the
    residual starts at the identity and walks away from it while pulses
    are collected the entire time.  At 10 km and 3e-3 C/s the endpoint
    amplitude is 0.095 while the mean over the run is 0.57 -- a factor of
    six, and the second one is what a rate ratio can possibly equal.

    Sampled at the same block midpoints `simulate_bb84_time_bin` uses,
    rather than as a continuous integral, so the prediction and the model
    agree about what "during the run" means rather than nearly agreeing.
    """
    if rate == 0.0:
        return 1.0
    f, U = _aligned(km, rate)
    scale = KEY_TRANSFER_S / (num_bits - 1)
    total = 0.0
    for i in range(blocks):
        lo = (i * num_bits) // blocks
        hi = ((i + 1) * num_bits) // blocks
        t_mid = 0.5 * (lo + max(hi - 1, lo)) * scale
        R = U @ f.at(t_mid).birefringence_matrix()
        total += abs(R[0, 0]) ** 2
    return total / blocks


def field_perturbation(km):
    """How hard birefringence hits the field at this distance.

    The common reference between this validator and
    `validate_duplinskiy_birefringence.py`: the same perturbation that the
    balanced topology absorbs to a bit-identical QBER is the one the
    polarisation chain turns into a 20+ pp swing.
    """
    J = _fibre(km).birefringence_matrix()
    if J is None:
        return 0.0
    E = np.array([[1.0 + 0j, 1.0 + 0j]]) / np.sqrt(2)
    return float(np.max(np.abs(np.transpose(J @ np.transpose(E)) - E))
                 / np.max(np.abs(E)))


def operator_section(failures):
    """The deterministic half: what SU(2) forces, before any pulse is sent.

    Every prediction the Monte Carlo is checked against is derived here,
    from the fibre's own operator, so the two are independent statements
    rather than the same arithmetic run twice.
    """
    print("\n  The operator, before any Monte Carlo")
    print("  " + "-" * 74)
    print(f"  {'km':>5} {'max|dE|/|E|':>12} {'||U00|-|U11||':>14} "
          f"{'|argU00+argU11|':>16} {'|U00|^2':>9} {'2argU11':>9} "
          f"{'-> QBER%':>9}")
    preds = {}
    for km in OP_DISTANCES:
        J = _fibre(km).birefringence_matrix()
        amp = abs(J[0, 0]) ** 2
        phi = 2.0 * np.angle(J[1, 1])
        q = (1.0 - math.cos(phi)) / 2.0
        preds[km] = (amp, phi, q)
        d_amp = abs(abs(J[0, 0]) - abs(J[1, 1]))
        d_arg = abs(np.angle(J[0, 0]) + np.angle(J[1, 1]))
        pert = field_perturbation(km)
        print(f"  {km:5d} {pert:12.4f} {d_amp:14.2e} {d_arg:16.2e} "
              f"{amp:9.5f} {phi:+9.4f} {100 * q:9.3f}")
        # SU(2) is the load-bearing fact.  If it ever failed the arms could
        # become unbalanced and every closed form below would be void.
        if d_amp > 1e-12 or d_arg > 1e-12:
            failures.append(
                f"{km} km: the Jones matrix is not SU(2) "
                f"(||U00|-|U11|| = {d_amp:.2e}, |argU00+argU11| = "
                f"{d_arg:.2e}); the amplitude/phase split this validator "
                f"rests on does not hold")
        if pert < 1e-3:
            failures.append(
                f"{km} km: birefringence barely moves the field "
                f"(max|dE|/|E| = {pert:.2e}), so any null below would be a "
                f"skipped code path rather than a result")

    print(f"\n  Bob aligns at t=0, then the fibre walks away from it: "
          f"R = U_comp @ J(t), at t = {KEY_TRANSFER_S:g} s")
    print(f"  {'km':>5} {'dT/dt (C/s)':>12} {'|R00|^2':>9} {'2argR11':>9} "
          f"{'-> QBER%':>9} {'dphi/dt (rad/s)':>17}")
    drift_pred = {}
    for km in OP_DISTANCES:
        for rate in (DRIFT_QUIET_C_S, DRIFT_LOUD_C_S):
            amp, phi = residual(km, rate)
            drift_pred[(km, rate)] = (amp, phi)
            print(f"  {km:5d} {rate:12.0e} {amp:9.5f} "
                  f"{phi:+9.4f} {100 * (1 - math.cos(phi)) / 2:9.3f} "
                  f"{phi / KEY_TRANSFER_S:17.3e}")
    print("    Gobby's measured ARM-length drift, for scale: 0.05 deg/s =")
    print("    8.727e-4 rad/s, so the quiet rate is the comparable one.")
    print("    Note the ordering: the amplitude is still within 1 % of unity")
    print("    while the phase has already moved, because a near-identity")
    print("    residual costs O(eps^2) in amplitude and O(eps) in phase.")
    print("    Fibre drift and arm drift are degenerate in their effect on")
    print("    QBER, which is why fibre drift is off by default in the")
    print("    replication: the paper's 3.3 % floor is already assigned in")
    print("    full to modulator bias plus arm drift, and a third term would")
    print("    count one measurement twice.")
    return preds, drift_pred


# One point on the drift curve.  A plain tuple had grown to six
# fields, three of which are amplitudes, which is how a figure ends
# up plotting the wrong one.
SweepRow = namedtuple('SweepRow',
                      'rate result ratio mean_amp end_amp sd')


def drift_sweep(ref_p, quick, failures):
    """Rate loss against drift rate, with Bob's servo holding the phase.

    The two drift rows show the paper's sentence at two points; this shows
    it as a curve.  With the phase held, what drift costs is bit rate, and
    the cost grows with the rate while the QBER does not follow it.

    Structured the way the rest of this file is: the operator gives the
    prediction, free, and the Monte Carlo is checked against it.  So the
    curve is not merely drawn, it is drawn twice by independent routes.

    Returns the rows, for the CSV and the figure.
    """
    print(f"\n  Drift sweep at {SWEEP_KM} km, servo re-locking every "
          f"{SERVO_INTERVAL_S:g} s")
    print("  " + "-" * 74)

    # The premise first.  The sweep only means anything while |R00|^2 falls
    # with the rate, so that is checked -- at every operator distance,
    # because it costs nothing -- rather than assumed from having looked
    # once.
    for km in OP_DISTANCES:
        amps = [residual(km, r)[0] for r in DRIFT_SWEEP_C_S]
        for (r0, a0), (r1, a1) in zip(zip(DRIFT_SWEEP_C_S, amps),
                                      list(zip(DRIFT_SWEEP_C_S, amps))[1:]):
            if a1 > a0 + 1e-12:
                failures.append(
                    f"{km} km: |R00|^2 rises from {a0:.6f} to {a1:.6f} "
                    f"between {r0:g} and {r1:g} C/s, so the swept range "
                    f"has run past where the residual rotation stays "
                    f"inside one turn.  A curve through a wrapped rotation "
                    f"reports where it happened to land, not what drift "
                    f"costs -- shorten DRIFT_SWEEP_C_S")

    n = budget(SWEEP_KM, 'polarisation_multiplexed', quick)
    print(f"  {'rate C/s':>10} {'sifted':>8} {'rate/ref':>9} "
          f"{'mean|R00|^2':>12} {'end|R00|^2':>11} {'QBER %':>17}")
    rows = []
    for rate in DRIFT_SWEEP_C_S:
        amp = mean_residual_amp(SWEEP_KM, rate, n)
        end_amp, _ = residual(SWEEP_KM, rate)
        if rate == 0.0:
            # Zero drift IS the reference: same budget, same seed, same
            # arithmetic.  Re-running it would spend a cell to recompute a
            # number already in hand.
            r = ref_p
        else:
            r = _run(SWEEP_KM, 'polarisation_multiplexed', n,
                     birefringence=True, compensate=True,
                     drift_temperature_rate_C_s=rate,
                     drift_blocks=DRIFT_BLOCKS,
                     phase_servo_interval_s=SERVO_INTERVAL_S)
        got = r['n_sifted'] / max(ref_p['n_sifted'], 1)
        sig = _sigma(r['qber'], r['n_sifted'])
        print(f"  {rate:10.0e} {r['n_sifted']:8d} {got:9.4f} {amp:12.5f} "
              f"{end_amp:11.5f} {100 * r['qber']:11.3f} +/-{100 * sig:5.3f}")

        # The rate must track the MEAN |R00|^2, which is the servo's whole
        # point: it turns a phase, so the amplitude loss is left exactly
        # where a unitary put it.  Same tolerance shape as the
        # uncompensated cell -- the 15 % floor is systematic, because a
        # cell with little signal has background as a larger share of its
        # sifted count.
        sd = got * math.sqrt(1.0 / max(r['n_sifted'], 1)
                             + 1.0 / max(ref_p['n_sifted'], 1))
        rows.append(SweepRow(rate, r, got, amp, end_amp, sd))
        if abs(got - amp) > max(4.0 * sd, 0.15 * amp):
            failures.append(
                f"{SWEEP_KM} km drift sweep at {rate:g} C/s: rate ratio "
                f"{got:.4f} +/- {sd:.4f} against a mean |R00|^2 of "
                f"{amp:.4f} (endpoint {end_amp:.4f}).  The servo corrects "
                f"phase only, so the amplitude the rotation took is "
                f"exactly what the rate should keep")

        # And the QBER must NOT follow.  This is the half of the paper's
        # sentence the servo exists to reproduce.
        if r['qber'] > ref_p['qber'] + 0.03:
            failures.append(
                f"{SWEEP_KM} km drift sweep at {rate:g} C/s: QBER rose to "
                f"{100 * r['qber']:.2f} % against a "
                f"{100 * ref_p['qber']:.2f} % reference.  Bob re-locks "
                f"every {SERVO_INTERVAL_S:g} s, so the fringe should hold "
                f"whatever the amplitude does")

    # Endpoints, not adjacent steps.  The first few predicted separations
    # are ~1e-4 in ratio against counting noise of order 1e-2, so adjacent
    # monotonicity is not resolvable and asserting it would fail at random.
    # The per-point checks above already pin the shape; this pins that
    # there IS a fall.
    lo, hi = rows[0].ratio, rows[-1].ratio
    if hi >= lo - 0.10:
        failures.append(
            f"{SWEEP_KM} km drift sweep: the rate did not fall across the "
            f"range -- {lo:.4f} at 0 C/s against {hi:.4f} at "
            f"{DRIFT_SWEEP_C_S[-1]:g} C/s.  With the phase held, bit rate "
            f"is the thing drift is supposed to cost")

    print("    Adjacent steps are reported, not asserted: the first few are")
    print("    ~1e-4 in ratio against counting noise of order 1e-2.")
    return rows


def run(quick=False):
    os.makedirs(OUT_DIR, exist_ok=True)
    failures = []
    print("=" * 78)
    print("Gobby: fibre impairments -- balanced null, polmux swing, and drift")
    print("=" * 78)
    print(f"  seed {SEED}, target {TARGET_SIFTED} sifted per Monte Carlo cell"
          + ("   [QUICK -- not quotable]" if quick else ""))

    preds, drift_pred = operator_section(failures)

    print("\n  Monte Carlo")
    print("  " + "-" * 74)
    print(f"  {'km':>4} {'topology':<10} {'configuration':<28} {'sifted':>7} "
          f"{'rate/ref':>9} {'QBER %':>18}")
    rows = []
    sweep_ref = None
    for km in MC_DISTANCES:
        amp_pred, _, q_pred = preds[km]
        cells = []

        # --- balanced: the physics null ------------------------------
        nb = budget(km, 'balanced', quick, target=NULL_SIFTED)
        ref_b = _run(km, 'balanced', nb)
        bir_b = _run(km, 'balanced', nb, birefringence=True, compensate=False)
        cells += [('balanced', REFERENCE, ref_b, ref_b),
                  ('balanced', 'birefringence, no align', bir_b, ref_b)]
        if not _same(bir_b, ref_b):
            failures.append(
                f"{km} km balanced: birefringence changed the result "
                f"(QBER {100 * ref_b['qber']:.3f} -> {100 * bir_b['qber']:.3f} "
                f"%).  Both time bins share one polarisation here, so a "
                f"rotation is common to the interfering pair and must cancel")

        # --- polarisation-multiplexed --------------------------------
        npm = budget(km, 'polarisation_multiplexed', quick)
        ref_p = _run(km, 'polarisation_multiplexed', npm)
        if km == SWEEP_KM:
            sweep_ref = ref_p
        aligned = _run(km, 'polarisation_multiplexed', npm,
                       birefringence=True, compensate=True)
        uncomp = _run(km, 'polarisation_multiplexed', npm,
                      birefringence=True, compensate=False)
        quiet = _run(km, 'polarisation_multiplexed', npm,
                     birefringence=True, compensate=True,
                     drift_temperature_rate_C_s=DRIFT_QUIET_C_S,
                     drift_blocks=DRIFT_BLOCKS)
        loud = _run(km, 'polarisation_multiplexed', npm,
                    birefringence=True, compensate=True,
                    drift_temperature_rate_C_s=DRIFT_LOUD_C_S,
                    drift_blocks=DRIFT_BLOCKS)
        servo = _run(km, 'polarisation_multiplexed', npm,
                     birefringence=True, compensate=True,
                     drift_temperature_rate_C_s=DRIFT_LOUD_C_S,
                     drift_blocks=DRIFT_BLOCKS,
                     phase_servo_interval_s=SERVO_INTERVAL_S)
        cells += [('polmux', REFERENCE, ref_p, ref_p),
                  ('polmux', 'birefringence, aligned', aligned, ref_p),
                  ('polmux', UNCOMP, uncomp, ref_p),
                  ('polmux', QUIET, quiet, ref_p),
                  ('polmux', LOUD, loud, ref_p),
                  ('polmux', SERVO, servo, ref_p)]

        if not _same(aligned, ref_p):
            failures.append(
                f"{km} km polmux: alignment did not return an exact null.  "
                f"U_comp = J^dagger gives U_comp @ J = I for any unitary, so "
                f"this is arithmetic and a deviation is a bug, not physics")

        # Uncompensated must match BOTH closed forms derived above.
        got_rate = uncomp['n_sifted'] / max(ref_p['n_sifted'], 1)
        # Tolerance from the counts rather than a flat fraction, so a smoke
        # run relaxes it honestly instead of failing on noise.  The 15 %
        # floor is systematic, not statistical: the uncompensated cell has
        # far less signal, so background clicks are a larger share of its
        # sifted count and push the ratio slightly above |U00|^2.
        sd_ratio = got_rate * math.sqrt(1.0 / max(uncomp['n_sifted'], 1)
                                        + 1.0 / max(ref_p['n_sifted'], 1))
        if abs(got_rate - amp_pred) > max(4.0 * sd_ratio, 0.15 * amp_pred):
            failures.append(
                f"{km} km polmux uncompensated: rate ratio {got_rate:.4f} "
                f"+/- {sd_ratio:.4f} against |U00|^2 = {amp_pred:.4f}.  The "
                f"amplitude term is common to both arms, so that ratio is "
                f"the prediction")
        s = _sigma(uncomp['qber'], uncomp['n_sifted'])
        if abs(uncomp['qber'] - q_pred) > 4.0 * s:
            failures.append(
                f"{km} km polmux uncompensated: QBER "
                f"{100 * uncomp['qber']:.2f} +/- {100 * s:.2f} % against "
                f"(1-cos 2argU11)/2 = {100 * q_pred:.2f} %, a "
                f"{abs(uncomp['qber'] - q_pred) / s:.1f} sigma miss")

        # Three drift rows, and together they are the paper's sentence and
        # its proviso.
        #
        #   quiet  the residual stays small, so nothing moves: the regime
        #          where "reduces the bit rate, but does not degrade the
        #          QBER" holds on its own.
        #   loud   the residual grows past a radian.  BOTH halves go: the
        #          rate falls with the amplitude and the QBER follows the
        #          phase.  This is the model without Bob's stretcher, and
        #          it contradicts the paper.
        #   servo  the same drift with Bob holding his operating point.
        #          The QBER comes back and the rate does NOT, which is the
        #          paper's ordering restored.
        #
        # The last comparison is the one that carries weight, and it only
        # carries it because loud is harsh enough for the rate loss to sit
        # outside the counting noise.
        q_rate = quiet['n_sifted'] / max(ref_p['n_sifted'], 1)
        l_rate = loud['n_sifted'] / max(ref_p['n_sifted'], 1)
        s_rate = servo['n_sifted'] / max(ref_p['n_sifted'], 1)
        if q_rate < 0.85:
            failures.append(
                f"{km} km polmux {QUIET}: the rate fell to {q_rate:.3f} of "
                f"reference.  A residual this close to the identity costs "
                f"O(eps^2) in amplitude and cannot cost that much")
        if l_rate > 0.85:
            failures.append(
                f"{km} km polmux {LOUD}: the rate only fell to {l_rate:.3f}, "
                f"against a predicted {drift_pred[(km, DRIFT_LOUD_C_S)][0]:.3f}. "
                f"Without a clear rate loss here the servo row below cannot "
                f"show that the servo leaves it alone")
        # The servo turns a phase, so it cannot return photons.  If the
        # rate recovers, it is correcting amplitude and is the wrong model.
        if s_rate > l_rate + 0.10:
            failures.append(
                f"{km} km polmux {SERVO}: the rate recovered from "
                f"{l_rate:.3f} to {s_rate:.3f}.  A phase correction cannot "
                f"do that -- the servo is touching the amplitude")
        if servo['qber'] > ref_p['qber'] + 0.03:
            failures.append(
                f"{km} km polmux {SERVO}: QBER stayed at "
                f"{100 * servo['qber']:.2f} % against a "
                f"{100 * ref_p['qber']:.2f} % reference.  Bob is re-locking "
                f"every {SERVO_INTERVAL_S:g} s, which should hold the fringe")
        if quiet['qber'] > ref_p['qber'] + 0.03:
            failures.append(
                f"{km} km polmux {QUIET}: QBER moved to "
                f"{100 * quiet['qber']:.2f} % at {DRIFT_QUIET_C_S:g} C/s, "
                f"where the accumulated phase only reaches "
                f"{drift_pred[(km, DRIFT_QUIET_C_S)][1]:+.3f} rad.  Gobby's "
                f"'does not degrade the QBER' should still hold in this "
                f"regime")
        if loud['qber'] <= ref_p['qber'] + 0.03:
            failures.append(
                f"{km} km polmux {LOUD}: QBER did not move "
                f"({100 * loud['qber']:.2f} % against a "
                f"{100 * ref_p['qber']:.2f} % reference) even though the "
                f"accumulated phase reaches "
                f"{drift_pred[(km, DRIFT_LOUD_C_S)][1]:+.3f} rad.  Bob "
                f"aligned at t=0, so the residual cannot still be the "
                f"identity")

        for topo, lab, r, ref in cells:
            print(f"  {km:4d} {topo:<10} {lab:<28} {r['n_sifted']:>7} "
                  f"{r['n_sifted'] / max(ref['n_sifted'], 1):>9.4f} "
                  f"{100 * r['qber']:>11.3f} +/-"
                  f"{100 * _sigma(r['qber'], r['n_sifted']):5.3f}")
            # The floor a cell must clear depends on what it claims: the
            # balanced cells assert bit-identity, which no budget sharpens.
            target = NULL_SIFTED if topo == 'balanced' else TARGET_SIFTED
            # A smoke run spends a twelfth of the pulses, so its floor has
            # to scale with the budget rather than sit at a fixed count.
            floor = max(20, target // 24) if quick else target // 2
            # The uncompensated cell is EXPECTED to lose most of its rate --
            # that is the measurement -- so it is exempt from the floor.
            if r['n_sifted'] < floor and lab != UNCOMP:
                failures.append(f"{km} km {topo} {lab}: only {r['n_sifted']} "
                                f"sifted, below the quotable standard")
            rows.append((km, topo, lab, r, ref))
        print(f"  {'':4} {'':10} {'predicted, uncompensated':<28} {'':>7} "
              f"{amp_pred:>9.4f} {100 * q_pred:>11.3f}")

    sweep = drift_sweep(sweep_ref, quick, failures)

    _write_csv(rows, preds, sweep, quick)
    _figure(rows, preds, sweep, quick)

    print()
    if failures:
        print("[FAIL]")
        for f in failures:
            print(f"  - {f}")
        return 1
    print("[PASS] balanced absorbs the rotation to a bit-identical QBER")
    print("[PASS] polmux alignment returns the exact arithmetic null")
    print("[PASS] uncompensated matches |U00|^2 in rate and "
          "(1-cos 2argU11)/2 in QBER")
    print("[PASS] quiet drift moves neither rate nor QBER")
    print("[PASS] loud drift costs both, which is the model WITHOUT Bob's")
    print("       stretcher and is not what the paper reports")
    print("[PASS] the phase servo returns the QBER and leaves the rate lost")
    print("       -- 'reduces the bit rate, but does not degrade the QBER'")
    print("[PASS] |R00|^2 falls monotonically across the swept range at "
          "every")
    print("       distance, so the sweep stays inside one turn of the "
          "residual")
    print(f"[PASS] swept {len(DRIFT_SWEEP_C_S)} drift rates with the servo "
          f"on: the rate tracks")
    print("       |R00|^2 at every one, and the QBER tracks none of them")
    return 0


def _stem(quick):
    """Smoke runs write to their own files, so an under-powered run can
    never replace a quotable artifact."""
    return 'val_gobby_impairments--quick' if quick else 'val_gobby_impairments'


def _write_csv(rows, preds, sweep, quick=False):
    path = os.path.join(OUT_DIR, _stem(quick) + '.csv')
    with open(path, 'w') as fh:
        fh.write("# validate_gobby_impairments.py"
                 + (" --quick  [NOT QUOTABLE]\n" if quick else "\n"))
        fh.write(f"# seed={SEED} target_sifted={TARGET_SIFTED} "
                 f"drift_quiet_C_s={DRIFT_QUIET_C_S:g} "
                 f"drift_loud_C_s={DRIFT_LOUD_C_S:g} "
                 f"servo_interval_s={SERVO_INTERVAL_S:g} "
                 f"run_duration_s={KEY_TRANSFER_S:g}\n")
        fh.write("# Link budget imported from val_gobby/validate_gobby.py.\n")
        fh.write("# predicted_rate = |U00|^2 and predicted_qber = "
                 "(1-cos(2*arg(U11)))/2, both derived from the fibre's own\n"
                 "# Jones matrix, which is SU(2).  They apply to the "
                 "uncompensated cell only.\n")
        fh.write(f"# 'servo sweep' rows are the drift-rate curve at "
                 f"{SWEEP_KM} km, servo on.  predicted_rate there is "
                 f"|R00|^2 for\n"
                 f"# R = U_comp @ J(t) AVERAGED over the run at the "
                 f"model's {DRIFT_BLOCKS} block midpoints -- not its "
                 f"endpoint, which is\n"
                 f"# smaller by a factor of six at the loudest rate "
                 f"because Bob aligns at t=0 and the residual grows from\n"
                 f"# the identity while pulses are collected throughout.  "
                 f"predicted_qber is EMPTY on purpose: the servo holds\n"
                 f"# the fringe, so the unserved phase term is not the "
                 f"prediction for those rows.  Swept rates: "
                 + ' '.join(f'{r:g}' for r in DRIFT_SWEEP_C_S) + "\n")
        fh.write("distance_km,topology,configuration,pulses,sifted,errors,"
                 "qber,sigma,rate_vs_reference,predicted_rate,"
                 "predicted_qber,field_perturbation\n")
        for km, topo, lab, r, ref in rows:
            amp, _, q = preds[km]
            pr = f"{amp:.6f}" if lab == UNCOMP else ""
            pq = f"{q:.6f}" if lab == UNCOMP else ""
            # The labels contain commas, so the field must be quoted or
            # every column after it shifts by one.
            fh.write(f"{km},{topo},\"{lab}\",{r['n_total']},{r['n_sifted']},"
                     f"{r['n_errors']},{r['qber']:.6f},"
                     f"{_sigma(r['qber'], r['n_sifted']):.6f},"
                     f"{r['n_sifted'] / max(ref['n_sifted'], 1):.6f},"
                     f"{pr},{pq},{field_perturbation(km):.6f}\n")
        # The drift sweep, same columns.  `predicted_rate` is |R00|^2 at
        # that rate; `predicted_qber` is left empty on purpose -- the
        # servo holds the fringe, so the unserved (1-cos(2argR11))/2 is
        # NOT the prediction for these rows and printing it would invite
        # exactly the wrong comparison.
        for row in sweep:
            rate, r = row.rate, row.result
            lab = f"servo sweep, drift {rate:g} C/s"
            fh.write(f"{SWEEP_KM},polarisation_multiplexed,\"{lab}\","
                     f"{r['n_total']},{r['n_sifted']},{r['n_errors']},"
                     f"{r['qber']:.6f},"
                     f"{_sigma(r['qber'], r['n_sifted']):.6f},"
                     f"{row.ratio:.6f},{row.mean_amp:.6f},,"
                     f"{field_perturbation(SWEEP_KM):.6f}\n")
    print(f"\n  CSV: {path}")


def _sweep_panel(ax, sweep):
    """The drift-rate curve: what the servo saves and what it cannot.

    Both halves of the paper's sentence in one axes.  Rate falls with the
    drift rate and follows the mean |R00|^2 the operator predicts; QBER,
    on the right axis, does not follow -- which is the servo working.

    The prediction is drawn at the swept rates and nowhere between them.
    A smooth interpolating curve through a finer grid would look better
    and would promise values that nothing here checked.

    Zero goes on a log axis at the left edge, labelled, because the
    alternative is either dropping the reference point or drawing the
    whole curve linearly and losing four decades.
    """
    nz = [r for r in DRIFT_SWEEP_C_S if r > 0]
    floor = min(nz) / 3.0
    x = [floor if w.rate == 0.0 else w.rate for w in sweep]

    # The prediction is drawn at the swept rates and nowhere else.  A
    # smooth interpolating curve would promise values between them that
    # nothing here checked.
    ax.plot(x, [w.mean_amp for w in sweep],
            '-', lw=1.6, alpha=0.6, color='tab:purple', zorder=1,
            label=r'mean $|R_{00}|^2$, from the operator')
    ax.errorbar(x, [w.ratio for w in sweep], yerr=[w.sd for w in sweep],
                fmt='o', ms=6, capsize=3, lw=1.0, color='tab:purple',
                markeredgecolor='0.25', zorder=3,
                label='sifted rate / reference')
    ax.set_xscale('log')
    ax.set_xlim(floor / 1.6, max(nz) * 1.6)
    ax.set_ylim(0, 1.15)
    ax.set_xlabel('fibre drift rate (C/s)')
    ax.set_ylabel('sifted rate / reference')
    ax.grid(True, alpha=0.3)
    # Name the fake zero rather than letting it read as 3e-6.
    ticks = [floor] + nz
    ax.set_xticks(ticks)
    ax.set_xticklabels(['0'] + [f'{r:g}' for r in nz], fontsize=7)

    q = ax.twinx()
    qs = [100 * w.result['qber'] for w in sweep]
    errs = [100 * _sigma(w.result['qber'], w.result['n_sifted'])
            for w in sweep]
    q.errorbar(x, qs, yerr=errs, fmt='s--', ms=4, lw=1.0, capsize=2,
               color='tab:olive', label='QBER, servo holding')
    q.set_ylabel('QBER (%)', color='tab:olive')
    q.tick_params(axis='y', labelcolor='tab:olive')
    q.set_ylim(0, max(6.0, 1.6 * max(qs + [1.0])))

    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = q.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=7, loc='lower left', framealpha=0.9)
    ax.set_title(f'Drift rate vs rate loss at {SWEEP_KM} km, servo on',
                 fontsize=10)


def _figure(rows, preds, sweep, quick=False):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("  matplotlib missing, skipping the figure")
        return

    colour = {REFERENCE: '0.55',
              'birefringence, no align': 'tab:green',
              'birefringence, aligned': 'tab:blue',
              UNCOMP: 'tab:red',
              QUIET: 'tab:olive',
              LOUD: 'tab:orange',
              SERVO: 'tab:purple'}
    kms = sorted({r[0] for r in rows})
    fig, axes = plt.subplots(1, 3, figsize=(17.0, 5.2))

    for ax, metric, ylab in ((axes[0], 'qber', 'QBER (%)'),
                             (axes[1], 'rate', 'sifted rate / reference')):
        x = np.arange(len(kms))
        w = 0.16
        for i, lab in enumerate(ORDER):
            vals, errs = [], []
            for km in kms:
                m = [r for r in rows if r[0] == km and r[2] == lab]
                if not m:
                    vals.append(np.nan)
                    errs.append(0.0)
                    continue
                _, _, _, r, ref = m[0]
                if metric == 'qber':
                    vals.append(100 * r['qber'])
                    errs.append(100 * _sigma(r['qber'], r['n_sifted']))
                else:
                    vals.append(r['n_sifted'] / max(ref['n_sifted'], 1))
                    errs.append(0.0)
            ax.bar(x + (i - 3.0) * w, vals, w, yerr=errs, capsize=3,
                   color=colour[lab], edgecolor='0.25', linewidth=0.6,
                   label=lab if ax is axes[0] else None)
        # The closed forms, drawn as targets rather than as fitted lines.
        for j, km in enumerate(kms):
            amp, _, q = preds[km]
            target = 100 * q if metric == 'qber' else amp
            # Centred on the UNCOMPENSATED bar, which is index 3 of six
            # drawn at x + (i - 2.5)*w, i.e. x + 0.5*w.  The prediction
            # applies to that cell alone, so it has to sit over it.
            ax.plot([x[j] - w * 0.5, x[j] + w * 0.5], [target, target],
                    color='k', lw=2.2, zorder=5,
                    label=('predicted from the Jones matrix'
                           if (ax is axes[0] and j == 0) else None))
        ax.set_xticks(x)
        ax.set_xticklabels([f'{k} km' for k in kms])
        ax.set_ylabel(ylab)
        ax.grid(True, axis='y', alpha=0.3)
        # Fixed limits so runs are comparable by eye, and so the tallest
        # bar is never the one the axes decide to clip.
        ax.set_ylim(0, 100 if metric == 'qber' else 1.15)

    axes[0].set_title('QBER: balanced absorbs the rotation, polmux does not',
                      fontsize=10)
    axes[1].set_title('Sifted rate: what a rotation actually costs',
                      fontsize=10)
    _sweep_panel(axes[2], sweep)
    # Figure-level legend below both panels.  Inside the axes it sat on top
    # of the tallest bar and its prediction marker, which are the two things
    # the reader most needs to see.
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, fontsize=8, ncol=4, loc='lower center',
               frameon=False, bbox_to_anchor=(0.5, -0.01))
    fig.suptitle(f'Fibre impairments in the Gobby chain (seed {SEED}, '
                 f'{KEY_TRANSFER_S:g} s transfer, servo every '
                 f'{SERVO_INTERVAL_S:g} s; bars at drift '
                 f'{DRIFT_QUIET_C_S:g} and {DRIFT_LOUD_C_S:g} C/s, '
                 f'sweep over {len(DRIFT_SWEEP_C_S)} rates to '
                 f'{max(DRIFT_SWEEP_C_S):g})'
                 + ('   [QUICK -- not quotable]' if quick else ''),
                 fontsize=11)
    fig.tight_layout(rect=(0, 0.09, 1, 1))
    png = os.path.join(OUT_DIR, _stem(quick) + '.png')
    fig.savefig(png, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  figure: {png}")


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--quick', action='store_true',
                    help='a twelfth of the pulses, for a smoke run; writes '
                         'to --quick names so it cannot replace the quotable '
                         'artifact')
    a = ap.parse_args()
    sys.exit(run(quick=a.quick))
