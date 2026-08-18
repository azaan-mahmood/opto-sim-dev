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
DRIFT_QUIET_C_S = 1e-4
DRIFT_LOUD_C_S = 1e-3

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
ORDER = (REFERENCE, 'birefringence, no align', 'birefringence, aligned',
         UNCOMP, QUIET, LOUD)


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
            f = _fibre(km, rate=rate)
            U = f.birefringence_matrix().conj().T
            R = U @ f.at(KEY_TRANSFER_S).birefringence_matrix()
            phi = 2.0 * np.angle(R[1, 1])
            drift_pred[(km, rate)] = (abs(R[0, 0]) ** 2, phi)
            print(f"  {km:5d} {rate:12.0e} {abs(R[0, 0]) ** 2:9.5f} "
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
        aligned = _run(km, 'polarisation_multiplexed', npm,
                       birefringence=True, compensate=True)
        uncomp = _run(km, 'polarisation_multiplexed', npm,
                      birefringence=True, compensate=False)
        quiet = _run(km, 'polarisation_multiplexed', npm,
                     birefringence=True, compensate=True,
                     drift_temperature_rate_C_s=DRIFT_QUIET_C_S,
                     drift_blocks=100)
        loud = _run(km, 'polarisation_multiplexed', npm,
                    birefringence=True, compensate=True,
                    drift_temperature_rate_C_s=DRIFT_LOUD_C_S,
                    drift_blocks=100)
        cells += [('polmux', REFERENCE, ref_p, ref_p),
                  ('polmux', 'birefringence, aligned', aligned, ref_p),
                  ('polmux', UNCOMP, uncomp, ref_p),
                  ('polmux', QUIET, quiet, ref_p),
                  ('polmux', LOUD, loud, ref_p)]

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

        # The paper's proviso, both sides of it.  Same impairment, same
        # alignment, two rates: quiet keeps the accumulated phase well under
        # a radian and the QBER holds, loud pushes it past one and the QBER
        # goes with it.  The rate survives in BOTH, which is the asymmetry
        # worth pinning -- a near-identity residual costs O(eps^2) in
        # amplitude against O(eps) in phase, so QBER is always what breaks
        # first.
        for lab, r, rate in ((QUIET, quiet, DRIFT_QUIET_C_S),
                             (LOUD, loud, DRIFT_LOUD_C_S)):
            amp_r, _ = drift_pred[(km, rate)]
            got = r['n_sifted'] / max(ref_p['n_sifted'], 1)
            if got < 0.7:
                failures.append(
                    f"{km} km polmux {lab}: the rate fell to {got:.3f} of "
                    f"reference against a predicted {amp_r:.3f}.  A residual "
                    f"this close to the identity cannot cost that much "
                    f"amplitude")
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

    _write_csv(rows, preds, quick)
    _figure(rows, preds, quick)

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
    print("[PASS] drift spares the rate at both rates, and degrades QBER")
    print("       only once the accumulated phase passes a radian -- the")
    print("       proviso in Gobby's own sentence, measured not assumed")
    return 0


def _stem(quick):
    """Smoke runs write to their own files, so an under-powered run can
    never replace a quotable artifact."""
    return 'val_gobby_impairments--quick' if quick else 'val_gobby_impairments'


def _write_csv(rows, preds, quick=False):
    path = os.path.join(OUT_DIR, _stem(quick) + '.csv')
    with open(path, 'w') as fh:
        fh.write("# validate_gobby_impairments.py"
                 + (" --quick  [NOT QUOTABLE]\n" if quick else "\n"))
        fh.write(f"# seed={SEED} target_sifted={TARGET_SIFTED} "
                 f"drift_quiet_C_s={DRIFT_QUIET_C_S:g} "
                 f"drift_loud_C_s={DRIFT_LOUD_C_S:g} "
                 f"run_duration_s={KEY_TRANSFER_S:g}\n")
        fh.write("# Link budget imported from val_gobby/validate_gobby.py.\n")
        fh.write("# predicted_rate = |U00|^2 and predicted_qber = "
                 "(1-cos(2*arg(U11)))/2, both derived from the fibre's own\n"
                 "# Jones matrix, which is SU(2).  They apply to the "
                 "uncompensated cell only.\n")
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
    print(f"\n  CSV: {path}")


def _figure(rows, preds, quick=False):
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
              LOUD: 'tab:orange'}
    kms = sorted({r[0] for r in rows})
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 5.2))

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
            ax.bar(x + (i - 2.5) * w, vals, w, yerr=errs, capsize=3,
                   color=colour[lab], edgecolor='0.25', linewidth=0.6,
                   label=lab if ax is axes[0] else None)
        # The closed forms, drawn as targets rather than as fitted lines.
        for j, km in enumerate(kms):
            amp, _, q = preds[km]
            target = 100 * q if metric == 'qber' else amp
            # Centred on the UNCOMPENSATED bar, which is index 3 of six
            # drawn at x + (i - 2.5)*w, i.e. x + 0.5*w.  The prediction
            # applies to that cell alone, so it has to sit over it.
            ax.plot([x[j], x[j] + w], [target, target],
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
    # Figure-level legend below both panels.  Inside the axes it sat on top
    # of the tallest bar and its prediction marker, which are the two things
    # the reader most needs to see.
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, fontsize=8, ncol=4, loc='lower center',
               frameon=False, bbox_to_anchor=(0.5, -0.01))
    fig.suptitle(f'Fibre impairments in the Gobby chain (seed {SEED}, '
                 f'drift {DRIFT_QUIET_C_S:g} and {DRIFT_LOUD_C_S:g} C/s over '
                 f'{KEY_TRANSFER_S:g} s)'
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
