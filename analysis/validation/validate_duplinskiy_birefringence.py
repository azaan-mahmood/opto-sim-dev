"""The birefringence swing: the other half of a Level-3 pair.

The **null** side is measured by `validate_gobby_impairments.py`:
birefringence perturbs the field hard, max|dE|/|E| of order 1-2, while the
QBER of the BALANCED topology stays bit-identical -- same qber, same
sifted count, same error count.  That is what a phase-encoded,
path-matched interferometer should do, because both time bins ride one
polarisation and the rotation is common to the interfering pair.

The qualifier matters.  The Gobby *replication* runs the
polarisation-multiplexed topology, where the arms leave on orthogonal
polarisations and a rotation does reach them -- as a common amplitude
|U00|^2 and a relative phase 2*arg(U11), the Jones matrix being SU(2).
So "the time-bin chain is blind to birefringence" is a statement about a
topology, not about an encoding.

But a null alone cannot validate an impairment model.  says so
itself, and names the successor:

    "It is not a validation of the impairment models at the QBER level.
     The observable is constitutionally blind to these three... Testing
     them against an observable requires an encoding they actually
     degrade: polarisation encoding, where an uncompensated SU(2) rotation
     maps directly onto the bit.  The Duplinskiy chain is that chain, and
     it carries a built-in control in its `compensate` flag.  A first
     probe at low statistics shows the swing is real -- compensated
     ~0-4.8 %, uncompensated 14-50 % across 10/50/100 km -- but the counts
     (2-50 sifted) are far too thin to quote."

This runs that probe properly.

What a Level-3 pair is
----------------------
The **same unchanged impairment model** must produce **opposite required
outcomes** in two chains:

    birefringence   ->  exact NULL in Gobby      (, done)
                    ->  large SWING in Duplinskiy (here)

A model that silently no-ops passes the null and fails the swing.  A model
that over-applies passes the swing and breaks the null.  Neither chain
alone can catch either failure, which is why the pair is worth more than
the sum of its halves.

Scope: birefringence only
-------------------------
also lists CD and PMD, and notes they were "hardcoded off"
in this chain.   exposed both -- but then showed they
are **inert here for a different reason**: this chain evaluates a single
time sample, and CD and PMD act on a time-resolved field.  So they cannot
be tested by this sweep, and claiming otherwise would repeat exactly the
mistake caught.  Birefringence is the one of the three that a
single-sample polarisation observable can see.

References
----------
[1] Duplinskiy, Ustimchik, Kanapin, Kurochkin & Kurochkin, Opt. Express
    25(23), 28886-28897 (2017).
"""
import argparse
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from src.channel import FiberRealization
from src.protocols.bb84_duplinskiy import simulate_bb84_duplinskiy

OUT_DIR = os.path.join(os.path.dirname(__file__), '..', 'val_duplinskiy')

SEED = 42
TARGET_SIFTED = 3000

# Measured sifted yield per pulse, used to budget each cell to the # standard.  Throughput is ~430k pulses/s after 's 8-outcome
# precompute, so even the 100 km cells cost about four minutes each --
# an earlier estimate of two hours used the pre- figure of 30k/s.
YIELD = {10: 1.72e-3, 50: 3.38e-4, 100: 3.25e-5}
DISTANCES = (10, 50, 100)


def budget(km, quick):
    n = int(TARGET_SIFTED / YIELD[km] * 1.10)
    return max(200_000, n // 12) if quick else n


def field_perturbation(km):
    """How hard birefringence actually hits the field at this distance.

    needed this on the Gobby side to show its null was
    physics rather than a config key being ignored.  The same number is
    worth having here, as the common reference: it is the *same*
    perturbation that this chain turns into a 20+ pp swing and Gobby
    turns into nothing.
    """
    fib = FiberRealization(L_m=km * 1000, temperature=25, bend_radius=None,
                           attenuation_factor=0.2, cd=False, pmd=False,
                           model='auto', seed=SEED)
    J = fib.birefringence_matrix()
    if J is None:
        return 0.0
    E = np.array([[1.0 + 0j, 1.0 + 0j]]) / np.sqrt(2)
    E_out = np.transpose(J @ np.transpose(E))
    return float(np.max(np.abs(E_out - E)) / np.max(np.abs(E)))


def run(quick=False):
    os.makedirs(OUT_DIR, exist_ok=True)
    failures = []
    print("=" * 74)
    print("Duplinskiy: the birefringence swing, against Gobby's null")
    print("=" * 74)
    print(f"  seed {SEED}, target {TARGET_SIFTED} sifted per cell "
          f"(standard)")
    print("  The null side is measured by "
          "analysis/validation/validate_gobby_impairments.py,")
    print("  not transcribed here.  What it establishes is bit-identity in")
    print("  the BALANCED topology -- same qber, same sifted, same errors --")
    print("  under a field perturbation of the same size this chain sees.")
    print("  A QBER value is the wrong thing to quote for it: the claim is")
    print("  exactness, and a number would invite comparison against error")
    print("  bars that do not apply to it.")

    print(f"\n  {'km':>5} {'pulses':>12} {'compensated':>22} "
          f"{'uncompensated':>22} {'swing':>20}")
    rows = []
    for km in DISTANCES:
        n = budget(km, quick)
        out = {}
        for comp in (True, False):
            r = simulate_bb84_duplinskiy(n, fiber_length=km, compensate=comp,
                                         seed=SEED)
            s, e, q = r['n_sifted'], r['n_errors'], r['qber']
            sig = math.sqrt(max(q * (1 - q), 1e-12) / s) if s else float('nan')
            out[comp] = (q, sig, s, e)
            if s < (250 if quick else TARGET_SIFTED):
                failures.append(f"{km} km, compensate={comp}: only {s} sifted; "
                                "below the quotable standard")
        (qc, sc, nc, _), (qu, su, nu, _) = out[True], out[False]
        d = qu - qc
        ds = math.hypot(sc, su)
        rows.append((km, n, qc, sc, nc, qu, su, nu, d, ds))
        print(f"  {km:5d} {n:12,} "
              f"{100 * qc:8.2f} +/-{100 * sc:5.2f} % ({nc:5d}) "
              f"{100 * qu:8.2f} +/-{100 * su:5.2f} % ({nu:5d}) "
              f"{100 * d:+8.2f} +/-{100 * ds:5.2f} pp")
        if d < 5 * ds:
            failures.append(f"{km} km: the swing is only {d / ds:.1f} sigma; "
                            "birefringence is not reaching the observable, so "
                            "the Level-3 pair is not demonstrated")

    print("\n  the same perturbation, measured at the field level")
    for km in DISTANCES:
        print(f"    {km:5d} km   max|dE|/|E| = {field_perturbation(km):.3f}")
    print("    -> this is what Gobby absorbs to bit-identical QBER and what")
    print("       this chain turns into the swing above. Same model, same")
    print("       magnitude, opposite required outcomes.")

    _figure(rows, quick)
    _write_csv(rows, quick)

    print()
    if failures:
        print("[FAIL]")
        for f in failures:
            print(f"  - {f}")
        return 1
    print("[PASS] birefringence swings the polarisation observable at every")
    print("       distance, completing the pair against null")
    return 0


def _stem(quick):
    """Smoke runs write to their own files.

    Sharing paths with the full run meant `--quick` silently replaced a
    quotable figure with an under-powered one, and nothing warned.  """
    return ('val_duplinskiy_birefringence--quick' if quick
            else 'val_duplinskiy_birefringence')


def _write_csv(rows, quick=False):
    path = os.path.join(OUT_DIR, _stem(quick) + '.csv')
    with open(path, 'w') as fh:
        fh.write("# Duplinskiy birefringence swing, "
                 "validate_duplinskiy_birefringence.py\n")
        fh.write(f"# seed={SEED} target_sifted={TARGET_SIFTED}\n")
        fh.write("distance_km,pulses,qber_compensated,sigma_compensated,"
                 "sifted_compensated,qber_uncompensated,sigma_uncompensated,"
                 "sifted_uncompensated\n")
        for (km, n, qc, sc, nc, qu, su, nu, _, _) in rows:
            fh.write(f"{km},{n},{qc:.6f},{sc:.6f},{nc},"
                     f"{qu:.6f},{su:.6f},{nu}\n")
    print(f"\n  CSV: {path}")


def _figure(rows, quick=False):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("  matplotlib missing, skipping the figure")
        return

    km = [r[0] for r in rows]
    qc = [100 * r[2] for r in rows]
    sc = [100 * r[3] for r in rows]
    qu = [100 * r[5] for r in rows]
    su = [100 * r[6] for r in rows]

    fig, ax = plt.subplots(figsize=(8.6, 5.2))
    x = np.arange(len(km))
    w = 0.34
    ax.bar(x - w / 2, qc, w, yerr=sc, capsize=4, color='tab:blue',
           edgecolor='0.25', linewidth=0.7, label='compensation on')
    ax.bar(x + w / 2, qu, w, yerr=su, capsize=4, color='tab:red',
           edgecolor='0.25', linewidth=0.7, label='compensation off')
    for i in x:
        ax.text(i - w / 2, qc[i] + sc[i] + 1.2, f'{qc[i]:.2f}', ha='center',
                fontsize=8)
        ax.text(i + w / 2, qu[i] + su[i] + 1.2, f'{qu[i]:.2f}', ha='center',
                fontsize=8)

    # No Gobby line here.  Its null is bit-identity rather than a value, so
    # drawing it as a horizontal QBER would misrepresent the claim -- and
    # the two chains do not share an error budget, so the heights would not
    # be comparable even when both are measured.  See
    # validate_gobby_impairments.py.
    ax.set_xticks(x)
    ax.set_xticklabels([f'{k} km' for k in km])
    ax.set_xlabel('fibre length')
    ax.set_ylabel('QBER (%)')
    ax.set_ylim(0, 100)
    ax.set_title('QBER with and without polarisation compensation\n'
                 'Duplinskiy chain, birefringence only', fontsize=11)
    ax.grid(True, axis='y', alpha=0.3)
    ax.legend(fontsize=9, loc='upper left')

    fig.tight_layout()
    png = os.path.join(OUT_DIR, _stem(quick) + '.png')
    fig.savefig(png, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  figure: {png}")


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--quick', action='store_true',
                    help='a twelfth of the pulses, for a smoke run')
    a = ap.parse_args()
    sys.exit(run(quick=a.quick))
