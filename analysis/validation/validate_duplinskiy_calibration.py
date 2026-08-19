"""The paper's Fig. 6: the eight-column calibration histogram.

Duplinskiy et al., Opt. Express 25(23), 28886 (2017), sec. 5:

    "Alice applies four different voltages to the phase modulator,
     corresponding to phase shifts of 0, pi/2, pi and 3pi/2.  Bob applies
     two voltages, corresponding to 0 and pi/2. ... We begin the
     calibration by amassing the statistics of detector clicks
     corresponding to each pair of voltages, receiving a histogram of
     eight columns (Fig. 6), as shown in the Table 1."

Those eight (Alice phase, Bob phase) pairs are exactly this chain's
8-outcome RESPONSE table, one for one, and the paper's Table 1 lists them
with the same basis and bit assignment.  So Fig. 6 is directly
reproducible, and it tests something QBER does not: the **calibration
criterion** the tuning algorithm actually runs on.

The criterion
-------------
The paper states which columns must be indistinguishable once the
polarisation controllers are set correctly:

    "the pulses that experienced the pairs of shifts (0, 0) and
     (pi/2, pi/2) should not be distinguishable.  This means that these
     two pulses will have the same statistics of clicks."

and then three more pairs, at relative shifts of pi/2, pi and 3pi/2.
Four pairs in total, and every one of them must match.

A mapping is needed, and it is not a fudge
------------------------------------------
This chain's outcome depends on the phase **sum** phi_A + phi_B, not the
difference, so the paper's pairs look broken until Bob's basis labels are
mapped as established: chain-X(pi/2) is the paper's
"linear"(0), chain-C(0) the paper's "circular"(pi/2).  With that mapping
every pair holds exactly.

That makes this an independent re-confirmation of register entry A4 --
the third route to it, after the state table and the Jones
algebra.  Here it comes from the calibration criterion instead.

References
----------
[1] Duplinskiy, Ustimchik, Kanapin, Kurochkin & Kurochkin, Opt. Express
    25(23), 28886-28897 (2017).  Fig. 6, Table 1, sec. 5 (Tuning).
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from src.channel import optics, FiberRealization
from src.channel.phase_modulator import PhaseModulator
from src.detectors.spad import spad

OUT_DIR = os.path.join(os.path.dirname(__file__), '..', 'val_duplinskiy')

# This script has no --quick flag; it takes --pulses, and 20000 is the
# full run.  Anything smaller is a reduced run and writes to its own file
# rather than replacing the committed figure.
FULL_PULSES = 20000


def _stem(reduced):
    return ('val_duplinskiy_calibration--quick' if reduced
            else 'val_duplinskiy_calibration')

SEED = 42
DISTANCE = 50

# SEED reaches `FiberRealization`, which seeds its own stream, but the
# detectors draw from the GLOBAL numpy state: `spad.detect` calls
# `np.random.random` for the dark count, the click and the afterpulse, and
# `np.random.exponential` for the afterpulse delay.  Without this line the
# click histogram -- which is the entire figure -- came out different on
# every run, so the committed artifact could not be reproduced and every
# harness pass left the repository dirty for no reason.
#
# Same fix and same placement as `validate_cwlaser.py`.
np.random.seed(SEED)

# The paper's Table 1, verbatim.  Pulse numbers 1-8, with C = circular and
# L = linear in its notation; this chain calls the linear basis X.
TABLE_1 = [
    # (Alice bit, Alice basis, Bob basis, sifted bit or None)
    (0, 'C', 'X', None),
    (1, 'X', 'X', 1),
    (1, 'C', 'X', None),
    (0, 'X', 'X', 0),
    (0, 'C', 'C', 0),
    (1, 'X', 'C', None),
    (1, 'C', 'C', 1),
    (0, 'X', 'C', None),
]

PHASE_NAME = {0: '0', 1: r'$\pi/2$', 2: r'$\pi$', 3: r'$3\pi/2$'}

# Alice's four voltages as phase indices (units of pi/2), and which
# (basis, bit) each corresponds to in this chain's encoding.
ALICE_PHASE = {(0, 'X'): 0, (0, 'C'): 1, (1, 'X'): 2, (1, 'C'): 3}

# Bob's two, and the mapping to the paper's labelling.
# The chain's X basis drives PM2 at pi/2, which is the paper's "linear" at
# 0; the chain's C basis drives it at 0, the paper's "circular" at pi/2.
BOB_PHASE = {'X': 1, 'C': 0}
BOB_PHASE_PAPER = {'X': 0, 'C': 1}

# The paper's four indistinguishability pairs, in ITS (Alice, Bob) phase
# indices.  Equivalent to "outcomes depend only on phi_A - phi_B".
PAPER_PAIRS = [((0, 0), (1, 1)), ((1, 0), (2, 1)),
               ((2, 0), (3, 1)), ((3, 0), (0, 1))]


def build_chain():
    pm_a = PhaseModulator(crystal_cut='X', modulation='DC')
    pm_b = PhaseModulator(crystal_cut='X', modulation='DC')
    fibre = FiberRealization(L_m=DISTANCE * 1000, temperature=25,
                             bend_radius=None, attenuation_factor=0.2,
                             cd=False, pmd=False, model='auto', seed=SEED)
    J = fibre.birefringence_matrix()
    return pm_a, pm_b, fibre, (None if J is None else J.conj().T)


def port_powers(a_bit, a_basis, b_basis, chain, mu=0.1, gate=20e-9):
    """The two analyser port powers for one (Alice, Bob) setting."""
    pm_a, pm_b, fibre, U = chain
    Vpi = pm_a.Vpi
    ppp = mu * (6.626e-34 * 3e8 / 1550e-9) / gate

    v_a = (Vpi / 2 if a_bit == 0 else 3 * Vpi / 2) if a_basis == 'C' else \
          (0 if a_bit == 0 else Vpi)
    E = np.sqrt(ppp / 2.0) * np.ones((1, 2), dtype=complex)
    E = pm_a.modulate(E_field=E, V=v_a)
    E = fibre.apply(E, dt=1e-7)
    if U is not None:
        E = np.transpose(U @ np.transpose(E))
    E = optics.voa(E, 2.0)
    E = pm_b.modulate(E_field=E, V=(0 if b_basis == 'C' else Vpi / 2))
    Ex, Ey = optics.circular_analyser(E)
    return float(np.mean(np.abs(Ex) ** 2)), float(np.mean(np.abs(Ey) ** 2))


def click_counts(n_pulses, chain):
    """Detector click counts for each of the eight settings.

    This is the histogram the paper's tuning algorithm accumulates: the
    same eight columns, on the same two detectors.
    """
    d1 = spad(wavelength=1550e-9, quantum_efficiency=0.10, dead_time=13e-6,
              dark_count_rate=15.0, afterpulse_prob=0.05, gate_width=20e-9)
    d2 = spad(wavelength=1550e-9, quantum_efficiency=0.10, dead_time=13e-6,
              dark_count_rate=15.0, afterpulse_prob=0.05, gate_width=20e-9)
    dt = 1.0 / 10e6

    counts = {}
    t = 0.0
    for k, (bit, ab, bb, _) in enumerate(TABLE_1):
        px, py = port_powers(bit, ab, bb, chain)
        c1 = c2 = 0
        for i in range(n_pulses):
            c1 += d1.detect(px, t)
            c2 += d2.detect(py, t)
            t += dt
        counts[k] = (c1, c2, px, py)
    return counts


def check_pairs(chain, failures):
    """The paper's four indistinguishability conditions, on port fractions."""
    print("\n  the paper's four indistinguishability pairs (sec. 5)")
    print("    a pair is indistinguishable when the two settings give the")
    print("    same port statistics.  Bob's labels are mapped per .")

    # Index the chain's eight settings by the PAPER's (Alice, Bob) phases.
    by_paper = {}
    for bit, ab, bb, _ in TABLE_1:
        pa = ALICE_PHASE[(bit, ab)]
        pb_paper = BOB_PHASE_PAPER[bb]
        px, py = port_powers(bit, ab, bb, chain)
        by_paper[(pa, pb_paper)] = px / (px + py)

    print(f"\n    {'pair (Alice, Bob) phases':<34}{'port-1 fraction':>22}   match")
    worst = 0.0
    for (p, q) in PAPER_PAIRS:
        f1, f2 = by_paper[p], by_paper[q]
        d = abs(f1 - f2)
        worst = max(worst, d)
        lab = (f"({PHASE_NAME[p[0]]},{PHASE_NAME[p[1]]}) <-> "
               f"({PHASE_NAME[q[0]]},{PHASE_NAME[q[1]]})")
        lab = lab.replace('$', '').replace('\\pi', 'pi')
        print(f"    {lab:<34}{f1:10.6f} / {f2:9.6f}   "
              f"{'yes' if d < 1e-9 else 'NO'}")
    print(f"\n    worst mismatch across the four pairs: {worst:.2e}")
    if worst > 1e-9:
        failures.append(f"the paper's indistinguishability pairs do not hold "
                        f"(worst {worst:.2e}); either the readout or the "
                        "basis mapping is wrong")
    else:
        print("    -> all four hold. Independent re-confirmation of A4, from")
        print("       the calibration criterion rather than the state table.")
    return by_paper


def run(n_pulses=20000):
    os.makedirs(OUT_DIR, exist_ok=True)
    failures = []
    print("=" * 72)
    print("Duplinskiy Fig. 6: the eight-column calibration histogram")
    print("=" * 72)
    print(f"  {DISTANCE} km, seed {SEED}, {n_pulses:,} pulses per column")

    chain = build_chain()
    counts = click_counts(n_pulses, chain)

    print(f"\n  {'#':>2}  {'Alice':>7} {'basis':>6} {'Bob':>5}  "
          f"{'sifted':>6}   {'D1':>7} {'D2':>7}")
    for k, (bit, ab, bb, sift) in enumerate(TABLE_1):
        c1, c2, _, _ = counts[k]
        print(f"  {k + 1:2d}  {bit:>7} {ab:>6} {bb:>5}  "
              f"{'-' if sift is None else sift:>6}   {c1:7d} {c2:7d}")
    print("    (Table 1's assignment; '-' = bases differ, sifted out)")

    # Assert Table 1's sifted-bit column rather than reading it off.  For a
    # matched basis the majority detector gives the bit; for a mismatched
    # one the two must be comparable, which is what makes the pulse
    # discardable rather than wrong.
    print("\n  does the histogram reproduce Table 1's sifted bits?")
    for k, (bit, ab, bb, sift) in enumerate(TABLE_1):
        c1, c2, _, _ = counts[k]
        if sift is not None:
            got = 0 if c1 > c2 else 1
            ok = got == sift
            print(f"    pulse {k + 1}: bases match, D1={c1} D2={c2} "
                  f"-> bit {got}, Table 1 says {sift}   "
                  f"{'ok' if ok else 'MISMATCH'}")
            if not ok:
                failures.append(f"pulse {k + 1} decodes to {got}, Table 1 "
                                f"says {sift}")
        else:
            tot = c1 + c2
            frac = c1 / tot if tot else 0.5
            # 5 sigma on a fair coin at this count.
            tol = 5 * (0.25 / max(tot, 1)) ** 0.5
            ok = abs(frac - 0.5) <= tol
            print(f"    pulse {k + 1}: bases differ, D1={c1} D2={c2} "
                  f"-> {100 * frac:.1f} % on D1 (expect 50)   "
                  f"{'ok' if ok else 'SKEWED'}")
            if not ok:
                failures.append(f"pulse {k + 1} has mismatched bases but "
                                f"splits {100 * frac:.1f}/{100 * (1 - frac):.1f}; "
                                "a discarded pulse should carry no bit "
                                "information")

    by_paper = check_pairs(chain, failures)
    _figure(counts, by_paper, n_pulses < FULL_PULSES)

    print()
    if failures:
        print("[FAIL]")
        for f in failures:
            print(f"  - {f}")
        return 1
    print("[PASS] the eight columns reproduce Table 1 and satisfy all four")
    print("       of the paper's indistinguishability conditions")
    return 0


def _figure(counts, by_paper, reduced=False):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("  matplotlib missing, skipping the figure")
        return

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6), sharey=True)
    xs = np.arange(len(TABLE_1))
    labels = [f"{k + 1}\n{ab}{bit}\n{bb}"
              for k, (bit, ab, bb, _) in enumerate(TABLE_1)]
    matched = [ab == bb for _, ab, bb, _ in TABLE_1]

    for d, (axis, title) in enumerate(zip(ax, ('detector 1', 'detector 2'))):
        vals = [counts[k][d] for k in xs]
        colours = ['tab:blue' if m else '0.75' for m in matched]
        axis.bar(xs, vals, color=colours, edgecolor='0.25', linewidth=0.7)
        for k in xs:
            axis.text(k, vals[k] + max(vals) * 0.02, str(vals[k]),
                      ha='center', fontsize=8)
        axis.set_xticks(xs)
        axis.set_xticklabels(labels, fontsize=8)
        axis.set_xlabel('pulse number / Alice basis+bit / Bob basis')
        axis.set_title(title, fontsize=10)
        axis.grid(True, axis='y', alpha=0.3)
    ax[0].set_ylabel('clicks')

    handles = [plt.Rectangle((0, 0), 1, 1, color='tab:blue'),
               plt.Rectangle((0, 0), 1, 1, color='0.75')]
    ax[0].legend(handles, ['bases match (sifted)', 'bases differ (discarded)'],
                 fontsize=8, loc='upper right')

    fig.suptitle('Duplinskiy chain: click statistics for the eight '
                 'calibration settings\n'
                 f'{DISTANCE} km, reproducing the paper\'s Fig. 6 and Table 1',
                 fontsize=11)
    fig.tight_layout()
    png = os.path.join(OUT_DIR, _stem(reduced) + '.png')
    fig.savefig(png, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"\n  figure: {png}")


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--pulses', type=int, default=FULL_PULSES,
                    help='pulses per column (default 20000)')
    a = ap.parse_args()
    sys.exit(run(n_pulses=a.pulses))
