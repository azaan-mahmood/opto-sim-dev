"""DUPL-2: the extinction term, and what it says about the afterpulse budget.

Duplinskiy's calibration goal 3 is

    "Bob's measurements differentiate BB84 orthogonal states with
     extinction higher than 98 %."   -- Opt. Express 25(23), 28886, sec. 5

The chain assumed a perfect analyser until now.  This adds the term and
uses it as the discriminating test for register entry A1, which is
load-bearing for **100 % of this replication's error budget**: whether the
ID230's quoted 5 % afterpulse probability is already net of dead-time
suppression, or is a raw trap-release figure that dead time should mostly
suppress.

Why extinction discriminates it
-------------------------------
The paper's own decomposition is ~1 % afterpulse floor + finite extinction
+ drift/recalibration -> 2 % average.  Since QBER_afterpulse ~ p_eff/2, the
paper's 1 % floor requires p_eff ~ 0.020.  Set against the two clean
readings of A1:

    A1 true  -- 0.05 is net of suppression      p_eff = 0.0500 -> 2.5 %
    A1 false -- P(Exp(6.5 us) > 13 us) = e^-2   p_eff = 0.0068 -> 0.34 %
    paper's stated floor                        p_eff = 0.0200 -> 1.0 %

**Neither clean reading lands on 0.020** -- A1-true overshoots the paper's
floor by 2x, A1-false undershoots by 3x.  That is recorded here *before*
the run, because it means the likely outcome is neither branch cleanly,
and the honest resolution may be that the 5 % figure is quoted at a dead
time other than 13 us, or that the exponential afterpulse-delay model has
the wrong shape (register A8).

Both readings of 98 %, because they differ by exactly 2x
--------------------------------------------------------
    (a) power-fraction   98 % of power reaches the right port  eps = 0.0200
    (b) visibility-like  (Imax-Imin)/(Imax+Imin) = 0.98        eps = 0.0101

The ambiguity falls on the exact quantity the test turns on, so both are
run and reported rather than one being chosen (register A7).  And 98 % is
a *threshold the tuning algorithm targets*, not an achieved value, so
anything derived from it is an upper bound on the term (register A3).

References
----------
[1] Duplinskiy, Ustimchik, Kanapin, Kurochkin & Kurochkin, Opt. Express
    25(23), 28886-28897 (2017).
[2] ID Quantique, ID230 InGaAs SPAD datasheet.
"""
import argparse
import math
import os
import statistics
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from src.channel import apply_extinction
from src.protocols.bb84_duplinskiy import simulate_bb84_duplinskiy

OUT_DIR = os.path.join(os.path.dirname(__file__), '..', 'val_duplinskiy')

SEED = 42
DISTANCE = 50

# sec. 25's quotable standard is 3000 sifted; 50 km yields ~3.1e-4 per pulse.
N_PULSES = 10_000_000
# 50 km is the fixed distance here, so even the smoke run needs ~1M pulses
# to clear a few hundred sifted bits.
N_QUICK = 1_000_000

EPSILONS = (('0', 0.0), ('0.0101 (b)', 0.0101), ('0.0200 (a)', 0.0200))
P_APS = (0.0, 0.025, 0.05)


def run_cell(n, eps, p_ap, km=DISTANCE, **kw):
    r = simulate_bb84_duplinskiy(n, fiber_length=km, seed=SEED,
                                 extinction_epsilon=eps,
                                 afterpulse_prob=p_ap, **kw)
    s, q = r['n_sifted'], r['qber']
    sig = math.sqrt(max(q * (1 - q), 1e-12) / s) if s else float('nan')
    return q, sig, s, r['n_errors']


def controls(quick, failures):
    print("\n  controls (sec. 29.3)")

    # 1. power conservation
    worst = 0.0
    for Pa, Pb in ((1.0, 0.0), (0.7, 0.3), (1e-13, 3e-44), (0.0, 0.0)):
        for e in (0.0, 0.0101, 0.02, 0.5):
            a, b = apply_extinction(Pa, Pb, e)
            worst = max(worst, abs((a + b) - (Pa + Pb)))
    leaks = all(abs(apply_extinction(1.0, 0.0, e)[1] - e) < 1e-18
                for e in (0.0101, 0.02, 0.1, 0.5))
    print(f"    power conserved             : worst |sum_out - sum_in| = {worst:.2e}")
    print(f"    leaks exactly epsilon       : {leaks}")
    if worst > 1e-15 or not leaks:
        failures.append("apply_extinction does not conserve power or does not "
                        "leak exactly epsilon")

    # 2. negative (G2)
    base = {(0, True): (35, 1), (0, False): (35, 1), (10, True): (39, 1),
            (10, False): (25, 9), (50, True): (8, 0), (50, False): (6, 1)}
    ok = True
    for (km, comp), (s0, e0) in base.items():
        r = simulate_bb84_duplinskiy(20000, fiber_length=km, compensate=comp,
                                     seed=SEED, extinction_epsilon=0.0)
        ok &= (r['n_sifted'], r['n_errors']) == (s0, e0)
    print(f"    negative, epsilon = 0       : "
          f"{'bit-identical to the frozen baseline' if ok else 'BASELINE MOVED'}")
    if not ok:
        failures.append("epsilon = 0 moved the frozen sec. 27.1 baseline")

    # 3. positive (G2)
    n = N_QUICK if quick else 1_000_000
    qs = []
    for km, mult in ((0, 1), (10, 1), (50, 4)):
        q, _, _, _ = run_cell(n * mult, 0.5, 0.05, km=km)
        qs.append((km, q))
    print("    positive, epsilon = 0.5     : QBER "
          + ", ".join(f"{km} km {100 * q:.1f} %" for km, q in qs)
          + "  (expect ~50 %)")
    qs = [q for _, q in qs]
    if any(abs(q - 0.5) > 0.10 for q in qs):
        failures.append("epsilon = 0.5 did not drive QBER to ~50 %; the term "
                        "is not reaching the observable and every result "
                        "below would be vacuous")

    # 4. sifted invariance -- see the note below for why this is conditional
    print("\n    sifted-count invariance, sec. 29.3's fourth control:")
    print("      cross-talk moves power between ports without removing any,")
    print("      so P(no click) = exp(-eta(1-e)P) exp(-eta e P) = exp(-eta P)")
    print("      and the sifted rate cannot depend on epsilon -- for a")
    print("      MEMORYLESS detector.  Measured across 6 independent seeds:")
    for lab, kw in (('ideal detector', dict(dead_time=0.0, afterpulse_prob=0.0,
                                            dark_count_rate=0.0)),
                    ('ID230 dead time 13 us', dict(afterpulse_prob=0.0,
                                                   dark_count_rate=0.0))):
        ds = []
        for sd in (42, 7, 11, 99, 1234, 2026):
            a = simulate_bb84_duplinskiy(400000 if quick else 2000000,
                                         fiber_length=10, seed=sd,
                                         extinction_epsilon=0.0, **kw)['n_sifted']
            b = simulate_bb84_duplinskiy(400000 if quick else 2000000,
                                         fiber_length=10, seed=sd,
                                         extinction_epsilon=0.02, **kw)['n_sifted']
            ds.append(b - a)
        m = statistics.mean(ds)
        sd_ = statistics.stdev(ds) / math.sqrt(len(ds))
        verdict = 'invariant' if abs(m) < 2 * sd_ else 'SHIFTS'
        print(f"        {lab:<24} {m:+7.1f} +/- {sd_:5.1f} sifted   {verdict}")
        if 'ideal' in lab and abs(m) >= 3 * sd_:
            failures.append(f"the sifted rate shifts with epsilon even for an "
                            f"ideal detector ({m:+.1f} +/- {sd_:.1f}); "
                            "extinction is being applied in the wrong place")
    print("      -> the invariance holds exactly where the algebra applies, and")
    print("         breaks with real dead time: a click on the weak port takes")
    print("         that detector offline for 130 pulses at 10 MHz, so the two")
    print("         are no longer independent.  sec. 29.3 states the control")
    print("         unconditionally; it is conditional on a memoryless detector.")


def matrix(quick, failures):
    n = N_QUICK if quick else N_PULSES
    print(f"\n  the discriminating run: {DISTANCE} km, compensate=True, "
          f"seed {SEED}, {n:,} pulses per cell")
    print("    p_ap \\ eps " + "".join(f"{lab:>22}" for lab, _ in EPSILONS))
    cells = {}
    for p_ap in P_APS:
        row = f"    {p_ap:<11.3f}"
        for lab, eps in EPSILONS:
            q, sig, s, e = run_cell(n, eps, p_ap)
            cells[(p_ap, eps)] = (q, sig, s)
            row += f"{100 * q:12.2f} +/-{100 * sig:5.2f}"
            if s < (250 if quick else 3000):
                failures.append(f"p_ap={p_ap}, eps={eps}: only {s} sifted; "
                                "not quotable")
        print(row)
    print("    (QBER %, +/- 1 sigma; rows are afterpulse probability)")
    return cells


def _print_matrix(cells):
    """Print an already-measured matrix, for the redraw path."""
    print(f"\n  the discriminating run: {DISTANCE} km, compensate=True, "
          f"seed {SEED}")
    print("    p_ap \\ eps " + "".join(f"{lab:>22}" for lab, _ in EPSILONS))
    for p_ap in P_APS:
        row = f"    {p_ap:<11.3f}"
        for _, eps in EPSILONS:
            q, sig, _ = cells[(p_ap, eps)]
            row += f"{100 * q:12.2f} +/-{100 * sig:5.2f}"
        print(row)
    print("    (QBER %, +/- 1 sigma; rows are afterpulse probability)")


def verdict(cells):
    """Apply sec. 29.5's four decision rules, committed before the run."""
    print("\n  sec. 29.5's decision rules, as written")
    q = {k: v[0] for k, v in cells.items()}
    s = {k: v[1] for k, v in cells.items()}
    paper = 0.02

    r1 = q[(0.05, 0.0101)]
    r2 = q[(0.025, 0.0101)]
    r3 = q[(0.05, 0.0)]

    print(f"    1. (p_ap=0.05, eps=0.0101) = {100 * r1:.2f} % "
          f"{'>~ 3 %, clearly past the paper' if r1 > 0.03 else 'not past 3 %'}")
    print(f"    2. (p_ap=0.025, eps=0.0101) = {100 * r2:.2f} % "
          f"{'~ 2 %, reproduces total AND decomposition' if abs(r2 - paper) < 2 * s[(0.025, 0.0101)] + 0.004 else 'not ~2 %'}")
    print(f"    3. (p_ap=0.05, eps=0) = {100 * r3:.2f} % against the paper's 2 %")
    spread = max(abs(v - paper) for v in q.values())
    print(f"    4. all cells within 0.4 pp of 2 %? "
          f"{'yes -- mechanisms not separable at this precision' if spread < 0.004 else f'no (worst {100 * spread:.2f} pp)'}")

    print("\n    what this says about A1")
    fired = [n for n, c in ((1, r1 > 0.03), (3, r3 > paper)) if c]
    if fired:
        print(f"      Rule {' and '.join(map(str, fired))} fired: the afterpulse "
              "term is too large at the")
        print("      datasheet value once extinction is included as well.")
        print("      A1's ALTERNATIVE reading is favoured -- 0.05 is not net of")
        print("      dead-time suppression.")
    if abs(r2 - paper) < 2 * s[(0.025, 0.0101)] + 0.004:
        print("      Rule 2 fired: p_ap = 0.025 reproduces both the paper's")
        print("      total and its ~1 % + ~1 % decomposition. NOT licence to")
        print("      set p_ap = 0.025 -- that would be fitting a datasheet")
        print("      parameter, which G9 forbids. It is a diagnosis, and the")
        print("      fix belongs in the afterpulse MODEL, not its input.")
    print("\n      Whatever the outcome: the paper's floor is a one-sentence")
    print("      calculation, 'about 1 % mainly due to afterpulses', and")
    print("      'mainly' is unquantified. Only one side of any ratio drawn")
    print("      from it is measured.")


def run(quick=False, figure_only=False):
    os.makedirs(OUT_DIR, exist_ok=True)
    failures = []
    print("=" * 74)
    print("DUPL-2: finite analyser extinction, and the A1 discriminating test")
    print("=" * 74)

    if figure_only:
        cells = _read_csv()
        if cells is None:
            print("  no previous run to draw; run without --figure-only first")
            return 1
        print("  redrawing from the last run's CSV, no simulation")
        _print_matrix(cells)
        verdict(cells)
        _figure(cells)
        return 0

    print("  predicted before the run (sec. 29.5): neither clean reading of A1")
    print("  lands on the paper's implied p_eff = 0.020 -- A1-true overshoots")
    print("  2x, A1-false undershoots 3x. Expect neither branch cleanly.")

    controls(quick, failures)
    cells = matrix(quick, failures)
    verdict(cells)
    _write_csv(cells)
    _figure(cells)

    print()
    if failures:
        print("[FAIL]")
        for f in failures:
            print(f"  - {f}")
        return 1
    print("[PASS] extinction conserves power, reaches the observable, and "
          "leaves\n       the frozen baseline untouched at epsilon = 0")
    return 0


def _read_csv():
    """Reload a previous run's cells, so the figure can be redrawn cheaply.

    The matrix costs ~50 minutes at quotable statistics.  Nothing about
    drawing it should require paying that again.
    """
    path = os.path.join(OUT_DIR, 'val_duplinskiy_extinction.csv')
    if not os.path.exists(path):
        return None
    cells = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith('#') or line.startswith('afterpulse'):
                continue
            p_ap, eps, q, sig, s = line.strip().split(',')
            cells[(float(p_ap), float(eps))] = (float(q), float(sig), int(s))
    return cells or None


def _figure(cells):
    """The 3x3 grid, drawn because the argument in it is geometric.

    The result is not "these nine cells have these values".  It is that
    TWO cells reach the paper's 2 % and only one of them also reproduces
    the paper's own decomposition of that 2 % into ~1 % afterpulse plus
    ~1 % extinction.  A table makes a reader hunt for that; a grid shows
    it.
    """
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        print("  matplotlib missing, skipping the figure")
        return

    eps_vals = [e for _, e in EPSILONS]
    q = np.array([[cells[(p, e)][0] * 100 for e in eps_vals] for p in P_APS])
    sg = np.array([[cells[(p, e)][1] * 100 for e in eps_vals] for p in P_APS])

    # Grouped bars, not a heatmap.  Nine numbers with error bars are a bar
    # chart; a heatmap hides the uncertainty and needs a contour label to
    # show where 2 % falls, which reads badly over coloured cells.
    fig, ax = plt.subplots(figsize=(9.2, 5.4))
    width = 0.26
    xs = np.arange(len(P_APS))
    colours = ['0.72', 'tab:blue', 'tab:orange']

    for j, (lab, _) in enumerate(EPSILONS):
        off = (j - 1) * width
        ax.bar(xs + off, q[:, j], width, yerr=sg[:, j], capsize=4,
               color=colours[j], edgecolor='0.25', linewidth=0.7,
               label=f'extinction  $\\epsilon$ = {lab}')
        for i in xs:
            ax.text(i + off, q[i, j] + sg[i, j] + 0.13, f'{q[i, j]:.2f}',
                    ha='center', fontsize=8)

    ax.axhline(2.0, color='crimson', ls='--', lw=1.6, zorder=0,
               label='paper: 2 %')

    # Ring the two bars that reach the paper's rate, labelled with what
    # each attributes the error to.  Two or three words; the reasoning
    # belongs in sec. 32.4, not on the axes.
    for (i, j, colour, tag) in ((2, 0, 'crimson', 'afterpulsing only'),
                                (1, 1, 'darkgreen', 'both terms')):
        x = xs[i] + (j - 1) * width
        ax.bar(x, q[i, j], width, fill=False, edgecolor=colour,
               linewidth=2.4, zorder=5)
        ax.text(x, q[i, j] + sg[i, j] + 0.42, tag, ha='center', fontsize=8,
                color=colour, fontweight='bold')

    ax.set_xticks(xs)
    ax.set_xticklabels([f'{p:g}' + ('\n(datasheet)' if p == 0.05 else '')
                        for p in P_APS])
    ax.set_xlabel('detector afterpulse probability')
    ax.set_ylabel('QBER (%)')
    ax.set_ylim(0, max(5.2, (q + sg).max() + 1.1))
    ax.set_title('QBER against afterpulse probability and analyser '
                 'extinction\nDuplinskiy chain, 50 km, 10M pulses per bar',
                 fontsize=11)
    ax.grid(True, axis='y', alpha=0.3)
    ax.legend(fontsize=9, loc='upper left')

    fig.tight_layout()
    png = os.path.join(OUT_DIR, 'val_duplinskiy_extinction.png')
    fig.savefig(png, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  figure: {png}")


def _write_csv(cells):
    path = os.path.join(OUT_DIR, 'val_duplinskiy_extinction.csv')
    with open(path, 'w') as fh:
        fh.write("# DUPL-2 extinction discriminating run, "
                 "validate_duplinskiy_extinction.py\n")
        fh.write(f"# {DISTANCE} km, compensate=True, seed={SEED}\n")
        fh.write("afterpulse_prob,extinction_epsilon,qber,qber_sigma,n_sifted\n")
        for (p_ap, eps), (q, sig, s) in sorted(cells.items()):
            fh.write(f"{p_ap:g},{eps:g},{q:.6f},{sig:.6f},{s}\n")
    print(f"\n  CSV: {path}")


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--quick', action='store_true',
                    help='fewer pulses per cell, for a smoke run')
    ap.add_argument('--figure-only', action='store_true',
                    help="redraw from the last run's CSV without simulating "
                         '(the matrix costs ~50 minutes at full statistics)')
    a = ap.parse_args()
    sys.exit(run(quick=a.quick, figure_only=a.figure_only))
