"""The urban link: Duplinskiy's Fig. 7, and the key rates nobody checked.

The paper's second experiment is a deployed 30 km urban line, and it is
fully specified -- more fully than the lab run, in fact, because it comes
with two ABSOLUTE key rates:

    "an environment experiment has been carried out for a 30 km urban line
     with high losses (13 dB) ... The residual illumination has been about
     200Hz. To reduce its impact, system has been upgraded to a 5 ns
     detection window ... the average photons per pulse has been lowered to
     0.02. Under these conditions, 106 bit/s sifted key rate has been
     obtained (265 bit/s raw key during key distribution regime). Figure 7
     illustrates QBER statistics vs. time during urban tests for 20 hours,
     the average value being 5.5%. ... about 20% has been spent for tuning.
     For this experiment the lower bound for the QBER is about 4,5% mostly
     caused by the noise from the telecommunicational channel, together
     with the low intensity of the key pulses."   -- Opt. Express 25(23), sec. 6

Three registered assumptions
----------------------------
1.  **alpha = 13/30 = 0.4333 dB/km is a LUMPED term** (guardrail G3), the
    same treatment `apply_extinction` gets.  A deployed urban span is
    ~0.2 dB/km of fibre plus splices, patch panels and connectors; the
    model carries one attenuation coefficient, so 13 dB over 30 km goes in
    as an effective per-km figure.  The alternative -- 65 km at
    0.2 dB/km -- also gives 13 dB, but then the birefringence realisation
    is that of 65 km rather than the 30 km the paper states.

2.  **215 Hz is not a dark count rate.**  15 Hz of it is the ID230's, and
    200 Hz is stray light through Bob's WDM filter at 1554.94 nm.  Same
    effect on the observable, uncorrelated counts inside the gate;
    different mechanism.  This is the naming defect sec. 19 raised against
    Gobby's `DCR`, recorded rather than renamed here.

3.  **13 dB is the LINE, not the link.**  Bob's ~2 dB is stated separately
    in the paper's sec. 4, so the total is 15 dB.  The alternative reading
    is not argued -- it is measured below and rejected by the paper's own
    key rate.

Reproduced, and reproduced-in-distribution
------------------------------------------
Fig. 7's *configuration* needs nothing the paper withholds.  What it does
withhold is a drift rate, and that governs the WIDTH of the scatter band,
not whether the run can be built.  So:

  - the mean and the physics are tested here;
  - the band width is REPORTED, never matched.  Block size is a
    presentation choice, so tuning it until the bands agree would be
    fitting (G9), and the drift rate is not in the paper and is not
    invented.

A 20 hour trace at 5 MHz is 3.6e11 pulses, about 11 days of compute.  It is
not run.  The paper's own figure shows no trend across its 20 hours, so a
stationary scatter band is its entire content, and the distribution is what
compares.

The claim this trace tests, stated before the run
-------------------------------------------------
Blocks are not obviously binomial.  Afterpulsing correlates consecutive
clicks, so the block-to-block variance may exceed the binomial prediction.

    Measure Var(block QBER) / [q(1-q)/N].  If it is 1, the scatter is pure
    counting statistics.  If it exceeds 1, afterpulse correlation inflates
    it, and part of Fig. 7's width needs no environmental drift at all.

At 120 blocks a ratio beyond about 1.26 is resolvable.  This is the same
afterpulse model registers A1 and A8 turn on, reached from a third
direction after sec. 32's discriminating grid.

References
----------
[1] Duplinskiy, Ustimchik, Kanapin, Kurochkin & Kurochkin, Opt. Express
    25(23), 28886-28897 (2017).
[2] ID Quantique, ID230 InGaAs SPAD datasheet.
"""
import argparse
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from src.protocols.bb84_duplinskiy import simulate_bb84_duplinskiy

OUT_DIR = os.path.join(os.path.dirname(__file__), '..', 'val_duplinskiy')

SEED = 42

# --- the paper's urban configuration, sec. 6 -----------------------------
KM = 30                    # "a 30 km urban line"
ALPHA = 13.0 / 30.0        # "high losses (13 dB)" -- lumped, assumption 1
BOB_LOSS = 2.0             # "about 2dB losses for the optical scheme of Bob's device", sec. 4
MU = 0.02                  # "the average photons per pulse has been lowered to 0.02"
GATE = 5e-9                # "upgraded to a 5 ns detection window"
REP = 5e6                  # "effective key generation frequency was reduced to 5 MHz"
STRAY = 215.0              # 200 Hz residual illumination + the ID230's 15 Hz
DUTY = 0.80                # "about 20% has been spent for tuning"

# the paper's own numbers, for comparison only
PAPER_QBER_MEAN = 0.055
PAPER_QBER_LOWER = 0.045
PAPER_BAND = (0.03, 0.09)
PAPER_RAW = 265.0
PAPER_SIFTED = 106.0

# Measured sifted yield at this configuration, used only to budget runs.
YIELD = 3.34e-5
TARGET_SIFTED = 3000       # sec. 25's quotable standard

N_BLOCKS = 120
BLOCK_SIFTED = 300         # -> sigma ~1.26 pp per block


def _base(**kw):
    cfg = dict(fiber_length=KM, alpha_dB=ALPHA, bob_loss_dB=BOB_LOSS, mu=MU,
               gate_width=GATE, rep_rate=REP, dark_count_rate=STRAY,
               seed=SEED)
    cfg.update(kw)
    return cfg


def _sigma(q, s):
    return math.sqrt(max(q * (1 - q), 1e-12) / s) if s else float('nan')


# --- Part A: the static point --------------------------------------------

ROWS = (
    ('as published',            {}),
    ('stray light off',         dict(dark_count_rate=15.0)),
    ('20 ns gate',              dict(gate_width=20e-9)),
    ('extinction 0.0101 (b)',   dict(extinction_epsilon=0.0101)),
    ('extinction 0.0200 (a)',   dict(extinction_epsilon=0.0200)),
)


def static_rows(n, failures, quick=False):
    standard = 200 if quick else TARGET_SIFTED
    print(f"\n  Part A -- the static point, {n:,} pulses per row")
    print(f"  {'':<24}{'QBER':>18}{'sifted':>10}{'errors':>9}")
    out = {}
    for label, kw in ROWS:
        r = simulate_bb84_duplinskiy(n, **_base(**kw))
        q, s, e = r['qber'], r['n_sifted'], r['n_errors']
        sig = _sigma(q, s)
        out[label] = (q, sig, s, e, r)
        print(f"  {label:<24}{100 * q:10.2f} +/-{100 * sig:5.2f} %{s:10,}{e:9,}")
        if s < standard:
            failures.append(f"row '{label}': only {s} sifted, below "
                            f"{standard}; not quotable")
    print(f"\n  the paper: lower bound about {100 * PAPER_QBER_LOWER:.1f} %, "
          f"average {100 * PAPER_QBER_MEAN:.1f} %")
    print("    'stray light off' tests the paper's ATTRIBUTION of the floor to")
    print("    channel noise; '20 ns gate' tests its REASON for the 5 ns window.")
    print("    Both extinction rows are reported, neither is chosen -- register")
    print("    A7's factor-of-two ambiguity is not resolved by this run.")
    return out


# --- Part B: the absolute key rates --------------------------------------

def key_rates(baseline, n_alt, failures):
    print("\n  Part B -- the absolute key rates")
    print(f"    the paper gives two: {PAPER_RAW:.0f} bit/s raw 'during key "
          f"distribution regime',")
    print(f"    and {PAPER_SIFTED:.0f} bit/s sifted.  They decompose exactly:")
    print(f"      {PAPER_SIFTED:.0f} / {PAPER_RAW:.0f} = "
          f"{PAPER_SIFTED / PAPER_RAW:.4f} = sifting (1/2) x duty ({DUTY:g})")
    print(f"      {PAPER_RAW:.0f} x 0.5 x {DUTY:g} = "
          f"{PAPER_RAW * 0.5 * DUTY:.1f}   (the paper says {PAPER_SIFTED:.0f})")
    print("    so the raw figure is measured while distributing and the sifted")
    print("    one carries the 20 % tuning downtime.  Ours, same convention:")

    r = baseline[4]
    raw = r['raw_key_rate']
    sifted = r['sifted_key_rate'] * DUTY
    f_raw, f_sift = raw / PAPER_RAW, sifted / PAPER_SIFTED
    print(f"\n      {'':<22}{'ours':>10}{'paper':>10}{'ratio':>9}")
    print(f"      {'raw clicks':<22}{raw:10.0f}{PAPER_RAW:10.0f}{f_raw:8.2f}x")
    print(f"      {'sifted (x duty)':<22}{sifted:10.0f}{PAPER_SIFTED:10.0f}"
          f"{f_sift:8.2f}x")
    common = 10 * math.log10(0.5 * (f_raw + f_sift))
    print(f"\n    one common factor on both, worth {common:.2f} dB of unmodelled")
    print("    loss -- about what a WDM filter costs, and the paper adds one at")
    print("    Bob for this experiment without quoting its insertion loss.")
    print("    Nothing here is fitted: mu, the losses, eta, the gate and the")
    print("    repetition rate are all stated values.")

    # The competing reading of "13 dB", settled by an independent number
    # rather than by preference.  Under an end-to-end reading Bob's 2 dB is
    # inside the 13, so the fibre carries only 11.
    alt = simulate_bb84_duplinskiy(n_alt, **_base(alpha_dB=11.0 / 30.0))
    f_alt = alt['raw_key_rate'] / PAPER_RAW
    print(f"\n    the competing reading of '13 dB' -- end-to-end, Bob included:")
    print(f"      raw {alt['raw_key_rate']:.0f} /s against {PAPER_RAW:.0f} "
          f"-> {f_alt:.2f}x")
    print("      a clean factor of two, not a plausible connector.  The line-only")
    print("      reading is selected by measurement, not by taste.")

    if f_raw > 2.0 or f_sift > 2.0:
        failures.append(f"the key rate is off by {max(f_raw, f_sift):.2f}x; "
                        "more than a factor of two means a missing mechanism, "
                        "not an unquoted connector")
    return raw, sifted, f_raw, f_sift, alt['raw_key_rate']


# --- Part C: Fig. 7 as a block trace -------------------------------------

def trace(n_blocks, block_sifted, failures):
    block_size = int(block_sifted / YIELD)
    n = block_size * n_blocks
    print(f"\n  Part C -- Fig. 7, as a block trace")
    print(f"    {n_blocks} blocks x {block_size:,} pulses = {n:,} pulses,")
    print(f"    {n / REP:.0f} s of link time at {REP / 1e6:g} MHz")
    print(f"    (the paper's 20 hours would be {20 * 3600 * REP:.1e} pulses, "
          "about 11 days)")

    r = simulate_bb84_duplinskiy(n, block_size=block_size, **_base())
    # Empty blocks are dropped only to keep the chi-square finite; at the
    # operating block size (~300 sifted) the filter never fires, and the
    # guard below fails the run if it ever starts to, because dropping
    # empty blocks would bias the distribution rather than tidy it.
    blocks = [(s, e) for (s, e) in r['blocks'] if s > 0]
    if len(blocks) < len(r['blocks']):
        failures.append(f"{len(r['blocks']) - len(blocks)} blocks carried no "
                        "sifted bits and were dropped; the block is too short "
                        "and the distribution is biased")
    if len(blocks) < n_blocks // 2:
        failures.append(f"only {len(blocks)} of {n_blocks} blocks carried any "
                        "sifted bits; the block is too short to plot")
        return None

    qs = [e / s for (s, e) in blocks]
    tot_s = sum(s for s, _ in blocks)
    tot_e = sum(e for _, e in blocks)
    q_bar = tot_e / tot_s
    mean_s = tot_s / len(blocks)

    qs_sorted = sorted(qs)
    lo = qs_sorted[int(0.05 * len(qs_sorted))]
    hi = qs_sorted[min(int(0.95 * len(qs_sorted)), len(qs_sorted) - 1)]

    print(f"\n    pooled QBER      {100 * q_bar:.2f} +/-{100 * _sigma(q_bar, tot_s):.2f} % "
          f"over {tot_s:,} sifted")
    print(f"    the paper        {100 * PAPER_QBER_MEAN:.1f} % average, "
          f"about {100 * PAPER_BAND[0]:.0f}-{100 * PAPER_BAND[1]:.0f} % across Fig. 7")
    print(f"    block spread     5th-95th percentile "
          f"{100 * lo:.2f}-{100 * hi:.2f} %, {int(mean_s)} sifted per block")

    # The pre-registered test.  Pearson chi-square against a common rate:
    # under binomial sampling X ~ chi2(B-1), so X/(B-1) is the variance
    # ratio and sqrt(2/(B-1)) its own standard error.
    dof = len(blocks) - 1
    X = sum((e - s * q_bar) ** 2 / (s * q_bar * (1 - q_bar))
            for (s, e) in blocks)
    ratio = X / dof
    ratio_sig = math.sqrt(2.0 / dof)
    print(f"\n    the pre-registered test: is the block scatter binomial?")
    print(f"      Var(measured) / Var(binomial) = {ratio:.3f} +/- {ratio_sig:.3f}"
          f"   ({dof + 1} blocks)")
    if ratio - 1.0 > 3 * ratio_sig:
        print(f"      -> INFLATED by {100 * (ratio - 1):.0f} %. Afterpulsing "
              "correlates consecutive")
        print("         clicks, so part of Fig. 7's width needs no drift at all.")
    elif 1.0 - ratio > 3 * ratio_sig:
        print("      -> SUPPRESSED below binomial, which nothing in the model")
        print("         predicts and is worth chasing before it is quoted.")
    else:
        print("      -> consistent with pure counting statistics. Afterpulse")
        print("         correlation does not survive to the block scale; the")
        print("         13 us dead time is far shorter than a block.")

    print(f"\n    what is NOT claimed: our per-block sigma is "
          f"{100 * _sigma(q_bar, mean_s):.2f} pp BY CHOICE of")
    print("    block size, which is presentation, not physics. The paper's band")
    print("    is wider; closing that gap would need its drift rate, which it")
    print("    never states and which is not invented here.")
    return blocks, q_bar, mean_s, ratio, ratio_sig, block_size


# --- reporting -----------------------------------------------------------

def _stem(quick):
    """Smoke runs write to their own files.

    They used to share the full run's paths, so `--quick` silently replaced
    a 1.08e9-pulse artifact with a 43M-pulse one, and `--figure-only` then
    redrew from the smoke data.  Nothing warned; the PNG simply became
    under-powered.  This is sec. 32.6's note about `val_duplinskiy/` holding
    regenerated artifacts, with teeth.
    """
    return 'val_duplinskiy_urban--quick' if quick else 'val_duplinskiy_urban'


def _write_csv(rows, tr, rates, quick=False):
    path = os.path.join(OUT_DIR, _stem(quick) + '.csv')
    with open(path, 'w') as fh:
        fh.write("# Duplinskiy urban link, validate_duplinskiy_urban.py\n")
        fh.write(f"# {KM} km, alpha={ALPHA:.4f} dB/km (13 dB lumped), "
                 f"bob={BOB_LOSS} dB, mu={MU}, gate={GATE:g}, rep={REP:g}, "
                 f"stray={STRAY} Hz, seed={SEED}\n")
        fh.write("section,label,qber,qber_sigma,n_sifted,n_errors\n")
        for label, (q, sig, s, e, _) in rows.items():
            fh.write(f"static,{label},{q:.6f},{sig:.6f},{s},{e}\n")
        if rates is not None:
            raw, sifted, f_raw, f_sift, alt = rates
            fh.write(f"rate,raw_per_s,{raw:.4f},,,\n")
            fh.write(f"rate,sifted_per_s,{sifted:.4f},,,\n")
            fh.write(f"rate,raw_ratio_to_paper,{f_raw:.4f},,,\n")
            fh.write(f"rate,sifted_ratio_to_paper,{f_sift:.4f},,,\n")
            fh.write(f"rate,raw_per_s_13dB_endtoend,{alt:.4f},,,\n")
        if tr is not None:
            blocks, q_bar, mean_s, ratio, ratio_sig, block_size = tr
            fh.write(f"trace,block_size_pulses,{block_size},,,\n")
            fh.write(f"trace,variance_ratio,{ratio:.6f},{ratio_sig:.6f},,\n")
            for i, (s, e) in enumerate(blocks):
                fh.write(f"block,{i},{e / s:.6f},,{s},{e}\n")
    print(f"\n  CSV: {path}")


def _read_csv(quick=False):
    path = os.path.join(OUT_DIR, _stem(quick) + '.csv')
    if not os.path.exists(path):
        return None
    rows, blocks, meta = {}, [], {}
    with open(path) as fh:
        for line in fh:
            if line.startswith('#') or line.startswith('section'):
                continue
            parts = line.rstrip('\n').split(',')
            kind, label = parts[0], parts[1]
            if kind == 'static':
                rows[label] = (float(parts[2]), float(parts[3]),
                               int(parts[4]), int(parts[5]), None)
            elif kind == 'block':
                blocks.append((int(parts[4]), int(parts[5])))
            elif kind in ('rate', 'trace'):
                meta[label] = (float(parts[2]),
                               float(parts[3]) if parts[3] else None)
    if not blocks:
        return None
    tot_s = sum(s for s, _ in blocks)
    tot_e = sum(e for _, e in blocks)
    tr = (blocks, tot_e / tot_s, tot_s / len(blocks),
          meta.get('variance_ratio', (float('nan'), None))[0],
          meta.get('variance_ratio', (float('nan'), 0.0))[1] or 0.0,
          int(meta.get('block_size_pulses', (0, None))[0]))
    return rows, tr


def _figure(tr, quick=False):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        print("  matplotlib missing, skipping the figure")
        return

    blocks, q_bar, mean_s, ratio, _, block_size = tr
    qs = np.array([100.0 * e / s for (s, e) in blocks])
    t = np.arange(len(qs)) * block_size / REP

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.4, 4.8),
                                   gridspec_kw={'width_ratios': [2, 1]})

    ax1.axhspan(100 * PAPER_BAND[0], 100 * PAPER_BAND[1], color='0.88',
                zorder=0, label='paper Fig. 7, 3-9 %')
    ax1.axhline(100 * PAPER_QBER_MEAN, color='crimson', ls='--', lw=1.5,
                label='paper average, 5.5 %')
    ax1.plot(t, qs, color='tab:blue', lw=0.9, marker='o', ms=2.6,
             label='simulated blocks')
    ax1.axhline(100 * q_bar, color='tab:blue', ls=':', lw=1.5,
                label=f'our average, {100 * q_bar:.2f} %')
    ax1.set_xlabel('elapsed link time (s)')
    ax1.set_ylabel('QBER (%)')
    ax1.set_ylim(0, max(10.0, qs.max() + 1.5))
    ax1.set_title(f'QBER per {block_size / REP:.1f} s block, 30 km urban '
                  f'line at 13 dB', fontsize=11)
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=8, loc='upper right', ncol=2)

    ax2.hist(qs, bins=18, color='tab:blue', alpha=0.75, edgecolor='0.3',
             density=True, label='simulated blocks')
    sig = 100 * math.sqrt(q_bar * (1 - q_bar) / mean_s)
    x = np.linspace(qs.min() - 2 * sig, qs.max() + 2 * sig, 400)
    ax2.plot(x, np.exp(-0.5 * ((x - 100 * q_bar) / sig) ** 2)
             / (sig * math.sqrt(2 * math.pi)), color='k', lw=1.6,
             label=f'binomial, sigma = {sig:.2f} pp')
    ax2.axvline(100 * PAPER_QBER_MEAN, color='crimson', ls='--', lw=1.5,
                label='paper average')
    ax2.set_xlabel('QBER (%)')
    ax2.set_ylabel('density')
    ax2.set_title(f'block distribution, variance ratio {ratio:.2f}',
                  fontsize=11)
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=8)

    fig.tight_layout()
    png = os.path.join(OUT_DIR, _stem(quick) + '.png')
    fig.savefig(png, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  figure: {png}")


def run(quick=False, figure_only=False):
    os.makedirs(OUT_DIR, exist_ok=True)
    failures = []
    print("=" * 74)
    print("Duplinskiy: the 30 km urban line, its key rates, and Fig. 7")
    print("=" * 74)
    print(f"  {KM} km at {ALPHA:.4f} dB/km (13 dB, lumped) + {BOB_LOSS} dB Bob,"
          f" mu {MU}, {GATE * 1e9:g} ns gate,")
    print(f"  {REP / 1e6:g} MHz, {STRAY:g} Hz uncorrelated counts "
          f"(200 stray + 15 dark), seed {SEED}")

    if figure_only:
        loaded = _read_csv(quick)
        if loaded is None:
            print("\n  no previous run to draw; run without --figure-only first")
            return 1
        print(f"\n  redrawing from {_stem(quick)}.csv, no simulation")
        _figure(loaded[1], quick)
        return 0

    print("\n  stated before the run: the block scatter is predicted to be")
    print("  binomial unless afterpulse correlation survives to the block")
    print("  scale. The verdict is read off the variance ratio, not chosen.")

    # the negative control every parameter in this project carries
    a = simulate_bb84_duplinskiy(200_000, seed=SEED)
    b = simulate_bb84_duplinskiy(200_000, seed=SEED, block_size=20_000)
    same = ((a['n_sifted'], a['n_errors'], a['qber'])
            == (b['n_sifted'], b['n_errors'], b['qber']))
    partition = (sum(s for s, _ in b['blocks']),
                 sum(e for _, e in b['blocks'])) == (b['n_sifted'],
                                                     b['n_errors'])
    print(f"\n  control: block_size leaves the frozen baseline "
          f"{'bit-identical' if same else 'MOVED'}, "
          f"blocks partition the run: {partition}")
    if not same:
        failures.append("block_size moved the frozen sec. 27.1 baseline; it "
                        "must consume no RNG")
    if not partition:
        failures.append("the blocks do not sum to the run totals")

    n_static = int(TARGET_SIFTED / YIELD * 1.15)
    if quick:
        n_static //= 12
    rows = static_rows(n_static, failures, quick=quick)
    rates = key_rates(rows['as published'], n_static // 4, failures)
    tr = trace(N_BLOCKS if not quick else 24,
               BLOCK_SIFTED if not quick else 60, failures)

    _write_csv(rows, tr, rates, quick)
    if tr is not None:
        _figure(tr, quick)

    print()
    if failures:
        print("[FAIL]")
        for f in failures:
            print(f"  - {f}")
        return 1
    print("[PASS] the urban configuration runs at quotable statistics, both of")
    print("       the paper's absolute key rates reproduce to one common")
    print("       factor, and Fig. 7's distribution is measured")
    return 0


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--quick', action='store_true',
                    help='a twelfth of the pulses, for a smoke run')
    ap.add_argument('--figure-only', action='store_true',
                    help="redraw from the last run's CSV without simulating")
    a = ap.parse_args()
    sys.exit(run(quick=a.quick, figure_only=a.figure_only))
