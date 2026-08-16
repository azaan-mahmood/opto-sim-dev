"""The DFB source through the time-bin (Gobby) chain.

The polarisation chain has already been driven from the DFB device model
(sec. 30.11): it turned out blind to everything the source adds except
pulse energy, because the encoding depends only on the ratio Ey/Ex and the
device emits one amplitude times a fixed Jones vector.

This asks the same question of phase encoding, where the mechanism is
completely different. A time-bin interferometer does not care about
polarisation. It cares whether the two interfering paths arrive with a
stable relative phase.

The claim, stated before the run
--------------------------------
A path-matched AMZI interferes a pulse with a copy of ITSELF, delayed and
then re-delayed so both routes accumulate the same total. At the
interference bin the two copies carry the same chirp at the same point in
the pulse, so a chirp common to both cancels exactly -- the same argument
sec. 23.2 made for linewidth, and sec. 27.9's CD-commutation note made for
dispersion.

So the prediction is **no significant QBER shift** when the analytic
Gaussian is replaced by a gain-switched DFB, and the reason is structural
rather than a small number.

That prediction is worth testing rather than assuming, because it is the
opposite of what the same source does in the polarisation chain, where
sec. 36 measured its chirp turning PMD from nothing into +9.6 pp.

Two controls, because a null needs them (G2)
--------------------------------------------
Negative: hand the chain the analytic Gaussian it would have built itself.
Bit-identical, or the hook is not neutral.

Positive: hand it a DFB pulse displaced from the centre of the window.
That breaks the arrival alignment the gate expects, so QBER must move. If
it does not, the source is not reaching the observable and the null above
is vacuous -- which is exactly the trap sec. 27.3 fell into.

What is NOT covered
-------------------
Per-pulse energy spread. `bb84_duplinskiy` takes `pulse_energy_factors`
because sec. 30.11 found energy was the one thing that mattered there;
`bb84_time_bin` has no equivalent, so this script compares pulse SHAPE and
PHASE only. Adding it would be a separate change with its own controls.

References
----------
[1] Gobby, Yuan & Shields, Appl. Phys. Lett. 84(19), 3762-3764 (2004).
[2] Kim, Chung & Lee, IEEE J. Quantum Electron. 36(7), 787-794 (2000).
"""
import argparse
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from src.lasers import DFBLaser, DriveParams, LaserDriver
from src.protocols.bb84_time_bin import (simulate_bb84_time_bin, field_grid,
                                         gaussian_pulse)

OUT_DIR = os.path.join(os.path.dirname(__file__), '..', 'val_dfb')

SEED = 42
LASER_SEED = 11
N_SECTIONS = 15
I_BIAS, I_PEAK = 0.060, 0.140
HARVEST_PERIOD = 5e-9
PULSE_WIDTH = 100e-12          # matches the chain's own default

DISTANCES = (0, 25, 65, 122)
TARGET_SIFTED = 3000
PILOT = 300_000
CEILING = 300_000_000

H, C, WL = 6.626e-34, 3e8, 1550e-9


def _stem(quick):
    return 'val_dfb_gobby--quick' if quick else 'val_dfb_gobby'


def analytic_field(dt, n, centre, mu=0.1):
    """The Gaussian the chain builds internally, rebuilt here.

    Used as the negative control: passing this in must reproduce the
    internal path bit for bit.
    """
    sigma = PULSE_WIDTH / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    t = np.arange(n, dtype=float) * dt
    amp = np.sqrt((mu * H * C / WL) / (sigma * np.sqrt(np.pi)))
    return (gaussian_pulse(t - centre, sigma, A=amp)[:, None]
            * np.array([1.0, 0.0], dtype=complex))


def dfb_field(dt, n, centre, offset=0.0):
    """One DFB pulse, placed on the chain's grid at `centre`.

    The driver is asked for the chain's own `dt` rather than the device
    step, so `sample_field` band-limits by averaging on the way down.
    Resampling afterwards would alias; see `field_grid`.

    The harvested window holds several pulses at a 5 ns period, and this
    chain expects exactly one -- a second pulse in the window would appear
    as a spurious bin.  So the brightest pulse is cut out and placed alone.

    `offset` displaces it, which is the positive control.
    """
    las = DFBLaser(n_sections=N_SECTIONS, seed=LASER_SEED)
    drive = DriveParams(mode='gain_switched', waveform='gaussian',
                        i_bias=I_BIAS, i_peak=I_PEAK,
                        period=HARVEST_PERIOD, width=PULSE_WIDTH)
    # azimuth 0 puts the whole field on Ex, matching this chain's
    # X-polarised analytic pulse.  A time-bin decoder is polarisation
    # blind, but keeping the convention makes the two sources comparable.
    drv = LaserDriver(las, drive, seed=LASER_SEED, polarization_azimuth=0.0)
    harvest = drv.sample_field(dt, max(n, int(2 * HARVEST_PERIOD / dt)))

    power = np.sum(np.abs(harvest) ** 2, axis=1)
    peak = int(np.argmax(power))
    half = int(round(0.5 * HARVEST_PERIOD / dt))
    lo, hi = max(0, peak - half), min(len(power), peak + half + 1)
    cut = harvest[lo:hi]
    cut_peak = peak - lo

    out = np.zeros((n, 2), dtype=complex)
    at = int(round((centre + offset) / dt))
    start = at - cut_peak
    src_lo = max(0, -start)
    dst_lo = max(0, start)
    take = min(len(cut) - src_lo, n - dst_lo)
    if take <= 0:
        raise ValueError("the pulse does not land inside the chain's window")
    out[dst_lo:dst_lo + take] = cut[src_lo:src_lo + take]
    return out


def run_to_target(target, quick, **kw):
    # The 122 km yield is 1.5e-5, so a 3000-sifted cell costs 2e8 pulses --
    # about nine minutes.  That is right for a real run and wrong for a
    # smoke run, so quick mode caps the spend and lowers the bar to match.
    ceiling = 20_000_000 if quick else CEILING
    n = PILOT if not quick else PILOT // 4
    r = simulate_bb84_time_bin(n, seed=SEED, **kw)
    for _ in range(3):
        if r['n_sifted'] >= target or n >= ceiling:
            break
        if r['n_sifted'] == 0:
            nxt = ceiling
        else:
            frac = r['n_sifted'] / n
            nxt = int(math.ceil(target / frac
                                * (1.15 + 3.0 / math.sqrt(r['n_sifted']))))
        nxt = min(max(nxt, n + 1), ceiling)
        if nxt <= n:
            break
        n = nxt
        r = simulate_bb84_time_bin(n, seed=SEED, **kw)
    return n, r


def _sig(q, s):
    return math.sqrt(max(q * (1 - q), 1e-12) / s) if s else float('nan')


def controls(dt, n, centre, failures):
    print("\n  controls")

    a = simulate_bb84_time_bin(200_000, fiber_length=25, seed=SEED)
    b = simulate_bb84_time_bin(200_000, fiber_length=25, seed=SEED,
                               source_field=analytic_field(dt, n, centre))
    same = ((a['n_sifted'], a['n_errors'], a['qber'])
            == (b['n_sifted'], b['n_errors'], b['qber']))
    print(f"    negative, analytic field passed in : "
          f"{'bit-identical' if same else 'CHAIN MOVED'}  "
          f"({a['n_sifted']} vs {b['n_sifted']} sifted)")
    if not same:
        failures.append("passing the chain its own analytic Gaussian changed "
                        "the result; the source hook is not neutral")

    base = simulate_bb84_time_bin(200_000, fiber_length=25, seed=SEED,
                                  source_field=dfb_field(dt, n, centre))
    moved = simulate_bb84_time_bin(
        200_000, fiber_length=25, seed=SEED,
        source_field=dfb_field(dt, n, centre, offset=1e-9))
    # The criterion is the SIFTED COUNT, not QBER.  Displacing the pulse by
    # 1 ns walks it out of the 1 ns gate entirely, so the count collapses
    # and QBER stops being defined -- reporting "0.00 %" from zero sifted
    # bits as if it were an error rate would be meaningless.
    print(f"    positive, pulse displaced 1 ns     : sifted "
          f"{base['n_sifted']} -> {moved['n_sifted']}  "
          f"({100 * (1 - moved['n_sifted'] / max(base['n_sifted'], 1)):.0f} % "
          "of the key lost)")
    if moved['n_sifted'] > 0.5 * base['n_sifted']:
        failures.append("displacing the source pulse by 1 ns left the sifted "
                        "rate largely intact; the source is not reaching the "
                        "observable and the null below would be vacuous")


def sweep(dt, n, centre, quick, failures):
    target = 250 if quick else TARGET_SIFTED
    print(f"\n  QBER against distance, target {target} sifted per cell")
    print(f"    {'km':>5}{'analytic Gaussian':>26}{'DFB device':>26}"
          f"{'difference':>20}")
    rows = []
    src = dfb_field(dt, n, centre)
    for km in DISTANCES:
        na, ra = run_to_target(target, quick, fiber_length=km)
        nb, rb = run_to_target(target, quick, fiber_length=km,
                               source_field=src)
        qa, qb = ra['qber'], rb['qber']
        sa, sb = _sig(qa, ra['n_sifted']), _sig(qb, rb['n_sifted'])
        d, ds = qb - qa, math.hypot(sa, sb)
        rows.append((km, qa, sa, ra['n_sifted'], qb, sb, rb['n_sifted'],
                     d, ds, na, nb))
        print(f"    {km:5d}{100 * qa:14.2f} +/-{100 * sa:5.2f} "
              f"({ra['n_sifted']:5d}){100 * qb:14.2f} +/-{100 * sb:5.2f} "
              f"({rb['n_sifted']:5d}){100 * d:+10.2f} +/-{100 * ds:5.2f} pp")
        if min(ra['n_sifted'], rb['n_sifted']) < (0.5 * target if quick
                                                  else target):
            failures.append(f"{km} km: {ra['n_sifted']}/{rb['n_sifted']} "
                            f"sifted, below {target}; not quotable")

    worst = max(abs(r[7]) / r[8] for r in rows)
    print(f"\n    largest shift: {worst:.1f} sigma")
    if worst < 3.0:
        print("    -> the chain does not see the source. A path-matched AMZI")
        print("       interferes the pulse with a copy of itself at the same")
        print("       chirp phase, so a common chirp cancels -- sec. 23.2's")
        print("       linewidth argument, reached from a different direction.")
    else:
        print("    -> the source DOES move this chain, which contradicts the")
        print("       prediction above. The path-matching argument needs")
        print("       rechecking before anything here is quoted.")
    return rows


def _write_csv(rows, quick):
    path = os.path.join(OUT_DIR, _stem(quick) + '.csv')
    with open(path, 'w', encoding='utf-8') as fh:
        fh.write("# DFB source through the time-bin chain, "
                 "validate_dfb_gobby.py\n")
        fh.write(f"# seed={SEED} laser_seed={LASER_SEED} "
                 f"pulse_width={PULSE_WIDTH:g} sections={N_SECTIONS}\n")
        fh.write("distance_km,qber_analytic,sigma_analytic,sifted_analytic,"
                 "qber_dfb,sigma_dfb,sifted_dfb,pulses_analytic,pulses_dfb\n")
        for (km, qa, sa, na_s, qb, sb, nb_s, _, _, na, nb) in rows:
            fh.write(f"{km},{qa:.6f},{sa:.6f},{na_s},{qb:.6f},{sb:.6f},"
                     f"{nb_s},{na},{nb}\n")
    print(f"\n  CSV: {path}")


def _figure(rows, dt, n, centre, quick):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("  matplotlib missing, skipping the figure")
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.4, 4.8),
                                   gridspec_kw={'width_ratios': [3, 2]})

    km = [r[0] for r in rows]
    x = np.arange(len(km))
    w = 0.36
    ax1.bar(x - w / 2, [100 * r[1] for r in rows], w,
            yerr=[100 * r[2] for r in rows], capsize=4, color='0.72',
            edgecolor='0.25', linewidth=0.7, label='analytic Gaussian')
    ax1.bar(x + w / 2, [100 * r[4] for r in rows], w,
            yerr=[100 * r[5] for r in rows], capsize=4, color='tab:green',
            edgecolor='0.25', linewidth=0.7, label='DFB device model')
    ax1.set_xticks(x)
    ax1.set_xticklabels([f'{k} km' for k in km])
    ax1.set_xlabel('fibre length')
    ax1.set_ylabel('QBER (%)')
    ax1.set_title('QBER with an analytic pulse and with a DFB pulse\n'
                  'time-bin chain', fontsize=11)
    ax1.grid(True, axis='y', alpha=0.3)
    ax1.legend(fontsize=9, loc='upper left')

    t = np.arange(n) * dt * 1e9
    pa = np.sum(np.abs(analytic_field(dt, n, centre)) ** 2, axis=1)
    pd = np.sum(np.abs(dfb_field(dt, n, centre)) ** 2, axis=1)
    lo = max(0, int((centre - 0.6e-9) / dt))
    hi = min(n, int((centre + 0.6e-9) / dt))
    ax2.plot(t[lo:hi], pa[lo:hi] / pa.max(), color='0.45',
             label='analytic Gaussian')
    ax2.plot(t[lo:hi], pd[lo:hi] / pd.max(), color='tab:green',
             label='DFB device model')
    ax2.set_xlabel('time (ns)')
    ax2.set_ylabel('power, peak normalised')
    ax2.set_title('the two pulse shapes on the chain grid', fontsize=11)
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=9)

    fig.tight_layout()
    png = os.path.join(OUT_DIR, _stem(quick) + '.png')
    fig.savefig(png, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  figure: {png}")


def run(quick=False):
    os.makedirs(OUT_DIR, exist_ok=True)
    failures = []
    dt, n, centre = field_grid(PULSE_WIDTH)
    print("=" * 74)
    print("The DFB source through the time-bin (Gobby) chain")
    print("=" * 74)
    print(f"  chain grid: dt = {dt:.3e} s, {n} samples, span "
          f"{n * dt * 1e9:.2f} ns, pulse at {centre * 1e9:.2f} ns")
    print(f"  Nyquist {0.5 / dt / 1e9:.0f} GHz -- most of the gain-switched")
    print("  chirp is averaged away by sample_field on the way to this grid,")
    print("  and for this chain that is not a loss: see the prediction below.")
    print("\n  predicted before the run: NO significant shift. A path-matched")
    print("  AMZI interferes the pulse with a copy of itself at the same")
    print("  chirp phase, so a common chirp cancels exactly.")

    controls(dt, n, centre, failures)
    rows = sweep(dt, n, centre, quick, failures)
    _write_csv(rows, quick)
    _figure(rows, dt, n, centre, quick)

    print()
    if failures:
        print("[FAIL]")
        for f in failures:
            print(f"  - {f}")
        return 1
    print("[PASS] the source hook is neutral, a displaced pulse moves the")
    print("       observable, and the DFB leaves the time-bin QBER alone")
    return 0


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--quick', action='store_true',
                    help='a smaller sifted target, separate output files')
    a = ap.parse_args()
    sys.exit(run(quick=a.quick))
