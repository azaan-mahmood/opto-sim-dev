"""Impairment scenarios for the polarisation-encoded (Duplinskiy) chain.

The counterpart to `val_system_scenarios.py`, which does this for the
time-bin chain.  Every row is an explicit, self-contained configuration run
at a recorded seed and pulse count, and the script refuses to emit a table
whose rows cannot support their conclusions.

No new physics.  Every number here has already been measured by one of the
six validators; this puts them in one place, at one seed, with one
statistical standard, so the rows can be compared to each other rather than
across scripts.

    validate_dfb_duplinskiy.py          the DFB source through the chain
    validate_duplinskiy_drift.py        temperature and bend mismatch
    validate_duplinskiy_extinction.py   analyser extinction, A1
    validate_duplinskiy_calibration.py  Fig. 6 and Table 1
    validate_duplinskiy_birefringence.py  the compensation swing
    validate_duplinskiy_urban.py        the 30 km urban line
    validate_duplinskiy_dispersion.py   CD and PMD

Two things to read carefully
---------------------------
**The source changes between rows.** Rows 1-4 and 8-13 use the flat
analytic field this chain has always used, one time sample of diagonal
polarisation.  Rows 5-7 use the DFB device model, because CD and PMD are
frequency-domain operators and cannot act on one sample.  Row 5 is the
control for rows 6 and 7: it is the same source with no dispersion, so the
comparison isolates the impairment rather than the source swap.

**One row is at a different distance.** The urban row is the paper's second
experiment, a 30 km deployed line at 13 dB, and it is marked as such.
Everything else is the 50 km spool the paper's lab measurement used.

References
----------
[1] Duplinskiy, Ustimchik, Kanapin, Kurochkin & Kurochkin, Opt. Express
    25(23), 28886-28897 (2017).
[2] ID Quantique, ID230 InGaAs SPAD datasheet.
"""
import argparse
import math
import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from src.lasers import DFBLaser, DriveParams, LaserDriver
from src.protocols.bb84_duplinskiy import simulate_bb84_duplinskiy

OUT = os.path.join(os.path.dirname(__file__), 'val_duplinskiy')

DISTANCE_DEFAULT = 50.0
TARGET_SIFTED_DEFAULT = 3000
MIN_SIFTED = 2400          # the write guard, deliberately below the target
PILOT_BITS = 400_000
CEILING_DEFAULT = 400_000_000

# The DFB operating point, matching validate_duplinskiy_dispersion.py.
LASER_SEED = 11
N_SECTIONS = 15
I_BIAS, I_PEAK = 0.060, 0.140
HARVEST_PERIOD = 5e-9
PULSE_WIDTH = 150e-12
N_SAMPLES = 8192

_SOURCE = {}


def dfb_source():
    """The DFB field on the device grid, built once and reused.

    On the device grid deliberately: `sample_field` decimates by averaging,
    which band-limits, and PMD acts through the source bandwidth -- so a
    coarser grid understates the impairment systematically rather than
    noisily (measured in validate_duplinskiy_dispersion.py).
    """
    if 'field' not in _SOURCE:
        las = DFBLaser(n_sections=N_SECTIONS, seed=LASER_SEED)
        drive = DriveParams(mode='gain_switched', waveform='gaussian',
                            i_bias=I_BIAS, i_peak=I_PEAK,
                            period=HARVEST_PERIOD, width=PULSE_WIDTH)
        drv = LaserDriver(las, drive, seed=LASER_SEED,
                          polarization_azimuth=np.pi / 4)
        _SOURCE['field'] = drv.sample_field(las.dt, N_SAMPLES)
        _SOURCE['dt'] = las.dt
    return dict(source_field=_SOURCE['field'], source_dt=_SOURCE['dt'])


# --- The rows.  Each config is complete; nothing is inherited. -----------
#
# `compensate=True` is the paper's operating mode: the channel is
# quasi-static, so Bob inverts its Jones matrix.  Rows that leave it on and
# change only the fibre state are EXACTLY inert by construction, which is
# why the drift rows set a calibration mismatch instead.
SCENARIOS = [
    dict(name='No impairments',
         config=dict(alpha_dB=0.0, model='auto', compensate=True,
                     afterpulse_prob=0.0, dark_count_rate=0.0),
         note='floor: sifting loss and double clicks only'),
    dict(name='Attenuation, birefringence compensated',
         config=dict(compensate=True, afterpulse_prob=0.0,
                     dark_count_rate=0.0),
         note='birefringence cannot be switched off from this chain, and '
              'with compensation on it is exactly inert: U_comp = J^dagger'),
    dict(name='Birefringence uncompensated',
         config=dict(compensate=False, afterpulse_prob=0.0,
                     dark_count_rate=0.0),
         note='positive control; a fixed SU(2) scrambles the encoding'),
    dict(name='DFB source, no dispersion',
         config=dict(compensate=True, afterpulse_prob=0.0,
                     dark_count_rate=0.0),
         source='dfb',
         note='control for the two rows below; same source, no CD or PMD'),
    dict(name='+ Chromatic dispersion',
         config=dict(compensate=True, afterpulse_prob=0.0,
                     dark_count_rate=0.0, cd=True),
         source='dfb',
         note='structural null: H(f) multiplies both components equally'),
    dict(name='+ PMD',
         config=dict(compensate=True, afterpulse_prob=0.0,
                     dark_count_rate=0.0, pmd=True),
         source='dfb',
         note="opposite phase ramps; the paper's sec. 4 chirp mechanism"),
    dict(name='+ Afterpulsing 5 %',
         config=dict(compensate=True, afterpulse_prob=0.05,
                     dark_count_rate=0.0),
         note='ID230 datasheet value; register A1'),
    dict(name='+ Dark counts 15 Hz',
         config=dict(compensate=True, afterpulse_prob=0.05,
                     dark_count_rate=15.0),
         note='ID230 datasheet value'),
    dict(name='+ Extinction 0.0101',
         config=dict(compensate=True, afterpulse_prob=0.05,
                     dark_count_rate=15.0, extinction_epsilon=0.0101),
         note="visibility-like reading of the paper's 98 %; register A7"),
    dict(name='+ Extinction 0.0200',
         config=dict(compensate=True, afterpulse_prob=0.05,
                     dark_count_rate=15.0, extinction_epsilon=0.0200),
         note='power-fraction reading of the same 98 %; register A7'),
    dict(name='+ Drift 0.05 C after calibration',
         config=dict(compensate=True, afterpulse_prob=0.05,
                     dark_count_rate=15.0, temperature=25.05,
                     calibration_temperature=25.0),
         note='past the 0.02 C tolerance measured at 50 km'),
    dict(name='+ Bent to 0.3 m after calibration',
         config=dict(compensate=True, afterpulse_prob=0.05,
                     dark_count_rate=15.0, bend_radius=0.3,
                     calibration_bend_radius=None),
         note='Ulrich 0.135*(r/R)^2, calibrated straight then disturbed'),
    dict(name="Paper's lab configuration",
         config=dict(compensate=True, afterpulse_prob=0.05,
                     dark_count_rate=15.0, mu=0.1, gate_width=20e-9,
                     rep_rate=10e6),
         note='50 km spool, 10 dB, 20 ns gate; the paper reports 2 %. '
              'DELIBERATELY identical to the dark-count row: the ladder '
              'reaches the paper\'s configuration exactly at that rung, and '
              'spelling it out is worth one repeated line'),
    dict(name="Paper's urban line",
         config=dict(compensate=True, afterpulse_prob=0.05,
                     dark_count_rate=215.0, mu=0.02, gate_width=5e-9,
                     rep_rate=5e6, alpha_dB=13.0 / 30.0, bob_loss_dB=2.0),
         distance_km=30.0,
         note='13 dB over 30 km, 200 Hz stray light; the paper reports 5.5 %'),
]


def qber_err(q, n):
    return math.sqrt(max(q * (1 - q), 1e-12) / n) if n else float('nan')


def build(scenario, distance, seed):
    cfg = dict(scenario['config'])
    cfg.setdefault('fiber_length', scenario.get('distance_km', distance))
    cfg['seed'] = seed
    if scenario.get('source') == 'dfb':
        cfg.update(dfb_source())
    return cfg


def run_to_target(config, target_sifted, ceiling, pilot_bits=PILOT_BITS,
                  max_retries=3):
    """Grow the pulse count until `target_sifted` sifted bits are collected.

    Same contract as `val_system_scenarios.run_to_target`: pilot, scale,
    re-run from the same seed, and retry on undershoot, because a one-shot
    extrapolation from a noisy pilot lands short about half the time.
    """
    n = int(min(pilot_bits, ceiling))
    r = simulate_bb84_duplinskiy(n, **config)
    for _ in range(max_retries):
        if r['n_sifted'] >= target_sifted or n >= ceiling:
            break
        if r['n_sifted'] == 0:
            n_next = ceiling
        else:
            frac = r['n_sifted'] / n
            headroom = 1.15 + 3.0 / math.sqrt(max(r['n_sifted'], 1))
            n_next = int(math.ceil(target_sifted / frac * headroom))
        n_next = int(min(max(n_next, n + 1), ceiling))
        if n_next <= n:
            break
        n = n_next
        r = simulate_bb84_duplinskiy(n, **config)
    return n, r


def check_statistical_power(rows, min_sifted, allow_underpowered):
    weak = [(n, s) for n, s in rows if s < min_sifted]
    if not weak:
        return
    detail = '; '.join(f"{n}: {s} sifted" for n, s in weak)
    if allow_underpowered:
        print(f"\n  !! WARNING: {len(weak)} row(s) below {min_sifted} sifted "
              f"({detail}).")
        print("     Emitting anyway because --allow-underpowered was given.")
        print("     This table is a smoke run. Do NOT cite it as a result.")
        return
    raise RuntimeError(
        f"Refusing to write an under-powered table: {len(weak)} row(s) have "
        f"fewer than {min_sifted} sifted bits ({detail}).\n"
        f"Rows that differ only in which impairment is enabled come out "
        f"bit-identical at low counts, and the table then cannot "
        f"distinguish 'this impairment does nothing' from 'this impairment "
        f"was never applied'.\n"
        f"Fix: raise --target-sifted or --ceiling. To emit a smoke run "
        f"anyway, pass --allow-underpowered."
    )


def commit_hash():
    try:
        return subprocess.run(
            ['git', 'rev-parse', '--short', 'HEAD'], capture_output=True,
            text=True, cwd=os.path.dirname(__file__)).stdout.strip() or 'unknown'
    except Exception:
        return 'unknown'


def _write_csv(rows, seed, distance, target, commit, stem):
    path = os.path.join(OUT, stem + '.csv')
    with open(path, 'w', encoding='utf-8') as f:
        f.write("# Impairment scenarios, polarisation-encoded (Duplinskiy) "
                "chain\n")
        f.write(f"# script: analysis/val_duplinskiy_scenarios.py  "
                f"seed: {seed}  distance: {distance:g} km  "
                f"target sifted: {target}  commit: {commit}\n")
        f.write("# chain: source -> PM1 -> FiberRealization -> compensation "
                "-> VOA -> PM2 -> circular analyser -> 2x SPAD\n")
        f.write("# qber_err: binomial sqrt(q(1-q)/n_sifted)\n")
        f.write("scenario,distance_km,qber,qber_err,n_sifted,n_errors,"
                "pulses,source,note,config\n")
        for sc, cfg, r, err, used in rows:
            src = 'DFB' if sc.get('source') == 'dfb' else 'analytic'
            shown = {k: v for k, v in cfg.items()
                     if k not in ('source_field', 'source_dt', 'seed')}
            f.write(f"\"{sc['name']}\",{cfg['fiber_length']:g},"
                    f"{r['qber']:.6f},{err:.6f},{r['n_sifted']},"
                    f"{r['n_errors']},{used},{src},\"{sc['note']}\","
                    f"\"{shown}\"\n")
    print(f"\nSaved: {path}")


def _write_markdown(rows, seed, distance, target, commit, stem):
    """A markdown table, not LaTeX.

    `val_system_scenarios.py` emits .tex for the manuscript.  This one does
    not: it feeds the working record, and a manuscript table is the
    author's to place.
    """
    path = os.path.join(OUT, stem + '.md')
    # Explicit utf-8: the default codec on this platform is cp1252 and
    # silently mangles the +/- sign into a replacement character.  Same
    # class of defect as .
    with open(path, 'w', encoding='utf-8') as f:
        f.write("# Impairment scenarios, polarisation-encoded chain\n\n")
        f.write(f"Generated by `analysis/val_duplinskiy_scenarios.py` at "
                f"commit `{commit}`, seed {seed}, {distance:g} km unless "
                f"the row says otherwise, every row grown to at least "
                f"{target} sifted bits.\n\n")
        f.write("| scenario | km | source | QBER | sifted | pulses | note |\n")
        f.write("|---|---|---|---|---|---|---|\n")
        for sc, cfg, r, err, used in rows:
            src = 'DFB' if sc.get('source') == 'dfb' else 'analytic'
            f.write(f"| {sc['name']} | {cfg['fiber_length']:g} | {src} | "
                    f"{100 * r['qber']:.2f} ± {100 * err:.2f} % | "
                    f"{r['n_sifted']:,} | {used:,} | {sc['note']} |\n")
    print(f"Saved: {path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--distance', type=float, default=DISTANCE_DEFAULT,
                    help='fibre length in km for every row that does not '
                         'set its own (default 50, the paper\'s spool)')
    ap.add_argument('--target-sifted', type=int,
                    default=TARGET_SIFTED_DEFAULT,
                    help='grow each row until this many sifted bits')
    ap.add_argument('--ceiling', type=int, default=CEILING_DEFAULT,
                    help='hard cap on pulses per row')
    ap.add_argument('--allow-underpowered', action='store_true',
                    help='emit even if rows fall below the write guard')
    ap.add_argument('--quick', action='store_true',
                    help='small target, separate output files, for a smoke run')
    a = ap.parse_args()

    if a.quick:
        a.target_sifted = 250
        a.ceiling = min(a.ceiling, 20_000_000)
    stem = ('val_duplinskiy_scenarios--quick' if a.quick
            else f'val_duplinskiy_scenarios--seed{a.seed}')

    os.makedirs(OUT, exist_ok=True)
    commit = commit_hash()
    print("=" * 78)
    print("Impairment scenarios, polarisation-encoded (Duplinskiy) chain")
    print("=" * 78)
    print(f"  seed {a.seed}, {a.distance:g} km unless the row says "
          f"otherwise, commit {commit}")
    print(f"  every row grown to >= {a.target_sifted} sifted bits "
          f"(ceiling {a.ceiling:,} pulses)")

    print(f"\n  {'scenario':<40}{'km':>5}{'src':>10}{'QBER':>18}"
          f"{'sifted':>9}{'pulses':>13}")
    rows = []
    for sc in SCENARIOS:
        cfg = build(sc, a.distance, a.seed)
        used, r = run_to_target(cfg, a.target_sifted, a.ceiling)
        err = qber_err(r['qber'], r['n_sifted'])
        rows.append((sc, cfg, r, err, used))
        src = 'DFB' if sc.get('source') == 'dfb' else 'analytic'
        print(f"  {sc['name']:<40}{cfg['fiber_length']:5.0f}{src:>10}"
              f"{100 * r['qber']:11.2f} +/-{100 * err:5.2f}"
              f"{r['n_sifted']:9,}{used:13,}")

    check_statistical_power([(sc['name'], r['n_sifted'])
                             for sc, _, r, _, _ in rows],
                            250 if a.quick else MIN_SIFTED,
                            a.allow_underpowered)

    _write_csv(rows, a.seed, a.distance, a.target_sifted, commit, stem)
    _write_markdown(rows, a.seed, a.distance, a.target_sifted, commit, stem)
    print(f"\nTotal pulses: {sum(u for _, _, _, _, u in rows):,}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
