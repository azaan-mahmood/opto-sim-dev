"""Impairment scenarios for the polarisation-encoded (Duplinskiy) chain.

Every row is an explicit, self-contained configuration run at a recorded
seed and pulse count, and the script refuses to emit a table whose rows
cannot support their conclusions.

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

What is checked
---------------
Introducing no physics limits what this script may honestly assert, and
`check_claims` stays inside that limit: every check is one row's own
`note` restated as arithmetic.  Four are exact -- the compensated rows
carry zero errors, the DFB source swap is inert without dispersion, CD is
a structural null, and the paper's lab row is bit-identical to the
dark-count row it deliberately repeats.  Four more assert that a row is
not inert, by its error count and never by its size: the uncompensated
control, PMD, and the two calibration mismatches all have magnitudes
that belong to the fibre realization rather than to the model, and move
with `--seed`.  Two comparisons need statistical resolution a smoke
budget does not have, and are skipped on any reduced
run -- `--quick` or `--allow-underpowered` -- rather than asserted and
left to fail at random.  That is the same predicate which names the
output file, so a run cannot be quotable in its filename and unquotable
in its checks, or the reverse.

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

# What the paper reports for the two configurations rows 13 and 14
# reproduce.  Printed for comparison and deliberately NOT asserted; see
# `check_claims` for why the ordering is checked and the offset is not.
PAPER_LAB_QBER = 0.02
PAPER_URBAN_QBER = 0.055
# validate_duplinskiy_urban.py's artifact, read rather than copied.
URBAN_VALIDATOR_CSV = 'val_duplinskiy_urban.csv'
URBAN_VALIDATOR_ROW = 'static,as published,'

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
    dict(key='null', name='No impairments',
         config=dict(alpha_dB=0.0, model='auto', compensate=True,
                     afterpulse_prob=0.0, dark_count_rate=0.0),
         note='floor: sifting loss and double clicks only'),
    dict(key='compensated', name='Attenuation, birefringence compensated',
         config=dict(compensate=True, afterpulse_prob=0.0,
                     dark_count_rate=0.0),
         note='birefringence cannot be switched off from this chain, and '
              'with compensation on it is exactly inert: U_comp = J^dagger'),
    dict(key='uncompensated', name='Birefringence uncompensated',
         config=dict(compensate=False, afterpulse_prob=0.0,
                     dark_count_rate=0.0),
         note='positive control; a fixed SU(2) scrambles the encoding'),
    dict(key='dfb', name='DFB source, no dispersion',
         config=dict(compensate=True, afterpulse_prob=0.0,
                     dark_count_rate=0.0),
         source='dfb',
         note='control for the two rows below; same source, no CD or PMD'),
    dict(key='cd', name='+ Chromatic dispersion',
         config=dict(compensate=True, afterpulse_prob=0.0,
                     dark_count_rate=0.0, cd=True),
         source='dfb',
         note='structural null: H(f) multiplies both components equally'),
    dict(key='pmd', name='+ PMD',
         config=dict(compensate=True, afterpulse_prob=0.0,
                     dark_count_rate=0.0, pmd=True),
         source='dfb',
         note="opposite phase ramps; the paper's sec. 4 chirp mechanism"),
    dict(key='afterpulse', name='+ Afterpulsing 5 %',
         config=dict(compensate=True, afterpulse_prob=0.05,
                     dark_count_rate=0.0),
         note='ID230 datasheet value; register A1'),
    dict(key='dark', name='+ Dark counts 15 Hz',
         config=dict(compensate=True, afterpulse_prob=0.05,
                     dark_count_rate=15.0),
         note='ID230 datasheet value'),
    dict(key='ext_visibility', name='+ Extinction 0.0101',
         config=dict(compensate=True, afterpulse_prob=0.05,
                     dark_count_rate=15.0, extinction_epsilon=0.0101),
         note="visibility-like reading of the paper's 98 %; register A7"),
    dict(key='ext_power', name='+ Extinction 0.0200',
         config=dict(compensate=True, afterpulse_prob=0.05,
                     dark_count_rate=15.0, extinction_epsilon=0.0200),
         note='power-fraction reading of the same 98 %; register A7'),
    dict(key='drift', name='+ Drift 0.05 C after calibration',
         config=dict(compensate=True, afterpulse_prob=0.05,
                     dark_count_rate=15.0, temperature=25.05,
                     calibration_temperature=25.0),
         note='past the 0.02 C tolerance measured at 50 km'),
    dict(key='bend', name='+ Bent to 0.3 m after calibration',
         config=dict(compensate=True, afterpulse_prob=0.05,
                     dark_count_rate=15.0, bend_radius=0.3,
                     calibration_bend_radius=None),
         note='Ulrich 0.135*(r/R)^2, calibrated straight then disturbed'),
    dict(key='lab', name="Paper's lab configuration",
         config=dict(compensate=True, afterpulse_prob=0.05,
                     dark_count_rate=15.0, mu=0.1, gate_width=20e-9,
                     rep_rate=10e6),
         note='50 km spool, 10 dB, 20 ns gate; the paper reports 2 %. '
              'DELIBERATELY identical to the dark-count row: the ladder '
              'reaches the paper\'s configuration exactly at that rung, and '
              'spelling it out is worth one repeated line'),
    dict(key='urban', name="Paper's urban line",
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

    Pilot, scale, re-run from the same seed, and retry on undershoot,
    because a one-shot extrapolation from a noisy pilot lands short about
    half the time.  `validate_gobby.run_to_target` does the same thing for
    the time-bin chain.
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


def _identical(a, b):
    """Whether two rows agree in every count, not merely in QBER.

    QBER alone would pass two rows differing in numerator and denominator
    by the same ratio, which is what a changed budget looks like.
    """
    return (a[2]['n_sifted'] == b[2]['n_sifted']
            and a[2]['n_errors'] == b[2]['n_errors']
            and a[4] == b[4])


def check_claims(rows, reduced):
    """Assert what the rows already claim.  0 on pass, 1 on fail.

    Nothing here is a new criterion.  Every check is one row's `note`
    restated as arithmetic, which is the only honest kind of check this
    script can make: it introduces no physics, so it can verify only that
    the numbers it collects still stand in the relations its own
    annotations assert.

    Budget is respected rather than hoped over.  The exact checks -- the
    nulls, the bit-identities, and the four "this row is not inert" error
    counts -- hold at any pulse count, because an equality does not get
    sharper with statistics.  The two that need statistical resolution
    are skipped on a reduced run and say so: a check that fails at random
    is worse than no check.

    No check is a magnitude, because the magnitudes here are not the
    model's.  0.05 C of calibration mismatch costs 35.9 % QBER at seed 42
    and 4.0 % at seed 7: how far the Jones matrix turns per degree is a
    property of the fibre draw.  Any threshold separating those two is a
    fitted number wearing a physical name.

    `reduced` is the predicate that also names the output file, and not
    `--quick` -- `--allow-underpowered` waives the write guard at any
    target, so it too makes a run unquotable.  A run cannot be quotable
    in its filename and unquotable in its checks, or the reverse.

    The gate is the budget and never the observed gap.  A ladder that had
    genuinely stopped climbing would show a gap near zero, hence no
    resolution, hence a skip -- and the check would never fire on the one
    case it exists for.
    """
    by = {sc['key']: (sc, cfg, r, err, used)
          for sc, cfg, r, err, used in rows}
    fail, ok, skipped = [], [], []

    # -- exact: compensation opens no error channel at all ---------------
    dirty = [(k, by[k][2]['n_errors']) for k in
             ('null', 'compensated', 'dfb', 'cd') if by[k][2]['n_errors']]
    if dirty:
        fail.append(
            "compensation is not exactly inert: "
            + "; ".join(f"{by[k][0]['name']} has {n} error(s)"
                        for k, n in dirty)
            + ". U_comp = J^dagger leaves no error channel open, so a "
              "non-zero count here is the inverse failing to reach the "
              "field rather than a small residual")
    else:
        ok.append("compensation is exactly inert -- four rows at zero "
                  "errors, not merely at a small QBER")

    # -- exact: the source swap changes nothing without dispersion -------
    if not _identical(by['dfb'], by['compensated']):
        fail.append(
            f"the DFB source is not inert without dispersion: "
            f"{by['dfb'][2]['n_sifted']} sifted from {by['dfb'][4]:,} "
            f"pulses against {by['compensated'][2]['n_sifted']} from "
            f"{by['compensated'][4]:,} on the analytic source. With cd and "
            f"pmd both off the fibre is frequency-independent, so only the "
            f"pulse energy reaches the detector and the two sources must "
            f"agree. If they do not, the CD and PMD rows below are "
            f"measuring the source swap as well as the impairment")
    else:
        ok.append("the DFB source swap is inert without dispersion, so the "
                  "CD and PMD rows isolate the impairment")

    # -- exact: CD is a structural null ----------------------------------
    if not _identical(by['cd'], by['dfb']):
        fail.append(
            f"chromatic dispersion is no longer a structural null: "
            f"{by['cd'][2]['n_sifted']} sifted, {by['cd'][2]['n_errors']} "
            f"errors against {by['dfb'][2]['n_sifted']} and "
            f"{by['dfb'][2]['n_errors']} on the same source with CD off. "
            f"H(f) is a scalar on this chain and multiplies both "
            f"polarisation components equally, so it cannot rotate the "
            f"encoding")
    else:
        ok.append("chromatic dispersion is a structural null -- H(f) "
                  "multiplies both components equally")

    # -- exact: the paper's lab row IS the dark-count row ----------------
    if not _identical(by['lab'], by['dark']):
        fail.append(
            f"the paper's lab configuration no longer reproduces the "
            f"dark-count row ({by['lab'][2]['n_sifted']} sifted / "
            f"{by['lab'][2]['n_errors']} errors against "
            f"{by['dark'][2]['n_sifted']} / {by['dark'][2]['n_errors']}). "
            f"The lab row spells out mu=0.1, gate_width=20 ns and "
            f"rep_rate=10 MHz where the dark-count row relies on those "
            f"being the defaults of simulate_bb84_duplinskiy, so one of "
            f"the three defaults has moved")
    else:
        ok.append("the paper's lab configuration is bit-identical to the "
                  "dark-count row, so the three defaults it spells out "
                  "still are the defaults")

    # -- structural: the rows that must not be inert ---------------------
    #
    # Asserted as "produces errors at all", never as a magnitude.  Those
    # magnitudes belong to the fibre realization and move with `--seed`:
    # 0.05 C of calibration mismatch costs 35.9 % at seed 42 and 4.0 % at
    # seed 7, both correct.  What holds for every draw is that an exact
    # inverse is exactly inert and anything else is not, so that is what
    # is checked; the magnitudes are printed instead.
    not_inert = (
        ('uncompensated', 'the uncompensated positive control',
         'a fixed SU(2) with no compensation has to reach the encoding; '
         'if it does not, every compensated row above it is passing '
         'vacuously'),
        ('pmd', 'PMD',
         'CD and PMD differ in kind -- CD is one scalar H(f), and is a '
         'null two checks above; PMD puts opposite phase ramps on the two '
         'components -- so PMD must not come out a null as well'),
        ('drift', '0.05 C of drift after calibration',
         'compensation calibrated at one fibre state and run at another '
         'is no longer the exact inverse, so it cannot be inert'),
        ('bend', 'a bend to 0.3 m after calibration',
         'as above, calibrated straight and then disturbed'),
    )
    for key, what, why in not_inert:
        r = by[key][2]
        if r['n_errors'] == 0:
            fail.append(f"{what} produced zero errors in {r['n_sifted']} "
                        f"sifted bits, making it exactly as inert as the "
                        f"compensated null. {why}")
        else:
            ok.append(f"{what} is not inert -- {r['n_errors']} errors in "
                      f"{r['n_sifted']} sifted ({100 * r['qber']:.1f} %)")

    # -- resolution-limited: skipped when the budget cannot see it -------
    lo, hi = by['afterpulse'], by['ext_power']
    ladder_gap = hi[2]['qber'] - lo[2]['qber']
    ladder_sigma = math.hypot(lo[3], hi[3])
    urban, lab = by['urban'], by['lab']
    urban_gap = urban[2]['qber'] - lab[2]['qber']
    urban_sigma = math.hypot(urban[3], lab[3])

    if reduced:
        skipped.append(
            f"the impairment ladder ({100 * ladder_gap:+.2f} pp end to "
            f"end) and the urban-over-lab ordering "
            f"({100 * urban_gap:+.2f} pp). One sigma at this budget is "
            f"{100 * ladder_sigma:.2f} pp and {100 * urban_sigma:.2f} pp, "
            f"so neither is resolvable and asserting them here would fail "
            f"at random")
    else:
        if ladder_gap <= 3.0 * ladder_sigma:
            fail.append(
                f"the impairment ladder does not climb: afterpulsing alone "
                f"reads {100 * lo[2]['qber']:.2f} % and the far end of the "
                f"ladder {100 * hi[2]['qber']:.2f} %, a gap of "
                f"{100 * ladder_gap:+.2f} pp against "
                f"{100 * ladder_sigma:.2f} pp of combined sigma. Each rung "
                f"adds an impairment and keeps the ones below it")
        else:
            ok.append(f"the impairment ladder climbs end to end, "
                      f"{100 * lo[2]['qber']:.2f} % to "
                      f"{100 * hi[2]['qber']:.2f} % "
                      f"({ladder_gap / ladder_sigma:.1f} sigma)")

        if urban_gap <= 2.0 * urban_sigma:
            fail.append(
                f"the urban line is not above the lab spool: "
                f"{100 * urban[2]['qber']:.2f} % against "
                f"{100 * lab[2]['qber']:.2f} %, a gap of "
                f"{100 * urban_gap:+.2f} pp against "
                f"{100 * urban_sigma:.2f} pp of combined sigma. The paper "
                f"reports {100 * PAPER_URBAN_QBER:.1f} % over "
                f"{100 * PAPER_LAB_QBER:.0f} %, and the ORDERING is what "
                f"this chain is asked to reproduce")
        else:
            ok.append(f"the urban line sits above the lab spool, "
                      f"{100 * urban[2]['qber']:.2f} % against "
                      f"{100 * lab[2]['qber']:.2f} % "
                      f"({urban_gap / urban_sigma:.1f} sigma), the "
                      f"ordering the paper reports")

    _report_unasserted(by, urban, lab)

    print()
    if fail:
        print("[FAIL]")
        for line in fail:
            print(f"  - {line}")
        return 1
    for line in ok:
        print(f"[PASS] {line}")
    for line in skipped:
        print(f"  SKIPPED at this budget: {line}")
    return 0


def urban_validator_reading():
    """What `validate_duplinskiy_urban.py` measured for the same line.

    Read out of its committed artifact rather than copied into a constant
    here.  A number copied by hand goes out of date when the other script
    changes, and nothing shows that it has.  The claim being made is
    precisely that the two scripts agree, so it should be checked against
    the other script rather than against a copy of what it once said.

    Returns `(None, None)` if the artifact is absent, because a missing
    file means there is no comparison to report, not that it failed.
    """
    path = os.path.join(OUT, URBAN_VALIDATOR_CSV)
    try:
        with open(path, encoding='utf-8') as f:
            for line in f:
                if line.startswith(URBAN_VALIDATOR_ROW):
                    cols = line.split(',')
                    return float(cols[2]), float(cols[3])
    except (OSError, ValueError, IndexError):
        pass
    return None, None


def _report_unasserted(by, urban, lab):
    """The two comparisons worth printing and not worth gating on.

    The offset from the paper, because `validate_duplinskiy_urban.py`
    already carries it and this script introduces no physics with which to
    resolve it.  And the ladder's adjacent rungs, because their steps are
    smaller than this budget resolves -- which is the reason the ladder is
    asserted end to end rather than rung by rung, and is worth showing
    rather than asserting in a docstring.
    """
    print()
    print("  Reported, not asserted:")
    print(f"    urban {100 * urban[2]['qber']:.2f} +/-{100 * urban[3]:.2f} "
          f"% against the paper's {100 * PAPER_URBAN_QBER:.1f} %; "
          f"lab {100 * lab[2]['qber']:.2f} +/-{100 * lab[3]:.2f} % against "
          f"{100 * PAPER_LAB_QBER:.0f} %.")
    q_ref, s_ref = urban_validator_reading()
    if q_ref is None:
        print(f"    {URBAN_VALIDATOR_CSV} is not on disk, so the "
              f"cross-script comparison is unavailable this run.")
    else:
        agrees = abs(urban[2]['qber'] - q_ref) <= 2 * math.hypot(urban[3],
                                                                 s_ref)
        print(f"    The urban offset is carried by "
              f"validate_duplinskiy_urban.py too, which measures "
              f"{100 * q_ref:.2f} +/-{100 * s_ref:.2f} % for the")
        print(f"    same configuration -- the two scripts "
              f"{'agree' if agrees else 'DISAGREE'} with each other, which "
              f"is what this table exists to show. The gap to")
        print("    the paper is a parameter offset and is not this "
              "script's to resolve.")
    print()
    print("    Magnitudes that belong to this fibre realization, not to "
          "the model. These move with --seed and are")
    print("    checked only for being non-zero; 0.05 C of mismatch cost "
          "35.9 % at seed 42 and 4.0 % at seed 7:")
    for key in ('uncompensated', 'pmd', 'drift', 'bend'):
        row = by[key]
        print(f"      {row[0]['name']:<40} {100 * row[2]['qber']:6.2f} "
              f"+/-{100 * row[3]:4.2f} %")
    print()
    print("    Ladder steps. Adjacent rungs are NOT asserted -- see the "
          "sigma column:")
    rungs = ['afterpulse', 'dark', 'ext_visibility', 'ext_power']
    for a_key, b_key in zip(rungs, rungs[1:]):
        a_row, b_row = by[a_key], by[b_key]
        d = b_row[2]['qber'] - a_row[2]['qber']
        s = math.hypot(a_row[3], b_row[3])
        print(f"      {a_row[0]['name']:<22} -> {b_row[0]['name']:<22} "
              f"{100 * d:+6.2f} pp  ({d / s:5.1f} sigma)")


def commit_hash():
    try:
        return subprocess.run(
            ['git', 'rev-parse', '--short', 'HEAD'], capture_output=True,
            text=True, cwd=os.path.dirname(__file__)).stdout.strip() or 'unknown'
    except Exception:
        return 'unknown'


def _stem(base, reduced):
    """Smoke runs write to their own files.

    Same shape as `validate_gobby._stem`, and for the same reason: a run
    that opts out of the statistical-power guard is by definition not
    quotable, so it must not be able to replace an artifact that is.

    `reduced` is driven by both flags that make a run unquotable.
    `--quick` lowers the target; `--allow-underpowered` waives the guard
    at whatever target was asked for.  Either is enough.

    The marker goes last, before the extension: `.gitignore` matches
    `*--quick.csv` and `*--quick.md`, so a name like
    `val_duplinskiy_scenarios--quick--seed42.csv` would fall outside the
    ignore rules.
    """
    return f'{base}--quick' if reduced else base


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

    This feeds the working record rather than a manuscript, and a
    manuscript table is the author's to place.
    """
    path = os.path.join(OUT, stem + '.md')
    # Explicit utf-8: the default codec on this platform is cp1252 and
    # silently mangles the +/- sign into a replacement character.
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
                    help='emit even if rows fall below the write guard. '
                         'Writes to the --quick stem, because a run that '
                         'waives the guard is not quotable and must not be '
                         'able to replace a table that is')
    ap.add_argument('--quick', action='store_true',
                    help='small target, separate output files, for a smoke run')
    a = ap.parse_args()

    if a.quick:
        a.target_sifted = 250
        a.ceiling = min(a.ceiling, 20_000_000)
    stem = _stem(f'val_duplinskiy_scenarios--seed{a.seed}',
                 a.quick or a.allow_underpowered)

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

    # Write before checking, so a run that fails a check still leaves the
    # table that failed it.  Same order as validate_gobby_impairments.
    _write_csv(rows, a.seed, a.distance, a.target_sifted, commit, stem)
    _write_markdown(rows, a.seed, a.distance, a.target_sifted, commit, stem)
    print(f"\nTotal pulses: {sum(u for _, _, _, _, u in rows):,}")
    return check_claims(rows, a.quick or a.allow_underpowered)


if __name__ == '__main__':
    sys.exit(main())
