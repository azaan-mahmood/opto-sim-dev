#!/usr/bin/env python3
"""Run every validator with a single seed and report a pass/fail summary.

Checks that each script exits 0 and emits its expected output file, and
exits non-zero if any fails.

Usage:
    python run_all.py                      # quick where it is safe
    python run_all.py --full               # publication budget throughout
    python run_all.py --seed 123           # custom seed
    python run_all.py --skip DFB-dupl      # skip by name (repeatable)
    python run_all.py --list               # show the roster and exit

Gobby runs only under ``--full``; see the ``'blocked'`` note below for
why it has no honest cheap setting.

Why the roster is not a flat list of paths
------------------------------------------
The validators do not share a command-line interface, and pretending they
do is how this harness fell behind in the first place.  Three differences
matter, and each is recorded per entry rather than assumed:

``seeded``
    Only the eight original channel validators take ``--seed``.  The DFB
    and Duplinskiy ones do not, and passing it to them is an argparse
    error and an exit code of 2 -- a harness failure indistinguishable, in
    the summary, from a physics failure.

``quick``
    ``'safe'``   -- takes ``--quick`` and writes smoke output to its own
                    ``--quick`` files, so a fast run cannot replace a
                    quotable artifact.  The harness uses it by default.
    ``'shared'`` -- takes ``--quick`` but writes to the SAME paths as the
                    full run.  The harness never passes it.  Running these
                    quick would silently swap a published figure for an
                    under-powered one, which is the exact defect commit
                    2360ad0 fixed elsewhere.  They therefore always run at
                    full budget and cost about eleven minutes between
                    them; giving them the ``_stem(quick)`` guard the other
                    seven already have would let them run quick too.
    ``'blocked'``-- has no honest cheap setting at all, so it runs only
                    under ``--full`` and is reported SKIP otherwise, with
                    the reason, rather than quietly left out.  Gobby is
                    the only one.  Its statistical-power guard refuses to
                    write a table with any row under 500 sifted bits, and
                    since GOBBY-1 corrected the link budget -- eta 0.10 to
                    eta_Bob 0.045, alpha 0.182 to 0.2, cutting the signal
                    3.68x -- clearing that at 122 km needs of order 1e8
                    pulses at that distance alone.  Below it the run
                    fails; above it the run is slow; and because
                    ``validate_gobby`` writes ``val_gobby_table.tex``
                    unconditionally, any middle setting would replace the
                    published 3000-sifted table with a weaker one.  The
                    ``_stem(quick)`` guard would fix this too.
    ``None``     -- no quick mode.  ``validate_duplinskiy_calibration``
                    reaches full at its default ``--pulses 20000``, and
                    ``validate_dfb_reflection`` takes seven seconds.

``output``
    Checked after the run.  ``{seed}`` fills in for the seeded validators;
    ``{q}`` becomes ``--quick`` only when the harness actually passed
    ``--quick``, so a quick run looks for the quick artifact rather than
    passing on a stale full one left behind by an earlier invocation.
"""

import argparse
import os
import subprocess
import sys
import time
from collections import namedtuple

BASE = os.path.dirname(os.path.abspath(__file__))

Validator = namedtuple('Validator', 'name script output seeded quick')


def _v(name, rel_script, rel_output, seeded=False, quick=None):
    return Validator(name,
                     os.path.join(BASE, *rel_script.split('/')),
                     os.path.join(BASE, *rel_output.split('/')),
                     seeded, quick)


VALIDATORS = [
    # --- channel and component validators (original eight) ---------------
    _v('CD', 'analysis/validation/validate_cd.py',
       'analysis/val_cd/val_cd--seed{seed}_table.csv', seeded=True),
    _v('PMD', 'analysis/validation/validate_pmd.py',
       'analysis/val_pmd/val_pmd--seed{seed}_table.csv', seeded=True),
    _v('Attenuation', 'analysis/validation/validate_attenuation.py',
       'analysis/val_attenuation/val_attenuation--seed{seed}_table.csv',
       seeded=True),
    _v('Birefringence', 'analysis/validation/validate_birefringence.py',
       'analysis/val_birefringence/val_birefringence--seed{seed}_table.csv',
       seeded=True),
    _v('APD', 'analysis/validation/validate_apd.py',
       'analysis/val_apd/val_apd--seed{seed}_table.csv', seeded=True),
    _v('CWLaser', 'analysis/validation/validate_cwlaser.py',
       'analysis/val_cwlaser/val_cwlaser--seed{seed}_table.csv', seeded=True),
    _v('MZM', 'analysis/validation/validate_mzm.py',
       'analysis/val_mzm/val_mzm--seed{seed}_table.csv', seeded=True),
    _v('Gobby', 'analysis/val_gobby/validate_gobby.py',
       'analysis/val_gobby/val_gobby_table.tex', seeded=True,
       quick='blocked'),

    # --- DFB device ------------------------------------------------------
    _v('DFB-reflection', 'analysis/validation/validate_dfb_reflection.py',
       'analysis/val_dfb/val_dfb_reflection--N15.csv'),
    _v('DFB-drive', 'analysis/validation/validate_dfb_drive.py',
       'analysis/val_dfb/val_dfb_drive.png', quick='shared'),
    _v('DFB-gobby', 'analysis/validation/validate_dfb_gobby.py',
       'analysis/val_dfb/val_dfb_gobby{q}.csv', quick='safe'),
    _v('DFB-dupl', 'analysis/validation/validate_dfb_duplinskiy.py',
       'analysis/val_dfb/val_dfb_duplinskiy_poincare.png', quick='shared'),

    # --- Duplinskiy polarisation chain -----------------------------------
    _v('DUPL-birefringence',
       'analysis/validation/validate_duplinskiy_birefringence.py',
       'analysis/val_duplinskiy/val_duplinskiy_birefringence{q}.csv',
       quick='safe'),
    _v('DUPL-calibration',
       'analysis/validation/validate_duplinskiy_calibration.py',
       'analysis/val_duplinskiy/val_duplinskiy_calibration.png'),
    _v('DUPL-dispersion',
       'analysis/validation/validate_duplinskiy_dispersion.py',
       'analysis/val_duplinskiy/val_duplinskiy_dispersion{q}.csv',
       quick='safe'),
    _v('DUPL-drift', 'analysis/validation/validate_duplinskiy_drift.py',
       'analysis/val_duplinskiy/val_duplinskiy_drift{q}.csv', quick='safe'),
    _v('DUPL-extinction',
       'analysis/validation/validate_duplinskiy_extinction.py',
       'analysis/val_duplinskiy/val_duplinskiy_extinction{q}.csv',
       quick='safe'),
    _v('DUPL-urban', 'analysis/validation/validate_duplinskiy_urban.py',
       'analysis/val_duplinskiy/val_duplinskiy_urban{q}.csv', quick='safe'),
]

parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--seed', type=int, default=42)
parser.add_argument('--gobby-target-sifted', type=int, default=3000,
                    help='sifted bits per distance point for the Gobby '
                         'validator under --full (default 3000, the '
                         'published budget). A flat --bits count is no '
                         'longer usable: since GOBBY-1 the sifted fraction '
                         'spans ~200x across the sweep, so any flat budget '
                         'either starves 122 km or wastes hours on 0 km.')
parser.add_argument('--full', action='store_true',
                    help='run the quick-capable validators at publication '
                         'budget instead of smoke budget')
parser.add_argument('--skip', action='append', default=[],
                    help='validator name to skip (repeatable, '
                         'case-insensitive)')
parser.add_argument('--list', action='store_true',
                    help='print the roster and exit without running anything')
args = parser.parse_args()

if args.list:
    print(f"{'name':22s} {'seed':5s} {'quick':8s} script")
    for v in VALIDATORS:
        print(f"{v.name:22s} {'yes' if v.seeded else '-':5s} "
              f"{v.quick or '-':8s} {os.path.relpath(v.script, BASE)}")
    print(f"\n{len(VALIDATORS)} validators")
    sys.exit(0)

# Fail loudly on a typo rather than silently running everything.
known = {v.name.lower() for v in VALIDATORS}
skip = {name.lower() for name in args.skip}
unknown = skip - known
if unknown:
    print(f"ERROR: --skip names not in the roster: {', '.join(sorted(unknown))}",
          file=sys.stderr)
    print(f"Known: {', '.join(v.name for v in VALIDATORS)}", file=sys.stderr)
    sys.exit(2)

# A missing script is a harness bug, not a validation failure. Say so before
# spending an hour discovering it one entry at a time.
missing = [v.name for v in VALIDATORS if not os.path.exists(v.script)]
if missing:
    print(f"ERROR: script(s) not found for: {', '.join(missing)}",
          file=sys.stderr)
    sys.exit(2)

mode = 'full' if args.full else 'quick where safe'
print('=' * 64)
print(f"  Validation Harness -- {len(VALIDATORS)} validators, "
      f"seed = {args.seed}, {mode}")
print('=' * 64)
print()

start = time.time()
results = []

for v in VALIDATORS:
    if v.name.lower() in skip:
        print(f"[{v.name}] SKIPPED (--skip)")
        results.append((v.name, 'SKIP', 0.0))
        continue

    # 'blocked' validators have no honest cheap setting, so they are
    # reported as a visible SKIP with the reason rather than silently
    # omitted or run at a budget that would degrade a published artifact.
    if v.quick == 'blocked' and not args.full:
        print(f"[{v.name}] SKIPPED (no honest cheap budget; use --full)")
        results.append((v.name, 'SKIP', 0.0))
        continue

    cmd = [sys.executable, v.script]
    if v.seeded:
        cmd += ['--seed', str(args.seed)]
    if v.name == 'Gobby':
        cmd += ['--target-sifted', str(args.gobby_target_sifted)]

    # Only 'safe' validators may be run quick: the 'shared' ones would
    # overwrite a quotable artifact with a smoke one.
    used_quick = (v.quick == 'safe' and not args.full)
    if used_quick:
        cmd += ['--quick']

    t0 = time.time()
    label = v.name + (' [quick]' if used_quick else '')
    print(f"[{label}] {'-' * max(0, 52 - len(label))}")

    try:
        ret = subprocess.run(cmd, text=True, capture_output=True, timeout=7200)
    except subprocess.TimeoutExpired:
        print(f"[{v.name}] FAILED (timeout after 2 h)")
        results.append((v.name, 'FAIL', time.time() - t0))
        continue

    elapsed = time.time() - t0
    stdout = ret.stdout or ''
    n_pass = stdout.count('[PASS]')
    n_fail_marker = stdout.count('[FAIL]')

    rc_ok = ret.returncode == 0
    out_file = v.output.format(seed=args.seed,
                               q='--quick' if used_quick else '')
    out_ok = os.path.exists(out_file)

    if not rc_ok:
        print(f"[{v.name}] FAILED (rc={ret.returncode})")
        if ret.stderr:
            print(ret.stderr[-2000:])
    if not out_ok:
        print(f"[{v.name}] FAILED (missing output {out_file})")
    # PMD reports through its table CSV rather than [PASS] lines.
    if n_pass == 0 and v.name != 'PMD':
        print(f"[{v.name}] WARNING: no [PASS] markers on stdout (rc=0)")

    status = 'PASS' if (rc_ok and out_ok) else 'FAIL'
    results.append((v.name, status, elapsed))
    marks = f"{n_pass} [PASS]" + (f", {n_fail_marker} [FAIL]"
                                  if n_fail_marker else "")
    print(f"[{v.name}] {status} in {elapsed:.1f}s ({marks}, "
          f"output: {os.path.basename(out_file)})")
    print()

total = time.time() - start
n_fail = sum(1 for _, s, _ in results if s == 'FAIL')
n_pass_v = sum(1 for _, s, _ in results if s == 'PASS')
n_skip = sum(1 for _, s, _ in results if s == 'SKIP')

print('=' * 64)
print(f"  Summary -- seed = {args.seed}, {mode}")
print('=' * 64)
for name, status, elapsed in sorted(results, key=lambda r: -r[2]):
    print(f"  {status:6s}  {name:22s}  {elapsed:7.1f}s")
print(f"  {'-' * 42}")
print(f"  {n_pass_v} passed, {n_fail} failed, {n_skip} skipped"
      f"{'':10s}{total:7.1f}s")
print('=' * 64)

sys.exit(1 if n_fail else 0)
