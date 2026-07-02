#!/usr/bin/env python3
"""
Run all Tier-1 channel validation scripts with a single seed.

Usage:
    python run_all.py                     # seed = 42
    python run_all.py --seed 123          # custom seed
"""

import subprocess, sys, argparse, os, time

BASE = os.path.dirname(os.path.abspath(__file__))

parser = argparse.ArgumentParser()
parser.add_argument('--seed', type=int, default=42)
args = parser.parse_args()

scripts = [
    ('CD',          os.path.join(BASE, 'analysis', 'validation', 'validate_cd.py')),
    ('PMD',         os.path.join(BASE, 'analysis', 'validation', 'validate_pmd.py')),
    ('Attenuation', os.path.join(BASE, 'analysis', 'validation', 'validate_attenuation.py')),
    ('Birefringence', os.path.join(BASE, 'analysis', 'validation', 'validate_birefringence.py')),
]

print(f"{'='*60}")
print(f"  Tier-1 Channel Validation -- seed = {args.seed}")
print(f"{'='*60}\n")

start = time.time()
results = []

for name, path in scripts:
    t0 = time.time()
    print(f"[{name}] {'-'*50}")
    ret = subprocess.run(
        [sys.executable, path, '--seed', str(args.seed)],
        capture_output=False, text=True
    )
    elapsed = time.time() - t0
    ok = ret.returncode == 0
    results.append((name, ok, elapsed))
    if not ok:
        print(f"[{name}] FAILED (rc={ret.returncode})")
    print()

total = time.time() - start

print(f"{'='*60}")
print(f"  Summary -- seed = {args.seed}")
print(f"{'='*60}")
for name, ok, elapsed in results:
    status = 'PASS' if ok else 'FAIL'
    print(f"  {status:6s}  {name:15s}  {elapsed:.1f}s")
print(f"  {'-'*30}")
print(f"  TOTAL                    {total:.1f}s")
print(f"{'='*60}")
