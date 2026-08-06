"""
System-Level Impairment Scenarios (BLOCK-2: Table 11 generator)
===============================================================
Generates the system-level impairment table for the time-bin BB84 chain:
each row is an *explicit, self-contained* impairment configuration, run
at a fixed distance with a recorded seed, and written to both a CSV and a
LaTeX table.  Every row's exact config dict is printed and stored, and the
LaTeX caption carries the generating script, seed, pulse count, and commit
hash — the reproducibility contract BLOCK-2 demands (the old
`paperwork/tables/val_system_table.tex` was hand-written: no script, no
seed, no CSV).

Note: `paperwork/` was deleted in an earlier session, so the table is
written next to the other `val_system/` artifacts; the manuscript author
pastes it back into `paperwork/tables/` when the paper is restored.

Chain (identical to `analysis/val_system.py`, which this imports):
  CWLaser(1 MHz, Y-pol) -> MZM carve (30 ps) -> encoder AMZI
    -> FiberRealization(birefringence + CD + PMD + attenuation)
    -> decoder AMZI (visibility V) -> 2x SPAD (ID230), Gobby defaults.

Physics to expect (established by ARCH-1, panels A-E):
  - time-bin encoding is immune to slow birefringence: rows that only
    toggle birefringence/CD/PMD move the sifted *rate*, not the QBER;
  - QBER is moved by decoder visibility V (e_opt = (1-V)/2) and by dark
    counts relative to the attenuated signal.
At 100 km (22 dB) with mu = 0.1 the attenuation-on rows are
sample-limited at a few M pulses; the table reports sifted counts and the
binomial error bar so that is explicit, and `--bits`/`--distance` let the
author regenerate at any budget.

Usage
-----
    python analysis/val_system_scenarios.py            # defaults below
    python analysis/val_system_scenarios.py --bits 50000000 --distance 75
"""
import argparse
import os
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from analysis.val_system import simulate_point, qber_err, VIS_DEFAULT

# --- Explicit impairment configurations (one full dict per row) --------
# Every row is self-contained: no hidden defaults.  Values that differ
# from the chain defaults are spelled out; `distance_km` is injected at
# runtime from --distance.
SCENARIOS = [
    dict(
        name='No impairments',
        config=dict(dispersion=False, birefringence=False, attenuation=False,
                    visibility=1.0, dcr=0.0,
                    note='floor: phase jitter + afterpulse + double-click'),
    ),
    dict(
        name='Attenuation only',
        config=dict(dispersion=False, birefringence=False, attenuation=True,
                    visibility=1.0, dcr=0.0,
                    note='Keiser 10^(-alpha L/10); no other impairment'),
    ),
    dict(
        name='+ CD',
        config=dict(dispersion=False, cd=True, pmd=False,
                    birefringence=False, attenuation=True,
                    visibility=1.0, dcr=0.0,
                    note='Agrawal H(Omega) = exp(-j beta2 Omega^2 L / 2)'),
    ),
    dict(
        name='+ PMD',
        config=dict(dispersion=False, cd=False, pmd=True,
                    birefringence=False, attenuation=True,
                    visibility=1.0, dcr=0.0,
                    note='Maxwellian DGD, Razavi'),
    ),
    dict(
        name='+ Birefringence',
        config=dict(dispersion=False, birefringence=True, attenuation=True,
                    visibility=1.0, dcr=0.0,
                    note='time-bin immune to quasi-static SU(2) rotation'),
    ),
    dict(
        name='+ Visibility V=0.934',
        config=dict(dispersion=False, birefringence=True, attenuation=True,
                    visibility=VIS_DEFAULT, dcr=0.0,
                    note='Gobby floor e_opt = (1-V)/2 = 3.3 %'),
    ),
    dict(
        name='+ Dark counts 10 kHz',
        config=dict(dispersion=False, birefringence=True, attenuation=True,
                    visibility=1.0, dcr=10000.0,
                    note='detector-noise floor (ID230 spec is 15 Hz)'),
    ),
    dict(
        name='Full chain',
        config=dict(dispersion=True, birefringence=True, attenuation=True,
                    visibility=VIS_DEFAULT, dcr=15.0,
                    note='Gobby-style defaults: CD+PMD+biref+V=0.934'),
    ),
]

BASE = dict(pulse_sigma=30e-12, mu=0.1)


def commit_hash():
    try:
        return subprocess.run(
            ['git', 'rev-parse', '--short', 'HEAD'],
            capture_output=True, text=True, cwd=os.path.dirname(__file__)
        ).stdout.strip() or 'unknown'
    except Exception:
        return 'unknown'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--bits', type=int, default=2000000,
                        help='Pulses per scenario (default 2M; 28 us/pulse '
                             '=> ~2 min per 5M at 100 km)')
    parser.add_argument('--distance', type=float, default=100.0,
                        help='Fibre length in km (default 100)')
    args = parser.parse_args()
    SEED, NUM_BITS, DIST = args.seed, args.bits, args.distance
    COMMIT = commit_hash()

    OUT = os.path.join(os.path.dirname(__file__), '..', 'val_system')
    os.makedirs(OUT, exist_ok=True)

    print(f"Impairment scenarios: {len(SCENARIOS)} configs x {NUM_BITS} "
          f"pulses, {DIST} km, seed {SEED}", flush=True)

    rows = []
    for s in SCENARIOS:
        config = dict(BASE)
        config.update(s['config'])
        config['fiber_length'] = DIST
        config['seed'] = SEED
        sim_config = {k: v for k, v in config.items() if k != 'note'}
        print(f"\n[{s['name']}]", flush=True)
        print(f"  config = {config}", flush=True)
        res = simulate_point(num_bits=NUM_BITS, **sim_config)
        q = res['qber']
        err = qber_err(q, res['n_sifted'])
        print(f"  QBER = {q*100:6.2f} +/- {err*100:5.2f}%  "
              f"sifted = {res['n_sifted']} ({res['n_sifted']/NUM_BITS*100:.3f}%)"
              f"  errors = {res['n_errors']}", flush=True)
        rows.append((s['name'], config, res, err))

    # --- CSV ---
    csv_path = os.path.join(OUT, f'val_system_scenarios--seed{SEED}.csv')
    with open(csv_path, 'w') as f:
        f.write(f"# Impairment scenarios, time-bin BB84 chain "
                f"(BLOCK-2, replaces hand-written Table 11)\n")
        f.write(f"# script: analysis/val_system_scenarios.py  seed: {SEED}  "
                f"bits: {NUM_BITS}  distance: {DIST} km  commit: {COMMIT}\n")
        f.write(f"# chain: CWLaser -> MZM carve -> encoder AMZI -> "
                f"FiberRealization -> decoder AMZI -> 2x SPAD\n")
        f.write(f"# qber_err: binomial sqrt(q(1-q)/n_sifted)\n")
        f.write("Scenario,QBER_fraction,QBER_err,Sifted_bits,Sifted_fraction,"
                "Errors,Config\n")
        for name, config, res, err in rows:
            f.write(f"{name},{res['qber']:.6f},{err:.6f},{res['n_sifted']},"
                    f"{res['n_sifted']/NUM_BITS:.6f},{res['n_errors']},"
                    f"{config}\n")
    print(f"\nSaved: {csv_path}")

    # --- LaTeX table ---
    tex_path = os.path.join(OUT, f'val_system_scenarios--seed{SEED}.tex')
    with open(tex_path, 'w') as f:
        f.write(f"% Generated by analysis/val_system_scenarios.py "
                f"--seed {SEED} --bits {NUM_BITS} --distance {DIST}\n")
        f.write(f"% Commit: {COMMIT} -- regenerate to match a new commit.\n")
        f.write("\\begin{table}[ht]\n\\centering\n")
        f.write("\\caption{System-level impairment scenarios for the "
                f"time-bin BB84 chain at {DIST:g} km "
                f"({NUM_BITS:,} pulses per scenario, seed {SEED}; "
                f"generated by \\texttt{{analysis/val\\_system\\_scenarios.py}} "
                f"@ \\texttt{{{COMMIT}}}). "
                f"QBER error bars are the binomial "
                f"$\\sqrt{{q(1-q)/n_\\mathrm{{sifted}}}}$. "
                f"Scenarios that only toggle birefringence/CD/PMD move the "
                f"sifted rate, not the QBER: time-bin encoding is immune to "
                f"the quasi-static fibre rotation, and CD/PMD broadening "
                f"stays inside the 1 ns gate.}}\n")
        f.write("\\begin{tabular}{lrrrr}\n\\hline\n")
        f.write("Scenario & QBER (\\%) & Sifted & Sifted (\\%) & Config "
                "\\\\\n\\hline\n")
        for name, config, res, err in rows:
            cfg = '; '.join(f"{k}={v}" for k, v in sorted(config.items())
                            if k not in ('seed', 'name') and v is not True)
            cfg = cfg.replace('_', '\\_').replace('%', '\\%')
            f.write(f"{name} & ${res['qber']*100:.2f} \\pm {err*100:.2f}$ & "
                    f"{res['n_sifted']} & {res['n_sifted']/NUM_BITS*100:.3f} "
                    f"& {cfg} \\\\\n")
        f.write("\\hline\n\\end{tabular}\n\\end{table}\n")
    print(f"Saved: {tex_path}")


if __name__ == "__main__":
    main()
