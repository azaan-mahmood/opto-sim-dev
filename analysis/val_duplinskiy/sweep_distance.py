"""Distance sweep for Duplinskiy et al. BB84 replication.

Generates QBER vs. distance curve matching paper Figure 6.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from src.protocols.bb84_duplinskiy import simulate_bb84_duplinskiy


def sweep_distance(distances, num_bits=100000, mu=0.1, seed=42):
    """Run QBER vs distance sweep."""
    results = []
    for d in distances:
        r = simulate_bb84_duplinskiy(
            num_bits=num_bits, fiber_length=d,
            mu=mu, seed=seed, verbose=False)
        results.append(r)
        print(f"  {d:6.1f} km  loss={r['total_loss_dB']:5.1f} dB  "
              f"sifted={r['n_sifted']:5d}  errors={r['n_errors']:4d}  "
              f"QBER={r['qber']*100:6.2f}%")
    return results


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Duplinskiy distance sweep')
    parser.add_argument('--bits', type=int, default=100000,
                        help='Pulses per distance point (default 100k)')
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()

    distances = [0, 10, 20, 30, 40, 50, 60, 75, 100]

    print("Duplinskiy et al. BB84 — Distance sweep")
    print(f"  mu=0.1, seed={args.seed}, bits={args.bits:,} per point")
    print(f"  {'Distance':>8}  {'Loss':>6}  {'Sifted':>7}  {'Errors':>6}  {'QBER':>7}")
    print("-" * 55)
    results = sweep_distance(distances, num_bits=args.bits, seed=args.seed)

    qber_pct = [r['qber'] * 100 for r in results]
    losses = [r['total_loss_dB'] for r in results]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Panel (a): QBER vs distance
    ax1.plot(distances, qber_pct, 'bo-', markersize=6, linewidth=1.5)
    ax1.set_xlabel('Fiber length (km)')
    ax1.set_ylabel('QBER (%)')
    ax1.set_title('Duplinskiy et al. BB84 -- QBER vs Distance\n(0 km: 2.6% vs paper 2%; longer distances include uncompensated birefringence)',
                  fontsize=10)
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(0, max(qber_pct) * 1.15)

    # Validation marker: 0 km back-to-back matches paper
    if 0 in distances:
        idx = distances.index(0)
        ax1.plot(0, qber_pct[idx], 'r*', markersize=15, zorder=5,
                 label=f'{qber_pct[idx]:.1f}% @ 0 km (paper: ~2%)')
        ax1.legend(fontsize=8)

    # Panel (b): Total loss vs distance
    ax2.plot(distances, losses, 'rs-', markersize=6, linewidth=1.5)
    ax2.set_xlabel('Fiber length (km)')
    ax2.set_ylabel('Total loss (dB)')
    ax2.set_title('Channel Loss vs Fiber Length')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    outdir = os.path.join(os.path.dirname(__file__), '..', '..', 'analysis', 'val_duplinskiy')
    os.makedirs(outdir, exist_ok=True)
    outpath = os.path.join(outdir, 'qber_vs_distance--seed42.png')
    plt.savefig(outpath, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"\nSaved: {outpath}")
