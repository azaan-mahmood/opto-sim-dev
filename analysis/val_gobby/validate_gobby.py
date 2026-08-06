"""Replicate Gobby, Yuan & Shields (2004) — QBER vs distance.

References
----------
[1] Gobby, C., Yuan, Z. L., & Shields, A. J. (2004). Quantum key
    distribution over 122 km of standard telecom fiber. Appl. Phys.
    Lett. 84(19), 3762-3764.
"""
import argparse
import numpy as np
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from src.protocols.bb84_time_bin import simulate_bb84_time_bin

# Gobby paper data points (from Fig 3)
GOBBY_DIST_KM = np.array([4.4, 65.0, 101.0, 122.0])
GOBBY_QBER = np.array([3.3, 3.3, 6.0, 8.9])  # percent

# Gobby parameters
LAM = 1550e-9
ALPHA_dB = 0.182  # SMF-28 at 1550 nm
MU = 0.1
ETA = 0.10
DCR = 15.0         # SPAD dark count rate (Hz)
GATE_WIDTH = 1e-9  # SPAD gate (s)
P_DARK = DCR * GATE_WIDTH  # dark count probability per gate
QBER_OPT = 0.033    # optical misalignment baseline
# OPEN-1: decoder visibility actually passed to the Monte Carlo (was a
# silent default of 1.0).  V = 0.934 corresponds to Gobby's stated 3.3%
# short-range QBER floor via e_opt = (1 - V)/2 [1]; verified in the 2nd
# pass: V = 0.934, no noise -> 3.2961% at 1.5M pulses.
VISIBILITY = 0.934


def model_qber(dist_km):
    """QBER from this work's analytic model.

    QBER = (P_dark / 2) / (mu * eta * T_link + P_dark) + QBER_opt

    where T_link = 10^(-alpha * L / 10).  The factor 1/2 accounts for the
    fact that a dark count lands in the correct detector half the time.

    Returns QBER in percent.
    """
    T_link = 10.0 ** (-ALPHA_dB * dist_km / 10.0)
    p_signal = MU * ETA * T_link  # small-signal approximation
    qber = (P_DARK / 2.0) / (p_signal + P_DARK) + QBER_OPT
    return qber * 100.0


def gobby_measured_qber(dist_km):
    """Gobby et al. (2004) measured QBER, interpolated from Fig 3 data.

    Returns the published value exactly at the four measured distances
    (4.4, 65, 101, 122 km); clamped at the endpoints elsewhere.
    """
    return np.interp(dist_km, GOBBY_DIST_KM, GOBBY_QBER)


def simulate_qber(dist_km, num_bits, seed=42, verbose=False,
                  visibility=VISIBILITY):
    """Run Monte Carlo at a given distance."""
    results = simulate_bb84_time_bin(
        num_bits=num_bits,
        fiber_length=dist_km,
        mu=MU, wavelength=LAM,
        repetition_rate=2.5e6, pulse_width=100e-12,
        spad_eta=ETA, dark_count_rate=DCR,
        afterpulse_prob=0.05, dead_time=13e-6,
        gate_width=GATE_WIDTH,
        visibility=visibility,      # OPEN-1: e_opt = (1-V)/2 = 3.3%
        phase_error=0.0,
        seed=seed, verbose=verbose,
    )
    return results


def run_validation(num_bits=200000, seed=42, distances=None,
                   visibility=VISIBILITY):
    if distances is None:
        distances = [0, 4, 10, 20, 40, 65, 80, 100, 122]
    print("=" * 60)
    print("Gobby et al. 2004 — QBER vs Distance Validation")
    print(f"  Alpha: {ALPHA_dB} dB/km, Mu: {MU}, Eta: {ETA}")
    print(f"  Pulses per point: {num_bits}, Seed: {seed}, "
          f"Visibility: {visibility}")
    print("=" * 60)

    sim_dist = []
    sim_qber = []
    sim_sifted = []

    for d in distances:
        print(f"\n--- Distance: {d} km ---")
        r = simulate_qber(d, num_bits, seed=seed, verbose=True,
                          visibility=visibility)
        sim_dist.append(d)
        sim_qber.append(r['qber'] * 100.0)
        sim_sifted.append(r['n_sifted'])
        if r['n_sifted'] > 0:
            print(f"  QBER: {r['qber']*100:.2f}%  (Gobby measured: {gobby_measured_qber(d):.1f}%, "
                  f"this work analytic: {model_qber(d):.1f}%)")

    sim_dist = np.array(sim_dist)
    sim_qber = np.array(sim_qber)
    sim_sifted = np.array(sim_sifted)

    # Analytical curve
    dense_dist = np.linspace(0, 130, 131)
    analytic = model_qber(dense_dist)

    # Save results table
    table_path = os.path.join(os.path.dirname(__file__), 'val_gobby_table.tex')
    with open(table_path, 'w') as f:
        f.write(r"\begin{tabular}{rccccc}" + "\n")
        f.write(r"  Distance & Pulses & Sifted & QBER & This work & Gobby et al. \\" + "\n")
        f.write(r"  (km)     &        & bits   & (\%) & analytic (\%) & measured (\%) \\" + "\n")
        f.write(r"\hline" + "\n")
        for d, q, s in zip(sim_dist, sim_qber, sim_sifted):
            g_model = model_qber(d)
            g_meas = gobby_measured_qber(d)
            f.write(f"  {d:.0f} & {num_bits} & {s} & {q:.2f} & {g_model:.1f} & {g_meas:.1f} \\\\\n")
        f.write(r"\end{tabular}" + "\n")
    print(f"\nTable saved to {table_path}")

    # Plot
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(dense_dist, analytic, 'b-', label=f'This work, analytic (QBER_opt={QBER_OPT*100:.0f}%)')
        ax.plot(sim_dist, sim_qber, 'ro-', label=f'Monte Carlo ({num_bits} pulses)')
        ax.plot(GOBBY_DIST_KM, GOBBY_QBER, 'gs', label='Gobby paper (Fig 3)')
        ax.set_xlabel('Distance (km)')
        ax.set_ylabel('QBER (%)')
        ax.set_title('Time-bin BB84 — QBER vs Distance')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.axhline(y=11.0, color='r', linestyle='--', alpha=0.5, label='BB84 threshold (11%)')
        ax.set_xlim(0, 130)
        ax.set_ylim(0, 50)

        png_path = os.path.join(os.path.dirname(__file__), f'val_gobby--seed{seed}.png')
        fig.savefig(png_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"Figure saved to {png_path}")
    except ImportError:
        print("matplotlib not available — skipping figure")

    print("\nDone.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Gobby et al. 2004 QBER validation")
    parser.add_argument('--bits', type=int, default=200000)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--visibility', type=float, default=VISIBILITY,
                        help='Decoder visibility (default 0.934 = Gobby '
                             '3.3%% floor; use 1.0 for the perfect-decoder '
                             'control)')
    args = parser.parse_args()
    run_validation(num_bits=args.bits, seed=args.seed,
                   visibility=args.visibility)
