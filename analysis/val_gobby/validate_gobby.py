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
DCR = 15.0         # SPAD dark count rate (Hz) -- ID230 spec, see caveat below
GATE_WIDTH = 1e-9  # SPAD gate (s)
P_DARK = DCR * GATE_WIDTH  # dark count probability per gate (legacy export)

# NOTE on DCR.  15 Hz is the ID230 datasheet figure -- a 2020-era detector
# being used to replicate a 2004 experiment.  At 15 Hz the dark
# contribution is negligible at *every* distance in this sweep
# (P_dark/p_signal = 2.5e-4 even at 122 km), so the simulated QBER comes
# out flat and cannot reproduce Gobby's rise from 3.3% to 8.9%.  Fitting
# the 122 km endpoint against the Monte Carlo gives DCR ~ 1,788 Hz, which
# reproduces the whole sweep to a mean residual of 0.36 pp.  That is ~119x
# this spec and should be justified against Gobby's actual detector
# performance before being used in the manuscript.  Pass --dcr to override.
# OPEN-1: decoder visibility actually passed to the Monte Carlo (was a
# silent default of 1.0).  V = 0.934 corresponds to Gobby's stated 3.3%
# short-range QBER floor via e_opt = (1 - V)/2 [1]; verified in the 2nd
# pass: V = 0.934, no noise -> 3.2961% at 1.5M pulses.
VISIBILITY = 0.934

# --- OPEN-2: statistical power guards ------------------------------------
# A flat pulse budget is badly matched to this sweep: the sifted fraction
# runs from ~2.3e-3 at 0 km to ~1.4e-5 at 122 km (a factor of 164), so a
# flat count over-samples the short end and starves the long end -- which
# is the only end anyone scrutinises.  Flat 10M leaves 122 km with ~140
# sifted bits (sigma = 2.34 pp) against a 3.78 pp effect: 1.6 sigma, i.e.
# unable to distinguish the V = 0.934 prediction (8.33 %) from the
# perfect-decoder control (4.55 %).
#
# MIN_SIFTED is the 3-sigma separation threshold for that comparison
# (~500 sifted bits).  Rows below it cannot support a claim, so the table
# writer refuses to emit them unless explicitly overridden.
# MIN_SIFTED sits below TARGET_SIFTED_DEFAULT on purpose -- see the same
# note in val_system_scenarios.py.  Equal values make the run a knife-edge
# where any undershoot discards the whole sweep at the final write.
MIN_SIFTED = 500
TARGET_SIFTED_DEFAULT = 3000      # sigma ~ 0.5 pp at QBER ~8 %
PILOT_BITS = 200_000              # cheap probe to measure the sifted rate
CEILING_DEFAULT = 500_000_000     # hard cap so a bad estimate cannot run away


def model_qber(dist_km, dcr=DCR, visibility=VISIBILITY, afterpulse_prob=0.05):
    """First-order analytic QBER for the time-bin chain, in percent.

    Every click channel is weighted into a single ratio rather than added
    as a standalone offset:

        c_sig  = mu * eta * T_link            signal click probability
        c_dark = 2 * DCR * t_gate             BOTH detectors can dark-count
        c_prim = c_sig + c_dark               primary clicks
        c_tot  = c_prim * (1 + a)             plus afterpulses

        errors = c_sig*(1-V)/2                wrong port, finite visibility
               + c_dark/2                     dark count in the wrong detector
               + a*c_prim/2                   afterpulses land at random

        QBER   = errors / c_tot

    Correctness of the limits: as c_sig -> 0 this tends to 1/2 (a
    dark-count-dominated link is a coin flip), and as c_dark -> 0 it tends
    to [(1-V)/2 + a/2]/(1+a), the misalignment-plus-afterpulse floor.

    ------------------------------------------------------------------
    THIS IS AN APPROXIMATION.  THE MONTE CARLO IS AUTHORITATIVE.
    ------------------------------------------------------------------
    Do NOT fit detector parameters against this function.  The 8th pass
    did exactly that and it produced DCR = 7,593 Hz, which made the Monte
    Carlo overshoot to 17.5 % at 122 km against Gobby's 8.9 %.  Fitting
    the same endpoint against the MC directly gave 1,788 Hz and reproduced
    the sweep to a mean residual of 0.36 pp.

    Accuracy against the MC (V = 0.9792, seed 42, >=2500 sifted/point):

        DCR = 1,788 Hz   mean |residual| 0.59 pp,  max 3.03 pp (122 km)
        DCR = 7,593 Hz   mean |residual| 1.22 pp,  max 4.73 pp (122 km)

    The superseded form -- `(P_dark/2)/(p_signal + P_dark) + QBER_opt` --
    scored 2.86 pp and 4.34 pp on the same data.  It had three defects:
    it counted only one detector's dark rate, halved the dark error term
    that should be whole, and bolted the misalignment on as an unweighted
    additive constant instead of a share of the clicks.

    Residual error is concentrated at the dark-dominated end, where dead
    time, afterpulse chaining off dark clicks, and double-click tie-breaks
    all interact in ways this closed form does not attempt to capture.
    Treat it as a guide curve for figures, not as a predictive model.
    """
    T_link = 10.0 ** (-ALPHA_dB * dist_km / 10.0)
    c_sig = MU * ETA * T_link          # small-signal approximation
    c_dark = 2.0 * dcr * GATE_WIDTH    # two detectors
    c_prim = c_sig + c_dark

    errors = (c_sig * (1.0 - visibility) / 2.0
              + c_dark / 2.0
              + afterpulse_prob * c_prim / 2.0)
    c_tot = c_prim * (1.0 + afterpulse_prob)
    return errors / c_tot * 100.0


def gobby_measured_qber(dist_km):
    """Gobby et al. (2004) measured QBER, interpolated from Fig 3 data.

    Returns the published value exactly at the four measured distances
    (4.4, 65, 101, 122 km); clamped at the endpoints elsewhere.
    """
    return np.interp(dist_km, GOBBY_DIST_KM, GOBBY_QBER)


def simulate_qber(dist_km, num_bits, seed=42, verbose=False,
                  visibility=VISIBILITY, dcr=DCR, afterpulse_prob=0.05):
    """Run Monte Carlo at a given distance.

    `dcr` and `afterpulse_prob` are exposed (rather than pinned to the
    ID230 module constants) because the detector parameters are a *fitted
    input* to this replication, not a given: ID230 is a 2020-era part and
    Gobby is a 2004 experiment.  See the DCR discussion in OPEN-2.
    """
    results = simulate_bb84_time_bin(
        num_bits=num_bits,
        fiber_length=dist_km,
        mu=MU, wavelength=LAM,
        repetition_rate=2.5e6, pulse_width=100e-12,
        spad_eta=ETA, dark_count_rate=dcr,
        afterpulse_prob=afterpulse_prob, dead_time=13e-6,
        gate_width=GATE_WIDTH,
        visibility=visibility,      # OPEN-1: e_opt = (1-V)/2
        phase_error=0.0,
        seed=seed, verbose=verbose,
    )
    return results


def qber_err_pp(qber_frac, n_sifted):
    """Binomial s.d. of a QBER estimate, in percentage points."""
    if n_sifted <= 0:
        return float('inf')
    return np.sqrt(qber_frac * (1.0 - qber_frac) / n_sifted) * 100.0


def run_to_target(dist_km, target_sifted, ceiling, seed, visibility,
                  pilot_bits=PILOT_BITS, dcr=DCR):
    """Grow the pulse count at one distance until `target_sifted` sifted
    bits are collected, or `ceiling` pulses are spent.

    Runs a cheap pilot to measure this distance's sifted fraction, then
    scales to the required count and re-runs once from the same seed (so
    the result stays reproducible from the reported pulse count alone --
    it is a single deterministic run, not an accumulation).

    Returns (pulses_used, results_dict).
    """
    pilot_n = int(min(pilot_bits, ceiling))
    r = simulate_qber(dist_km, pilot_n, seed=seed, visibility=visibility,
                      dcr=dcr)

    n = pilot_n
    for _ in range(3):
        if r['n_sifted'] >= target_sifted or n >= ceiling:
            break
        if r['n_sifted'] == 0:
            n_next = ceiling    # no signal in the pilot: spend the budget
        else:
            # Re-estimate from the largest run so far. Headroom covers the
            # sampling error on the sifted count (~1/sqrt(n_sifted)); a
            # flat 15% undershoots about half the time, and an undershoot
            # that trips the write guard throws the whole sweep away.
            frac = r['n_sifted'] / n
            headroom = 1.15 + 3.0 / np.sqrt(max(r['n_sifted'], 1))
            n_next = int(np.ceil(target_sifted / frac * headroom))
        n_next = int(min(max(n_next, n + 1), ceiling))
        if n_next <= n:
            break
        n = n_next
        r = simulate_qber(dist_km, n, seed=seed, visibility=visibility,
                          dcr=dcr)

    return n, r


def check_statistical_power(rows, min_sifted, allow_underpowered):
    """Refuse to emit a table whose rows cannot support a claim (OPEN-2).

    `rows` is a sequence of (distance_km, n_sifted).  Raises unless every
    row clears `min_sifted`, or the caller explicitly opts out.
    """
    weak = [(d, s) for d, s, in rows if s < min_sifted]
    if not weak:
        return
    detail = ', '.join(f"{d:g} km: {s} sifted" for d, s in weak)
    if allow_underpowered:
        print(f"\n  !! WARNING: {len(weak)} row(s) below {min_sifted} sifted "
              f"bits ({detail}).")
        print("     Emitting anyway because --allow-underpowered was given.")
        print("     This table is a smoke run. Do NOT cite it as a result.")
        return
    raise RuntimeError(
        f"Refusing to write an under-powered table: {len(weak)} row(s) have "
        f"fewer than {min_sifted} sifted bits ({detail}).\n"
        f"A row below {min_sifted} cannot separate the V=0.934 prediction "
        f"from the perfect-decoder control at 3 sigma, so the artifact would "
        f"look like a result without being one (see OPEN-2 in "
        f"opto-sim-issues-and-fixes.md).\n"
        f"Fix: re-run with --target-sifted {TARGET_SIFTED_DEFAULT} (grows the "
        f"pulse count per distance), or raise --bits.\n"
        f"To emit anyway for a smoke test, pass --allow-underpowered."
    )


def run_validation(num_bits=200000, seed=42, distances=None,
                   visibility=VISIBILITY, target_sifted=None,
                   ceiling=CEILING_DEFAULT, allow_underpowered=False,
                   min_sifted=MIN_SIFTED, dcr=DCR):
    if distances is None:
        distances = [0, 4, 10, 20, 40, 65, 80, 100, 122]
    print("=" * 60)
    print("Gobby et al. 2004 — QBER vs Distance Validation")
    print(f"  Alpha: {ALPHA_dB} dB/km, Mu: {MU}, Eta: {ETA}")
    if target_sifted:
        print(f"  Budget: --target-sifted {target_sifted} per point "
              f"(ceiling {ceiling:,} pulses)")
    else:
        print(f"  Budget: flat {num_bits:,} pulses per point")
    print(f"  Seed: {seed}, Visibility: {visibility}")
    print("=" * 60)

    sim_dist = []
    sim_qber = []
    sim_sifted = []
    sim_pulses = []

    for d in distances:
        print(f"\n--- Distance: {d} km ---")
        if target_sifted:
            used, r = run_to_target(d, target_sifted, ceiling, seed,
                                    visibility, dcr=dcr)
        else:
            used, r = num_bits, simulate_qber(d, num_bits, seed=seed,
                                              verbose=True,
                                              visibility=visibility, dcr=dcr)
        sim_dist.append(d)
        sim_qber.append(r['qber'] * 100.0)
        sim_sifted.append(r['n_sifted'])
        sim_pulses.append(used)
        if r['n_sifted'] > 0:
            sd = qber_err_pp(r['qber'], r['n_sifted'])
            print(f"  QBER: {r['qber']*100:.2f} +/- {sd:.2f}%  "
                  f"({used:,} pulses, {r['n_sifted']} sifted)")
            print(f"        (Gobby measured: {gobby_measured_qber(d):.1f}%, "
                  f"this work analytic: {model_qber(d, dcr, visibility):.1f}%)")
        else:
            print(f"  No sifted bits in {used:,} pulses.")

    sim_dist = np.array(sim_dist)
    sim_qber = np.array(sim_qber)
    sim_sifted = np.array(sim_sifted)
    sim_pulses = np.array(sim_pulses)

    # OPEN-2 guard: refuse to emit an artifact that looks like a result
    # but cannot support one.  Checked before any file is written.
    check_statistical_power(list(zip(sim_dist, sim_sifted)), min_sifted,
                            allow_underpowered)

    # Analytical curve
    dense_dist = np.linspace(0, 130, 131)
    analytic = model_qber(dense_dist, dcr, visibility)

    # Save results table
    table_path = os.path.join(os.path.dirname(__file__), 'val_gobby_table.tex')
    with open(table_path, 'w') as f:
        f.write(f"% Generated by analysis/val_gobby/validate_gobby.py "
                f"--seed {seed} --visibility {visibility}"
                + (f" --target-sifted {target_sifted}\n" if target_sifted
                   else f" --bits {num_bits}\n"))
        f.write("% QBER error bars are binomial sqrt(q(1-q)/n_sifted).\n")
        f.write("%\n% !! DO NOT PUBLISH THIS TABLE YET -- see GOBBY-1 (section 18)\n"
                "%    in opto-sim-issues-and-fixes.md.  The chain is currently\n"
                "%    mis-parameterised against the source paper: alpha should be\n"
                "%    0.2 dB/km (not 0.182), eta should be eta_Bob = 0.045\n"
                "%    (including Bob's 5 dB apparatus loss, not the bare detector\n"
                "%    efficiency), the gate is 3.5 ns (not 1 ns), and the error\n"
                "%    count rate should come from Gobby's measured\n"
                "%    P_e = 8.5e-7/clock -- of which the majority is clock-laser\n"
                "%    stray light, which this chain does not model.  Visibility is\n"
                "%    also an OUTPUT, V = S/(S + 2*P_e), not an input.\n"
                "%    Regenerate once GOBBY-1 section 18.5 is implemented.\n%\n")
        f.write(r"\begin{tabular}{rccccc}" + "\n")
        f.write(r"  Distance & Pulses & Sifted & QBER & This work & Gobby et al. \\" + "\n")
        f.write(r"  (km)     &        & bits   & (\%) & analytic (\%) & measured (\%) \\" + "\n")
        f.write(r"\hline" + "\n")
        for d, q, s, n in zip(sim_dist, sim_qber, sim_sifted, sim_pulses):
            g_model = model_qber(d, dcr, visibility)
            g_meas = gobby_measured_qber(d)
            sd = qber_err_pp(q / 100.0, s)
            f.write(f"  {d:.0f} & {n:,} & {s} & ${q:.2f} \\pm {sd:.2f}$ "
                    f"& {g_model:.1f} & {g_meas:.1f} \\\\\n")
        f.write(r"\end{tabular}" + "\n")
    print(f"\nTable saved to {table_path}")

    # CSV alongside the .tex so the numbers are machine-readable
    csv_path = os.path.join(os.path.dirname(__file__),
                            f'val_gobby--seed{seed}.csv')
    with open(csv_path, 'w') as f:
        f.write(f"# validate_gobby.py seed={seed} visibility={visibility} "
                f"{'target_sifted=' + str(target_sifted) if target_sifted else 'bits=' + str(num_bits)}\n")
        f.write("distance_km,pulses,sifted_bits,qber_pct,qber_err_pp,"
                "analytic_pct,gobby_measured_pct\n")
        for d, q, s, n in zip(sim_dist, sim_qber, sim_sifted, sim_pulses):
            f.write(f"{d:g},{n},{s},{q:.4f},{qber_err_pp(q/100.0, s):.4f},"
                    f"{model_qber(d, dcr, visibility):.4f},"
                    f"{gobby_measured_qber(d):.4f}\n")
    print(f"CSV saved to {csv_path}")

    # Plot
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(dense_dist, analytic, 'b-', label=f'This work, analytic (V={visibility:g}, DCR={dcr:g} Hz)')
        mc_label = (f'Monte Carlo (>={target_sifted} sifted/pt)'
                    if target_sifted else f'Monte Carlo ({num_bits:,} pulses)')
        yerr = [qber_err_pp(q / 100.0, s)
                for q, s in zip(sim_qber, sim_sifted)]
        ax.errorbar(sim_dist, sim_qber, yerr=yerr, fmt='ro-', capsize=3,
                    label=mc_label)
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
    parser.add_argument('--bits', type=int, default=200000,
                        help='Flat pulses per point. Ignored when '
                             '--target-sifted is given.')
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--visibility', type=float, default=VISIBILITY,
                        help='Decoder visibility (default 0.934 = Gobby '
                             '3.3%% floor; use 1.0 for the perfect-decoder '
                             'control)')
    parser.add_argument('--target-sifted', type=int, default=None,
                        metavar='N',
                        help=f'Grow the pulse count per distance until N '
                             f'sifted bits are collected, instead of a flat '
                             f'--bits budget. Recommended: '
                             f'{TARGET_SIFTED_DEFAULT} (sigma ~0.5 pp). A '
                             f'flat budget starves the long distances, which '
                             f'are the contested ones -- see OPEN-2.')
    parser.add_argument('--dcr', type=float, default=DCR,
                        help=f'SPAD dark count rate in Hz (default {DCR:g}, '
                             f'the ID230 spec). At the spec value the dark '
                             f'contribution is negligible at every distance '
                             f'and the sweep comes out flat; ~1788 Hz '
                             f'reproduces Gobby (MC-fitted at 122 km).')
    parser.add_argument('--ceiling', type=int, default=CEILING_DEFAULT,
                        help='Hard cap on pulses per point under '
                             '--target-sifted (default 500M)')
    parser.add_argument('--min-sifted', type=int, default=MIN_SIFTED,
                        help=f'Refuse to write tables if any row has fewer '
                             f'sifted bits than this (default {MIN_SIFTED})')
    parser.add_argument('--allow-underpowered', action='store_true',
                        help='Emit tables even when rows fall below '
                             '--min-sifted. For smoke runs only; the '
                             'artifact must not be cited.')
    args = parser.parse_args()
    try:
        run_validation(num_bits=args.bits, seed=args.seed,
                       visibility=args.visibility,
                       target_sifted=args.target_sifted,
                       ceiling=args.ceiling,
                       min_sifted=args.min_sifted,
                       allow_underpowered=args.allow_underpowered,
                       dcr=args.dcr)
    except RuntimeError as exc:
        print(f"\nERROR: {exc}", file=sys.stderr)
        sys.exit(1)
