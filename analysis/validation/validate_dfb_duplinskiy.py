"""The Duplinskiy polarisation chain driven by the DFB device model.

Until now `bb84_duplinskiy.py` injected a flat analytic field -- one time
sample of `sqrt(power_per_pulse/2) * ones((1,2))`, no laser at all.  This
runs the same chain from `src/lasers/dfblaser.py` instead and measures
what changes.

The wiring is the paper's own (Duplinskiy et al., Opt. Express 25(23),
28886, 2017, Fig. 1 and Sec. 2):

    "Alice produces linearly polarized optical pulses using a 1550 nm
     laser source.  The subsequent polarization controller (PC 1) is
     configured in such a way that the amplitudes of the field along the
     ordinary and extraordinary axes of the crystal inside the modulator
     (PM 1) are equal."

The DFB emits one TE mode, which is linearly polarised, and PC1 is a Jones
rotation to 45 degrees -- `polarization_azimuth = pi/4` on the driver.
The paper's 10 MHz repetition sits far above the device's ~210 ps floor.

What this establishes
---------------------
1. **Stokes vectors through every stage**, for all four BB84 states.
2. **A control**: the 8-outcome response table must be unchanged by the
   swap.  Measured worst relative difference 2.6e-15.  If that ever
   breaks, the QBER comparison is no longer controlled and means nothing.
3. **QBER** at 0/10/50 km for the flat field and for the DFB at both ends
   of its usable drive window.

Why the polarisation state cannot move
--------------------------------------
`LaserDriver.sample_field` returns one complex amplitude times a fixed
Jones vector, so **both components carry the same amplitude**.  Normalised
Stokes parameters depend only on the ratio Ey/Ex, so the source's RIN,
phase noise and chirp divide out exactly -- measured DOP = 1.000000 with
the states exact.  A polarisation-encoding chain is blind to everything a
real source adds *except pulse energy*.  This is the same cancellation
that removes linewidth in a matched AMZI, reached independently from the
other side.

So the prediction, on record before the run: **a null**, with any residual
coming only through per-pulse energy.

The pulse width is not in the paper
-----------------------------------
Sec. 2, 6 and 7 give the 10 MHz repetition, the 20 ns detection window,
the ID230 detector figures, 0.1 photons per pulse and 2 % QBER -- but
never the optical pulse duration.  It is load-bearing here, because it
sets the per-pulse energy spread, which is the one thing the DFB
contributes.  The device bounds it: below ~100 ps it barely lases (53 %
shot-to-shot), above ~250 ps the optical pulse follows the drive and it is
no longer gain switching.  Within that window the spread runs 1.1 % to
11.1 %, so the comparison is run at both ends.  If they agree, the
unstated width is registered but not load-bearing for the result.

References
----------
[1] Duplinskiy, Ustimchik, Kanapin, Kurochkin & Kurochkin, "Low loss QKD
    optical scheme for fast polarization encoding", Opt. Express 25(23),
    28886-28897 (2017).
[2] Kim, Chung & Lee, IEEE J. Quantum Electron. 36(7), 787-794 (2000).
[3] Collett, E., "Field Guide to Polarization", SPIE Press, 2005, Ch. 2.
"""
import argparse
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from src.channel import optics, FiberRealization
from src.channel.phase_modulator import PhaseModulator
from src.lasers import DFBLaser, DriveParams, LaserDriver
from src.protocols.bb84_duplinskiy import simulate_bb84_duplinskiy
from src.visualization import compute_stokes_parameters

OUT_DIR = os.path.join(os.path.dirname(__file__), '..', 'val_dfb')


def _stem(quick):
    """Smoke runs write to their own files.

    Sharing paths with the full run meant `--quick` silently replaced a
    quotable figure with an under-powered one, and nothing warned: the PNG
    simply became worse.  `.gitignore` excludes the `--quick` names, so a
    smoke artifact cannot reach the repository either.
    """
    return ('val_dfb_duplinskiy_poincare--quick' if quick
            else 'val_dfb_duplinskiy_poincare')

SEED = 42
LASER_SEED = 11
N_SECTIONS = 15

# Gain-switched drive.  60 -> 200 mA is not from the paper either, which
# describes no laser beyond "1550 nm laser source"; it is the operating
# point characterised in validate_dfb_drive.py.
I_BIAS, I_PEAK = 0.060, 0.140

# Harvest period.  Must be >= 5 ns: at 2 ns the device has not recovered
# between pulses and the mean energy comes out 17 % high, where 5 and
# 20 ns agree to 3 %.  The paper runs at 100 ns, which would cost ~11 s of
# wall time per pulse to simulate directly.
HARVEST_PERIOD = 5e-9
HARVEST_WINDOW = 150e-9
SETTLE = 40e-9

# The two ends of the usable gain-switching window (see the docstring).
WIDTH_NARROW = 100e-12    # ~11 % per-pulse energy spread
WIDTH_WIDE = 200e-12      # ~1 % per-pulse energy spread

STATES = (('X', 0, 'D'), ('X', 1, 'A'), ('C', 0, 'R'), ('C', 1, 'L'))


def dfb_source(width, azimuth=np.pi / 4, dt=5e-12, n=4096):
    """DFB field at PM1's input, and its per-pulse energy factors.

    The azimuth is the paper's PC1: it makes the amplitudes along the two
    crystal axes equal.  Returns (field, factors, diagnostics).
    """
    las = DFBLaser(n_sections=N_SECTIONS, seed=LASER_SEED)
    drive = DriveParams(mode='gain_switched', waveform='gaussian',
                        i_bias=I_BIAS, i_peak=I_PEAK,
                        period=HARVEST_PERIOD, width=width)
    drv = LaserDriver(las, drive, seed=LASER_SEED, polarization_azimuth=azimuth)
    field = drv.sample_field(dt, n)

    # Per-pulse energies, taken from the device timebase so the pulse is
    # resolved rather than decimated.
    res = drv.run(t_end=SETTLE + HARVEST_WINDOW, record_every=1)
    m = res.t >= SETTLE
    t, P = res.t[m] - SETTLE, res.P_right[m]
    energies = []
    for k in range(1, int(t[-1] // HARVEST_PERIOD)):
        w = (t >= k * HARVEST_PERIOD) & (t < (k + 1) * HARVEST_PERIOD)
        if w.sum() > 8:
            energies.append(P[w].sum() * las.dt)
    energies = np.array(energies)
    factors = energies / energies.mean()
    diag = dict(n_pulses=len(factors), mean_J=float(energies.mean()),
                spread=float(100 * factors.std()), dt_dev=las.dt)
    return field, factors, diag


def _fmt(S, dop, psi, chi):
    return (f"[{S[0]:.3f}, {S[1]:+.6f}, {S[2]:+.6f}, {S[3]:+.6f}]  "
            f"DOP={dop:.6f}  psi={math.degrees(psi):+7.2f} deg  "
            f"chi={math.degrees(chi):+7.2f} deg")


def stokes(E):
    S, (psi, chi) = compute_stokes_parameters(E)
    dop = float(np.sqrt(S[1] ** 2 + S[2] ** 2 + S[3] ** 2))
    return S, dop, psi, chi


def stokes_section(field, fiber_km, failures):
    """Track all four BB84 states through every stage. Returns plot data."""
    pm_a = PhaseModulator(crystal_cut='X', modulation='DC')
    pm_b = PhaseModulator(crystal_cut='X', modulation='DC')
    Vpi = pm_a.Vpi
    fibre = FiberRealization(L_m=fiber_km * 1000, temperature=25,
                             bend_radius=None, attenuation_factor=0.2,
                             cd=False, pmd=False, model='auto', seed=SEED)
    J = fibre.birefringence_matrix()
    U = None if J is None else J.conj().T

    print(f"\n  Stokes through the chain, {fiber_km} km "
          f"(states from the paper's Eqs. 4-5)")
    S, dop, psi, chi = stokes(field)
    print(f"    source, after PC1 (45 deg): {_fmt(S, dop, psi, chi)}")
    if abs(S[2] - 1.0) > 1e-9 or dop < 1 - 1e-9:
        failures.append("the field entering PM1 is not D with DOP 1; PC1 is "
                        "not delivering equal amplitudes on the crystal axes")

    stages = {'after PM1': [], 'after fibre': [], 'after PC2': [], 'at PBS': []}
    for basis, bit, name in STATES:
        v_a = (Vpi / 2 if bit == 0 else 3 * Vpi / 2) if basis == 'C' else \
              (0 if bit == 0 else Vpi)
        E1 = pm_a.modulate(E_field=field, V=v_a)
        E2 = fibre.apply(E1, dt=1e-7)
        E3 = np.transpose(U @ np.transpose(E2)) if U is not None else E2
        E4 = pm_b.modulate(E_field=optics.voa(E3, 2.0),
                           V=(0 if basis == 'C' else Vpi / 2))
        print(f"    {basis}{bit} = {name}")
        for lab, E in (('after PM1', E1), ('after fibre', E2),
                       ('after PC2', E3), ('at PBS', E4)):
            S, dop, psi, chi = stokes(E)
            stages[lab].append((name, S))
            print(f"      {lab:12s} {_fmt(S, dop, psi, chi)}")
            if dop < 1 - 1e-9:
                failures.append(f"{name} lost its degree of polarisation at "
                                f"'{lab}' (DOP={dop:.6f}); a unitary chain "
                                "cannot depolarise")

    # The four encoded states must be exact.
    want = {'D': (0, +1, 0), 'A': (0, -1, 0), 'R': (0, 0, +1), 'L': (0, 0, -1)}
    for name, S in stages['after PM1']:
        tgt = want[name]
        err = max(abs(S[1] - tgt[0]), abs(S[2] - tgt[1]), abs(S[3] - tgt[2]))
        if err > 1e-9:
            failures.append(f"state {name} after PM1 is off its target Stokes "
                            f"vector by {err:.2e}")
    return stages


def response_control(field, failures):
    """The deterministic chain must be unchanged by the source swap.

    Checked through the observable rather than by reaching into the
    protocol's table: with `pulse_energy_factors` left off, swapping the
    flat field for the DFB must give bit-identical sifted and error
    counts.  If it does not, the QBER comparison below is not a controlled
    experiment and nothing in it means anything.
    """
    print("\n  control: the source swap must change nothing deterministic")
    clean = True
    for L in (0, 10, 50):
        ra = simulate_bb84_duplinskiy(20000, fiber_length=L, seed=SEED)
        rb = simulate_bb84_duplinskiy(20000, fiber_length=L, seed=SEED,
                                      source_field=field)
        same = (ra['n_sifted'], ra['n_errors']) == (rb['n_sifted'], rb['n_errors'])
        clean &= same
        print(f"    {L:3d} km  flat {ra['n_sifted']}/{ra['n_errors']}   "
              f"DFB {rb['n_sifted']}/{rb['n_errors']}   "
              f"{'identical' if same else 'DIFFERENT'}")
    if clean:
        print("    -> only the per-pulse energy sequence can move the result.")
    else:
        failures.append("the DFB source changed the deterministic response; "
                        "the QBER comparison below is not controlled")


def qber_section(sources, distances, n_pulses, failures):
    print("\n  QBER, flat field against the DFB at both ends of its drive window")
    print("    source          km   sifted  errors     QBER +/- 1 sigma      vs flat")
    base = {}
    for label, field, factors in sources:
        for L in distances:
            n = n_pulses[L]
            r = simulate_bb84_duplinskiy(
                n, fiber_length=L, seed=SEED,
                source_field=field, pulse_energy_factors=factors)
            s, e = r['n_sifted'], r['n_errors']
            q = r['qber']
            sig = math.sqrt(max(q * (1 - q), 1e-12) / s) if s else float('nan')
            if label == 'flat':
                base[L] = (q, sig)
                delta = ''
            else:
                q0, s0 = base[L]
                d = q - q0
                ds = math.sqrt(sig ** 2 + s0 ** 2)
                delta = f"  {d * 100:+.3f} +/- {ds * 100:.3f} pp"
            print(f"    {label:14s} {L:3d} {s:8d} {e:7d}   "
                  f"{q * 100:6.3f} +/- {sig * 100:.3f} %{delta}")
            if s == 0:
                failures.append(f"{label} at {L} km produced no sifted bits")
    return base


def run(quick=False, do_qber=True):
    os.makedirs(OUT_DIR, exist_ok=True)
    failures = []
    print("=" * 72)
    print("Duplinskiy polarisation chain, driven by the DFB device model")
    print("=" * 72)

    print("\n  source: gain-switched DFB, "
          f"{I_BIAS * 1e3:.0f} -> {(I_BIAS + I_PEAK) * 1e3:.0f} mA, "
          f"PC1 at 45 deg")
    sources = []
    fields = {}
    for label, width in (('DFB 200ps', WIDTH_WIDE), ('DFB 100ps', WIDTH_NARROW)):
        field, factors, diag = dfb_source(width)
        fields[label] = field
        sources.append((label, field, factors))
        print(f"    {label}: {diag['n_pulses']} pulses harvested at "
              f"{HARVEST_PERIOD * 1e9:.0f} ns, mean {diag['mean_J']:.3e} J, "
              f"energy spread {diag['spread']:.2f} %")

    stages = stokes_section(fields['DFB 200ps'], 50, failures)
    response_control(fields['DFB 200ps'], failures)

    if do_qber:
        distances = (0, 10) if quick else (0, 10, 50)
        # Sized for >= 3000 sifted bits per cell, the standard used
        # here for a quotable QBER.  Measured yields: ~2.4e-3 sifted
        # per pulse at 0 km, ~1.7e-3 at 10 km, ~4e-4 at 50 km.
        n_pulses = {0: 200_000 if quick else 1_400_000,
                    10: 300_000 if quick else 2_000_000,
                    50: 400_000 if quick else 10_000_000}
        all_sources = [('flat', None, None)] + sources
        qber_section(all_sources, distances, n_pulses, failures)

    _figure(stages, quick)

    print()
    if failures:
        print("[FAIL]")
        for f in failures:
            print(f"  - {f}")
        return 1
    print("[PASS] the four BB84 states are exact with DOP 1, and the source")
    print("       swap changes nothing deterministic in the chain")
    return 0


def _figure(stages, quick=False):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("  matplotlib missing, skipping the figure")
        return

    fig = plt.figure(figsize=(13, 4.6))
    order = ['after PM1', 'after fibre', 'after PC2']
    titles = ['encoded at Alice (PM1)', 'after 50 km of fibre',
              "after Bob's compensation (PC2)"]
    colours = {'D': 'tab:blue', 'A': 'tab:orange',
               'R': 'tab:green', 'L': 'tab:red'}

    for k, (lab, title) in enumerate(zip(order, titles)):
        ax = fig.add_subplot(1, 3, k + 1, projection='3d')
        u, v = np.mgrid[0:2 * np.pi:60j, 0:np.pi:30j]
        ax.plot_wireframe(np.cos(u) * np.sin(v), np.sin(u) * np.sin(v),
                          np.cos(v), color='0.85', linewidth=0.4)
        for name, S in stages[lab]:
            ax.scatter(S[1], S[2], S[3], s=55, color=colours[name],
                       depthshade=False, label=name)
            ax.plot([0, S[1]], [0, S[2]], [0, S[3]], color=colours[name],
                    lw=1.0, alpha=0.6)
        ax.set_xlabel('S1'); ax.set_ylabel('S2'); ax.set_zlabel('S3')
        ax.set_xlim(-1, 1); ax.set_ylim(-1, 1); ax.set_zlim(-1, 1)
        ax.set_title(title, fontsize=10)
        if k == 0:
            ax.legend(fontsize=8, loc='upper left')

    fig.suptitle('BB84 polarisation states on the Poincare sphere, '
                 'DFB-driven Duplinskiy chain', fontsize=11)
    fig.tight_layout()
    png = os.path.join(OUT_DIR, _stem(quick) + '.png')
    fig.savefig(png, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"\n  figure: {png}")


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--quick', action='store_true',
                    help='fewer pulses and two distances, for a smoke run')
    ap.add_argument('--no-qber', action='store_true',
                    help='Stokes and the control only; skip the QBER runs')
    a = ap.parse_args()
    sys.exit(run(quick=a.quick, do_qber=not a.no_qber))
