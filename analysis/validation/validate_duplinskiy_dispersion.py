"""CD and PMD in a polarisation chain: one is inert by algebra, one is not.

found both impairments exposed but inert in this chain, and
blamed the single time sample: CD and PMD are frequency-domain operators,
and the FFT of one sample is the sample. That was right, and it left the
pair untestable until a source could supply a time-resolved field. The DFB
device model can, so this runs them.

Two separate findings, and they must not be confused with each other.

CD is inert for a reason that has nothing to do with sampling
------------------------------------------------------------
`apply_cd` multiplies BOTH field components by the same transfer function
H(f) -- chromatic dispersion is polarisation-independent. Everything
upstream in this chain keeps the two components proportional:

  - the source emits one complex amplitude times a fixed Jones vector,
    so Ey(t) = c*Ex(t) exactly (measured: the ratio is constant to 4e-17);
  - a phase modulator multiplies one component by a constant phase;
  - the fibre's birefringence is a constant 2x2 matrix.

A constant ratio in time is a constant ratio in frequency, so H(f) scales
both equally and Ey'(t) = c'*Ex'(t) still. The normalised Stokes vector
depends only on that ratio, so it cannot move. CD would therefore stay
inert here at ANY sampling rate, with any source that emits a single
polarisation state. This is the same cancellation as linewidth
argument and source-blindness result.

It is a structural null, not a measurement, and reporting it as evidence
that the CD model works would be wrong. The positive control below is the
point: CD perturbs the field enormously -- max|dE|/|E| above 1.0 -- while
moving the port powers by 2e-16.

PMD is live, and the source's chirp is why
------------------------------------------
`_apply_pmd_fixed` gives the two components OPPOSITE phase ramps, so it
breaks the proportionality that CD preserves. How much it breaks depends on
how much bandwidth the source has, because the ramp is `omega*dgd/2`.

This is the paper's own mechanism, stated in its sec. 4:

    "for pulses wider than a nanosecond the PMD itself does not produce
     significant error rate. However, power modulation of a semiconductor
     laser diode, which is used to generate the pulses, results in phase
     variation with time, i.e. in chirping of these pulses. This effect
     combined with the shift between the two orthogonal polarization
     components caused by the PMD leads to a significant degree of
     polarization degradation as SOP changes dramatically within a single
     pulse."

A gain-switched DFB chirps. So the chain should show exactly what the paper
describes: PMD alone is small, PMD plus a chirped source is not.

Note Bob's compensation cannot help. It inverts a fixed SU(2) matrix, and
PMD is frequency-dependent, so no static rotation undoes it.

Why the time base needed fixing first
-------------------------------------
The chain passed `dt=1/rep_rate` to `FiberRealization.apply` -- 100 ns, the
spacing BETWEEN pulses, where the operators need the spacing between
samples WITHIN one. Measured below: at 100 ns the PMD perturbation reads
2.8e-06 instead of 0.55, five orders down. That null looks exact and means
nothing, which is failure mode in a second costume.

References
----------
[1] Duplinskiy, Ustimchik, Kanapin, Kurochkin & Kurochkin, Opt. Express
    25(23), 28886-28897 (2017).
[2] Kim, Chung & Lee, IEEE J. Quantum Electron. 36(7), 787-794 (2000).
[3] Agrawal, G. P., "Nonlinear Fiber Optics", 5th ed., Academic Press,
    2013, sec. 2.4 (chromatic dispersion) and sec. 1.2.3 (PMD).
"""
import argparse
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from src.channel import FiberRealization, optics
from src.lasers import DFBLaser, DriveParams, LaserDriver
from src.protocols.bb84_duplinskiy import simulate_bb84_duplinskiy

OUT_DIR = os.path.join(os.path.dirname(__file__), '..', 'val_duplinskiy')

SEED = 42
LASER_SEED = 11
N_SECTIONS = 15

# The gain-switched operating point characterised in validate_dfb_drive.py.
# 150 ps sits mid-window between 100 and 200 ps ends (A9).
I_BIAS, I_PEAK = 0.060, 0.140
HARVEST_PERIOD = 5e-9
PULSE_WIDTH = 150e-12
N_SAMPLES = 8192

DISTANCES = (10, 50, 100)
# Measured sifted yield per pulse at mu = 0.1, from sweep.
YIELD = {10: 1.72e-3, 50: 3.38e-4, 100: 3.25e-5}
TARGET_SIFTED = 3000

CASES = (('no dispersion', {}), ('CD only', dict(cd=True)),
         ('PMD only', dict(pmd=True)))


def _stem(quick):
    return ('val_duplinskiy_dispersion--quick' if quick
            else 'val_duplinskiy_dispersion')


def source(dt=None, n=N_SAMPLES):
    """The DFB field at PM1's input, on the device's own time grid.

    The grid matters here in a way it does not for the polarisation-only
    validators, which evaluate a single time sample and so have no
    bandwidth to lose.  `sample_field` decimates by averaging, which
    band-limits, and PMD's effect scales with the source bandwidth -- so
    a coarse grid does not merely add noise, it systematically
    understates the impairment.  Quantified in `grid_scan`.
    """
    las = DFBLaser(n_sections=N_SECTIONS, seed=LASER_SEED)
    drive = DriveParams(mode='gain_switched', waveform='gaussian',
                        i_bias=I_BIAS, i_peak=I_PEAK,
                        period=HARVEST_PERIOD, width=PULSE_WIDTH)
    drv = LaserDriver(las, drive, seed=LASER_SEED,
                      polarization_azimuth=np.pi / 4)
    dt = las.dt if dt is None else dt
    return drv.sample_field(dt, n), dt, las.dt


def _dop(F):
    """Degree of polarisation of the time-averaged coherency matrix.

    DOP = 1 means every time sample carries the same polarisation state.
    Below 1 means the state varies across the pulse, which is precisely
    what the paper says PMD plus chirp does.
    """
    Ex, Ey = F[:, 0], F[:, 1]
    s0 = np.mean(np.abs(Ex) ** 2 + np.abs(Ey) ** 2)
    s1 = np.mean(np.abs(Ex) ** 2 - np.abs(Ey) ** 2)
    s2 = 2 * np.mean(np.real(Ex * np.conj(Ey)))
    s3 = -2 * np.mean(np.imag(Ex * np.conj(Ey)))
    return float(np.sqrt(s1 ** 2 + s2 ** 2 + s3 ** 2) / s0)


def _ports(F):
    Ex, Ey = optics.circular_analyser(F)
    return float(np.mean(np.abs(Ex) ** 2)), float(np.mean(np.abs(Ey) ** 2))


def _fibre(km, **kw):
    cfg = dict(birefringence=False, cd=False, pmd=False, attenuation=False)
    cfg.update(kw)
    return FiberRealization(L_m=km * 1000, seed=SEED, **cfg)


def field_level(E, dt, failures):
    print("\n  the field level: what each operator does before detection")

    ratio = E[:, 1] / E[:, 0]
    spread = float(np.std(ratio))
    print(f"    source Ey/Ex is constant to {spread:.2e}  "
          "-- one amplitude, one Jones vector")
    if spread > 1e-12:
        failures.append(f"the source's Ey/Ex varies by {spread:.2e}; the CD "
                        "argument below assumes a common envelope")

    print(f"\n    {'':<10}{'max|dE|/|E|':>14}{'Ey/Ex drift':>14}"
          f"{'port powers':>14}{'DOP':>10}")
    rows = {}
    for km in DISTANCES:
        for lab, kw in (('CD', dict(cd=True)), ('PMD', dict(pmd=True))):
            fib = _fibre(km, **kw)
            out = fib.apply(E, dt=dt)
            d = float(np.max(np.abs(out - E)) / np.max(np.abs(E)))
            # Only where Ex carries real amplitude.  In the dark gaps
            # between pulses both components are near zero and their ratio
            # is numerically meaningless, which would put a large number in
            # a column whose whole purpose is to show a small one.
            m = np.abs(out[:, 0]) > 1e-3 * np.max(np.abs(out[:, 0]))
            drift = float(np.max(np.abs(out[m, 1] / out[m, 0] - ratio[m])))
            p0, p1 = _ports(E), _ports(out)
            move = max(abs(p1[0] - p0[0]), abs(p1[1] - p0[1])) / max(p0)
            dop = _dop(out)
            rows[(km, lab)] = (d, drift, move, dop, fib._dgd)
            print(f"    {km:3d} km {lab:<4}{d:14.4f}{drift:14.2e}"
                  f"{move:14.2e}{dop:10.6f}")

    print("\n    CD: the field is perturbed hard and the observable does not")
    print("        move at all. That is algebra, not a working model --")
    print("        H(f) multiplies both components equally.")
    print("    PMD: opposite phase ramps, so the state varies across the")
    print("         pulse and DOP falls. The paper's sec. 4 mechanism.")

    for km in DISTANCES:
        d, drift, move, dop, _ = rows[(km, 'CD')]
        if d < 0.1:
            failures.append(f"CD at {km} km barely perturbs the field "
                            f"({d:.2e}); the null below would be vacuous")
        if move > 1e-12:
            failures.append(f"CD at {km} km moved the port powers by "
                            f"{move:.2e}; it should be exactly inert")
        d, _, move, dop, dgd = rows[(km, 'PMD')]
        if dop > 0.999:
            failures.append(f"PMD at {km} km left DOP at {dop:.6f}; it is not "
                            "reaching the observable")
    return rows


def false_null(E, dt, failures):
    print("\n  the wrong time base, and the null it produces")
    fib = _fibre(100, pmd=True)
    got = {}
    for lab, d_t in (('sample interval', dt), ('pulse period, 10 MHz', 1e-7)):
        out = fib.apply(E, dt=d_t)
        got[lab] = float(np.max(np.abs(out - E)) / np.max(np.abs(E)))
        print(f"    dt = {d_t:.3e} s  ({lab:<21}) "
              f"max|dE|/|E| = {got[lab]:.3e}")
    print("    -> the chain used to pass the pulse period. Five orders of")
    print("       magnitude of impairment vanish, and nothing errors.")
    print("       `source_dt` is now required when cd or pmd is on and the")
    print("       field has more than one sample.")
    if got['pulse period, 10 MHz'] > 1e-3:
        failures.append("the pulse period no longer suppresses PMD; the "
                        "argument for requiring source_dt needs rechecking")


def grid_scan(failures):
    """How much PMD you see depends on how well the source is resolved."""
    print("\n  resolving the source: PMD scales with the bandwidth kept")
    _, _, dt_dev = source()
    print(f"    device step {dt_dev:.4e} s, Nyquist {0.5 / dt_dev / 1e9:.0f} GHz")
    print(f"    {'sample dt':>12}{'Nyquist':>12}{'DOP at 100 km':>16}")
    rows = []
    for mult in (1, 2, 4, 8):
        dt = dt_dev * mult
        E, _, _ = source(dt=dt, n=N_SAMPLES // mult)
        out = _fibre(100, pmd=True).apply(E, dt=dt)
        rows.append((dt, _dop(out)))
        print(f"    {dt:12.3e}{0.5 / dt / 1e9:10.0f} GHz{_dop(out):16.6f}")
    print("    -> a coarser grid averages the chirp away and understates the")
    print("       impairment. This is not numerical noise, it is a systematic")
    print("       loss of the bandwidth PMD acts on, so the device grid is")
    print("       the honest choice and every number above uses it.")
    if rows[-1][1] < rows[0][1]:
        failures.append("coarsening the grid increased the measured "
                        "depolarisation; the grid argument is backwards")
    return rows


def qber_sweep(E, dt, quick, failures):
    print("\n  QBER, DFB source, compensate=True, seed 42")
    print(f"    {'':<16}" + "".join(f"{lab:>22}" for lab, _ in CASES))
    cells = {}
    for km in DISTANCES:
        n = int(TARGET_SIFTED / YIELD[km] * 1.15)
        if quick:
            n //= 12
        row = f"    {km:3d} km {n:9,}"
        for lab, kw in CASES:
            r = simulate_bb84_duplinskiy(n, fiber_length=km, seed=SEED,
                                         source_field=E, source_dt=dt, **kw)
            q, s = r['qber'], r['n_sifted']
            sig = math.sqrt(max(q * (1 - q), 1e-12) / s) if s else float('nan')
            cells[(km, lab)] = (q, sig, s)
            row += f"{100 * q:12.2f} +/-{100 * sig:5.2f}"
            if s < (250 if quick else TARGET_SIFTED):
                failures.append(f"{km} km, {lab}: only {s} sifted; not quotable")
        print(row)
    print("    (QBER %, +/- 1 sigma)")

    print("\n    CD against no dispersion: bit-identical at every distance?")
    for km in DISTANCES:
        a, b = cells[(km, 'no dispersion')], cells[(km, 'CD only')]
        same = (a[0], a[2]) == (b[0], b[2])
        print(f"      {km:3d} km  {'yes' if same else 'NO'}  "
              f"({a[2]} vs {b[2]} sifted, {100 * a[0]:.4f} vs {100 * b[0]:.4f} %)")
        if not same:
            failures.append(f"CD moved QBER at {km} km; the algebra says it "
                            "cannot")

    print("\n    PMD against no dispersion:")
    for km in DISTANCES:
        a, b = cells[(km, 'no dispersion')], cells[(km, 'PMD only')]
        d = b[0] - a[0]
        ds = math.hypot(a[1], b[1])
        print(f"      {km:3d} km  {100 * a[0]:6.2f} % -> {100 * b[0]:6.2f} %   "
              f"{100 * d:+6.2f} +/-{100 * ds:5.2f} pp   {d / ds:5.1f} sigma")
        # Quick mode carries ~290 sifted per cell against ~3400, so the
        # same swing lands near 2 sigma there.  A threshold a smoke run
        # cannot meet turns the smoke run into a permanent failure, which
        # is why the floor scales with the budget rather than sitting at
        # the quotable value.
        floor = 1.5 if quick else 3.0
        if d < floor * ds:
            failures.append(f"PMD at {km} km moved QBER by only "
                            f"{d / ds:.1f} sigma; the swing is not established")
    return cells


def _write_csv(rows, cells, grid, quick):
    path = os.path.join(OUT_DIR, _stem(quick) + '.csv')
    with open(path, 'w') as fh:
        fh.write("# CD and PMD in the Duplinskiy chain, "
                 "validate_duplinskiy_dispersion.py\n")
        fh.write(f"# DFB source, {PULSE_WIDTH * 1e12:g} ps gain-switched, "
                 f"device grid, seed={SEED} laser_seed={LASER_SEED}\n")
        fh.write("section,distance_km,case,value,extra\n")
        for (km, lab), (d, drift, move, dop, dgd) in rows.items():
            fh.write(f"field,{km},{lab}_perturbation,{d:.6e},\n")
            fh.write(f"field,{km},{lab}_port_move,{move:.6e},\n")
            fh.write(f"field,{km},{lab}_dop,{dop:.6f},{dgd:.6e}\n")
        for (km, lab), (q, sig, s) in cells.items():
            fh.write(f"qber,{km},{lab},{q:.6f},{sig:.6f}\n")
        for dt, dop in grid:
            fh.write(f"grid,100,dop,{dop:.6f},{dt:.6e}\n")
    print(f"\n  CSV: {path}")


def _figure(rows, cells, grid, quick):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("  matplotlib missing, skipping the figure")
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.4, 4.8),
                                   gridspec_kw={'width_ratios': [3, 2]})

    x = np.arange(len(DISTANCES))
    w = 0.26
    colours = ('0.72', 'tab:blue', 'tab:red')
    for j, (lab, _) in enumerate(CASES):
        q = [100 * cells[(km, lab)][0] for km in DISTANCES]
        e = [100 * cells[(km, lab)][1] for km in DISTANCES]
        ax1.bar(x + (j - 1) * w, q, w, yerr=e, capsize=4, color=colours[j],
                edgecolor='0.25', linewidth=0.7, label=lab)
        for i in x:
            ax1.text(i + (j - 1) * w, q[i] + e[i] + 0.3, f'{q[i]:.2f}',
                     ha='center', fontsize=8)
    ax1.set_xticks(x)
    ax1.set_xticklabels([f'{km} km' for km in DISTANCES])
    ax1.set_xlabel('fibre length')
    ax1.set_ylabel('QBER (%)')
    ax1.set_title('QBER with chromatic dispersion and with PMD\n'
                  'DFB source, polarisation compensation on', fontsize=11)
    ax1.grid(True, axis='y', alpha=0.3)
    ax1.legend(fontsize=9, loc='upper left')

    dops_pmd = [rows[(km, 'PMD')][3] for km in DISTANCES]
    dops_cd = [rows[(km, 'CD')][3] for km in DISTANCES]
    ax2.plot(DISTANCES, dops_pmd, 'o-', color='tab:red', label='PMD')
    ax2.plot(DISTANCES, dops_cd, 's--', color='tab:blue', label='CD')
    ax2.axhline(1.0, color='0.4', lw=1.0)
    ax2.set_xlabel('fibre length (km)')
    ax2.set_ylabel('degree of polarisation')
    ax2.set_ylim(0.4, 1.05)
    ax2.set_title('polarisation held across the pulse', fontsize=11)
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
    print("=" * 74)
    print("CD and PMD in a polarisation chain, with a time-resolved source")
    print("=" * 74)

    E, dt, dt_dev = source()
    print(f"  DFB gain-switched, {PULSE_WIDTH * 1e12:g} ps, "
          f"{N_SAMPLES} samples at dt = {dt:.4e} s")
    print(f"  span {N_SAMPLES * dt * 1e9:.2f} ns, "
          f"Nyquist {0.5 / dt / 1e9:.0f} GHz")
    print("\n  stated before the run: CD must be EXACTLY inert here, because")
    print("  it multiplies both components by the same H(f) and the source")
    print("  keeps them proportional. PMD must not be, because it gives them")
    print("  opposite phase ramps. Verdicts are read off the numbers.")

    rows = field_level(E, dt, failures)
    false_null(E, dt, failures)
    grid = grid_scan(failures)
    cells = qber_sweep(E, dt, quick, failures)

    _write_csv(rows, cells, grid, quick)
    _figure(rows, cells, grid, quick)

    print()
    if failures:
        print("[FAIL]")
        for f in failures:
            print(f"  - {f}")
        return 1
    print("[PASS] CD is inert by algebra with a positive control on the field;")
    print("       PMD swings the observable, completing a second pair against")
    print("       null in the time-bin chain")
    return 0


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--quick', action='store_true',
                    help='a twelfth of the pulses, for a smoke run')
    a = ap.parse_args()
    sys.exit(run(quick=a.quick))
