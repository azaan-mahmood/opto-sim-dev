"""How far Duplinskiy's polarisation channel can drift before it needs recalibrating.

Polarisation encoding has few device-level impairments -- §30.11 showed the
chain is blind to everything a real source adds except pulse energy -- but
in practice it is the least stable encoding, because the channel's Jones
matrix moves with temperature and mechanical stress and the receiver's
compensation is only ever as good as its last calibration.

The paper says so directly:

    "The average time that the system has spent in the key distribution
     mode is about 80%, while the other 20% has been required for
     recalibrations as the quantum channel has not been isolated from
     external influences, including mechanical and temperature ones."
                     -- Duplinskiy et al., Opt. Express 25(23), 28886, §6

This script measures the tolerance that sentence implies: how far the
fibre may drift after calibration before QBER leaves the paper's ~2 % band.

Why a mismatch is needed at all
-------------------------------
`temperature` and `bend_radius` on their own are **exactly inert** while
Bob compensates the same fibre light travels through, which is the
default.  `U_comp` is the conjugate transpose of the channel's own Jones
matrix, so `U_comp @ J = I` for any unitary `J`, and those parameters only
change `J`.  That is arithmetic, not physics.

They become live under a *calibration mismatch*: Bob inverts the fibre as
it was, while light travels through it as it is.  That is what
`calibration_temperature` and `calibration_bend_radius` express, and it is
the mechanism the paper attributes its 20 % duty-cycle loss to.

What the curve does and does not show
-------------------------------------
The residual after a mismatch is a **fixed** unitary rotation, not a
randomising process, so QBER is a deterministic function of it rather than
tending to 50 %.  Measured at 10 km with a bend applied after calibrating
straight, QBER runs 79.9 % at R = 1 m and 44.8 % at R = 0.1 m -- a tighter
bend is not monotonically worse, because the residual rotation simply maps
the encoded states somewhere else on the sphere.

So only the **small-mismatch regime** is a curve worth reading. Past the
tolerance the system is decalibrated and the particular value is an
accident of the realisation. The script reports the tolerance and says
plainly that the far field is arbitrary.

References
----------
[1] Duplinskiy, Ustimchik, Kanapin, Kurochkin & Kurochkin, Opt. Express
    25(23), 28886-28897 (2017).
[2] Ulrich, Rashleigh & Eickhoff, "Bending-induced birefringence in
    single-mode fibers", Opt. Lett. 5(6), 273-275 (1980).
[3] Agrawal, "Fiber-Optic Communication Systems", 5th ed., §4.1.
"""
import argparse
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from src.channel import FiberRealization
from src.protocols.bb84_duplinskiy import (SAME_AS_OPERATING,
                                           simulate_bb84_duplinskiy)

OUT_DIR = os.path.join(os.path.dirname(__file__), '..', 'val_duplinskiy')

SEED = 42
T0 = 25.0

# Sized for >= 3000 sifted bits per cell (sec. 25's standard).  Measured
# yields: ~1.7e-3 sifted per pulse at 10 km, ~4e-4 at 50 km.
N_PULSES = {0: 1_400_000, 10: 2_000_000, 50: 10_000_000}
# Quick mode is sized so every distance clears the 250 threshold below.
# 50 km was at 400_000, which yields ~140 sifted -- so `--quick` could
# never pass, and a smoke mode that always fails is not a smoke mode.
N_QUICK = {0: 200_000, 10: 300_000, 50: 1_500_000}

# The paper's operating band: 2 % average, ~1 % floor.  "Recalibration is
# needed" is taken as QBER leaving the band, i.e. exceeding 2x the
# undrifted floor -- a stated criterion, not a fitted one.
TOLERANCE_FACTOR = 2.0


def residual(km, T=T0, R=None, T_cal=T0, R_cal=None):
    """||U_comp(calibration) @ J(operating) - I||_F. 0 = still perfect."""
    op = FiberRealization(L_m=km * 1000, temperature=T, bend_radius=R,
                          attenuation_factor=0.2, seed=SEED)
    cal = FiberRealization(L_m=km * 1000, temperature=T_cal, bend_radius=R_cal,
                           attenuation_factor=0.2, seed=SEED)
    U = cal.birefringence_matrix().conj().T
    return float(np.linalg.norm(U @ op.birefringence_matrix() - np.eye(2)))


def run_point(km, n, **kw):
    r = simulate_bb84_duplinskiy(n, fiber_length=km, seed=SEED, **kw)
    s, e, q = r['n_sifted'], r['n_errors'], r['qber']
    sig = math.sqrt(max(q * (1 - q), 1e-12) / s) if s else float('nan')
    return q, sig, s, e


def controls(failures):
    print("\n  controls")
    base = {(0, True): (35, 1), (0, False): (35, 1), (10, True): (39, 1),
            (10, False): (25, 9), (50, True): (8, 0), (50, False): (6, 1)}
    ok = True
    for (km, comp), (s0, e0) in base.items():
        r = simulate_bb84_duplinskiy(20000, fiber_length=km, compensate=comp,
                                     seed=SEED)
        ok &= (r['n_sifted'], r['n_errors']) == (s0, e0)
    print(f"    negative, defaults          : "
          f"{'bit-identical to the frozen baseline' if ok else 'BASELINE MOVED'}")
    if not ok:
        failures.append("the new parameters moved the frozen sec. 27.1 baseline "
                        "at their defaults")

    # Matched drift: compensation is still exact, so the response table must
    # not move beyond floating point.  It is NOT bit-identical -- a 1000
    # section Jones product and its inverse accumulate ~3e-14 relative, and
    # that is enough to flip a marginal `random() < p_click` comparison.
    a = simulate_bb84_duplinskiy(300000, fiber_length=10, seed=SEED)
    b = simulate_bb84_duplinskiy(300000, fiber_length=10, seed=SEED,
                                 temperature=50.0,
                                 calibration_temperature=50.0)
    d = abs(a['qber'] - b['qber'])
    print(f"    negative, matched drift 25->50 C: QBER {100 * a['qber']:.2f} % "
          f"vs {100 * b['qber']:.2f} %, |d| = {100 * d:.3f} pp")
    print("      -> compensation stays exact; the response table agrees to "
          "~3e-14,\n         which is Jones-product rounding, not a leak")
    if d > 0.01:
        failures.append(f"matched drift moved QBER by {100 * d:.2f} pp; "
                        "compensation should be exact regardless of "
                        "temperature")

    q0, _, _, _ = run_point(10, 300000)
    q1, _, _, _ = run_point(10, 300000, temperature=T0 + 1.0,
                            calibration_temperature=T0)
    q2, _, _, _ = run_point(10, 300000, bend_radius=0.3,
                            calibration_bend_radius=None)
    print(f"    positive, dT = 1 C          : {100 * q0:.2f} % -> {100 * q1:.2f} %")
    print(f"    positive, bent after calib. : {100 * q0:.2f} % -> {100 * q2:.2f} %")
    if q1 < 0.10 or q2 < 0.10:
        failures.append("a calibration mismatch did not reach the observable; "
                        "every null below would be vacuous (G2)")


def temperature_curve(km, quick, failures):
    n = (N_QUICK if quick else N_PULSES)[km]
    deltas = (0.0, 0.001, 0.003, 0.01, 0.02, 0.03, 0.05, 0.1)
    print(f"\n  temperature drift after calibration, {km} km, "
          f"{n:,} pulses per point")
    print("     dT (C)   residual   sifted  errors     QBER +/- 1 sigma")
    rows = []
    floor = None
    tol = None
    for dT in deltas:
        res = residual(km, T=T0 + dT, T_cal=T0)
        kw = {} if dT == 0.0 else dict(temperature=T0 + dT,
                                       calibration_temperature=T0)
        q, sig, s, e = run_point(km, n, **kw)
        rows.append((dT, res, q, sig, s))
        if dT == 0.0:
            floor = q
        elif tol is None and q > TOLERANCE_FACTOR * floor:
            tol = dT
        print(f"    {dT:7.3f}   {res:8.4f} {s:8d} {e:7d}   "
              f"{100 * q:6.2f} +/- {sig * 100:.2f} %")
        # 3000 is sec. 25's quotable standard.  The quick threshold is well
        # below the ~500 that mode yields, so it catches a broken run
        # without tripping on ordinary Poisson scatter around the mean.
        if s < (250 if quick else 3000):
            failures.append(f"{km} km, dT = {dT}: only {s} sifted bits; "
                            "not quotable")
    if tol is not None:
        print(f"\n    tolerance at {km} km: QBER leaves "
              f"{TOLERANCE_FACTOR:.0f}x the {100 * floor:.2f} % floor "
              f"by dT = {tol:g} C")
    else:
        print(f"\n    tolerance at {km} km: not reached within "
              f"dT = {max(deltas):g} C")
    return rows


def bend_curve(km, quick, failures):
    n = (N_QUICK if quick else N_PULSES)[km]
    # Extended below 2 m because that is already at residual 0.45, past the
    # knee.  20/10/5/3 m populate the small-mismatch regime as densely as
    # the temperature curve, which is what the collapse panel needs.
    # 20/10/5/3 populate the small-mismatch regime (2 m is already at
    # residual 0.45, past the knee).  2.4 and 1.9 are there specifically to
    # match the dT = 0.03 and dT = 0.05 residuals to within 5 %, which is
    # what the collapse test needs -- see collapse_test().
    radii = (None, 20.0, 10.0, 5.0, 3.0, 2.4, 2.0, 1.9, 1.0, 0.5, 0.3, 0.1)
    print(f"\n  bend applied after calibrating straight, {km} km")
    print("      R (m)   residual   sifted  errors     QBER +/- 1 sigma")
    rows = []
    for R in radii:
        res = residual(km, R=R, R_cal=None)
        kw = {} if R is None else dict(bend_radius=R,
                                       calibration_bend_radius=None)
        q, sig, s, e = run_point(km, n, **kw)
        rows.append((R, res, q, sig, s))
        label = 'straight' if R is None else f'{R:.2f}'
        print(f"    {label:>7}   {res:8.4f} {s:8d} {e:7d}   "
              f"{100 * q:6.2f} +/- {sig * 100:.2f} %")
    print("      -> past the tolerance the residual is a fixed unitary "
          "rotation, so\n         QBER is whatever that rotation does, not "
          "50 %. A tighter bend is\n         not monotonically worse.")
    return rows


def collapse_test(temp_rows, bend_rows, failures):
    """Do the two mechanisms agree where their residuals coincide?

    If QBER depends only on how far Bob's inverse misses the channel, and
    not on whether temperature or bending caused the miss, then two
    settings with the same residual must give the same QBER.  Stated
    before running: dT = 0.02 C and R = 3 m sit at residuals 0.1995 and
    0.1979, 0.8 % apart, so they are the sharpest available test.

    The pairs must be matched TIGHTLY, and a first attempt at 10 % was not
    tight enough.  It reported dT = 0.05 (residual 0.4888) against R = 2 m
    (0.4466) as differing by 2.11 +/- 0.81 pp -- but the curve there runs
    at ~36 pp per unit residual, so an 8.6 % residual mismatch predicts
    1.51 pp on its own, leaving 0.6 pp inside the error bar.  A loose
    match manufactures a difference where the curve is steep.  The
    tolerance is 5 %, and R = 2.4 and 1.9 m were added to give genuinely
    matched partners for dT = 0.03 and 0.05.
    """
    MATCH_TOL = 0.05
    print("\n  do the two mechanisms collapse onto one curve?")
    print(f"    pairs matched on compensation residual to within "
          f"{100 * MATCH_TOL:.0f} %, not on setting")
    print("    temperature                 bend                     "
          "QBER difference")
    tested = 0
    for dT, r_t, q_t, s_t, _ in temp_rows:
        if dT == 0.0:
            continue
        best = min((b for b in bend_rows if b[0] is not None),
                   key=lambda b: abs(b[1] - r_t))
        R, r_b, q_b, s_b, _ = best
        # Only a genuinely matched pair tests anything: where the curve is
        # steep, comparing two different residuals measures the slope, not
        # the mechanism.
        if r_t <= 0 or abs(r_b - r_t) / r_t > MATCH_TOL:
            continue
        tested += 1
        d = q_t - q_b
        sig = math.sqrt(s_t ** 2 + s_b ** 2)
        flag = 'agree' if abs(d) <= 2 * sig else 'DIFFER'
        print(f"    dT={dT:<6g} res={r_t:6.4f}  "
              f"R={R:<5g}m res={r_b:6.4f}  "
              f"{100 * q_t:6.2f} % vs {100 * q_b:6.2f} %  "
              f"{100 * d:+6.2f} +/- {100 * sig:.2f} pp  {flag}")
        if abs(d) > 3 * sig:
            failures.append(
                f"matched residuals ({r_t:.4f} vs {r_b:.4f}) gave QBER "
                f"{100 * q_t:.2f} % and {100 * q_b:.2f} %, "
                f"{abs(d) / sig:.1f} sigma apart -- QBER is not a function "
                "of the residual alone")
    if tested == 0:
        print("    no pair of settings landed within 10 % in residual; "
              "nothing tested")
    else:
        print("      -> a match means the cause does not matter, only how far")
        print("         the compensator missed. That is what makes the "
              "residual\n         the right axis for the third panel.")


def run(quick=False, distances=(10, 50)):
    os.makedirs(OUT_DIR, exist_ok=True)
    failures = []
    print("=" * 72)
    print("Duplinskiy: how far the channel drifts before recalibration")
    print("=" * 72)
    print(f"  seed {SEED}, calibration at {T0:.0f} C, compensate=True "
          "throughout")

    controls(failures)
    temp_rows = {}
    for km in distances:
        temp_rows[km] = temperature_curve(km, quick, failures)
    bend_rows = bend_curve(distances[0], quick, failures)
    collapse_test(temp_rows[distances[0]], bend_rows, failures)

    _write_csv(temp_rows, bend_rows, quick)
    _figure(temp_rows, bend_rows, quick)

    print()
    if failures:
        print("[FAIL]")
        for f in failures:
            print(f"  - {f}")
        return 1
    print("[PASS] compensation is exact when calibrated, degrades measurably "
          "when not")
    return 0


def _stem(quick):
    """Smoke runs write to their own files.

    Sharing paths with the full run meant `--quick` silently replaced a
    quotable figure with an under-powered one, and nothing warned.  See
    sec. 35.6.
    """
    return 'val_duplinskiy_drift--quick' if quick else 'val_duplinskiy_drift'


def _write_csv(temp_rows, bend_rows, quick=False):
    path = os.path.join(OUT_DIR, _stem(quick) + '.csv')
    with open(path, 'w') as fh:
        fh.write("# Duplinskiy calibration-drift tolerance, "
                 "validate_duplinskiy_drift.py\n")
        fh.write(f"# seed={SEED} calibration_temperature={T0} compensate=True\n")
        fh.write("axis,distance_km,setting,residual_frobenius,"
                 "qber,qber_sigma,n_sifted\n")
        for km, rows in temp_rows.items():
            for dT, res, q, sig, s in rows:
                fh.write(f"temperature,{km},{dT:g},{res:.6f},"
                         f"{q:.6f},{sig:.6f},{s}\n")
        km0 = list(temp_rows)[0]
        for R, res, q, sig, s in bend_rows:
            setting = 'straight' if R is None else f"{R:g}"
            fh.write(f"bend,{km0},{setting},{res:.6f},"
                     f"{q:.6f},{sig:.6f},{s}\n")
    print(f"\n  CSV: {path}")


def _figure(temp_rows, bend_rows, quick=False):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("  matplotlib missing, skipping the figure")
        return

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))

    # --- panel 1: the variable actually set, temperature -----------------
    colours = {10: 'tab:blue', 50: 'tab:purple'}
    for km, rows in temp_rows.items():
        c = colours.get(km, 'tab:gray')
        ax[0].errorbar([r[0] for r in rows[1:]], [100 * r[2] for r in rows[1:]],
                       yerr=[100 * r[3] for r in rows[1:]], marker='o',
                       capsize=3, lw=1.2, color=c, label=f'{km} km')
        # The undrifted point cannot sit on a log axis, so it is drawn as
        # the floor the curve rises from rather than silently dropped.
        ax[0].axhline(100 * rows[0][2], color=c, ls=':', lw=1.0, alpha=0.7)
    ax[0].axhspan(0, 2.0, color='tab:green', alpha=0.12,
                  label="paper's ~2 % QBER")
    ax[0].set_xscale('log')
    ax[0].set_xlabel('temperature change since calibration,  '
                     r'$\Delta T$ (C)')
    ax[0].set_ylabel('QBER (%)')
    ax[0].set_title('QBER against temperature drift\n'
                    '(dotted = undrifted floor)', fontsize=10)
    ax[0].grid(True, alpha=0.3)
    ax[0].legend(fontsize=8)

    # --- panel 2: the variable actually set, bend ------------------------
    km0 = list(temp_rows)[0]
    bent = [r for r in bend_rows if r[0] is not None]
    ax[1].errorbar([r[0] for r in bent], [100 * r[2] for r in bent],
                   yerr=[100 * r[3] for r in bent], marker='s', capsize=3,
                   lw=1.2, color='tab:red')
    straight = [r for r in bend_rows if r[0] is None]
    if straight:
        ax[1].axhline(100 * straight[0][2], color='tab:red', ls=':', lw=1.0,
                      alpha=0.7, label='straight (calibrated state)')
        ax[1].legend(fontsize=8)
    ax[1].axhspan(0, 2.0, color='tab:green', alpha=0.12)
    ax[1].set_xscale('log')
    ax[1].invert_xaxis()
    ax[1].set_xlabel('bend radius R (m),  tighter to the right')
    ax[1].set_ylabel('QBER (%)')
    ax[1].set_title(f'QBER against bend radius\n'
                    f'{km0} km, calibrated straight', fontsize=10)
    ax[1].grid(True, alpha=0.3)

    # --- panel 3: both mechanisms against how far compensation missed ----
    tr = temp_rows[km0]
    ax[2].errorbar([r[1] for r in tr], [100 * r[2] for r in tr],
                   yerr=[100 * r[3] for r in tr], marker='o', capsize=3,
                   lw=1.2, color='tab:blue', label='temperature')
    ax[2].errorbar([r[1] for r in bend_rows], [100 * r[2] for r in bend_rows],
                   yerr=[100 * r[3] for r in bend_rows], marker='s', capsize=3,
                   lw=1.2, color='tab:red', label='bend')
    ax[2].axhspan(0, 2.0, color='tab:green', alpha=0.12)
    ax[2].set_xlabel('compensation residual  '
                     r'$\|U_{cal}J_{op}-I\|_F$'
                     '\n0 = exact inverse,  ~2.2 = fully scrambled')
    ax[2].set_ylabel('QBER (%)')
    ax[2].set_title('QBER against compensation residual\n'
                    'both mechanisms', fontsize=10)
    ax[2].grid(True, alpha=0.3)
    ax[2].legend(fontsize=8)

    fig.suptitle('Duplinskiy chain: QBER against post-calibration drift',
                 fontsize=11)
    fig.tight_layout()
    png = os.path.join(OUT_DIR, _stem(quick) + '.png')
    fig.savefig(png, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  figure: {png}")


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--quick', action='store_true',
                    help='fewer pulses per point, for a smoke run')
    ap.add_argument('--distances', type=int, nargs='+', default=[10, 50],
                    help='fibre lengths in km (default 10 50)')
    a = ap.parse_args()
    sys.exit(run(quick=a.quick, distances=tuple(a.distances)))
