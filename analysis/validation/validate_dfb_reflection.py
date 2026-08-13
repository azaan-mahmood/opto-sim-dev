"""Reflection spectrum of the DFB grating (Kim, Chung & Lee 2000, Fig. 4).

Reproduces the paper's reflection spectrum: power reflectivity and phase
against wavelength across the stopband.

This is deliberately built out of the model's OWN split-step section
operator rather than a textbook grating formula, so it tests the
implementation instead of testing an equation.  The derivation is written
out in full in dfblaser-notes.md; the short version follows.

One section maps the waves entering it to the waves leaving it, with
p = exp((G - j*delta)*dz), c00 = sech(gamma*dz), c01 = j*tanh(gamma*dz):

    F(i+1) = p * (c00 * F(i) + c01 * R(i+1))
    R(i)   = p * (c01 * F(i) + c00 * R(i+1))

That is a scattering form and does not cascade, because its inputs and
outputs sit at different planes.  Solving the second line for R(i+1) and
substituting into the first gives a transfer matrix on (F, R) at one
plane:

        [  p/c00      c01/c00   ]
    T = [                       ]
        [ -c01/c00    1/(p*c00) ]

which uses c00^2 - c01^2 = sech^2 + tanh^2 = 1.  det(T) = 1 exactly, and
p cancels out of it, so the identity holds under gain or loss.

The grating is uniform, so the device is T**N.  Illuminating from the left
with nothing incident from the right, R(N) = 0, the second row gives

    r = -T21 / T22

Validation target: at zero detuning the exact uniform-grating peak power
reflectivity is tanh^2(kappa*L).  The section operator must reproduce it.

References
----------
[1] Kim, Chung & Lee, "An Efficient Split-Step Time-Domain Dynamic
    Modeling of DFB/DBR Laser Diodes", IEEE J. Quantum Electron. 36(7),
    787-794 (2000).  Fig. 4 (reflection spectrum), Fig. 5 (convergence in
    kappa*dz), Sec. II ("15 subsections or more are enough" at 600 um).
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from src.lasers.dfblaser import Laser

OUT_DIR = os.path.join(os.path.dirname(__file__), '..', 'val_dfb')

# Exact uniform-grating zeros satisfy sqrt(delta^2 - kappa^2)*L = m*pi.
MAX_ZERO_ORDER = 8


def reflection(lam_m, n_sections, lossy=False):
    """Complex reflection coefficient r(lambda) for the grating.

    Returns (laser, delta, r).  `lossy=False` sets the field gain to zero,
    i.e. a passive grating, which is what makes the nulls true zeros; with
    the waveguide loss on they fill in.
    """
    las = Laser(n_sections=n_sections)
    c00, c01 = las._coupling()          # constant in wavelength
    dz, n = las.dz, las.n
    g = -las.alpha / 2.0 if lossy else 0.0

    delta = (2.0 * np.pi / lam_m) * las.n_eff_0 - np.pi / las.bragg_condition
    r = np.empty_like(lam_m, dtype=complex)
    for k, d in enumerate(delta):
        p = np.exp((g - 1j * d) * dz)
        t = np.array([[p / c00, c01 / c00],
                      [-c01 / c00, 1.0 / (p * c00)]], dtype=complex)
        tt = np.linalg.matrix_power(t, n)
        r[k] = -tt[1, 0] / tt[1, 1]
    return las, delta, r


def analytic_zeros(las, lam_lo, lam_hi, lam_ref=1550e-9):
    """Zero-reflection wavelengths, from the coupled-mode condition.

    Computed independently of the model so they are a real check on it.
    """
    kappa, length = las.kappa, las.grating_length
    dlam_dd = lam_ref ** 2 / (2.0 * np.pi * las.n_eff_0)
    out = []
    for m in range(1, MAX_ZERO_ORDER):
        d_m = np.sqrt(kappa ** 2 + (m * np.pi / length) ** 2)
        for sign in (-1.0, +1.0):
            lam0 = lam_ref - sign * d_m * dlam_dd
            if lam_lo <= lam0 <= lam_hi:
                out.append(lam0)
    return np.sort(np.array(out))


def run(n_sections=15, n_ref=200, lam_lo=1548e-9, lam_hi=1552e-9, points=4001):
    os.makedirs(OUT_DIR, exist_ok=True)
    lam = np.linspace(lam_lo, lam_hi, points)

    las, _, r = reflection(lam, n_sections)
    _, _, r_ref = reflection(lam, n_ref)
    _, _, r_lossy = reflection(lam, n_sections, lossy=True)

    R, R_ref = np.abs(r) ** 2, np.abs(r_ref) ** 2
    kappa, length = las.kappa, las.grating_length
    exact = float(np.tanh(kappa * length) ** 2)

    print("=" * 62)
    print("DFB grating reflection spectrum (Kim et al. 2000, Fig. 4)")
    print("=" * 62)
    print(f"  kappa = {kappa / 100:.0f} cm^-1   L = {length * 1e6:.0f} um   "
          f"kappa*L = {kappa * length:.2f}")
    print(f"  kappa*dz at N={n_sections}: {las.kappa_dz:.4f}")
    print(f"\n  peak power reflectivity")
    print(f"    N = {n_sections:<4d}          {R.max():.5f}")
    print(f"    N = {n_ref:<4d}          {R_ref.max():.5f}")
    print(f"    tanh^2(kappa*L)     {exact:.5f}   <- exact")
    err = abs(R.max() - exact)
    print(f"    |error| at N={n_sections}: {err:.2e}")

    zeros = analytic_zeros(las, lam_lo, lam_hi)
    print(f"\n  analytic zero-reflection wavelengths ({len(zeros)} in window):")
    for z in zeros:
        i = int(np.argmin(np.abs(lam - z)))
        print(f"    {z * 1e9:9.4f} nm   model |r|^2 = {R[i]:.3e}")
    print("  (not identically zero because the sweep is on a finite grid "
          "and\n   the nulls are narrower than the sample spacing)")

    phase = np.unwrap(np.angle(r))
    print(f"\n  unwrapped phase: {phase.min():.3f} to {phase.max():.3f} rad")

    csv_path = os.path.join(OUT_DIR, f'val_dfb_reflection--N{n_sections}.csv')
    with open(csv_path, 'w') as f:
        f.write(f"# DFB grating reflection spectrum, "
                f"validate_dfb_reflection.py\n")
        f.write(f"# kappa={kappa:.1f} m^-1 L={length:.6e} m N={n_sections} "
                f"kappa_dz={las.kappa_dz:.6f}\n")
        f.write(f"# peak |r|^2 {R.max():.6f} vs exact tanh^2(kappa L) "
                f"{exact:.6f}\n")
        f.write("wavelength_nm,detuning_m^-1,power_reflectivity,"
                "phase_rad,power_reflectivity_lossy\n")
        d = (2.0 * np.pi / lam) * las.n_eff_0 - np.pi / las.bragg_condition
        for k in range(len(lam)):
            f.write(f"{lam[k] * 1e9:.6f},{d[k]:.6e},{R[k]:.6e},"
                    f"{phase[k]:.6f},{np.abs(r_lossy[k]) ** 2:.6e}\n")
    print(f"\n  CSV: {csv_path}")

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("  matplotlib missing, skipping the figure")
        return R.max(), exact

    nm = lam * 1e9
    fig, ax1 = plt.subplots(figsize=(9, 5.5))
    ax1.plot(nm, R_ref, color='0.75', lw=3, label=f'N = {n_ref} (converged)')
    ax1.plot(nm, R, color='tab:blue', lw=1.4,
             label=f"N = {n_sections} (paper's choice)")
    ax1.plot(nm, np.abs(r_lossy) ** 2, color='tab:green', lw=1.0, ls=':',
             label=r'N = %d, with $\alpha$ = %.0f cm$^{-1}$'
                   % (n_sections, las.alpha / 100))
    ax1.set_xlabel('wavelength (nm)')
    ax1.set_ylabel(r'power reflectivity  $|r|^2$')
    ax1.set_ylim(0, 1)
    ax1.set_xlim(nm[0], nm[-1])
    ax1.grid(True, alpha=0.3)

    ax2 = ax1.twinx()
    ax2.plot(nm, phase, color='tab:red', lw=1.0, alpha=0.85, label='phase')
    ax2.set_ylabel('phase of $r$ (rad)', color='tab:red')
    ax2.set_ylim(-np.pi, 2 * np.pi)
    ax2.set_yticks([-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi,
                    3 * np.pi / 2, 2 * np.pi])
    ax2.set_yticklabels([r'$-\pi$', r'$-\pi/2$', '0', r'$\pi/2$', r'$\pi$',
                         r'$3\pi/2$', r'$2\pi$'])
    ax2.tick_params(axis='y', colors='tab:red')

    for z in zeros:
        ax1.annotate('', xy=(z * 1e9, 0.0), xytext=(z * 1e9, 0.13),
                     arrowprops=dict(arrowstyle='-|>', color='k', lw=1.1))

    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2, loc='upper right', fontsize=8, framealpha=0.9)
    ax1.set_title('DFB grating reflection spectrum  '
                  rf'($\kappa$ = {kappa / 100:.0f} cm$^{{-1}}$, '
                  rf'L = {length * 1e6:.0f} $\mu$m, '
                  rf'$\kappa L$ = {kappa * length:.1f})'
                  '\narrows mark the analytic zero-reflection wavelengths')
    fig.tight_layout()
    png = os.path.join(OUT_DIR, f'val_dfb_reflection--N{n_sections}.png')
    fig.savefig(png, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  figure: {png}")
    return R.max(), exact


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--sections', type=int, default=15,
                    help="sections in the device under test (default 15, "
                         "the paper's recommendation for 600 um)")
    ap.add_argument('--reference-sections', type=int, default=200,
                    help='sections for the converged reference curve')
    ap.add_argument('--points', type=int, default=4001,
                    help='wavelength samples across the window')
    ap.add_argument('--tolerance', type=float, default=1e-4,
                    help='allowed |peak - tanh^2(kappa L)| before failing')
    a = ap.parse_args()

    peak, exact = run(n_sections=a.sections, n_ref=a.reference_sections,
                      points=a.points)
    if abs(peak - exact) > a.tolerance:
        print(f"\n[FAIL] peak reflectivity off the exact value by "
              f"{abs(peak - exact):.2e} (tolerance {a.tolerance:g})")
        sys.exit(1)
    print(f"\n[PASS] peak reflectivity matches tanh^2(kappa*L) to "
          f"{abs(peak - exact):.2e}")
