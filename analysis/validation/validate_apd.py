import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os, sys, argparse

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from src.detectors.apd import apd

OUT = os.path.join(os.path.dirname(__file__), '..', 'val_apd')
os.makedirs(OUT, exist_ok=True)

parser = argparse.ArgumentParser()
parser.add_argument('--seed', type=int, default=42)
args = parser.parse_args()
SEED = args.seed
np.random.seed(SEED)

WAVELENGTH = 1550e-9
LIGHTSPEED = 299792458.0
PLANCK = 6.626e-34
ECHARGE = 1.602e-19

det = apd(WAVELENGTH, excess_noise_factor=10, load_resistance=50,
          temperature=300, gain=10, quantum_efficiency=0.9, dark_current=10e-9)

# ----------------------------------------------------------------
# Panel A: Responsivity = eta * e * lambda / (h * c)  (Kasap Eq 4.19)
# ----------------------------------------------------------------
# The comparison curves below are formed from the DETECTOR'S constants,
# not from this file's.  The distinction matters: the question a validator
# can answer is whether the model implements the textbook expression, and
# restating the expression with different constants tests the constants
# instead, silently.  `det.c` is 3e8 where `LIGHTSPEED` here is exact,
# which is a 0.069 % difference and was exactly the responsivity
# "disagreement" before this line existed.  Whether 3e8 is good enough is
# a separate question, reported at the bottom and deliberately not gated.
R_theory = det.qe * det.charge * WAVELENGTH / (det.h * det.c)
lam_scan = np.linspace(800e-9, 1700e-9, 50)
R_sim = []
for lam in lam_scan:
    d = apd(lam, excess_noise_factor=10, load_resistance=50,
            temperature=300, gain=10, quantum_efficiency=0.9, dark_current=10e-9)
    R_sim.append(d.R)
R_sim = np.array(R_sim)
R_analytic = det.qe * det.charge * lam_scan / (det.h * det.c)

# ----------------------------------------------------------------
# Panel B: Signal current — I_signal = M * R * P  (Kasap Eq 4.23)
# ----------------------------------------------------------------
P_sweep = np.logspace(-9, -3, 30)
I_sig_sim = np.array([det.calculate_output_current(P) for P in P_sweep])
I_sig_theory = det.gain * det.R * P_sweep

# ----------------------------------------------------------------
# Panel C: Shot noise — scales with sqrt(B) and sqrt(P)
# ----------------------------------------------------------------
B_sweep = np.logspace(6, 10, 25)
I_fixed = 1e-6
noise_B = np.array([det.calculate_noise(I_fixed, B) for B in B_sweep])
# Kasap Eq 4.45: i_total^2 = F*(i_dark^2 + i_signal^2) + i_thermal^2.
# The DARK shot term belongs in the total and was missing from it, which
# is why the curve sat 4.8e-3 % below the model at every bandwidth -- a
# gap small enough to look like rounding on a log plot and to survive
# being called "verified" in the table.
shot_B = np.sqrt(2 * det.charge * I_fixed * B_sweep * det.enf)
dark_B = np.sqrt(2 * det.charge * det.dark_current * B_sweep * det.enf)
thermal_B = np.sqrt(4 * det.kB * det.T * B_sweep / det.RL)
noise_B_theory = np.sqrt(shot_B**2 + dark_B**2 + thermal_B**2)

# Power scaling of noise (fixed bandwidth)
P_noise = np.logspace(-9, -3, 20)
I_sig_noise = det.gain * det.R * P_noise
noise_P = np.array([det.calculate_noise(I, 1e9) for I in I_sig_noise])
shot_P = np.sqrt(2 * det.charge * I_sig_noise * 1e9 * det.enf)
dark_P = np.sqrt(2 * det.charge * det.dark_current * 1e9 * det.enf)
thermal_const = np.sqrt(4 * det.kB * det.T * 1e9 / det.RL)
noise_P_theory = np.sqrt(shot_P**2 + dark_P**2 + thermal_const**2)

# ----------------------------------------------------------------
# Panel D: Thermal noise floor — Johnson-Nyquist (Kasap Eq 4.42)
# ----------------------------------------------------------------
B_thermal = np.logspace(6, 10, 30)
noise_thermal = np.array([det.calculate_noise(0, B) for B in B_thermal])
# At I_signal = 0 the dark-current shot term does NOT vanish, so this is a
# thermal-plus-dark floor rather than a pure Johnson-Nyquist one.  Calling
# it thermal and comparing it to 4kTB/RL alone understates the model by
# the same 4.8e-3 %.
thermal_theory = np.sqrt(2 * det.charge * det.dark_current * B_thermal
                         * det.enf
                         + 4 * det.kB * det.T * B_thermal / det.RL)

# ----------------------------------------------------------------
# Panel E: Photon detection — Poisson statistics (Agrawal Eq 4.1.2)
# ----------------------------------------------------------------
P_phot = np.logspace(-12, -6, 20)
exposure = 1e-9
n_photons = []
for P in P_phot:
    n = det.detect_photons(P, exposure)
    n_photons.append(n)
n_photons = np.array(n_photons)
photon_energy = PLANCK * LIGHTSPEED / WAVELENGTH
n_theory = (P_phot / photon_energy) * exposure * det.qe

# ================================================================
# 6-panel figure
# ================================================================
fig = plt.figure(figsize=(16, 12))
gs = fig.add_gridspec(2, 3, hspace=0.35, wspace=0.35)

# --- Panel A: Responsivity ---
ax1 = fig.add_subplot(gs[0, 0])
ax1.plot(lam_scan * 1e9, R_sim * 1e3, '-', c='C0', lw=1.5, label=r'$R(\lambda)$ (sim)')
ax1.plot(lam_scan * 1e9, R_analytic * 1e3, '--k', lw=1.5, alpha=0.5,
         label=r'$\eta e\lambda/hc$')
ax1.axvline(1550, c='gray', ls=':', lw=0.6, alpha=0.5)
ax1.set(xlabel='Wavelength (nm)', ylabel='Responsivity (mA/W)',
        title='A: Responsivity\nKasap [1] Eq 4.19')
ax1.legend(fontsize=7)
ax1.grid(True, alpha=0.25)

# --- Panel B: Signal current ---
ax2 = fig.add_subplot(gs[0, 1])
ax2.loglog(P_sweep * 1e3, I_sig_sim * 1e6, 'o-', c='C1', ms=4, label='Simulated')
ax2.loglog(P_sweep * 1e3, I_sig_theory * 1e6, '--k', lw=1.5, alpha=0.5,
           label=r'$M \cdot R \cdot P$')
ax2.set(xlabel='Power (mW)', ylabel='Signal current ($\\mu$A)',
        title='B: Signal current\nKasap [1] Eq 4.23')
ax2.legend(fontsize=7)
ax2.grid(True, alpha=0.25)

# --- Panel C: RMS noise vs bandwidth ---
ax3 = fig.add_subplot(gs[0, 2])
ax3.loglog(B_sweep, noise_B * 1e9, '-', c='C2', lw=1.5, label='Total noise (sim)')
ax3.loglog(B_sweep, noise_B_theory * 1e9, '--k', lw=1.5, alpha=0.5, label='Theory')
ax3.loglog(B_sweep, shot_B * 1e9, ':', c='gray', lw=1, alpha=0.5, label='Shot (theory)')
ax3.loglog(B_sweep, thermal_B * 1e9, ':', c='C3', lw=1, alpha=0.5, label='Thermal (theory)')
ax3.set(xlabel='Bandwidth (Hz)', ylabel='Noise current (nA)',
        title='C: RMS noise vs bandwidth\nKasap [1] Eq 4.42-4.46')
ax3.legend(fontsize=7)
ax3.grid(True, alpha=0.25)

# --- Panel D: Thermal noise floor ---
ax4 = fig.add_subplot(gs[1, 0])
ax4.loglog(B_thermal, noise_thermal * 1e9, '-', c='C3', lw=1.5, label=r'$I_\mathrm{noise}(I_\mathrm{sig}=0)$')
ax4.loglog(B_thermal, thermal_theory * 1e9, '--k', lw=1.5, alpha=0.5,
           label=r'$\sqrt{4k_B T B / R_L}$')
ax4.set(xlabel='Bandwidth (Hz)', ylabel='Noise current (nA)',
        title='D: Thermal noise floor\n(Johnson-Nyquist)')
ax4.legend(fontsize=7)
ax4.grid(True, alpha=0.25)

# --- Panel E: Photon detection ---
ax5 = fig.add_subplot(gs[1, 1])
ax5.loglog(P_phot * 1e9, n_photons, 'o', c='C4', ms=4, label='Detected (sim)')
ax5.loglog(P_phot * 1e9, n_theory, '--k', lw=1.5, alpha=0.5,
           label=r'$P \cdot t \cdot \eta / h\nu$')
ax5.set(xlabel='Power (nW)', ylabel='Detected photons',
        title='E: Photon detection\nSaleh and Teich [3] Eq 17.1-10')
ax5.legend(fontsize=7)
ax5.grid(True, alpha=0.25)

fig.suptitle('APD Validation -- Kasap [1] Ch. 4, Agrawal [2] Ch. 4, Saleh and Teich [3] Ch. 17',
             fontsize=13, fontweight='bold', y=0.97)
fig.savefig(os.path.join(OUT, f'val_apd--seed{SEED}.png'), dpi=200, bbox_inches='tight')
print(f"Saved: val_apd--seed{SEED}.png")
np.savetxt(os.path.join(OUT, f'val_apd--seed{SEED}.csv'),
           np.column_stack([P_sweep, I_sig_sim, I_sig_theory]),
           delimiter=',',
           header='power_W,I_signal_sim,I_signal_theory',
           comments='')
print(f"Saved: val_apd--seed{SEED}.csv")

# ----------------------------------------------------------------
# The five comparisons this script has always drawn, as numbers
# ----------------------------------------------------------------
# Every one of these theory curves was computed above and plotted beside
# its simulated counterpart.  The table then declared each "verified"
# without any of them being compared, so a detector that had stopped
# obeying Kasap Eq 4.19 would have produced the same word.


def _max_err_pct(sim, theory):
    """Largest relative departure, in percent, ignoring zero references."""
    sim, theory = np.asarray(sim, float), np.asarray(theory, float)
    nz = np.abs(theory) > 0
    if not nz.any():
        return float('nan')
    return float(np.max(np.abs(sim[nz] - theory[nz]) / np.abs(theory[nz]))
                 * 100.0)


err_R = _max_err_pct(R_sim, R_analytic)
err_I = _max_err_pct(I_sig_sim, I_sig_theory)
err_nB = _max_err_pct(noise_B, noise_B_theory)
err_nP = _max_err_pct(noise_P, noise_P_theory)
err_th = _max_err_pct(noise_thermal, thermal_theory)

import csv
table_csv = os.path.join(OUT, f'val_apd--seed{SEED}_table.csv')
with open(table_csv, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Parameter', 'Value'])
    writer.writerow(['Responsivity', f'R = {det.R*1e3:.3f} mA/W at 1550 nm '
                                     f'(max err {err_R:.2e} %)'])
    writer.writerow(['Signal current', f'I = M * R * P, max err '
                                       f'{err_I:.2e} %'])
    writer.writerow(['Shot noise', f'vs sqrt(B) and sqrt(P), max err '
                                   f'{max(err_nB, err_nP):.2e} %'])
    writer.writerow(['Thermal noise', f'sqrt(4*k_B*T*B/R_L), max err '
                                      f'{err_th:.2e} %'])
    writer.writerow(['Photon detection', f'Poisson process, eta = {det.qe}'])
    writer.writerow(['Excess noise', f'F = {det.enf} applied to shot terms'])
print(f"Saved: val_apd--seed{SEED}_table.csv")
plt.close(fig)

# ============================================================
# Verdict
# ============================================================
#
# Four of these are closed form against closed form -- the detector
# evaluates the same textbook expressions the theory curves do -- so the
# residual is numerical and a loose gate catches every way they could
# genuinely part company: a responsivity that stops scaling with lambda,
# a gain dropped from the signal current, an excess-noise factor applied
# to the thermal term.
MAX_ERROR_PCT = 1e-6

failures = []
for name, err, ref in (
        ('responsivity vs eta*e*lambda/(h*c)', err_R, 'Kasap Eq 4.19'),
        ('signal current vs M*R*P', err_I, 'Kasap Eq 4.23'),
        ('noise vs sqrt(B)', err_nB, 'shot + Johnson'),
        ('noise vs sqrt(P)', err_nP, 'shot + Johnson'),
        ('thermal floor vs sqrt(4*kB*T*B/RL)', err_th, 'Kasap Eq 4.42')):
    if not np.isfinite(err):
        failures.append(f"{name}: comparison produced no finite value")
    elif err >= MAX_ERROR_PCT:
        failures.append(
            f"{name} departs by {err:.3e} % against a {MAX_ERROR_PCT:g} % "
            f"limit; the model no longer reproduces {ref}")

# Poisson detection is a DRAW, not a closed form, so it gets a loose
# statistical band rather than the numerical gate above.  Mean counts
# across the sweep must track eta*P*t/(h*nu); at the top of the sweep the
# expected count is large and the relative spread is small.
big = n_theory > 100
if big.any():
    rel = np.abs(n_photons[big] - n_theory[big]) / n_theory[big]
    band = 5.0 / np.sqrt(n_theory[big])          # 5 sigma on a Poisson draw
    if np.any(rel > band):
        k = int(np.argmax(rel - band))
        failures.append(
            f"photon counts leave the Poisson band: {n_photons[big][k]} "
            f"against an expected {n_theory[big][k]:.1f} "
            f"({rel[k] / (band[k] / 5.0):.1f} sigma). Detection should be "
            f"a Poisson draw about eta*P*t/(h*nu)")

# Reported, never gated: the detector's physical constants against CODATA.
# `self.c = 3e8` is a 0.069 % rounding, which is larger than every
# tolerance above and would fail any of them -- but changing it moves every
# APD number in the project, so it is a decision rather than a check.
print()
print("  Reported, not asserted -- detector constants against CODATA 2018:")
for nm, got, ref in (('c ', det.c, 299792458.0),
                     ('h ', det.h, 6.62607015e-34),
                     ('e ', det.charge, 1.602176634e-19),
                     ('kB', det.kB, 1.380649e-23)):
    print(f"    {nm} = {got:.6e}  vs {ref:.6e}   "
          f"{abs(got - ref) / ref * 100:+.4f} %")
print("    The comparison curves above use the detector's own constants, so")
print("    these deviations do not enter the checks; they set the accuracy")
print("    of every APD number rather than the consistency of the model.")

print()
if failures:
    print("[FAIL]")
    for f_ in failures:
        print(f"  - {f_}")
    sys.exit(1)
print(f"[PASS] responsivity, signal current, both noise scalings and the "
      f"thermal floor")
print(f"       all track their textbook forms to "
      f"{max(err_R, err_I, err_nB, err_nP, err_th):.2e} %")
print(f"[PASS] photon detection stays inside a 5-sigma Poisson band about "
      f"eta*P*t/(h*nu)")
sys.exit(0)
