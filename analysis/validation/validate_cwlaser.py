import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats as sp_stats
import os, sys, argparse

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from src.lasers.cwlaser import CWLaser
from src.visualization import eye_diagram

OUT = os.path.join(os.path.dirname(__file__), '..', 'val_cwlaser')
os.makedirs(OUT, exist_ok=True)

parser = argparse.ArgumentParser()
parser.add_argument('--seed', type=int, default=42)
args = parser.parse_args()
SEED = args.seed
np.random.seed(SEED)

LIGHTSPEED = 299792458.0
PLANCK = 6.626e-34
WAVELENGTH = 1550e-9

# ----------------------------------------------------------------
# Panel A: Power calibration — mean(|E|^2) vs power_dbm
# ----------------------------------------------------------------
p_dbm_vals = np.linspace(-20, 10, 31)
P_meas = []
for p in p_dbm_vals:
    las = CWLaser(WAVELENGTH, power_dbm=p, linewidth=1e6, rin_density=-200)
    E = las.sample_field(1e-12, 2000)
    P_meas.append(np.mean(np.sum(np.abs(E)**2, axis=1)))
P_meas = np.array(P_meas)
P_theory = 10.0 ** (p_dbm_vals / 10.0) * 1e-3
p_err = np.abs(P_meas - P_theory) / P_theory * 100

# ----------------------------------------------------------------
# Panel B: Phase noise variance vs time (Wiener process, Henry [1])
# ----------------------------------------------------------------
las_ph = CWLaser(WAVELENGTH, power_dbm=0, linewidth=100e6, rin_density=-200)
dt_ph = 1e-12
N_ph = 50000
t_ph = np.arange(N_ph) * dt_ph * 1e9  # ns
phi = np.unwrap(np.angle(las_ph.sample_field(dt_ph, N_ph)[:, 0]))
phi_var = np.var(phi)
D_phi = 2.0 * np.pi * 100e6
theory_var = D_phi * t_ph[-1] / 1e9  # t in seconds

# Phase increment distribution
inc = np.diff(np.unwrap(np.angle(las_ph.sample_field(dt_ph, 20000)[:, 0])))
inc_theory_std = np.sqrt(D_phi * dt_ph)

# ----------------------------------------------------------------
# Panel C: Phase increment PDF (Gaussian, Henry [1])
# ----------------------------------------------------------------
inc_bins = 60
inc_hist, inc_edges = np.histogram(inc, bins=inc_bins, density=True)
inc_centers = 0.5 * (inc_edges[:-1] + inc_edges[1:])
inc_gauss = sp_stats.norm.pdf(inc_centers, 0, inc_theory_std)

# ----------------------------------------------------------------
# Panel D: RIN PSD — comparison to Coldren [2] Eq 5.3.38
# ----------------------------------------------------------------
las_rin = CWLaser(WAVELENGTH, power_dbm=0, linewidth=1e6, rin_density=-130,
                  relaxation_frequency=5e9, damping_rate=1.88e10)
dt_rin = 5e-12
N_rin = 2**15
E_rin = las_rin.sample_field(dt_rin, N_rin)
P_rin = np.sum(np.abs(E_rin)**2, axis=1)
P_mean = P_rin.mean()
P_norm = (P_rin - P_mean) / P_mean
fs_rin = 1.0 / dt_rin
freq_rin = np.fft.rfftfreq(N_rin, dt_rin) * 1e-9
PSD = np.abs(np.fft.rfft(P_norm))**2 / (N_rin * fs_rin)

omega_R = 2.0 * np.pi * 5e9
gamma = 1.88e10
omega_th = 2.0 * np.pi * freq_rin * 1e9
num = gamma**2 + omega_th**2
den = (omega_R**2 - omega_th**2)**2 + (gamma * omega_th)**2
H_sq = np.divide(num, den, out=np.zeros_like(num), where=den > 0)
H_sq /= H_sq[0]
RIN_lin = 10.0 ** (-130.0 / 10.0)
theory_psd = H_sq * RIN_lin

# ----------------------------------------------------------------
# Panel E: Polarization Jones vector — azimuth scan
# ----------------------------------------------------------------
psi_vals = np.linspace(0, np.pi, 50)
Ex_mag = []
Ey_mag = []
for psi in psi_vals:
    las_pol = CWLaser(WAVELENGTH, power_dbm=0, linewidth=1e6,
                      polarization_azimuth=psi, polarization_ellipticity=0)
    pol = las_pol._polarization_vector()
    Ex_mag.append(np.abs(pol[0]))
    Ey_mag.append(np.abs(pol[1]))
Ex_mag = np.array(Ex_mag)
Ey_mag = np.array(Ey_mag)
norm_check = np.sqrt(Ex_mag**2 + Ey_mag**2)

# The two claims the table has always made about polarisation, as
# numbers.  A linear state at azimuth psi is (cos psi, sin psi), and
# the Jones vector carries direction only -- power is set separately,
# so a non-unit norm would double-count one of them.
norm_dev = float(np.max(np.abs(norm_check - 1.0)))
azimuth_dev = float(max(
    np.max(np.abs(Ex_mag - np.abs(np.cos(psi_vals)))),
    np.max(np.abs(Ey_mag - np.abs(np.sin(psi_vals))))))

# ================================================================
# 6-panel figure
# ================================================================
fig = plt.figure(figsize=(16, 12))
gs = fig.add_gridspec(2, 3, hspace=0.35, wspace=0.35)

# --- Panel A: Power calibration ---
ax1 = fig.add_subplot(gs[0, 0])
ax1.plot(p_dbm_vals, P_theory * 1e3, '-k', lw=1.5, alpha=0.6,
         label=r'Theory: $P = 10^{P_\mathrm{dBm}/10}$ \si{mW}')
ax1.plot(p_dbm_vals, P_meas * 1e3, 'o', c='C0', ms=4, label='measured')
ax1.set(xlabel='Power (dBm)', ylabel='Power (mW)',
        title='A: Power calibration\n(Henry [1], Coldren [2])')
ax1.legend(fontsize=7)
ax1.grid(True, alpha=0.25)

# --- Panel B: Phase noise variance vs time ---
ax2 = fig.add_subplot(gs[0, 1])
ax2.plot(t_ph, np.cumsum(np.diff(phi, prepend=0)**2),
         '-', c='C1', lw=1, label=r'$\sum \Delta\phi^2$ (simulated)')
ax2.plot(t_ph, D_phi * t_ph / 1e9, '--k', lw=1.5, alpha=0.6,
         label=r'$D_\phi \cdot t$ (theory)')
ax2.set(xlabel='Time (ns)', ylabel=r'$\sum \Delta\phi^2$ (rad$^2$)',
        title='B: Phase diffusion (Wiener process)\nHenry [1] Eq~18')
ax2.legend(fontsize=7)
ax2.grid(True, alpha=0.25)

# --- Panel C: Phase increment distribution ---
ax3 = fig.add_subplot(gs[0, 2])
ax3.bar(inc_centers, inc_hist, width=inc_centers[1]-inc_centers[0],
        fc='C2', alpha=0.5, label='Increments')
ax3.plot(inc_centers, inc_gauss, '-k', lw=1.5, label='Gaussian fit')
ax3.set(xlabel=r'$\Delta\phi$ (rad)', ylabel='Probability density',
        title='C: Phase increment PDF\nHenry [1]')
ax3.legend(fontsize=7)
ax3.grid(True, alpha=0.25)

# --- Panel D: RIN PSD ---
ax4 = fig.add_subplot(gs[1, 0])
ax4.loglog(freq_rin, PSD, '-', c='C3', lw=0.8, alpha=0.7, label='Simulated PSD')
ax4.loglog(freq_rin, theory_psd, '--k', lw=1.5, label='Coldren Eq 5.3.38')
ax4.axvline(5, c='gray', ls=':', lw=0.6, alpha=0.5, label=r'$f_{RO} = 5$ GHz')
ax4.set(xlabel='Frequency (GHz)', ylabel='RIN PSD (Hz$^{-1}$)',
        title='D: RIN spectral density\nColdren [2] Eq 5.3.38')
ax4.legend(fontsize=7)
ax4.grid(True, alpha=0.25)

# --- Panel E: Polarization Jones vector ---
ax5 = fig.add_subplot(gs[1, 1])
ax5.plot(psi_vals, Ex_mag, '-', c='C0', lw=1.5, label=r'$|E_x|$')
ax5.plot(psi_vals, Ey_mag, '-', c='C1', lw=1.5, label=r'$|E_y|$')
ax5.plot(psi_vals, norm_check, '--k', lw=1, alpha=0.5, label=r'$\sqrt{|E_x|^2 + |E_y|^2}$')
ax5.set(xlabel='Azimuth $\\psi$ (rad)', ylabel='Magnitude',
        title='E: Polarization Jones vector\n(Yariv [3] Ch. 6)')
ax5.legend(fontsize=7)
ax5.grid(True, alpha=0.25)

# --- Panel F: NRZ-OOK eye through a real MZM ---
# The grid was always 2x3 and the last cell was empty.  A source is only
# half-characterised by its spectra: the eye is what a receiver actually
# sees, and it puts the phase noise and RIN of the panels above into the
# units a link budget is written in.
#
# `rin_density` is set well above the -200 dB/Hz used for the calibration
# panels, which is deliberately negligible: at that level the rails are
# flat and the eye shows nothing.  -130 dB/Hz is an ordinary DFB and makes
# the amplitude noise visible as rail thickness.
#
# Seeded twice over, and it needs both: `np.random.seed(SEED)` at the top
# of this file pins the laser's own noise, and `seed=SEED` pins the bit
# pattern, which is drawn from a separate generator.
ax6 = fig.add_subplot(gs[1, 2])
las_eye = CWLaser(WAVELENGTH, power_dbm=0, linewidth=1e6, rin_density=-130,
                  polarization_azimuth=np.pi / 4)
eye_diagram(las_eye, bitrate=10e9, n_bits=128, spb=64, ax=ax6, seed=SEED,
            title='F: NRZ-OOK eye, 10 Gbaud\nMZM at switching voltage')
ax6.grid(True, alpha=0.25)

fig.suptitle('CW Laser Validation — Henry [1], Coldren [2], Yariv [3]',
             fontsize=13, fontweight='bold', y=0.97)
fig.savefig(os.path.join(OUT, f'val_cwlaser--seed{SEED}.png'), dpi=200, bbox_inches='tight')
print(f"Saved: val_cwlaser--seed{SEED}.png")
np.savetxt(os.path.join(OUT, f'val_cwlaser--seed{SEED}.csv'),
           np.column_stack([p_dbm_vals, P_theory, P_meas, p_err]),
           delimiter=',',
           header='power_dbm,P_theory_W,P_meas_W,error_pct',
           comments='')
print(f"Saved: val_cwlaser--seed{SEED}.csv")

import csv
table_csv = os.path.join(OUT, f'val_cwlaser--seed{SEED}_table.csv')
with open(table_csv, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Parameter', 'Value'])
    writer.writerow(['Power calibration', f'max error = {p_err.max():.2e} %'])
    writer.writerow(['Phase diffusion coeff', f'D_phi = {D_phi:.2e} rad^2/s'])
    writer.writerow(['Phase increments', f'Gaussian, std = {inc_theory_std:.2e} rad'])
    writer.writerow(['RIN resonance', 'f_RO = 5 GHz, damping = 1.88e10 rad/s'])
    writer.writerow(['Polarization',
                     f'unit norm to {norm_dev:.2e}; azimuth follows '
                     f'(cos psi, sin psi) to {azimuth_dev:.2e}'])
print(f"Saved: val_cwlaser--seed{SEED}_table.csv")
plt.close(fig)

# ============================================================
# Verdict
# ============================================================
#
# Three of these are deterministic and get tight gates.  The two that are
# DRAWS -- phase diffusion and the increment width -- get statistical
# bands, because a Wiener process is only required to land near its
# variance, not on it.
#
# Phase variance is the loose one on purpose.  var(phi) over a single
# realization of a Wiener process is itself a random variable with a
# spread comparable to its mean, so the check is an order-of-magnitude
# band: it catches diffusion switched off, or running at the wrong D_phi,
# and does not pretend to more.
# Power is NOT a closed form here: the mean is taken over 2000 samples of
# a laser carrying RIN, so it scatters.  At -200 dB/Hz over a 5e11 Hz
# Nyquist band the per-sample relative RMS is sqrt(1e-20*5e11) = 7.1e-5,
# so the standard error of a 2000-sample mean is 1.6e-6, or 1.6e-4 %.
# 1e-2 % sits ~60x above that and far below any calibration error.
POWER_TOL_PCT = 1e-2
POL_TOL = 1e-9            # Jones vector algebra, deterministic
INC_STD_TOL_FRAC = 0.05   # sample std of ~20k Gaussian increments
# Quadratic variation of N increments has relative SD sqrt(2/N); at
# N = 50000 that is 0.63 %, so 5 % is ~8 sigma.
QV_TOL_FRAC = 0.05

failures = []

if p_err.max() >= POWER_TOL_PCT:
    worst = int(np.argmax(p_err))
    failures.append(
        f"optical power departs from 10^(dBm/10) mW by {p_err.max():.3e} % "
        f"at {p_dbm_vals[worst]:.1f} dBm ({P_meas[worst]:.6e} W against "
        f"{P_theory[worst]:.6e} W)")

if norm_dev >= POL_TOL:
    failures.append(
        f"the polarisation Jones vector is not unit norm (worst deviation "
        f"{norm_dev:.3e}); it carries direction, and power is set "
        f"separately, so a non-unit norm double-counts one of them")

# A linear state at azimuth psi is (cos psi, sin psi).  The table claimed
# "linear/circular verified" without this ever being compared.
if azimuth_dev >= POL_TOL:
    failures.append(
        f"a linear state at azimuth psi is (cos psi, sin psi); the scan "
        f"departs from it by {azimuth_dev:.3e}")

inc_std_dev = abs(inc.std() - inc_theory_std) / inc_theory_std
if inc_std_dev >= INC_STD_TOL_FRAC:
    failures.append(
        f"phase increments have std {inc.std():.4e} rad against "
        f"sqrt(D_phi*dt) = {inc_theory_std:.4e} ({inc_std_dev:.2%} out). "
        f"A Wiener process puts the linewidth entirely in this width")

# The quantity that equals D_phi*t is the QUADRATIC VARIATION -- the sum
# of squared increments -- which is exactly what panel B plots.  var(phi)
# across the trace is a different statistic entirely (a Brownian path
# scattered about its own time-average), and comparing it to D_phi*t is a
# category error rather than a loose check.
qv = float(np.sum(np.diff(phi, prepend=0.0) ** 2))
qv_dev = abs(qv - theory_var) / theory_var if theory_var else float('inf')
if qv_dev >= QV_TOL_FRAC:
    failures.append(
        f"phase diffusion is off: the quadratic variation is {qv:.4f} "
        f"rad^2 against D_phi*t = {theory_var:.4f} ({qv_dev:.2%} out). "
        f"A Wiener process accumulates variance linearly at D_phi = "
        f"2*pi*linewidth")

print()
if failures:
    print("[FAIL]")
    for f_ in failures:
        print(f"  - {f_}")
    sys.exit(1)
print(f"[PASS] optical power tracks 10^(dBm/10) mW to {p_err.max():.2e} % "
      f"over {p_dbm_vals.min():.0f} to {p_dbm_vals.max():.0f} dBm")
print(f"[PASS] the polarisation vector is unit norm to {norm_dev:.2e} and "
      f"follows (cos psi, sin psi) to {azimuth_dev:.2e}")
print(f"[PASS] phase increments are Gaussian with std sqrt(D_phi*dt) "
      f"({inc_std_dev:.2%} out), and the quadratic variation returns "
      f"D_phi*t to {qv_dev:.2%}")
sys.exit(0)
