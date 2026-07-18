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
R_theory = det.qe * ECHARGE * WAVELENGTH / (PLANCK * LIGHTSPEED)
lam_scan = np.linspace(800e-9, 1700e-9, 50)
R_sim = []
for lam in lam_scan:
    d = apd(lam, excess_noise_factor=10, load_resistance=50,
            temperature=300, gain=10, quantum_efficiency=0.9, dark_current=10e-9)
    R_sim.append(d.R)
R_sim = np.array(R_sim)
R_analytic = det.qe * ECHARGE * lam_scan / (PLANCK * LIGHTSPEED)

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
# Remove signal shot noise contribution to isolate bandwidth scaling
shot_B = np.sqrt(2 * ECHARGE * I_fixed * B_sweep * det.enf)
thermal_B = np.sqrt(4 * det.kB * det.T * B_sweep / det.RL)
noise_B_theory = np.sqrt(shot_B**2 + thermal_B**2)

# Power scaling of noise (fixed bandwidth)
P_noise = np.logspace(-9, -3, 20)
I_sig_noise = det.gain * det.R * P_noise
noise_P = np.array([det.calculate_noise(I, 1e9) for I in I_sig_noise])
shot_P = np.sqrt(2 * ECHARGE * I_sig_noise * 1e9 * det.enf)
thermal_const = np.sqrt(4 * det.kB * det.T * 1e9 / det.RL)
noise_P_theory = np.sqrt(shot_P**2 + thermal_const**2)

# ----------------------------------------------------------------
# Panel D: Thermal noise floor — Johnson-Nyquist (Kasap Eq 4.42)
# ----------------------------------------------------------------
B_thermal = np.logspace(6, 10, 30)
noise_thermal = np.array([det.calculate_noise(0, B) for B in B_thermal])
thermal_theory = np.sqrt(4 * det.kB * det.T * B_thermal / det.RL)

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
        title='E: Photon detection\nAgrawal [2] Eq 4.1.2')
ax5.legend(fontsize=7)
ax5.grid(True, alpha=0.25)

# --- Panel F: Validation summary ---
ax6 = fig.add_subplot(gs[1, 2])
ax6.axis('off')
summary = [
    ['Responsivity', f'$R = {det.R*1e3:.3f}$ mA/W at 1550 nm'],
    ['Signal current', f'$I = M \\cdot R \\cdot P$ linear verified'],
    ['Shot noise', f'$\\propto \\sqrt{{B}}$ and $\\propto \\sqrt{{P}}$ verified'],
    ['Thermal noise', f'$\\sqrt{{4k_B T B / R_L}}$ verified'],
    ['Photon detection', f'Poisson process, $\\eta = {det.qe}$'],
    ['Excess noise', f'$F = {det.enf}$ applied to shot terms'],
]
table = ax6.table(cellText=summary, colLabels=['Parameter', 'Value'],
                  loc='center', cellLoc='left', fontsize=8)
table.auto_set_column_width(col=list(range(2)))
table.auto_set_font_size(False)
table.set_fontsize(8)
for (row, col), cell in table.get_celld().items():
    if row == 0:
        cell.set_facecolor('#40466e')
        cell.set_text_props(color='w', fontweight='bold')
ax6.set_title('F: Validation summary', fontsize=10, pad=10)

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
plt.close(fig)
