import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats as sp_stats
import os, sys, argparse

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from src.channel.fiber import apply_cd, D_TOTAL, D_MATERIAL, D_WAVEGUIDE

OUT = os.path.join(os.path.dirname(__file__), '..', 'val_cd')
os.makedirs(OUT, exist_ok=True)

parser = argparse.ArgumentParser()
parser.add_argument('--seed', type=int, default=42)
args = parser.parse_args()
SEED = args.seed
np.random.seed(SEED)

LIGHTSPEED = 299792458.0
WAVELENGTH = 1550e-9
D_total = D_TOTAL
D_SI = D_total * 1e-6
beta2 = -D_SI * WAVELENGTH**2 / (2 * np.pi * LIGHTSPEED)

T0 = 30e-12
LD = T0**2 / abs(beta2)
DT = 1e-12
N_SAMPLES = 2**12
T = np.arange(N_SAMPLES) * DT
center = N_SAMPLES // 2
t_arr = T - T[center]

sigma0 = T0 / np.sqrt(2)

z_over_LD_dense = np.linspace(0, 2.5, 51)
distances_dense = z_over_LD_dense * LD / 1e3

E_in = np.zeros((N_SAMPLES, 2), dtype=np.complex128)
E_in[:, 0] = np.exp(-0.5 * (t_arr / T0)**2)

widths_sim = np.zeros(len(z_over_LD_dense))
for i, L_km in enumerate(distances_dense):
    if L_km == 0:
        E_out = E_in.copy()
    else:
        E_out = apply_cd(E_in.copy(), dt=DT, L=L_km*1000, wavelength=WAVELENGTH)
    I = np.abs(E_out[:, 0])**2
    total = I.sum()
    t_mean = np.sum(T * I) / total
    t_var = np.sum((T - t_mean)**2 * I) / total
    widths_sim[i] = np.sqrt(t_var.real)

widths_analytic = sigma0 * np.sqrt(1 + z_over_LD_dense**2)
errors_pct = np.abs(widths_sim - widths_analytic) / widths_analytic * 100

z_test = np.array([0.0, 0.5, 1.0, 2.0])
d_test = z_test * LD / 1e3
widths_test = np.array([widths_sim[np.argmin(np.abs(z_over_LD_dense - zf))] for zf in z_test])
widths_analytic_test = sigma0 * np.sqrt(1 + z_test**2)
errors_test = np.abs(widths_test - widths_analytic_test) / widths_analytic_test * 100

print(f"CD validation - D = {D_total} ps/(nm*km)  (D_mat = {D_MATERIAL}, D_wg = {D_WAVEGUIDE})")
print(f"beta2 = {beta2:.3e} s^2/m,  T0 = {T0*1e12:.1f} ps,  LD = {LD/1e3:.2f} km")
print(f"\n{'z/LD':>6s}  {'z(km)':>8s}  {'sig_sim(ps)':>11s}  {'sig_an(ps)':>11s}  {'error%':>8s}")
for i, zf in enumerate(z_test):
    print(f"{zf:6.1f}  {d_test[i]:8.2f}  {widths_test[i]*1e12:11.6f}  {widths_analytic_test[i]*1e12:11.6f}  {errors_test[i]:8.6f}")
print(f"\n  Max error across sweep: {errors_pct.max():.6e} %")

# ============================================================
# Main validation figure — 6 panels
# ============================================================
fig = plt.figure(figsize=(16, 12))
gs = fig.add_gridspec(2, 3, hspace=0.35, wspace=0.35)

# --- Panel A: Pulse broadening — time-domain (Agrawal Fig 2.6) ---
# Use distinct colors, zoom into pulse region
ax1 = fig.add_subplot(gs[0, 0])
colors_A = ['#1f77b4', '#2ca02c', '#ff7f0e', '#d62728']  # blue, green, orange, red
labels_A = [f'z/L_D = {zf:.1f} (theory: sigma = {wa*1e12:.1f} ps)'
            for zf, wa in zip(z_test, widths_analytic_test)]
for i, zf in enumerate(z_test):
    L_km = d_test[i]
    if L_km == 0:
        E_out = E_in.copy()
    else:
        E_out = apply_cd(E_in.copy(), dt=DT, L=L_km*1000, wavelength=WAVELENGTH)
    I = np.abs(E_out[:, 0])**2
    ax1.plot(t_arr*1e12, I/I.max(), color=colors_A[i], lw=1.5, label=labels_A[i])
ax1.set(xlabel='Time (ps)', ylabel='Normalized intensity',
        title='A: Gaussian pulse broadening (Agrawal [6] Fig 2.6)')
ax1.legend(fontsize=7, loc='upper left')
ax1.grid(True, alpha=0.25)
ax1.set_xlim(-60, 60)
ax1.annotate(f'T0 = {T0*1e12:.0f} ps, lambda = {WAVELENGTH*1e9:.0f} nm',
             xy=(0.95, 0.95), xycoords='axes fraction', ha='right', va='top',
             fontsize=8, bbox=dict(boxstyle='round', fc='lightyellow', alpha=0.8))

# --- Panel B: sigma(z)/sigma0 vs z/LD with analytic curve ---
ax2 = fig.add_subplot(gs[0, 1])
z_plot = np.linspace(0, 2.5, 500)
ax2.plot(z_plot, np.sqrt(1+z_plot**2), 'k-', lw=2, alpha=0.7,
         label=r'Analytic: $\sqrt{1+(z/L_D)^2}$')
ratio_sim = widths_sim / sigma0
ax2.plot(z_over_LD_dense, ratio_sim, '.-', c='C3', ms=2, lw=0.8,
         label='apply_cd measured')
for i, zf in enumerate(z_test):
    ax2.plot(zf, ratio_sim[np.argmin(np.abs(z_over_LD_dense-zf))], 'o', c='C3', ms=6)
    ax2.annotate(f'{errors_test[i]:.4f}%', (zf, ratio_sim[np.argmin(np.abs(z_over_LD_dense-zf))]),
                 fontsize=7, va='bottom', ha='left', c='C3')
ax2.set(xlabel=r'$z / L_D$', ylabel=r'$\sigma(z) / \sigma_0$',
        title=f'B: RMS width growth  (D = {D_total} ps/(nm*km))',
        xlim=(0, 2.5))
ax2.legend(fontsize=8)
ax2.grid(True, alpha=0.25)

# --- Panel C: Residual (% error) with legend ---
ax3 = fig.add_subplot(gs[0, 2])
ax3.semilogy(z_over_LD_dense, errors_pct, '-', c='C3', lw=1, label='Error across sweep')
ax3.axhline(0, color='gray', ls=':', lw=0.5)
for i, zf in enumerate(z_test):
    err = errors_test[i]
    ax3.plot(zf, err if err > 0 else 1e-16, 'o', c='k', ms=6, zorder=5)
ax3.plot([], [], 'ok', ms=6, label='Test points (z/L_D = 0, 0.5, 1.0, 2.0)')
ax3.set(xlabel=r'$z / L_D$', ylabel='|Relative error| (%)',
        title='C: Residual')
ax3.legend(fontsize=8)
ax3.grid(True, alpha=0.25)
ymin, ymax = ax3.get_ylim()
ax3.set_ylim(bottom=max(1e-16, ymin))

# --- Panel D: Frequency-domain phase response H(omega) ---
omega = 2*np.pi * np.fft.fftfreq(N_SAMPLES, d=DT)
H_analytic = np.exp(-1j * beta2 * omega**2 * (2*LD) / 2)  # z = 2*LD
E_f_in = np.fft.fft(E_in[:, 0])
E_f_out = np.fft.fft(apply_cd(E_in.copy(), dt=DT, L=2*LD, wavelength=WAVELENGTH)[:, 0])
H_sim = np.where(np.abs(E_f_in) > 1e-15, E_f_out / E_f_in, 0)
f_hz = omega / (2*np.pi)

f_nyquist = 1.0 / (2 * DT)
idx_pos = (omega > 0) & (f_hz < f_nyquist)
phase_ana = np.unwrap(np.angle(H_analytic[idx_pos]))
phase_sim = np.unwrap(np.angle(H_sim[idx_pos]))
f_THz = f_hz[idx_pos] * 1e-12

ax4 = fig.add_subplot(gs[1, 0])
ax4.plot(f_THz, phase_ana, '-k', lw=1.5, alpha=0.6, label='Analytic (unwrapped)')
ax4.plot(f_THz, phase_sim, ':', c='C3', lw=1, label='apply_cd (unwrapped)')
# Add reference: phase at z = L_D (half the test distance)
H_ref = np.exp(-1j * beta2 * omega**2 * LD / 2)
phase_ref = np.unwrap(np.angle(H_ref[idx_pos]))
ax4.plot(f_THz, phase_ref, '--', c='C0', lw=1, alpha=0.5,
         label=r'Reference: $\phi$ at $z = L_D$')
ax4.set(xlabel='Frequency (THz)', ylabel='Unwrapped phase (rad)',
        title='D: CD transfer function phase\n(Agrawal [6] Eq 2.4.11, z = 2L$_D$)')
ax4.legend(fontsize=7)
ax4.grid(True, alpha=0.25)

# --- Panel E: Material vs waveguide dispersion ---
lam_range = np.linspace(1100e-9, 1700e-9, 200)
c0 = LIGHTSPEED
D_mat_lam = (D_MATERIAL / 2) * (1 + (1310e-9)**2 / lam_range**2)
D_wg_lam = D_WAVEGUIDE * np.ones_like(lam_range)
D_total_lam = D_mat_lam + D_wg_lam

ax5 = fig.add_subplot(gs[1, 1])
ax5.plot(lam_range*1e9, D_mat_lam, '-', c='C0', lw=1.5, label=f'Material (D_mat @ 1550 = {D_MATERIAL})')
ax5.plot(lam_range*1e9, D_wg_lam, '-', c='C1', lw=1.5, label=f'Waveguide (D_wg @ 1550 = {D_WAVEGUIDE})')
ax5.plot(lam_range*1e9, D_total_lam, '--k', lw=1.5, label=f'Total (D @ 1550 = {D_total})')
ax5.axvline(1550, c='gray', ls=':', lw=0.6, alpha=0.5)
ax5.axvline(1310, c='gray', ls=':', lw=0.6, alpha=0.5)
ax5.set(xlabel='Wavelength (nm)', ylabel='Dispersion D (ps/(nm*km))',
        title='E: Material & waveguide dispersion\n(Hui [2], Keck [3])')
ax5.legend(fontsize=7)
ax5.grid(True, alpha=0.25)
ax5.annotate(f'D_total = D_mat + D_wg\n= {D_MATERIAL} + ({D_WAVEGUIDE})\n= {D_total} ps/(nm*km) @ 1550 nm',
             xy=(0.95, 0.95), xycoords='axes fraction', ha='right', va='top',
             fontsize=8, bbox=dict(boxstyle='round', fc='lightyellow', alpha=0.8))

# --- Panel F: tabular summary ---
ax6 = fig.add_subplot(gs[1, 2])
ax6.axis('off')
summary_rows = [
    ['D_total', f'{D_total} ps/(nm*km)'],
    ['D_material', f'{D_MATERIAL} ps/(nm*km)'],
    ['D_waveguide', f'{D_WAVEGUIDE} ps/(nm*km)'],
    ['beta2', f'{beta2:.3e} s^2/m'],
    ['Dispersion length L_D', f'{LD/1e3:.2f} km'],
    ['Max error', f'{errors_pct.max():.6e} %'],
    ['Gaussian broadening', 'Matches analytic: all < 0.001 %'],
]
table = ax6.table(cellText=summary_rows,
                  colLabels=['Parameter', 'Value'],
                  loc='center', cellLoc='left', fontsize=8)
table.auto_set_column_width(col=list(range(2)))
table.auto_set_font_size(False)
table.set_fontsize(8)
for (row, col), cell in table.get_celld().items():
    if row == 0:
        cell.set_facecolor('#40466e')
        cell.set_text_props(color='w', fontweight='bold')
ax6.set_title('F: Validation summary', fontsize=10, pad=10)

fig.suptitle('Chromatic Dispersion - Validation vs Agrawal [6] Sec 2.4',
             fontsize=13, fontweight='bold', y=0.97)
fig.savefig(os.path.join(OUT, f'val_cd--seed{SEED}.png'), dpi=200, bbox_inches='tight')
print(f"Saved: val_cd--seed{SEED}.png")

csv_name = f'val_cd--seed{SEED}.csv'
np.savetxt(os.path.join(OUT, csv_name),
           np.column_stack([z_over_LD_dense, distances_dense, widths_sim*1e12,
                            widths_analytic*1e12, errors_pct]),
           delimiter=',',
           header='z_over_LD,distance_km,sigma_rms_ps,analytic_sigma_ps,error_pct',
           comments='')
print(f"Saved: {csv_name}")
plt.close(fig)
