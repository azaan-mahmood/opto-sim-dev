import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os, sys, argparse
from scipy import stats as sp_stats

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from src.channel.fiber import apply_attenuation

OUT = os.path.join(os.path.dirname(__file__), '..', 'val_attenuation')
os.makedirs(OUT, exist_ok=True)

parser = argparse.ArgumentParser()
parser.add_argument('--seed', type=int, default=42)
args = parser.parse_args()
SEED = args.seed
np.random.seed(SEED)

N_SAMPLES = 50000
E_in = np.zeros((N_SAMPLES, 2), dtype=complex)
E_in[:, 0] = np.sqrt(0.5) * (np.random.randn(N_SAMPLES) + 1j * np.random.randn(N_SAMPLES))
P_in = np.mean(np.abs(E_in)**2)

ALPHA_NOMINAL = 0.182
distances_km = np.linspace(0, 200, 81)
P_out_theory = P_in * 10 ** (-ALPHA_NOMINAL * distances_km / 10.0)
P_out_meas = np.array([
    np.mean(np.abs(apply_attenuation(E_in.copy(), L_km, ALPHA_NOMINAL))**2)
    for L_km in distances_km
])
errors_pct = np.abs(P_out_meas - P_out_theory) / P_out_theory * 100
loss_dB = -10 * np.log10(P_out_meas / P_in)

slope, intercept, r_value, p_value, std_err = sp_stats.linregress(
    distances_km, loss_dB)
alpha_fitted = slope
r2 = r_value**2
print(f"Attenuation validation - alpha = {ALPHA_NOMINAL} dB/km")
print(f"  Fitted alpha = {alpha_fitted:.6f} dB/km  (R^2 = {r2:.10f})")
print(f"  Max |error| = {errors_pct.max():.6e} %")

fig = plt.figure(figsize=(14, 9))
gs = fig.add_gridspec(2, 3, hspace=0.35, wspace=0.40)

# --- Panel A: Exponential power decay (Keiser Fig 3.2 style) ---
ax1 = fig.add_subplot(gs[0, 0])
ax1.semilogy(distances_km, P_out_theory, 'k-', lw=2, label='Keiser Eq 3.6')
ax1.semilogy(distances_km, P_out_meas, 'o', ms=3, c='C3', label='apply_attenuation')
ax1.set(xlabel='Distance (km)', ylabel='Optical power (a.u.)',
        title='A: Power decay (Keiser [1] Eq 3.6)')
ax1.legend(fontsize=8, loc='upper right')
ax1.grid(True, alpha=0.25)
ax1.annotate(f'alpha = {ALPHA_NOMINAL} dB/km\nSMF-28 @ 1550 nm',
             xy=(0.95, 0.95), xycoords='axes fraction', ha='right', va='top',
             fontsize=8, bbox=dict(boxstyle='round', fc='lightyellow', alpha=0.8))

# --- Panel B: dB loss vs distance with linear fit ---
ax2 = fig.add_subplot(gs[0, 1])
ax2.plot(distances_km, loss_dB, 'o', ms=3, c='C0', label='Measured loss')
L_dense = np.linspace(0, 200, 400)
ax2.plot(L_dense, ALPHA_NOMINAL * L_dense, '--k', lw=1, label=f'Nominal α = {ALPHA_NOMINAL}')
ax2.plot(L_dense, alpha_fitted * L_dense, '-', c='C3', lw=1.5,
         label=f'Fit α = {alpha_fitted:.5f} dB/km')
ax2.set(xlabel='Distance (km)', ylabel='Loss (dB)',
        title=f'B: Loss vs distance  (R² = {r2:.10f})')
ax2.legend(fontsize=8)
ax2.grid(True, alpha=0.25)

# --- Panel C: Residual ---
ax3 = fig.add_subplot(gs[0, 2])
ax3.plot(distances_km, errors_pct, '-', c='C3', lw=1)
ax3.axhline(0, color='gray', ls=':', lw=0.5)
ax3.set(xlabel='Distance (km)', ylabel='|Relative error| (%)',
        title='C: Residual')
ax3.set_yscale('log')
ax3.grid(True, alpha=0.25)
ymin, ymax = ax3.get_ylim()
if ymin < 1e-16:
    ax3.set_ylim(bottom=1e-16)

# --- Panel D: Attenuation coefficient at multiple λ (literature ref) ---
lam_nm = np.array([800, 1000, 1200, 1310, 1385, 1450, 1550, 1625])
typical_dB_per_km = np.array([2.5, 1.8, 0.5, 0.34, 0.55, 0.28, 0.182, 0.20])
ax4 = fig.add_subplot(gs[1, 0])
ax4.semilogy(lam_nm, typical_dB_per_km, 'o-', c='gray', lw=1.5, zorder=1)
ax4.axvline(1550, c='C3', ls='--', lw=0.8, alpha=0.6)
ax4.axhline(ALPHA_NOMINAL, c='C3', ls='--', lw=0.8, alpha=0.6)
ax4.annotate('1550 nm', xy=(1550, 0.182), xytext=(1400, 0.5),
             fontsize=7, ha='center',
             arrowprops=dict(arrowstyle='->', lw=0.8))
ax4.set(xlabel='Wavelength (nm)', ylabel='Attenuation (dB/km)',
        title='D: SMF-28 attenuation spectrum (Keiser Fig 3.2)\nLiterature reference',
        xlim=(750, 1650))
ax4.grid(True, alpha=0.25)
ax4.text(0.98, 0.98, 'Rayleigh\nscattering\nlimit',
         transform=ax4.transAxes, fontsize=7, ha='right', va='top',
         bbox=dict(boxstyle='round', fc='lightblue', alpha=0.3))

# --- Panel E: Attenuation coefficient accuracy at multiple distances ---
# Extract alpha from each segment using paired ratios
alpha_extracted = -10 / distances_km[1:] * np.log10(P_out_meas[1:] / P_in)
# Use scatter plot with mean line
ax5 = fig.add_subplot(gs[1, 1])
ax5.plot(distances_km[1:], alpha_extracted, 'o-', ms=3, c='C0', lw=0.8,
         label=f'Extracted alpha (mean = {alpha_extracted.mean():.6f})')
ax5.axhline(ALPHA_NOMINAL, c='C3', ls='--', lw=1, label=f'Nominal alpha = {ALPHA_NOMINAL}')
ax5.set(xlabel='Distance (km)', ylabel='Extracted alpha (dB/km)',
        title=f'E: Alpha consistency across distance')
ax5.legend(fontsize=8)
ax5.grid(True, alpha=0.25)
ax5.set_ylim(ALPHA_NOMINAL - 0.01*ALPHA_NOMINAL, ALPHA_NOMINAL + 0.01*ALPHA_NOMINAL)

# --- Panel F: tabular summary ---
ax6 = fig.add_subplot(gs[1, 2])
ax6.axis('off')
# Instead of a table, show a concise summary of the validation result
test_points = [0, 20, 50, 100, 150, 200]
row_data = []
for L in test_points:
    idx = np.argmin(np.abs(distances_km - L))
    row_data.append([f'{L}', f'{P_out_meas[idx]:.8e}',
                     f'{P_out_theory[idx]:.8e}',
                     f'{errors_pct[idx]:.6e}'])
col_labels = ['L (km)', 'P_meas (a.u.)', 'P_theory (a.u.)', 'Error (%)']
table = ax6.table(cellText=row_data, colLabels=col_labels,
                  loc='center', cellLoc='center', fontsize=7)
table.auto_set_column_width(col=list(range(len(col_labels))))
table.auto_set_font_size(False)
table.set_fontsize(7)
ax6.set_title('F: Validation summary', fontsize=10, pad=10)

fig.suptitle('Fiber Attenuation - Validation vs Keiser [1] Eq 3.6',
             fontsize=13, fontweight='bold', y=0.98)
fig.savefig(os.path.join(OUT, f'val_attenuation--seed{SEED}.png'), dpi=200, bbox_inches='tight')
print(f"Saved: val_attenuation--seed{SEED}.png")

csv_name = f'val_attenuation--seed{SEED}.csv'
np.savetxt(os.path.join(OUT, csv_name),
           np.column_stack([distances_km, P_out_meas, P_out_theory, errors_pct,
                            loss_dB, alpha_fitted * distances_km]),
           delimiter=',',
           header='distance_km,power_meas,power_theory,error_pct,loss_dB,loss_fit_dB',
           comments='')
print(f"Saved: {csv_name}")
plt.close(fig)
