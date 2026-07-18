import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats as sp_stats
import os, sys, argparse

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from src.channel.fiber import apply_birefringence

OUT = os.path.join(os.path.dirname(__file__), '..', 'val_birefringence')
os.makedirs(OUT, exist_ok=True)

parser = argparse.ArgumentParser()
parser.add_argument('--seed', type=int, default=42)
args = parser.parse_args()
SEED = args.seed
np.random.seed(SEED)

WAVELENGTH = 1550e-9
r_fiber = 62.5e-6
T0_C = 25.0
biref_T0 = 0.87e-5
temp_coeff = -5e-7
bend_factor = 0.135

def dbeta(biref, lam=WAVELENGTH):
    return 2 * np.pi * biref / lam

# --- Self-consistency checks ---
def test_power_conservation():
    E = np.random.randn(1000, 2) + 1j * np.random.randn(1000, 2)
    P_in = np.mean(np.abs(E)**2)
    for L_m in [1, 10, 100, 1000]:
        E_out = apply_birefringence(E.copy(), L_m, wavelength=WAVELENGTH)
        P_out = np.mean(np.abs(E_out)**2)
        assert abs(P_out - P_in) / P_in < 1e-12
    print("  [PASS] Power conservation (unitary Jones matrix)")

def test_phase_ratio():
    E = np.ones((2, 2), dtype=complex)
    J0 = [apply_birefringence(E.copy(), L, wavelength=WAVELENGTH)[0, 0] for L in [1, 2]]
    ratio = J0[1] * np.conj(J0[0])
    expected = np.exp(1j * dbeta(biref_T0) * 0.5)
    assert np.abs(ratio - expected) < 1e-12
    print("  [PASS] Phase ratio: complex error < 1e-12")

def test_temperature_dependence():
    E = np.ones((100, 2), dtype=complex)
    Js = [apply_birefringence(E.copy(), 1000, wavelength=WAVELENGTH, temperature=T)[0, 0]
          for T in [0, 25, 50]]
    assert not np.allclose(Js[0], Js[1]) and not np.allclose(Js[1], Js[2])
    print("  [PASS] Temperature sensitivity detected")

def test_wavelength_dependence():
    E = np.ones((100, 2), dtype=complex)
    Js = [apply_birefringence(E.copy(), 1000, wavelength=lam)[0, 0]
          for lam in [1310e-9, 1550e-9]]
    assert not np.allclose(Js[0], Js[1])
    print("  [PASS] Wavelength dependence detected")

print("Birefringence validation via apply_birefringence")
test_power_conservation()
test_phase_ratio()
test_temperature_dependence()
test_wavelength_dependence()

# ============================================================
# Main validation figure — 6 panels
# Panels A/B/C: theory and model overlaid (short L, no wrap)
# ============================================================
fig = plt.figure(figsize=(14, 10))
gs = fig.add_gridspec(2, 3, hspace=0.35, wspace=0.35)

# --- Panel A: Phase vs length (short L, no wrapping) ---
L_short = np.linspace(0, 0.5, 100)  # max phase ~9 rad, no wrap
E_in = np.ones((1, 2), dtype=complex)
dbeta_val = dbeta(biref_T0)

J_sim_A = np.array([apply_birefringence(E_in.copy(), L, wavelength=WAVELENGTH)[0, 0]
                     for L in L_short])
J_ana_A = np.exp(1j * dbeta_val * L_short / 2)
phase_sim_A = np.angle(J_sim_A)
phase_ana_A = np.angle(J_ana_A)

ax1 = fig.add_subplot(gs[0, 0])
ax1.plot(L_short, phase_ana_A, '-', c='k', lw=2, label='Theory: dbeta*L/2')
ax1.plot(L_short[::5], phase_sim_A[::5], 'o', ms=4, c='C3', label='apply_birefringence')
ax1.set(xlabel='Fibre length (m)', ylabel='Phase of J[0,0] (rad)',
        title='A: Phase vs length (Agrawal [6] Eq 4.1.2)')
ax1.legend(fontsize=8)
ax1.grid(True, alpha=0.25)
ax1.annotate(f'delta_beta = {dbeta_val:.1f} rad/m\nL_B = {2*np.pi/dbeta_val*1e3:.1f} mm\n'
             f'Phase error max = {np.abs(phase_sim_A - phase_ana_A).max():.2e} rad',
             xy=(0.05, 0.95), xycoords='axes fraction', ha='left', va='top',
             fontsize=8, bbox=dict(boxstyle='round', fc='lightyellow', alpha=0.8))

# --- Panel B: Temperature sensitivity (short L, no wrapping) ---
temp_C = np.linspace(-20, 60, 33)
L_temp = 0.1  # max phase ~6 rad, no wrap
E_in = np.ones((1, 2), dtype=complex)

J_sim_B = np.array([apply_birefringence(E_in.copy(), L_temp, wavelength=WAVELENGTH,
                                          temperature=T)[0, 0]
                     for T in temp_C])
J_ana_B = np.array([np.exp(1j * dbeta(biref_T0 + temp_coeff * (T - T0_C)) * L_temp / 2)
                     for T in temp_C])
phase_sim_B = np.angle(J_sim_B)
phase_ana_B = np.angle(J_ana_B)

ax2 = fig.add_subplot(gs[0, 1])
ax2.plot(temp_C, phase_ana_B, '-', c='k', lw=2, label='Theory: dbeta(T)*L/2')
ax2.plot(temp_C[::3], phase_sim_B[::3], 'o', ms=4, c='C1', label='apply_birefringence')
ax2.set(xlabel='Temperature (C)', ylabel='Phase of J[0,0] (rad)',
        title='B: Temperature sensitivity')
ax2.legend(fontsize=8)
ax2.grid(True, alpha=0.25)
ax2.annotate(f'temp_coeff = -5e-7 /C\nPhase error max = {np.abs(phase_sim_B - phase_ana_B).max():.2e} rad',
             xy=(0.05, 0.05), xycoords='axes fraction', ha='left', va='bottom',
             fontsize=8, bbox=dict(boxstyle='round', fc='lightyellow', alpha=0.8))

# --- Panel C: Bend radius sweep (short L, no wrapping) ---
bend_radii = np.logspace(np.log10(0.002), np.log10(0.02), 40)
L_bend = 0.01  # max phase ~2.9 rad at R=2mm, no wrap
E_in = np.ones((1, 2), dtype=complex)

J_sim_C = np.array([apply_birefringence(E_in.copy(), L_bend, wavelength=WAVELENGTH,
                                          temperature=T0_C, bend_radius=R)[0, 0]
                     for R in bend_radii])
dn_bend = bend_factor * (r_fiber / bend_radii)**2
J_ana_C = np.exp(1j * dbeta(biref_T0 + dn_bend) * L_bend / 2)
phase_sim_C = np.angle(J_sim_C)
phase_ana_C = np.angle(J_ana_C)

ax3 = fig.add_subplot(gs[0, 2])
ax3.plot(bend_radii*1e3, phase_ana_C, '-', c='k', lw=2, label='Theory: Ulrich [7]')
ax3.plot(bend_radii[::3]*1e3, phase_sim_C[::3], 'o', ms=4, c='C2', label='apply_birefringence')
ax3.set(xlabel='Bend radius (mm)', ylabel='Phase of J[0,0] (rad)',
        title='C: Bend-induced birefringence\n(Ulrich [7] Eq 1)')
ax3.legend(fontsize=8)
ax3.grid(True, alpha=0.25)
ax3.annotate(f'Phase error max = {np.abs(phase_sim_C - phase_ana_C).max():.2e} rad',
             xy=(0.95, 0.95), xycoords='axes fraction', ha='right', va='top',
             fontsize=9, fontweight='bold',
             bbox=dict(boxstyle='round', fc='lightgreen', alpha=0.6))

# --- Panel D: Extracted vs expected delta_n ---
L_ref = 10.0
delta_L = 0.001
E_in = np.ones((1, 2), dtype=complex)

J_ref = np.array([apply_birefringence(E_in.copy(), L_ref, wavelength=WAVELENGTH,
                                       temperature=T0_C, bend_radius=R)[0, 0]
                   for R in bend_radii])
J_del = np.array([apply_birefringence(E_in.copy(), L_ref + delta_L, wavelength=WAVELENGTH,
                                       temperature=T0_C, bend_radius=R)[0, 0]
                   for R in bend_radii])

dn_extracted = np.angle(J_del * np.conj(J_ref)) * WAVELENGTH / (np.pi * delta_L)
dn_expected = biref_T0 + dn_bend

ax4 = fig.add_subplot(gs[1, 0])
ax4.loglog(bend_radii*1e3, dn_extracted, 'o', ms=4, c='C2', label='Extracted delta_n')
ax4.loglog(bend_radii*1e3, dn_expected, '--k', lw=1.5, label='Ulrich [7] theory')
ax4.set(xlabel='Bend radius (mm)', ylabel='Total delta_n',
        title='D: delta_n vs bend radius\n(Ulrich [7] Eq 1)')
ax4.legend(fontsize=7)
ax4.grid(True, alpha=0.25)
dn_err = np.abs(dn_extracted - dn_expected) / dn_expected * 100
ax4.annotate(f'max error: {dn_err.max():.4f}%',
             xy=(0.05, 0.95), xycoords='axes fraction', ha='left', va='top',
             fontsize=9, fontweight='bold',
             bbox=dict(boxstyle='round', fc='lightgreen', alpha=0.6))

# --- Panel E: Beat length vs wavelength ---
lam_range = np.linspace(800e-9, 1700e-9, 50)
L_B_vals = lam_range / biref_T0

ax5 = fig.add_subplot(gs[1, 1])
ax5.plot(lam_range*1e9, L_B_vals*1e3, '-', c='C4', lw=1.5)
ax5.axvline(1550, c='gray', ls='--', lw=0.6, alpha=0.5)
ax5.axvline(1310, c='gray', ls='--', lw=0.6, alpha=0.5)
ax5.set(xlabel='Wavelength (nm)', ylabel='Beat length L_B (mm)',
        title='E: Beat length vs wavelength\nL_B = lambda / delta_n (Agrawal [6])')
ax5.grid(True, alpha=0.25)
L_B_1550 = 1550e-9 / biref_T0 * 1e3
L_B_1310 = 1310e-9 / biref_T0 * 1e3
ax5.annotate(f'@ 1550 nm: L_B = {L_B_1550:.2f} mm\n'
             f'@ 1310 nm: L_B = {L_B_1310:.2f} mm',
             xy=(0.05, 0.95), xycoords='axes fraction', ha='left', va='top',
             fontsize=8, bbox=dict(boxstyle='round', fc='lightyellow', alpha=0.8))

# --- Panel F: Validation summary ---
ax6 = fig.add_subplot(gs[1, 2])
ax6.axis('off')
summary = [
    ['Power conservation', 'Unitary Jones', 'err < 1e-12'],
    ['Phase vs length', f'dbeta = {dbeta_val:.1f} rad/m',
     f'phase err = {np.abs(phase_sim_A - phase_ana_A).max():.2e}'],
    ['Temperature', f'L = {L_temp} m',
     f'phase err = {np.abs(phase_sim_B - phase_ana_B).max():.2e}'],
    ['Bend (Ulrich [7])', f'L = {L_bend} m',
     f'phase err = {np.abs(phase_sim_C - phase_ana_C).max():.2e}'],
    ['Wavelength', '1310 vs 1550 nm differ', 'PASS'],
]
table = ax6.table(cellText=summary,
                  colLabels=['Test', 'Metric', 'Result'],
                  loc='center', cellLoc='center', fontsize=7)
table.auto_set_column_width(col=list(range(3)))
table.auto_set_font_size(False)
table.set_fontsize(7)
for (row, col), cell in table.get_celld().items():
    if row == 0:
        cell.set_facecolor('#40466e')
        cell.set_text_props(color='w', fontweight='bold')
ax6.set_title('F: Validation summary', fontsize=10, pad=10)

fig.suptitle('Fiber Birefringence - Validation vs Agrawal [6] Eq 4.1.2 & Ulrich [7]',
             fontsize=13, fontweight='bold', y=0.98)
fig.savefig(os.path.join(OUT, f'val_birefringence--seed{SEED}.png'), dpi=200, bbox_inches='tight')
print(f"Saved: val_birefringence--seed{SEED}.png")

csv_name = f'val_birefringence--seed{SEED}.csv'
np.savetxt(os.path.join(OUT, csv_name),
           np.column_stack([L_short, phase_sim_A, phase_ana_A]),
           delimiter=',', header='length_m,phase_sim_rad,phase_analytic_rad', comments='')
print(f"Saved: {csv_name}")
plt.close(fig)
