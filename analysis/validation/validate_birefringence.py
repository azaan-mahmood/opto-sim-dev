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

def dbeta_from_biref(biref, lam=WAVELENGTH):
    return 2 * np.pi * biref / lam

def jones_error(J_sim, J_ana):
    return np.abs(J_sim - J_ana)

# --- Self-consistency checks ---
def test_power_conservation():
    E = np.random.randn(1000, 2) + 1j * np.random.randn(1000, 2)
    P_in = np.mean(np.abs(E)**2)
    for L_m in [1, 10, 100, 1000]:
        E_out = apply_birefringence(E.copy(), L_m, wavelength=WAVELENGTH)
        P_out = np.mean(np.abs(E_out)**2)
        err = abs(P_out - P_in) / P_in
        assert err < 1e-12, f"Power not conserved at L={L_m}m"
    print("  [PASS] Power conservation (unitary Jones matrix)")

def test_phase_ratio():
    E = np.ones((2, 2), dtype=complex)
    J0 = [apply_birefringence(E.copy(), L, wavelength=WAVELENGTH)[0, 0] for L in [1, 2]]
    ratio = J0[1] * np.conj(J0[0])
    expected = np.exp(1j * dbeta_from_biref(biref_T0) * 0.5)
    err = np.abs(ratio - expected)
    assert err < 1e-12, f"Phase ratio error: {err:.2e}"
    print(f"  [PASS] Phase ratio: complex error = {err:.2e}")

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
# Main validation figure — 6 panels (all error-based, no unwrap)
# ============================================================
fig = plt.figure(figsize=(14, 10))
gs = fig.add_gridspec(2, 3, hspace=0.35, wspace=0.35)

# --- Panel A: Jones matrix error vs length ---
L_m_dense = np.logspace(0, 4, 200)
E_in = np.ones((1, 2), dtype=complex)
dbeta_val = dbeta_from_biref(biref_T0)

J_sim_arr = np.array([
    apply_birefringence(E_in.copy(), L, wavelength=WAVELENGTH)[0, 0]
    for L in L_m_dense
])
J_ana_arr = np.exp(1j * dbeta_val * L_m_dense / 2)
err_A = jones_error(J_sim_arr, J_ana_arr)

ax1 = fig.add_subplot(gs[0, 0])
ax1.plot(L_m_dense, err_A, '-', c='C0', lw=1.5)
ax1.set(xlabel='Fibre length (m)', ylabel='|J_sim - J_analytic|',
        title='A: Jones matrix error vs length\n(Agrawal [6] Eq 4.1.2)',
        ylim=(-1e-16, 2e-16))
ax1.axhline(0, color='gray', ls=':', lw=0.5)
ax1.grid(True, alpha=0.25)
ax1.annotate(f'max error = {err_A.max():.2e}\ndelta_beta = {dbeta_val:.1f} rad/m\n'
             f'L_B = {2*np.pi/dbeta_val*1e3:.1f} mm',
             xy=(0.95, 0.05), xycoords='axes fraction', ha='right', va='bottom',
             fontsize=8, bbox=dict(boxstyle='round', fc='lightyellow', alpha=0.8))

# --- Panel B: Temperature — Jones matrix error at each T ---
# Compare J_sim(T) with J_analytic(T) directly (no ratio, no unwrap)
temp_C = np.linspace(-20, 60, 33)
E_in = np.ones((1, 2), dtype=complex)
Js_T = np.array([
    apply_birefringence(E_in.copy(), 1000, wavelength=WAVELENGTH, temperature=T)[0, 0]
    for T in temp_C
])

# Analytic Jones element at each temperature
J_ana_T = np.array([
    np.exp(1j * dbeta_from_biref(biref_T0 + temp_coeff * (T - T0_C)) * 1000 / 2)
    for T in temp_C
])

temp_err = jones_error(Js_T, J_ana_T)

ax2 = fig.add_subplot(gs[0, 1])
ax2.plot(temp_C, temp_err, 'o-', ms=3, c='C1', lw=1)
ax2.axhline(0, color='gray', ls=':', lw=0.5)
ax2.set(xlabel='Temperature (C)', ylabel='|J_sim - J_analytic|',
        title='B: Temperature sensitivity\n(Jones matrix error)',
        ylim=(-1e-16, 2e-16))
ax2.grid(True, alpha=0.25)
ax2.annotate(f'max error = {temp_err.max():.2e}\ntemp_coeff = -5e-7 /C',
             xy=(0.05, 0.05), xycoords='axes fraction', ha='left', va='bottom',
             fontsize=8, bbox=dict(boxstyle='round', fc='lightyellow', alpha=0.8))

# --- Panel C: Bend sweep — Jones matrix error ---
bend_radii = np.logspace(np.log10(0.002), np.log10(0.02), 40)
L_bend = 10.0
E_in = np.ones((1, 2), dtype=complex)

J_sim_B = np.array([
    apply_birefringence(E_in.copy(), L_bend, wavelength=WAVELENGTH,
                        temperature=T0_C, bend_radius=R)[0, 0]
    for R in bend_radii
])
dn_bend = bend_factor * (r_fiber / bend_radii)**2
J_ana_B = np.exp(1j * dbeta_from_biref(biref_T0 + dn_bend) * L_bend / 2)
err_C = jones_error(J_sim_B, J_ana_B)

ax3 = fig.add_subplot(gs[0, 2])
ax3.plot(bend_radii*1e3, err_C, 'o-', ms=3, c='C2', lw=1.2)
ax3.axhline(0, color='gray', ls=':', lw=0.5)
ax3.set(xlabel='Bend radius (mm)', ylabel='|J_sim - J_analytic|',
        title='C: Bend-induced birefringence\n(Ulrich [7] Eq 1)',
        ylim=(-1e-16, 2e-16))
ax3.grid(True, alpha=0.25)
ax3.annotate(f'max error = {err_C.max():.2e}',
             xy=(0.95, 0.95), xycoords='axes fraction', ha='right', va='top',
             fontsize=9, fontweight='bold',
             bbox=dict(boxstyle='round', fc='lightgreen', alpha=0.6))

# --- Panel D: Extracted vs expected delta_n using tiny delta_L ---
# Use delta_L = 0.001m (1mm) so phase diff < pi even at R=2mm
# At R=2mm: dbeta = 571 rad/m, delta_L=0.001 => dphi = 571*0.001/2 = 0.286 rad < pi
L_ref = 10.0
delta_L = 0.001
dn_total_expected = biref_T0 + dn_bend

J_ref = np.array([
    apply_birefringence(E_in.copy(), L_ref, wavelength=WAVELENGTH,
                        temperature=T0_C, bend_radius=R)[0, 0]
    for R in bend_radii
])
J_del = np.array([
    apply_birefringence(E_in.copy(), L_ref + delta_L, wavelength=WAVELENGTH,
                        temperature=T0_C, bend_radius=R)[0, 0]
    for R in bend_radii
])

ratio_J = J_del * np.conj(J_ref)
dphi_extracted = np.angle(ratio_J)
dn_extracted = dphi_extracted * WAVELENGTH / (np.pi * delta_L)

ax4 = fig.add_subplot(gs[1, 0])
ax4.loglog(bend_radii*1e3, dn_extracted, 'o', ms=4, c='C2', label='Extracted delta_n')
ax4.loglog(bend_radii*1e3, dn_total_expected, '--k', lw=1.5, label='Ulrich [7] theory')
ax4.set(xlabel='Bend radius (mm)', ylabel='Total delta_n',
        title='D: delta_n vs bend radius\n(Ulrich [7] Eq 1)')
ax4.legend(fontsize=7)
ax4.grid(True, alpha=0.25)
dn_err = np.abs(dn_extracted - dn_total_expected) / dn_total_expected * 100
ax4.annotate(f'max error: {dn_err.max():.4f}%',
             xy=(0.05, 0.95), xycoords='axes fraction', ha='left', va='top',
             fontsize=9, fontweight='bold',
             bbox=dict(boxstyle='round', fc='lightgreen', alpha=0.6))

# --- Panel E: Beat length L_B vs wavelength ---
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
     f'max err = {err_A.max():.2e}'],
    ['Temperature', f'dphi/dT per step',
     f'max err = {temp_err.max():.2e}'],
    ['Bend (Ulrich [7])', f'max bend err = {err_C.max():.2e}',
     f'dn err = {dn_err.max():.4f}%'],
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
           np.column_stack([L_m_dense, err_A]),
           delimiter=',', header='length_m,jones_error', comments='')
print(f"Saved: {csv_name}")
plt.close(fig)
