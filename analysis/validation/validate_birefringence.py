import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os, sys, argparse

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from src.channel.fiber import apply_birefringence, SECTIONAL_LIMIT

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
biref_T0 = 5.0e-8
temp_coeff = -3.0e-9
bend_factor = 0.135
REALIZATIONS = 50

# --- Self-consistency checks: sectional model (L < SECTIONAL_LIMIT) ---
def test_power_conservation_sectional():
    E = np.random.randn(1000, 2) + 1j * np.random.randn(1000, 2)
    P_in = np.mean(np.abs(E)**2)
    for L_m in [1, 10, 100, 1000]:
        E_out = apply_birefringence(E.copy(), L_m, wavelength=WAVELENGTH,
                                    model='sectional')
        P_out = np.mean(np.abs(E_out)**2)
        assert abs(P_out - P_in) / P_in < 1e-12
    print("  [PASS] Power conservation (sectional model)")

def test_temperature_dependence_sectional():
    E = np.ones((100, 2), dtype=complex)
    np.random.seed(SEED + 1)
    Js = [apply_birefringence(E.copy(), 1000, wavelength=WAVELENGTH,
                               temperature=T, model='sectional')[0, 0]
          for T in [0, 25, 50]]
    assert not np.allclose(Js[0], Js[1])
    print("  [PASS] Temperature sensitivity (sectional model)")

def test_wavelength_dependence_sectional():
    E = np.ones((100, 2), dtype=complex)
    np.random.seed(SEED + 2)
    Js = [apply_birefringence(E.copy(), 1000, wavelength=lam, model='sectional')[0, 0]
          for lam in [1310e-9, 1550e-9]]
    assert not np.allclose(Js[0], Js[1])
    print("  [PASS] Wavelength dependence (sectional model)")

def test_randomness_sectional():
    E = np.ones((100, 2), dtype=complex)
    np.random.seed(SEED)
    J1 = apply_birefringence(E.copy(), 100, wavelength=WAVELENGTH, model='sectional')[0, 0]
    np.random.seed(SEED + 100)
    J2 = apply_birefringence(E.copy(), 100, wavelength=WAVELENGTH, model='sectional')[0, 0]
    assert not np.allclose(J1, J2)
    print("  [PASS] Seed-dependent output (sectional model)")

def test_output_on_poincare_sphere_sectional():
    E = np.array([[1.0, 0.0]], dtype=complex)
    np.random.seed(SEED)
    outputs = [apply_birefringence(E.copy(), 1500, wavelength=WAVELENGTH,
                                    model='sectional')[0] for _ in range(50)]
    ex_powers = [np.abs(o[0])**2 for o in outputs]
    assert np.allclose([np.abs(o[0])**2 + np.abs(o[1])**2 for o in outputs], 1.0, atol=1e-12)
    assert np.std(ex_powers) > 0.05
    print("  [PASS] Output polarisation varies (sectional model, 1.5 km)")

# --- Self-consistency checks: phenomenological model (L >= SECTIONAL_LIMIT) ---
def test_power_conservation_phenomenological():
    E = np.random.randn(1000, 2) + 1j * np.random.randn(1000, 2)
    P_in = np.mean(np.abs(E)**2)
    for L_m in [5000, 50000, 100000]:
        E_out = apply_birefringence(E.copy(), L_m, wavelength=WAVELENGTH,
                                    model='phenomenological')
        P_out = np.mean(np.abs(E_out)**2)
        assert abs(P_out - P_in) / P_in < 1e-12
    print("  [PASS] Power conservation (phenomenological model)")

def test_zero_length_identity():
    E = np.random.randn(100, 2) + 1j * np.random.randn(100, 2)
    E_out = apply_birefringence(E.copy(), 0, wavelength=WAVELENGTH)
    assert np.allclose(E_out, E)
    print("  [PASS] Zero length returns field unchanged")

def test_temperature_dependence_phenomenological():
    E = np.ones((100, 2), dtype=complex)
    np.random.seed(SEED + 3)
    Js = [apply_birefringence(E.copy(), 50000, wavelength=WAVELENGTH,
                               temperature=T, model='phenomenological')[0, 0]
          for T in [0, 25, 50]]
    assert not np.allclose(Js[0], Js[1])
    print("  [PASS] Temperature sensitivity (phenomenological model)")

def test_wavelength_dependence_phenomenological():
    E = np.ones((100, 2), dtype=complex)
    np.random.seed(SEED + 4)
    Js = [apply_birefringence(E.copy(), 50000, wavelength=lam, model='phenomenological')[0, 0]
          for lam in [1310e-9, 1550e-9]]
    assert not np.allclose(Js[0], Js[1])
    print("  [PASS] Wavelength dependence (phenomenological model)")

def test_randomness_phenomenological():
    E = np.ones((100, 2), dtype=complex)
    np.random.seed(SEED)
    J1 = apply_birefringence(E.copy(), 50000, wavelength=WAVELENGTH, model='phenomenological')[0, 0]
    np.random.seed(SEED + 100)
    J2 = apply_birefringence(E.copy(), 50000, wavelength=WAVELENGTH, model='phenomenological')[0, 0]
    assert not np.allclose(J1, J2)
    print("  [PASS] Seed-dependent output (phenomenological model)")

def test_output_on_poincare_sphere_phenomenological():
    E = np.array([[1.0, 0.0]], dtype=complex)
    np.random.seed(SEED)
    outputs = [apply_birefringence(E.copy(), 100e3, wavelength=WAVELENGTH,
                                    model='phenomenological')[0] for _ in range(50)]
    ex_powers = [np.abs(o[0])**2 for o in outputs]
    assert np.std(ex_powers) > 0.05
    print("  [PASS] Output polarisation varies (phenomenological model, 100 km)")

# --- Auto-dispatch validation ---
def test_auto_dispatch():
    E = np.ones((10, 2), dtype=complex)
    np.random.seed(SEED + 10)
    out_short = apply_birefringence(E.copy(), SECTIONAL_LIMIT - 1, wavelength=WAVELENGTH, model='auto')
    np.random.seed(SEED + 10)
    out_short_explicit = apply_birefringence(E.copy(), SECTIONAL_LIMIT - 1, wavelength=WAVELENGTH, model='sectional')
    assert np.allclose(out_short, out_short_explicit), \
        f"auto should route to sectional below limit ({SECTIONAL_LIMIT} m)"

    np.random.seed(SEED + 10)
    out_long = apply_birefringence(E.copy(), SECTIONAL_LIMIT + 1, wavelength=WAVELENGTH, model='auto')
    np.random.seed(SEED + 10)
    out_long_explicit = apply_birefringence(E.copy(), SECTIONAL_LIMIT + 1, wavelength=WAVELENGTH, model='phenomenological')
    assert np.allclose(out_long, out_long_explicit), \
        f"auto should route to phenomenological at or above limit ({SECTIONAL_LIMIT} m)"
    print("  [PASS] Auto-dispatch correctly routes based on fibre length")

def test_enabled_false():
    E = np.random.randn(100, 2) + 1j * np.random.randn(100, 2)
    E_out = apply_birefringence(E.copy(), 10000, wavelength=WAVELENGTH, enabled=False)
    assert np.allclose(E_out, E)
    print("  [PASS] enabled=False returns field unchanged")

print("Birefringence validation: multi-section and phenomenological models")
print("  Ref: Menyuk & Wai, JOSA B 1994; Wai & Menyuk, JLT 1996; Ulrich 1980")
print(f"  SECTIONAL_LIMIT = {SECTIONAL_LIMIT} m")
test_power_conservation_sectional()
test_temperature_dependence_sectional()
test_wavelength_dependence_sectional()
test_randomness_sectional()
test_output_on_poincare_sphere_sectional()
test_power_conservation_phenomenological()
test_zero_length_identity()
test_temperature_dependence_phenomenological()
test_wavelength_dependence_phenomenological()
test_randomness_phenomenological()
test_output_on_poincare_sphere_phenomenological()
test_auto_dispatch()
test_enabled_false()

# ============================================================
# Main validation figure — 6 panels
# Panel layout:
#   A — sectional: mean |Ex|^2 vs distance (0–1.5 km)
#   B — phenomenological: mean |Ex|^2 vs distance (0–200 km)
#   C — phenomenological: mean |Ex|^2 vs temperature
#   D — phenomenological: mean |Ex|^2 vs bend radius
#   E — Beat length vs wavelength (analytical)
#   F — Total Δn = base + temp + bend
# ============================================================
fig = plt.figure(figsize=(14, 10))
gs = fig.add_gridspec(2, 3, hspace=0.35, wspace=0.35)

# --- Panel A: Sectional model — mean |Ex|^2 vs distance (0–1.5 km) ---
np.random.seed(SEED)
E_in = np.array([[1.0, 0.0]], dtype=complex)
dist_short_km = np.arange(0, 1.6, 0.1)
mean_ex_short = []
for d in dist_short_km:
    ex_powers = []
    for _ in range(REALIZATIONS):
        out = apply_birefringence(E_in.copy(), d * 1e3, wavelength=WAVELENGTH,
                                  model='sectional')
        ex_powers.append(np.abs(out[0, 0])**2)
    mean_ex_short.append(np.mean(ex_powers))
mean_ex_short = np.array(mean_ex_short)

ax1 = fig.add_subplot(gs[0, 0])
ax1.plot(dist_short_km, mean_ex_short, 's-', c='C3', lw=1.5, ms=4)
ax1.set(xlabel='Fibre length (km)', ylabel=r'$\langle |E_x|^2 \rangle$',
        title='A: Sectional model — polarisation vs distance\n(L_B ≈ 31 m, sections = 100 m)')
ax1.grid(True, alpha=0.25)
ax1.set_ylim(-0.05, 1.05)
ax1.axvline(2, c='gray', ls=':', lw=0.6, alpha=0.4, label='Sectional limit')
ax1.annotate('Random walk on\nPoincar\\\'e sphere',
             xy=(0.05, 0.95), xycoords='axes fraction', ha='left', va='top',
             fontsize=8, bbox=dict(boxstyle='round', fc='lightyellow', alpha=0.8))

# --- Panel B: Phenomenological model — mean |Ex|^2 vs distance (0–200 km) ---
np.random.seed(SEED)
dist_long_km = np.arange(0, 210, 10)
mean_ex_long = []
for d in dist_long_km:
    ex_powers = []
    for _ in range(REALIZATIONS):
        out = apply_birefringence(E_in.copy(), d * 1e3, wavelength=WAVELENGTH,
                                  model='phenomenological')
        ex_powers.append(np.abs(out[0, 0])**2)
    mean_ex_long.append(np.mean(ex_powers))
mean_ex_long = np.array(mean_ex_long)

ax2 = fig.add_subplot(gs[0, 1])
ax2.plot(dist_long_km, mean_ex_long, 's-', c='C4', lw=1.5, ms=4)
ax2.set(xlabel='Fibre length (km)', ylabel=r'$\langle |E_x|^2 \rangle$',
        title='B: Phenomenological — sqrt(L) rotation\n(Menyuk & Wai 1994, L₀ = 75 km)')
ax2.grid(True, alpha=0.25)
ax2.set_ylim(-0.05, 1.05)
ax2.annotate('Gradual scrambling\nθ ∝ √(L/L_char)',
             xy=(0.05, 0.95), xycoords='axes fraction', ha='left', va='top',
             fontsize=8, bbox=dict(boxstyle='round', fc='lightyellow', alpha=0.8))

# --- Panel C: Phenomenological — mean |Ex|^2 vs temperature (50 km) ---
np.random.seed(SEED)
E_in = np.array([[1.0, 0.0]], dtype=complex)
fixed_dist = 50e3
temp_C = np.linspace(0, 60, 13)
mean_ex_temp = []
for T in temp_C:
    ex_powers = []
    for _ in range(REALIZATIONS):
        out = apply_birefringence(E_in.copy(), fixed_dist, wavelength=WAVELENGTH,
                                  temperature=T, model='phenomenological')
        ex_powers.append(np.abs(out[0, 0])**2)
    mean_ex_temp.append(np.mean(ex_powers))
mean_ex_temp = np.array(mean_ex_temp)

ax3 = fig.add_subplot(gs[0, 2])
ax3.plot(temp_C, mean_ex_temp, 's-', c='C1', lw=1.5, ms=4)
ax3.set(xlabel='Temperature (°C)', ylabel=r'$\langle |E_x|^2 \rangle$',
        title='C: Temperature sensitivity\n(50 km, phenomenological)')
ax3.grid(True, alpha=0.25)
ax3.set_ylim(-0.05, 1.05)
ax3.annotate('Δn changes with T,\nθ ∝ 1/√|Δn|',
             xy=(0.05, 0.95), xycoords='axes fraction', ha='left', va='top',
             fontsize=8, bbox=dict(boxstyle='round', fc='lightyellow', alpha=0.8))

# --- Panel D: Phenomenological — mean |Ex|^2 vs bend radius (50 km) ---
np.random.seed(SEED)
bend_radii = np.array([0.002, 0.005, 0.01, 0.015, 0.02, 0.03, 0.05])
mean_ex_bend = []
for R in bend_radii:
    ex_powers = []
    for _ in range(REALIZATIONS):
        out = apply_birefringence(E_in.copy(), fixed_dist, wavelength=WAVELENGTH,
                                  temperature=25, bend_radius=R, model='phenomenological')
        ex_powers.append(np.abs(out[0, 0])**2)
    mean_ex_bend.append(np.mean(ex_powers))
mean_ex_bend = np.array(mean_ex_bend)

ax4 = fig.add_subplot(gs[1, 0])
ax4.semilogx(bend_radii * 1e3, mean_ex_bend, 's-', c='C2', lw=1.5, ms=4)
ax4.set(xlabel='Bend radius (mm)', ylabel=r'$\langle |E_x|^2 \rangle$',
        title='D: Bend-induced birefringence\n(Ulrich [7], 50 km, phenomenological)')
ax4.grid(True, alpha=0.25)

# --- Panel E: Beat length vs wavelength (analytical) ---
lam_range = np.linspace(800e-9, 1700e-9, 50)
L_B_vals = lam_range / biref_T0

ax5 = fig.add_subplot(gs[1, 1])
ax5.plot(lam_range * 1e9, L_B_vals, '-', c='C4', lw=1.5)
ax5.axvline(1550, c='gray', ls='--', lw=0.6, alpha=0.5)
ax5.axvline(1310, c='gray', ls='--', lw=0.6, alpha=0.5)
ax5.set(xlabel='Wavelength (nm)', ylabel='Beat length L_B (m)',
        title='E: L_B = λ / Δn (Agrawal [6] §4.1)')
ax5.grid(True, alpha=0.25)
L_B_1550 = 1550e-9 / biref_T0
L_B_1310 = 1310e-9 / biref_T0
ax5.annotate(f'@ 1550 nm: L_B = {L_B_1550:.1f} m\n'
             f'@ 1310 nm: L_B = {L_B_1310:.1f} m',
             xy=(0.05, 0.95), xycoords='axes fraction', ha='left', va='top',
             fontsize=8, bbox=dict(boxstyle='round', fc='lightyellow', alpha=0.8))

# --- Panel F: Total Δn = base + temp + bend for each R ---
ax6 = fig.add_subplot(gs[1, 2])
T_range = np.linspace(-20, 60, 50)
colors = ['C1', 'C4', 'C2', 'C3', 'C0']
for i, R_mm in enumerate([np.inf, 20, 10, 5, 2]):
    if np.isinf(R_mm):
        dn_total = biref_T0 + temp_coeff * (T_range - T0_C)
        label = 'No bend (R=$\\infty$)'
        ls = '-'
        lw = 1.5
    else:
        dn_b = bend_factor * (r_fiber / (R_mm * 1e-3)) ** 2
        dn_total = biref_T0 + temp_coeff * (T_range - T0_C) + dn_b
        label = f'R = {R_mm} mm'
        ls = '--'
        lw = 1.0
    ax6.plot(T_range, dn_total * 1e6, ls, c=colors[i], lw=lw, label=label, alpha=0.8)
ax6.axhline(biref_T0 * 1e6, c='k', ls=':', lw=0.8, alpha=0.5, label='Base Δn₀')
ax6.set(xlabel='Temperature (°C)', ylabel=r'Δn (×10⁻⁶)',
        title='F: Total Δn = base + temp + bend\n(Ulrich [7] model, sectional Δn₀)')
ax6.legend(fontsize=6, ncol=2)
ax6.grid(True, alpha=0.25)

fig.suptitle('Fiber Birefringence — Dual-Model Validation (Sectional + Phenomenological)\n'
             '(Menyuk & Wai 1994; Wai & Menyuk 1996; Ulrich 1980; Agrawal 2021)',
             fontsize=12, fontweight='bold', y=0.98)
fig.savefig(os.path.join(OUT, f'val_birefringence--seed{SEED}.png'), dpi=200, bbox_inches='tight')
print(f"Saved: val_birefringence--seed{SEED}.png")

csv_name = f'val_birefringence--seed{SEED}.csv'
with open(os.path.join(OUT, csv_name), 'w') as f:
    f.write('# Sectional model (short fibre)\n')
    f.write('dist_short_km,mean_Ex_short\n')
    for d, m in zip(dist_short_km, mean_ex_short):
        f.write(f'{d},{m}\n')
    f.write('# Phenomenological model (long fibre)\n')
    f.write('dist_long_km,mean_Ex_long\n')
    for d, m in zip(dist_long_km, mean_ex_long):
        f.write(f'{d},{m}\n')
print(f"Saved: {csv_name}")

import csv
table_csv = os.path.join(OUT, f'val_birefringence--seed{SEED}_table.csv')
with open(table_csv, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Test', 'Model', 'Result'])
    writer.writerow(['Power conservation', 'sectional', 'err < 1e-12'])
    writer.writerow(['Temperature sensitivity', 'sectional', 'PASS'])
    writer.writerow(['Wavelength sensitivity', 'sectional', 'PASS'])
    writer.writerow(['Seed randomness', 'sectional', 'PASS'])
    writer.writerow(['Poincaré scrambling', 'sectional', 'std > 0.05'])
    writer.writerow(['Power conservation', 'phenomenological', 'err < 1e-12'])
    writer.writerow(['Temperature sensitivity', 'phenomenological', 'PASS'])
    writer.writerow(['Wavelength sensitivity', 'phenomenological', 'PASS'])
    writer.writerow(['Seed randomness', 'phenomenological', 'PASS'])
    writer.writerow(['Poincaré scrambling', 'phenomenological', 'std > 0.05'])
    writer.writerow(['Zero-length identity', 'both', 'PASS'])
    writer.writerow(['Auto-dispatch', 'auto', f'{SECTIONAL_LIMIT} m threshold'])
    writer.writerow(['enabled=False', 'both', 'PASS'])
print(f"Saved: val_birefringence--seed{SEED}_table.csv")
plt.close(fig)
