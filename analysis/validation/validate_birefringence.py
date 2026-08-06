import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os, sys, argparse

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from src.channel.fiber import apply_birefringence
from src.visualization.stokes import compute_stokes_parameters

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

# --- Self-consistency checks moved to the test suite (REPRO-4) ---
# The 13 internal invariants (power conservation, temperature/wavelength
# sensitivity, seed dependence, auto-dispatch, enabled=False, ...) used to
# live here, printed as [PASS] lines. They are unit-test material, not
# literature validation: they now live in tests/test_fiber.py
# (TestBirefringenceSelfConsistency, 13 tests). The literature comparisons
# are TestBirefringenceDepolarization (Menyuk & Wai [10] p^N law) and
# TestUlrichBendLaw (Ulrich [7] Eq. 1) in the same file.

# --- Model note (PHYS-5, 5th pass) ---
# The former phenomenological model (single SU(2) rotation with fitted
# theta = min(pi, sqrt(L/L_char)*pi/2)) was removed: the multi-section
# model is quasi-static, converges to uniform SU(2) within a few hundred
# metres, and costs ~0.16 ms/apply at 122 km. 'auto' == 'sectional' at all
# lengths. Temperature/bend dependence is only visible in the
# single-section regime (L ~ metres); beyond ~1 km the ensemble mean is
# already the uniform-SU(2) value 1/3 + (2/3)cos^2(theta/2) -> 1/2.

print("Birefringence validation: multi-section model (all lengths)")
print("  Ref: Menyuk & Wai, JOSA B 1994; Wai & Menyuk, JLT 1996; Ulrich 1980")
print("  Self-consistency checks: see tests/test_fiber.py "
      "(TestBirefringenceSelfConsistency)")

# ============================================================
# Main validation figure — 6 panels
# Panel layout:
#   A — multi-section: mean |Ex|^2 vs distance (0–1.5 km)
#   B — multi-section: mean |Ex|^2 vs distance (0–200 km, ensemble)
#   C — multi-section: mean |Ex|^2 vs temperature (2 m, single section)
#   D — multi-section: mean |Ex|^2 vs bend radius (2 m, single section)
#   E — Beat length vs wavelength (analytical)
#   F — Total Δn = base + temp + bend
# ============================================================
fig = plt.figure(figsize=(14, 10))
gs = fig.add_gridspec(2, 3, hspace=0.35, wspace=0.35)

# --- Panel A: mean |Ex|^2 vs distance (0–1.5 km) ---
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
        title='A: Multi-section — polarisation vs distance\n(L_B ≈ 31 m, correlation length = 50 m)')
ax1.grid(True, alpha=0.25)
ax1.set_ylim(-0.05, 1.05)
ax1.annotate('Random walk on\nPoincar\\\'e sphere',
             xy=(0.05, 0.95), xycoords='axes fraction', ha='left', va='top',
             fontsize=8, bbox=dict(boxstyle='round', fc='lightyellow', alpha=0.8))

# --- Panel B: mean |Ex|^2 vs distance (0–200 km, ensemble) ---
# With L_c = 50 m and per-section retardance delta ~ 10 rad, the ensemble
# mean reaches its uniform-SU(2) value (1/2) within a few hundred metres —
# the honest long-distance behaviour of the single model (PHYS-5).
np.random.seed(SEED)
dist_long_km = np.arange(0, 210, 10)
mean_ex_long = []
for d in dist_long_km:
    ex_powers = []
    for _ in range(REALIZATIONS):
        out = apply_birefringence(E_in.copy(), d * 1e3, wavelength=WAVELENGTH,
                                  model='sectional')
        ex_powers.append(np.abs(out[0, 0])**2)
    mean_ex_long.append(np.mean(ex_powers))
mean_ex_long = np.array(mean_ex_long)

ax2 = fig.add_subplot(gs[0, 1])
ax2.plot(dist_long_km, mean_ex_long, 's-', c='C4', lw=1.5, ms=4)
ax2.set(xlabel='Fibre length (km)', ylabel=r'$\langle |E_x|^2 \rangle$',
        title='B: Multi-section — long-distance ensemble\n(0–200 km)')
ax2.grid(True, alpha=0.25)
ax2.set_ylim(-0.05, 1.05)
ax2.annotate('Uniform SU(2)\nplateau: |E_x|² → ½\nwithin ~200 m',
             xy=(0.05, 0.95), xycoords='axes fraction', ha='left', va='top',
             fontsize=8, bbox=dict(boxstyle='round', fc='lightyellow', alpha=0.8))

# --- Panel C: mean |Ex|^2 vs temperature (2 m, single section) ---
# At 2 m the fibre is one correlation cell: delta(T) = 2*pi*|delta_n(T)|*L/lambda
# is in the visible regime (delta ~ 0.3–0.5 rad), so the temperature
# coefficient -3.0e-9 / C shows up. Beyond ~1 km the ensemble is saturated
# and temperature is invisible in the mean.
np.random.seed(SEED)
E_in = np.array([[1.0, 0.0]], dtype=complex)
fixed_dist = 2.0  # m — single correlation cell
temp_C = np.linspace(0, 60, 13)
mean_ex_temp = []
for T in temp_C:
    ex_powers = []
    for _ in range(REALIZATIONS):
        out = apply_birefringence(E_in.copy(), fixed_dist, wavelength=WAVELENGTH,
                                  temperature=T, model='sectional')
        ex_powers.append(np.abs(out[0, 0])**2)
    mean_ex_temp.append(np.mean(ex_powers))
mean_ex_temp = np.array(mean_ex_temp)

ax3 = fig.add_subplot(gs[0, 2])
ax3.plot(temp_C, mean_ex_temp, 's-', c='C1', lw=1.5, ms=4)
ax3.set(xlabel='Temperature (°C)', ylabel=r'$\langle |E_x|^2 \rangle$',
        title='C: Temperature sensitivity\n(2 m, single correlation cell)')
ax3.grid(True, alpha=0.25)
ax3.set_ylim(-0.05, 1.05)
ax3.annotate('Δn(T) = 5e-8 − 3e-9·(T−25)\nvisible where δ(T) < π',
             xy=(0.05, 0.95), xycoords='axes fraction', ha='left', va='top',
             fontsize=8, bbox=dict(boxstyle='round', fc='lightyellow', alpha=0.8))

# --- Panel D: mean |Ex|^2 vs bend radius (2 m, single section) ---
np.random.seed(SEED)
bend_radii = np.array([0.001, 0.0015, 0.002, 0.003, 0.005, 0.01, 0.02])
mean_ex_bend = []
for R in bend_radii:
    ex_powers = []
    for _ in range(REALIZATIONS):
        out = apply_birefringence(E_in.copy(), fixed_dist, wavelength=WAVELENGTH,
                                  temperature=25, bend_radius=R, model='sectional')
        ex_powers.append(np.abs(out[0, 0])**2)
    mean_ex_bend.append(np.mean(ex_powers))
mean_ex_bend = np.array(mean_ex_bend)

ax4 = fig.add_subplot(gs[1, 0])
ax4.semilogx(bend_radii * 1e3, mean_ex_bend, 's-', c='C2', lw=1.5, ms=4)
ax4.set(xlabel='Bend radius (mm)', ylabel=r'$\langle |E_x|^2 \rangle$',
        title='D: Bend-induced birefringence\n(Ulrich [7], 2 m, single correlation cell)')
ax4.grid(True, alpha=0.25)
ax4.annotate('Δn_bend = 0.135·(r/R)²\nwraps mod 2π below ~1.5 mm',
             xy=(0.05, 0.95), xycoords='axes fraction', ha='left', va='top',
             fontsize=8, bbox=dict(boxstyle='round', fc='lightyellow', alpha=0.8))

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

fig.suptitle('Fiber Birefringence — Multi-Section Model Validation (single model, all lengths)\n'
             '(Menyuk & Wai 1994; Wai & Menyuk 1996; Ulrich 1980; Agrawal 2021)',
             fontsize=12, fontweight='bold', y=0.98)
fig.savefig(os.path.join(OUT, f'val_birefringence--seed{SEED}.png'), dpi=200, bbox_inches='tight')
print(f"Saved: val_birefringence--seed{SEED}.png")

csv_name = f'val_birefringence--seed{SEED}.csv'
with open(os.path.join(OUT, csv_name), 'w') as f:
    f.write('# Multi-section model (short fibre)\n')
    f.write('dist_short_km,mean_Ex_short\n')
    for d, m in zip(dist_short_km, mean_ex_short):
        f.write(f'{d},{m}\n')
    f.write('# Multi-section model (long fibre)\n')
    f.write('dist_long_km,mean_Ex_long\n')
    for d, m in zip(dist_long_km, mean_ex_long):
        f.write(f'{d},{m}\n')
print(f"Saved: {csv_name}")

import csv
table_csv = os.path.join(OUT, f'val_birefringence--seed{SEED}_table.csv')
with open(table_csv, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Test', 'Model', 'Result'])
    writer.writerow(['Panel figure (A-F)', 'sectional', f'saved val_birefringence--seed{SEED}.png'])
    writer.writerow(['Poincaré convergence', 'sectional', '|mean(S)|: first -> last'])
    writer.writerow(['Self-consistency (13)', 'sectional', 'moved to tests/test_fiber.py'])
print(f"Saved: {table_csv}")
plt.close(fig)

# ============================================================
# Convergence figure: Poincaré sphere scatter plots for the
# multi-section (sectional) model at increasing distances,
# showing the transition from ordered → uniform SU(2).
#
# With the physical correlation length L_c = 50 m (PHYS-4) exceeding the
# beat length L_B ~ 31 m, a fibre shorter than L_c is a *single* random-axis
# section whose retardation phase (alpha = 2*pi*|Δn|*L/lambda) already
# wraps past 2*pi within a few metres — so |mean(S)| oscillates with L
# rather than decaying smoothly until L is a few correlation lengths.
# The first distance below is chosen well under one beat length, where
# that oscillation hasn't kicked in yet, to still show a coherent state.
# ============================================================
from mpl_toolkits.mplot3d import Axes3D

np.random.seed(SEED)
E_in = np.array([[1.0, 0.0]], dtype=complex)
N_REAL = 500
conv_distances = [0.002, 0.05, 1.0, 10.0]  # km

fig2 = plt.figure(figsize=(12, 10))
gs2 = fig2.add_gridspec(2, 2, hspace=0.30, wspace=0.30)
titles = [
    f'A: L = {conv_distances[0]*1e3:.0f} m — single correlation cell, clustered',
    f'B: L = {conv_distances[1]*1e3:.0f} m — one correlation length, partial scrambling',
    f'C: L = {conv_distances[2]} km — near-uniform',
    f'D: L = {conv_distances[3]} km — uniform SU(2) (Haar)',
]

for idx, (d_km, title) in enumerate(zip(conv_distances, titles)):
    ax = fig2.add_subplot(gs2[idx // 2, idx % 2], projection='3d')
    S_vals = []
    for _ in range(N_REAL):
        out = apply_birefringence(E_in.copy(), d_km * 1e3, wavelength=WAVELENGTH,
                                  model='sectional')
        (S0, S1, S2, S3), _ = compute_stokes_parameters(out)
        S_vals.append([S1, S2, S3])
    S_vals = np.array(S_vals)
    mean_S = np.mean(S_vals, axis=0)
    ax.scatter(S_vals[:, 0], S_vals[:, 1], S_vals[:, 2],
               s=4, alpha=0.5, c='C3', marker='.')
    ax.scatter([mean_S[0]], [mean_S[1]], [mean_S[2]],
               s=80, c='k', marker='o', label=f'Mean = ({mean_S[0]:.3f}, {mean_S[1]:.3f}, {mean_S[2]:.3f})')
    ax.set_xlim(-1.1, 1.1); ax.set_ylim(-1.1, 1.1); ax.set_zlim(-1.1, 1.1)
    ax.set_xlabel('S1'); ax.set_ylabel('S2'); ax.set_zlabel('S3')
    ax.set_title(title, fontsize=9, fontweight='bold')
    ax.legend(fontsize=7, loc='upper left')
    ax.set_box_aspect([1, 1, 1])

fig2.suptitle('Convergence of Multi-Section Model to Uniform SU(2) on Poincare Sphere\n'
              f'({N_REAL} realisations per distance, input Ex=1, Ey=0)',
              fontsize=12, fontweight='bold', y=0.98)
fig2.savefig(os.path.join(OUT, f'val_birefringence_poincare--seed{SEED}.png'),
             dpi=200, bbox_inches='tight')
print(f"Saved: val_birefringence_poincare--seed{SEED}.png")
plt.close(fig2)

# Quantitative convergence check: mean Stokes vector → (0,0,0) as L increases
conv_results = []
for d_km in conv_distances:
    np.random.seed(SEED)
    S_all = []
    for _ in range(N_REAL):
        out = apply_birefringence(E_in.copy(), d_km * 1e3, wavelength=WAVELENGTH,
                                  model='sectional')
        (_, S1, S2, S3), _ = compute_stokes_parameters(out)
        S_all.append([S1, S2, S3])
    S_all = np.array(S_all)
    mean_norm = np.linalg.norm(np.mean(S_all, axis=0))
    conv_results.append(mean_norm)
    print(f"  L={d_km:6.3f} km: |mean(S)| = {mean_norm:.4f}  (0 -> uniform, 1 -> single state)")
assert conv_results[0] > 0.5, f"Sub-beat-length fibre should have high mean Stokes: {conv_results[0]}"
assert conv_results[-1] < 0.10, f"Long fibre should have near-zero mean Stokes: {conv_results[-1]}"
print(f"  [PASS] Poincare sphere convergence: |mean(S)| = {conv_results[0]:.3f} -> {conv_results[-1]:.3f}")
