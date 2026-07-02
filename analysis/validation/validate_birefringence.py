import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os, argparse

OUT = os.path.join(os.path.dirname(__file__), '..', 'val_biref')

parser = argparse.ArgumentParser()
parser.add_argument('--seed', type=int, default=42)
args = parser.parse_args()
SEED = args.seed
np.random.seed(SEED)

WAVELENGTH = 1550e-9   # m
FIBER_RADIUS = 62.5e-6 # m, SMF-28 cladding radius

# --- Strain-optic coefficients for silica (p_ij: Pockels coefficients) ---
# Ref: Xu & Stroud, "Stress-optic coefficient of silica fiber"
N_SILICA = 1.444
P11 = 0.121
P12 = 0.270
C_BEND = 0.5 * N_SILICA**3 * np.sqrt((P11 - P12)**2 + (P11 - 3*P12)**2 / 4.0)
# For a linear birefringence model, the simpler form is:
# Δn_bend = -0.5 · n³ · (p₁₁ - p₁₂) · (r/R)²
# The magnitude of bend-induced birefringence:
BEND_COEFF = 0.5 * N_SILICA**3 * abs(P11 - P12) * FIBER_RADIUS**2

print(f"Bend birefringence coefficient: dn_bend = {BEND_COEFF:.3e} / R^2  (R in m)")
print(f"  At R = 1 cm: dn_bend = {BEND_COEFF/0.01**2:.3e}")
print(f"  At R = 5 mm: dn_bend = {BEND_COEFF/0.005**2:.3e}")
print()

# --- Bend radii (mm) ---
R_mm = np.logspace(np.log10(5), np.log10(100), 200)
R_m = R_mm * 1e-3

# --- Intrinsic birefringence values for the curves ---
# Yuan [4] Fig 1 uses Δn₀ = 0.5, 1.0, 2.0 × 10⁻⁵
delta_n0_vals = [0.5e-5, 1.0e-5, 2.0e-5]
labels = [r'$\Delta n_0 = 0.5\times 10^{-5}$',
          r'$\Delta n_0 = 1.0\times 10^{-5}$',
          r'$\Delta n_0 = 2.0\times 10^{-5}$']

# --- Current model from fiber.py for reference ---
# fiber.py uses: Δn = Δn₀ + bend_effect_factor × num_bends
# where bend_effect_factor = 2.4e-4 and num_bends is a count (not R-dependent)
# For comparison, we show the "bend effect" as a horizontal line at the
# beat length corresponding to bend_effect_factor / revolution
# (assuming 1 bend at R gives some Δn contribution)
print("Current fiber.py model: bend_effect_factor = 2.4e-4 per bend")
print("  This is a fixed scalar, NOT a function of bend radius R.")
print("  To reproduce Yuan Fig 1, the bend contribution must scale as 1/R^2.")
print()

# --- Reproduce Yuan Fig 1: L_B vs R for various Δn₀ ---
fig, ax = plt.subplots(figsize=(8, 5))

for dn0, lab in zip(delta_n0_vals, labels):
    # Total birefringence = intrinsic + bend-induced
    dn_bend = BEND_COEFF / R_m**2
    dn_total = dn0 + dn_bend
    L_B = WAVELENGTH / dn_total      # beat length (m)
    ax.plot(R_mm, L_B * 100, linewidth=1.5, label=lab)

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Bend radius R (mm)')
ax.set_ylabel('Beat length $L_B$ (cm)')
ax.set_title('Birefringence validation — Yuan [4] Fig 1')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3, which='both')
ax.set_xlim(5, 100)

# Annotate key scaling regimes
ax.annotate('$L_B \\propto R$ (bend-dominated)', xy=(8, 0.4), fontsize=8,
            rotation=35, color='gray', alpha=0.7)
ax.annotate('plateau (intrinsic $\\Delta n_0$ dominates)', xy=(30, 5), fontsize=8,
            color='gray', alpha=0.7)

fig.tight_layout()
fname = f'val_biref--seed{SEED}.png'
fig.savefig(os.path.join(OUT, fname), dpi=150, bbox_inches='tight')
print(f"Saved: {fname}")

# --- CSV ---
csv_name = f'val_biref--seed{SEED}.csv'
header = 'R_mm,dn0_0.5e-5_LB_cm,dn0_1.0e-5_LB_cm,dn0_2.0e-5_LB_cm'
lb_vals = []
for dn0 in delta_n0_vals:
    dn_total = dn0 + BEND_COEFF / R_m**2
    lb_vals.append(WAVELENGTH / dn_total * 100)
np.savetxt(os.path.join(OUT, csv_name),
           np.column_stack([R_mm] + lb_vals),
           delimiter=',', header=header, comments='')
print(f"Saved: {csv_name}")
plt.close(fig)
