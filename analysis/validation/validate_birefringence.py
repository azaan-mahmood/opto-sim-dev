import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os, sys, argparse

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from src.channel.fiber import cable

OUT = os.path.join(os.path.dirname(__file__), '..', 'val_biref')
os.makedirs(OUT, exist_ok=True)

parser = argparse.ArgumentParser()
parser.add_argument('--seed', type=int, default=42)
args = parser.parse_args()
SEED = args.seed
np.random.seed(SEED)

WAVELENGTH = 1550e-9
N_SAMPLES = 1024
DT = 1e-12
E_in = np.zeros((N_SAMPLES, 2), dtype=np.complex128)
E_in[:, 0] = 1.0 + 0j
E_in[:, 1] = 1.0 + 0j

# fiber.py coefficients
DN0 = 0.87e-5
TCOEFF = -5e-7
R_FIBER = 62.5e-6
BEND_COEFF = 0.135

def measure_delta_n(L_m, temp, bend_R):
    E_out = cable(L_m/1000, E_in.copy(), dt=DT, wavelength=WAVELENGTH,
                  dispersion=False, attenuation_factor=0.0,
                  temperature=temp, bend_radius=bend_R)
    # Relative phase between Ex and Ey (symmetric Jones formulation):
    #   dphi = dbeta·L = 2π·L·Δn/λ
    dphi = np.angle(np.mean(E_out[:, 0])) - np.angle(np.mean(E_out[:, 1]))
    dphi = np.arctan2(np.sin(dphi), np.cos(dphi))  # wrap to [-pi, pi]
    return dphi * WAVELENGTH / (2 * np.pi * L_m)

# --- 1. Base birefringence: sweep length at T=25, no bends ---
L_vals_m = np.array([1e-3, 2e-3, 5e-3, 1e-2, 2e-2])
dn_base = np.array([measure_delta_n(L, 25.0, None) for L in L_vals_m])
coeff_L = np.polyfit(L_vals_m, dn_base, 0)  # constant fit (should = DN0)
print(f"Base birefringence: measured = {np.mean(dn_base):.3e}, fiber.py = {DN0:.3e}")
print(f"  Error: {abs((np.mean(dn_base)-DN0)/DN0)*100:.4f}%")

# --- 2. Temperature: L = 1 cm, sweep 0-50 C ---
L_T = 1e-2
T_vals = np.array([0, 10, 25, 40, 50])
dn_T = np.array([measure_delta_n(L_T, T, None) for T in T_vals])
coeff_T = np.polyfit(T_vals - 25, dn_T, 1)
print(f"\nTemperature coeff: measured = {coeff_T[0]:.3e} /C, fiber.py = {TCOEFF:.3e} /C")
temp_err = abs((coeff_T[0] - TCOEFF) / TCOEFF) * 100
print(f"  Error: {temp_err:.2f}%")

# --- 3. Bend birefringence: L = 0.1 mm, sweep bend_radius ---
# Δn_bend = 0.135·(r_fiber/R)²  (Ulrich [7], Smith [8], Shibata [9]).
# Short fiber (0.1 mm) prevents phase wrapping for tight bends.
L_B = 1e-4
R_vals = np.array([0.002, 0.003, 0.005, 0.01, 0.02])  # 2 mm – 2 cm
dn_with_bend = np.array([measure_delta_n(L_B, 25.0, R) for R in R_vals])
dn_base_at_LB = measure_delta_n(L_B, 25.0, None)  # subtract base
dn_bend_measured = dn_with_bend - dn_base_at_LB

# Fit measured Δn_bend vs (r_f/R)² — slope should = BEND_COEFF
inv_R2_scale = (R_FIBER / R_vals) ** 2
slope, intercept = np.polyfit(inv_R2_scale, dn_bend_measured, 1)
dn_bend_expected = BEND_COEFF * inv_R2_scale
print(f"\nBend coeff (slope of dn vs (r_f/R)^2): measured = {slope:.4e}, fiber.py = {BEND_COEFF:.4e}")
bend_err = abs((slope - BEND_COEFF) / BEND_COEFF) * 100
print(f"  Error: {bend_err:.4f}%")
print(f"  Intercept: {intercept:.3e} (should be ~0)")

# --- Plot ---
fig, axes = plt.subplots(1, 3, figsize=(12, 4))

# Base dn
axes[0].axhline(DN0, c='gray', ls='--', label=f'fiber.py: {DN0:.3e}')
axes[0].plot(L_vals_m*1e3, dn_base, 'o-', c='C0', label='cable() measured')
axes[0].set(xlabel='Length (mm)', ylabel='Delta n',
            title=f'Base birefringence (T=25, no bends)')
axes[0].legend(); axes[0].grid(True, alpha=0.3)

# Temperature
dn_T_expected = DN0 + TCOEFF * (T_vals - 25)
axes[1].plot(T_vals, dn_T_expected, 's--', c='C1', ms=4, label=f'fiber.py: {TCOEFF:.1e}/C')
axes[1].plot(T_vals, dn_T, 'o-', c='C0', ms=5, label=f'cable() fit: {coeff_T[0]:.2e}/C')
axes[1].set(xlabel='Temperature (C)', ylabel='Delta n',
            title=f'Temperature (L = {L_T*1e3:.0f} mm)')
axes[1].legend(); axes[1].grid(True, alpha=0.3)

# Bends
axes[2].plot(R_vals*1e3, dn_bend_expected, 's--', c='C1', ms=4,
             label=f'$0.135\\cdot(r_f/R)^2$')
axes[2].plot(R_vals*1e3, dn_bend_measured, 'o-', c='C0', ms=5,
             label=f'cable() fit; slope={slope:.3e}')
axes[2].set(xlabel='Bend radius (mm)', ylabel='Delta n bend contribution',
            title=f'Bend birefringence (L = {L_B*1e3:.0f} mm)')
axes[2].legend(); axes[2].grid(True, alpha=0.3)

fig.tight_layout()
fname = f'val_biref--seed{SEED}.png'
fig.savefig(os.path.join(OUT, fname), dpi=150)
print(f"\nSaved: {fname}")

csv_name = f'val_biref--seed{SEED}.csv'
np.savetxt(os.path.join(OUT, csv_name),
           np.column_stack([L_vals_m*1e3, dn_base,
                           T_vals, dn_T,
                           R_vals*1e3, dn_bend_measured, dn_bend_expected]),
           delimiter=',',
           header='L_mm,dn_base,T_C,dn_T,R_mm,dn_bend_meas,dn_bend_exp', comments='')
print(f"Saved: {csv_name}")
plt.close(fig)
