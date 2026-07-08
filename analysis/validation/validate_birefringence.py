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
BEND_EFF = 2.4e-4

def measure_delta_n(L_m, temp, n_bends):
    E_out = cable(L_m/1000, E_in.copy(), dt=DT, wavelength=WAVELENGTH,
                  dispersion=False, attenuation_factor=0.0,
                  temperature=temp, num_bends=n_bends)
    # Relative phase between Ex and Ey (symmetric Jones formulation):
    #   dphi = dbeta·L = 2π·L·Δn/λ
    dphi = np.angle(np.mean(E_out[:, 0])) - np.angle(np.mean(E_out[:, 1]))
    dphi = np.arctan2(np.sin(dphi), np.cos(dphi))  # wrap to [-pi, pi]
    return dphi * WAVELENGTH / (2 * np.pi * L_m)

# --- 1. Base birefringence: sweep length at T=25, no bends ---
L_vals_m = np.array([1e-3, 2e-3, 5e-3, 1e-2, 2e-2])
dn_base = np.array([measure_delta_n(L, 25.0, 0) for L in L_vals_m])
coeff_L = np.polyfit(L_vals_m, dn_base, 0)  # constant fit (should = DN0)
print(f"Base birefringence: measured = {np.mean(dn_base):.3e}, fiber.py = {DN0:.3e}")
print(f"  Error: {abs((np.mean(dn_base)-DN0)/DN0)*100:.4f}%")

# --- 2. Temperature: L = 1 cm, sweep 0-50 C ---
L_T = 1e-2
T_vals = np.array([0, 10, 25, 40, 50])
dn_T = np.array([measure_delta_n(L_T, T, 0) for T in T_vals])
coeff_T = np.polyfit(T_vals - 25, dn_T, 1)
print(f"\nTemperature coeff: measured = {coeff_T[0]:.3e} /C, fiber.py = {TCOEFF:.3e} /C")
temp_err = abs((coeff_T[0] - TCOEFF) / TCOEFF) * 100
print(f"  Error: {temp_err:.2f}%")

# --- 3. Bends: L = 0.1 mm, sweep 0-5 bends ---
L_B = 1e-4
bend_vals = np.array([0, 1, 2, 3, 5])
dn_B = np.array([measure_delta_n(L_B, 25.0, nb) for nb in bend_vals])
coeff_B = np.polyfit(bend_vals, dn_B, 1)
dn_B_expected = DN0 + BEND_EFF * bend_vals
print(f"\nBend coeff: measured = {coeff_B[0]:.3e} /bend, fiber.py = {BEND_EFF:.3e} /bend")
bend_err = abs((coeff_B[0] - BEND_EFF) / BEND_EFF) * 100
print(f"  Error: {bend_err:.2f}%")

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
dn_B_fit = coeff_B[0] * bend_vals + coeff_B[1]
axes[2].plot(bend_vals, dn_B_expected, 's--', c='C1', ms=4, label=f'fiber.py: {BEND_EFF:.1e}/bend')
axes[2].plot(bend_vals, dn_B, 'o-', c='C0', ms=5, label=f'cable() fit: {coeff_B[0]:.2e}/bend')
axes[2].set(xlabel='Num bends', ylabel='Delta n',
            title=f'Bends (L = {L_B*1e3:.3f} mm)')
axes[2].legend(); axes[2].grid(True, alpha=0.3)

fig.tight_layout()
fname = f'val_biref--seed{SEED}.png'
fig.savefig(os.path.join(OUT, fname), dpi=150)
print(f"\nSaved: {fname}")

csv_name = f'val_biref--seed{SEED}.csv'
np.savetxt(os.path.join(OUT, csv_name),
           np.column_stack([L_vals_m*1e3, dn_base,
                           T_vals, dn_T,
                           bend_vals, dn_B]),
           delimiter=',',
           header='L_mm,dn_base,T_C,dn_T,num_bends,dn_bends', comments='')
print(f"Saved: {csv_name}")
plt.close(fig)
