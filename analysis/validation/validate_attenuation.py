import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os, sys, argparse

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from src.channel.fiber import cable

OUT = os.path.join(os.path.dirname(__file__), '..', 'val_att')

parser = argparse.ArgumentParser()
parser.add_argument('--seed', type=int, default=42)
args = parser.parse_args()
SEED = args.seed
np.random.seed(SEED)

WAVELENGTH = 1550e-9
ATT_FACTOR = 0.182       # dB/km, SMF-28 at 1550 nm
DISTANCES = np.linspace(0, 200, 41)

# --- generate a CW field ---
N_SAMPLES = 1024
DT = 1e-12
E_in = np.ones((N_SAMPLES, 2), dtype=np.complex128) * np.sqrt(0.5)  # 0.5 W total

# --- sweep ---
ratios = []
for L in DISTANCES:
    if L == 0:
        power_out = np.mean(np.abs(E_in)**2)
    else:
        E_out = cable(
            fiber_length=L,
            E=E_in.copy(),
            dt=DT,
            wavelength=WAVELENGTH,
            dispersion=False,
            attenuation_factor=ATT_FACTOR,
            temperature=25.0,
            num_bends=0,
        )
        power_out = np.mean(np.abs(E_out)**2)
    ratios.append(power_out / np.mean(np.abs(E_in)**2))

ratios = np.array(ratios)

# --- analytic ---
analytic = 10.0 ** (-ATT_FACTOR * DISTANCES / 10.0)
errors = np.abs(ratios - analytic) / analytic * 100

# --- dB scale ---
ratios_db = 10 * np.log10(ratios)
analytic_db = -ATT_FACTOR * DISTANCES

print(f"Attenuation validation: {ATT_FACTOR} dB/km at {WAVELENGTH*1e9:.0f} nm")
print(f"{'L (km)':>8s}  {'P_out/P_in':>12s}  {'analytic':>12s}  {'error(%)':>8s}  {'dB':>8s}")
for i, L in enumerate(DISTANCES[::5]):
    print(f"{L:8.1f}  {ratios[i]:12.6f}  {analytic[i]:12.6f}  {errors[i]:8.4f}  {ratios_db[i]:8.4f}")

# --- plot ---
fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))

# Linear scale
axes[0].plot(DISTANCES, ratios, 'o-', color='C0', markersize=3, label='Simulation')
axes[0].plot(DISTANCES, analytic, '--', color='C1', linewidth=1.5, label='$10^{-\\alpha L/10}$')
axes[0].set_xlabel('Distance (km)')
axes[0].set_ylabel('$P_{out} / P_{in}$')
axes[0].set_title(f'Attenuation — SMF-28 at {WAVELENGTH*1e9:.0f} nm ({ATT_FACTOR} dB/km)')
axes[0].legend(fontsize=9)
axes[0].grid(True, alpha=0.3)

# dB scale
axes[1].plot(DISTANCES, ratios_db, 'o-', color='C0', markersize=3, label='Simulation')
axes[1].plot(DISTANCES, analytic_db, '--', color='C1', linewidth=1.5, label=f'{-ATT_FACTOR:.3f} dB/km × L')
axes[1].set_xlabel('Distance (km)')
axes[1].set_ylabel('Attenuation (dB)')
axes[1].set_title('Log scale — linear fit')
axes[1].legend(fontsize=9)
axes[1].grid(True, alpha=0.3)

fig.tight_layout()
fname = f'val_att--seed{SEED}.png'
fig.savefig(os.path.join(OUT, fname), dpi=150, bbox_inches='tight')
print(f"\nSaved: {fname}")

csv_name = f'val_att--seed{SEED}.csv'
header = 'distance_km,power_ratio,analytic_ratio,error_pct'
np.savetxt(os.path.join(OUT, csv_name),
           np.column_stack([DISTANCES, ratios, analytic, errors]),
           delimiter=',', header=header, comments='')
print(f"Saved: {csv_name}")
plt.close(fig)
