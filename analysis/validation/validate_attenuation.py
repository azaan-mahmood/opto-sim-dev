import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os, sys, argparse

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from src.channel.fiber import apply_attenuation

OUT = os.path.join(os.path.dirname(__file__), '..', 'val_attenuation')
os.makedirs(OUT, exist_ok=True)

parser = argparse.ArgumentParser()
parser.add_argument('--seed', type=int, default=42)
args = parser.parse_args()
SEED = args.seed
np.random.seed(SEED)

N_SAMPLES = 10000
E_in = np.zeros((N_SAMPLES, 2), dtype=complex)
E_in[:, 0] = np.sqrt(0.5) * (np.random.randn(N_SAMPLES) + 1j * np.random.randn(N_SAMPLES))
P_in = np.mean(np.abs(E_in)**2)

distances_km = np.linspace(0, 200, 41)
attenuation_dB_per_km = 0.182
P_out_theory = P_in * 10 ** (-attenuation_dB_per_km * distances_km / 10.0)
P_out_meas = np.array([
    np.mean(np.abs(apply_attenuation(E_in.copy(), L_km, attenuation_dB_per_km))**2)
    for L_km in distances_km
])

errors = np.abs(P_out_meas - P_out_theory) / P_out_theory * 100

print(f"Attenuation validation via apply_attenuation — alpha = {attenuation_dB_per_km} dB/km")
print(f"{'Distance(km)':>12s}  {'P_meas':>12s}  {'P_theory':>12s}  {'error%':>8s}")
for i, L in enumerate(distances_km):
    print(f"{L:12.1f}  {P_out_meas[i]:12.8f}  {P_out_theory[i]:12.8f}  {errors[i]:8.6f}")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))
ax1.semilogy(distances_km, P_out_theory, 'k-', lw=1.5, label=f'Theory: $10^{{-\\alpha L/10}}$')
ax1.semilogy(distances_km, P_out_meas, 'o', ms=4, c='C3', label='apply_attenuation')
ax1.set(xlabel='Distance (km)', ylabel='Output Power (a.u.)',
        title='Attenuation validation (apply_attenuation)')
ax1.legend(fontsize=8); ax1.grid(True, alpha=0.3)

ax2.plot(distances_km, errors, 'o-', ms=3, c='C3')
ax2.set(xlabel='Distance (km)', ylabel='|Error| (%)', title='Percent error')
ax2.grid(True, alpha=0.3)

fig.tight_layout()
fname = f'val_attenuation--seed{SEED}.png'
fig.savefig(os.path.join(OUT, fname), dpi=150)
print(f"\nSaved: {fname}")

csv_name = f'val_attenuation--seed{SEED}.csv'
np.savetxt(os.path.join(OUT, csv_name),
           np.column_stack([distances_km, P_out_meas, P_out_theory, errors]),
           delimiter=',', header='distance_km,power_meas,power_theory,error_pct', comments='')
print(f"Saved: {csv_name}")
plt.close(fig)
