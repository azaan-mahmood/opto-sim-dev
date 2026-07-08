import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os, sys, argparse

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from src.channel.fiber import cable

OUT = os.path.join(os.path.dirname(__file__), '..', 'val_att')
os.makedirs(OUT, exist_ok=True)

parser = argparse.ArgumentParser()
parser.add_argument('--seed', type=int, default=42)
args = parser.parse_args()
SEED = args.seed
np.random.seed(SEED)

WAVELENGTH = 1550e-9
ATT_FACTOR = 0.182
DISTANCES = np.linspace(0, 200, 41)
N_SAMPLES = 1024
DT = 1e-12

E_in = (np.ones((N_SAMPLES, 2), dtype=np.complex128) * np.sqrt(0.5))
power_in = np.mean(np.abs(E_in)**2)

ratios = []
for L in DISTANCES:
    if L == 0:
        ratios.append(1.0)
    else:
        E_out = cable(L, E_in.copy(), dt=DT, wavelength=WAVELENGTH,
                      dispersion=False, attenuation_factor=ATT_FACTOR,
                      temperature=25.0, bend_radius=None)
        ratios.append(np.mean(np.abs(E_out)**2) / power_in)

ratios = np.array(ratios)
analytic = 10.0 ** (-ATT_FACTOR * DISTANCES / 10.0)
errors = np.abs(ratios - analytic) / analytic * 100

print(f"Attenuation via cable() — {ATT_FACTOR} dB/km SMF-28 at 1550 nm")
print(f"{'L(km)':>7s}  {'P_out/P_in':>11s}  {'analytic':>11s}  {'error%':>7s}")
for i in range(0, len(DISTANCES), 5):
    print(f"{DISTANCES[i]:7.1f}  {ratios[i]:11.6f}  {analytic[i]:11.6f}  {errors[i]:7.4f}")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))
ax1.plot(DISTANCES, ratios, 'o-', c='C0', ms=3, label='cable()')
ax1.plot(DISTANCES, analytic, '--', c='C1', lw=1.5, label=r'$10^{-\alpha L/10}$')
ax1.set(xlabel='Distance (km)', ylabel=r'$P_{out}/P_{in}$',
        title=f'Attenuation — {ATT_FACTOR} dB/km @ {WAVELENGTH*1e9:.0f} nm')
ax1.legend(); ax1.grid(True, alpha=0.3)

ax2.plot(DISTANCES, 10*np.log10(ratios), 'o-', c='C0', ms=3, label='cable()')
ax2.plot(DISTANCES, -ATT_FACTOR*DISTANCES, '--', c='C1', lw=1.5, label=f'{-ATT_FACTOR:.3f} dB/km')
ax2.set(xlabel='Distance (km)', ylabel='Attenuation (dB)',
        title='Log scale — linear fit')
ax2.legend(); ax2.grid(True, alpha=0.3)

fig.tight_layout()
fname = f'val_att--seed{SEED}.png'
fig.savefig(os.path.join(OUT, fname), dpi=150)
print(f"\nSaved: {fname}")

csv_name = f'val_att--seed{SEED}.csv'
np.savetxt(os.path.join(OUT, csv_name),
           np.column_stack([DISTANCES, ratios, analytic, errors]),
           delimiter=',', header='distance_km,power_ratio,analytic_ratio,error_pct', comments='')
print(f"Saved: {csv_name}")
plt.close(fig)
