import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os, sys, argparse

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from src.channel.fiber import cable

OUT = os.path.join(os.path.dirname(__file__), '..', 'val_cd')
os.makedirs(OUT, exist_ok=True)

parser = argparse.ArgumentParser()
parser.add_argument('--seed', type=int, default=42)
args = parser.parse_args()
SEED = args.seed
np.random.seed(SEED)

from src.channel.fiber import D_TOTAL

LIGHTSPEED = 299792458.0
WAVELENGTH = 1550e-9
D_total = D_TOTAL             # ps/(nm·km), imported dynamically from fiber.py
D_SI = D_total * 1e-6
beta2 = -D_SI * WAVELENGTH**2 / (2 * np.pi * LIGHTSPEED)

T0 = 30e-12
LD = T0**2 / abs(beta2)
z_over_LD = np.array([0.0, 0.5, 1.0, 2.0])
DT = 1e-12
N_SAMPLES = 2**12
T = np.arange(N_SAMPLES) * DT
center = N_SAMPLES // 2
t_arr = T - T[center]

print(f"beta2 = {beta2:.3e} s^2/m  (D = {D_total} ps/(nm km))")
print(f"T0 = {T0*1e12:.1f} ps,  LD = {LD/1e3:.2f} km\n")

# --- Generate Gaussian pulse and propagate through cable() ---
E_in = np.zeros((N_SAMPLES, 2), dtype=np.complex128)
E_in[:, 0] = np.exp(-0.5 * (t_arr / T0)**2)

distances_km = z_over_LD * LD / 1e3
sigma0 = T0 / np.sqrt(2)
widths = []

for L_km in distances_km:
    if L_km == 0:
        E_out = E_in.copy()
    else:
        E_out = cable(L_km, E_in.copy(), dt=DT, wavelength=WAVELENGTH,
                      dispersion=True, pm_dispersion=0.0,
                      attenuation_factor=0.0, temperature=25.0, num_bends=0)
    I = np.abs(E_out[:, 0])**2
    total = I.sum()
    t_mean = np.sum(T * I) / total
    t_var = np.sum((T - t_mean)**2 * I) / total
    widths.append(np.sqrt(t_var.real))

widths = np.array(widths)
analytic = sigma0 * np.sqrt(1 + z_over_LD**2)
ratio_m = widths / sigma0
ratio_a = analytic / sigma0
errors = np.abs(ratio_m - ratio_a) / ratio_a * 100

print(f"{'z/LD':>6s}  {'z(km)':>8s}  {'sigma(ps)':>11s}  {'analytic':>11s}  {'error%':>7s}")
for i, zf in enumerate(z_over_LD):
    print(f"{zf:6.1f}  {distances_km[i]:8.2f}  {widths[i]*1e12:11.6f}  {analytic[i]*1e12:11.6f}  {errors[i]:7.4f}")

# --- Plot ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))
colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(z_over_LD)))
for i, L_km in enumerate(distances_km):
    E_out = cable(L_km, E_in.copy(), dt=DT, wavelength=WAVELENGTH,
                  dispersion=True, pm_dispersion=0.0,
                  attenuation_factor=0.0, temperature=25.0, num_bends=0) if L_km > 0 else E_in
    I = np.abs(E_out[:, 0])**2
    ax1.plot(T*1e12, I/I.max(), color=colors[i],
             label=f'z/L_D = {z_over_LD[i]:.1f}')
ax1.set(xlabel='Time (ps)', ylabel='Normalized intensity',
        title='Gaussian pulse broadening (cable() output)')
ax1.legend(fontsize=8); ax1.grid(True, alpha=0.3)

z_plot = np.linspace(0, 2.5, 200)
ax2.plot(z_plot, np.sqrt(1+z_plot**2), 'k-', lw=1.5,
         label=r'Analytic $\sqrt{1+(z/L_D)^2}$')
ax2.plot(z_over_LD, ratio_m, 'o', c='C3', ms=6, label='cable() measured')
for zf, err, rm in zip(z_over_LD, errors, ratio_m):
    ax2.annotate(f'{err:.4f}%', (zf, rm), fontsize=7, va='bottom', ha='left', c='C3')
ax2.set(xlabel=r'$z/L_D$', ylabel=r'$\sigma(z)/\sigma_0$',
        title=f'CD validation (D={D_total} ps/(nm·km), T0={T0*1e12:.0f} ps)',
        xlim=(0, 2.5))
ax2.legend(fontsize=8); ax2.grid(True, alpha=0.3)

fig.tight_layout()
fname = f'val_cd--seed{SEED}.png'
fig.savefig(os.path.join(OUT, fname), dpi=150)
print(f"\nSaved: {fname}")

csv_name = f'val_cd--seed{SEED}.csv'
np.savetxt(os.path.join(OUT, csv_name),
           np.column_stack([z_over_LD, distances_km, widths*1e12, analytic*1e12, errors]),
           delimiter=',', header='z_over_LD,distance_km,sigma_rms_ps,analytic_sigma_ps,error_pct', comments='')
print(f"Saved: {csv_name}")
plt.close(fig)
