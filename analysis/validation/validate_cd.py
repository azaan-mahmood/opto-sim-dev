import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os, sys, argparse

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from src.channel.fiber import cable

OUT = os.path.join(os.path.dirname(__file__), '..', 'val_cd')

parser = argparse.ArgumentParser()
parser.add_argument('--seed', type=int, default=42)
args = parser.parse_args()
SEED = args.seed
rng = np.random.default_rng(SEED)
np.random.seed(SEED)

LIGHTSPEED = 299792458.0
WAVELENGTH = 1550e-9

# --- fibre dispersion parameters (must match fiber.py) ---
D_material = 17.0       # ps/(nm·km)
D_waveguide = -3.0      # ps/(nm·km)
D_total = D_material + D_waveguide    # 14 ps/(nm·km)
D_SI = D_total * 1e-6                # s/m²
beta2 = -D_SI * WAVELENGTH**2 / (2 * np.pi * LIGHTSPEED)   # s²/m

# --- Gaussian pulse parameters ---
T0 = 30e-12                        # 1/e half-width of intensity (s)
LD = T0**2 / abs(beta2)            # dispersion length (m)
z_over_LD = np.array([0.0, 0.5, 1.0, 2.0])
distances = z_over_LD * LD         # m
DT = 1e-12
N_SAMPLES = 2**12                  # plenty of margin
T = np.arange(N_SAMPLES) * DT
center_idx = N_SAMPLES // 2

print(f"beta2 = {beta2:.3e} s^2/m  (D = {D_total} ps/(nm km))")
print(f"T0 = {T0*1e12:.1f} ps,  LD = {LD/1e3:.2f} km")
print()

# --- generate input Gaussian pulse (complex envelope) ---
t_arr = T - T[center_idx]
pulse_envelope = np.exp(-0.5 * (t_arr / T0)**2)        # amplitude: exp(-t²/(2T0²))
E_in = np.zeros((N_SAMPLES, 2), dtype=np.complex128)
E_in[:, 0] = pulse_envelope + 0j                        # Ex only

# --- propagate at each distance ---
widths = []
for z_frac in z_over_LD:
    L_m = z_frac * LD
    if L_m == 0:
        E_out = E_in.copy()
    else:
        E_out = cable(
            fiber_length=L_m / 1e3,
            E=E_in,
            dt=DT,
            wavelength=WAVELENGTH,
            dispersion=True,
            pm_dispersion=0.0,          # no PMD
            attenuation_factor=0.0,      # no loss
            temperature=25.0,
            num_bends=0,
        )

    intensity = np.abs(E_out[:, 0])**2
    # Gaussian fit: intensity = A * exp(-(t - t0)^2 / (2 * sigma^2))
    # Use the moments of the intensity for a robust estimate:
    total = intensity.sum()
    if total == 0:
        widths.append(0.0)
        continue
    t_mean = np.sum(T * intensity) / total
    t_var = np.sum((T - t_mean)**2 * intensity) / total
    sigma_rms = np.sqrt(t_var.real)        # RMS width
    # For a Gaussian pulse, sigma_rms = T0 / sqrt(2), and
    # the analytic broadened RMS width = sqrt(T0²/2 + z²*|β₂|²/T0²)
    # which simplifies to sigma_rms = (T0/sqrt(2)) * sqrt(1 + (z/LD)²)
    widths.append(sigma_rms)

widths = np.array(widths)

# --- analytic prediction ---
sigma0 = T0 / np.sqrt(2)                   # initial RMS width
analytic = sigma0 * np.sqrt(1 + z_over_LD**2)
ratio_measured = widths / sigma0
ratio_analytic = analytic / sigma0
errors = np.abs(ratio_measured - ratio_analytic) / ratio_analytic * 100

print(f"{'z/LD':>6s}  {'z (km)':>8s}  {'sigma_rms(ps)':>15s}  {'analytic(ps)':>15s}  {'error(%)':>8s}")
for i, zf in enumerate(z_over_LD):
    print(f"{zf:6.1f}  {distances[i]/1e3:8.2f}  {widths[i]*1e12:15.6f}  {analytic[i]*1e12:15.6f}  {errors[i]:8.4f}")

# --- plot ---
fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))

# Left: pulse shapes at each distance
for i, zf in enumerate(z_over_LD):
    if distances[i] == 0:
        L_km = 0
        intensity = np.abs(E_in[:, 0])**2
        label = 'z = 0'
    else:
        L_km = distances[i] / 1e3
        E_out = cable(
            fiber_length=L_km,
            E=E_in.copy(),
            dt=DT,
            wavelength=WAVELENGTH,
            dispersion=True,
            pm_dispersion=0.0,
            attenuation_factor=0.0,
            temperature=25.0,
            num_bends=0,
        )
        intensity = np.abs(E_out[:, 0])**2
        label = f'z/L_D = {zf:.1f}'
    axes[0].plot(T * 1e12, intensity / intensity.max(), label=label)
axes[0].set_xlabel('Time (ps)')
axes[0].set_ylabel('Normalized intensity')
axes[0].set_title('Gaussian pulse broadening (CD only)')
axes[0].legend(fontsize=8)
axes[0].grid(True, alpha=0.3)

# Right: RMS width ratio vs z/L_D
z_plot = np.linspace(0, 2.5, 200)
analytic_plot = sigma0 * np.sqrt(1 + z_plot**2)
axes[1].plot(z_plot, analytic_plot / sigma0, 'k-', linewidth=1.5, label='Analytic: $\\sqrt{1 + (z/L_D)^2}$')
axes[1].plot(z_over_LD, ratio_measured, 'o', color='C3', markersize=6, label='Simulation')
for i, zf in enumerate(z_over_LD):
    axes[1].annotate(f'  {errors[i]:.3f}%', (zf, ratio_measured[i]),
                     fontsize=7, va='bottom', color='C3')
axes[1].set_xlabel('$z / L_D$')
axes[1].set_ylabel('$\\sigma(z) / \\sigma_0$')
axes[1].set_title(f'CD validation — D = {D_total} ps/(nm·km), $T_0$ = {T0*1e12:.0f} ps')
axes[1].legend(fontsize=8)
axes[1].grid(True, alpha=0.3)
axes[1].set_xlim(0, 2.5)

fig.tight_layout()
fname = f'val_cd--seed{SEED}.png'
fig.savefig(os.path.join(OUT, fname), dpi=150, bbox_inches='tight')
print(f"\nSaved: {fname}")

# --- CSV ---
csv_name = f'val_cd--seed{SEED}.csv'
header = 'z_over_LD,distance_km,sigma_rms_ps,analytic_sigma_ps,error_pct'
np.savetxt(os.path.join(OUT, csv_name),
           np.column_stack([z_over_LD, distances/1e3, widths*1e12, analytic*1e12, errors]),
           delimiter=',', header=header, comments='')
print(f"Saved: {csv_name}")
plt.close(fig)
