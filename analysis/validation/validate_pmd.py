import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import os, argparse

OUT = os.path.join(os.path.dirname(__file__), '..', 'val_pmd')

parser = argparse.ArgumentParser()
parser.add_argument('--seed', type=int, default=42)
args = parser.parse_args()
SEED = args.seed
np.random.seed(SEED)

N_REALIZATIONS = 20000
PM_DISPERSION = 0.1e-12        # s / sqrt(m)
L_KM = 50
L_M = L_KM * 1e3

# --- replicate fiber.py DGD generation ---
pmd_sd = PM_DISPERSION * np.sqrt(L_M)   # seconds (RMS DGD)
scale_rayleigh = pmd_sd * np.sqrt(2.0 / np.pi)

dgd_samples = np.random.rayleigh(scale_rayleigh, size=N_REALIZATIONS)

# --- statistics ---
mean_dgd = np.mean(dgd_samples)
rms_dgd = np.sqrt(np.mean(dgd_samples**2))
target_mean = pmd_sd                      # the code aims for mean = pmd_sd
target_rms = pmd_sd                       # the code claims RMS = pmd_sd
# Actual Rayleigh(scale=σ) has: mean = σ·√(π/2), variance = (2-π/2)·σ²
# With scale = pmd_sd·√(2/π): mean = pmd_sd·√(2/π)·√(π/2) = pmd_sd ✓
# RMS² = variance + mean² = (2-π/2)·σ² + σ²·π/2 = 2·σ²
# RMS = σ·√2 = pmd_sd·√(2/π)·√2 = pmd_sd·2/√π ≈ 1.128·pmd_sd
rayleigh_mean = scale_rayleigh * np.sqrt(np.pi / 2)
rayleigh_rms = scale_rayleigh * np.sqrt(2)

# --- Maxwellian comparison ---
# Maxwell(scale=a) is the distribution of sqrt(X1²+X2²+X3²) where Xi ~ N(0,a²)
# mean = 2a·√(2/π), RMS = a·√3
# For PMD, we want mean = pmd_sd (RMS DGD coefficient)
# So a = pmd_sd / √3, which gives mean = 2·pmd_sd/√3·√(2/π) ≈ 0.921·pmd_sd
# and RMS = pmd_sd
maxwell_a = pmd_sd / np.sqrt(3)
maxwell_mean = 2 * maxwell_a * np.sqrt(2 / np.pi)
maxwell_rms = maxwell_a * np.sqrt(3)

print(f"PMD DGD distribution validation ({N_REALIZATIONS} realizations)")
print(f"  PMD coeff = {PM_DISPERSION*1e12:.2f} ps/sqrt(m)")
print(f"  L = {L_KM} km")
print(f"  Target RMS DGD (pmd_sd) = {pmd_sd*1e12:.4f} ps")
print()
print(f"Current model (Rayleigh):")
print(f"  scale = {scale_rayleigh*1e12:.4f} ps")
print(f"  Measured mean DGD = {mean_dgd*1e12:.4f} ps  (expected = {rayleigh_mean*1e12:.4f})")
print(f"  Measured RMS DGD  = {rms_dgd*1e12:.4f} ps  (expected = {rayleigh_rms*1e12:.4f})")
print(f"  RMS/mean = {rms_dgd/mean_dgd:.4f}  (Rayleigh theory: {np.sqrt(4/np.pi):.4f})")
print()
print(f"Maxwellian (PMD theory):")
print(f"  scale a = {maxwell_a*1e12:.4f} ps  (a = RMS/sqrt(3))")
print(f"  Expected mean = {maxwell_mean*1e12:.4f} ps")
print(f"  Expected RMS  = {maxwell_rms*1e12:.4f} ps")
print(f"  RMS/mean = {maxwell_rms/maxwell_mean:.4f}  (Maxwell theory: {np.sqrt(3/(8/np.pi)):.4f})")

# --- fit both distributions to measured data ---
rayleigh_fit_scale = np.sqrt(np.mean(dgd_samples**2) / 2)   # MLE for Rayleigh
rayleigh_fit = stats.rayleigh(scale=rayleigh_fit_scale)

# Maxwell fit (method of moments using RMS)
maxwell_fit_a = rms_dgd / np.sqrt(3)
maxwell_fit = stats.maxwell(scale=maxwell_fit_a)

# KS tests
ks_rayleigh = stats.kstest(dgd_samples, 'rayleigh', args=(rayleigh_fit_scale,))
ks_maxwell = stats.kstest(dgd_samples, 'maxwell', args=(maxwell_fit_a,))
print()
print(f"KS test vs fitted Rayleigh:  statistic={ks_rayleigh.statistic:.6f}, p-value={ks_rayleigh.pvalue:.4f}")
print(f"KS test vs fitted Maxwell:   statistic={ks_maxwell.statistic:.6f}, p-value={ks_maxwell.pvalue:.4f}")

# --- plot ---
fig, ax = plt.subplots(figsize=(8, 5))
bins = np.linspace(0, max(dgd_samples)*1.05, 80)
ax.hist(dgd_samples * 1e12, bins=bins, density=True, alpha=0.6, color='C0',
        label=f'Simulated ({N_REALIZATIONS} draws)')

dgd_grid = np.linspace(0, max(dgd_samples)*1.05, 500)
pdf_rayleigh = stats.rayleigh.pdf(dgd_grid, scale=scale_rayleigh)
pdf_maxwell = stats.maxwell.pdf(dgd_grid, scale=maxwell_a)
pdf_rayleigh_fit = stats.rayleigh.pdf(dgd_grid, scale=rayleigh_fit_scale)
pdf_maxwell_fit = stats.maxwell.pdf(dgd_grid, scale=maxwell_fit_a)

ax.plot(dgd_grid * 1e12, pdf_rayleigh * 1e12, '--', color='C1', linewidth=1.5,
        label=f'Rayleigh (params from code): mean={rayleigh_mean*1e12:.2f} ps')
ax.plot(dgd_grid * 1e12, pdf_maxwell * 1e12, '--', color='C2', linewidth=1.5,
        label=f'Maxwell (PMD theory): mean={maxwell_mean*1e12:.2f} ps')
ax.plot(dgd_grid * 1e12, pdf_rayleigh_fit * 1e12, '-', color='C1', linewidth=1,
        label='Rayleigh (fitted)')
ax.plot(dgd_grid * 1e12, pdf_maxwell_fit * 1e12, '-', color='C2', linewidth=1,
        label='Maxwell (fitted)')

ax.axvline(mean_dgd * 1e12, color='C0', linestyle=':', alpha=0.7, label=f'Mean = {mean_dgd*1e12:.3f} ps')
ax.axvline(target_mean * 1e12, color='gray', linestyle=':', alpha=0.5, label=f'Target RMS = {target_mean*1e12:.3f} ps')

ax.set_xlabel('DGD (ps)')
ax.set_ylabel('Probability density')
ax.set_title(f'PMD DGD distribution — {L_KM} km, PMD coeff = {PM_DISPERSION*1e12:.2f} ps/√km')
ax.legend(fontsize=7, loc='upper right')
ax.grid(True, alpha=0.3)

fig.tight_layout()
fname = f'val_pmd_dgd--seed{SEED}.png'
fig.savefig(os.path.join(OUT, fname), dpi=150, bbox_inches='tight')
print(f"\nSaved: {fname}")

# --- CSV ---
csv_name = f'val_pmd_dgd--seed{SEED}.csv'
header = 'dgd_ps'
np.savetxt(os.path.join(OUT, csv_name), dgd_samples.reshape(-1, 1) * 1e12,
           delimiter=',', header=header, comments='')
print(f"Saved: {csv_name}")
plt.close(fig)
