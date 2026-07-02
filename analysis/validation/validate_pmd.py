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

# --- replicate fiber.py DGD generation (updated: Maxwellian) ---
pmd_sd = PM_DISPERSION * np.sqrt(L_M)          # seconds, RMS DGD
maxwell_scale = pmd_sd / np.sqrt(3)
dgd_samples = stats.maxwell.rvs(scale=maxwell_scale, size=N_REALIZATIONS)

# --- statistics ---
mean_dgd = np.mean(dgd_samples)
rms_dgd = np.sqrt(np.mean(dgd_samples**2))
# Maxwell(scale=a): mean = 2a·√(2/π), RMS = a·√3
# With a = pmd_sd/√3: RMS = pmd_sd ✓, mean = 2·pmd_sd/√3·√(2/π) ≈ 0.921·pmd_sd
expected_mean = 2 * maxwell_scale * np.sqrt(2 / np.pi)
expected_rms = maxwell_scale * np.sqrt(3)

print(f"PMD DGD distribution validation ({N_REALIZATIONS} realizations)")
print(f"  PMD coeff = {PM_DISPERSION*1e12:.2f} ps/sqrt(m)")
print(f"  L = {L_KM} km")
print(f"  Target RMS DGD (pmd_sd) = {pmd_sd*1e12:.4f} ps")
print(f"  Maxwell scale a = {maxwell_scale*1e12:.4f} ps  (a = pmd_sd/sqrt(3))")
print()
print(f"Current model (Maxwellian):")
print(f"  Measured mean DGD = {mean_dgd*1e12:.4f} ps  (expected = {expected_mean*1e12:.4f})")
print(f"  Measured RMS DGD  = {rms_dgd*1e12:.4f} ps  (expected = {expected_rms*1e12:.4f})")
print(f"  RMS/mean = {rms_dgd/mean_dgd:.4f}  (Maxwell theory: {np.sqrt(3/(8/np.pi)):.4f})")

# --- KS test ---
ks = stats.kstest(dgd_samples, lambda x: stats.maxwell.cdf(x, scale=maxwell_scale))
print()
print(f"KS test vs Maxwell:  statistic={ks.statistic:.6f}, p-value={ks.pvalue:.4f}")

# --- plot ---
fig, ax = plt.subplots(figsize=(8, 5))
bins = np.linspace(0, max(dgd_samples)*1.05, 80)
ax.hist(dgd_samples * 1e12, bins=bins, density=True, alpha=0.6, color='C0',
        label=f'Simulated ({N_REALIZATIONS} draws)')

dgd_grid = np.linspace(0, max(dgd_samples)*1.05, 500)
pdf = stats.maxwell.pdf(dgd_grid, scale=maxwell_scale)
ax.plot(dgd_grid * 1e12, pdf * 1e12, '-', color='C2', linewidth=1.5,
        label=f'Maxwellian: mean={expected_mean*1e12:.2f} ps, RMS={expected_rms*1e12:.2f} ps')

ax.axvline(mean_dgd * 1e12, color='C0', linestyle=':', alpha=0.7, label=f'Sim mean = {mean_dgd*1e12:.3f} ps')
ax.axvline(pmd_sd * 1e12, color='gray', linestyle=':', alpha=0.5, label=f'Target RMS = {pmd_sd*1e12:.3f} ps')

ax.set_xlabel('DGD (ps)')
ax.set_ylabel('Probability density')
ax.set_title(f'PMD DGD distribution — {L_KM} km, PMD coeff = {PM_DISPERSION*1e12:.2f} ps/√km')
ax.legend(fontsize=8, loc='upper right')
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
