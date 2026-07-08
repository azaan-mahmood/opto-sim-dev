import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import os, sys, argparse

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from src.channel.fiber import cable, _dgd_sampled

OUT = os.path.join(os.path.dirname(__file__), '..', 'val_pmd')
os.makedirs(OUT, exist_ok=True)

parser = argparse.ArgumentParser()
parser.add_argument('--seed', type=int, default=42)
parser.add_argument('--n-realizations', type=int, default=5000)
args = parser.parse_args()
SEED = args.seed
N_REAL = args.n_realizations
np.random.seed(SEED)

L_KM = 100
PMD_COEFF = 0.1e-12
DT = 0.5e-12
N_SAMPLES = 4096
T = np.arange(N_SAMPLES) * DT
center = N_SAMPLES // 2

PULSE_SIGMA = 5e-12
t_arr = T - T[center]
E_in = np.zeros((N_SAMPLES, 2), dtype=np.complex128)
E_in[:, 0] = np.exp(-0.5 * (t_arr / PULSE_SIGMA)**2)
E_in[:, 1] = E_in[:, 0].copy()

pmd_sd = PMD_COEFF * np.sqrt(L_KM * 1e3)
maxwell_scale = pmd_sd / np.sqrt(3)

print(f"PMD validation via cable() — {N_REAL} realizations, {L_KM} km")
print(f"  Expected RMS DGD = {pmd_sd*1e12:.2f} ps")
print(f"  Expected mean DGD = {2*maxwell_scale*np.sqrt(2/np.pi)*1e12:.2f} ps")
print()

# Clear any prior recorded DGDs, then call cable() — each call appends its DGD
_dgd_sampled.clear()
for i in range(N_REAL):
    E_out = cable(L_KM, E_in.copy(), dt=DT, wavelength=1550e-9,
                  dispersion=True, pm_dispersion=PMD_COEFF,
                  attenuation_factor=0.0, temperature=25.0, bend_radius=None)
    if (i + 1) % 1000 == 0:
        print(f"  {i+1}/{N_REAL} complete")

dgds = np.array(_dgd_sampled) * 1e12  # convert to ps, exact values from cable()
mean_dgd = np.mean(dgds)
rms_dgd = np.sqrt(np.mean(dgds**2))
expected_mean = 2 * maxwell_scale * np.sqrt(2 / np.pi) * 1e12
expected_rms = pmd_sd * 1e12

print(f"\nResults:")
print(f"  Mean DGD = {mean_dgd:.3f} ps  (expected = {expected_mean:.3f} ps)")
print(f"  RMS  DGD = {rms_dgd:.3f} ps  (expected = {expected_rms:.3f} ps)")

ks = stats.kstest(dgds, lambda x: stats.maxwell.cdf(x * 1e-12, scale=maxwell_scale))
print(f"  KS test vs Maxwell:  D={ks.statistic:.5f}, p={ks.pvalue:.4f}")

# --- Plot ---
fig, ax = plt.subplots(figsize=(8, 5))
bins = np.linspace(0, max(dgds)*1.05, 60)
ax.hist(dgds, bins=bins, density=True, alpha=0.6, color='C0',
        label=f'cable() extracted DGD ({len(dgds)} valid)')

dgd_ps = np.linspace(0, max(dgds)*1.05, 500)
pdf = stats.maxwell.pdf(dgd_ps * 1e-12, scale=maxwell_scale) * 1e-12
ax.plot(dgd_ps, pdf, '-', c='C2', lw=1.5,
        label=f'Maxwellian (RMS = {expected_rms:.1f} ps)')
ax.axvline(mean_dgd, c='C0', ls=':', alpha=0.7, label=f'Mean = {mean_dgd:.2f} ps')
ax.axvline(expected_rms, c='gray', ls=':', alpha=0.5, label=f'Target RMS = {expected_rms:.1f} ps')

ax.set(xlabel='DGD (ps)', ylabel='Probability density',
       title=f'PMD: DGD extracted from cable() output — {L_KM} km')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

fig.tight_layout()
fname = f'val_pmd_dgd--seed{SEED}.png'
fig.savefig(os.path.join(OUT, fname), dpi=150)
print(f"\nSaved: {fname}")

csv_name = f'val_pmd_dgd--seed{SEED}.csv'
np.savetxt(os.path.join(OUT, csv_name), dgds, delimiter=',',
           header='dgd_ps_extracted_from_cable', comments='')
print(f"Saved: {csv_name}")
plt.close(fig)
