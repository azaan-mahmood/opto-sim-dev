import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import os, sys, argparse

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from src.channel.fiber import apply_pmd

OUT = os.path.join(os.path.dirname(__file__), '..', 'val_pmd')
os.makedirs(OUT, exist_ok=True)

parser = argparse.ArgumentParser()
parser.add_argument('--seed', type=int, default=42)
parser.add_argument('--n-realizations', type=int, default=10000)
args = parser.parse_args()
SEED = args.seed
N_REAL = args.n_realizations
np.random.seed(SEED)

L_KM = 100
PMD_COEFF = 0.1  # ps/sqrt(km) — Corning SMF-28 Ultra spec <= 0.1 ps/sqrt(km)

# Gates for the verdict at the bottom of this file.  Stated here so the
# "Expected" column of the committed table and the checks are the same
# numbers rather than two opinions about them.
KS_P_FLOOR = 1e-3        # Maxwell shape; permissive, a true sample rarely
                         # falls below even at 10k draws
COEFF_TOL_FRAC = 0.02    # fitted vs nominal ps/sqrt(km), 2000 draws per
                         # point across seven lengths
MIN_R2 = 0.999           # RMS DGD linear in sqrt(L)
DT = 0.5e-12
N_SAMPLES = 4096
T = np.arange(N_SAMPLES) * DT
center = N_SAMPLES // 2

PULSE_SIGMA = 5e-12
t_arr = T - T[center]
E_in = np.zeros((N_SAMPLES, 2), dtype=np.complex128)
E_in[:, 0] = np.exp(-0.5 * (t_arr / PULSE_SIGMA)**2)
E_in[:, 1] = E_in[:, 0].copy()

pmd_sd = PMD_COEFF * 1e-12 * np.sqrt(L_KM)   # seconds
maxwell_scale = pmd_sd / np.sqrt(3)

print(f"PMD validation via apply_pmd - {N_REAL} realizations, {L_KM} km")
print(f"  Expected RMS DGD = {pmd_sd*1e12:.3f} ps")
print(f"  Expected mean DGD = {2*maxwell_scale*np.sqrt(2/np.pi)*1e12:.3f} ps")

dgds = []
for i in range(N_REAL):
    _, dgd = apply_pmd(E_in.copy(), dt=DT, L=L_KM*1000, pmd_coeff_ps_sqrt_km=PMD_COEFF)
    dgds.append(dgd)
    if (i + 1) % 2000 == 0:
        print(f"  {i+1}/{N_REAL} complete")

dgds_ps = np.array(dgds) * 1e12
mean_dgd = np.mean(dgds_ps)
rms_dgd = np.sqrt(np.mean(dgds_ps**2))
expected_mean = 2 * maxwell_scale * np.sqrt(2 / np.pi) * 1e12
expected_rms = pmd_sd * 1e12

ks = stats.kstest(dgds, lambda x: stats.maxwell.cdf(x, scale=maxwell_scale))
print(f"\nResults:")
print(f"  Mean DGD = {mean_dgd:.3f} ps  (expected = {expected_mean:.3f} ps)")
print(f"  RMS  DGD = {rms_dgd:.3f} ps  (expected = {expected_rms:.3f} ps)")
print(f"  KS test vs Maxwell:  D={ks.statistic:.5f}, p={ks.pvalue:.4f}")

# ---- DGD vs sqrt(L) sweep ----
L_sweep_km = np.array([10, 25, 50, 75, 100, 150, 200])
N_PER_L = 2000
rms_dgd_sweep = np.zeros(len(L_sweep_km))
for j, Lk in enumerate(L_sweep_km):
    dgd_list = []
    for _ in range(N_PER_L):
        _, dgd = apply_pmd(E_in.copy(), dt=DT, L=Lk*1000, pmd_coeff_ps_sqrt_km=PMD_COEFF)
        dgd_list.append(dgd)
    rms_dgd_sweep[j] = np.sqrt(np.mean(np.array(dgd_list)**2)) * 1e12

sqrt_L = np.sqrt(L_sweep_km * 1e3)
slope_pmd, int_pmd, r2_pmd, p_pmd, se_pmd = stats.linregress(sqrt_L, rms_dgd_sweep)
print(f"\nDGD vs sqrt(L) sweep:")
print(f"  Fitted PMD coeff = {slope_pmd*np.sqrt(1e3):.5f} ps/sqrt(km) (nominal = {PMD_COEFF:.5f} ps/sqrt(km))")
print(f"  R^2 = {r2_pmd:.6f}")

# ============================================================
# Main validation figure — 6 panels
# ============================================================
fig = plt.figure(figsize=(16, 12))
gs = fig.add_gridspec(2, 3, hspace=0.35, wspace=0.35)

# --- Panel A: DGD histogram + Maxwellian fit (Razavi Fig 2.11) ---
ax1 = fig.add_subplot(gs[0, 0])
bins = np.linspace(0, max(dgds_ps)*1.05, 50)
ax1.hist(dgds_ps, bins=bins, density=True, alpha=0.6, color='C0',
         label=f'apply_pmd ({len(dgds_ps)} realizations)')

dgd_plot = np.linspace(0, max(dgds_ps)*1.05, 500)
pdf_theory = stats.maxwell.pdf(dgd_plot * 1e-12, scale=maxwell_scale) * 1e-12
ax1.plot(dgd_plot, pdf_theory, '-', c='k', lw=1.5,
         label=f'Maxwellian (RMS = {expected_rms:.1f} ps)')
ax1.axvline(mean_dgd, c='C0', ls=':', lw=1, alpha=0.7,
            label=f'Measured mean = {mean_dgd:.2f} ps')
ax1.axvline(expected_rms, c='gray', ls='--', lw=1, alpha=0.5,
            label=f'Expected RMS = {expected_rms:.1f} ps')
ax1.set(xlabel='DGD (ps)', ylabel='Probability density',
        title='A: DGD distribution\n(Razavi [5] Fig 2.11)')
ax1.legend(fontsize=7)
ax1.grid(True, alpha=0.25)
ax1.annotate(f'KS stat D = {ks.statistic:.5f}\np-value = {ks.pvalue:.4f}',
             xy=(0.95, 0.95), xycoords='axes fraction', ha='right', va='top',
             fontsize=8, bbox=dict(boxstyle='round', fc='lightyellow', alpha=0.8))

# --- Panel B: QQ plot for Maxwellian ---
ax2 = fig.add_subplot(gs[0, 1])
dgds_sorted = np.sort(dgds_ps)
theoretical_quantiles = stats.maxwell.ppf(
    np.linspace(1/(N_REAL+1), N_REAL/(N_REAL+1), N_REAL),
    scale=maxwell_scale) * 1e12
ax2.plot(theoretical_quantiles, dgds_sorted, '.', ms=1, c='C0', alpha=0.5)
q_min = min(theoretical_quantiles.min(), dgds_sorted.min())
q_max = max(theoretical_quantiles.max(), dgds_sorted.max())
ax2.plot([q_min, q_max], [q_min, q_max], '--k', lw=1, label='y = x (ideal)')
ax2.set(xlabel='Theoretical Maxwellian quantiles (ps)',
        ylabel='Sample quantiles (ps)',
        title='B: QQ-plot (Maxwellian fit)')
ax2.legend(fontsize=8)
ax2.grid(True, alpha=0.25)
ax2.axis('equal')
corr = np.corrcoef(theoretical_quantiles, dgds_sorted)[0, 1]
ax2.annotate(f'r = {corr:.6f}', xy=(0.05, 0.95), xycoords='axes fraction',
             ha='left', va='top', fontsize=8,
             bbox=dict(boxstyle='round', fc='lightyellow', alpha=0.8))

# --- Panel C: Residual (observed - expected) ---
ax3 = fig.add_subplot(gs[0, 2])
n_bins = 30
counts, bin_edges = np.histogram(dgds_ps, bins=n_bins, density=True)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
pdf_at_bins = stats.maxwell.pdf(bin_centers * 1e-12, scale=maxwell_scale) * 1e-12
residuals = counts - pdf_at_bins
ax3.bar(bin_centers, residuals, width=bin_centers[1]-bin_centers[0],
        alpha=0.6, color='C3')
ax3.axhline(0, c='k', lw=0.5)
ax3.set(xlabel='DGD (ps)', ylabel='Residual (data − theory)',
        title='C: Histogram residual')
ax3.grid(True, alpha=0.25)

# --- Panel D: DGD vs sqrt(L) (PMD scaling, Razavi [5] Sec 2.5) ---
sqrt_km = np.sqrt(L_sweep_km)
slope_ps_per_sqrtkm = slope_pmd * np.sqrt(1e3)
nominal_ps_per_sqrtkm = PMD_COEFF

ax4 = fig.add_subplot(gs[1, 0])
ax4.plot(sqrt_km, rms_dgd_sweep, 'o-', ms=5, c='C2', lw=1.5,
         label='RMS DGD from apply_pmd')
sqrt_km_dense = np.linspace(0, max(sqrt_km), 200)
ax4.plot(sqrt_km_dense, slope_ps_per_sqrtkm * sqrt_km_dense, '--k', lw=1, alpha=0.6,
         label=f'Fit: DGD_rms = {slope_ps_per_sqrtkm:.4f} * sqrt(L) (ps)')
ax4.plot(sqrt_km_dense, nominal_ps_per_sqrtkm * sqrt_km_dense, ':', c='gray', lw=1, alpha=0.5,
         label=f'Nominal: DGD_rms = {nominal_ps_per_sqrtkm:.4f} * sqrt(L) (ps)')
ax4.set(xlabel='sqrt(L) (sqrt(km))', ylabel='RMS DGD (ps)',
        title=f'D: PMD scaling DGD ~ sqrt(L)\n(Razavi [5] Sec 2.5)')
ax4.legend(fontsize=7)
ax4.grid(True, alpha=0.25)
ax4.annotate(f'R^2 = {r2_pmd:.6f}',
             xy=(0.95, 0.05), xycoords='axes fraction', ha='right', va='bottom',
             fontsize=8, bbox=dict(boxstyle='round', fc='lightyellow', alpha=0.8))

# --- Panel E: PMD coeff vs distance (consistency) ---
pmd_coeff_extracted = rms_dgd_sweep / np.sqrt(L_sweep_km)
ax5 = fig.add_subplot(gs[1, 1])
ax5.plot(L_sweep_km, pmd_coeff_extracted, 'o-', ms=4, c='C4', lw=1.2)
ax5.axhline(PMD_COEFF, c='gray', ls='--', lw=1, alpha=0.6,
            label=f'Nominal PMD coeff = {PMD_COEFF:.3f} ps/√km')
ax5.set(xlabel='Fibre length (km)',
        ylabel='Extracted PMD coeff (ps/√km)',
        title='E: PMD coefficient consistency',
        ylim=(0, max(pmd_coeff_extracted)*1.3))
ax5.legend(fontsize=8)
ax5.grid(True, alpha=0.25)

fig.suptitle('Polarization Mode Dispersion - Validation vs Razavi [5] Sec 2.5',
             fontsize=13, fontweight='bold', y=0.97)
fig.savefig(os.path.join(OUT, f'val_pmd_dgd--seed{SEED}.png'), dpi=200, bbox_inches='tight')
print(f"Saved: val_pmd_dgd--seed{SEED}.png")

csv_name = f'val_pmd_dgd--seed{SEED}.csv'
np.savetxt(os.path.join(OUT, csv_name), dgds_ps, delimiter=',',
           header='dgd_ps_from_apply_pmd', comments='')
print(f"Saved: {csv_name}")

import csv
table_csv = os.path.join(OUT, f'val_pmd--seed{SEED}_table.csv')
with open(table_csv, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Parameter', 'Value', 'Expected'])
    writer.writerow(['Mean DGD (ps)', f'{mean_dgd:.3f}', f'{expected_mean:.3f}'])
    writer.writerow(['RMS DGD (ps)', f'{rms_dgd:.3f}', f'{expected_rms:.3f}'])
    writer.writerow(['KS D-stat', f'{ks.statistic:.5f}', '< 0.05 typical'])
    writer.writerow(['KS p-value', f'{ks.pvalue:.4f}',
                     f'> {KS_P_FLOOR:g} (not rejected)'])
    writer.writerow(['PMD coeff (ps/sqrt(km))', f'{slope_ps_per_sqrtkm:.4f}',
                     f'{nominal_ps_per_sqrtkm:.4f}'])
    # r_value, squared.  `linregress` returns r as its third value; this
    # row has been labelled R2 while carrying r.  Both sit near 1 here, so
    # the mislabel never showed.
    writer.writerow(['DGD sqrt(L) R2', f'{r2_pmd ** 2:.6f}', '~ 1.0'])
print(f"Saved: val_pmd--seed{SEED}_table.csv")
plt.close(fig)

# ============================================================
# Verdict
# ============================================================
#
# This one is a Monte Carlo, not a closed form, so the tolerances are
# derived from the sampling distribution rather than chosen.  Maxwell with
# scale `a` has
#
#     mean = 2a*sqrt(2/pi),  var = a^2*(3pi-8)/pi,  m2 = 3a^2,  m4 = 15a^4
#
# so the standard error of the sample mean is sd/sqrt(N), and the relative
# standard error of the sample RMS is sqrt((m4/m2^2 - 1)/(4N)) = 1/sqrt(6N).
# Both therefore tighten with `--n-realizations` instead of being fixed at
# whatever one run happened to produce.
SIGMA_GATE = 5.0

sd_maxwell = maxwell_scale * np.sqrt((3 * np.pi - 8) / np.pi) * 1e12
se_mean = sd_maxwell / np.sqrt(N_REAL)
se_rms_rel = 1.0 / np.sqrt(6.0 * N_REAL)

failures = []

dev_mean = abs(mean_dgd - expected_mean)
if dev_mean > SIGMA_GATE * se_mean:
    failures.append(
        f"mean DGD {mean_dgd:.4f} ps against an expected "
        f"{expected_mean:.4f} ps -- {dev_mean / se_mean:.1f} sigma on a "
        f"standard error of {se_mean:.4f} ps from {N_REAL} realizations")

dev_rms_rel = abs(rms_dgd - expected_rms) / expected_rms
if dev_rms_rel > SIGMA_GATE * se_rms_rel:
    failures.append(
        f"RMS DGD {rms_dgd:.4f} ps against an expected {expected_rms:.4f} "
        f"ps -- {dev_rms_rel / se_rms_rel:.1f} sigma. RMS DGD is what the "
        f"PMD coefficient is defined by, so this is the number the spec "
        f"names")

if ks.pvalue < KS_P_FLOOR:
    failures.append(
        f"the DGD sample is not Maxwellian: KS D = {ks.statistic:.5f}, "
        f"p = {ks.pvalue:.2e} against a {KS_P_FLOOR:g} floor. A random "
        f"birefringence walk gives a Maxwell DGD; another shape means the "
        f"walk is not random")

coeff_dev = abs(slope_ps_per_sqrtkm - nominal_ps_per_sqrtkm) \
    / nominal_ps_per_sqrtkm
if coeff_dev > COEFF_TOL_FRAC:
    failures.append(
        f"the fitted PMD coefficient is {slope_ps_per_sqrtkm:.5f} "
        f"ps/sqrt(km) against the {nominal_ps_per_sqrtkm:.5f} asked for, "
        f"a relative error of {coeff_dev:.3%}. The sweep runs over the "
        f"same model, so it has to return its own input")

if r2_pmd ** 2 < MIN_R2:
    failures.append(
        f"RMS DGD is not linear in sqrt(L) (R^2 = {r2_pmd ** 2:.6f}). "
        f"PMD accumulates as a random walk, which is what makes the "
        f"coefficient a per-sqrt(km) number at all")

print()
if failures:
    print("[FAIL]")
    for f_ in failures:
        print(f"  - {f_}")
    sys.exit(1)
print(f"[PASS] mean and RMS DGD sit within {SIGMA_GATE:g} sigma of Maxwell "
      f"over {N_REAL} realizations")
print(f"[PASS] the DGD sample is Maxwellian (KS p = {ks.pvalue:.3f})")
print(f"[PASS] RMS DGD grows as sqrt(L) (R^2 = {r2_pmd ** 2:.6f}) and the "
      f"fit returns {slope_ps_per_sqrtkm:.5f} ps/sqrt(km)")
sys.exit(0)
