import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os, sys, argparse

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from src.channel.mzm import MZM
from src.channel.phase_modulator import PhaseModulator

OUT = os.path.join(os.path.dirname(__file__), '..', 'val_mzm')
os.makedirs(OUT, exist_ok=True)

parser = argparse.ArgumentParser()
parser.add_argument('--seed', type=int, default=42)
args = parser.parse_args()
SEED = args.seed
np.random.seed(SEED)

# Common test field (Ey-only for X-cut, Ex-only for Y-cut)
field_ey = np.array([0.0, 1.0], dtype=complex)
field_ex = np.array([1.0, 0.0], dtype=complex)
field_45 = np.array([1.0, 1.0], dtype=complex) / np.sqrt(2)

def tile_field(E_scalar, N):
    """Tile a (2,) field to (N, 2) for array modulation."""
    return np.tile(E_scalar, (N, 1))

mzm_pp = MZM(mode='push-pull')
mzm_sd = MZM(mode='single-drive')
V_pi = mzm_pp.V_pi

# ----------------------------------------------------------------
# Panel A: Transfer function — cos^2 characteristic
# ----------------------------------------------------------------
V_scan = np.linspace(0, 2 * V_pi, 200)
N_scan = len(V_scan)
E_pp_arr = mzm_pp.modulate(tile_field(field_ey, N_scan), V_scan)
E_sd_arr = mzm_sd.modulate(tile_field(field_ey, N_scan), V_scan)
P_pp = np.sum(np.abs(E_pp_arr)**2, axis=1)
P_sd = np.sum(np.abs(E_sd_arr)**2, axis=1)
P_in = 1.0
theory_pp = P_in * np.cos(np.pi * V_scan / (2 * V_pi))**2
theory_sd = theory_pp  # same intensity, different phase

# Test points for push-pull
idx_null = np.argmin(np.abs(V_scan - V_pi))
idx_peak = np.argmin(np.abs(V_scan - 0.0))
idx_quad = np.argmin(np.abs(V_scan - V_pi / 2))

# ----------------------------------------------------------------
# Panel B: Phase response — push-pull vs single-drive chirp
# ----------------------------------------------------------------
E_pp_ph = mzm_pp.modulate(tile_field(field_ey, N_scan), V_scan)
E_sd_ph = mzm_sd.modulate(tile_field(field_ey, N_scan), V_scan)

# Chirp phase: the exp(j*pi*V/(2*V_pi)) part only.
# E_out for single-drive = cos(pi*V/(2*V_pi)) * exp(j*pi*V/(2*V_pi))
# np.angle(E_out) mixes the cosine sign with the chirp.
# Remove the cosine sign: arg(cos) = 0 when cos > 0, pi when cos < 0.
phi_sd_raw = np.angle(E_sd_ph[:, 1])
cos_val = np.cos(np.pi * V_scan / (2 * V_pi))
cos_sign = np.where(cos_val < 0, np.pi, 0.0)
phi_sd_chirp = np.unwrap(phi_sd_raw - cos_sign)

# Mask null crossings where amplitude is near zero (phase undefined)
threshold = 0.01
amp_sd = np.abs(E_sd_ph[:, 1])
mask_sd = amp_sd > threshold * np.max(amp_sd)

# Push-pull analytic: zero chirp (cos is real, no frequency deviation)
# Single-drive analytic: chirp = pi * V / (2 * V_pi)
theory_sd_phase = np.pi * V_scan / (2 * V_pi)

# ----------------------------------------------------------------
# Panel C: Extinction ratio effect
# ----------------------------------------------------------------
V_scan2 = np.linspace(0, 2 * V_pi, 200)
figures = []
legends = []
N_scan2 = len(V_scan2)
for er_db in [None, 20, 10]:
    mzm_er = MZM(mode='push-pull', extinction_ratio_db=er_db)
    E_er = mzm_er.modulate(tile_field(field_ey, N_scan2), V_scan2)
    P_vals = np.sum(np.abs(E_er)**2, axis=1) / np.sum(np.abs(E_er)**2, axis=1).max()
    if er_db is None:
        figures.append(P_vals)
        legends.append('Ideal (infinite ER)')
    else:
        figures.append(P_vals)
        legends.append(f'ER = {er_db} dB')

er_null = [f[np.argmin(np.abs(V_scan2 - V_pi))] for f in figures]

# ----------------------------------------------------------------
# Panel D: Insertion loss
# ----------------------------------------------------------------
il_vals_db = np.linspace(0, 6, 7)
P_il_meas = []
P_il_theory = []
for il in il_vals_db:
    mzm_il = MZM(mode='push-pull', insertion_loss_db=il)
    E = mzm_il.modulate(field_ey, 0.0)
    P_il_meas.append(np.sum(np.abs(E)**2))
    P_il_theory.append(10.0 ** (-il / 10.0) * 1.0)

# ----------------------------------------------------------------
# Panel E: Crystal cut — X-cut modulates Ey, Y-cut modulates Ex
# ----------------------------------------------------------------
pm_x = PhaseModulator(crystal_cut='X')
pm_y = PhaseModulator(crystal_cut='Y')
mzm_x = MZM(pm=pm_x)
mzm_y = MZM(pm=pm_y)

V_pi_x = mzm_x.V_pi
V_pi_y = mzm_y.V_pi
V_test = np.linspace(0, 2 * V_pi_x, 100)

# X-cut: field with both Ex and Ey
E_out_x = np.array([mzm_x.modulate(field_45, V) for V in V_test])
Px_x = np.abs(E_out_x[:, 0])**2  # Ex power after X-cut (should be constant)
Py_x = np.abs(E_out_x[:, 1])**2  # Ey power after X-cut (should be modulated)

# Y-cut
E_out_y = np.array([mzm_y.modulate(field_45, V) for V in np.linspace(0, 2 * V_pi_y, 100)])
Px_y = np.abs(E_out_y[:, 0])**2  # Ex power after Y-cut (should be modulated)
Py_y = np.abs(E_out_y[:, 1])**2  # Ey power after Y-cut (should be constant)

# ================================================================
# 6-panel figure
# ================================================================
fig = plt.figure(figsize=(16, 12))
gs = fig.add_gridspec(2, 3, hspace=0.35, wspace=0.35)

# --- Panel A: Transfer function ---
ax1 = fig.add_subplot(gs[0, 0])
ax1.plot(V_scan / V_pi, P_pp, '-', c='C0', lw=1.5, label='Push-pull (sim)')
ax1.plot(V_scan / V_pi, theory_pp, '--k', lw=1.5, alpha=0.5, label=r'$\cos^2(\pi V/2V_\pi)$')
ax1.plot(V_scan / V_pi, P_sd, ':', c='C1', lw=1.5, label='Single-drive (sim)')
ax1.axvline(0.5, c='gray', ls=':', lw=0.6, alpha=0.5)
ax1.axvline(1.0, c='gray', ls=':', lw=0.6, alpha=0.5)
ax1.annotate('Quadrature\n(V$_{\\pi}$/2)\nP = 0.5', xy=(0.5, 0.5), fontsize=7,
             ha='center', va='bottom')
ax1.annotate('Null\n(V$_{\\pi}$)\nP = 0', xy=(1.0, 0.05), fontsize=7,
             ha='center', va='top')
ax1.annotate('Peak\nP = 1.0', xy=(0.0, 1.0), fontsize=7,
             ha='left', va='bottom', xytext=(0.15, 0.92),
             arrowprops=dict(arrowstyle='->', color='gray', lw=0.8))
ax1.set(xlabel=r'$V / V_\pi$', ylabel='Output power (a.u.)',
        title='A: Transfer function\nAgrawal [1] Sec 4.2')
ax1.legend(fontsize=7)
ax1.grid(True, alpha=0.25)

# --- Panel B: Phase response (chirp) ---
ax2 = fig.add_subplot(gs[0, 1])
ax2.axhline(0, c='C0', lw=1.5, label='Push-pull (theory)')
ax2.plot(V_scan[mask_sd] / V_pi, phi_sd_chirp[mask_sd], '-', c='C1', lw=1.5,
         label='Single-drive (sim)')
ax2.plot(V_scan / V_pi, theory_sd_phase, '--k', lw=1.5, alpha=0.5,
         label=r'$\pi V / 2V_\pi$ (theory)')
ax2.set(xlabel=r'$V / V_\pi$', ylabel='Chirp phase (rad)',
        title='B: Phase response (chirp)\nKoyama and Iga [2]')
ax2.legend(fontsize=7)
ax2.grid(True, alpha=0.25)

# --- Panel C: Extinction ratio ---
ax3 = fig.add_subplot(gs[0, 2])
for Pv, leg in zip(figures, legends):
    ax3.plot(V_scan2 / V_pi, Pv, lw=1.5, label=leg)
ax3.set(xlabel=r'$V / V_\pi$', ylabel='Normalized output',
        title='C: Extinction ratio effect')
ax3.legend(fontsize=7)
ax3.grid(True, alpha=0.25)

# --- Panel D: Insertion loss ---
ax4 = fig.add_subplot(gs[1, 0])
ax4.plot(il_vals_db, P_il_meas, 'o-', c='C3', ms=5, label='Measured')
ax4.plot(il_vals_db, P_il_theory, '--k', lw=1.5, alpha=0.5,
         label=r'$10^{-IL/10}$')
ax4.set(xlabel='Insertion loss (dB)', ylabel='Output power (a.u.)',
        title='D: Insertion loss')
ax4.legend(fontsize=7)
ax4.grid(True, alpha=0.25)

# --- Panel E: Crystal cut ---
ax5 = fig.add_subplot(gs[1, 1])
ax5.plot(V_test / V_pi_x, Px_x, '-', c='C0', lw=1.5, label='X-cut: |Ex|$^2$ (pass)')
ax5.plot(V_test / V_pi_x, Py_x, '--', c='C0', lw=1.5, label='X-cut: |Ey|$^2$ (mod)')
V_y_norm = np.linspace(0, 2 * V_pi_y, 100) / V_pi_y
ax5.plot(V_y_norm, Px_y, '-', c='C1', lw=1.5, alpha=0.7, label='Y-cut: |Ex|$^2$ (mod)')
ax5.plot(V_y_norm, Py_y, '--', c='C1', lw=1.5, alpha=0.7, label='Y-cut: |Ey|$^2$ (pass)')
ax5.set(xlabel=r'$V / V_\pi$', ylabel='Power (a.u.)',
        title='E: Crystal cut selection\nWeis and Gaylord [3]')
ax5.legend(fontsize=7)
ax5.grid(True, alpha=0.25)

fig.suptitle('MZM Validation -- Agrawal [1] Sec 4.2, Koyama and Iga [2], Weis and Gaylord [3]',
             fontsize=13, fontweight='bold', y=0.97)
fig.savefig(os.path.join(OUT, f'val_mzm--seed{SEED}.png'), dpi=200, bbox_inches='tight')
print(f"Saved: val_mzm--seed{SEED}.png")
csv_data = np.column_stack([V_scan / V_pi, P_pp, theory_pp, P_sd,
                            np.zeros_like(V_scan), phi_sd_chirp, theory_sd_phase])
np.savetxt(os.path.join(OUT, f'val_mzm--seed{SEED}.csv'),
           csv_data, delimiter=',',
           header='V_over_Vpi,P_pushpull,P_theory,P_singledrive,phi_pp(chirp),phi_sd(chirp),phi_sd_theory',
           comments='')
print(f"Saved: val_mzm--seed{SEED}.csv")

# ----------------------------------------------------------------
# The claims the table has always made, as numbers
# ----------------------------------------------------------------
# Every one of these was drawn against its analytic curve and then
# declared "verified" or given a bare ideal value, with nothing comparing
# them.  A modulator that had lost its cos^2 would have printed the same
# table.
err_transfer = float(np.max(np.abs(P_pp - theory_pp)))


# The peak, quadrature and null are evaluated AT their voltages rather
# than at the nearest sweep point.  V_scan is 200 points across
# [0, 2*V_pi], so it contains neither V_pi/2 nor V_pi: the nearest samples
# sit at 0.5025*V_pi and 0.9950*V_pi, where cos^2 is 0.49605 and 6.2e-5.
# The table reported those as "P_out = 0.5" and "P_out = 0", which were
# claims about points the sweep never visited.
def _power_at(volts):
    return float(np.sum(np.abs(mzm_pp.modulate(field_ey, volts)) ** 2))


P_peak = _power_at(0.0)
P_quad = _power_at(V_pi / 2.0)
P_null = _power_at(V_pi)

# Push-pull is chirp-free by construction: both arms move oppositely, so
# the transfer is real and the phase can only be 0 or pi.  Measured on the
# same sign-corrected, null-masked phase the single-drive panel uses.
phi_pp_raw = np.angle(E_pp_ph[:, 1])
_cos_val = np.cos(np.pi * V_scan / (2 * V_pi))
_cos_sign = np.where(_cos_val < 0, np.pi, 0.0)
phi_pp_chirp = np.unwrap(phi_pp_raw - _cos_sign)
amp_pp = np.abs(E_pp_ph[:, 1])
mask_pp = amp_pp > threshold * np.max(amp_pp)
chirp_pp_max = float(np.max(np.abs(phi_pp_chirp[mask_pp])))

err_chirp = float(np.max(np.abs(phi_sd_chirp[mask_sd]
                                - theory_sd_phase[mask_sd])))
err_il = float(np.max(np.abs(np.array(P_il_meas)
                             - np.array(P_il_theory))))

# Crystal cut: the unmodulated component must stay flat across the sweep.
# Peak-to-peak relative to its own mean, so "flat" is scale-free.
def _flatness(p):
    p = np.asarray(p, float)
    return float((p.max() - p.min()) / max(p.mean(), 1e-300))


leak_x = _flatness(Px_x)     # X-cut must leave Ex alone
leak_y = _flatness(Py_y)     # Y-cut must leave Ey alone

import csv
table_csv = os.path.join(OUT, f'val_mzm--seed{SEED}_table.csv')
with open(table_csv, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Parameter', 'Value'])
    writer.writerow(['Transfer function',
                     f'cos^2(pi*V/(2*V_pi)), max err {err_transfer:.2e}'])
    writer.writerow(['Peak (V=0)', f'P_out = {P_peak:.8f}'])
    writer.writerow(['Quadrature (V=V_pi/2)', f'P_out = {P_quad:.8f}'])
    writer.writerow(['Null (V=V_pi)', f'P_out = {P_null:.3e}'])
    writer.writerow(['Chirp (push-pull)',
                     f'max |phi| = {chirp_pp_max:.2e} rad'])
    writer.writerow(['Chirp (single-drive)',
                     f'phi = pi*V/(2*V_pi), max err {err_chirp:.2e} rad'])
    writer.writerow(['ER degradation', f'10 dB ER: null depth = {er_null[2]:.3f}'])
    writer.writerow(['Insertion loss',
                     f'10^(-IL/10), max err {err_il:.2e}'])
    writer.writerow(['Crystal cut',
                     f'X-cut leaves Ex flat to {leak_x:.2e}; '
                     f'Y-cut leaves Ey flat to {leak_y:.2e}'])
print(f"Saved: val_mzm--seed{SEED}_table.csv")
plt.close(fig)

# ============================================================
# Verdict
# ============================================================
#
# All closed form against closed form, so the gates sit far above the
# numerical residual and catch the gross failures instead: a transfer that
# has stopped being cos^2, a push-pull arm that has acquired chirp, an
# insertion loss applied in amplitude rather than power, a crystal cut
# modulating the wrong component.
TOL_POWER = 1e-9        # absolute, on a transfer normalised to 1
TOL_CHIRP_RAD = 1e-9    # absolute, radians
TOL_FLAT = 1e-9         # fractional peak-to-peak on the untouched axis

failures = []

if err_transfer >= TOL_POWER:
    failures.append(
        f"push-pull transfer departs from cos^2(pi*V/2V_pi) by "
        f"{err_transfer:.3e}, past {TOL_POWER:g}")

for label, got, want in (('peak (V=0)', P_peak, 1.0),
                         ('quadrature (V=V_pi/2)', P_quad, 0.5),
                         ('null (V=V_pi)', P_null, 0.0)):
    if abs(got - want) >= 1e-6:
        failures.append(
            f"{label}: P_out = {got:.8f} against {want}, which is where "
            f"cos^2 puts it")

if chirp_pp_max >= TOL_CHIRP_RAD:
    failures.append(
        f"push-pull is not chirp-free: max |phi| = {chirp_pp_max:.3e} rad. "
        f"Both arms move oppositely, so the transfer is real and the phase "
        f"can only be 0 or pi")

if err_chirp >= TOL_CHIRP_RAD:
    failures.append(
        f"single-drive chirp departs from pi*V/(2*V_pi) by "
        f"{err_chirp:.3e} rad, past {TOL_CHIRP_RAD:g}")

if err_il >= TOL_POWER:
    failures.append(
        f"insertion loss departs from 10^(-IL/10) by {err_il:.3e}. Applied "
        f"in amplitude instead of power it would be out by a square")

if leak_x >= TOL_FLAT:
    failures.append(
        f"an X-cut modulator disturbed Ex by {leak_x:.3e} peak-to-peak; "
        f"X-cut drives Ey and must leave Ex alone")
if leak_y >= TOL_FLAT:
    failures.append(
        f"a Y-cut modulator disturbed Ey by {leak_y:.3e} peak-to-peak; "
        f"Y-cut drives Ex and must leave Ey alone")

print()
if failures:
    print("[FAIL]")
    for f_ in failures:
        print(f"  - {f_}")
    sys.exit(1)
print(f"[PASS] push-pull transfer is cos^2(pi*V/2V_pi) to "
      f"{err_transfer:.2e}, with peak/quadrature/null on their exact values")
print(f"[PASS] push-pull carries no chirp ({chirp_pp_max:.2e} rad) while "
      f"single-drive follows pi*V/2V_pi to {err_chirp:.2e} rad")
print(f"[PASS] insertion loss scales as 10^(-IL/10) to {err_il:.2e}")
print(f"[PASS] each crystal cut leaves the other axis flat "
      f"({max(leak_x, leak_y):.2e} peak-to-peak)")
sys.exit(0)
