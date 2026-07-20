"""
System-Level BB84 QKD Demonstration
====================================
Demonstrates that independently validated components compose correctly by
simulating BB84 QKD under combined impairments (CD, PMD, birefringence,
attenuation) with realistic laser, modulator, and detector models.

Uses MZM-carved Gaussian pulses (matching bb84_test_dispersion.py) to
provide sufficient spectral bandwidth for CD and PMD to induce measurable
QBER.

Produces a 5-panel figure:
  A: QBER vs fibre distance (0--1000 km)
  B: QBER vs pulse width (5--50 ps) at the critical distance (75 km)
  C: QBER vs fibre temperature (0--60 C) at 75 km
  D: QBER vs PMD coefficient (0.0--0.3 ps/sqrt(km)) at 75 km
  E: QBER vs bend radius (2 mm -- 5 cm) at 75 km

All panels use seeded RNG for reproducibility.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os, sys, time, argparse, random

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from src.lasers.cwlaser import CWLaser
from src.channel import propagate, optics
from src.channel.phase_modulator import PhaseModulator
from src.channel.mzm import MZM
from src.detectors.apd import apd

OUT = os.path.join(os.path.dirname(__file__), '..', 'val_system')
os.makedirs(OUT, exist_ok=True)

parser = argparse.ArgumentParser()
parser.add_argument('--seed', type=int, default=42)
parser.add_argument('--bits', type=int, default=3000,
                    help='Bits per QBER measurement')
args = parser.parse_args()
SEED = args.seed
NUM_BITS = args.bits

# ── Effective bit rate ──────────────────────────────────────────────
# Each bit occupies n_samples * dt = 4000 * 1e-12 = 4 ns → 250 MHz.
# This is well within the demonstrated range of high-speed BB84:
#   Takesue et al. (2007, Nature Photonics 1, 343):  1.6 GHz, 200 km
#   Dixon et al. (2008, APL 93, 131107):              GHz clock
#   Namekata et al. (2010, Opt. Express 18, 17237):   1.5 GHz, 85 km
BIT_RATE_MHZ = 250


def simulate_bb84_full(num_bits, fiber_length=50, pulse_sigma=30e-12,
                       dispersion=True, temperature=25, bend_radius=None,
                       pm_dispersion=0.1e-12, seed=None, n_samples=4000):
    """Full BB84 with MZM-carved pulses, all impairments."""
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)

    alice_bits, alice_bases = [], []
    bob_bits, bob_bases = [], []

    detector = apd(wavelength=1550e-9, quantum_efficiency=0.9, gain=10,
                   excess_noise_factor=10, load_resistance=50, temperature=25)

    alice_laser = CWLaser(
        wavelength=1550e-9, polarization_azimuth=np.pi,
        polarization_ellipticity=np.pi / 2, power_dbm=10,
        linewidth=1e6, rin_density=-140,
    )

    pm_alice = PhaseModulator(crystal_cut='X', modulation="DC")
    pm_bob = PhaseModulator(crystal_cut='X', modulation="DC")
    Vpi = pm_alice.Vpi

    dt = 1e-12
    t = np.arange(n_samples) * dt
    t_center = n_samples * dt / 2

    mzm = MZM(mode='push-pull', bias_voltage=Vpi)
    V_pulse = Vpi * np.exp(-0.5 * ((t - t_center) / pulse_sigma) ** 2)

    for _ in range(num_bits):
        alice_basis = random.choice(['C', 'X'])
        alice_bit = random.randint(0, 1)

        if alice_basis == 'C':
            phase_alice = Vpi / 2 if alice_bit == 0 else 3 * Vpi / 2
        else:
            phase_alice = 0 if alice_bit == 0 else Vpi

        E = alice_laser.sample_field(dt=dt, n_samples=n_samples)
        E = mzm.modulate(E_in=E, V=V_pulse)
        E = optics.polarizer(E, '45')
        E = pm_alice.modulate(E_field=E, V=phase_alice)

        E = propagate(
            fiber_length=fiber_length, E=E, dt=dt, dispersion=dispersion,
            attenuation_factor=0.182, temperature=temperature,
            bend_radius=bend_radius, pm_dispersion=pm_dispersion,
        )

        bob_basis = random.choice(['C', 'X'])
        phase_bob = 0 if bob_basis == 'C' else Vpi / 2
        E = pm_bob.modulate(E_field=E, V=phase_bob)

        Ex, Ey = optics.pbs(E)
        I_x = detector.output(E=Ex, bandwidth=1e6)
        I_y = detector.output(E=Ey, bandwidth=1e6)

        noise_floor = detector.calculate_noise(0, 1e6)
        threshold = 3 * noise_floor

        if I_x > threshold or I_y > threshold:
            bob_bit = 0 if I_x > I_y else 1
        else:
            bob_bit = random.randint(0, 1)

        alice_bits.append(alice_bit)
        alice_bases.append(alice_basis)
        bob_bits.append(bob_bit)
        bob_bases.append(bob_basis)

    sifted_indices = [i for i in range(num_bits) if alice_bases[i] == bob_bases[i]]
    sifted_alice = [alice_bits[i] for i in sifted_indices]
    sifted_bob = [bob_bits[i] for i in sifted_indices]

    qber = 0.0
    if len(sifted_alice) > 0:
        errors = sum(sifted_alice[i] != sifted_bob[i] for i in range(len(sifted_alice)))
        qber = errors / len(sifted_alice)

    return qber


def sweep(pname, vals, fixed):
    base = dict(fiber_length=50, pulse_sigma=30e-12, dispersion=True)
    qs = []
    for v in vals:
        kw = {**base, **fixed, pname: v}
        q = simulate_bb84_full(NUM_BITS, seed=SEED, **kw)
        qs.append(q)
        print(f"  {pname}={v:<12}  QBER={q*100:5.1f}%")
    return np.array(qs)


# ── Panel A: QBER vs distance ─────────────────────────────────────
print("Panel A: QBER vs distance (30 ps pulses, 250 MHz effective bit rate)")
# Fine grid up to 200 km, coarse grid beyond (saturation + dark-count regime)
distances = np.concatenate([np.arange(0, 201, 10),
                             np.arange(250, 1001, 50)])
qber_on = sweep('fiber_length', distances, dict(dispersion=True))
qber_off = sweep('fiber_length', distances, dict(dispersion=False))

# ── Panel B: QBER vs pulse width (at "critical" distance 75 km) ──
print("\nPanel B: QBER vs pulse width (75 km)")
CRIT_DIST = 75
pws_ps = np.array([5, 7, 10, 15, 20, 30, 40, 50])
qber_pulse = sweep('pulse_sigma', pws_ps * 1e-12,
                    dict(fiber_length=CRIT_DIST))

# ── Panel C: QBER vs temperature (at 75 km) ──────────────────────
print("\nPanel C: QBER vs temperature (75 km)")
temps = np.array([0, 10, 20, 25, 30, 40, 50, 60])
qber_temp = sweep('temperature', temps, dict(fiber_length=CRIT_DIST))

# ── Panel D: QBER vs PMD coeff (at 75 km) ────────────────────────
print("\nPanel D: QBER vs PMD coefficient (75 km)")
pmds = np.array([0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3])
qber_pmd = sweep('pm_dispersion', pmds * 1e-12, dict(fiber_length=CRIT_DIST))

# ── Panel E: QBER vs bend radius (at 75 km) ──────────────────────
print("\nPanel E: QBER vs bend radius (75 km)")
bends = np.array([0.002, 0.005, 0.01, 0.015, 0.02, 0.03, 0.05])
qber_bend = sweep('bend_radius', bends, dict(fiber_length=CRIT_DIST))

# ── Save CSV ──────────────────────────────────────────────────────
csv_path = os.path.join(OUT, f'val_system--seed{SEED}.csv')
with open(csv_path, 'w') as f:
    f.write("Panel,Parameter,Value,QBER_fraction\n")
    for d, q in zip(distances, qber_on):
        f.write(f"A,distance_km,{d},{q:.6f}\n")
    for d, q in zip(distances, qber_off):
        f.write(f"A,distance_km_no_disp,{d},{q:.6f}\n")
    for pw, q in zip(pws_ps, qber_pulse):
        f.write(f"B,pulse_width_ps,{pw},{q:.6f}\n")
    for t, q in zip(temps, qber_temp):
        f.write(f"C,temperature_C,{t},{q:.6f}\n")
    for p, q in zip(pmds, qber_pmd):
        f.write(f"D,pmd_coeff_ps_sqrt_km,{p},{q:.6f}\n")
    for r, q in zip(bends, qber_bend):
        f.write(f"E,bend_radius_m,{r},{q:.6f}\n")
print(f"Saved: {csv_path}")

# ── Figure ────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 3, figsize=(18, 12))
axes = axes.flatten()

ax1 = axes[0]
ax1.plot(distances, qber_on * 100, 's-', color='C3', lw=1.5, ms=5,
         label='All impairments')
ax1.plot(distances, qber_off * 100, 'o-', color='C0', lw=1.5, ms=5,
         label='Attenuation only')
ax1.axvline(CRIT_DIST, color='grey', ls=':', lw=1, alpha=0.5,
            label=f'{CRIT_DIST} km (Panels B--E)')
ax1.set_xlabel('Fibre length (km)')
ax1.set_ylabel('QBER (%)')
ax1.set_title('A: QBER vs Distance (30 ps pulse, 0\u20131000 km)', fontweight='bold')
ax1.legend(fontsize=8)
ax1.grid(True, alpha=0.3)
ax1.set_ylim(-5, 105)

ax2 = axes[1]
ax2.plot(pws_ps, qber_pulse * 100, 'D-', color='C2', lw=1.5, ms=5)
ax2.set_xlabel('Pulse $\\sigma$ (ps)')
ax2.set_ylabel('QBER (%)')
ax2.set_title(f'B: QBER vs Pulse Width ({CRIT_DIST} km)', fontweight='bold')
ax2.grid(True, alpha=0.3)
ax2.set_ylim(-5, 105)

ax3 = axes[2]
ax3.plot(temps, qber_temp * 100, '^-', color='C1', lw=1.5, ms=5)
ax3.set_xlabel('Temperature ($^\\circ$C)')
ax3.set_ylabel('QBER (%)')
ax3.set_title(f'C: QBER vs Temperature ({CRIT_DIST} km)', fontweight='bold')
ax3.grid(True, alpha=0.3)
ax3.set_ylim(-5, 105)

ax4 = axes[3]
ax4.plot(pmds, qber_pmd * 100, 'v-', color='C4', lw=1.5, ms=5)
ax4.set_xlabel('PMD coeff. (ps/$\\sqrt{\\mathrm{km}}$)')
ax4.set_ylabel('QBER (%)')
ax4.set_title(f'D: QBER vs PMD Coeff. ({CRIT_DIST} km)', fontweight='bold')
ax4.grid(True, alpha=0.3)
ax4.set_ylim(-5, 105)

ax5 = axes[4]
ax5.plot(bends * 1000, qber_bend * 100, 'o-', color='C5', lw=1.5, ms=5)
ax5.set_xlabel('Bend radius (mm)')
ax5.set_ylabel('QBER (%)')
ax5.set_title(f'E: QBER vs Bend Radius ({CRIT_DIST} km)', fontweight='bold')
ax5.grid(True, alpha=0.3)
ax5.invert_xaxis()
ax5.set_ylim(-5, 105)

axes[5].set_visible(False)

fig.suptitle(
    'System-Level BB84 Demonstration -- Independent Components, Combined Impairments\n'
    r'CWLaser (Henry 1982) $\rightarrow$ MZM Pulse Carving $\rightarrow$ '
    r'Fibre (CD+PMD+Biref+Att) $\rightarrow$ APD (Kasap 2013)'
    rf'   |   {BIT_RATE_MHZ} MHz eff. bit rate (Takesue 2007, Dixon 2008)',
    fontsize=11, fontweight='bold', y=0.98
)

fig_path = os.path.join(OUT, f'val_system--seed{SEED}.png')
fig.savefig(fig_path, dpi=200, bbox_inches='tight')
print(f"Saved: {fig_path}")
plt.close(fig)
