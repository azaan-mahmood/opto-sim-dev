"""
System-Level Time-Bin BB84 QKD Demonstration
=============================================
Demonstrates that the independently validated components compose correctly
by simulating BB84 QKD under combined impairments (birefringence, CD, PMD,
attenuation) with realistic laser, modulator, interferometer, and SPAD
models — rebuilt on the time-bin (phase-encoded) path per ARCH-1.

Signal chain (per bit):
  CWLaser → MZM (Gaussian carve) → AsymmetricMZI(encoder, phi_A)
    → FiberRealization (birefringence + CD + PMD + attenuation)
    → AsymmetricMZI(decoder, phi_B) → 2x SPAD (PHYS-3 sifting)

Why time-bin encoding?  The old Section-5 demo used polarization encoding
between Ex/Ey measured with a classical linear APD at 10 mW and 1 MHz
bandwidth — fully classical, ~13 orders above the single-photon regime,
and outside the manuscript's own -20..0 dBm validity range (ARCH-1 in
opto-sim-issues-and-fixes.md).  The rebuilt demo is genuine QKD: 0.1
photons/pulse, gated SPADs, basis sifting with non-detections discarded.
A second consequence is physical: time-bin phase encoding is immune to
slow birefringence — the two bins pass through the *same* (quasi-static,
LTI) fibre operator, so the interference phase is preserved regardless of
the channel's random SU(2) rotation.  The impairments therefore enter
through pulse broadening vs. the SPAD gate window (CD/PMD) and through
attenuation vs. dark counts, not through polarization scrambling.

Panel set (rebuilt):
  A: QBER vs fibre distance, all impairments vs attenuation-only
  B: QBER vs pulse width (CD/PMD broadening vs 1 ns gate)
  C: QBER vs decoder visibility V (e_opt = (1-V)/2, Gobby floor)
  D: QBER vs mean photon number mu (dark-count vs multi-photon tradeoff)
  E: QBER vs SPAD dark count rate (detector noise floor at 122 km)

Linearity / quasi-static shortcut
---------------------------------
The fibre and AMZIs are linear and time-invariant, and the impairments
are quasi-static (drawn once per point).  The response of the 8 (phi_A,
phi_B) encoding states is therefore computed exactly once per parameter
point: the early and late basis pulses are propagated separately through
the fibre and decoder, and per bit the gate powers follow

    P_c(phi_A, phi_B, theta) = g0_c + 2 Re[S_c . exp(-j delta)],
    delta = phi_B - phi_A - theta,

where theta ~ N(0, 2*pi*dnu*T_delay) is the laser phase jitter across
the AMZI differential delay (Wiener phase noise, Lorentzian line).  This
is identical to naive per-bit field propagation up to floating-point
rounding (F(E1 + c*E2) = F(E1) + c*F(E2) exactly for a linear channel),
and it makes each parameter point seconds rather than hours.  Detection
statistics (photon Poisson, dark counts, afterpulses, dead time) are
still simulated bit by bit, which is where the Monte Carlo lives.

References
----------
[1] Gobby, C., Yuan, Z. L. & Shields, A. J., "Quantum key distribution
    over 122 km of standard telecom fiber", Appl. Phys. Lett. 84(19),
    3762-3764, 2004.  (AMZI delay 5.8 ns; e_opt = (1-V)/2; V = 0.934
    reproduces the 3.3% short-range floor.)
[2] Bennett, C. H. & Brassard, G., Proc. IEEE 1984 — BB84 sifting.
[3] ID Quantique, ID230 InGaAs SPAD datasheet — 10% QE, 13 us dead time,
    15 Hz DCR, 5% afterpulsing.
[4] Agrawal, G., "Nonlinear Fiber Optics" 5th ed. — CD operator
    H(Omega) = exp(-j*beta2*Omega^2*L/2).
[5] Menyuk & Wai, JLT 12(2), 298-307, 1994 — fibre birefringence model.
[6] Keiser, G., "Optical Fiber Communications" — 10^(-alpha*L/10).

All panels use seeded RNG for reproducibility and report sifted-bit
counts; error bars are the binomial s.d. sqrt(q(1-q)/n_sifted).
"""
import argparse
import os
import random
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from src.lasers.cwlaser import CWLaser
from src.channel import FiberRealization
from src.channel.interferometer import AsymmetricMZI
from src.channel.mzm import MZM
from src.channel.phase_modulator import PhaseModulator
from src.detectors.spad import spad

H = 6.626e-34
C0 = 2.99792458e8
LAM = 1550e-9
ALPHA_DB = 0.182           # SMF-28 at 1550 nm [6]
DELAY = 5.8e-9             # Gobby AMZI differential delay [1]
REP_RATE = 2.5e6           # pulses/s (Gobby)
DT = 1e-12                 # sampling interval
SIGMA_MAX = 50e-12         # widest pulse in panel B
N_SAMPLES = int(np.ceil((2.0 * DELAY + 5.0 * SIGMA_MAX) / DT))
GATE_WIDTH = 1e-9          # SPAD gate (matches validate_gobby)
ETA = 0.10                 # SPAD quantum efficiency [3]
AFTERPULSE = 0.05          # [3]
DEAD_TIME = 13e-6          # [3]
VIS_DEFAULT = 0.934        # e_opt = (1-V)/2 = 3.3% [1]
LINE_WIDTH = 1e6           # CWLaser linewidth (Hz)
THETA_SIGMA = np.sqrt(2.0 * np.pi * LINE_WIDTH * DELAY)


def simulate_point(num_bits, seed, fiber_length=75, pulse_sigma=30e-12,
                   dispersion=True, visibility=VIS_DEFAULT, mu=0.1,
                   dcr=15.0, temperature=25, bend_radius=None,
                   pmd_coeff_ps_sqrt_km=0.1, birefringence=True,
                   attenuation=True, cd=None, pmd=None,
                   gate_width=None, delay=None):
    """Run the time-bin BB84 Monte Carlo at one parameter point.

    `gate_width` (s, default 1 ns) sets both the SPAD gate and the
    integration window.

    `delay` (s, default `DELAY` = 5.8 ns) is the AMZI differential delay,
    i.e. the time-bin separation.  It exists for the OPEN-3 *code-path
    check*: at the Gobby delay, chromatic dispersion cannot produce bin
    crosstalk (CD and the AMZI are both LTI and therefore commute, so CD
    cannot move the constructive/destructive port ratio at all; and even
    ignoring that, spill across 5.8 ns would need z ~ 5,674 km).  Shrinking
    the delay to ~200 ps brings the crosstalk threshold down to z ~ 191 km,
    where CD *does* measurably move the QBER — which demonstrates the CD
    code path is live rather than silently skipped.  That configuration is
    a code-path check, **not** a physical scenario.

    Passing `delay=None` reproduces the module constants exactly
    (bit-identical output); any other value recomputes the sample window
    and the laser phase-jitter scale locally.

    Returns dict with qber (fraction), n_sifted, n_errors, plus the
    per-point configuration.  See module docstring for the linearity
    shortcut that makes each point cheap without changing the physics.

    Impairment toggles (used by `val_system_scenarios.py`, BLOCK-2):
    `dispersion` is the legacy alias enabling both CD and PMD; `cd` /
    `pmd` override it when not None.  `birefringence` and `attenuation`
    toggle the other two chain impairments.
    """
    random.seed(seed)
    np.random.seed(seed)

    photon_energy = H * C0 / LAM

    # Default path uses the module constants verbatim so output stays
    # bit-identical; a custom delay resizes the window and rescales the
    # Wiener phase-jitter sigma (which grows as sqrt(delay)).
    if delay is None:
        delay_s = DELAY
        n_samples = N_SAMPLES
        theta_sigma = THETA_SIGMA
    else:
        delay_s = float(delay)
        n_samples = int(np.ceil((2.0 * delay_s + 5.0 * SIGMA_MAX) / DT))
        theta_sigma = np.sqrt(2.0 * np.pi * LINE_WIDTH * delay_s)

    delay_samples = int(delay_s / DT)

    t = np.arange(n_samples) * DT
    pulse_center = delay_s / 2.0

    # --- Source: CW laser + MZM-carved Gaussian pulse ---
    # The MZM is X-cut: it modulates Ey only.  The source must therefore
    # be Y-polarized, otherwise the unmodulated Ex component passes
    # through as a CW floor that destroys the single-photon calibration
    # (this is why the old classical demo could run at 10 mW without
    # noticing).  Time-bin encoding is polarization-agnostic, so any
    # input SOP is fine physically.
    laser = CWLaser(wavelength=LAM, power_dbm=-10, linewidth=LINE_WIDTH,
                    rin_density=-140, polarization_azimuth=np.pi / 2,
                    polarization_ellipticity=0.0)
    pm = PhaseModulator(crystal_cut='X', modulation='DC')
    Vpi = pm.Vpi
    mzm = MZM(mode='push-pull', bias_voltage=Vpi)
    # Small-signal carve at quadrature: envelope proportional to V_pulse
    V_pulse = Vpi * 0.3 * np.exp(-0.5 * ((t - pulse_center) / pulse_sigma) ** 2)
    E_cw = laser.sample_field(dt=DT, n_samples=n_samples)
    E_carved = mzm.modulate(E_in=E_cw, V=V_pulse)
    # Calibrate pulse energy to mu photons (project convention: calibrate
    # once at the source output)
    energy = float(np.sum(np.abs(E_carved) ** 2) * DT)
    E_carved *= np.sqrt(mu * photon_energy / energy)

    # --- Encoder: 50:50 split, delayed arm ---
    E_early = E_carved / np.sqrt(2.0)
    E_late = np.roll(E_carved / np.sqrt(2.0), delay_samples, axis=0)
    E_late[:delay_samples] = 0.0

    # --- Fibre (quasi-static; one realization per point) ---
    fibre = FiberRealization(
        L_m=fiber_length * 1000, temperature=temperature,
        bend_radius=bend_radius,
        pmd_coeff_ps_sqrt_km=pmd_coeff_ps_sqrt_km,
        attenuation_factor=ALPHA_DB,
        birefringence=birefringence,
        attenuation=attenuation,
        cd=dispersion if cd is None else cd,
        pmd=dispersion if pmd is None else pmd,
        seed=seed)
    F_E = fibre.apply(E_early, dt=DT)
    F_L = fibre.apply(E_late, dt=DT)

    # --- Decoder: finite-visibility combiner (r^2 + s^2 = 1, 2rs = V) ---
    v = float(visibility)
    r = np.sqrt(0.5 + 0.5 * np.sqrt(max(0.0, 1.0 - v * v)))
    s = np.sqrt(0.5 - 0.5 * np.sqrt(max(0.0, 1.0 - v * v)))
    T1_c = s * np.roll(F_E / np.sqrt(2.0), delay_samples, axis=0)
    T1_c[:delay_samples] = 0.0
    T2_c = r * F_L / np.sqrt(2.0)
    T1_d = -T1_c
    T2_d = T2_c

    # Gate-window quadratic-form coefficients (see module docstring)
    gate_width = GATE_WIDTH if gate_width is None else float(gate_width)
    idx = int(pulse_center / DT) + delay_samples
    gate_half = int(gate_width / DT / 2)
    w = slice(idx - gate_half, idx + gate_half + 1)
    n_win = max(1, gate_half * 2 + 1)
    g0_c = float(np.sum(np.abs(T1_c[w]) ** 2 + np.abs(T2_c[w]) ** 2) / n_win)
    S_c = complex(np.sum(np.conj(T1_c[w]) * T2_c[w]) / n_win)
    g0_d = float(np.sum(np.abs(T1_d[w]) ** 2 + np.abs(T2_d[w]) ** 2) / n_win)
    S_d = complex(np.sum(np.conj(T1_d[w]) * T2_d[w]) / n_win)

    # --- SPADs (ID230 specs [3]) ---
    spd_c = spad(wavelength=LAM, quantum_efficiency=ETA, dead_time=DEAD_TIME,
                 dark_count_rate=dcr, afterpulse_prob=AFTERPULSE,
                 gate_width=gate_width)
    spd_d = spad(wavelength=LAM, quantum_efficiency=ETA, dead_time=DEAD_TIME,
                 dark_count_rate=dcr, afterpulse_prob=AFTERPULSE,
                 gate_width=gate_width)

    # Sifting is accumulated inline rather than into per-pulse lists.
    # The lists were O(num_bits) in memory (five of them), which is fine at
    # the 1e6 pulses this script used to run but reaches multiple GB at the
    # ~1e8 pulses `--target-sifted` needs at 100 km -- the original
    # 8th-pass scenario run was killed by exactly that.  Counting inline is
    # O(1) and touches no RNG draw, so results are bit-identical.
    n_sifted = 0
    n_errors = 0

    for i in range(num_bits):
        # --- Alice: encode phase phi_A (X: 0/pi, Y: pi/2 / 3pi/2) ---
        basis_a = np.random.randint(0, 2)
        bit_a = np.random.randint(0, 2)
        phi_A = basis_a * np.pi / 2.0 + bit_a * np.pi

        # --- Bob: basis selection (X: 0, Y: pi/2) ---
        basis_b = np.random.randint(0, 2)
        phi_B = basis_b * np.pi / 2.0

        # Laser phase jitter across the AMZI delay (Wiener, per bit)
        theta = np.random.normal(0.0, theta_sigma)
        delta = phi_B - phi_A - theta
        e_delta = np.exp(-1j * delta)
        P_c = g0_c + 2.0 * np.real(S_c * e_delta)
        P_d = g0_d + 2.0 * np.real(S_d * e_delta)

        t_gate = i / REP_RATE
        click_c = spd_c.detect(P_c, t_gate)
        click_d = spd_d.detect(P_d, t_gate)

        # --- Bob's bit (PHYS-3: no click is a non-detection, not a bit) ---
        if click_c and not click_d:
            bob_bit = 0
        elif click_d and not click_c:
            bob_bit = 1
        elif click_c and click_d:
            bob_bit = np.random.randint(0, 2)
        else:
            bob_bit = -1

        # --- Sifting (PHYS-3): same basis AND a click ---
        if basis_a == basis_b and (click_c or click_d):
            n_sifted += 1
            if bit_a != bob_bit:
                n_errors += 1

    qber = n_errors / n_sifted if n_sifted > 0 else 0.0

    return {'qber': qber, 'n_sifted': n_sifted, 'n_errors': n_errors}


def sweep(pname, vals, fixed, num_bits, seed):
    """Run a parameter sweep; returns (values, qber, sifted)."""
    qs, sfs = [], []
    for v in vals:
        kw = {**fixed, pname: v}
        r = simulate_point(num_bits, seed=seed, **kw)
        qs.append(r['qber'])
        sfs.append(r['n_sifted'])
        print(f"  {pname}={v:<12}  QBER={r['qber']*100:5.1f}%  "
              f"sifted={r['n_sifted']}  errors={r['n_errors']}", flush=True)
    return np.array(vals), np.array(qs), np.array(sfs)


def qber_err(q, n):
    """Binomial s.d. of the QBER estimate (fraction)."""
    n = np.maximum(n, 1)
    return np.sqrt(q * (1.0 - q) / n)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--bits', type=int, default=1000000,
                        help='Pulses per QBER measurement (default 1M)')
    args = parser.parse_args()
    SEED = args.seed
    NUM_BITS = args.bits

    OUT = os.path.join(os.path.dirname(__file__), '..', 'val_system')
    os.makedirs(OUT, exist_ok=True)

    print(f"System-level time-bin BB84 demo — {NUM_BITS} pulses/point, "
          f"seed {SEED}")
    print(f"  theta_sigma = {THETA_SIGMA*1e3:.1f} mrad (laser jitter "
          f"across {DELAY*1e9:.1f} ns delay)")

    # ── Panel A: QBER vs distance ───────────────────────────────────
    print("\nPanel A: QBER vs distance (30 ps pulses, all impairments "
          "vs attenuation only)")
    distances = np.array([0, 4, 10, 20, 40, 65, 80, 100, 122])
    d_on, q_on, s_on = sweep('fiber_length', distances,
                             dict(pulse_sigma=30e-12, dispersion=True),
                             NUM_BITS, SEED)
    d_off, q_off, s_off = sweep('fiber_length', distances,
                                dict(pulse_sigma=30e-12, dispersion=False),
                                NUM_BITS, SEED)

    # ── Panel B: QBER vs pulse width (at 75 km) ─────────────────────
    print("\nPanel B: QBER vs pulse width (75 km)")
    CRIT_DIST = 75
    sigmas_ps = np.array([5, 7, 10, 15, 20, 30, 40, 50])
    b_vals, q_b, s_b = sweep('pulse_sigma', sigmas_ps * 1e-12,
                             dict(fiber_length=CRIT_DIST), NUM_BITS, SEED)

    # ── Panel C: QBER vs decoder visibility (at 75 km) ──────────────
    print("\nPanel C: QBER vs decoder visibility (75 km)")
    c_vals, q_c, s_c = sweep('visibility',
                             np.array([0.90, 0.92, 0.934, 0.95, 0.97, 0.99,
                                       1.00]),
                             dict(fiber_length=CRIT_DIST), NUM_BITS, SEED)

    # ── Panel D: QBER vs mean photon number (at 75 km) ──────────────
    print("\nPanel D: QBER vs mean photon number (75 km)")
    d_vals, q_d, s_d = sweep('mu', np.array([0.02, 0.05, 0.1, 0.2, 0.5, 1.0,
                                             2.0]),
                             dict(fiber_length=CRIT_DIST), NUM_BITS, SEED)

    # ── Panel E: QBER vs dark count rate (at 122 km) ────────────────
    # The ID230 spec (15 Hz) contributes <0.1% QBER even at 122 km, so
    # the sweep runs to 10 kHz to expose the detector noise floor.
    print("\nPanel E: QBER vs dark count rate (122 km, 0-10 kHz)")
    e_vals, q_e, s_e = sweep('dcr', np.array([0.0, 100.0, 300.0, 1000.0,
                                              3000.0, 10000.0]),
                             dict(fiber_length=122), NUM_BITS, SEED)

    # ── Save CSV ────────────────────────────────────────────────────
    csv_path = os.path.join(OUT, f'val_system--seed{SEED}.csv')
    with open(csv_path, 'w') as f:
        f.write(f"# System-level time-bin BB84 demo (ARCH-1 rebuild), "
                f"seed {SEED}, {NUM_BITS} pulses/point\n")
        f.write("# Chain: CWLaser(1 MHz) -> MZM carve -> encoder AMZI "
                "-> FiberRealization -> decoder AMZI -> 2x SPAD\n")
        f.write("# Linearity shortcut: channel response computed once per "
                "point; detection Monte Carlo per bit.\n")
        f.write("# Defaults: V=0.934, mu=0.1, DCR=15 Hz, eta=0.10, gate "
                "1 ns, alpha=0.182 dB/km, delay=5.8 ns, rep=2.5 MHz\n")
        f.write("# Floor: visibility 3.3% + phase jitter ~0.8% + afterpulse "
                "~1.5% + double-click ~0.5% (~6.4% at 0 km)\n")
        f.write("# Panel E sweeps DCR to 10 kHz (ID230 spec 15 Hz is "
                "<0.1% even at 122 km)\n")
        f.write("Panel,Parameter,Value,QBER_fraction,Sifted_bits\n")
        for x, q, s in zip(distances, q_on, s_on):
            f.write(f"A_allimp,distance_km,{x},{q:.6f},{s}\n")
        for x, q, s in zip(distances, q_off, s_off):
            f.write(f"A_attonly,distance_km,{x},{q:.6f},{s}\n")
        for x, q, s in zip(sigmas_ps, q_b, s_b):
            f.write(f"B,pulse_sigma_ps,{x},{q:.6f},{s}\n")
        for x, q, s in zip(c_vals, q_c, s_c):
            f.write(f"C,visibility,{x},{q:.6f},{s}\n")
        for x, q, s in zip(d_vals, q_d, s_d):
            f.write(f"D,mu_photons,{x},{q:.6f},{s}\n")
        for x, q, s in zip(e_vals, q_e, s_e):
            f.write(f"E,dcr_hz,{x},{q:.6f},{s}\n")
    print(f"Saved: {csv_path}")

    # ── Figure ──────────────────────────────────────────────────────
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    axes = axes.flatten()

    ax1 = axes[0]
    ax1.errorbar(distances, q_on * 100, yerr=qber_err(q_on, s_on) * 100,
                 fmt='s-', color='C3', lw=1.5, ms=5, capsize=3,
                 label='All impairments')
    ax1.errorbar(distances, q_off * 100, yerr=qber_err(q_off, s_off) * 100,
                 fmt='o-', color='C0', lw=1.5, ms=5, capsize=3,
                 label='Attenuation only')
    ax1.axvline(CRIT_DIST, color='grey', ls=':', lw=1, alpha=0.5,
                label=f'{CRIT_DIST} km (Panels B--D)')
    ax1.set_xlabel('Fibre length (km)')
    ax1.set_ylabel('QBER (%)')
    ax1.set_title('A: QBER vs Distance (30 ps pulse)', fontweight='bold')
    ax1.legend(fontsize=8)
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(-1, 30)

    ax2 = axes[1]
    ax2.errorbar(sigmas_ps, q_b * 100, yerr=qber_err(q_b, s_b) * 100,
                 fmt='D-', color='C2', lw=1.5, ms=5, capsize=3)
    ax2.set_xlabel('Pulse $\\sigma$ (ps)')
    ax2.set_ylabel('QBER (%)')
    ax2.set_title(f'B: QBER vs Pulse Width ({CRIT_DIST} km)',
                  fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(-1, 20)

    ax3 = axes[2]
    ax3.errorbar(c_vals, q_c * 100, yerr=qber_err(q_c, s_c) * 100,
                 fmt='^-', color='C1', lw=1.5, ms=5, capsize=3)
    ax3.axvline(VIS_DEFAULT, color='grey', ls=':', lw=1, alpha=0.5,
                label='Gobby V = 0.934')
    ax3.set_xlabel('Decoder visibility $V$')
    ax3.set_ylabel('QBER (%)')
    ax3.set_title(f'C: QBER vs Visibility ({CRIT_DIST} km)',
                  fontweight='bold')
    ax3.legend(fontsize=8)
    ax3.grid(True, alpha=0.3)
    ax3.set_ylim(-1, 20)

    ax4 = axes[3]
    ax4.errorbar(d_vals, q_d * 100, yerr=qber_err(q_d, s_d) * 100,
                 fmt='v-', color='C4', lw=1.5, ms=5, capsize=3)
    ax4.set_xlabel('Mean photons per pulse $\\mu$')
    ax4.set_ylabel('QBER (%)')
    ax4.set_title(f'D: QBER vs $\\mu$ ({CRIT_DIST} km)', fontweight='bold')
    ax4.grid(True, alpha=0.3)
    ax4.set_ylim(-1, 30)

    ax5 = axes[4]
    ax5.errorbar(e_vals, q_e * 100, yerr=qber_err(q_e, s_e) * 100,
                 fmt='o-', color='C5', lw=1.5, ms=5, capsize=3)
    ax5.set_xlabel('Dark count rate (Hz)')
    ax5.set_ylabel('QBER (%)')
    ax5.set_title('E: QBER vs DCR (122 km)', fontweight='bold')
    ax5.grid(True, alpha=0.3)
    ax5.set_ylim(-1, 60)

    ax6 = axes[5]
    ax6.axis('off')
    notes = (
        "Chain: CWLaser (1 MHz) $\\rightarrow$ MZM carve\n"
        "$\\rightarrow$ encoder AMZI $\\rightarrow$ FiberRealization\n"
        "$\\rightarrow$ decoder AMZI $\\rightarrow$ 2$\\times$ SPAD\n\n"
        "Defaults: $\\mu$ = 0.1, $\\eta$ = 10 %, gate 1 ns,\n"
        "DCR 15 Hz, afterpulse 5 %, dead time 13 $\\mu$s,\n"
        "$\\alpha$ = 0.182 dB/km, delay 5.8 ns, 2.5 MHz\n\n"
        "Time-bin encoding is immune to slow birefringence:\n"
        "both bins traverse the same quasi-static fibre operator,\n"
        "so the interference phase survives CD, PMD and the\n"
        "channel's random SU(2) rotation (impairments enter the\n"
        "chain but affect rate, not QBER).\n\n"
        "QBER floor: visibility (1-V)/2 = 3.3 % + laser phase\n"
        "jitter across the delay ~0.8 % + afterpulse ~1.5 %\n"
        "+ double-click ~0.5 %  ->  ~6.4 % at 0 km.\n\n"
        "Error bars: binomial $\\sqrt{q(1-q)/n_\\mathrm{sifted}}$;\n"
        "long-range points are sample-limited at 1M pulses\n"
        "(the 10M-pulse Gobby table is the precise anchor)."
    )
    ax6.text(0.02, 0.98, notes, transform=ax6.transAxes,
             fontsize=9, va='top', family='monospace')

    fig.suptitle(
        'System-Level Time-Bin BB84 QKD -- Independent Components, '
        'Combined Impairments\n'
        r'CWLaser $\rightarrow$ MZM Carve $\rightarrow$ Encoder AMZI '
        r'$\rightarrow$ Fibre (Biref+CD+PMD+Att) '
        r'$\rightarrow$ Decoder AMZI $\rightarrow$ 2$\times$ SPAD'
        rf'   |   seed {SEED}, {NUM_BITS} pulses/point',
        fontsize=11, fontweight='bold', y=0.98
    )

    fig_path = os.path.join(OUT, f'val_system--seed{SEED}.png')
    fig.savefig(fig_path, dpi=200, bbox_inches='tight')
    print(f"Saved: {fig_path}")
    plt.close(fig)


if __name__ == "__main__":
    main()
