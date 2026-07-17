import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os, sys, argparse

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from src.channel.fiber import apply_birefringence

OUT = os.path.join(os.path.dirname(__file__), '..', 'val_birefringence')
os.makedirs(OUT, exist_ok=True)

parser = argparse.ArgumentParser()
parser.add_argument('--seed', type=int, default=42)
args = parser.parse_args()
SEED = args.seed
np.random.seed(SEED)

WAVELENGTH = 1550e-9

def test_power_conservation():
    E = np.random.randn(1000, 2) + 1j * np.random.randn(1000, 2)
    P_in = np.mean(np.abs(E)**2)
    for L_m in [1, 10, 100, 1000]:
        E_out = apply_birefringence(E.copy(), L_m, wavelength=WAVELENGTH)
        P_out = np.mean(np.abs(E_out)**2)
        err = abs(P_out - P_in) / P_in
        assert err < 1e-12, f"Power not conserved at L={L_m}m: err={err:.2e}"
    print("  [PASS] Power conservation (unitary Jones matrix)")

def test_phase_shift_scales_with_length():
    E = np.ones((100, 2), dtype=complex)
    E[:, 1] = 0
    # Dense sampling so np.unwrap can resolve wraparound
    Ls = np.arange(1, 101)
    phases = []
    for L_m in Ls:
        E_out = apply_birefringence(E.copy(), L_m, wavelength=WAVELENGTH)
        phases.append(np.angle(E_out[0, 0]))
    phases = np.unwrap(np.array(phases))
    dphi_10_1 = phases[9] - phases[0]
    dphi_100_10 = phases[99] - phases[9]
    assert dphi_10_1 != 0 and dphi_100_10 != 0, "Phase shift must be non-zero"
    ratio = dphi_100_10 / dphi_10_1
    assert 8.0 < abs(ratio) < 12.0, f"Expected |ratio| ≈ 10, got {ratio:.2f}"
    print("  [PASS] Phase shift scales linearly with length")

def test_temperature_dependence():
    E = np.ones((100, 2), dtype=complex)
    phases = []
    for T in [0, 25, 50]:
        E_out = apply_birefringence(E.copy(), 1000, wavelength=WAVELENGTH, temperature=T)
        phases.append(np.angle(E_out[0, 0]))
    assert not np.allclose(phases[0], phases[1]), "0 and 25 should differ"
    assert not np.allclose(phases[1], phases[2]), "25 and 50 should differ"
    print("  [PASS] Temperature sensitivity detected")

def test_wavelength_dependence():
    E = np.ones((100, 2), dtype=complex)
    phases = []
    for lam in [1310e-9, 1550e-9]:
        E_out = apply_birefringence(E.copy(), 1000, wavelength=lam)
        phases.append(np.angle(E_out[0, 0]))
    assert not np.allclose(phases[0], phases[1]), "1310 and 1550 nm should differ"
    print("  [PASS] Wavelength dependence detected")

print("Birefringence validation via apply_birefringence")
test_power_conservation()
test_phase_shift_scales_with_length()
test_temperature_dependence()
test_wavelength_dependence()

L_m_range = np.logspace(0, 4, 50)
E_in = np.ones((1, 2), dtype=complex)
phases = []
for L_m in L_m_range:
    E_out = apply_birefringence(E_in.copy(), L_m, wavelength=WAVELENGTH)
    phases.append(np.angle(E_out[0, 0]))
phases = np.unwrap(np.array(phases))  # now continuous

fig, ax1 = plt.subplots(figsize=(8, 4))
ax1.semilogx(L_m_range, phases, '.-', c='C0', ms=4)
ax1.set(xlabel='Fibre length (m)', ylabel='Phase shift of Ex (rad)',
        title=f'Birefringence: phase vs length (λ = {WAVELENGTH*1e9:.0f} nm)')
ax1.grid(True, alpha=0.3)

fig.tight_layout()
fname = f'val_birefringence--seed{SEED}.png'
fig.savefig(os.path.join(OUT, fname), dpi=150)
print(f"\nSaved: {fname}")

csv_name = f'val_birefringence--seed{SEED}.csv'
np.savetxt(os.path.join(OUT, csv_name),
           np.column_stack([L_m_range, phases]),
           delimiter=',', header='length_m,phase_shift_rad', comments='')
print(f"Saved: {csv_name}")
plt.close(fig)
