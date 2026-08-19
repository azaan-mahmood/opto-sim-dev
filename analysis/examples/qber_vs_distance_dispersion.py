"""RETIRED to analysis/examples/. Kept because it runs and because the
reasoning in it is worth reading, not because anything depends on it.
Nothing in src/ or run_all.py imports from this directory.

What it was: a QBER-versus-distance sweep with dispersion on, through
bb84_test_dispersion.

What replaced it: analysis/validation/validate_duplinskiy_dispersion.py,
which asks the same question of a chain that replicates a real
experiment. Judged not worth rewiring on its own terms -- the sweep is
a special case of what the protocol validators already cover.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os, sys, time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from src.protocols.bb84_test_dispersion import simulate_bb84_dispersion

OUT = os.path.join(os.path.dirname(__file__))
SEED = 42
NUM_BITS = 300
PULSE_SIGMA = 5e-12
DISTANCES = np.arange(10, 201, 10)

def sweep_qber(dispersion):
    qbers = []
    for L in DISTANCES:
        t0 = time.time()
        qber = simulate_bb84_dispersion(
            NUM_BITS, fiber_length=L, pulse_sigma=PULSE_SIGMA,
            dispersion=dispersion, seed=SEED
        )
        elapsed = time.time() - t0
        print(f"  {L:3d} km  QBER={qber*100:5.1f}%  ({elapsed:.1f}s)")
        qbers.append(qber)
    return qbers

print(f"Sweeping QBER vs distance ({NUM_BITS} bits, {PULSE_SIGMA*1e12:.0f} ps pulse)")
print("--- With dispersion ---")
qber_disp = sweep_qber(True)
print("--- Without dispersion ---")
qber_none = sweep_qber(False)

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(DISTANCES, np.array(qber_disp) * 100, 's-', color='C3', linewidth=1.5,
        markersize=5, label=f'Dispersion ON ({PULSE_SIGMA*1e12:.0f} ps pulse)')
ax.plot(DISTANCES, np.array(qber_none) * 100, 'o-', color='C0', linewidth=1.5,
        markersize=5, label='Dispersion OFF')
ax.set_xlabel('Fiber length (km)')
ax.set_ylabel('Sifted QBER (%)')
ax.set_title(f'BB84 QBER vs Distance — MZM-carved {PULSE_SIGMA*1e12:.0f} ps Gaussian pulse')
ax.legend()
ax.grid(True, alpha=0.3)
ax.set_ylim(-1, 55)

path = os.path.join(OUT, 'qber_vs_distance_dispersion.png')
fig.savefig(path, dpi=150, bbox_inches='tight')
print(f"\nSaved: {path}")
plt.close(fig)
