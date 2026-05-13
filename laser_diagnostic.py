"""
Diagnostic script to evaluate the physics of the SolidStateLaser rate equations.
Prints key values at each step to identify inconsistencies.
"""
import numpy as np
from scipy.integrate import solve_ivp
from src.lasers.sslaser import SolidStateLaser

laser = SolidStateLaser(
    wavelength=1550e-9, polarization_azimuth=np.pi,
    polarization_ellipticity=np.pi / 2, power_dbm=0, frequency=1e6
)

print("=" * 55)
print("  LASER PARAMETER AUDIT")
print("=" * 55)
print(f"  Total atoms          N0  = {laser.N0:.3e}")
print(f"  Lower level tau      t1  = {laser.tau1:.3e} s")
print(f"  Upper level tau      t2  = {laser.tau2:.3e} s")
print(f"  Photon lifetime      tc  = {laser.tau_c:.3e} s")
print(f"  Gain damping         a   = {laser.alpha}")
print(f"  Pump rate            Rp  = {laser.Rp:.3e} s^-1")
print(f"  Sigma (abs/em)       sig = {laser.sigma12:.3e} m^2")
print(f"  Speed of light       c   = {laser.c:.3e} m/s")
print(f"  Input power (dbm)    Pdbm= {laser.power_dbm} dBm")
print(f"  Input power (mW)     Pmw = {laser.power_mw:.3e} mW")
print(f"  Frequency used       f   = {laser.frequency:.3e} Hz")
print(f"  Initial photon den.  I0  = {laser.I_0:.3e}")
print()

# Check initial conditions
N1_0 = 2e23
N2_0 = 0
print("=" * 55)
print("  INITIAL CONDITIONS AUDIT")
print("=" * 55)
print(f"  N1_0 = {N1_0:.3e}  (lower level initial pop)")
print(f"  N2_0 = {N2_0:.3e}  (upper level initial pop)")
print(f"  I_0  = {laser.I_0:.3e}  (initial photon density)")
print()
print(f"  [!] N0 (total atoms)      = {laser.N0:.3e}")
print(f"  [!] N1_0 >> N0?           = {N1_0 > laser.N0}  (N1_0={N1_0:.2e} vs N0={laser.N0:.2e})")
print(f"  [!] N1_0 + N2_0 vs N0    : {N1_0+N2_0:.2e} vs {laser.N0:.2e} --> ratio = {(N1_0+N2_0)/laser.N0:.1f}x")
print()

# Check gain vs loss at t=0
gain_at_t0 = laser.c * laser.sigma12 * (N2_0 - laser.g2byg1 * N1_0)
loss_at_t0 = laser.alpha / laser.tau_c
print("=" * 55)
print("  GAIN/LOSS AUDIT at t=0")
print("=" * 55)
print(f"  Gain coeff  g(t=0)  = c*sig*(N2-N1) = {gain_at_t0:.3e} s^-1")
print(f"  Loss coeff  l       = a/tc          = {loss_at_t0:.3e} s^-1")
print(f"  Net gain at t=0     = {gain_at_t0 - loss_at_t0:.3e} s^-1")
print(f"  [!] Gain < 0 at t=0? {gain_at_t0 < 0} -- means photons are ABSORBED, not amplified!")
print()

# Run ODE and show snapshots
N1_0_corrected = laser.N0 * 0.999  # physically consistent: almost all atoms in ground state
N2_0_corrected = laser.N0 * 0.001
print("=" * 55)
print("  ODE SOLUTION SNAPSHOTS (Original N1_0=2e23)")
print("=" * 55)
y0_orig = [N1_0, N2_0, laser.I_0]
sol = solve_ivp(laser.rate, [0, 1e-6], y0_orig,
                t_eval=np.linspace(0, 1e-6, 1000), method='BDF', rtol=1e-6, atol=1e-9)
for idx in [0, 100, 250, 500, 999]:
    print(f"  t={sol.t[idx]*1e9:6.1f} ns | N1={sol.y[0][idx]:.3e} | N2={sol.y[1][idx]:.3e} | I={sol.y[2][idx]:.3e}")

print()
print("=" * 55)
print("  ODE SOLUTION SNAPSHOTS (Corrected N1_0=N0*0.999)")
print("=" * 55)
y0_fixed = [N1_0_corrected, N2_0_corrected, laser.I_0]
sol_fixed = solve_ivp(laser.rate, [0, 1e-6], y0_fixed,
                      t_eval=np.linspace(0, 1e-6, 1000), method='BDF', rtol=1e-6, atol=1e-9)
for idx in [0, 100, 250, 500, 999]:
    print(f"  t={sol_fixed.t[idx]*1e9:6.1f} ns | N1={sol_fixed.y[0][idx]:.3e} | N2={sol_fixed.y[1][idx]:.3e} | I={sol_fixed.y[2][idx]:.3e}")

print()
print("  Diagnostic complete.")
