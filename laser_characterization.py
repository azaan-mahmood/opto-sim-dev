import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from src.lasers.sslaser import SolidStateLaser
import os

# Ensure analysis directory exists
os.makedirs('analysis', exist_ok=True)


def get_raw_sol(laser):
    """Re-run rate equations to get full solution object (N1, N2, I over time)."""
    # Enforce physical initial conditions as established in sslaser.py
    N1_0 = laser.N0
    N2_0 = 0
    y0 = [N1_0, N2_0, laser.I_0]
    t_span = [0, 20e-3] # 20ms window to capture atomic lifetimes
    t_eval = np.linspace(0, 20e-3, 2000)
    sol = solve_ivp(
        laser.rate, t_span, y0,
        t_eval=t_eval, method='BDF', rtol=1e-6, atol=1e-9
    )
    return sol


# ─────────────────────────────────────────────────────────
# Graph 1: Population Dynamics - N1(t), N2(t), I(t)
# ─────────────────────────────────────────────────────────
def plot_population_dynamics():
    print("Graph 1: Plotting Population Dynamics...")
    laser = SolidStateLaser(
        wavelength=1550e-9, polarization_azimuth=np.pi,
        polarization_ellipticity=np.pi / 2, power_dbm=0, frequency=1e6,
        noise_std=1e6 # Real Langevin noise integrated into the ODE
    )
    sol = get_raw_sol(laser)
    t = sol.t * 1e3  # Convert to milliseconds

    fig, axs = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
    fig.suptitle("Laser Rate Equation: Population Dynamics", fontsize=14, fontweight='bold')

    axs[0].plot(t, sol.y[0], color='steelblue', linewidth=2)
    axs[0].set_ylabel("$N_1(t)$ (atoms)")
    axs[0].set_title("Lower Energy Level Population $N_1$")
    axs[0].grid(True, alpha=0.4)

    axs[1].plot(t, sol.y[1], color='darkorange', linewidth=2)
    axs[1].set_ylabel("$N_2(t)$ (atoms)")
    axs[1].set_title("Upper Energy Level Population $N_2$")
    axs[1].grid(True, alpha=0.4)

    axs[2].plot(t, np.abs(sol.y[2]), color='mediumseagreen', linewidth=2)
    axs[2].set_ylabel("$I(t)$ (photons/m³)")
    axs[2].set_xlabel("Time (ms)")
    axs[2].set_title("Intracavity Photon Density $I$ (with Physical Langevin Noise)")
    axs[2].grid(True, alpha=0.4)

    plt.tight_layout()
    plt.savefig('analysis/graph1_population_dynamics.png', dpi=300)
    plt.close()
    print("  -> Saved: analysis/graph1_population_dynamics.png")

def plot_combined_populations():
    print("Graph 6: Plotting Combined Population Dynamics...")
    laser = SolidStateLaser(
        wavelength=1550e-9, polarization_azimuth=np.pi,
        polarization_ellipticity=np.pi / 2, power_dbm=0, frequency=1e6
    )
    sol = get_raw_sol(laser)
    t = sol.t * 1e3 # ms
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(t, sol.y[0], label='Lower Level ($N_1$)', color='steelblue', alpha=0.8)
    ax.plot(t, sol.y[1], label='Upper Level ($N_2$)', color='darkorange', alpha=0.8)
    
    # Shade the inter-region
    ax.fill_between(t, sol.y[0], sol.y[1], color='gray', alpha=0.2, label='Inversion Region')
    
    # Highlight crossing point (Transparency Threshold)
    crossing_idx = np.where(sol.y[1] > sol.y[0])[0]
    if len(crossing_idx) > 0:
        ax.axvline(t[crossing_idx[0]], color='red', linestyle='--', label='Transparency Point')
        
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Population (atoms)")
    ax.set_title("Combined Population Dynamics: Inversion Analysis", fontsize=14)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('analysis/graph6_combined_populations.png', dpi=300)
    plt.close()
    print("  -> Saved: analysis/graph6_combined_populations.png")


# ─────────────────────────────────────────────────────────
# Graph 2: Population Inversion ΔN(t)
# ─────────────────────────────────────────────────────────
def plot_population_inversion():
    print("Graph 2: Plotting Population Inversion...")
    laser = SolidStateLaser(
        wavelength=1550e-9, polarization_azimuth=np.pi,
        polarization_ellipticity=np.pi / 2, power_dbm=0, frequency=1e6,
        noise_std=1e6
    )
    sol = get_raw_sol(laser)
    t = sol.t * 1e3 # ms
    g2g1 = laser.g2byg1
    delta_N = sol.y[1] - g2g1 * sol.y[0]  # N2 - (g2/g1)*N1

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(t, delta_N, color='crimson', linewidth=2, label='$\\Delta N = N_2 - (g_2/g_1) N_1$')
    ax.axhline(0, color='black', linestyle='--', linewidth=1, label='Threshold (ΔN = 0)')
    ax.fill_between(t, delta_N, 0, where=(delta_N > 0), alpha=0.2, color='green', label='Population Inversion (Gain)')
    ax.fill_between(t, delta_N, 0, where=(delta_N < 0), alpha=0.2, color='red', label='Absorption (Loss)')
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("$\\Delta N$ (atoms)")
    ax.set_title("Population Inversion $\\Delta N(t)$ — Lasing Threshold Analysis", fontsize=13)
    ax.legend()
    ax.grid(True, alpha=0.4)

    plt.tight_layout()
    plt.savefig('analysis/graph2_population_inversion.png', dpi=300)
    plt.close()
    print("  -> Saved: analysis/graph2_population_inversion.png")


# ─────────────────────────────────────────────────────────
# Graph 3: L-I Curve — Output Power vs Pump Rate
# ─────────────────────────────────────────────────────────
def plot_li_curve():
    print("Graph 3: Plotting L-I Curve (Output Power vs Pump Rate)...")
    # Use a base laser to get defaults, then sweep Rp manually
    base = SolidStateLaser(
        wavelength=1550e-9, polarization_azimuth=np.pi,
        polarization_ellipticity=np.pi / 2, power_dbm=0, frequency=1e6,
        noise_std=1e6
    )
    tau2 = base.tau2
    Rp_threshold = base.g2byg1 / tau2

    Rp_values = np.linspace(0.1 * Rp_threshold, 5 * Rp_threshold, 40)
    power_outputs = []

    for Rp_val in Rp_values:
        # Temporarily patch Rp and re-run ODE
        laser = SolidStateLaser(
            wavelength=1550e-9, polarization_azimuth=np.pi,
            polarization_ellipticity=np.pi / 2, power_dbm=0, frequency=1e6,
            noise_std=1e6
        )
        laser.Rp = Rp_val
        laser.power_out, _ = laser.out_pow()
        power_outputs.append(laser.power_out * 1e3)  # mW

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(Rp_values / Rp_threshold, power_outputs, color='royalblue', linewidth=2, marker='o', markersize=4)
    ax.axvline(x=1.0, color='red', linestyle='--', label='Threshold ($R_p = R_{p,th}$)')
    ax.set_xlabel("Normalized Pump Rate $R_p / R_{p,th}$")
    ax.set_ylabel("Output Power (mW)")
    ax.set_title("L-I Curve: Output Power vs. Normalized Pump Rate", fontsize=13)
    ax.legend()
    ax.grid(True, alpha=0.4)

    plt.tight_layout()
    plt.savefig('analysis/graph3_li_curve.png', dpi=300)
    plt.close()
    print("  -> Saved: analysis/graph3_li_curve.png")


# ─────────────────────────────────────────────────────────
# Graph 4: Polarization Ellipse / Lissajous
# ─────────────────────────────────────────────────────────
def plot_polarization_ellipse():
    print("Graph 4: Plotting Polarization Ellipse...")
    configs = [
        (0, 0, "Linear Horizontal (φ=0, χ=0)"),
        (np.pi/4, np.pi/4, "Linear 45° (φ=π/4, χ=π/4)"),
        (np.pi/2, 0, "Linear Vertical (φ=π/2, χ=0)"),
        (np.pi/4, -np.pi/4, "Right Elliptical (φ=π/4, χ=-π/4)"),
        (np.pi/2, np.pi/2, "Left Circular (φ=π/2, χ=π/2)"),
    ]
    colors = ['steelblue', 'darkorange', 'mediumseagreen', 'crimson', 'purple']

    fig, axs = plt.subplots(1, len(configs), figsize=(16, 4))
    fig.suptitle("Polarization States (Ex vs Ey Lissajous)", fontsize=13, fontweight='bold')

    for ax, (phi, chi, label), color in zip(axs, configs, colors):
        laser = SolidStateLaser(
            wavelength=1550e-9, polarization_azimuth=phi,
            polarization_ellipticity=chi, power_dbm=0, frequency=1e6,
            noise_std=1e6
        )
        E = laser.get_electric_field(normalize=False, over_period=True)
        Ex_real = np.real(E[:, 0])
        Ey_real = np.real(E[:, 1])

        ax.plot(Ex_real, Ey_real, color=color, linewidth=2)
        ax.set_aspect('equal')
        ax.set_title(label, fontsize=8)
        ax.set_xlabel("$Re(E_x)$")
        ax.set_ylabel("$Re(E_y)$")
        ax.grid(True, alpha=0.4)

    plt.tight_layout()
    plt.savefig('analysis/graph4_polarization_ellipse.png', dpi=300)
    plt.close()
    print("  -> Saved: analysis/graph4_polarization_ellipse.png")


# ─────────────────────────────────────────────────────────
# Graph 5: Output Power vs Input Power (dBm sweep)
# ─────────────────────────────────────────────────────────
def plot_power_transfer():
    print("Graph 5: Plotting Output Power vs Input Power...")
    dbm_values = np.linspace(-20, 10, 30)
    power_out_mw = []

    for dbm in dbm_values:
        laser = SolidStateLaser(
            wavelength=1550e-9, polarization_azimuth=np.pi,
            polarization_ellipticity=np.pi / 2, power_dbm=dbm, frequency=1e6,
            noise_std=1e6
        )
        power_out_mw.append(laser.power_out * 1e3)  # convert to mW

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(dbm_values, power_out_mw, color='darkorange', linewidth=2, marker='o', markersize=4)
    ax.set_xlabel("Input Power (dBm)")
    ax.set_ylabel("Output Power (mW)")
    ax.set_title("Power Transfer: Output Power vs. Input Power (dBm Sweep)", fontsize=13)
    ax.grid(True, alpha=0.4)

    plt.tight_layout()
    plt.savefig('analysis/graph5_power_transfer.png', dpi=300)
    plt.close()
    print("  -> Saved: analysis/graph5_power_transfer.png")

# ─────────────────────────────────────────────────────────
# Graph 6: Combined Population Dynamics (Inversion)
# ─────────────────────────────────────────────────────────
def plot_combined_populations():
    print("Graph 6: Plotting Combined Population Dynamics...")
    laser = SolidStateLaser(
        wavelength=1550e-9, polarization_azimuth=np.pi,
        polarization_ellipticity=np.pi / 2, power_dbm=0, frequency=1e6,
        noise_std=1e6 # Physically integrated Langevin noise
    )
    sol = get_raw_sol(laser)
    t = sol.t * 1e3 # ms
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(t, sol.y[0], label='Lower Level ($N_1$)', color='steelblue', alpha=0.9, linewidth=1.5)
    ax.plot(t, sol.y[1], label='Upper Level ($N_2$)', color='darkorange', alpha=0.9, linewidth=1.5)
    
    # Shade the inter-region to show inversion dynamics
    ax.fill_between(t, sol.y[0], sol.y[1], where=(sol.y[1] > sol.y[0]), 
                    color='green', alpha=0.15, label='Population Inversion')
    ax.fill_between(t, sol.y[0], sol.y[1], where=(sol.y[1] <= sol.y[0]), 
                    color='red', alpha=0.05, label='Absorption Regime')
    
    # Highlight crossing point (Transparency Threshold)
    crossing_idx = np.where(sol.y[1] > sol.y[0])[0]
    if len(crossing_idx) > 0:
        ax.axvline(t[crossing_idx[0]], color='red', linestyle='--', linewidth=1, label='Transparency Threshold')
        
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Population (atoms/m³)")
    ax.set_title("Combined Population Dynamics: Ground State vs. Excited State", fontsize=14, fontweight='bold')
    ax.legend(loc='best', frameon=True, shadow=True)
    ax.grid(True, linestyle=':', alpha=0.6)
    
    plt.tight_layout()
    plt.savefig('analysis/graph6_combined_populations.png', dpi=300)
    plt.close()
    print("  -> Saved: analysis/graph6_combined_populations.png")


if __name__ == "__main__":
    print("=" * 55)
    print("  Laser Characterization - SolidStateLaser (Er:Yb)")
    print("=" * 55)
    plot_population_dynamics()
    plot_population_inversion()
    plot_li_curve()
    plot_polarization_ellipse()
    plot_power_transfer()
    plot_combined_populations()
    print("=" * 55)
    print("  All graphs saved to: analysis/")
    print("=" * 55)
