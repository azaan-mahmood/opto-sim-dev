import numpy as np
import matplotlib.pyplot as plt
from src.lasers import NdYAGLaser

def characterize_ndyag():
    print("=======================================================")
    print("  Nd:YAG 4-Level Laser Characterization")
    print("=======================================================")
    
    # Initialize Laser
    laser = NdYAGLaser(noise_std=1e6)
    
    # Run Simulation
    t_span = [0, 500e-6] # 500us to see startup spikes
    power, sol = laser.out_pow(t_span=t_span)
    
    t = sol.t * 1e6 # Convert to microseconds
    N2 = sol.y[0]
    I = sol.y[1]
    
    # Plotting
    fig, axs = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    fig.suptitle("Nd:YAG 4-Level Dynamics: Inversion & Photons", fontsize=14, fontweight='bold')
    
    # 1. Population Inversion (N2)
    # In a 4-level system, N2 is the inversion because N1 ~ 0.
    axs[0].plot(t, N2, color='darkorange', linewidth=2)
    axs[0].set_ylabel("Inversion $N_2$ (atoms/m³)")
    axs[0].set_title("Population Inversion Dynamics (Always Positive)")
    axs[0].grid(True, alpha=0.4)
    
    # 2. Photon Density (I)
    axs[1].plot(t, I, color='mediumseagreen', linewidth=2)
    axs[1].set_ylabel("Photon Density $I$ (photons/m³)")
    axs[1].set_xlabel("Time (µs)")
    axs[1].set_title("Laser Output (Showing Relaxation Oscillations)")
    axs[1].grid(True, alpha=0.4)
    
    plt.tight_layout()
    import os
    if not os.path.exists('analysis'):
        os.makedirs('analysis')
    plt.savefig('analysis/ndyag_dynamics.png', dpi=300)
    print("  -> Saved: analysis/ndyag_dynamics.png")
    
    # Print steady state inversion
    ss_inversion = N2[-1]
    print(f"  Steady State Inversion: {ss_inversion:.2e} atoms/m³")
    print("=======================================================")

if __name__ == "__main__":
    characterize_ndyag()
