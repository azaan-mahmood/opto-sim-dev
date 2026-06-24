import numpy as np
import matplotlib.pyplot as plt
import os
from src.deprecated.sslaser import SolidStateLaser
from src.channel.phase_modulator import PhaseModulator
from src.channel import optics

def compare_dc_rf():
    print("Running DC vs RF Modulation comparison...")
    
    # 1. Generate Optical Pulse
    laser = SolidStateLaser(
        wavelength=1550e-9,
        polarization_azimuth=np.pi,
        polarization_ellipticity=np.pi / 2,
        power_dbm=0, 
        frequency=1e6
    )
    E_in = laser.get_electric_field(normalize=False, over_period=True)
    t = np.linspace(0, 2 * np.pi / laser.frequency, 1000)
    
    # Align polarization to 45 degrees
    E_in = optics.polarizer(E_in, '45')
    
    # Calculate initial power
    power_in = np.abs(E_in[:, 0])**2 + np.abs(E_in[:, 1])**2
    peak_idx = np.argmax(power_in)
    t_peak = t[peak_idx]
    
    # 2. Setup Modulators
    pm_dc = PhaseModulator(crystal_cut='X', modulation="DC")
    pm_rf = PhaseModulator(crystal_cut='X', modulation="RF")
    Vpi = pm_dc.Vpi
    
    # 3. DC Modulation (Constant Vpi)
    E_dc = pm_dc.modulate(E_in, Vpi)
    
    # 4. RF Modulation (Gaussian RF pulse)
    # The RF pulse is centered on the optical peak, but is narrower
    rf_width = (t[-1] / 10) / 2 # Narrow electrical pulse
    V_rf = Vpi * np.exp(-((t - t_peak)**2) / (2 * rf_width**2))
    E_rf = pm_rf.modulate(E_in, V_rf)
    
    # 5. Interference Measurement
    # A pi phase shift on the Y axis flips 45 deg polarization to -45 deg.
    # Passing through a -45 deg polarizer reveals how much of the pulse achieved the full phase shift.
    E_dc_out = optics.polarizer(E_dc, '-45')
    E_rf_out = optics.polarizer(E_rf, '-45')
    
    power_dc_out = np.abs(E_dc_out[:, 0])**2 + np.abs(E_dc_out[:, 1])**2
    power_rf_out = np.abs(E_rf_out[:, 0])**2 + np.abs(E_rf_out[:, 1])**2
    
    # Plotting
    fig, axs = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
    
    # Plot 1: Optical Input Power
    axs[0].plot(t, power_in, 'k-', linewidth=2, label="Input Optical Power (45°)")
    axs[0].set_ylabel("Power")
    axs[0].set_title("Input Optical Pulse (from SolidStateLaser)")
    axs[0].legend()
    axs[0].grid(True)
    
    # Plot 2: Applied Voltages
    axs[1].plot(t, np.ones_like(t) * Vpi, 'b--', linewidth=2, label="DC Voltage (V = Vπ)")
    axs[1].plot(t, V_rf, 'r-', linewidth=2, label="RF Voltage (Gaussian)")
    axs[1].set_ylabel("Voltage (V)")
    axs[1].set_title("Applied Modulation Voltage")
    axs[1].legend()
    axs[1].grid(True)
    
    # Plot 3: Output Power after -45 polarizer
    axs[2].plot(t, power_dc_out, 'b--', linewidth=2, label="DC Output Power (-45°)")
    axs[2].plot(t, power_rf_out, 'r-', linewidth=2, label="RF Output Power (-45°)")
    axs[2].fill_between(t, power_dc_out, power_rf_out, color='gray', alpha=0.3, label="Information Loss / Distortion")
    axs[2].set_xlabel("Time (s)")
    axs[2].set_ylabel("Power")
    axs[2].set_title("Modulated Output (Interfered via -45° Polarizer)")
    axs[2].legend()
    axs[2].grid(True)
    
    plt.tight_layout()
    
    # Ensure analysis dir exists
    if not os.path.exists('analysis'):
        os.makedirs('analysis')
        
    save_path = 'analysis/dc_vs_rf_modulation.png'
    plt.savefig(save_path, dpi=300)
    print(f"Comparison graph saved as '{save_path}'")

if __name__ == "__main__":
    compare_dc_rf()
