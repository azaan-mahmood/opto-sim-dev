import numpy as np
import matplotlib.pyplot as plt
from src.protocols.bb84_high_bitrate import simulate_bb84_high_bitrate

def run_qber_vs_bitrate_experiment():
    print("Running QBER vs Bitrate (Bandwidth) Experiment...")
    
    # Test bandwidths from 10 MHz to 10 GHz
    bandwidths = np.logspace(7, 10, num=10)
    qbers = []
    
    num_bits = 500
    # Use a fixed distance where signal is attenuated but not entirely dead
    fixed_distance = 100 
    
    for bw in bandwidths:
        print(f"Simulating for Bandwidth: {bw/1e6:.1f} MHz...")
        try:
            qber = simulate_bb84_high_bitrate(num_bits=num_bits, fiber_length=fixed_distance, bandwidth=bw, show_pol=False)
            qbers.append(qber)
            print(f"  -> QBER: {qber*100:.2f}%")
        except Exception as e:
            print(f"  -> Error: {e}")
            qbers.append(np.nan)
            
    # Plotting the results
    plt.figure(figsize=(10, 6))
    plt.plot(bandwidths / 1e6, [q * 100 for q in qbers], marker='o', linestyle='-', color='g', linewidth=2)
    plt.xscale('log')
    plt.title('Quantum Bit Error Rate (QBER) vs. Detector Bandwidth')
    plt.xlabel('Detector Bandwidth / Bitrate (MHz)')
    plt.ylabel('QBER (%)')
    plt.grid(True)
    plt.axhline(y=11.0, color='r', linestyle='--', label='Theoretical Abort Threshold (11%)')
    plt.legend()
    
    # Save the figure
    plt.savefig('qber_vs_bitrate.png', dpi=300)
    print("Experiment completed. Plot saved as 'qber_vs_bitrate.png'.")

if __name__ == "__main__":
    run_qber_vs_bitrate_experiment()
