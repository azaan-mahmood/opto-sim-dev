import numpy as np
import matplotlib.pyplot as plt
from src.protocols.examples.bb84_ideal import simulate_bb84

def run_qber_vs_distance_experiment():
    print("Running QBER vs Distance Experiment...")
    distances = np.arange(10, 201, 20)  # From 10 km to 200 km
    qbers = []
    
    num_bits = 500  # Number of bits per simulation point
    
    for dist in distances:
        print(f"Simulating for distance: {dist} km...")
        try:
            # Note: Set dispersion=True in cable() within bb84_ideal.py to see real effects, 
            # currently it defaults to dispersion=False in bb84_ideal.py
            results = simulate_bb84(num_bits=num_bits, fiber_length=dist, show_pol=False)
            qber = results[7] # QBER is the 8th item returned
            qbers.append(qber)
            print(f"  -> QBER: {qber*100:.2f}%")
        except Exception as e:
            print(f"  -> Error at {dist} km: {e}")
            qbers.append(np.nan)
            
    # Plotting the results
    plt.figure(figsize=(10, 6))
    plt.plot(distances, [q * 100 for q in qbers], marker='o', linestyle='-', color='b', linewidth=2)
    plt.title('Quantum Bit Error Rate (QBER) vs. Fiber Distance')
    plt.xlabel('Fiber Length (km)')
    plt.ylabel('QBER (%)')
    plt.grid(True)
    plt.axhline(y=11.0, color='r', linestyle='--', label='Theoretical Abort Threshold (11%)')
    plt.legend()
    
    # Save the figure
    plt.savefig('qber_vs_distance.png', dpi=300)
    print("Experiment completed. Plot saved as 'qber_vs_distance.png'.")

if __name__ == "__main__":
    run_qber_vs_distance_experiment()
