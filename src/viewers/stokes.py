import numpy as np
import matplotlib.pyplot as plt

def compute_stokes_parameters(E):
    """
    Calculate the Stokes parameters from the electric field over one period.

    Parameters:
    E (numpy array): The electric field array with shape (N, 2), where N is the number of time points,
                     and the columns represent the Ex and Ey components.

    Returns:
    S0, S1, S2, S3 (float): The Stokes parameters.
    """
    Ex = E[:, 0]
    Ey = E[:, 1]

    S0 = np.real(np.mean(Ex*np.conj(Ex) + Ey*np.conj(Ey)))
    S1 = np.real(np.mean(Ex*np.conj(Ex) - Ey*np.conj(Ey))) / S0
    S2 = 2 * np.real(np.mean(Ex * np.conj(Ey))) / S0
    S3 = -2 * np.imag(np.mean(Ex * np.conj(Ey))) / S0

    # Avoid division by zero
    if S0 == 0:
        raise ValueError("S0 is zero, cannot compute Psi and Chi")

    # Calculate Psi and Chi
    psi = 0.5 * np.arctan2(S2, S1)
    chi = 0.5 * np.arcsin(S3/S0)
    # Output results
    # print(f"S0 = {S0:.3f}\nS1 = {S1:.3f}\nS2 = {S2:.3f}\nS3 = {S3:.3f}")
    # print(f"Psi (polarization azimuth) = {np.rad2deg(psi):.3f}°")
    # print(f"Chi (polarization ellipticity) = {np.rad2deg(chi):.3f}°")
    return S0, S1, S2, S3


def poincare(s1, s2, s3):
    """Plot the polarization point on the Poincaré sphere."""
    # Create a sphere
    u, v = np.mgrid[0:2*np.pi:100j, 0:np.pi:50j]
    x = np.cos(u) * np.sin(v)
    y = np.sin(u) * np.sin(v)
    z = np.cos(v)

    # Plot
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(x, y, z, color='lightblue', alpha=0.3, linewidth=0)

    # Axes
    ax.quiver(0, 0, 0, 1, 0, 0, color='r', label='S1', arrow_length_ratio=0.1)
    ax.quiver(0, 0, 0, 0, 1, 0, color='g', label='S2', arrow_length_ratio=0.1)
    ax.quiver(0, 0, 0, 0, 0, 1, color='b', label='S3', arrow_length_ratio=0.1)

    # Plot the Stokes vector (normalized)
    ax.scatter(s1, s2, s3, color='k', s=100, label='Polarization State')

    # Labels and formatting
    ax.set_xlabel('S1')
    ax.set_ylabel('S2')
    ax.set_zlabel('S3')
    ax.set_title("Poincaré Sphere")
    ax.legend()
    ax.set_box_aspect([1, 1, 1])
    plt.tight_layout()
    plt.show()