import numpy as np
import matplotlib.pyplot as plt

# --- LITERATURE SOURCES ---
# [1] Collett, E., "Field Guide to Polarization", SPIE Press, 2005, Ch. 2.
#     Stokes parameters definitions and the Poincaré sphere representation.
# [2] Hecht, E., "Optics", 4th ed., Addison-Wesley, 2002, Ch. 8.
#     Polarization ellipsometry and Stokes parameter formalism.
# [3] Born, M. & Wolf, E., "Principles of Optics", 7th ed., Cambridge, 1999.
#     Section 1.4: The polarization ellipse and Stokes parameters.


def compute_stokes_parameters(E: np.ndarray) -> tuple:
    """
    Calculate the Stokes parameters from the electric field over one period.

    Parameters:
    E (numpy array): The electric field array with shape (N, 2), where N is the number of time points,
                     and the columns represent the Ex and Ey components.

    Returns:
    tuple: A tuple containing the Stokes parameters [S0, S1, S2, S3] and the polarization ellipse parameters [psi, chi].
    """
    Ex = E[:, 0]
    Ey = E[:, 1]

    # Stokes parameters from time-averaged field correlations (Collett [1])
    S0 = np.real(np.mean(Ex*np.conj(Ex) + Ey*np.conj(Ey)))

    if S0 == 0:
        raise ValueError("S0 is zero — no optical power present")

    # Normalized Stokes parameters (Collett [1], Eq. 2.12–2.15)
    S1 = np.real(np.mean(Ex*np.conj(Ex) - Ey*np.conj(Ey))) / S0
    S2 = 2 * np.real(np.mean(Ex * np.conj(Ey))) / S0
    # Collett [1] Eq 2.15: S3 = 2 * Im(<Ex* conj(Ey)>), sign gives handedness
    S3 = -2 * np.imag(np.mean(Ex * np.conj(Ey))) / S0
    S0 = 1.0  # Normalize S0 to 1 for pure states, so S1, S2, S3 are relative values on the Poincaré sphere.    
    # Polarization ellipse parameters (Collett [1], Eq. 2.28–2.29)
    psi = 0.5 * np.arctan2(S2, S1)   #— orientation angle
    chi = 0.5 * np.arcsin(S3)        #— ellipticity angle (S3 already normalized)
    # Note: Clip to [-1, 1] guards against floating-point noise in S3
    # that can exceed unity and cause arcsin to return NaN (Born & Wolf [3]).

    return [S0, S1, S2, S3], [psi, chi]


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


def cos_sim(v1: np.ndarray, v2: np.ndarray) -> np.ndarray:
    """
    Returns the cosine similarity index between two vectors.
    Parameters:
        v1 (np.ndarrary): Vector 1, must be a np.column_stack type.
        v2 (np.ndarrary): Vector 2, must be a np.column_stack type.
    Raises:
        No Error
    """
    num = np.dot(v1, v2.T)
    mag_a = np.linalg.norm(v1)
    mag_b = np.linalg.norm(v2)
    sim = num/(mag_a*mag_b)
    return sim