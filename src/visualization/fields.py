import numpy as np
import matplotlib.pyplot as plt

# --- LITERATURE SOURCES ---
# [1] Hecht, E., "Optics", 4th ed., Addison-Wesley, 2002, Ch. 8.
#     Electric field representation of polarized light.


def plot_field(E, frequency, title)->None:
    Ex = np.real(E[:, 0])
    Ey = np.real(E[:, 1])

    E_combined = np.sqrt(Ex ** 2 + Ey ** 2)

    plt.figure()
    t = np.linspace(0, 2 * np.pi / frequency, 1000)
    plt.subplot(3, 1, 1)
    plt.plot(t, Ex, label="Ex")
    plt.ylabel("Ex")
    plt.grid(True)

    plt.subplot(3, 1, 2)
    plt.plot(t, Ey, label="Ey", color="orange")
    plt.ylabel("Ey")
    plt.grid(True)

    plt.subplot(3, 1, 3)
    plt.plot(t, E_combined, label="E", color="Red")
    plt.xlabel("Time")
    plt.ylabel("E")
    plt.grid(True)

    plt.suptitle(title)
    plt.tight_layout()
    plt.show()
