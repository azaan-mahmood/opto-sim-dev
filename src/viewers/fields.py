def plot_field(E, frequency, title)->None:
    Ex = np.real(E[:, 0])
    Ey = np.real(E[:, 1])

    E_combined = np.sqrt(Ex ** 2 + Ey ** 2)

    plt.figure()
    t = np.linspace(0, 2 * np.pi / frequency, 1000)
    plt.subplot(3, 1, 1)
    plt.plot(t, Ex, label="Ex")
    plt.title(title)
    plt.xlabel("Time")
    plt.ylabel("Ex")
    plt.grid(True)

    plt.subplot(3, 1, 2)
    plt.plot(t, Ey, label="Ey", color="orange")
    plt.xlabel("Time")
    plt.ylabel("Ey")
    plt.grid(True)

    plt.subplot(3, 1, 3)
    plt.plot(t, E_combined, label="E", color="Red")
    plt.xlabel("Time")
    plt.ylabel("E")
    plt.grid(True)

    plt.tight_layout()
    plt.show()
