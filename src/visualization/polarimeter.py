import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.patches import Ellipse
import cmath


def polarimeter(E:np.ndarray, component:str='r', title:str="Polarization Ellipse")->None:
    """
    Plot the polarization ellipse with an arrow indicating handedness.

    Parameters:
    E = The electric field array with shape (N, 2), where N is the number of time points,
        and the columns represent the Ex and Ey components.
    component = Accepts a character. 'r' = real, 'i' = imaginary, 'c' = both real and imaginary
    """
    Ex = E[:, 0]  # Electric field component Ex
    Ey = E[:, 1]  # Electric field component Ey

    plt.figure()

    # Plot based on selected component
    if component == 'r':
        plt.plot(Ex.real, Ey.real, 'r-', label='Real')  # Plot real parts
        X, Y = Ex.real, Ey.real
    elif component == 'i':
        plt.plot(Ex.imag, Ey.imag, 'b-', label='Imaginary')  # Plot imaginary parts
        X, Y = Ex.imag, Ey.imag
    elif component == 'c':
        plt.plot(Ex.real, Ey.real, 'r-', label='Real')  # Plot real parts
        plt.plot(Ex.imag, Ey.imag, 'b-', label='Imaginary')  # Plot imaginary parts
        X, Y = Ex.real, Ey.real
    else:
        raise Exception("Invalid Component")

    # Determine arrow direction
    idx = len(X) // 3  # Pick an index for arrow placement
    dx = X[idx + 1] - X[idx]  # Change in x
    dy = Y[idx + 1] - Y[idx]  # Change in y

    # Place an arrow in the direction of rotation
    plt.quiver(X[idx], Y[idx], dx, dy, angles='xy', scale_units='xy', scale=0.5, color='red')

    # Plot formatting
    plt.title(title)
    plt.xlabel("Ex")
    plt.ylabel("Ey")
    plt.axis('equal')
    plt.grid(True)
    plt.legend()
    plt.show()
