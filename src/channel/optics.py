import numpy as np


def coupler_split(power, E, ratio=0.5):
    if 0 <= ratio <= 1:
        port_power = power * ratio
        tap_power = power * (1-ratio)
        jones = np.array([
            [1, 1],
            [1, -1]
        ])
        port_E = E
        tap_E = jones @ np.transpose(E)
        return port_power, port_E, tap_power, np.transpose(tap_E)
    else:
        raise Exception("Incorrect Ratio")


def coupler_combine(power_port, port_E, power_tap, tap_E, out_ports=1):
    if out_ports==1:
        E_out = (port_E + tap_E)
        pout = power_port+power_tap
        return pout, E_out
    elif out_ports==2:
        E_out1 = (port_E + 1j * tap_E)
        E_out2 = (1j * port_E + tap_E)
        return power_port, E_out1, power_tap, E_out2
    else:
        raise Exception("Incorrect Number of Ports. Argument accepted is 1 or 2.")

def halfwave(E, theta=0, rotation=True):
    if not rotation:
        half_matrix = np.exp(-1j*np.pi/2)*np.array([
            [1, 0],
            [0, -1]
        ])
        E_half = half_matrix @ np.transpose(E)
        return np.transpose(E_half)
    elif rotation:
        theta = np.radians(theta)
        half_matrix = np.exp(-1j*np.pi/2)*np.array([
            [np.cos(theta)**2 - np.sin(theta)**2, 2*np.cos(theta)*np.sin(theta)],
            [2*np.cos(theta)*np.sin(theta), -np.cos(theta)**2 + np.sin(theta)**2]
        ])

        E_half = half_matrix @ np.transpose(E)
        return np.transpose(E_half)
    else:
        raise Exception("Incorrect Rotation of waveplate")


def quarterwave(E, theta=0, rotation=True):
    if not rotation:
        quarter_matrix = np.exp(-1j*np.pi/4)*np.array([
            [1, 0],
            [0, -1j]
        ])
        E_quarter = quarter_matrix @ np.transpose(E)
        return np.transpose(E_quarter)
    elif rotation:
        theta = np.radians(theta)
        quarter_matrix = np.exp(-1j*np.pi/4)*np.array([
            [np.cos(theta) ** 2 + 1j * np.sin(theta) ** 2, (1 - 1j) * np.sin(theta) * np.cos(theta)],
            [(1 - 1j) * np.sin(theta) * np.cos(theta), np.sin(theta) ** 2 + 1j * np.cos(theta) ** 2]
        ])
        E_quarter = quarter_matrix @ np.transpose(E)
        return np.transpose(E_quarter)
    else:
        raise Exception("Incorrect Rotation of waveplate")


def polarization_rotator(E, rot_angle):
    theta = np.radians(rot_angle)
    rot_matrix = np.array([
        [np.cos(theta), np.sin(theta)],
        [-np.sin(theta), np.cos(theta)]
    ])

    E_rot = rot_matrix @ np.transpose(E)
    return np.transpose(E_rot)

def polarization_controller(E, qwp1=0, hwp=0,qwp2=0):
    E = quarterwave(E, qwp1)
    E = halfwave(E, hwp)
    E = quarterwave(E, qwp2)
    return E


def voa(E, attenuation_dB):
    """Variable optical attenuator.

    Parameters
    ----------
    E : ndarray (N, 2) — complex-envelope field.
    attenuation_dB : float — attenuation in dB (positive = loss).

    Returns
    -------
    ndarray (N, 2) — attenuated field.
    """
    atten_lin = 10.0 ** (-attenuation_dB / 10.0)
    return E * np.sqrt(atten_lin)


def polarizer(E, polarization):
    if polarization == 'H':
        J_H = np.array([
            [0.99, 0],
            [0, 0.01]
        ])
        return np.transpose(J_H @ np.transpose(E))
    elif polarization == 'V':
        J_V = np.array([
            [0.001, 0],
            [0, 0.99]
        ])
        return np.transpose(J_V @ np.transpose(E))
    elif polarization == '45':
        J_45 = np.array([
            [1 / 2, 1 / 2],
            [1 / 2, 1 / 2]
        ])
        # J_45 = 0.5*np.array([
        #     [1, 1],
        #     [1, 1]
        # ])
        return np.transpose(J_45 @ np.transpose(E))
    elif polarization == '-45':
        J_i45 = np.array([
            [1 / 2, -1 / 2],
            [-1 / 2, 1 / 2]
        ])
        # J_i45 = np.array([
        #     [1, -1],
        #     [-1, 1]
        # ])
        return np.transpose(J_i45 @ np.transpose(E))
    else:
        raise Exception("Wrong Value")

#Version 1
# def beam_splitter(E):
#     # E = E * np.sqrt(power)
#     J_x =(1/np.sqrt(2)) * np.array([[0.95, 0],
#                                     [0, 0.05]])
#
#     J_y = (1/np.sqrt(2)) * np.array([[0.05, 0],
#                                     [0, 0.95]])
#
#     # Calculate the transmitted and reflected components
#     E_x = np.transpose(J_x @ np.transpose(E))
#     E_y = np.transpose(J_y @ np.transpose(E))
#     E_y = halfwave(E_y, rotation=False)
#
#     # Calculate the power in each component
#     P_x = np.linalg.norm(E_x) ** 2
#     P_y = np.linalg.norm(E_y) ** 2
#     print(P_x, P_y)
#     # return P_x, E_x, P_y, E_y
#     return E_x, E_y


def hadamard(E):
    # E is assumed to be an (N,2) array.
    J = (1 / np.sqrt(2)) * np.array([[1, 1],
                                     [1, -1]])
    E_out = E @ J.T  # Apply transformation along the last axis
    return E_out

def circular_analyser(E):
    """Circular-basis analyser: quarter-wave plate followed by a
    polarising beam splitter (QWP+PBS), projecting onto the L/R circular
    basis rather than H/V.

    This is *not* a polarising beam splitter (PHYS-6 in
    opto-sim-issues-and-fixes.md — this function used to be misnamed
    `pbs`). A true PBS projects onto H and V (`diag(1,0)` / `diag(0,1)`,
    see `pbs()` below) and is blind to the relative phase between Ex and
    Ey. This function instead converts that relative phase into an
    intensity imbalance between its two output ports:

        E_x = (Ex - i*Ey) / sqrt(2),   E_y = (-i*Ex + Ey) / sqrt(2)

    which is what phase-encoded schemes (e.g. Section 5's 45°-polarised,
    phase-modulated BB84) need for detection — a true H/V PBS would give
    a fixed 50/50 split regardless of phase and discriminate nothing.
    Any caller relying on this behaviour must call `circular_analyser()`
    explicitly, not `pbs()`.

    Parameters
    ----------
    E : ndarray (N, 2) — complex envelope [Ex, Ey].

    Returns
    -------
    E_x, E_y : ndarray (N,) — the two output-port fields.
    """
    J = (1 / np.sqrt(2)) * np.array([[1, -1j],
                                    [-1j, 1]])
    E_out = E @ J.T  # Apply transformation along the last axis
    E_x = E_out[:, 0]
    E_y = E_out[:, 1]
    return E_x, E_y


def pbs(E):
    """Polarising beam splitter: projects onto H and V.

    A true PBS transmits the H component to one port and reflects the V
    component to the other, and is blind to the relative phase between Ex
    and Ey — unlike `circular_analyser()` (formerly misnamed `pbs`, see
    PHYS-6 in opto-sim-issues-and-fixes.md), which projects onto the
    circular basis and *does* depend on that phase.

    Parameters
    ----------
    E : ndarray (N, 2) — complex envelope [Ex, Ey].

    Returns
    -------
    E_x, E_y : ndarray (N,) — the H-port and V-port fields.
    """
    E_x = E[:, 0]
    E_y = E[:, 1]
    return E_x, E_y
#
# Version 2
# def beam_splitter(E, power):
#     # Scale electric field by the square root of power
#     # E = E * np.sqrt(power)
#
#     # Jones matrices for transmitted and reflected components
#     J_t = (1 / np.sqrt(2)) * np.array([[1, 0],
#                                        [0, 1]])
#     J_r = (1 / np.sqrt(2)) * np.array([[1, 0],
#                                        [0, -1]])
#
#     # Calculate transmitted (E_x) and reflected (E_y) fields
#     E_x = np.transpose(J_t @ np.transpose(E))  # Transmitted
#     E_y = np.transpose(J_r @ np.transpose(E))  # Reflected
#
#     # Optional: apply a half-wave plate or other transformations to E_y
#     # E_y = halfwave(E_y, rotation=False) if needed
#
#     # Calculate the power in each component
#     P_x = np.linalg.norm(E_x) ** 2
#     P_y = np.linalg.norm(E_y) ** 2
#
#     return P_x, E_x, P_y, E_y

# Version 3
# def beam_splitter(E):
#     E_x_normalized = E[:, 0]
#     E_y_normalized = E[:, 1]
#     E_x_normalized = np.asarray(E_x_normalized, dtype=complex)
#     E_y_normalized = np.asarray(E_y_normalized, dtype=complex)
#     print(f"E_x_normalized[0]: {E_x_normalized[0]} (Type: {type(E_x_normalized[0])})")
#     print(f"E_y_normalized[0]: {E_y_normalized[0]} (Type: {type(E_y_normalized[0])})")
#     theta = np.angle(E_y_normalized[0]) - np.angle(E_x_normalized[0])
#     print(f'Theta = {theta}')
#
#     T = np.cos(theta)**2  # Transmission coefficient
#     R = np.sin(theta)**2  # Reflection coefficient
#     print(f'T = {T}, R = {R}')
#     # T = 0.5
#     # R = 0.5
#     T_matrix = np.array([[T, 0],
#                          [0, R]])
#
#     E_out = np.transpose(T_matrix @ np.transpose(E))
#     E_x = E_out[:, 0]
#     E_y = E_out[:, 1]
#     P_H = np.linalg.norm(E_x) ** 2
#     P_V = np.linalg.norm(E_y) ** 2
#     return P_H, E_out[:, 0], P_V, E_out[:, 1]

def beam_combiner(P_x, E_x, P_y, E_y, normalized=True):
    if normalized and np.linalg.norm(E_x) ** 2 < 1 and np.linalg.norm(E_y) ** 2 < 1:
        E = E_x + E_y
        pout = P_x+P_y
        return E, pout
    elif not normalized:
        E = E_x + E_y
        P_out = np.linalg.norm(E) ** 2
        return E, P_out
