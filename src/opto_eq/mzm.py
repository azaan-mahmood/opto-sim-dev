import numpy as np


class MZM:
    """
    Mach-Zehnder Modulator in push-pull configuration.

    Models the interferometric conversion of a phase modulation into an
    amplitude (intensity) modulation.  The input field is split into two
    arms, each arm receives a phase shift proportional to the applied
    voltage, and the arms recombine:

        E_out = 0.5 * E_in * [exp(j * phi_1) + exp(j * phi_2)]

    For push-pull operation (phi_1 = -phi_2 = pi * V / V_pi):
        E_out = E_in * cos(pi * V / V_pi)

    where V_pi is the single-arm voltage for a pi phase shift.
    The switching voltage (null-to-null) is V_pi / 2.

    Parameters
    ----------
    V_pi : float
        Voltage required for a pi phase shift in one arm (Volts).
        The intensity goes from maximum at V=0 to minimum at V=V_pi/2.
        Default 5.0 V.
    bias_voltage : float
        DC bias voltage applied equally to both arms (Volts).
        Adds a common-mode phase exp(j * pi * V_bias / V_pi).
        Default 0.0 (null bias — maximum extinction at V=0).

    References
    ----------
    [1] Agrawal, G. P., "Fiber-Optic Communication Systems", 4th ed.,
        Wiley, 2010, §4.2: External Modulation and Mach-Zehnder Modulators.
    """
    def __init__(self, V_pi: float = 5.0, bias_voltage: float = 0.0):
        self.V_pi = V_pi
        self.bias_voltage = bias_voltage

    def modulate(self, E_in: np.ndarray, V_signal: np.ndarray | float
                 ) -> np.ndarray:
        """
        Apply push-pull modulation to an optical field.

        Parameters
        ----------
        E_in : ndarray, shape (N, 2)
            Input optical field [Ex, Ey] at N time samples.
        V_signal : ndarray, shape (N,) or float
            Applied signal voltage at each sample.

        Returns
        -------
        ndarray, shape (N, 2)
            Modulated optical field.
        """
        phi_bias = np.pi * self.bias_voltage / self.V_pi
        phi_signal = np.pi * np.asarray(V_signal) / self.V_pi
        tf = np.cos(phi_signal) * np.exp(1j * phi_bias)
        return E_in * tf[:, np.newaxis]

    @property
    def switching_voltage(self) -> float:
        """Voltage for a full ON-to-OFF transition (null)."""
        return self.V_pi / 2.0
