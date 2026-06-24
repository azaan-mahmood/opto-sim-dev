import numpy as np
from .phase_modulator import PhaseModulator


class MZM:
    """
    Mach-Zehnder Modulator built from PhaseModulator-defined arm phase shifts.

    Physical model: an MZI where the input field is split into two arms
    (Y-branch), each arm may contain a PhaseModulator that imparts a voltage-
    dependent phase via the Pockels effect, and the arms recombine at a
    second Y-branch to produce interference.  The real Y-branch has finite
    extinction ratio from imperfect 50:50 splitting; insertion loss is a
    separate power-scaling parameter.

    Two electrode configurations are supported:

      *push-pull* (default)
          Both arms are modulated with opposite voltages.  No residual
          frequency chirp because the common-mode phase cancels.

          Transfer function per polarisation component:
              phi_1  =  pi * (V/2 + V_bias/2) / V_pi
              phi_2  =  pi * (-V/2 - V_bias/2) / V_pi
              E_out  =  sqrt(IL) * E_in *
                        cos( pi * (V + V_bias) / (2 * V_pi) )

      *single-drive*
          Only arm 1 is modulated; arm 2 is a passive reference waveguide.
          Residual frequency chirp appears as exp(j * pi * V / (2 * V_pi)).

          Transfer function per polarisation component:
              phi_1  =  pi * (V + V_bias) / V_pi
              phi_2  =  0
              E_out  =  sqrt(IL) * E_in * exp(j * pi * (V + V_bias) / (2 * V_pi)) *
                        cos( pi * (V + V_bias) / (2 * V_pi) )

    The crystal cut (X or Y) of the internal PhaseModulator determines
    which transverse field component (Ey or Ex) acquires the phase shift;
    the orthogonal component passes through unchanged.  This is physically
    correct for LiNbO3 modulators where the RF field is aligned to either
    the extraordinary or ordinary axis.

    Parameters
    ----------
    mode : {'push-pull', 'single-drive'}, default 'push-pull'
        Electrode configuration.
    pm : PhaseModulator or None
        PhaseModulator used by the modulated arm(s).  If None (default), a
        default X-cut DC PhaseModulator is created from crystal parameters.
        In push-pull mode, the same PhaseModulator type is used for both
        arms; in single-drive mode only arm 1 uses it.
    insertion_loss_db : float or None
        Total device insertion loss in dB.  None (default) means ideal
        (0 dB, no insertion loss).
    extinction_ratio_db : float or None
        Extinction ratio in dB, set by the Y-branch splitting imbalance.
        None (default) means ideal (infinite extinction, perfect 50:50
        splitting).
    bias_voltage : float
        Static DC bias voltage that shifts the operating point.  V_bias = 0
        biases the MZM at minimum transmission (null); V_bias = V_pi/2
        biases at quadrature (50 % transmission).  Default 0.0.

    References
    ----------
    [1]  Agrawal, G. P., "Fiber-Optic Communication Systems", 4th ed.,
         Wiley, 2010, §4.2: External Modulation and Mach-Zehnder Modulators.
    [2]  Koyama, M. and Iga, K., "Frequency chirping in external modulators",
         J. Lightwave Technol. 6(1), 87–93 (1988).
         --- single-drive vs. push-pull chirp comparison.
    [3]  LiNbO3 Pockels coefficients from Weis and Gaylord, "Lithium Niobate:
         Summary of Physical Properties and Crystal Structure", Appl. Phys. A
         37, 191–203 (1985).
    """

    def __init__(self, mode='push-pull', pm=None,
                 insertion_loss_db=None, extinction_ratio_db=None,
                 bias_voltage=0.0):
        if pm is None:
            pm = PhaseModulator(crystal_cut='X', modulation='RF')
        self._pm = pm
        self._mode = mode
        self._Vpi = pm.Vpi

        # Insertion loss
        if insertion_loss_db is None:
            self._il_lin = 1.0
        else:
            self._il_lin = 10.0 ** (-insertion_loss_db / 10.0)

        # Extinction ratio → Y-branch splitting imbalance
        if extinction_ratio_db is None:
            self._r = 0.5
        else:
            delta = 10.0 ** (-extinction_ratio_db / 20.0)
            self._r = 0.5 * (1.0 + delta)
        self._one_minus_r = 1.0 - self._r

        self._bias_voltage = bias_voltage

    # ── properties ──────────────────────────────────────────────────────

    @property
    def V_pi(self):
        """Half-wave voltage (voltage for ON-to-OFF transition)."""
        return self._Vpi

    @property
    def switching_voltage(self):
        """Alias for V_pi — voltage required for full ON-to-OFF extinction."""
        return self._Vpi

    @property
    def mode(self):
        return self._mode

    @mode.setter
    def mode(self, value):
        if value not in ('push-pull', 'single-drive'):
            raise ValueError(f"mode must be 'push-pull' or 'single-drive', got {value!r}")
        self._mode = value

    # ── modulation ──────────────────────────────────────────────────────

    def modulate(self, E_in, V):
        """
        Apply the modulation voltage to the optical field.

        Parameters
        ----------
        E_in : ndarray, shape (..., 2)
            Input optical field [Ex, Ey].  May be a single sample (2,) or
            a sequence of N samples (N, 2).
        V : float or ndarray, shape (N,)
            Modulation voltage(s).  A scalar is applied uniformly to all
            samples; a 1-D array is applied sample-by-sample.

        Returns
        -------
        ndarray, shape (..., 2)
            Modulated optical field.
        """
        V_pi = self._Vpi
        V = np.asarray(V, dtype=float)
        Vb = self._bias_voltage

        if self._mode == 'push-pull':
            V1 = V / 2.0 + Vb / 2.0
            V2 = -V / 2.0 - Vb / 2.0
        else:
            V1 = V + Vb
            V2 = 0.0

        phi1 = np.pi * V1 / V_pi
        phi2 = np.pi * V2 / V_pi

        # Y-branch split
        E_arm1 = E_in * np.sqrt(self._r)
        E_arm2 = E_in * np.sqrt(self._one_minus_r)

        # Phase modulation in each arm (crystal cut determines Ex vs Ey)
        E_arm1 = self._apply_phase(E_arm1, phi1)
        E_arm2 = self._apply_phase(E_arm2, phi2)

        # Y-branch recombination + insertion loss
        E_out = np.sqrt(self._il_lin) * (
            np.sqrt(self._r) * E_arm1 + np.sqrt(self._one_minus_r) * E_arm2
        )
        return E_out

    # ── helpers ─────────────────────────────────────────────────────────

    def _apply_phase(self, E, phi):
        """Apply phase shift ``phi`` to the polarisation component
        determined by the internal PhaseModulator's crystal cut."""
        E_out = E.copy()
        if self._pm.crystal_cut == 'X':
            E_out[..., 1] *= np.exp(1j * phi)
        else:
            E_out[..., 0] *= np.exp(1j * phi)
        return E_out

    def __repr__(self):
        return (f"MZM(mode={self._mode!r}, "
                f"V_pi={self._Vpi:.3f} V, "
                f"cut={self._pm.crystal_cut})")
