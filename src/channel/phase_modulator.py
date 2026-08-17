import numpy as np
from typing import Dict

class PhaseModulator:
    def __init__(self, crystal_cut='X', modulation="DC", params=None,
                 phase_error_rad=0.0, phase_noise_rad=0.0, rng=None,
                 bias_offset_v=0.0) -> None:
        """
        Initialize a Phase Modulator using the Pockels effect in LiNbO3.

        Parameters:
            crystal_cut (str): 'X' or 'Y'. Determines modulation axis.
            modulation (str): 'DC' or 'RF'. Only 'DC' is implemented currently.
            params (dict, optional): Dictionary of modulator parameters.
                Accepts keys: 'wavelength', 'n_o', 'n_e', 'd', 'L', 'Gamma', 'r13', 'r22'.
            phase_error_rad (float): STATIC phase offset added to every
                modulation, in radians (default 0). Models a calibration
                offset -- a systematically mis-set Vpi.
            phase_noise_rad (float): per-call Gaussian phase jitter, sigma
                in radians (default 0). Models drive-voltage noise, which
                is random pulse to pulse -- see "Modulation error" below.
            rng (numpy.random.Generator, optional): source for the jitter.
                Defaults to the global numpy.random state.
            bias_offset_v (float): STATIC bias error expressed as a drive
                voltage (default 0). Converted through this modulator's own
                V_pi as phi = pi*V/V_pi, so it is the same mechanism as
                `phase_error_rad` stated in the units the hardware is
                actually set in. Supplying both raises.

        Modulation error
        ----------------
        A modulator's bias is never set exactly. Two knobs describe the
        error, and for a given magnitude they give the same QBER floor to
        within about 1%:

            static offset d:   QBER = (1 - cos d) / 2
            Gaussian jitter s: QBER = (1 - exp(-s^2/2)) / 2

        `bias_offset_v` is the static form, in volts. On the default
        crystal parameters V_pi = 3.8826 V, so a 20.93 deg offset is
        452 mV. `phase_noise_rad` is the jitter form, in radians: drive
        voltage noise resampled per pulse, default 0.

        Which one applies is a statement about the hardware. A bias set
        once and left is static; noise on the drive electronics is
        per-pulse random. Drive electronics typically deliver 0.1 to 1% of
        V_pi, which is 0.0002 to 0.025% QBER, so jitter is a small term
        unless the modulator is unusually noisy.

        Interferometer arm-length drift is a separate mechanism and lives
        on `AsymmetricMZI.phase_drift_rad_s`, not here.

        NOT MODELLED: LiNbO3 DC bias drift (charge migration, defect-bound
        charge, the LN/SiO2 RC loop). Its relaxation runs minutes to days,
        so over runs of this length it is indistinguishable from a constant
        offset and is already absorbed by `bias_offset_v`. Implement it when
        a protocol simulates tens of minutes or longer.

        Raises:
            RuntimeError: If unknown crystal_cut, modulation, or parameter key is provided.

        Important Notes:
            - RF modulation is not yet implemented.
            - Based on Z-propagation with Ex, Ey transverse field simulation.
            - X-Cut means Y surface, so phase applied to Ey.
            - Y-Cut means X surface, so phase applied to Ex.
            - X-Cut uses r13 coefficient (Ey modulation), Y-Cut uses r22 (Ex modulation).

        V_π formula:
            V_π = λ·d / (2·n_o³·r·Γ·L)

            This is the voltage needed by the *push-pull MZM* to achieve a
            π relative phase shift between its arms. A single arm alone would
            require 2× this voltage for π phase shift.

        References
        ----------
        [1] Weis, R. S. & Gaylord, T. K., "Lithium Niobate: Summary of
            Physical Properties and Crystal Structure", Appl. Phys. A 37,
            191-203, 1985.  — LiNbO3 Pockels coefficients r13, r22.
        [2] Alferness, R. C., "Titanium-diffused lithium niobate waveguide
            devices", in Guided-Wave Optoelectronics, T. Tamir (ed.),
            Springer, 1988, Ch. 4.  — V_π formula for LiNbO3 modulators.
        """
        self.default_params = {
            "wavelength": 1550e-9,   # Wavelength (in nm)
            "n_o": 2.2,     # Ordinary refractive index of LiNbO3
            "n_e": 2.14,     # Extraordinary refractive index of LiNbO3, only used for non Z-propagating.
            "d": 24.133e-6,       # Separation of the Electrodes, also called Se sometimes
            "L": 5.588e-2,       # Length of LiNbO3 film
            "Gamma": 0.8,   # Physical/Spatial Overlap factor between electrode and waveguide
            "r13": 10.12e-12,     # For X-Cut Config
            "r22": 3.4e-12     # For Y-Cut Config
        }
        if params is None:
            params = self.default_params
        else:
            params = {**self.default_params, **params}

        for key in params:
            if key not in self.default_params:
                raise RuntimeError(f"Unidentified Parameter: {key}")

        if crystal_cut not in ["X", "Y"]:
            raise RuntimeError(f"Unidentified Crystal Cut: {crystal_cut}")
        else:
            self.crystal_cut = crystal_cut

        # --- MODIFIED: 2026-05-06 16:02:19 ---
        if modulation not in ["DC", "RF"]:
            raise RuntimeError(f"Unidentified Modulation Type: {modulation}")
        else:
            self.modulation = modulation
        # --- END MODIFIED ---

        for key, value in params.items():
            setattr(self, key, value)

        self.__Vpi = self.get_vpi()

        # A bias error is one mechanism; accepting it twice in two unit
        # systems would silently sum them.
        if bias_offset_v != 0.0 and phase_error_rad != 0.0:
            raise RuntimeError(
                "bias_offset_v and phase_error_rad describe the same static "
                "bias error in different units; supply one, not both.")
        if bias_offset_v != 0.0:
            phase_error_rad = np.pi * bias_offset_v / self.__Vpi

        self.bias_offset_v = float(bias_offset_v)
        self.phase_error_rad = float(phase_error_rad)
        self.phase_noise_rad = float(phase_noise_rad)
        self._rng = rng

    def _phase_imperfection(self):
        """Static offset plus a fresh jitter draw (0.0 when both are off)."""
        extra = self.phase_error_rad
        if self.phase_noise_rad > 0.0:
            rng = self._rng if self._rng is not None else np.random
            extra += rng.normal(0.0, self.phase_noise_rad)
        return extra

    def __repr__(self):
        return (f"PhaseModulator(cut={self.crystal_cut}, "
                f"modulation={self.modulation}, Vpi={self.Vpi:.2f} V, "
                f"phase_error={self.phase_error_rad:.4f} rad, "
                f"phase_noise={self.phase_noise_rad:.4f} rad)")

    def get_vpi(self)->float:
        """
        Gets half-angle voltage of the modulator
        :return: Voltage Vπ
        """
        if self.crystal_cut == "X":
            bot = (2 * self.n_o ** 3 * self.r13 * self.Gamma * self.L)
            if bot == 0:
                raise ZeroDivisionError("Invalid Parameters: Denominator of Vπ calculation is Zero!")
            return (self.wavelength * self.d) / bot
        else:
            bot = (2 * self.n_o ** 3 * self.r22 * self.Gamma * self.L)
            if bot == 0:
                raise ZeroDivisionError("Invalid Parameters: Denominator of Vπ calculation is Zero!")
            return (self.wavelength * self.d) / bot

    def get_phi(self, V):
        """
        Gets the phi value of the modulator
        :param V: Modulation voltage (in Volts)
        :return: Phase shift φ (in radians)
        """
        if np.size(V) == 1:
            return (np.pi * V) / self.__Vpi

    @property
    def Vpi(self)->float:
        return self.__Vpi

    def modulate(self, E_field: np.ndarray, V:float)->np.ndarray:
        """
        Modulate the E_field depending on Modulation Type
        :param E_field:
        :param V: Modulation DC Voltage
        :return: Modulated E_Field as numpy array of (N, 2)
        """
        if E_field.shape[-1] != 2:
            raise ValueError("E_field must be a 2D array with shape (N, 2)")

        if self.modulation == "DC" and np.size(V) == 1:
            phi = (np.pi * V) / self.__Vpi + self._phase_imperfection()
            if self.crystal_cut == "X":
                pm = np.array([
                    [1, 0],
                    [0, np.exp(1j * phi)]
                ])
                return np.transpose(pm @ np.transpose(E_field))
            else:
                pm = np.array([
                    [np.exp(1j * phi), 0],
                    [0, 1]
                ])
                return np.transpose(pm @ np.transpose(E_field))
        # --- MODIFIED: 2026-05-06 16:02:19 ---
        elif self.modulation == "RF" and np.size(V) > 1:
            if len(V) != E_field.shape[0]:
                raise ValueError("For RF modulation, V array length must match E_field time steps (N).")
            
            # Calculate time-varying phase shift
            phi = ((np.pi * V) / self.__Vpi  # phi is now an array of size N
                   + self._phase_imperfection())
            
            # Initialize output array
            E_out = np.zeros_like(E_field, dtype=complex)
            
            if self.crystal_cut == "X":
                # Phase shift applied to Ey
                E_out[:, 0] = E_field[:, 0]
                E_out[:, 1] = E_field[:, 1] * np.exp(1j * phi)
            else: # Y-Cut
                # Phase shift applied to Ex
                E_out[:, 0] = E_field[:, 0] * np.exp(1j * phi)
                E_out[:, 1] = E_field[:, 1]
                
            return E_out
        # --- END MODIFIED ---
