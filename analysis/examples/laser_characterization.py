"""RETIRED to analysis/examples/. Kept because it runs and because the
reasoning in it is worth reading, not because anything depends on it.
Nothing in src/ or run_all.py imports from this directory.

What it was: an eight-panel CWLaser dashboard -- power calibration,
phase noise, RIN, polarisation, and eye diagrams.

What replaced it: analysis/validation/validate_cwlaser.py covers the
spectra with pass/fail criteria, which this file never had.

The eye diagrams were the exception and were kept. They moved to
src/visualization/eye.py, which this file now imports, and are drawn by
validate_cwlaser.py and validate_dfb_drive.py. Porting them turned up
three faults that had always been there: no seed, so every run differed;
no polarisation control, so the X-cut MZM was modulating a component the
source had not put its light on; and an unfiltered drive, so the edges
were one sample wide and there was no eye opening at all.
"""
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import os, sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from src.lasers.cwlaser import CWLaser
from src.visualization.stokes import compute_stokes_parameters, poincare
from src.channel.mzm import MZM
# One implementation, and it is not this file's any more: the eye was
# ported to src/visualization/ when this script retired, and gained a
# seed, a polarisation controller and a band-limited drive there.
from src.visualization import eye_diagram

OUT = os.path.join(os.path.dirname(__file__))


def _lorentzian(f, f0, hwhm, S0):
    """Lorentzian lineshape: S(f) = S0 * (hwhm/pi) / ((f-f0)^2 + hwhm^2)"""
    return S0 * (hwhm / np.pi) / ((f - f0)**2 + hwhm**2)


# --------------------------------------------─
#  1. Power & field convention
# --------------------------------------------─
def plot_power_convention(laser, ax):
    """Verify mean(|E|^2) = P_W via laser.sample_field()."""
    E = laser.sample_field(1e-12, 5000)
    P_meas = np.mean(np.sum(np.abs(E)**2, axis=1))
    P_exp  = laser._power_w
    err_pc  = 100 * abs(P_meas - P_exp) / P_exp
    ax.bar(['Expected', 'Measured'], [P_exp, P_meas],
           color=['steelblue', 'tomato'], width=0.5)
    ax.set_ylabel('Power (W)')
    ax.set_title(f'Field power convention\n(error = {err_pc:.2f} %)')
    ax.yaxis.major.formatter.set_powerlimits((0, 4))
    return P_meas, P_exp


# --------------------------------------------─
#  2. Optical spectrum / linewidth
# --------------------------------------------─
def plot_optical_spectrum(laser, ax):
    """
    Complex-envelope PSD via Welch on the Ex component of sample_field().
    Fit a Lorentzian to extract the FWHM linewidth.
    """
    fs = 200e9
    dt = 1.0 / fs
    T  = 5e-6
    N  = int(T / dt)

    E = laser.sample_field(dt, N)
    Ex = E[:, 0] - np.mean(E[:, 0])

    f, Pxx = signal.welch(Ex, fs=fs, nperseg=min(N//8, 16384),
                          return_onesided=True, scaling='density')

    peak_idx = np.argmax(Pxx)
    f_peak = f[peak_idx]

    fit_mask = np.abs(f - f_peak) < 1e9
    f_fit = f[fit_mask]
    P_fit = Pxx[fit_mask]

    if len(f_fit) < 10:
        ax.text(0.5, 0.5, 'Insufficient resolution\nfor linewidth fit',
                transform=ax.transAxes, ha='center')
        return

    from scipy.optimize import curve_fit
    def lor(f, hwhm, S0):
        return _lorentzian(f, f_peak, hwhm, S0)

    try:
        popt, _ = curve_fit(lor, f_fit, P_fit,
                            p0=[laser.linewidth / 2, np.max(P_fit) * np.pi * laser.linewidth / 2],
                            bounds=([1e3, 1e-20], [1e12, 1e10]))
        hwhm_fit = popt[0]
        fwhm_fit = 2 * hwhm_fit
    except Exception:
        fwhm_fit = np.nan

    ax.plot(f * 1e-9, Pxx, lw=0.8, label='Measured spectrum')
    if not np.isnan(fwhm_fit):
        f_smooth = np.linspace(f[peak_idx] - 1e9, f[peak_idx] + 1e9, 500)
        ax.plot(f_smooth * 1e-9, lor(f_smooth, *popt), '--r',
                label=f'Lorentzian fit\nFWHM = {fwhm_fit/1e6:.2f} MHz')

    ax.axvline(laser.linewidth / 2 * 1e-9, color='grey', ls=':',
               label=f'Specified \u00bd-\u0394\u03bd = {laser.linewidth/2e6:.2f} MHz')

    ax.set_xlabel('Frequency offset (GHz)')
    ax.set_ylabel('PSD (a.u.)')
    ax.set_title('Optical spectrum (complex envelope)')
    ax.legend(fontsize=8)
    ax.set_yscale('log')
    ax.set_xlim(-1, 1)
    return fwhm_fit


# --------------------------------------------─
#  3. RIN spectrum — extracted from field power
# --------------------------------------------─
def theoretical_rin_psd(f, f_RO, gamma, RIN_0):
    """Coldren Eq 5.3.38 normalised to RIN_0 at DC."""
    omega   = 2 * np.pi * f
    omega_R = 2 * np.pi * f_RO
    num     = gamma**2 + omega**2
    den     = (omega_R**2 - omega**2)**2 + (gamma * omega)**2
    h2      = (num / den) / (gamma**2 / omega_R**4)
    return RIN_0 * h2


def plot_rin_spectrum(laser, ax):
    """RIN PSD via sample_field() power time-series vs. Coldren theory."""
    fs = 100e9
    dt = 1.0 / fs
    T  = 10e-6
    N  = int(T / dt)

    E = laser.sample_field(dt, N)
    power = np.sum(np.abs(E)**2, axis=1)
    rin = power / laser._power_w - 1.0

    f, Pxx = signal.welch(rin, fs=fs, nperseg=min(N//16, 8192),
                          return_onesided=True, scaling='density')

    ax.plot(f * 1e-9, 10 * np.log10(Pxx), lw=0.8, label='Measured (from field)')

    f_th = np.logspace(7, 10, 500)
    RIN_th = theoretical_rin_psd(f_th, laser._relaxation_freq,
                                 laser._damping_rate, laser._rin_linear)
    ax.plot(f_th * 1e-9, 10 * np.log10(RIN_th), '--k', lw=1.2,
            label='Coldren Eq 5.3.38')

    ax.axvline(laser._relaxation_freq * 1e-9, color='grey', ls=':',
               label=f'$f_{{RO}}$ = {laser._relaxation_freq/1e9:.1f} GHz')
    ax.axhline(laser.rin_density, color='grey', ls='--',
               label=f'DC RIN = {laser.rin_density:.0f} dB/Hz')

    ax.set_xlabel('Frequency (GHz)')
    ax.set_ylabel('RIN (dB/Hz)')
    ax.set_title('Relative Intensity Noise')
    ax.legend(fontsize=8)
    ax.set_xscale('log')
    ax.set_xlim(1e-1, f[-1] * 1e-9)


# --------------------------------------------─
#  4. Phase noise — extracted from field phase
# --------------------------------------------─
def plot_phase_noise(laser, ax):
    """
    Structure function D_phi(tau) = <[phi(t+tau) - phi(t)]^2>
    extracted from the complex envelope via sample_field().
    For a Wiener process D_phi(tau) = 2*pi*linewidth * tau.
    """
    fs = 100e9
    dt = 1.0 / fs
    T  = 2e-6
    N  = int(T / dt)

    E = laser.sample_field(dt, N)
    phi = np.unwrap(np.angle(E[:, 0]))

    max_tau = min(N // 4, 5000)
    lags = np.arange(1, max_tau)
    D    = np.zeros_like(lags, dtype=float)

    for i, tau in enumerate(lags):
        diff = phi[tau:] - phi[:-tau]
        D[i] = np.mean(diff**2)

    tau_vals = lags * dt

    ax.plot(tau_vals * 1e9, D, lw=0.8, label='Measured (from field)')
    ax.plot(tau_vals * 1e9, 2 * np.pi * laser.linewidth * tau_vals,
            '--k', lw=1.2, label=r'$2\pi \Delta\nu \cdot \tau$')

    ax.set_xlabel('Lag $\\tau$ (ns)')
    ax.set_ylabel('$D_\\phi(\\tau)$ (rad$^2$)')
    ax.set_title('Phase noise structure function')
    ax.legend(fontsize=8)


# --------------------------------------------─
#  5. Polarisation (Stokes)
# --------------------------------------------─
def plot_stokes_summary(laser, ax):
    """Display Stokes parameters as text."""
    E = laser.instantaneous_field(over_period=True, normalize=False)
    [S0, S1, S2, S3], [psi, chi] = compute_stokes_parameters(E)
    ax.axis('off')
    text = (f'Stokes parameters:\n'
            f'S0 = {S0:.4f}\n'
            f'S1 = {S1:.4f}\n'
            f'S2 = {S2:.4f}\n'
            f'S3 = {S3:.4f}\n\n'
            f'Polarisation ellipse:\n'
            f'psi = {np.rad2deg(psi):.1f} deg\n'
            f'chi = {np.rad2deg(chi):.1f} deg')
    ax.text(0.5, 0.5, text, transform=ax.transAxes, ha='center', va='center',
            fontsize=10, family='monospace')
    return S1, S2, S3


# --------------------------------------------─
#  Main
# --------------------------------------------─
def main():
    laser = CWLaser(
        wavelength=1550e-9,
        power_dbm=0,
        linewidth=1e6,
        rin_density=-130,
        polarization_azimuth=np.pi / 4,
        polarization_ellipticity=0.0,
        relaxation_frequency=5e9,
        damping_rate=1.88e10,
    )

    print(laser)
    print(f"  P_0 = {laser._power_w * 1e3:.4f} mW  ({laser.power_dbm:.0f} dBm)")

    fig = plt.figure(figsize=(14, 10))

    ax1  = fig.add_subplot(3, 3, 1)
    plot_power_convention(laser, ax1)

    ax2  = fig.add_subplot(3, 3, 2)
    plot_optical_spectrum(laser, ax2)

    ax3  = fig.add_subplot(3, 3, 3)
    plot_phase_noise(laser, ax3)

    ax4  = fig.add_subplot(3, 3, 4)
    plot_rin_spectrum(laser, ax4)

    ax5  = fig.add_subplot(3, 3, 5)
    S1, S2, S3 = plot_stokes_summary(laser, ax5)

    ax6  = fig.add_subplot(3, 3, 6)
    eye_diagram(laser, bitrate=5e9,  n_bits=128, spb=64, ax=ax6,
                title='5 Gbaud eye')

    ax7  = fig.add_subplot(3, 3, 7)
    eye_diagram(laser, bitrate=10e9, n_bits=128, spb=32, ax=ax7,
                title='10 Gbaud eye')

    ax8  = fig.add_subplot(3, 3, 8)
    eye_diagram(laser, bitrate=25e9, n_bits=128, spb=32, ax=ax8,
                title='25 Gbaud eye')

    ax9  = fig.add_subplot(3, 3, 9)
    ax9.axis('off')

    plt.tight_layout()
    path = os.path.join(OUT, 'laser_characterization.png')
    fig.savefig(path, dpi=150)
    print(f'\nSaved -> {path}')
    plt.close()

    poincare(S1, S2, S3)
    path3 = os.path.join(OUT, 'poincare_sphere.png')
    plt.gcf().savefig(path3, dpi=150)
    print(f'Saved -> {path3}')
    plt.close()

    fig2, axes = plt.subplots(1, 3, figsize=(15, 4))
    for ax, br, lbl in zip(axes,
                           [5e9, 10e9, 25e9],
                           ['5 Gbaud', '10 Gbaud', '25 Gbaud']):
        eye_diagram(laser, bitrate=br, n_bits=512, spb=64, ax=ax,
                    title=lbl)
    plt.tight_layout()
    path2 = os.path.join(OUT, 'eye_diagrams.png')
    fig2.savefig(path2, dpi=150)
    print(f'Saved -> {path2}')
    plt.close()


if __name__ == '__main__':
    main()
