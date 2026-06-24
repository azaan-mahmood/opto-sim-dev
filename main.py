from src.deprecated import sslaser as laser
from src.visualization.stokes import compute_stokes_parameters, poincare, cos_sim
from src.visualization import polarimeter
from src.channel import polarizer, PhaseModulator, cable
import numpy as np

source = laser.SolidStateLaser(
    wavelength=1550e-9,  # Laser wavelength
    polarization_azimuth=np.pi,  # 45° polarization
    polarization_ellipticity=np.pi/4,
    power_dbm=-5,  # arbitrary power unit
    frequency=5e6
)

phase_modulator = PhaseModulator()
Vpi = phase_modulator.Vpi
E = source.get_electric_field(normalize=False, over_period=True)
stokes, polarization = compute_stokes_parameters(E)
S0, S1, S2, S3 = stokes
psi, chi = polarization
print(f"-----------------Stokes of Source--------------------")
print(f"S0 = {S0:.3f}\nS1 = {S1:.3f}\nS2 = {S2:.3f}\nS3 = {S3:.3f}")
print(f"Psi (polarization azimuth) = {np.rad2deg(psi):.3f}°")
print(f"Chi (polarization ellipticity) = {np.rad2deg(chi):.3f}°")
print(f"------------------Stokes after 45---------------------")
E = polarizer(E, polarization="45")
stokes, polarization = compute_stokes_parameters(E)
S0, S1, S2, S3 = stokes
psi, chi = polarization
print(f"S0 = {S0:.3f}\nS1 = {S1:.3f}\nS2 = {S2:.3f}\nS3 = {S3:.3f}")
print(f"Psi (polarization azimuth) = {np.rad2deg(psi):.3f}°")
print(f"Chi (polarization ellipticity) = {np.rad2deg(chi):.3f}°")
# polarimeter(E, title="Source at 45 LP")
# poincare(S1, S2, S3)
print(f"------------------Stokes after First Phase Modulator---------------------")
E = phase_modulator.modulate(E_field=E, V=0)
stokes, polarization = compute_stokes_parameters(E)
S0, S1, S2, S3 = stokes
psi, chi = polarization
stokes_fpm = np.column_stack([S0, S1, S2, S3])
print(f"S0 = {S0:.3f}\nS1 = {S1:.3f}\nS2 = {S2:.3f}\nS3 = {S3:.3f}")
print(f"Psi (polarization azimuth) = {np.rad2deg(psi):.3f}°")
print(f"Chi (polarization ellipticity) = {np.rad2deg(chi):.3f}°")
polarimeter(E, title=f"First PM V = {0:.2f}")
print(f"------------------Stokes after QC---------------------")
# NOTE: dispersion=True requires the complex envelope (sample_field()).
# The over-period field below has only ~5 fs of data — CD on this has
# no modulation bandwidth to disperse.  Pass a placeholder dt for now.
E = cable(100, E, dt=1e-12, dispersion=True)
stokes, polarization = compute_stokes_parameters(E)
S0, S1, S2, S3 = stokes
psi, chi = polarization
stokes_QC = np.column_stack([S0, S1, S2, S3])
print(f"S0 = {S0:.3f}\nS1 = {S1:.3f}\nS2 = {S2:.3f}\nS3 = {S3:.3f}")
print(f"Psi (polarization azimuth) = {np.rad2deg(psi):.3f}°")
print(f"Chi (polarization ellipticity) = {np.rad2deg(chi):.3f}°")
polarimeter(E, title="After QC")
print(f"Similarity Cos Index = {cos_sim(stokes_fpm, stokes_QC)}")
print(f"------------------Stokes after Second Phase Modulator---------------------")
E = phase_modulator.modulate(E_field=E, V=0)
stokes, polarization = compute_stokes_parameters(E)
S0, S1, S2, S3 = stokes
psi, chi = polarization
print(f"S0 = {S0:.3f}\nS1 = {S1:.3f}\nS2 = {S2:.3f}\nS3 = {S3:.3f}")
print(f"Psi (polarization azimuth) = {np.rad2deg(psi):.3f}°")
print(f"Chi (polarization ellipticity) = {np.rad2deg(chi):.3f}°")
polarimeter(E, title=f"Second PM V = {0:.2f}")
