## opto-sim: A Physics-Based Optical Simulation Framework for QKD Systems

**opto-sim** is a modular headless Python library for simulating real-world optical systems at a backend level. It models components like lasers, LiNbO₃ phase modulators, beam splitters, fiber cables, and APDs, using **numerical methods** and **Stokes vector polarization analysis**. Think of it as a lightweight, scriptable backend akin to what powers tools like HFSS, Lumerical, or SIwave—minus the GUI.

**Work in Progress**  
This project is under active development and not production-ready, but reflects:

- A strong interest in engineering simulation
- Hands-on application of numerical methods to physical systems
- A **backend-first** approach suitable for scripting, automation, and scientific use
- Emphasis on **physics-first modeling**, internal APIs, and extendability.

## Project Overview

This simulator aims to mimic the behavior of physical optical components used in QKD experiments—specifically BB84 protocol setups—through mathematical models. It’s intended as a learning tool and experimental platform for simulating realistic noise conditions, signal degradation, and quantum encoding/decoding via optics.

### Current Features

- **Laser Source**  
  - Solid-state lasers modeled after Nd:YAG and Er:YAG systems.

- **Phase Modulators**  
  - Supports X-cut and Y-cut, Z-propagating configurations.

- **Avalanche Photodiodes (APDs)**  
  - Models include shot noise, thermal noise, and random Gaussian/Poisson noise.

- **Optical Components**  
  - Basic tools like beam splitters, couplers, etc

- **Visualization**  
  - Polarimetry
  - Calculation of **Stokes parameters** for polarization state analysis.
  - Poincare Sphere

- **Fiber Cable**
  - Simulates fiber transmission, optional dispersion modeling

### Mathematical Modeling

All models aim to closely reflect published experimental setups and component behaviors found in academic literature. While accuracy is a priority, simplifications are made when necessary to keep computations tractable.

## Technical Highlights
- Physics-first architecture: Built from experimental papers and known equations
- Modular backend API: Easily compose system-level simulations
- Pure Python + NumPy stack
- Runs in Jupyter/Colab with reproducible outputs and examples
- Output includes:
  - Stokes vector breakdown
  - Cosine similarity across transmission
  - Realistic phase and polarization transitions

## Sample Use Case
### Code
```python
from src.opto_eq import PhaseModulator
from src.lasers import sslaser as laser
from src.viewers.stokes import compute_stokes_parameters
from src.viewers import polarimeter
source = laser.SolidStateLaser(
    wavelength=1550e-9,  # Laser wavelength
    polarization_azimuth=np.pi,  # 45° polarization
    polarization_ellipticity=np.pi/4,
    power_dbm=-5,  # arbitrary power unit
    frequency=5e6
)
E_field = source.get_electric_field(normalize=False, over_period=True)
S0, S1, S2, S3 = compute_stokes_parameters(E_field)
print(f"S0 = {S0:.3f}\nS1 = {S1:.3f}\nS2 = {S2:.3f}\nS3 = {S3:.3f}")
polarimeter(E_field, title=f"Laser Output")
pm = PhaseModulator(crystal_cut='X', modulation='DC')
E_modulated = pm.modulate(E_field, V=3.3)
S0, S1, S2, S3 = compute_stokes_parameters(E_modulated)
print(f"S0 = {S0:.3f}\nS1 = {S1:.3f}\nS2 = {S2:.3f}\nS3 = {S3:.3f}")
polarimeter(E_modulated, title=f"Phase Modulator Output")
```
### Output
```output
S0 = 0.910
S1 = 0.000
S2 = -0.707
S3 = -0.707
```
![Laser Output](https://github.com/azaan-mahmood/opto-sim/blob/main/src/common/images/Laser_Out.png "Laser Output")
```output
S0 = 0.910
S1 = 0.000
S2 = 0.951
S3 = 0.309
```
![PM Output](https://github.com/azaan-mahmood/opto-sim/blob/main/src/common/images/PM_Out.png "Phase Modulator Output")

## Dependencies

- `numpy`
- `scipy`
- `matplotlib`

Install dependencies with:

```bash
pip install -r requirements.txt
```
## Licensing and Usage
🚫 This simulator is intended for demonstration and research purposes only. Please contact me for permission before using in personal research/projects or commercial projects.

## Contact ##
- Feel free to reach out via GitHub issues or messages if you’re curious about the project and want to collaborate.
