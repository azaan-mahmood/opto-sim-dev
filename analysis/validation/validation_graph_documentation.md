# Fiber Channel Validation — Graph Documentation

This document explains every panel in the four validation figures produced by `analysis/validation/`.
Each impairment is validated against its source literature. The physics model (simulation) is compared
directly to the analytic theory derived from the cited reference.

---

## 1. Fiber Attenuation (`val_attenuation--seed42.png`)

**Source literature:** Keiser, G., "Optical Fiber Communications", 5th ed., McGraw-Hill, 2015. Equation 3.6.

**Theory:** Optical power decays exponentially with fiber length:
```
P(L) = P_0 * 10^(-alpha * L / 10)
```
where alpha is the attenuation coefficient in dB/km. For SMF-28 at 1550 nm, alpha = 0.182 dB/km.

**Model:** `apply_attenuation(E, L_km, attenuation_factor)` multiplies the complex envelope field by
`sqrt(10^(-alpha * L / 10))`, preserving the field phase while reducing amplitude.

### Panel A: Power decay (Keiser Eq 3.6)

- **Black line:** Analytic theory — `P(L) = P_0 * 10^(-0.182*L/10)` plotted on a logarithmic y-axis.
- **Red dots:** Measured output power from `apply_attenuation` at each distance.
- **What the graph states:** The red dots lie exactly on the black line across 5 orders of magnitude (10^0 to 10^-4). The simulation reproduces the exponential decay predicted by Keiser Eq 3.6 with no deviation.

### Panel B: Loss vs distance

- **Blue dots:** Measured loss in dB at each distance: `-10*log10(P_out/P_in)`.
- **Black dashed line:** Nominal alpha = 0.182 dB/km (linear in dB).
- **Red solid line:** Linear fit to measured data: alpha_fit = 0.182000 dB/km (R^2 = 1.0).
- **What the graph states:** The dB loss is perfectly linear with distance, confirming that the attenuation coefficient is constant across all distances. The fitted alpha matches the nominal value exactly.

### Panel C: Residual

- **Red line:** Percentage error `|P_measured - P_theory| / P_theory * 100` at each distance.
- **What the graph states:** The maximum relative error is approximately 3.2e-14 %, which is floating-point precision. The simulation is exact to machine accuracy.

### Panel D: SMF-28 attenuation spectrum (Keiser Fig 3.2)

- **Gray line with markers:** Typical SMF-28 attenuation at standard wavelengths (800, 1000, 1200, 1310, 1385, 1450, 1550, 1625 nm). This is the literature reference curve.
- **Red dashed lines:** Markers at 1550 nm (the operating wavelength) and 0.182 dB/km (the nominal attenuation).
- **What the graph states:** The operating point (1550 nm, 0.182 dB/km) is correctly positioned on the standard SMF-28 attenuation spectrum. The attenuation decreases from 2.5 dB/km at 800 nm to 0.182 dB/km at 1550 nm, dominated by Rayleigh scattering at shorter wavelengths.

### Panel E: Alpha consistency across distance

- **Blue dots:** Attenuation coefficient extracted independently at each distance using `-10*log10(P_out/P_in) / L`.
- **Red dashed line:** Nominal alpha = 0.182 dB/km.
- **What the graph states:** The extracted alpha is constant at 0.182 dB/km for all distances from 5 to 200 km. There is no distance-dependent drift, confirming that the model applies Keiser Eq 3.6 uniformly.

### Panel F: Validation summary

- **Table:** Tabulated values of P_measured, P_theory, and error at 0, 20, 50, 100, 150, and 200 km. All errors are at floating-point precision (10^-14 %).

---

## 2. Fiber Birefringence (`val_birefringence--seed42.png`)

**Source literature:**
- Agrawal, G. P., "Fiber-Optic Communication Systems", 5th ed., Wiley, 2021. Equation 4.1.2.
- Ulrich, R. et al., "Bending-induced birefringence in single-mode fibers", Opt. Lett., vol. 5, no. 6, pp. 273-275, 1980.

**Theory:** Birefringence introduces a Jones matrix:
```
J = diag(exp(j * dbeta * L / 2), exp(-j * dbeta * L / 2))
```
where `dbeta = 2 * pi * delta_n / lambda`. The birefringence `delta_n` depends on temperature
and bend radius per the Ulrich model: `delta_n = delta_n_T0 + temp_coeff * (T - 25) + 0.135 * (r/R)^2`.

**Model:** `apply_birefringence(E, L, wavelength, temperature, bend_radius)` constructs the Jones matrix
and applies it to the field via matrix multiplication.

### Panel A: Phase vs length (Agrawal Eq 4.1.2)

- **Black line:** Analytic phase `phi = dbeta * L / 2` computed from the known birefringence `delta_n = 0.87e-5` at 25 C and 1550 nm.
- **Red dots:** Phase of the simulated Jones matrix element `J[0,0]` measured from `apply_birefringence`.
- **What the graph states:** The red dots lie exactly on the black line across the full length range (0 to 0.5 m). The phase increases linearly with length at 35.3 rad/m. The phase error between simulation and theory is 0.00e+00 rad (exact match). The beat length is 178.2 mm.

### Panel B: Temperature sensitivity

- **Black line:** Analytic phase `phi = dbeta(T) * L / 2` at each temperature, where `dbeta(T) = 2*pi*(delta_n_T0 + temp_coeff*(T-25)) / lambda`.
- **Orange dots:** Phase of simulated `J[0,0]` at each temperature, measured from `apply_birefringence`.
- **What the graph states:** The phase decreases linearly with temperature from -20 C to 60 C. The red dots track the black line exactly. The phase error is 0.00e+00 rad. The temperature coefficient is -5e-7 per degree C, which means the birefringence decreases as temperature increases.

### Panel C: Bend-induced birefringence (Ulrich Eq 1)

- **Black line:** Analytic phase `phi = dbeta(R) * L / 2` at each bend radius, where `dbeta(R)` includes the Ulrich bend term `0.135 * (r_fiber/R)^2`.
- **Green dots:** Phase of simulated `J[0,0]` at each bend radius.
- **What the graph states:** The phase increases sharply as the bend radius decreases (tighter bends cause more birefringence). The green dots lie exactly on the black line. The phase error is 0.00e+00 rad. At R = 2 mm, the phase is approximately 2.9 rad; at R = 20 mm, it is approximately 0.2 rad. This confirms the Ulrich model is implemented correctly.

### Panel D: delta_n vs bend radius (Ulrich Eq 1)

- **Green dots:** Birefringence `delta_n` extracted from the simulated Jones matrix using the ratio method (two closely spaced fiber lengths, taking the phase difference).
- **Black dashed line:** Ulrich theory: `delta_n = delta_n_T0 + 0.135 * (r_fiber/R)^2`.
- **What the graph states:** The extracted delta_n values fall exactly on the Ulrich theory curve across the full bend radius range (2 mm to 20 mm). The maximum extraction error is 0.0000%. This validates that the bend-induced birefringence follows the Ulrich model.

### Panel E: Beat length vs wavelength

- **Purple line:** Beat length `L_B = lambda / delta_n` as a function of wavelength, computed from the intrinsic birefringence `delta_n = 0.87e-5`.
- **Gray dashed lines:** Reference markers at 1550 nm and 1310 nm.
- **What the graph states:** The beat length increases linearly with wavelength. At 1550 nm, L_B = 178.16 mm; at 1310 nm, L_B = 150.57 mm. These values are consistent with standard single-mode fiber (long beat length indicates weak birefringence, as expected).

### Panel F: Validation summary

- **Table:** All five tests (power conservation, phase vs length, temperature, bend, wavelength) pass with zero error.

---

## 3. Chromatic Dispersion (`val_cd--seed42.png`)

**Source literature:** Agrawal, G. P., "Fiber-Optic Communication Systems", 5th ed., Wiley, 2021. Section 2.4, Eq 2.4.11, Fig 2.6.

**Theory:** A Gaussian pulse with initial width T0 broadens as it propagates:
```
sigma(z) = T0/sqrt(2) * sqrt(1 + (z/L_D)^2)
```
where `L_D = T0^2 / |beta2|` is the dispersion length, and `beta2 = -D * lambda^2 / (2*pi*c)` is the
group-velocity dispersion parameter. The CD transfer function is `H(omega) = exp(-j * beta2 * omega^2 * L / 2)`.

**Model:** `apply_cd(E, dt, L, wavelength)` applies the CD transfer function in the frequency domain via FFT.

### Panel A: Gaussian pulse broadening (Agrawal Fig 2.6)

- **Four curves:** Normalized pulse intensity `|E(t)|^2 / max` at z/L_D = 0.0 (blue), 0.5 (green), 1.0 (orange), and 2.0 (red).
- **Legend:** Shows the analytic RMS width (sigma) for each distance: 21.2 ps, 23.7 ps, 30.0 ps, 47.4 ps.
- **What the graph states:** The Gaussian pulse broadens as it propagates. At z/L_D = 0, the pulse is at its initial width (sigma = 21.2 ps). At z/L_D = 2.0, the pulse has broadened to sigma = 47.4 ps (2.24x the initial width). The broadening is symmetric and preserves the Gaussian shape, consistent with Agrawal Fig 2.6.

### Panel B: RMS width growth (D = 14.0 ps/(nm*km))

- **Black line:** Analytic curve `sigma(z)/sigma_0 = sqrt(1 + (z/L_D)^2)`.
- **Red dots with line:** Measured RMS width ratio from `apply_cd` at each distance.
- **Red percentage labels:** Relative error at z/L_D = 0.0, 0.5, 1.0, and 2.0 (all 0.0000%).
- **What the graph states:** The measured width ratio follows the analytic curve exactly. The RMS width grows from 1.0 (at z=0) to approximately 2.24 (at z = 2*L_D), matching the theoretical prediction. The maximum error is 3.0e-14 %, which is floating-point precision.

### Panel C: Residual

- **Red line:** Percentage error across the full distance sweep.
- **Black dots:** Test points at z/L_D = 0.0, 0.5, 1.0, and 2.0 (as labeled in the legend).
- **What the graph states:** The maximum relative error is approximately 3.0e-14 %, which is machine precision. The black dots at the test points show zero error. The simulation matches the analytic Gaussian broadening formula exactly.

### Panel D: CD transfer function phase (Agrawal Eq 2.4.11, z = 2L_D)

- **Black line:** Unwrapped analytic phase `phi(f) = -beta2 * (2*pi*f)^2 * L / 2` at z = 2*L_D.
- **Red dotted line:** Unwrapped phase of the simulated transfer function `H_sim = FFT(E_out) / FFT(E_in)`.
- **Blue dashed line:** Reference phase at z = L_D (half the test distance), shown for comparison.
- **What the graph states:** The analytic and simulated phases overlap perfectly, confirming that `apply_cd` implements the correct frequency-domain transfer function. The phase is quadratic in frequency (parabolic shape), which is the signature of chromatic dispersion. At the Nyquist frequency (0.5 THz), the total accumulated phase is approximately 1200 rad. The reference line at z = L_D shows approximately half the phase, as expected from the linear scaling `phi ~ L`.

### Panel E: Material and waveguide dispersion (Hui [2], Keck [3])

- **Blue line:** Material dispersion `D_mat(lambda)` computed from the approximate Sellmeier model.
- **Orange line:** Waveguide dispersion `D_wg = -3.0 ps/(nm*km)` (constant across wavelength).
- **Black dashed line:** Total dispersion `D_total = D_mat + D_wg`.
- **Gray dashed lines:** Reference markers at 1550 nm and 1310 nm.
- **What the graph states:** Material dispersion dominates and decreases with wavelength. Waveguide dispersion is small and approximately constant. At 1550 nm, D_total = D_mat + D_wg = 17.0 + (-3.0) = 14.0 ps/(nm*km), which is the standard SMF-28 dispersion value. The total dispersion crosses zero near 1310 nm, consistent with the zero-dispersion wavelength of standard fiber.

### Panel F: Validation summary

- **Table:** D_total = 14.0 ps/(nm*km), D_material = 17.0, D_waveguide = -3.0, beta2 = -1.786e-26 s^2/m, L_D = 50.40 km, max error = 3.0e-14 %. Gaussian broadening matches analytic to < 0.001 %.

---

## 4. Polarization Mode Dispersion (`val_pmd_dgd--seed42.png`)

**Source literature:** Razavi, B., "Design of Integrated Circuits for Optical Communications", 2nd ed., Wiley, 2012. Section 2.5, Fig 2.11.

**Theory:** PMD causes the differential group delay (DGD) between two principal states of polarization
to follow a Maxwellian distribution:
```
p(DGD) = (2/sqrt(pi)) * (DGD/sigma_m^3) * exp(-DGD^2/sigma_m^2)
```
where `sigma_m = D_pmd * sqrt(L) / sqrt(3)` is the Maxwellian scale parameter, and `D_pmd` is the
PMD coefficient. The RMS DGD equals `D_pmd * sqrt(L)`.

**Model:** `apply_pmd(E, dt, L, pm_dispersion)` samples a Maxwellian DGD and applies a frequency-domain
phase shift between the two polarization components.

### Panel A: DGD distribution (Razavi Fig 2.11)

- **Blue histogram:** Distribution of DGD values from 5000 realizations of `apply_pmd` at L = 100 km.
- **Black line:** Analytic Maxwellian PDF with RMS DGD = 31.6 ps.
- **Gray dashed line:** Expected RMS DGD = 31.6 ps.
- **Blue dotted line:** Measured mean DGD = 29.20 ps.
- **What the graph states:** The simulated DGD histogram matches the analytic Maxwellian distribution. The KS test p-value is 0.45, which means the null hypothesis (that the data comes from a Maxwellian) cannot be rejected. The measured mean (29.20 ps) and RMS (31.74 ps) match the expected values (29.14 ps and 31.62 ps) within statistical fluctuation.

### Panel B: QQ-plot (Maxwellian fit)

- **Blue dots:** Sample quantiles of the simulated DGD vs theoretical Maxwellian quantiles.
- **Black dashed line:** The ideal y = x line (perfect agreement).
- **What the graph states:** The blue dots fall on the y = x line with correlation r = 0.9998. This confirms that the simulated DGD distribution is Maxwellian across the full range (0 to 80 ps).

### Panel C: Histogram residual

- **Red bars:** Difference between the measured histogram density and the analytic Maxwellian PDF at each bin.
- **What the graph states:** The residuals are randomly distributed around zero with no systematic pattern. The maximum residual is approximately 0.003, which is consistent with statistical noise from 5000 realizations.

### Panel D: PMD scaling DGD ~ sqrt(L) (Razavi Sec 2.5)

- **Green dots with line:** RMS DGD measured from `apply_pmd` at fiber lengths of 10, 25, 50, 75, 100, 150, and 200 km.
- **Black dashed line:** Linear fit: DGD_rms = 3.1662 * sqrt(L) (ps).
- **Gray dotted line:** Nominal prediction: DGD_rms = 3.1623 * sqrt(L) (ps).
- **What the graph states:** The RMS DGD scales linearly with sqrt(L), confirming the fundamental PMD scaling law. The fitted PMD coefficient is 3.1662 ps/sqrt(km), which matches the nominal value of 3.1623 ps/sqrt(km) (0.1 ps/sqrt(m)). The R^2 = 0.99993, indicating a near-perfect linear relationship.

### Panel E: PMD coefficient consistency

- **Purple dots with line:** PMD coefficient extracted at each fiber length: `D_pmd = DGD_rms / sqrt(L)`.
- **Gray dashed line:** Nominal PMD coefficient = 0.100 ps/sqrt(km).
- **What the graph states:** The extracted PMD coefficient is approximately constant across all fiber lengths (10 to 200 km), confirming that the model applies PMD consistently regardless of distance. The slight variation is due to statistical fluctuation in the Maxwellian sampling.

### Panel F: Validation summary

- **Table:** Mean DGD = 29.205 ps (expected 29.135), RMS DGD = 31.736 ps (expected 31.623), KS p-value = 0.45 (not rejected), PMD coefficient = 3.1662 ps/sqrt(km) (expected 3.1623), DGD-sqrt(L) R^2 = 0.99993.
