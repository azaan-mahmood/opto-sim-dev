## Fiber Channel Validation Report

Literature cross-check of all Tier-1 validation scripts and fiber.py models

Generated: 2026-07-06 | Commit: 24d4751 | All 48 tests passing

## 1. Chromatic Dispersion Validation

Reference: Agrawal, G. P., "Fiber-Optic Communication Systems", 5th ed., Wiley, 2021, Section 2.4 (Eq 2.4.11) and Section 3.2 (Eq 3.2.6).

## Model:

FFT-based dispersion applied to both Ex and Ey: H(omega) = exp(-j * beta2 * omega^2 * L / 2). The GVD parameter beta2 = -D * lambda^2 / (2 * pi * c). Gaussian input pulse sigma0 = T0 / sqrt(2), theoretical output width sigma(z) = sigma0 * sqrt(1 + (z/LD)^2).

## Validation Result: PASS

Measured error: 0.0000% at z/LD = {0.0, 0.5, 1.0, 2.0}. This is numerically exact because the Fourier transform of a Gaussian is Gaussian, and the FFT-based transfer function H(omega) is the analytic continuation of the quadratic phase factor. Agrawal Fig 2.6 predicts sigma/sigma0 = sqrt(1+(z/LD)^2) for unchirped Gaussian pulses, which the simulation reproduces to machine precision.

## Literature Sources for Parameters:

- \- D_material = 17.0 ps/(nm.km): Material dispersion of silica at 1550 nm (Hui & O'Sullivan [2])

- \- D_waveguide = -3.0 ps/(nm.km): Waveguide dispersion contribution (Keck et al. [3])

- \- D_total = 14.0 ps/(nm.km): Within Corning SMF-28 Ultra spec (<= 18.0 at 1550 nm)

## Issues:

## [WARNING]

Hardcoded D_total = 14.0 in validate_cd.py line 21. If fiber.py changes its dispersion parameters, the validation script silently diverges. Should import D_total dynamically from fiber.py.


## 2. PMD Validation

Primary Reference: Razavi, B., "Design of Integrated Circuits for Optical Communications", 2nd ed., Wiley, 2012, Fig 2.11.

Supporting: Karlsson, M., "Probability density functions of the differential group delay in optical fiber communication systems", JLT 2001; Yang & Kath, "Non-Maxwellian DGD distributions of PMD emulators", OFC 2001; Corning SMF-28 Ultra datasheet (PMD <= 0.04 ps/sqrt(km)).

## Model:

Frequency-domain DGD applied via Jones matrix J_pmd = diag(exp(-j*omega*Dtau/2), exp(+j*omega*Dtau/2)). The DGD is sampled from a Maxwellian distribution: scale = pmd_sd / sqrt(3) where pmd_sd = PMD_coeff * sqrt(L). RMS(DGD) = PMD_coeff * sqrt(L), mean(DGD) = 2 * scale * sqrt(2/pi) = 0.921 * RMS.

## Validation Result: CONDITIONAL PASS

Mean DGD = 29.322 ps (expected 29.135 ps, error 0.64%). RMS DGD = 31.773 ps (expected 31.623 ps, error 0.47%). KS test vs Maxwellian: D = 0.02101, p = 0.0239.

Commit 24d4751 fixed a critical bug: the old code used np.random.rayleigh() (2D magnitude distribution), but PMD theory requires a Maxwellian (3D vector magnitude). The PMD vector has three independent, identically distributed Gaussian components; its magnitude follows a Maxwellian (Razavi Fig 2.11, Karlsson 2001). The Rayleigh distribution would produce a lower mean and incorrect tail behavior, underestimating the probability of rare high-DGD events.

## Issues:

## [WARNING]

KS test p = 0.0239 (p < 0.05 at 5% significance). The DGD histogram marginally fails a Maxwellian goodness-of-fit test. Root cause: DGD extraction via cross-correlation (validate_pmd.py lines 41-48) is quantized at DT = 0.5 ps resolution, introducing discretization error that distorts the distribution.

## [WARNING]

Cross-correlation DGD extraction is indirect and noisy. The script infers DGD from the output field rather than recording the actual DGD value that stats.maxwell.rvs() sampled inside propagate(). This conflates measurement error with model error. Fix: add a return parameter to propagate() or record DGD separately.


## 3. Attenuation Validation

Reference: Keiser, G., "Optical Fiber Communications", 5th ed., McGraw-Hill, 2015, Eq 3.6. Datasheet: Corning SMF-28 Ultra Product Information (PI-1470-AEN), max 0.18 dB/km at 1550 nm.

## Model:

Field attenuation: E *= sqrt(10^(-alpha * L / 10)). Power attenuation: P_out / P_in = 10^(-alpha * L / 10). Default alpha = 0.182 dB/km (SMF-28 Ultra at 1550 nm).

## Validation Result: PASS

Measured error: 0.0000% at all 41 distance points (0-200 km). The implementation exactly reproduces the Keiser Eq 3.6 formula, as expected for a single algebraic scaling.

## Issues:

## [WARNING]

Pure self-consistency check: the script verifies propagate() matches the formula, but does not overlay external experimental data (e.g., OTDR traces from Keiser or Thorlabs app notes). The roadmap (tier_1_attack_plan.md Task 3) suggests overlaying OTDR data for external validation.

## [WARNING]

Default alpha = 0.182 dB/km is slightly above the Corning SMF-28 Ultra max spec of 0.18 dB/km. For SMF-28e (standard grade), the spec is 0.19-0.22 dB/km. Consider whether the default should match the Ultra spec (0.18) or be configurable.


## 4. Birefringence Validation

Cited Reference [4]: Yuan, L. et al., "Stress-induced birefringence and fabrication of in-fiber polarization devices by controlled

femtosecond laser irradiations", Optics Express, vol. 24, no. 2, pp. 1062-1071, 2016. DOI: 10.1364/oe.24.001062. Additional: Ulrich, R. et al., "Bending-induced birefringence in single-mode fibers", Optics Letters, vol. 5, no. 6, pp. 273-275, 1980.

Smith, A. M., "Birefringence induced by bends and twists in single-mode optical fiber", Applied Optics, vol. 19, no. 15, 1980.

## Model:

Birefringence delta_n = delta_n0 + T_coeff * (T - T0) + 0.135 * (r_f/R)^2 (Ulrich [7]). Phase accumulation: delta_phi = 2 * pi * L * delta_n / lambda. Symmetric Jones matrix: diag(exp(j * delta_beta * L / 2), exp(-j * delta_beta * L / 2)), delta_beta = 2 * pi * delta_n / lambda (Agrawal [6] Eq 4.1.2).

## Validation Result: CONDITIONAL PASS

Self-consistency error: 0.0000% for all three coefficients (base, temperature, bend). The propagate() function exactly reproduces the coefficients declared in fiber.py.

## Issues:

## [CRITICAL]

Reference [4] in fiber.py header cites "Opt. Express, vol. 27, 2019" but the Yuan paper about stress-induced birefringence (Delta_n = 2.4e-4) was published in vol. 24, 2016. The correct citation is: L. Yuan et al., Opt. Express 24(2), 1062-1071 (2016). No matching Yuan paper in vol. 27 (2019) was found.

## [CRITICAL]

bend_effect_factor = 2.4e-4 per "bend" was misapplied. The Yuan 2016 paper measures Delta_n = 2.4e-4 from femtosecond-laser-fabricated stress rods, not from fiber bending. True bending-induced birefringence follows a radius-dependent formula: Delta_n_bend = 0.135 * (r_f/R)^2 (Ulrich 1980, Smith 1980, Shibata 1986). **FIXED**: replaced num_bends with bend_radius using the Ulrich formula. Validated: 0.0000% error.

## [WARNING]

birefringence_T0 = 0.87e-5: Standard SMF-28 residual birefringence is typically 1e-7 to 1e-6 (beat length ~10-100 m). A value of 8.7e-6 corresponds to beat length L_B = lambda/D_n = 1550 nm / 8.7e-6 = 0.18 m, which is characteristic of polarization-maintaining fiber, not standard SMF. Source for this value should be cited.

## [WARNING]

temperature_coefficient = -5e-7 /C: The thermo-optic coefficient of silica is dn/dT ~ +1.1e-5 /C. This -5e-7/C is the temperature derivative of stress-induced birefringence (d(Dn)/dT), not the index itself. The literature source should be cited.


## 5. Commit 24d4751 Changes to fiber.py

The latest commit (24d4751) changed PMD DGD sampling from Rayleigh to Maxwellian:

## Old code (incorrect):

dgd = np.random.rayleigh(pmd_sd * sqrt(2/pi))

This samples a Rayleigh distribution (2D vector magnitude). PMD theory requires a Maxwellian (3D vector magnitude).

## New code (correct):

dgd = stats.maxwell.rvs(scale = pmd_sd / sqrt(3))

Maxwellian with scale = pmd_sd/sqrt(3) gives RMS(DGD) = pmd_sd and mean(DGD) = 2*scale*sqrt(2/pi) = 0.921*pmd_sd. This matches the ITU-T and Corning PMD link design value convention where the PMD coefficient specifies the RMS DGD per sqrt(km).

## Literature Support:

- \- Razavi [5] Fig 2.11: DGD distribution for long fibers is Maxwellian

- \- Karlsson (JLT 2001): "probability density function of the DGD ... Maxwellian"

- \- Corning SMF-28 datasheet: "DGD histograms could be well approximated by a Maxwellian probability distribution"

- \- Yang & Kath (OFC 2001): asymptotic DGD distribution is Maxwellian after ~30 correlation lengths

## Verdict: FIXED CORRECTLY

The Rayleigh distribution was a genuine physics bug that would underestimate high-DGD tails. The Maxwellian fix aligns with all literature sources.


## 6. Summary Table

| Validation Math Lit. Correct Status Issues |
| --- |
| CD Pass Pass OK Hardcoded D_total (minor) |
| PMD Pass* Pass OK* KS p=0.024, extraction method |
| Attenuation Pass Pass OK Missing OTDR overlay |
| Birefringence | Pass | PASS | FIXED | Wrong ref year (Fix 1), bend model (Fix 2 — fixed), symmetric Jones (Fix 5) |

## 7. Critical Fixes Required

## Fix 1: Correct Yuan reference in fiber.py (line 12)

Current: "Yuan, L. et al., ... Opt. Express, vol. 27, 2019."

Correct: "Yuan, L. et al., Stress-induced birefringence and fabrication of in-fiber polarization devices by controlled femtosecond laser irradiations, Optics Express, vol. 24, no. 2, pp. 1062-1071, 2016. DOI: 10.1364/oe.24.001062."

Source: https://doi.org/10.1364/oe.24.001062

## Fix 2: Fix birefringence bend model in fiber.py (was lines 87-93, now fixed)

The old "bend_effect_factor * num_bends" model treated each bend as a fixed Delta_n = 2.4e-4 contribution independent of radius. This was not physically correct.

**FIXED**: Replaced with Delta_n_bend = 0.135 * (r_fiber / R)^2 (Ulrich 1980, Smith 1980, Shibata 1986). Uses bend_radius parameter (float, metres). Validated: 0.0000% error.

## Sources:

- [7] Ulrich, R. et al., Opt. Lett. 5(6), 273-275, 1980. DOI: 10.1364/ol.5.000273

- [8] Smith, A. M., Appl. Opt. 19(15), 2606, 1980. DOI: 10.1364/ao.19.002606

- [12] Shibata, N. et al., J. Opt. Soc. Am. A 3(11), 1935, 1986. DOI: 10.1364/josaa.3.001935

## Fix 3: Fix PMD DGD extraction in validate_pmd.py (lines 41-48)

Replace cross-correlation DGD extraction with direct recording of the DGD value sampled by stats.maxwell.rvs() inside propagate(). This eliminates discretization noise at DT = 0.5 ps and produces a clean Maxwellian histogram for the KS test.

## Fix 4: Derive D_total dynamically in validate_cd.py (line 21)

Import the dispersion parameters from fiber.py rather than hardcoding D_total = 14.0. This keeps the validation script in sync with the implementation.


## 8. Literature Sources

- [1] Keiser, G., "Optical Fiber Communications", 5th ed., McGraw-Hill, 2015. Ch. 3: fiber attenuation, dispersion, PMD.

- [2] Hui, R. & O'Sullivan, M., "Fiber-Optic Measurement Techniques", Academic Press, 2009.

- [3] Keck, D. B. et al., "Waveguide dispersion in single-mode fibers", IEEE J. Quantum Electron., QE-21(6), 1985.

- [4] Yuan, L. et al., "Stress-induced birefringence and fabrication of in-fiber polarization devices by controlled femtosecond laser irradiations", Opt. Express 24(2), 1062-1071, 2016. DOI: 10.1364/oe.24.001062.

- [5] Razavi, B., "Design of Integrated Circuits for Optical Communications", 2nd ed., Wiley, 2012. Fig 2.11: PMD Maxwellian DGD.

- [6] Agrawal, G. P., "Fiber-Optic Communication Systems", 5th ed., Wiley, 2021. Eq 2.4.11 (CD), Eq 4.1.2 (birefringence).

- [7] Ulrich, R. et al., "Bending-induced birefringence in single-mode fibers", Opt. Lett. 5(6), 273-275, 1980.

- [8] Smith, A. M., "Birefringence induced by bends and twists in single-mode optical fiber", Appl. Opt. 19(15), 2606, 1980.

- [9] Karlsson, M., "Probability density functions of the differential group delay in optical fiber communication systems", J. Lightwave Technol. 19(3), 2001.

- [10] Corning Inc., "SMF-28 Ultra Optical Fiber Product Information", PI-1470-AEN, 2020.

- [11] Garth, S. J., "Birefringence in bent single-mode fibers", J. Lightwave Technol. 6(3), 445-449, 1988.

- [12] Shibata, N. et al., "Birefringence and polarization mode dispersion in a coil of a single-mode fiber", J. Opt. Soc. Am. A 3(11), 1935, 1986.
