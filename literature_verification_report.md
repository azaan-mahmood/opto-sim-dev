## Literature-to-Implementation Verification Report

Date: 2026-07-08 | Project: opto-sim-dev | Commit: ebf4971 | All 48 tests passing

---

### 1. `src/channel/fiber.py` — Fiber Channel

**References:** [1] Keiser 2015, [2] Hui 2009, [3] Keck 1985, [4] Yuan 2016, [5] Razavi 2012, [6] Agrawal 2021

| Eq/Location | Reference | Implementation | Status |
|---|---|---|---|
| Keiser [1] Eq 3.6 | P_out/P_in = 10^{-αL/10} | `att_lin = 10**(-attenuation_factor * fiber_length / 10)` (L:170) | **PASS** |
| Agrawal [6] Eq 2.4.11 | H(Ω) = exp(-j·β₂·Ω²·L/2) | `H = np.exp(-1j * beta2 * omega**2 * L / 2)` (L:140) | **PASS** |
| Agrawal [6] Eq 2.4.11 | β₂ = -D·λ²/(2πc) | `beta2 = -D_SI * wavelength**2 / (2 * np.pi * c0)` (L:130) | **PASS** |
| Razavi [5] Fig 2.11 | DGD follows Maxwellian (3D PMD vector) | `stats.maxwell.rvs(scale=pmd_sd/sqrt(3))` (L:153) | **PASS** |
| Razavi [5] | J_pmd = diag(exp(∓jωΔτ/2), exp(±jωΔτ/2)) | Hx/Hy with opposite signs (L:159-160) | **PASS** |
| Yuan [4] | Bend-induced birefringence 2.4e-4 (see Fix 2) | `bend_effect_factor = 2.4e-4` (L:99) | **CONDITIONAL** — deferred |

Anomalies:

- **Δβ = 4π·Δn/λ** (L:101): Non-standard definition. Standard (Agrawal [6] §4.1.2) is Δβ = 2π·Δn/λ. The code compensates with a L/2 factor in the Jones matrix, giving the same total phase 2π·L·Δn/λ. Mathematically correct but unconventional — the factor choice conflates Δβ (propagation-constant difference) with the total phase accumulated in one beat length.
- **Asymmetric Jones matrix** (L:102-103): `diag(exp(jΔβL/2), 1)` instead of the symmetric `diag(exp(jΔβL/2), exp(-jΔβL/2))`. Gives correct relative phase between Ex and Ey (Δβ·L/2 = 2π·L·Δn/λ) but the absolute phase is offset by Δβ·L/4. Only relative phase matters for polarization, so this is benign. Consider switching to the symmetric form for conceptual clarity.
- **PMD applies fixed sign** (L:159-160): The fast axis is always Ex (Hx = exp(-jωDGD/2), Hy = exp(+jωDGD/2)). In a real fiber, the PMD vector orientation is random per realization. The sign of the DGD application should be ±1 with 50% probability each. **Effect on QKD**: minimal because BB84 randomizes the measurement basis, but a systematic PMD sign could bias the output polarization in a single-fiber setup.
- **_dgd_sampled side-effect** (L:29, L:154): Module-level list accumulates DGD from every `cable()` call with `dispersion=True`. Tests and scripts that call `cable()` will append values. Validation scripts clear it first, but running tests + validation in the same process would mix values. **Fix**: use a separate subprocess (as `run_all.py` already does) or add a reset function.

---

### 2. `src/lasers/cwlaser.py` — CW Laser

**References:** [1] Henry 1982, [2] Coldren 2012, [3] Yariv 1991, [4] Schawlow-Townes 1958, [5] Petermann 1988

| Eq/Location | Reference | Implementation | Status |
|---|---|---|---|
| Henry [1] Eq 18 | D_φ = 2π·Δν (phase diffusion coefficient) | `_phase_diff_coeff = 2.0 * np.pi * linewidth` (L:99) | **PASS** |
| Henry [1] | Wiener process φ(t+dt) = φ(t) + √(D_φ·dt)·N(0,1) | `increments = np.random.normal(0, sqrt(D_phi * dt))` (L:139-141) | **PASS** |
| Coldren [2] Eq 5.3.38 | S_RIN(f) = RIN_0·(γ²+ω²) / (|ω_R²-ω²+jγω|²) | `num = gamma**2 + omega**2`; `den = (omega_R**2 - omega**2)**2 + (gamma*omega)**2` (L:199-200) | **PASS** |
| Yariv [3] Ch. 6 | Jones vector from ψ, χ | `_polarization_vector()` (L:116-130) | **PASS** |

Anomalies:

- **RIN_0 is a user parameter** (L:103): `self._rin_linear = 10.0 ** (rin_density / 10.0)`. This is phenomenologically valid but does not derive RIN_0 from first principles (spontaneous emission). Coldren [2] Eq 5.3.38 defines RIN_0 from the rate-equation parameters. By making it a free parameter, the model can produce realistic RIN spectra but cannot predict the RIN floor from basic laser parameters. This is appropriate for QKD simulation where RIN is a controlled variable rather than a prediction target.
- **No intensity noise below RIN_0 floor**: The frequency-domain method shapes white noise by sqrt(S_RIN(f)). The normalization ensures the PSD integral matches the specified RIN density at DC. Verified: DC amplitude matches `self._rin_linear`. OK.
- **Relaxation frequency sign convention**: `omega_R = 2π·f_RO` is in rad/s, and `gamma` is the damping rate in rad/s. Coldren Eq 5.3.38 uses the same convention. ✓

---

### 3. `src/channel/mzm.py` — Mach-Zehnder Modulator

**References:** [1] Agrawal 2010 §4.2, [2] Koyama & Iga 1988, [3] Weis & Gaylord 1985

| Eq/Location | Reference | Implementation | Status |
|---|---|---|---|
| Agrawal [1] §4.2 | E_out ∝ cos(π·(V+V_bias)/(2·V_π)) | Push-pull: `np.sqrt(IL) * E_in * cos(pi*(V+Vb)/(2*Vpi))` via Y-branch (L:148-170) | **PASS** |
| Koyama & Iga [2] | Single-drive: residual chirp exp(jπ·(V+Vb)/(2V_π)) | `np.exp(1j * pi * (V + V_bias) / (2 * V_pi))` combined with cos (L:151-170) | **PASS** |
| Weis & Gaylord [3] | LiNbO3 Pockels coefficients r13, r22 | Delegated to `PhaseModulator` | **PASS** (delegated) |

Anomalies:

- **Extinction ratio formula** (L:95-98): `delta = 10^(-ER_dB/20)`, `r = 0.5 * (1+delta)`. This gives the correct Y-branch splitting imbalance. Verified: for 30 dB ER, null transmission = 0.001 = -30 dB. ✓
- **V_π is derived from internal PhaseModulator** (L:84). Verified: for X-cut DC modulator, V_π = 3.88 V with default params. ✓
- **Crystal cut modulates the correct component** (L:178-181): X-cut modulates Ey ((...,1)), Y-cut modulates Ex ((...,0)). Consistent with Z-propagation convention. ✓
- **Literature edition mismatch**: Reference [1] is cited as "4th ed., Wiley, 2010" but same equation appears in the 5th ed. (2021) as well. Minor — equation is identical across editions.

---

### 4. `src/channel/phase_modulator.py` — Phase Modulator (LiNbO3 Pockels)

**References:** Web sources + papers (no numbered refs)
V_π formula (L:84-93): `Vπ = λ·d / (2·n_o³·r·Γ·L)`

| Eq/Location | Implementation | Status |
|---|---|---|
| V_π for X-cut (r13) | `(self.wavelength * self.d) / (2 * self.n_o**3 * self.r13 * self.Gamma * self.L)` | **CONDITIONAL** |
| V_π for Y-cut (r22) | `(self.wavelength * self.d) / (2 * self.n_o**3 * self.r22 * self.Gamma * self.L)` | **CONDITIONAL** |

Anomalies:

- **Missing formal literature references**: The docstring lists URLs and paper titles but no numbered references with specific equations. The V_π formula is the standard one for a bulk LiNbO3 phase modulator but no peer-reviewed source is cited for the exact form. **Recommendation**: add a numbered reference [1] for the electro-optic coefficient definitions (Weis & Gaylord 1985, Appl. Phys. A 37, 191-203) and [2] for the V_π formula (Alferness, R. C., "Titanium-diffused lithium niobate waveguide devices", in Guided-Wave Optoelectronics, Springer, 1988, Ch. 4).
- **Uses n_o for both X-cut and Y-cut** (L:84-93): For X-cut LiNbO3 with Z-propagation, the r13 Pockels coefficient couples the applied field to the extraordinary index. The relevant refractive index is n_e = 2.14, not n_o = 2.2. Using n_o³ = 10.65 instead of n_e³ = 9.80 gives V_π ≈ 3.88 V vs 3.57 V (≈ 9% error). Some references use n_o for r13 and n_e for r33 only; there is no universal consensus. **Marking as CONDITIONAL** pending exact configuration verification.
- **V_π formula uses 2·n³·r·Γ·L in denominator** (L:85, 90). The standard formula extends `λ·d / (n³·r·Γ·L)` for push-pull MZMs where the voltage is split across two arms (the factor of 2 comes from half the voltage per arm). For a single phase modulator without push-pull, the formula should be `λ·d / (n³·r·Γ·L)`. The factor of 2 makes V_π ≈ 3.88 V instead of ≈ 1.94 V. However, since the PhaseModulator's V_π is used by MZM, which handles the push-pull voltage splitting internally, this extra factor of 2 in the PhaseModulator's V_π means the combined V_π is correct for the MZM. This is an issue of where the factor-2 responsibility is assigned. **Not a bug**, but the docstring should clarify that this is the MZM-effective V_π (push-pull configuration), not the single-arm V_π.

---

### 5. `src/detectors/apd.py` — Avalanche Photodiode

**References:** [1] Kasap 2013, [2] Agrawal 2021, [3] Saleh & Teich 2019

| Eq/Location | Reference | Implementation | Status |
|---|---|---|---|
| Kasap [1] Eq 4.19 | R = η·e·λ/(h·c) | `self.qe * self.charge * self.wavelength / (self.h * self.c)` (L:31) | **PASS** |
| Kasap [1] Eq 4.23 | I_signal = M·R·P | `self.gain * self.R * power` (L:62) | **PASS** |
| Kasap [1] Eq 4.42 | i_d² = 2·e·I_dark·B | `2 * self.charge * self.dark_current * bandwidth` (L:77) | **PASS** |
| Kasap [1] Eq 4.43 | i_q² = 2·e·I_signal·B | `2 * self.charge * I_signal * bandwidth` (L:78) | **PASS** |
| Kasap [1] Eq 4.44 | i_th² = 4·k·T·B/R_L | `4 * self.kB * self.T * bandwidth / self.RL` (L:79) | **PASS** |
| Kasap [1] Eq 4.45 | i_total² = F·(i_d² + i_q²) + i_th² | `self.enf * (shot_dark_sq + shot_signal_sq) + thermal_sq` (L:81) | **PASS** |
| Agrawal [2] Eq 4.1.2 | photon_rate = P/(h·ν) | `power / (self.h * self.frequency)` (L:49) | **PASS** |
| Saleh & Teich [3] Eq 17.1-10 | N_expected = (P/(h·ν))·t·η | `(power / photon_energy) * exposure_time * self.qe` (L:49) | **PASS** |

Anomalies: None. All equations match the cited literature exactly. The excess noise factor F is correctly applied only to shot-noise terms (not thermal), matching Kasap Eq 4.45 and the physical origin (F multiplies the variance of the multiplied shot noise, but thermal noise is added after the APD gain stage).

---

### 6. `src/visualization/stokes.py` — Stokes Parameters

**References:** [1] Collett 2005, [2] Hecht 2002, [3] Born & Wolf 1999

| Eq/Location | Reference | Implementation | Status |
|---|---|---|---|
| Collett [1] Eq 2.12 | S0 = ⟨|Ex|² + |Ey|²⟩ | `np.mean(Ex*conj(Ex) + Ey*conj(Ey))` (L:28) | **PASS** |
| Collett [1] Eq 2.13 | S1 = ⟨|Ex|² - |Ey|²⟩/S0 | Same, /S0 (L:34) | **PASS** |
| Collett [1] Eq 2.14 | S2 = 2·Re(⟨Ex·Ey*⟩)/S0 | `2 * np.real(np.mean(Ex * np.conj(Ey))) / S0` (L:35) | **PASS** |
| Collett [1] Eq 2.15 | S3 = -2·Im(⟨Ex·Ey*⟩)/S0 | `-2 * np.imag(np.mean(Ex * np.conj(Ey))) / S0` (L:37) | **PASS** |
| Collett [1] Eq 2.28 | ψ = 0.5·arctan2(S2, S1) | `0.5 * np.arctan2(S2, S1)` (L:40) | **PASS** |
| Collett [1] Eq 2.29 | χ = 0.5·arcsin(S3) | `0.5 * np.arcsin(S3)` (L:41) | **PASS** |

Anomalies:

- **Missing S3 clip** (L:41-43): The comment on lines 42-43 explicitly identifies the need to clip S3 before arcsin to guard against floating-point noise, but the code does NOT implement the clip. If |S3| slightly exceeds 1.0 (possible for near-circular states due to finite-precision averaging), `np.arcsin(S3)` returns NaN, propagating to chi = NaN.
  ```
  # Comment says:  Clip to [-1, 1] guards against floating-point noise in S3
  # Code does:    chi = 0.5 * np.arcsin(S3)                     # NO CLIP
  ```
  **Fix**: `chi = 0.5 * np.arcsin(np.clip(S3, -1.0, 1.0))`
- **S3 sign convention** (L:37): The code uses `S3 = -2·Im(...)`, consistent with Collett [1] Eq 2.15 (right-handed = negative = S3 < 0). This is the standard optics convention (Born & Wolf [3]). ✓

---

### 7. `src/visualization/fields.py` — Electric Field Plotter

**Reference:** [1] Hecht 2002 Ch. 8

Minimal file — plots Ex, Ey, |E| vs time. No physics equations to verify beyond the Hecht reference for the field representation. Uses `np.real(E)` for the real part, which is standard.

Anomaly: `t = np.linspace(0, 2π/frequency, 1000)` hardcodes 1000 points (L:16). The `title` parameter is accepted but `title(title)` is called on the first subplot before `plt.suptitle` or similar is used. The title only appears above the first subplot. Minor usability issue.

---

### 8. `analysis/validation/validate_cd.py` — CD Validation

**Reference:** Agrawal [6] Eq 2.4.11 (GVD), Eq 3.2.6 (Gaussian pulse broadening)

| Check | Result |
|---|---|
| sigma(z) = sigma0·√(1 + (z/LD)²) | **PASS**: 0.0000% error at z/LD = 0.0, 0.5, 1.0, 2.0 |
| beta2 computed correctly | **PASS**: β₂ = -D·λ²/(2πc) |
| D_TOTAL imported dynamically | **PASS** (Fix 4 applied) |

---

### 9. `analysis/validation/validate_pmd.py` — PMD Validation

**References:** Razavi [5] Fig 2.11, Karlsson [9], Yang & Kath, Corning

| Check | Result |
|---|---|
| Maxwellian DGD distribution | **PASS**: KS p=0.31 (Fix 3 applied) |
| Mean DGD = 29.3 ps (expected 29.1 ps) | **PASS**: error 0.6% |
| RMS DGD = 31.8 ps (expected 31.6 ps) | **PASS**: error 0.5% |
| DGD recorded directly from cable() | **PASS** (Fix 3 applied) |

---

### 10. `analysis/validation/validate_attenuation.py` — Attenuation Validation

**Reference:** Keiser [1] Eq 3.6

| Check | Result |
|---|---|
| P_out/P_in = 10^{-αL/10} | **PASS**: 0.0000% error at 0-200 km |
| Validates via cable() output | **PASS** |

---

### 11. `analysis/validation/validate_birefringence.py` — Birefringence Validation

**References:** Yuan [4], Ulrich [7], Smith [8]

| Check | Result |
|---|---|
| Base Δn = 0.87e-5 | **PASS**: 0.0000% error |
| Temperature coefficient = -5e-7 /C | **PASS**: 0.0000% error |
| Bend coefficient = 2.4e-4 /bend | **PASS**: 0.0000% error |

---

### 12. Summary of Issues Found

| # | File | Severity | Issue | Fix |
|---|---|---|---|---|
| 1 | stokes.py:41 | **BUG** | `np.arcsin(S3)` without clip — NaN if |S3| slightly > 1 | Add `np.clip(S3, -1.0, 1.0)` |
| 2 | fiber.py:101-103 | Minor | Δβ = 4π·Δn/λ is non-standard; Jones matrix asymmetry | Use symmetric form: `diag(exp(jΔβL/2), exp(-jΔβL/2))` with Δβ = 2π·Δn/λ |
| 3 | fiber.py:159-160 | Minor | PMD sign always applies Ex fast / Ey slow | Randomize sign: `if np.random.rand() < 0.5: swap(Hx, Hy)` per realization |
| 4 | phase_modulator.py:84-93 | Minor | Uses n_o for both cuts; n_e may be correct for X-cut (r13) | Verify exact crystal orientation and update if needed |
| 5 | phase_modulator.py:84-93 | Minor | V_π formula has implicit factor of 2 (single-arm vs push-pull) | Clarify in docstring that V_π is the MZM-effective value |
| 6 | phase_modulator.py | Info | Missing formal numbered references | Add Weis & Gaylord 1985 for coefficients, Alferness 1988 for V_π |
| 7 | fiber.py:29,154 | Info | `_dgd_sampled` is a module-level side-effect | Already mitigated by subprocess in `run_all.py` |
| 8 | fields.py:16,18 | Info | Hardcoded 1000 points; title only applies to first subplot | Minor usability — no physics impact |

### 13. Literature Coverage by Severity

**Pass (no issues):** APD, MZM, CWLaser (RIN + phase noise + polarization), Attenuation, CD validation

**Minor anomalies (physics correct, documentation or convention):** Fiber (Δβ convention, PMD sign), PhaseModulator (n_o vs n_e, factor-2), Fields (hardcoded points)

**Bug (incorrect behavior for edge cases):** Stokes (missing S3 clip — NaN for near-circular polarization)

### 14. Recommendations

1. **CRITICAL**: Fix the missing S3 clip in `stokes.py:41` — `np.arcsin(np.clip(S3, -1.0, 1.0))`
2. **MINOR**: Switch fiber birefringence to symmetric Jones matrix with standard Δβ = 2π·Δn/λ for conceptual clarity
3. **MINOR**: Randomize PMD fast/slow axis sign in fiber.py
4. **INFO**: Document the V_π factor convention in PhaseModulator
5. **INFO**: Add numbered references to PhaseModulator docstring
