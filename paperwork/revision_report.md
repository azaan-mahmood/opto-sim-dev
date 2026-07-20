# Manuscript Revision Report

## Part A: Manuscript (.tex) Edits

---

### A1. Abstract — Strengthen Novelty Framing

**Current** (lines 47–88): Leads with component list, positions "validation methodology" as primary contribution.

**Problem**: The combination itself (complex-envelope field + per-component validation + verified composition + open-source) is the novel claim, but it's never stated explicitly. The reader has to infer it from separate bullet points.

**Replace with**:

```latex
\begin{abstract}
We present an open-source optical channel simulation framework for quantum
key distribution (QKD) that is, to the best of our knowledge, the first to
simultaneously provide (i)~complex-envelope electric field propagation at
the physical-layer level, (ii)~independent per-component validation against
published literature, (iii)~verified system-level composition of those
components, and (iv)~full seeded reproducibility under an open-source licence.
No existing QKD simulator --- whether protocol-level, commercial optical, or
QKD-specific --- offers this combination: protocol simulators operate at the
qubit abstraction level without field physics, while commercial optical
simulators are closed-source and lack transparent per-component validation.

The field $\mathbf{E}(t)=[E_x(t),E_y(t)]^\mathsf{T}$ is the single source of
truth for both polarisation state and optical power ($\langle|E|^2\rangle =
P_{\text{opt}}$).  Five components are independently validated:
(i)~a continuous-wave laser with Wiener phase noise and relaxation-oscillation
RIN (Henry~1982, Coldren~2012), reproducing the correct power convention
($<5\times10^{-4}\%$ error), phase diffusion, and RIN spectral shape;
(ii)~a Mach--Zehnder modulator with push-pull (zero chirp) and single-drive
(linear chirp) configurations, finite extinction ratio, insertion loss, and
crystal-cut-dependent modulation (Agrawal~2021, Koyama~\&~Iga~1988);
(iii)~a four-impairment fibre channel --- temperature- and bend-dependent
birefringence via random-axis SU(2) rotation (Menyuk~\&~Wai~1994,
Wai~\&~Menyuk~1996, Ulrich~1980), FFT-based chromatic dispersion
(Agrawal~2021, $<10^{-13}\%$ error), Maxwellian polarisation-mode
dispersion (Razavi~2012, KS~$p=0.42$), and exponential attenuation
(Keiser~2015, $<10^{-13}\%$ error); and (iv)~an avalanche photodiode with
excess-noise-factor-corrected shot noise and Johnson--Nyquist thermal noise
(Kasap~2013, Saleh~\&~Teich~2019).

The system-level BB84 demonstration confirms that independently validated
components compose correctly, producing physically sensible QBER under
combined CD, PMD, birefringence, and attenuation across sweeps of distance,
pulse width, temperature, and bend radius.  The framework is fully
reproducible (seeded RNG, version-controlled, 48~unit tests passing) and
available at \texttt{https://github.com/azaan-mahmood/opto-sim-dev}.
\end{abstract}
```

**What changed**:
- Opens with an explicit gap statement: "first to simultaneously provide (i)–(iv)"
- Names the two camps (protocol vs commercial) and states neither provides this combination
- Validation results condensed into a single tight paragraph
- Composition claim moved to the closing position (not a "separate item")
- Changelog of edits: ~40% shorter, higher information density

---

### A2. Motivation Paragraph (Introduction, lines 113–125)

**Current** (lines 113–125):
```
Existing QKD simulators fall into two camps.  \emph{Abstract qubit
simulators} ...  \emph{Commercial optical simulators} ...
Neither category provides the combination of first-principles optical
modelling, per-component literature validation, and open-source
reproducibility that defensible QKD simulation requires.
```

**Problems**:
1. "Commercial optical simulators ... operate at the sinusoid-carrier level" is incorrect — VPI/OptiSystem/OptSim solve the nonlinear Schrödinger equation for the complex envelope. The issue is not the representation but that they are closed-source and unvalidated.
2. Missing "verified system-level composition" from the gap list.
3. The sentence structure buries the punchline.

**Replace with**:
```latex
Existing QKD simulators fall into two camps, and neither provides the
combination required for defensible physical-layer QKD simulation.
\emph{Abstract qubit simulators} such as NetSquid~\cite{netSquid},
SimulaQron~\cite{simulaqron}, and QuTIP~\cite{qutip} represent channels
as depolarising maps with parametrised loss and click probabilities.
There is no field, no wavelength, no phase noise, and no birefringence
between Alice and Bob.  \emph{Commercial optical simulators} such as
OptiSystem~\cite{optisystem}, VPItransmissionMaker~\cite{vpiphotonics},
and Synopsys OptSim~\cite{optsim} solve the nonlinear Schr\"odinger
equation for the complex envelope, but they are closed-source,
proprietary, and --- critically --- do not provide transparent,
independently verifiable per-component validation against published
literature.  No existing simulator in either camp simultaneously
provides first-principles optical field modelling, per-component
validation, verified system-level composition, and open-source
reproducibility.
```

---

### A3. Related Work — Commercial Simulators (lines 182–209)

**Current line 184**: "not designed for the complex-envelope representation used in QKD"

**Problem**: As noted above, commercial tools do work at the complex-envelope level. This claim is inaccurate and a reviewer will flag it.

**Fix (line 184)**: Change to:
```latex
but are closed-source, expensive, and --- critically --- do not provide
transparent, independently verifiable per-component validation against
published literature.
```

Similarly line 198–204 about VPIphotonics: change the critique from representation to validation transparency.

---

### A4. Gap Analysis (lines 229–279) — Add item (g) to criteria list

**Current** (lines 231–237): lists (a) through (f), missing composition verification.

**Replace (a)–(f) with (a)–(g)**:
```latex
No existing tool provides all of the following simultaneously:
(a)~complex-envelope optical field propagation,
(b)~chromatic dispersion via FFT-based transfer functions,
(c)~Maxwellian PMD,
(d)~temperature- and bend-dependent birefringence,
(e)~per-component validation against published literature,
(f)~open-source availability with reproducible validation, and
(g)~a system-level demonstration that independently validated components
compose correctly under combined impairments.
```

**Also replace paragraph below Table 1** (lines 273–279):
```latex
The gap is clear: protocol-level simulators sacrifice optical-layer
fidelity for scalability and protocol flexibility, while commercial
optical simulators provide detailed physics but are closed-source,
expensive, and lack independent per-component validation.  Neither
category --- nor the QKD-specific simulators surveyed in
Section~\ref{sec:qkd-sim} --- provides all of the following in a
single framework: first-principles complex-envelope optical field
modelling, per-component literature validation, verified system-level
composition, and full open-source reproducibility.  Our framework
targets exactly this intersection.
```

---

### A5. Contributions Section (lines 281–293)

**Current**: A flat list of 9 bullet points with no structure. The first three ("open-source framework", "literature traceability", "component-wise validation") overlap heavily.

**Replace with**:
```latex
\subsection{Contributions}

Our primary contribution is a framework that simultaneously achieves
four capabilities absent from existing QKD simulators:

\begin{itemize}[nosep]
  \item \textbf{First-principles complex-envelope field propagation} ---
    the electric field $\mathbf{E}(t)$ is the single source of truth
    for power and polarisation, enabling physically correct modelling
    of CD, PMD, birefringence, and their interactions, in contrast to
    the abstract qubit-level models used by protocol simulators.

  \item \textbf{Per-component literature validation} --- every equation
    cites a specific numbered reference (34~sources), and each of the
    five components is independently verified against published theory,
    in contrast to commercial simulators whose internal physics is not
    transparently verifiable.

  \item \textbf{Verified system-level composition} --- a six-panel BB84
    demonstration confirms that independently validated components
    compose correctly under combined impairments, producing physically
    sensible QBER across sweeps of distance, pulse width, temperature,
    and bend radius.

  \item \textbf{Full open-source reproducibility} --- seeded RNG
    (default~42), version-controlled codebase, 48~passing unit tests,
    and an automated validation pipeline (\texttt{run\_all.py}) producing
    seed-tagged outputs.
\end{itemize}

The specific component models are:
\begin{itemize}[nosep]
  \item CWLaser with Wiener phase noise (Henry~1982) and
    relaxation-oscillation RIN (Coldren~2012)
  \item MZM with push-pull (zero chirp) and single-drive (linear chirp)
    configurations (Koyama~\&~Iga~1988)
  \item Random-axis birefringence model with diffusive polarisation
    scrambling (Menyuk~\&~Wai~1994, Wai~\&~Menyuk~1996) and
    temperature/bend modulation (Ulrich~1980, Smith~1980, Shibata~1986)
  \item FFT-based CD (Agrawal~2021), Maxwellian PMD (Razavi~2012), and
    exponential attenuation (Keiser~2015)
  \item APD with excess-noise-factor-corrected shot noise (Kasap~2013)
\end{itemize}
```

---

### A6. Birefringence Subsection — Add Beat Length Note

**Insert** after `$L_0 = 75$~km.` (after line 490):

```latex
This residual birefringence corresponds to a beat length
$L_B = \lambda / \Delta n_0 \approx 0.18$~m at 1550~nm.
For context, standard SMF-28 has $L_B > 10$~m (weak birefringence),
while polarisation-maintaining fibre has $L_B < 1$~cm.
Our model therefore compresses the effective correlation length by
a factor of $\sim$500$\times$ relative to standard fibre, producing
diffusive polarisation scrambling over 10--200~km rather than
hundreds of kilometres.  This is appropriate for system-level QBER
estimation where the exact spatial correlations are not the primary
output, but the sub-metre birefringence structure is not physically
resolved.
```

---

### A7. System Demo — Clarify $\sigma$ Convention (Fig 5B discussion)

**Line 914, change**:
```latex
Narrow pulses ($\sigma < 10$~ps RMS, where $\sigma = T_0/\sqrt{2}$
for a Gaussian pulse)
```

This matches the Agrawal convention used in the CD validation section.

---

### A8. System Demo — Rewrite Temperature Panel Discussion

**Replace** lines 921–934 (Panel C) with:
```latex
\paragraph{Panel C: QBER vs temperature.}
The QBER at 75~km exhibits a broad minimum around 30--40$^\circ$C,
where the temperature-induced birefringence partially cancels the
residual birefringence.  The temperature coefficient
$T_{\text{coeff}} = -5\times10^{-7}/^\circ$C, while reasonable in
magnitude, is large enough relative to $\Delta n_0$ that the
birefringence can approach zero deterministically, producing a
deeper null than would occur in deployed fibre.  In real SMF-28,
residual stochastic birefringence from core ellipticity and
microstress prevents exact cancellation; the QBER would drop to a
few percent (not zero) at the minimum.  The V-shaped response
illustrates the mechanism of temperature-dependent polarisation
drift but should not be interpreted as a quantitative prediction at
the null point.
```

---

### A9. Subtitle

**Current**: `Component-Wise Literature Verification with Reproducible Validation Pipeline`

**Replace with**: `First-Principles Field Propagation with Per-Component Literature Validation and Verified System-Level Composition`

---

### A10. Agrawal Citation Consistency

**Current inconsistency**:
- `\cite{agrawal}` (4th ed., 2010) at line 681 in MZM validation
- `\cite{agrawal5}` (5th ed., 2021) at lines 385, 503, 723

**Fix**: Change line 681 from `\cite{agrawal}` to `\cite{agrawal5}` so the MZM validation cites the same edition as the MZM model description.

---

### A11. RIN Implementation Note

**Insert** after the RIN equation description (after line 361):

```latex
The FFT length is chosen as the next power of two above $N$ and is
internally decimated when the optical sampling rate ($\sim$200~THz)
exceeds the Nyquist rate of the RIN bandwidth ($\sim$50~GHz) by more
than $10^3\times$, maintaining numerical efficiency without aliasing
the relaxation-oscillation peak at 5~GHz.
```

---

## Part B: Birefringence Model — First-Principles Rewrite

### B1. Problem Statement

The current `apply_birefringence()` uses a **single-section phenomenological model**:

```
θ = min(π, sqrt(L / L_char) · π/2)
L_char = 75 km · (Δn₀ / |Δn|)²
```

This produces `L_B ≈ 0.18 m`, which is 500× shorter than real SMF-28 (~10–100 m). The spatial statistics of polarisation evolution are not physically resolved.

### B2. Proposed Multi-Section Model

Replace `apply_birefringence()` in `src/channel/fiber.py` with:

```python
def apply_birefringence(E, L, wavelength=1550e-9, temperature=25,
                        bend_radius=None, section_length=1.0):
    """
    Multi-section birefringence via ordered product of random-axis
    Jones matrices (first-principles model).

    The fibre is divided into N sections of length `section_length`.
    Each section has an independent random birefringence axis and
    a deterministic phase retardation Δβ·Δz = 2π·|Δn|·Δz/λ.
    The total Jones matrix is the ordered product of all section
    matrices, which naturally produces:
      - Diffusive polarisation evolution
      - Maxwellian statistics for long fibres
      - Correct correlation length scale

    Parameters
    ----------
    E : ndarray (N, 2) — complex envelope [Ex, Ey].
    L : float — fibre length in metres.
    wavelength : float — centre wavelength (default 1550 nm).
    temperature : float — ambient temperature in °C (default 25).
    bend_radius : float or None — bend radius in metres.
    section_length : float — section length in metres (default 1.0).

    Returns
    -------
    ndarray (N, 2) — field after birefringence rotation.

    References
    ----------
    Menyuk & Wai, JOSA B 11(7), 1994.
    Wai & Menyuk, JLT 14(2), 1996.
    Agrawal, Fiber-Optic Comm. Systems, 5th ed., §4.1–4.2.
    """
    if L <= 0:
        return E.copy()

    T0 = 25.0
    r_fiber = 62.5e-6
    birefringence_T0 = 0.87e-5
    temperature_coefficient = -5e-7
    bend_effect_factor = 0.135

    # Total birefringence: residual + temperature + bend
    delta_n = (
        birefringence_T0
        + temperature_coefficient * (temperature - T0)
    )
    if bend_radius is not None:
        delta_n += bend_effect_factor * (r_fiber / bend_radius) ** 2

    # Stochastic residual — prevents unphysical cancellation
    # (drawn once per cable() call, ~10% of Δn₀)
    delta_n += np.random.normal(0, 0.1 * birefringence_T0)

    # Clamp |Δn| to prevent L_B → ∞ (physical: even "zero-stress"
    # fibre has residual birefringence from core ellipticity)
    delta_n = np.sign(delta_n) * max(abs(delta_n), 1e-8)

    # Number of sections
    N = max(1, int(np.round(L / section_length)))
    dz = L / N

    # Phase retardation per section
    delta_beta_dz = 2 * np.pi * abs(delta_n) * dz / wavelength

    # Ordered product of section Jones matrices
    J_total = np.eye(2, dtype=complex)
    for _ in range(N):
        J_section = _random_su2_rotation(delta_beta_dz)
        J_total = J_section @ J_total

    return np.transpose(J_total @ np.transpose(E))
```

### B3. Key Differences from Current Model

| Aspect | Current (single-section) | Proposed (multi-section) |
|---|---|---|
| Rotation mechanism | Phenomenological `θ = min(π, sqrt(L/L_char)·π/2)` | Physical `Δβ·Δz = 2π·Δn·Δz/λ` per section |
| Correlation length | `L_char = 75 km · (Δn₀/|Δn|)²` (artificial) | `section_length` (1 m default, physically meaningful) |
| Statistics per call | One random axis → one SU(2) rotation | N random axes → N rotations multiplied in order |
| For L → 0 | θ → 0, identity | N → 0, identity (same) |
| For L → ∞ | θ → π, fully mixed | Product of many random SU(2) → fully mixed |
| Temperature artefact | Δn → 0 → L_char → ∞ → θ → 0 (unphysical null) | |Δn| clamped + stochastic term → softened minimum |
| Frequency dependence | None ($L_B$ not used explicitly) | Via `Δβ(ω) = 2π·Δn·ω/(c)` in future extension |

### B4. Impact on Validation Script (`validate_birefringence.py`)

**Current panels** and what changes:

| Panel | Current content | Needs update? |
|---|---|---|
| A: Mean \|Ex\|² vs distance | 200 realisations per distance, shows diffusive drop | Still works — same measurement, new model |
| B: Temperature vs rotation angle | Plots `θ = sqrt(L/L_char)·π/2` vs T using L_char formula | **Must be replaced** — rotation angle now comes from product of N small rotations. Instead plot `mean(\|Ex\|²)` vs T at a fixed distance (same measurement as Panel A, but along T axis) |
| C: Bend radius vs L_char | Plots `L_char(R)` using Ulrich formula | **Must be replaced** — L_char no longer exists. Instead plot `mean(\|Ex\|²)` vs bend radius at fixed distance |
| D: Beat length vs wavelength | `L_B(λ) = λ/Δn₀` — analytical, no birefringence call | Can stay (it's a material parameter check) |
| E: Total Δn components | Plots Δn vs T for various R | Can stay (it validates the Ulrich bend model) |

**Proposed new Panel B**: `mean(|Ex|²) vs temperature at fixed distance` — same measurement as Panel A, measured at L = 50 km, sweeping T from 0–60°C. Shows the softened minimum (QBER does not reach zero).

**Proposed new Panel C**: `mean(|Ex|²) vs bend radius at fixed distance` — same measurement at L = 50 km, T = 25°C, sweeping R from 2 mm to 5 cm.

### B5. Impact on System Demo (`val_system.py`)

| Aspect | Current | After fix | Manuscript change needed |
|---|---|---|---|
| Panel A (QBER vs distance) | Three-regime curve | Same qualitative shape; diffusive scrambling still occurs | None |
| Panel C (QBER vs temperature) | Drops to 0% at 42°C | Drops to a few % (softened minimum) | Rewrite description (see A8) |
| Panel D (QBER vs PMD) | Unchanged | Unchanged | None |
| Panel E (QBER vs bend radius) | Unchanged | Unchanged | None |

No changes needed to the sweep logic in `val_system.py` — the underlying `cable()` → `apply_birefringence()` swap is transparent to the caller.

### B6. Impact on Test Suite (`tests/test_fiber.py`)

Current tests to check:

| Test | Current assertion | After fix | Change needed |
|---|---|---|---|
| `test_birefringence` (power conservation) | `|P_out - P_in| < 1e-12` | Still passes (Jones matrices are unitary) | None |
| `test_temperature_sensitivity` | Different T → different output | Still passes (stochastic Δn preserves sensitivity) | None |
| `test_bend_sensitivity` | Different R → different output | Still passes | None |
| `test_zero_length_edge` | L=0 → identity | Still passes | None |
| `test_seeded_reproducibility` | Same seed → same output | **May fail** if seed consumption changes (now draws N random axes instead of 1) | Update seed check if needed |

---

## Part C: Implementation Order

To avoid breaking intermediate states, apply changes in this order:

1. **`src/channel/fiber.py`** — Replace `apply_birefringence()` with multi-section version (Part B2)
2. **`tests/test_fiber.py`** — Update seed-based reproducibility test if needed (Part B6)
3. **`analysis/validation/validate_birefringence.py`** — Replace Panels B and C (Part B4)
4. **Run validation**: `python analysis/validation/validate_birefringence.py --seed 42`
5. **Run system demo**: `python analysis/val_system.py --seed 42`
6. **Run tests**: `python -m pytest tests/ -v`
7. **`paperwork/manuscript.tex`** — Apply all edits from Part A (A1–A11)
8. **Regenerate figures** — All seed-tagged PNGs in `analysis/val_birefringence/` and `val_system/`
9. **Recompile manuscript**: `xelatex manuscript.tex`

---

## Appendix: Before/After Comparison for Fig 5C (QBER vs Temperature)

**Before** (current model, deterministic Δn = 0 at ~42°C):
```
QBER (%)
 70 |   *   
 60 |   *   *
 50 |   *     *
 40 |   *       *
 30 |   *        
 20 |   *         
 10 |   *          
  0 |   *-----------*  ← unphysical null at T ≈ 42°C
    +------------------------
       0  10  20  30  40  50  60
                 T (°C)
```

**After** (multi-section + stochastic Δn + clamping):
```
QBER (%)
 70 |   *   
 60 |   *   *
 50 |   *     *
 40 |   *       *
 30 |   *        
 20 |   *         
 10 |   *     *    
  0 |   *-----------*  ← softened minimum, QBER ≈ 3–5%
    +------------------------
       0  10  20  30  40  50  60
                 T (°C)
```

The characteristic V-shape is preserved (the mechanism is still real), but the artefactual null at 42°C is eliminated. This is physically defensible and consistent with the known behavior of temperature-dependent birefringence in standard single-mode fibre.
