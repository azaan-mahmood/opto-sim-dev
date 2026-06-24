# Research Roadmap: opto-sim Physical-Layer QKD Simulator

A tiered plan with effort estimates and diminishing-returns checks.

---

## What makes this simulator different

Most QKD simulators commit one of two policies.

**Sin 1 — Abstract qubit layer**: NetSquid, SimulaQron, QuTIP-based simulators represent quantum states as vectors
|\psi⟩ and channels as depolarising maps. A "detector" clicks with some probability. There is no field, no wavelength,
no phase noise, no birefringence, no *physics* between Alice's laser and Bob's APD.

**Sin 2 — Numeric microwave**: Lumerical, RSoft, OptiSystem are exquisite for photonic circuit design but are
closed-source, priced per seat, and designed for sinusoidal carriers — not for QKD, where the complex envelope carries
the quantum information.

**opto-sim is neither**. It propagates the complex-envelope electric field (power-bearing, calibrated once at the laser)
through physically parametrised components: a CWLaser with Wiener phase noise and relaxation-oscillation RIN,
a PhaseModulator whose V_π comes from LiNbO₃ crystal geometry, a MZM built from two modulators in an interferometer,
a fibre whose birefringence *depends on temperature and bend radius* via analytic formulae, and an APD whose shot noise
carries the excess-noise factor from the gain process.

This permits research questions that the other simulators simply cannot address:

- How does a 10 °C temperature swing on an aerial fibre change the sifted key rate?
- Can the relaxation-oscillation peak of the RIN spectrum leak information to an eavesdropper?
- What is the joint effect of CD + PMD + polarisation drift on BB84 phase encoding at 10 Gbaud?

---

## Tier 0 — Reproducibility & testing infrastructure

| Effort | Deliverable | Paper |
|---|---|---|
| 1–2 weeks | pytest suite + seeded-reproducibility | None (scaffolding) |

**What to do:**
1. Pin a fixed RNG seed per test so every "random" fibre realisation is repeatable.
2. Write unit tests for each component (CWLaser power & linewidth, MZM transfer curve, CD pulse broadening, APD noise variance).
3. Add a `tests/` directory with `conftest.py` and a `Makefile` / `tox.ini`.
4. Add a `--seed` CLI argument to BB84 scripts.

**Why first:** Without seeded randomness, every figure in every paper will be unreproducible.
Every hour spent here saves days of debugging false QBER changes later.

**Stop when:**
- Results are reproducible (same seed → same figures).
- Unit tests catch obvious regressions.

**Not:**
- 95 % coverage.
- Enterprise CI/CD pipeline.
- Perfect software engineering.

---

## Tier 1 — Channel validation and characterisation

| Effort | Deliverable | Paper |
|---|---|---|
| 2–4 weeks | Validated fibre model against published data | Conference (OFC, CLEO, ECOC) |

**What makes this novel:** No QKD simulator validates chromatic dispersion, PMD, and
birefringence against separate experimental datasets. Doing so makes the simulator defensible
in peer review.

**Tasks:**
1. **CD validation** (3 d). Reproduce Agrawal Fig 2.6 (Gaussian pulse broadening at z/L_D = 0.5, 1.0, 2.0).
   Fit β_2 from SMF-28 datasheet (D = 17 ps/(nm·km)) and confirm width ratio within 0.1 %.
2. **PMD validation** (3 d). Generate 10^4 fibre realisations, plot DGD histogram against
   Maxwellian PDF (Razavi [5] Fig 2.11). Compute ⟨Δτ⟩ and confirm it matches
   the PMD coefficient × √L.
3. **Attenuation validation** (1 d). Plot received power vs distance for an SMF-28 link
   at 1550 nm (0.182 dB/km). Overlay OTDR data from Keiser [1] or a Thorlabs application note.
4. **Birefringence validation** (2 d). Reproduce Yuan [4] Fig 1: beat length L_B vs fibre
   bend radius for various Δn_0. Confirm L_B ∝ R / Δn_0 scaling.
5. **Reproducibility pass** (2 d). All scripts save PNG & CSV with (seed, params) in the
   filename. A reader can re-run with `python run_all.py --seed 42`.

**Effort breakdown:** 11 d coding, 3 d writing, 1 d figures = **3 weeks**.

**Stop when:** Each model reproduces its reference figure within reasonable tolerance.
Spending weeks to shave the last 0.01 % off the CD fit yields zero research value —
the question is whether the model is *physically correct*, not whether it is perfect.

---

## Tier 2 — Environmental sensitivity of QKD

| Effort | Deliverable | Paper |
|---|---|---|
| 4–8 weeks | Temperature/bend → QBER curves, polarisation drift analysis | Journal (Optics Express, IEEE JLT, EPJ QT) |

**What makes this novel:** The birefringence formula in `fiber.py` (Yuan [4]) makes temperature
and bend radius *explicit* parameters. Every other QKD simulator treats polarisation as a
random unitary. This model can answer: *would burying the fibre vs aerial deployment change
the QBER statistics measurably?*

**Tasks:**
1. **Temperature sweep** (1 w). 0–50 °C, 25 km spans. Plot Δn(T), beat length, and resulting
   QBER for BB84 at 1 Gbaud. Does the phase shift θ(T) cross a zero (all bits flip)? How fast?
2. **Bend sensitivity** (1 w). R = 10–100 mm, fixed T. Plot Δn(R) via Yuan Eq 8–9.
   Coupled with temperature: at what (T, R) combination does the fibre become visibly
   polarisation-maintaining (L_B < 1 m)?
3. **Time-varying channel** (2 w). Inject a slow temperature ramp (1 °C/min) and run repeated
   BB84 exchanges. Plot QBER as a time series. Does the reconciliation protocol recover?
   How large a temperature swing is tolerable before the key rate drops below zero?
4. **Aerial vs buried comparison** (1 w). Use published weather data (ASHRAE handbook) to bound
   diurnal temperature variation. Compare QBER variance for aerial (fast, large swings) vs
   buried (slow, damped) scenarios.
5. **Write-up** (1 w). Figures, discussion, comparison to literature.

**Effort breakdown:** 5 w coding + simulation, 1 w writing, 1 w review = **7 weeks**.

**Stop when:** The effect size is small. If aerial fibre causes only 0.2 % additional QBER,
that may be physically interesting but not publication-worthy. Abandon Tier 2 and move to
Tier 3.

The environmental story is strongest when it makes a concrete, falsifiable prediction
("QBER variance in aerial fibre is 3× larger than buried"). If the signal is weak, do not
force a paper.

---

## Tier 3 — Physical-layer security analysis

| Effort | Deliverable | Paper |
|---|---|---|
| 3–5 months | Side-channel quantification, finite-key analysis | High-impact journal (IEEE T-IFS, NPJ Quantum Inf, PRX Quantum) |

**Why it matters:** Publication value ∝ **Novelty × Credibility × Importance**.

- Validation papers (Tier 1–2): high credibility, lower novelty.
- Security side-channel papers: high novelty, potentially high importance.

The physical-layer approach lets you simulate attacks that abstract simulators cannot even
represent. **This is where the simulator earns its keep.** If opto-sim can reveal a measurable
leakage mechanism that existing abstractions miss, that is a much stronger story than yet
another QBER-versus-distance graph.

**Tasks:**
1. **Laser noise side-channel** (3 w). CWLaser's RIN resonance (Coldren Eq 5.3.38) produces
   correlated amplitude fluctuations. Can Eve synchronise to the relaxation peak and extract
   bit information from the amplitude envelope? Quantify leaked information in bits/symbol.
2. **Finite-key effects with realistic noise** (4 w). Derive SKR from simulated QBER using
   the canonical GLLP + decoy-state formulas. Show how finite-key bounds change when the
   underlying noise is correlated (RIN, 1/f phase noise) rather than i.i.d. This is a
   genuinely open question in QKD security proofs.
3. **Detector side-channel** (2 w). APD dead time and afterpulsing probabilities create
   detector-efficiency mismatches. Simulate time-shift attacks (Makarov 2009) in the
   physical-layer framework. Confirm the attack success rate vs QBER penalty.
4. **Phase modulator flaw analysis** (2 w). MZM drive voltage ripple
   (V = V_π + δ·sin(ωt)) creates a deterministic phase error pattern. Can Eve exploit this
   pattern to learn Alice's basis choice?
5. **Write-up** (3 w). A full journal paper.

**Effort breakdown:** 3 m coding, 1 m writing, 1 m review = **5 months**.

**Stop when:** You have one strong, publishable result where the physical-layer approach
provides a clear advantage over abstract models. Do not implement every attack in the
literature. Pick the one where opto-sim uniquely shines (likely the laser-noise side-channel
or the correlated-noise finite-key analysis). One strong result beats a catalogue of re-runs.

---

## Tier 4 — Full QKD system optimisation

| Effort | Deliverable | Paper |
|---|---|---|
| 6–12 months | Decoy-state protocol, SKR optimisation, open-source release | Nature Communications, PRX Quantum + code release |

**Why it matters:** This closes the loop: physical model → protocol → secure key rate →
optimisation. A well-documented open-source framework with a reproducible paper attached
would be unique.

**Tasks:**
1. **Decoy-state protocol** (2 m). Implement 3-intensity decoy-state BB84
   (Lo, Ma & Chen 2005). Optimise μ, ν, ω intensities for the physical laser model. Show
   how finite-intensity statistical fluctuations interact with actual RIN to change the
   secure bound.
2. **SKR optimisation** (2 m). Multi-parameter sweep: laser power, modulation depth, detector
   bandwidth, fibre length, temperature. Produce contour plots of SKR in (distance, power)
   space. Include detector dead-time limitation at high rates.
3. **Framework clean-up and documentation** (1 m). API docs, tutorial notebooks, a one-command
   demo. Publish to PyPI or as a GitHub release with a DOI.
4. **Write-up** (2 m). The definitive paper showing what a physical-layer QKD simulator can do.

**Effort breakdown:** 6 m coding, 2 m writing, 2 m review = **10 months**.

**Diminishing returns:** This tier is where feature creep is most dangerous. Resist the urge
to add every protocol (E91, MDI-QKD, TF-QKD), every detector (SNSPD, TES), and every
modulation format (PPM, DPSK, CSK). The goal is *one* complete system story, not a tool
chest. Adding more protocols past the first is a **strongly diminishing** activity.

---

## Summary

| Tier | Time | Impact | Risk |
|---|---|---|---|
| 0: Reproducibility | 1–2 w | Scaffolding | None |
| 1: Channel validation | 3 w | Conference | Low |
| 2: Environmental | 7 w | Journal | Medium |
| 3: Security | 5 m | High-impact J. | Medium–High |
| 4: Full system | 10 m | Flagship | High |

### Recommended route

1. **Tier 0** (2 w). Every tier needs it.
2. **Tier 1** (3 w). Low-risk calibration paper at OFC/CLEO.
3. **Tier 2** (7 w). If the environmental signal is strong, publish at Optics Express.
   If weak (≤0.2 % QBER difference), **abandon** and skip directly to Tier 3.
4. **Tier 3** (5 m). This is the highest impact per hour spent for this codebase —
   *no other simulator can do the laser-noise side-channel analysis*. Invest here.
5. **Tier 4** only if the project needs a flagship capstone. The marginal impact of the
   tenth contour plot is much smaller than the first novel attack.

### One-page version

> **Do this first:** Tier 0 + Tier 1 (5 weeks total).  
> **Then publish** a channel-validation conference paper.  
> **If Tier 2's signal is strong** (≥1 % QBER delta), publish an environmental journal paper.  
> **If not**, skip to Tier 3.  
> **Then do** the laser-noise side-channel (3 months).  
> **Stop there.** You will have two or three papers and a defensible, validated simulator —
> the diminishing returns past this point are steep.
