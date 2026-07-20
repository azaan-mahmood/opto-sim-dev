# Tier 1 — Channel Validation: Plan of Attack

## Goal
Validate every physical model in the fibre channel against published experimental/literature data, making the simulator defensible in peer review.

## References
- Agrawal [6] — *Fiber-Optic Communication Systems*, Wiley
- Razavi [5] — *An Introduction to Quantum Communications*, Cambridge
- Keiser [1] — *Optical Fiber Communications*, McGraw-Hill
- Yuan [4] — "Birefringence in bent optical fiber", JLT 2020

---

## Task 1 — CD Validation (`analysis/validate_cd.py`)

**Reference:** Agrawal [6] Fig 2.6, Eq 3.2.6

**Method:**
- Launch Gaussian pulse through MZM (σ = 5–30 ps)
- Pass through `propagate()` at z/L_D = 0.5, 1.0, 2.0
- Fit output |E|² to Gaussian, measure σ(z)
- Compare to analytic: σ(z) = σ₀·√(1 + (z/L_D)²)

**Target:** Width ratio error < 0.1 %

**Output:** `analysis/val_cd_{seed}.png`, `val_cd_{seed}.csv`

---

## Task 2 — PMD Validation (`analysis/validate_pmd.py`)

**Reference:** Razavi [5] Fig 2.11

**Current state:** `fiber.py` L145 uses `np.random.rayleigh` (2D magnitude)

**Method:**
1. Generate 10⁴ fibre realisations with PMD ON (CD/birefringence OFF)
2. Extract DGD per run (from applied Jones matrix)
3. Histogram vs **both** Rayleigh and Maxwell PDFs (let data decide)
4. Compute ⟨Δτ⟩, confirm ⟨Δτ⟩ = PMD_coeff · √L
5. If Rayleigh ≠ literature expectation → fix to Maxwell

**Output:** `analysis/val_pmd_dgd_{seed}.png`, `val_pmd_{seed}.csv`

---

## Task 3 — Attenuation Validation (`analysis/validate_attenuation.py`)

**Reference:** Keiser [1] Eq 3.6, SMF-28 datasheet (0.182 dB/km at 1550 nm)

**Method:**
- Sweep L = 0–200 km, measure P_out / P_in
- Overlay: ideal 10^{-αL/10}

**Output:** `analysis/val_att_{seed}.png`, `val_att_{seed}.csv`

---

## Task 4 — Birefringence Validation (`analysis/validate_birefringence.py`)

**Reference:** Yuan [4] Fig 1

**Method:**
- Sweep bend radius R = 5–100 mm at fixed T = 25°C
- Compute Δn(R, T) via Yuan Eq 8–9
- Plot beat length L_B = λ/Δn vs R for Δn₀ = {0.5, 1.0, 2.0}×10^{-5}
- Confirm L_B ∝ R / Δn₀ scaling

**Output:** `analysis/val_biref_{seed}.png`, `val_biref_{seed}.csv`

---

## Task 5 — Reproducibility Pass

**Affects:** All validation scripts + `run_all.py`

**Changes:**
- Every script saves PNG + CSV with `{name}--seed{seed}--{params}.{ext}`
- New `run_all.py` at project root executes all 4 validations with `--seed 42`

---

## Estimated Effort

| Task | Days | Deliverable |
|------|------|-------------|
| CD validation | 1 | `validate_cd.py` |
| PMD validation (+ fix) | 1.5 | `validate_pmd.py` + possible `fiber.py` edit |
| Attenuation validation | 0.5 | `validate_attenuation.py` |
| Birefringence validation | 1 | `validate_birefringence.py` |
| Reproducibility pass | 1 | `run_all.py`, filename conventions |
| **Total** | **5 d** | 4 validated models + run_all.py |

## Key Decisions
- PMD: Validate current Rayleigh sampling first, compare to Maxwellian literature — fix only if needed
- No `propagate()` API changes expected (unless PMD fix changes DGD sampling)
- All scripts accept `--seed` for reproducibility
