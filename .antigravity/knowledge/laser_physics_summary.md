# Knowledge Item: Laser Simulation Physics & Rate Equations
## Repository: opto-sim
## Date: 2026-05-14

### 1. Overview
This KI documents the fundamental physics corrections and modeling decisions made during the development of the `opto-sim` framework. It serves as a primary context source for the transition from Er:Yb models to stable 4-level Nd:YAG models.

### 2. Key Physics Decisions & Corrections

#### A. Er:Yb Quasi-3-Level Model (`sslaser.py`)
- **Status**: Research-grade but numerically "stiff".
- **Corrections**:
    - **Unit System**: $N_0$ (ion density) in literature (Rao et al.) is often cited in $cm^{-3}$ ($8.59 \times 10^{20}$). For SI consistency in the simulation, this must be used as $8.59 \times 10^{26} m^{-3}$.
    - **Threshold Logic**: Population inversion ($N_2 > N_1$) in this system requires a pump rate $R_p > 1/\tau_2$. Pumping at exactly $1/\tau_2$ only reaches the transparency limit, not lasing.
    - **Atom Conservation**: Enforced $N_1 = N_0 - N_2$ to ensure physical consistency.

#### B. Nd:YAG 4-Level Model (`ndyag.py`)
- **Status**: Recommended for stable algorithm verification (Polarization/Detection).
- **Physics**: 
    - Implemented a canonical 4-level system based on **Saleh & Teich** and **Koechner**.
    - **Inversion Advantage**: Since the lower laser level decays in nanoseconds, $N_{lower} \approx 0$, ensuring $N_{upper} - N_{lower}$ is always positive during pumping.
    - **Parameters**: $\lambda = 1064nm$, $\tau_2 = 230\mu s$, $\sigma = 2.8 \times 10^{-23} m^2$.

#### C. Stochastic Langevin Noise Model
- **Implementation**: Integrated directly into the $dI/dt$ derivative in the `rate()` function.
- **Stability Fix**: To maintain compatibility with adaptive-step solvers (BDF/RK45), noise is generated using a **deterministic time-seed** (`np.random.seed(int(t * C))`). This ensures the solver sees a consistent derivative at every probe point, preventing infinite step-size reduction.

### 3. File Map for Context
- `src/lasers/ndyag.py`: Most stable reference model.
- `ndyag_characterization.py`: Verification script showing relaxation oscillations.
- `analysis/ndyag_dynamics.png`: Visual proof of stable inversion.
- `repository/`: Contains the 22-category research library (Papers/Datasheets).

### 4. Next Steps & Open Threads
- **Polarization**: Integrate the Stokes/Jones calculus with the high-gain Nd:YAG model.
- **Modulation**: Implement phase modulation dynamics using the datasheets in `repository/Phase Modulator`.
- **Detection**: Use `repository/Photodiodes` to refine the APD/SPAD jitter models.
