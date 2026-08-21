# Changelog

All timestamps are local time (UTC+5).

---

## 2026-08-20 — six validators that could not fail, and a theory curve that was wrong

### Session: an audit, and the first finding from it

An audit of the repository turned up five findings. This lands the first
two, which turned out to be one thing.

| Change | Files | Rationale |
|---|---|---|
| Real pass/fail criteria | `validate_{apd,attenuation,cd,cwlaser,mzm,pmd}.py` | Six roster validators had no `[PASS]`, no `[FAIL]`, no non-zero exit. They computed, wrote a table, and exited 0 unconditionally. |
| Hardcoded verdicts replaced by measurements | the five committed component tables | Seven claims across four artifacts asserted agreement that nothing had checked. |
| APD theory curves corrected | `validate_apd.py` | The curves the model was being compared against were themselves wrong, in two ways. |
| Harness warning restored to one line | `run_all.py` | It was firing six times a run, which is how it went unnoticed. |

**Six validators reported PASS for running.** `validate_apd`,
`validate_attenuation`, `validate_cd`, `validate_cwlaser`,
`validate_mzm` and `validate_pmd` contained no verdict of any kind. The
harness reported PASS because the script had not crashed and its output
file existed. Every one of them already computed a comparison against an
analytic prediction and simply never asserted it, so none of the criteria
below is invented — each is a comparison the script was already drawing.

They now emit 2, 2, 2, 3, 3 and 4 `[PASS]` lines (CD, Attenuation, APD,
PMD, CWLaser, MZM) and exit non-zero on failure. Every one of the six had
its failure path exercised: four fired on their own during development,
and the remaining two were broken deliberately — alpha moved 0.182 to 0.2
under a theory curve left at 0.182 (56 % departure, caught), and the PMD
coefficient moved 0.10 to 0.12 under an expectation left at 0.10 (48
sigma on the mean, 50 on the RMS, and the Maxwell KS test collapsing to
p = 5e-260).

**The warning that should have caught this was itself the casualty.**
`run_all.py` exempted PMD by name from its "no `[PASS]` markers" warning.
By the time anyone looked, six validators were tripping it every run, so
the one signal that could have surfaced the problem had become noise. The
exemption is gone; the warning now fires exactly once, for
`validate_gobby`, which is genuinely the last roster entry with no
verdict — it can fail on the statistical-power guard but not on physics.
That is left firing on purpose rather than silenced.

**The APD validator was comparing the model against the wrong curves.**
Tightening its gates surfaced two defects, and in both the model was
right:

* Its "theory" responsivity used the exact speed of light while the
  detector uses `self.c = 3e8`. That 0.0692 % difference was *exactly*
  the disagreement observed.
* Its noise theory omitted the dark-current shot term that
  `calculate_noise` includes — a 4.8e-3 % gap, small enough to look like
  rounding on a log plot and to survive being called "verified".

Formed from the detector's own constants and with the dark term restored,
all five comparisons agree to **1.75e-14 %**. A comparison curve built
with different constants tests the constants, silently, and reports the
answer as if it had tested the formula.

Whether `c = 3e8` is good enough is a separate question and is now
*reported* rather than gated: the validator prints c +0.0692 %,
kB +0.0470 %, e +0.0110 % and h +0.0011 % against CODATA 2018. Changing
them moves every APD number in the project, so it is a decision, not a
check.

**Two more claims that were never true as stated.** `validate_mzm`
reported `P_out = 0.5` at quadrature and `0` at the null; `V_scan` is 200
points across `[0, 2*V_pi]` and contains neither voltage — the nearest
samples sit at 0.5025 and 0.9950 of `V_pi`, where cos^2 gives 0.49605 and
6.2e-5. The three points are now evaluated at the voltages being claimed,
and come back 1.00000000, 0.50000000 and 2.6e-32. `validate_pmd` wrote
`linregress`'s `r` under a heading reading `R2`; both sit near 1, so it
never showed.

**What the criteria are.** Where the comparison is closed form against
closed form the gates sit far above the numerical residual and catch
gross failures — a transfer that has stopped being cos^2, an insertion
loss applied in amplitude rather than power, a crystal cut modulating the
wrong component. Where it is a draw, the tolerance comes from the
sampling distribution rather than from a chosen number: PMD's mean and
RMS DGD are gated at five standard errors derived from the Maxwell
moments, so they tighten with `--n-realizations` instead of being fixed
at whatever one run produced.

Suite unchanged at 378 passed, 1 skipped.

---

## 2026-08-20 — the drift curve, and two amplitudes that are easy to confuse

### Session: drift rate against rate loss, servo on

| Change | Files | Rationale |
|---|---|---|
| Drift-rate sweep, servo on | `analysis/validation/validate_gobby_impairments.py` | The two drift rows show the paper's sentence at two points. This shows it as a curve, which is the one thing a table cannot. |
| `mean_residual_amp()` | same | The rate ratio is a time average, so its prediction has to be one too. |
| `residual()` / `_aligned()` | same | One definition of Bob's alignment convention, for the operator table and the sweep both. |
| Third figure panel | same | Rate falling while QBER does not, in one axes. |

This is the carve-out kept when the sensitivity study was retired: drift
rate against **rate loss** rather than against QBER, which nothing
measured. It lands where the servo and both drift constants already live,
as `DRIFT_SWEEP_C_S`, rather than as a new validator.

**The operator chose the range, and it is checked rather than trusted.**
`|R00|^2` falls with drift rate only while the residual rotation stays
inside about one turn. At 10 km, 3e-3 C/s gives 0.095 and 1e-2 gives
0.963, with the accumulated phase going +4.41 rad to -5.92 — the rotation
wraps and the amplitude comes back. A sweep through that would draw a
curve that rises again, and the rise would report where the rotation
happened to land, not what drift costs. That is the trap
`validate_duplinskiy_drift.py` already records for bend radius, where 1 m
is worse than 0.1 m. So the sweep stops at 3e-3, and monotonicity over
exactly the swept range is now a check, at all three operator distances.

**The endpoint amplitude is not the prediction, and using it failed
loudly.** The first draft predicted the rate ratio from `|R00|^2` at
t = 120 s and missed at 3e-3 by a factor of six — 0.5644 measured against
0.0953. Bob aligns at t=0, so the residual *starts* at the identity and
walks away from it while pulses are collected the whole time; what a rate
ratio can equal is the time average, not the final instant.
`mean_residual_amp()` now samples at the model's own 100 block midpoints
rather than as a continuous integral, so prediction and model agree about
what "during the run" means instead of nearly agreeing. It returns 0.572
at 10 km and 3e-3 C/s, against 0.603 +/- 0.027 measured — 1.1 sigma.

That number also reproduces a comment already sitting in the file, which
said the mean rate factor at 3e-3 is 0.57 at 10 km. Agreement with a
number written before the function existed is worth more than agreement
with one written after.

Both amplitudes are printed side by side for exactly that reason: they
differ sixfold and only one of them is a prediction.

**What is asserted.** Per point, the rate ratio against the mean
`|R00|^2` — a servo turns a phase, so the amplitude a unitary took is
precisely what the rate should keep. Per point, that the QBER does *not*
follow. Across the range, that the rate falls end to end. Adjacent steps
are reported with the measurement, not asserted: the first few predicted
separations are ~1e-4 in ratio against counting noise of order 1e-2.

`DRIFT_BLOCKS` is now a named constant passed to both the simulation and
the prediction, so the two cannot silently diverge on how often the fibre
is re-evaluated.

Full budget 12 min 53 s, all checks passing. Suite unchanged at 378
passed, 1 skipped.

---

## 2026-08-20 — a study closed rather than deferred again, and three small items

### Session: the sensitivity study is retired; grid duplication, partial blocks, FFT churn

| Change | Files | Rationale |
|---|---|---|
| `field_grid()` deduplicated | `src/protocols/bb84_time_bin.py` | Two copies of one derivation, and the shape check compares a caller's field against one of them. |
| Partial-block behaviour documented | `src/protocols/bb84_duplinskiy.py` | A trailing short block is dropped. Deliberate, and nothing said so. |
| FFT nondeterminism documented in the artifact | `analysis/validation/validate_dfb_drive.py`, `val_dfb_drive--cw_rin.csv` | The explanation lived in a local file nobody diffing the CSV would ever see. |
| Sensitivity study retired | — | Decided, not deferred a third time. See below. |

**The sensitivity study is closed.** It was carried across sessions as
pending work, deferred twice. It has no git history, no spec, and — as it
turns out — no entry in the working note either: it existed only as a
name. Rather than defer it again, it is now a decided entry alongside the
phenomenological birefringence model.

Its content is already covered four ways: `validate_gobby.py` for QBER vs
distance, `validate_gobby_impairments.py` for impairment-by-topology with
the servo, `validate_duplinskiy_drift.py` for drift before recalibration,
and `val_duplinskiy_scenarios.py` for the ladder. More decisively, the
physics argues against a broad sweep and that argument was already
written down in the drift validator: past the calibration tolerance the
residual is a fixed unitary rotation, measured at 79.9 % QBER for
R = 1 m against 44.8 % for R = 0.1 m — a tighter bend not monotonically
worse — so "the particular value is an accident of the realisation." A
wide sweep would largely plot realization noise on a smooth-looking axis.

One carve-out is kept and located rather than lost: drift rate against
*rate loss* with the servo on is genuinely unmeasured, and belongs as a
short list of rates in `validate_gobby_impairments.py`, where the servo
and both drift constants already are — not as a new 34-minute validator.

**The grid deduplication had to be bit-identical, and that was tested
rather than argued.** `field_grid()` returns `(dt, n_samples,
pulse_center)`; `simulate_bb84_time_bin` recomputed all three with the
same expressions, and the shape check for a caller-supplied source field
compared against its own copy. The second site now calls the first. The
reasoning that identical expressions in identical order give identical
floats is sound and is still not evidence, so it was checked by running:
five configurations — both topologies, a birefringent fibre at 10 km, a
non-default `pulse_width`/`delay` grid, and the supplied-source shape
path — all identical to the last digit of QBER.

**The FFT churn is now documented where it is met.** Regenerating
`val_dfb_drive--cw_rin.csv` moved 1108 of 8193 rows, max relative change
9.98e-7 — one unit in the last of the seven figures `%.6e` prints — while
the summary line, the frequency axis and `val_dfb_drive.png` came back
byte-identical. Three header lines now say so in the file itself, so the
next person to diff two runs reads the explanation instead of
rediscovering it. The nondeterminism itself is not chased; documenting it
is the whole fix.

Suite unchanged at 378 passed, 1 skipped.

---

## 2026-08-20 — the last unchecked artifact, and a threshold that was a fitted number

### Session: `val_duplinskiy_scenarios` becomes a checked validator

| Change | Files | Rationale |
|---|---|---|
| `DUPL-scenarios` joins the harness | `run_all.py` | It was the only script producing committed artifacts that no regression check covered. Since the four proof-of-concept scripts retired, that was a one-item list rather than a survey. |
| Ten checks added | `analysis/val_duplinskiy_scenarios.py` | A roster line alone would have bought a validator that could only fail by crashing. |
| `_stem` adopted | `analysis/val_duplinskiy_scenarios.py` | A live defect, not a fit issue. See below. |
| Cross-reference read, not copied | `analysis/val_duplinskiy_scenarios.py` | A hand-copied number from another validator can rot quietly while still looking authoritative. |
| Roster count and `seeded` note corrected | `run_all.py` | Nine validators take `--seed` now, not eight. 20 validators, not 19. |

**`--allow-underpowered` was overwriting the quotable table.** The stem
was driven by `--quick` while the power guard was waived by a different
flag, so the one flag that declares a run unquotable was the one the
filename ignored. Running it printed "Do NOT cite it as a result" and then
wrote straight over `val_duplinskiy_scenarios--seed42.csv`. This is
exactly the failure `validate_gobby._stem` was written to prevent, and its
docstring describes this scenario; scenarios had simply never adopted the
helper. It now does, with `reduced = quick or allow_underpowered`.

**A threshold tuned on one seed is a fitted number.** The first draft of
the checks asserted magnitudes — that a 0.05 C calibration mismatch costs
more than 15 % QBER, that the uncompensated control exceeds 25 %. Running
the checks at three seeds broke it: **0.05 C of mismatch costs 35.9 % at
seed 42 and 4.0 % at seed 7**, and neither run is wrong. How far the Jones
matrix turns per degree is a property of the fibre draw, not of the model.

Those four checks now assert error counts — an exact inverse is exactly
inert, so a mismatch must produce a non-zero count, and that holds for
every draw — and the magnitudes are printed under "reported, not
asserted". The tell for the general case is that the threshold could only
be chosen by looking at one run's output.

**What the ten checks are.** Four exact: the compensated rows carry zero
errors; the DFB source swap is inert without dispersion; CD is a
structural null; and the paper's lab row is bit-identical to the
dark-count row it deliberately repeats. That last one is the sharpest —
the lab row spells out `mu=0.1`, `gate_width=20e-9`, `rep_rate=10e6` where
the dark-count row relies on those being `simulate_bb84_duplinskiy`'s
defaults, so the check silently guards the three defaults. Confirmed by
moving `mu` and watching it fail with a readable message.

Four assert a row is not inert. Two are orderings — the impairment ladder
end to end (5.9 sigma) and urban above lab (3.1 sigma) — and both are
skipped on a reduced run, because one sigma at the smoke budget is about
1.2 pp against steps of 0.4 pp. The gate is the budget and never the
observed gap: a ladder that had genuinely stopped climbing would show a
gap near zero, hence no resolution, hence a skip, and the check would
never fire on the one case it exists for.

**The ladder's adjacent rungs are reported and not asserted.** Their steps
are +0.39, +0.56 and +1.36 pp at 1.2, 1.6 and 3.3 sigma, so two of the
three are not resolvable even at full budget. Asserting adjacency would
have been a check that fails at random; the end-to-end climb is asserted
instead, and each step is printed with its significance.

**The urban row's offset is not this script's to resolve.** It reads
3.28 +/- 0.25 % against the paper's 5.5 %. `validate_duplinskiy_urban.py`
measures 3.38 +/- 0.31 % for the same configuration and already carries
the gap; the two agree with each other, which is the cross-script
consistency this table exists to show. So the ordering is asserted and the
offset is printed — and the comparison number is read out of the urban
validator's committed CSV rather than copied into a constant here.

**Costs, measured.** Quick 4.6-6.0 min across three seeds, 26.0M pulses.
Full 41 min 35 s, 308.6M pulses, 124k pulses/s. The working note recorded
this script at "~9 min", wrong by about 4.6x. An extrapolation from the
quick run mid-session gave 77 min, wrong the other way, because the quick
run carries proportionally more pilot and DFB-build overhead — the
arithmetic from the chain's known throughput was right and the
extrapolation from a smaller run was not.

**The regenerated table is byte-identical to the committed one** apart
from the commit hash in its header. Nothing here touched physics, and the
artifact says so. Suite unchanged at 378 passed, 1 skipped;
`run_all.py --list` now shows 20 validators.

---

## 2026-08-19 — the servo Gobby's receiver actually has, and a clear-out

### Session: phase servo and piezo stretcher, four scripts retired, two false statements cut

| Change | Files | Rationale |
|---|---|---|
| Bob's phase servo | `src/protocols/bb84_time_bin.py` | The drift model contradicted the paper it replicates. Gobby report that drift costs bit rate and not QBER; ours did the reverse, because a residual rotation is SU(2) and so costs O(eps^2) in amplitude against O(eps) in phase. Their receiver holds its operating point with a piezo stretcher, and that component was missing. |
| `PiezoFibreStretcher` | `src/channel/piezo_stretcher.py` (NEW) | The servo was a scalar knob with no range, no drive voltage and no device. Every other piece of hardware here lives in `src/channel/` with a cited formula. |
| Eye diagram ported out of the retiring script | `src/visualization/eye.py` (NEW), `validate_cwlaser.py`, `validate_dfb_drive.py` | The one thing worth keeping from `laser_characterization.py`. Porting it exposed three faults it had always had. |
| Four proof-of-concept scripts retired | `analysis/examples/` (NEW), `git mv` of `val_system.py`, `val_system_scenarios.py`, `qber_vs_distance_dispersion.py`, `laser_characterization.py` | Never literature-validated: their numbers could be checked only against themselves. The protocol chains replaced them. |
| Two false statements cut from an artifact header | `analysis/val_gobby/validate_gobby.py` | Not merely verbose — wrong. See below. |
| Figure titles | `validate_gobby.py`, `validate_dfb_drive.py`, `validate_duplinskiy_urban.py` | A figure travels away from the script that made it. |
| Two reproducibility defects fixed | `validate_duplinskiy_calibration.py`, `val_cwlaser--seed42.png` | Found by running the harness, not by reading the code. |
| Full Gobby sweep re-run | `val_gobby--seed42.{csv,png}`, `val_gobby_table.tex` | The regression check for the whole fibre-and-drift workstream. |

**The servo, and why it costs nothing to run.** The closed form is
`P = g0 + 2*Re(S*exp(i*delta))`, so a fibre phase `theta` sends
`S -> S*exp(i*theta)` and shifts `arg(S)` by exactly `theta`. `arg(S)` is
the fringe phase, which is what a real lock-in servo measures, so the
error signal was already being computed once per block. Verified sharply
rather than assumed: with the fibre frozen, the servo is bit-identical to
cancelling `2*arg(R11)` by hand.

It is **phase only**, which is the entire point. Inverting the full Jones
matrix would null the amplitude too and the rate loss would vanish with
it. Measured at 10 km and 3e-3 C/s over the paper's 120 s transfer: the
sifted rate falls to 0.589 of reference with the servo on and 0.589 with
it off, while the QBER goes from 41.6 % to 0.054 % against a 0.095 %
baseline. Sample-and-hold gives the crossover — 0.00 %, 0.39 %, 4.09 % and
14.57 % at re-lock intervals of 1.2 s, 12 s, 60 s and never.

Not offered on `bb84_duplinskiy`: there the residual is a full SU(2)
acting on the encoding itself, so a phase-only correction has nothing to
correct.

**The stretcher is a specific part.** Thorlabs FVP155P, DOC-103641 Rev B,
supplied by the user after asking for a citable design. Every default is
one line of that datasheet — half-wave voltage < 20 V, 7 pi stroke at
150 V, 80 kHz resonance, 0.1 dB insertion loss, 0.15 % residual AM — and
the values appear in that file only. The physics under it is Butter &
Hocker, *Appl. Opt.* 17(18) 2867 (1978).

`voltage_for` wraps into one fringe before converting, which is what makes
a finite stroke sufficient: an interferometer's operating point is
periodic, so a 5 pi correction and a pi one are the same point. The demand
therefore never exceeds 2*v_pi = 40 V against a 150 V limit however far
the fibre has drifted.

Insertion loss is applied **by default** and suppressed explicitly in the
Gobby chain, where `ETA_BOB = 0.045` already folds in "5 dB of loss in
Bob's apparatus". Residual AM is carried and not applied: the datasheet
bounds its size without giving its dependence on drive, and applying it
would mean inventing that shape.

**Two statements in a committed header were false, not just long.** It
told the reader afterpulsing "is the candidate for the paper's
acknowledged third error source"; page 6 resolves that source itself as
interferometer imperfection, bounded as minor. And it recorded a
"structural difference" — that Gobby split 1.6:1 while "this encoder AMZI
splits 50:50" — attributing residual to it. The chain passes
`split_ratio=SPLIT_RATIO=1.6` and reproduces the paper's ratio exactly, so
the header was blaming a difference that does not exist. Headers down from
37 of 47 lines to 23 of 33 in the CSV, 46 of 60 to 27 of 41 in the .tex,
with the (m)/(i) marking intact.

**The sweep was the regression check, and it passed.** 1.07e9 pulses
across nine distances to 122 km, every row bit-identical to the previous
artifact: same pulse count, same sifted count, same QBER to four decimals.
`bb84_time_bin.py` had changed substantially — `FiberRealization`, the
drift clock, blocked coefficient extraction — and bit-identity had only
been shown at 0 and 50 km on 300k pulses.

**Three faults in the eye diagram, and only one was known.** It was
unseeded, so every run differed. The X-cut MZM modulates Ey while
`LaserDriver` emits entirely on Ex, so the DFB eye was a flat line with
the modulator doing nothing, and `CWLaser` at 45 degrees left half its
power unmodulated as a floor capping extinction near 3 dB. And the drive
was an ideal rectangle, so edges were one sample wide and there was no eye
opening at all. Fixed with a seed, a polarisation controller implemented
as a unitary, and the fourth-order Bessel-Thomson drive at 0.75*baud the
eye-mask standards specify.

**Two reproducibility defects, found by running the harness.** The CW
laser figure was stale and that one was self-inflicted: regenerated
mid-way through fixing the eye, then committed alongside the finished
code. The calibration figure was never reproducible at all — three runs,
three hashes — because `SEED` reached `FiberRealization` while
`spad.detect` drew from the global `numpy` state that nothing had seeded.

**A property of this chain worth not relearning.** Its RNG desynchronises
on any float-level change: one flipped detection alters how many draws
follow, and every later pulse gets different numbers. Measured with a
deliberate nudge on the servo phase, 1e-15 relative moves the sifted count
about 13 % with no systematic trend — a reshuffled sample, not a different
answer. So a feature can be bit-identical when off only if the default
path's arithmetic is *literally* untouched; computing the same value a
different way is enough to break it.

**Scope corrections worth recording.** The figure-subtitle work was
believed to cover 21 scripts, then 8 by a regex, and was actually 3. Both
earlier counts were false negatives from a pattern that required the
newline in the first quoted segment, where these titles are built from
concatenated literals. Checked individually in the end.

Suite 332 -> 378 passed, 1 skipped. `run_all.py` 19 validators, all 18
non-Gobby passing in 51 minutes.

---

## 2026-08-18 — the Gobby chain gets a real fibre, and the fibre gets a clock

### Session: `FiberRealization` in `bb84_time_bin`, drift on both protocols, a new validator

| Change | Files | Rationale |
|---|---|---|
| Fibre impairments in the time-bin chain | `src/protocols/bb84_time_bin.py` | It modelled 122 km of fibre as one scalar multiply while `bb84_duplinskiy` and `bb84_test_dispersion` both built a `FiberRealization`. The justification was that both interfering paths share the operator so birefringence divides out -- true for the *balanced* topology, but the Gobby replication runs the polarisation-multiplexed one, where the arms leave on orthogonal polarisations and a rotation reaches them. |
| `FiberRealization.at(t)` and `drift_temperature_rate_C_s` | `src/channel/fiber.py` | `_J` was built once at construction and never rebuilt, so nothing could express a fibre changing between the moment Bob calibrates against it and the moment light travels through it. `calibration_temperature` / `calibration_bend_radius` gave a fixed two-state mismatch; what was missing is a mismatch that *grows*. |
| Blocked coefficient extraction | `src/protocols/bb84_time_bin.py`, `src/protocols/bb84_duplinskiy.py` | The closed form and the response table are each exact for one fibre state. Runs are now cut into `drift_blocks` pieces with the fibre held still inside each. A count rather than a pulse size, so raising `num_bits` for statistical power does not change the drift resolution. |
| `run_duration` on the polarisation chain | `src/protocols/bb84_duplinskiy.py` | Without it the pulse budget sets the simulated experiment length, so asking for tighter error bars quietly asks for a longer experiment. Same bug and same fix as the time-bin chain already had. |
| New validator | `analysis/validation/validate_gobby_impairments.py`, `run_all.py` (19 validators) | Three different required outcomes from one impairment model in one chain, plus the drift crossover. |
| Transcribed Gobby literals removed | `analysis/validation/validate_duplinskiy_birefringence.py` | It printed "no impairments 2.91 %, + birefringence 2.92 %" as the Gobby null, with no generating script -- the numbers came from the retired `val_system.py`. It now points at the validator that measures it. |

**The mechanism, because it is not the obvious one.** The sectional Jones
matrix is exactly SU(2): `|U00| = |U11|` and `arg(U00) = -arg(U11)`, both
to the float64 floor at 10/50/122 km across seeds. A rotation therefore
does exactly two things in the polarisation-multiplexed topology -- scale
**both** arms by the same `|U00|`, and shift their relative phase by
`2*arg(U11)`. There is no arm imbalance available to collapse the fringe.

Both halves of the paper's sentence follow from that. "Polarisation drift
reduces the bit rate" is the common amplitude; "but does not degrade the
QBER" holds because the phase is calibrated out, which their Bob does with
the piezo-driven fibre stretcher in his long arm; and the proviso, "provided
that the signal rate is significantly higher than the intrinsic error
rate", is the background claiming a larger share of a reduced signal.

Measured uncompensated against those two closed forms: QBER 86.77 +/- 1.60 %
at 10 km against a predicted 88.42, and 34.04 +/- 1.50 % at 50 km against
32.86; rate ratios 0.1068 and 0.3512 against 0.1017 and 0.3501. Cancelling
`2*arg(U11)` with the *modulator bias* parameter alone restores the QBER
while leaving the rate loss untouched, which is what establishes the
degeneracy.

**Drift is off by default in the replication, deliberately.** In this
topology fibre drift is degenerate with the interferometer arm drift
already modelled -- the fibre's phase enters `delta` exactly as an
arm-length offset would -- and `_bias_for_aggregate` has already assigned
the paper's full 3.3 % floor to modulator bias plus arm drift. A third
term would count one measurement twice. The paper also rules the fibre out
independently: the floor is flat to 65 km while a fibre effect grows with
length, classical fringe visibility is 99.96 % measured over the full
122 km link, and polarisation is reported stable for over 30 minutes there.

**A correction to an earlier diagnosis.** The QBER movement under an
uncompensated fibre was attributed to "arm imbalance at Bob's PBS". That
is wrong -- a unitary cannot imbalance the arms -- and the comments
carrying it have been replaced. The movement is the phase term alone.

**Bit-identity, which everything rests on.** Impairments off and drift off
reproduce the previous commit exactly at 0 and 50 km on both topologies:
same `qber`, `n_sifted`, `n_errors`. A static fibre with alignment on is an
exact null, `U_comp @ J = I` for any unitary, and is reported as arithmetic
rather than as evidence about impairments. Suite 350 passed, 1 skipped.

---

## 2026-08-17 — every validator in the harness; docstrings document the code, not its history

### Session: `run_all.py` 8 validators → 18, `--quick` output guards, docstring cleanup

| Change | Files | Rationale |
|---|---|---|
| Harness covers all 18 validators | `run_all.py` | It ran 8 of the 18 on disk. The other 10 already emit `[FAIL]` or exit non-zero, so they were written as regression checks and nothing invoked them. |
| Per-validator CLI in the roster | `run_all.py` | The newer validators do not share a command line with the original eight. None of the 10 takes `--seed`, which the old harness passed unconditionally; for them that is an argparse error and exit code 2, indistinguishable in the summary from a physics failure. |
| `--quick` output guard on the last three | `validate_dfb_drive.py`, `validate_dfb_duplinskiy.py`, `validate_gobby.py` | They wrote to the same paths whether or not the run was a smoke run, so a cheap run replaced a quotable artifact with an under-powered one. The other seven already had `_stem(quick)`. |
| `--allow-underpowered` no longer overwrites the published table | `validate_gobby.py` | It wrote straight over `val_gobby_table.tex`. The statistical-power guard refusing to write was the only thing keeping a smoke table out of the repository, and that flag exists to switch the guard off. A run that opts out of the guard is by definition not quotable, so it now writes to `--quick` names. |
| `*--quick.tex` ignored | `.gitignore` | The smoke block covered png, csv and md; Gobby's table is LaTeX. |
| Docstrings: documentation, not history | 13 modules under `src/`, 11 under `analysis/`, `run_all.py` | They had become a record of how the code got here: hypotheses tested and rejected, measurements from those tests, values that changed and why, corrections to earlier versions. |
| Generic components no longer discuss replications | `phase_modulator.py`, `interferometer.py`, `optics.py`, `bb84_time_bin.py` | A hardware model may cite a paper for a formula or a number, but not explain itself in terms of whether one experiment used it. `phase_modulator` had 5 Gobby mentions, none in a citation block. |
| References to the issue log removed | 26 files | `opto-sim-issues-and-fixes.md` is untracked and never pushed, so anyone cloning the repository read citations to a document they cannot open. 138 lines. |
| Measured sectional-birefringence behaviour recorded | `src/channel/fiber.py` | It was not written down anywhere. |

**The rule applied to docstrings.** A docstring says what the thing is:
what it models, units, defaults, the physics that constrains it, and its
citations. It may carry a standing caution about correct use, such as
that a value is a device specification rather than a fitting knob, or
that two mechanisms must not both be applied because they double-count
one measured aggregate. The test for a sentence is whether it would still
be written if the code had been right the first time.

Replication scripts are the opposite case and keep their parameter
justifications. `AFTERPULSE_PROB = 0.0` in `validate_gobby.py` still
explains itself from the paper's own Fig. 3, its stated error
probability, and its list of three mechanisms. What went was the story of
the value having been 0.05 and being changed.

**Two harness findings, both from running it.** Gobby's
`--gobby-bits 200000` default predates the link-budget correction, which
cut the signal 3.68x; at that budget the sweep now yields 318 sifted bits
at 0 km and zero at 122 km, and the power guard correctly refuses to
write. Replaced by `--gobby-target-sifted`, since a flat budget spans a
200x sifted fraction across the sweep and either starves 122 km or wastes
hours at 0 km. Separately, threading `quick` into
`validate_dfb_duplinskiy` missed that the PNG write sits inside
`_figure()`, which did not take the argument: `NameError` at write time,
caught on the harness's first run. That file had no automated caller
before.

**Measured sectional birefringence** (500 realizations per length, Haar
SU(2) reference mean rotation angle 2.217 rad): the model is already
indistinguishable from a uniformly random SU(2) by about **ten metres**,
where the ensemble mean Stokes vector has fallen to the finite-sample
floor of 1/sqrt(500). Unitarity error is machine epsilon times the
section count, so **1.9e-10 at 50,000 km**, further than any cable on the
planet. The build is one-time per realization, 2 ms at 122 km, because
`FiberRealization` freezes the matrix at construction.

**Verification.** 332 passed, 1 skipped throughout. All 173 bracketed
citations preserved, along with the Duplinskiy paper's own section
references and the Coldren, Agrawal, Keiser and Wiley section numbers.
`run_all.py` 18 passed, 0 failed, 0 skipped, with no tracked artifact
modified by a quick run.

---

## 2026-08-14 — DFB device model: validated, merged to `main`, driving Duplinskiy

### Session: the `dfblaser_work` branch, 11 commits, fast-forwarded to `main`

| Change | Files | Rationale |
|---|---|---|
| **One TE mode, not two** | `src/lasers/dfblaser.py` | The device tracked two field pairs sharing every coefficient and one carrier reservoir, differing only in the noise seed — one mode duplicated, not TE and TM. Both lased and split the output evenly (PER **0.01 dB**), and with independent seeds the relative phase was random, so the source emitted *unpolarised* light. Kim's Eq. (1) has a single modal function. Pair deleted; total power moved 5.375 → 5.225 mW at 100 mA, because gain clamping fixes the total however many modes share it. |
| `Laser` → `DFBLaser`, exports wired | `src/lasers/__init__.py`, `laser_driver.py`, `analysis/validation/validate_dfb_reflection.py` | Matches `CWLaser` in the same package; `__init__` exported only `CWLaser` before. |
| **`LaserDriver.sample_field(dt, n)`** | `src/lasers/laser_driver.py` | `(n, 2)` complex `[Ex, Ey]`, mirroring `CWLaser.sample_field`. Discards a 40 ns settle, removes the lasing mode's **549–591 GHz** offset from the Bragg reference before resampling (a 2 ps grid would alias it), decimates by averaging, and places the TE axis with a Jones vector now shared with `CWLaser`. |
| `cw` / `gain_switched` drive modes | `src/lasers/laser_driver.py` | "pulsed" described the output rather than what is done to the device. RIN, chirp and phase statistics are documented as device *outputs*, not driver inputs. |
| Convergence guard → `UserWarning` | `src/lasers/dfblaser.py` | It is a construction-time configuration choice; `RuntimeWarning` is what numpy raises for overflow, so filtering that would have trapped this instead. |
| **Chirp measurement corrected** | `analysis/validation/validate_dfb_drive.py` | `np.gradient` of the unwrapped phase at the 0.4933 ps step reported **297 GHz** — noise. It alternated sign every sample and spanned 1123 GHz, above the device Nyquist of 1013 GHz. Replaced with the amplitude-weighted single-lag autocorrelation: **70 GHz**, stable across analysis windows. |
| Rayleigh test for pulse-phase randomness | `analysis/validation/validate_dfb_drive.py` | A fixed cut on \|⟨e^{iφ}⟩\| is wrong: it falls as ~0.886/√N, so a 12 ns run tripped 0.5 at 0.516 on 22 pulses where 30 ns gives 0.086 on 58. `Z = N R²` is Exp(1) under uniformity at any N. |
| **Resonance-scaling claim withdrawn** | `validate_dfb_drive.py`, §30.4 | The script reported the intensity-noise resonance failing sqrt(I−I_th). At most currents there is no resolvable line to scale — five peaks within a factor of 1.5 at 100 mA — so `argmax` hopped between near-equals; a centroid does not scale either (r = +0.11). Withdrawn as **unmeasurable**, not as a model failure. |
| Replaced by RIN vs power | `validate_dfb_drive.py --rin-scaling` | `RIN ~ P^-1.20`, r = −0.89 against the textbook 1/P, over an order of magnitude in power, nothing fitted. Absolute level **−141 to −156 dB/Hz**, inside the −140 to −160 range `CWLaser` has to be *told*. |
| **`source_field` / `pulse_energy_factors`** | `src/protocols/bb84_duplinskiy.py` | The chain injected a flat analytic field, no laser at all. Both default to the old behaviour and are **bit-identical** to the frozen §27.1 baseline; factors are cycled in index order, never sampled, so no RNG is consumed. |
| Stokes + Poincaré validation | `analysis/validation/validate_dfb_duplinskiy.py` (new) | Four BB84 states through every stage, the 8-outcome control, and QBER across three sources. |
| Two dead files preserved then removed | `src/lasers/dfblaser_v2.py`, `drive.py` | Neither was in git history *anywhere* — always untracked — so deletion would have been permanent. Committed first (`0ce7c46`) so removal is reversible. |
| Three questions closed by decision | §30.4 | Efficiency gap (parameter-level, only lever is `alpha`, which G9 forbids adjusting), gain-bandwidth filter (wider than the simulation band), facet bifurcation (**spatial hole burning** — carriers symmetric at 100 mA, 9.2 % asymmetric at 120; correct behaviour for a uniform grating with AR facets). |

**Key result — a polarisation-encoding chain is blind to everything a
source adds except pulse energy.** `sample_field` returns one complex
amplitude times a *fixed* Jones vector, so both components carry the same
amplitude, and normalised Stokes parameters — which depend only on Ey/Ex —
cannot see RIN, phase noise or chirp. Measured to nine decimals at 80, 100
and 200 ps drive, in CW mode, and identically for `CWLaser`. The four BB84
states come out exact with **DOP = 1.000000**: D [0,+1,0], A [0,−1,0],
R [0,0,+1], L [0,0,−1], matching the paper's Eqs. (4)–(5). Same
cancellation as §23.2 (linewidth) and §26.6 (CD/PMD/birefringence),
reached from a third direction.

QBER, ≥3122 sifted per cell, combined across 0/10/50 km: **+0.082 ± 0.202 pp**
(200 ps drive) and **+0.006 ± 0.200 pp** (100 ps). A null, predicted before
the run, bounded below **0.40 pp at 2σ**. Duplinskiy states no optical
pulse width — registered as **A9** — but the device bounds it to 100–250 ps
and both ends agree to 0.38σ, so the unstated parameter is not load-bearing
for the result.

**Still open on the DFB:** `tests/test_dfblaser.py` is parked behind a
module skip. It targeted an API that never existed and five of its eleven
tests exercised the gain-bandwidth filter, which is deliberately absent.
The three validation scripts cover the model end to end meanwhile.

332 tests pass, 1 skipped.

---

## 2026-08-11 — A4 settled by algebra: `circular_analyser` ≡ the paper's PC3+PBS readout

### Session: VALID-1, register entry A4 closed — zero simulation runs

| Change | Files | Rationale |
|---|---|---|
| **A4 register row → SETTLED (TRUE)** | `opto-sim-issues-and-fixes.md` §28.6 | QWP+PBS (ours) ≡ HWP(22.5°)+PBS (paper) on the four BB84 states, proven by Jones calculus, not asserted. |
| **Derivation written** | §28.6a (new) | Paper quotes (eqs. 4–5, PM2/PC3), both composite matrices, the state-for-state port table, the operator identity. |
| **Downstream consequences updated (G8)** | §27.5, §28.8 item 2, §29.2, §29.8 | "Evidently equivalent" → proven; "A4 unsettled" caveats → discharged with reasons revised. |
| **Numeric spot-check** | deterministic matrix products (no simulation) | `J_circ = diag(1,−i)·H·diag(1,−i)` exactly; `P(port\|s,Δφ) = P(port\|s,Δφ−π/2)` for all eight (state, phase) pairs. |

**Key result:** the equivalence is structural, not coincidental. The
circular analyser is the paper's HWP(22.5°) sandwiched between two λ/4
phase plates — exactly the phase the paper assigns to Bob's PM2. The chain
folds the λ/4 into a fixed element (QWP+PBS); the paper folds it into
`pm_bob`. Port statistics agree state-for-state: chain-X ≡ paper-"linear"
exactly, chain-C ≡ paper-"circular" up to a free port exchange. Matched →
deterministic port, mismatched → ½/½, in both implementations.

**Nothing changed:** `circular_analyser` untouched (live in four
protocols), no numerical result moved, no code edits. The "coincidental
for the tested states" alternative reading is refuted.

332 tests pass (unchanged).

---

## 2026-08-09 — the Duplinskiy chain brought to standard; Gobby closed out

### Session: DUPL-1 — a polarisation chain that can carry the impairment models

| Change | Files | Rationale |
|---|---|---|
| **Five O(N) per-pulse lists removed** | `src/protocols/bb84_duplinskiy.py` | Same defect the Gobby chain carried until §17.5; several GB at sweep-scale pulse counts. Bit-identical. |
| **8-outcome precompute (PERF-2 analogue)** | same | The response is deterministic given `(alice_basis, alice_bit, bob_basis)` because the fibre Jones matrix is sampled once per run (ROOT-1). **29,938 → 386,245 pulses/s, 12.9x**, bit-identical, negative control fires. |
| **Detector parameterised** | same | QE / dead time / DCR / afterpulse were literals. All four defaults are stated verbatim in the paper. |
| **`cd` / `pmd` exposed** | same | Were hardcoded `False`. Exposed — and found **inert**, see below. |
| **Gobby closeout** | `analysis/val_gobby/validate_gobby.py` | χ²/dof relabelled; four stale "still open" claims corrected against the code. |

**Key results:**

- **The 8-outcome table is exact, and was verified to be legitimate first.** No stage in the field chain consumes randomness — `pm.modulate`, `fibre.apply`, `optics.voa`, `circular_analyser` all leave both RNG streams untouched, and `fibre.apply` is repeatable. Without that, precomputing would shift the stream and silently change every downstream detector draw. Bit-identity holds at 0/10/50 km × compensate on/off; the negative control (swapped analyser outputs) moves the 10 km QBER 0.0274 → 0.9833, so the gate can fail.
- **`cd` and `pmd` are inert on this chain by construction.** It builds its field as `np.ones((1,2))` — a **single time sample**. Measured through the same fibre: `max|dE|/|E|` = **0.000e+00** at one sample, **5.945e-01** at 4096. They are not nulls, they are no-ops. A polarisation impairment table built today could honestly carry **birefringence only**.
- **We reproduce Duplinskiy's 2 % — by the wrong mechanism.** The paper names **three** contributions: a ~1 % afterpulse floor, a finite instrumental extinction (calibration goal 3, ">98 %"), and drift/recalibration bringing the average to 2 % over an 80 % duty cycle. This model has **one**: with afterpulsing off the QBER is exactly **0.00 %**, and at the datasheet `p_ap = 0.05` it produces **2.07 ± 0.40 %**. **Occurrence #7** of the recurring pattern — one mechanism doing three mechanisms' work and landing on the right total. Caught before it was written up as a successful replication.
- **The old "0.98 % vs 2 %" gap was underpowered noise.** At 400k pulses the 50 km point carries ~140 sifted bits (σ ≈ 0.71 pp); at 4e6 pulses it reads 2.07 ± 0.40 % on 1,258 sifted.
- **Two open questions answered from the source.** Bob's detection is PC2 → PM2 → PC3 → **PBS**, with PC3 acting "similar to a half-wave plate… rotating polarization by 45°"; this model uses QWP+PBS via `circular_analyser`. Structurally different, composite evidently equivalent for these states — **reported, not changed**, since `circular_analyser` is live in four protocols. And `U_comp` is computed once and applied unchanged, where the paper is explicit that φ₁/φ₃ "are not guaranteed to remain stable"; the full cycle (floor ~1 %, ceiling ~5 %, 80 % duty) is citable but **not implemented**, pending the afterpulse question.

**Corrections to earlier claims in this document, each named rather than replaced:**

- **χ²/dof = 1.95 is an upper bound, not an estimate.** It divides our residuals by *our* error bars alone, granting Gobby's published values zero uncertainty. A σ of **0.329 pp** on their side brings it to 1.03 — comparable to our own mean σ of 0.3375 pp, against values quoted to one decimal with no error bars. The earlier reading ("~1.4× what sampling noise explains, so a small residual systematic") is a claim about our model that the number does not support, and is **withdrawn**.
- **The afterpulse question does not touch Gobby.** An earlier statement said fixing it "would move the Gobby results committed at v0.3.1". It cannot: `AFTERPULSE_PROB = 0.0`, so the dead-time path is never reached. Impact is Duplinskiy-only, where it is 100 % of the error budget.
- **The "2x" afterpulse figure is soft.** It compares our measurement to a one-sentence calculation whose word "mainly" is unquantified — ~2x at 1.0 %, ~3x at 0.7 %. Only one side is measured. Recorded as a discrepancy to characterise, not a defect; the counter-argument (that `p_ap = 0.05` is a datasheet figure already net of dead-time suppression) is strong.
- **Two literature entries retracted after reading the articles**: the "aerial ~1 ms / underground tens of seconds" drift figures are not in the cited paper, and the DPS-QKD arXiv identifier was guessed and wrong. **Verified by reading**: Zelmon, Small & Jundt (1997) — at 1.55 µm n_o = 2.2111, n_e = 2.1376; our defaults are good to 0.5 % on the indices and 1.5 % on V_π, but Δn is off 18 % (a difference of rounded near-equal numbers) and is inert today.

**Measured costs, recorded so the scope decision is not re-litigated from memory:** one full field-chain evaluation is 0.048 ms at one sample, 0.357 ms at 2048; the 8-outcome table pays it **8 times per run**, so a time-resolved field costs **0.004 %** of a 78 s sweep row. **Compute is not the constraint on any outstanding work** — an earlier draft costing the time-resolved field as "medium" was wrong.

**Scope decided:** build extinction and the birefringence sweep (the latter would be the first QBER-level validation of any fibre impairment model in this simulator); skip diode chirp, modulator crystal PMD and the PC tuning algorithm — the first two are compensated in source, so a faithful model outputs zero. Recommended next protocol is **DPS-QKD**, which reuses `AsymmetricMZI` and is aligned with deployed schemes.

332 tests pass.

---

## 2026-08-09 — the mechanisms Gobby names, at the components that cause them

### Session: GOBBY-7e — the record corrected, and the last open items closed

| Change | Files | Rationale |
|---|---|---|
| **Pattern count corrected** | `opto-sim-issues-and-fixes.md` §25.4, `CHANGELOG.md` | The "fourth occurrence... now five" claim was wrong in both halves. |
| **Sweep variance statistics** | §25.5, `CHANGELOG.md` | 0.414 pp compared against an interpolated value; χ²/dof was never stated. |
| **65 km cross-check reading settled** | `analysis/val_gobby/validate_gobby.py` | The flattering interpretation was available and is rejected in writing. |
| **`coupler_combine` → real convention** | `src/channel/optics.py`, `tests/test_optics.py` | Deferred since §21.3; the module now claims the real form in its own header. |
| **OPEN-5 provenance** | `analysis/val_system_scenarios.py` | Header stamped a pulse count no run used. |
| **OPEN-3 re-run at power** | same, + `val_system/*` | The blocker was a timing estimate stale by 4×. |

**Key results:**

- **Two of my own statements corrected, both named rather than replaced.** The pattern list claimed five occurrences with two added here; in fact §19.5's item 3 *is* the afterpulse case, so this work **closed** it rather than adding it, and §17.2(a) and §18.4 had been missed. Six total, one added and one closed. And the sweep residual was 0.4325 pp, not 0.414 — the earlier figure compared 100 km against the interpolated 5.9 % instead of Gobby's raw 101 km measurement of 6.0 %, which also moves that residual to −0.69 pp (2.0σ).
- **χ²/dof = 1.95 now stated.** Mean signed residual +0.0625 pp (essentially unbiased), variance 0.3250 pp², RMS 0.4977 pp against a mean Monte-Carlo σ of 0.3375 pp. The scatter is ~1.4× what sampling noise alone explains, so a small residual systematic sits on top of it. The replication is good — unbiased, sub-half-pp, nothing beyond 2σ — but it is **not a pure-noise fit** and had been presented as though it were.
- **The 65 km miss is real, and the escape route was refused.** Reading the stated <0.4 % as covering the *dark term alone* would pass comfortably (0.184 % at 65 km, crossing only at 82 km), but the paper's sentence names "dark counts **and** stray light", so the bound covers `P_E` in full. Taking the narrower reading would have been choosing the interpretation that flatters the model; it is recorded as considered and rejected. The miss stands at 0.485 %, over by 21 % at the endpoint and holding to 60.8 km. The bound would allow `P_E ≤ 7.00e-07` against the **8.50e-07 the paper itself measures** — their stated `P_E` and their stated bound are mildly inconsistent within this model. Same direction and size as the χ²/dof systematic; plausibly the same thing. `P_E` unchanged.
- **`coupler_combine` was verified unreachable before being changed.** `interferometer.py` imports only `coupler_split`; `bb84_time_bin.py` only `pbs`, `pbc`, `voa`; nothing under `src/` names it outside its own definition. A convention flip is what inverted the Y basis in GOBBY-4, so this was established rather than assumed. Two tests added beyond the matrix assertions: real input → real output, and interference as **cos** not **sin**, the property whose absence produced a flat ~50 % QBER in GOBBY-2 §19.
- **The second stale timing estimate found in this project, again off several-fold in the same direction.** `val_system_scenarios.py` warned "~35.8 us/pulse... budget ~5 h for eight rows"; re-measured at **118,000 pulses/s (8.5 us/pulse)** with a 4.8e-5 sifted fraction, i.e. **~1.2 h for eight rows**. The Gobby sweep was likewise believed to cost 4.5 h and measured 0.9 h. Both notes corrected in place with an instruction to re-measure rather than trust them.
- **OPEN-3 closed, with a caveat that matters more than the closure.** Re-run at `--target-sifted 3000` (648e6 pulses, 100 km), rows now carry 3,589–3,915 sifted bits against the 86 that prompted the issue; σ ≈ 0.27 pp, so any effect above ~0.8 pp would show. The four middle rows are **still bit-identical**, which needed a different check than statistics: identity to the bit has two indistinguishable causes, real invariance or an impairment that never ran — the same trap as GOBBY-6's vacuous negative control. Resolved twice over: each impairment demonstrably perturbs the field while conserving power (birefringence `max|dE|/|E|` = **1.679**, CD **0.595**, PMD **0.010**), and the CD code-path check moves the sifted rate **+31.0 %, 11.7σ**. `Full chain` being bit-identical to `+ Visibility` is then *predicted*, not suspicious.
- **What that table is not.** It demonstrates that time-bin/phase encoding is immune to CD, PMD and quasi-static birefringence — which is what the scheme exists to do, and what OPEN-3 said was asserted but undemonstrated. It is **not** a QBER-level validation of the impairment models: the observable is constitutionally blind to all three, so nothing here can show they are *correct*, only that they are live and that the scheme is immune. That requires polarisation encoding, where an uncompensated SU(2) rotation maps onto the bit. Recorded as the successor, not claimed as done.
- 317 tests pass.

---

### Session: GOBBY-7d — the nine-point sweep, and the two defects it exposed

| Change | Files | Rationale |
|---|---|---|
| **`run_duration` — drift clock separated from the detector clock** | `src/protocols/bb84_time_bin.py` | Drift advanced as `pulse_idx / repetition_rate`, so the pulse budget chosen for *statistical power* silently set the simulated *experiment duration*. Detector clock keeps true `1/f` spacing (dead time and afterpulsing are defined against real elapsed time); only drift moves to the declared duration. `None` is bit-identical to the old behaviour. |
| **Bias solved jointly with drift** | `analysis/val_gobby/validate_gobby.py` | `arccos(1 - 2*E_MOD)` is the bias-**only** solution; applying drift on top counted it twice. `_bias_for_aggregate()` solves the time-average instead. |
| **Nine-point sweep re-run and artifacts regenerated** | `analysis/val_gobby/*` | OPEN-2 (P0) closed. |
| **Drift-clock and joint-solve invariants** | `tests/test_polmux_interferometer.py` | Accumulated drift at the final pulse is 6.00° at any budget from 1e3 to 1e9 pulses. |

**Key results:**

- **The sweep did its job by failing first.** Run at the corrected model it produced **13.52 % at 122 km against Gobby's 8.9 %**, and every point matched a drift-aware prediction to within statistics (122 km predicted 13.54 %). The 122 km point needs 1e9 pulses = **500 s at 2 MHz**, against the paper's stated *two-minute* transfer — accumulating 25° of drift where theirs accumulates 6, and inflating the effective modulation error from 3.31 % to 8.60 %. The whole 4.6 pp excess was that.
- **The flaw was introduced in GOBBY-7**, where a run-duration knob was considered during planning and dropped as overreach. Recorded as such, not presented as a discovery. **More pulses must mean a better estimate of the same experiment, not a longer experiment.**
- **A second defect surfaced once the duration was set to the paper's 120 s**: the 0 km floor read 4.30 % against 3.3 %. Gobby attribute the floor to bias inaccuracy "**as well as phase drift during the experiment**" — 3.3 % is the *aggregate of both*, so deriving the bias from the full 3.3 % and then adding drift double-counts. Solved jointly: **d0 = 17.864°**, drifting to 23.864°, time-averaging to exactly **3.300 %** — a ramp *centred on* the 20.93° the naive reading assigns.
- **The recurring pattern, counted correctly.** An earlier draft of this entry said "fifth occurrence" and credited this session with two. Both were wrong: §19.5's item 3 *is* the afterpulse case, so this session **closed** it rather than adding it, and two earlier instances (§17.2(a), §18.4) were missed. The list runs to **six** — 3.3 % floor paid twice by V=0.934 and afterpulsing (§17.2a); 2.43 % as "optical misalignment", actually afterpulsing (BLOCK-3); 1,788 Hz as "dark count rate", actually a lumped term (§18.1); the 3.3 % floor read as interferometer visibility when visibility is an *output* (§18.4); `afterpulse_prob = 0.05` standing in for modulation error (§19.5, closed here); and the bias absorbing drift's share (§25.4, new here). `MU_EFF = 0.0793` is deliberately excluded — fitting to a target is a different failure mode. The recurrence is the point: it is the characteristic failure mode of replicating an aggregate whose components are *named but not individually quantified*.
- **The sweep, measured**, against Gobby's four real measurements (§18.9 established the other five rows compare against nearest-matched values and are not evidence):

  | km | this work | Gobby | residual | sigma |
  |---|---|---|---|---|
  | 4.4 | 3.25 ± 0.29 % | 3.3 % | −0.05 | 0.2σ |
  | 65 | 3.69 ± 0.28 % | 3.3 % | +0.39 | 1.4σ |
  | 101 | 5.31 ± 0.35 % | 6.0 % | −0.69 | 2.0σ |
  | **122** | **9.50 ± 0.43 %** | **8.9 %** | **+0.60** | **1.4σ** |

  **Mean |residual| 0.4325 pp**, mean signed residual +0.0625 pp (essentially unbiased), **variance 0.3250 pp²** (sd 0.5701), RMS 0.4977 pp, mean MC σ 0.3375 pp, **χ²/dof 1.95**. *Correction:* an earlier statement gave 0.414 pp and a −0.62 residual at 101 km; that compared against the interpolated 5.9 % rather than Gobby's raw 6.0 % measurement.
- **χ²/dof = 1.95 is stated, not buried.** The scatter is ~1.4× what Monte-Carlo statistics alone explain, so a small residual systematic sits on top of the sampling noise. The replication is unbiased in sign and sub-half-pp in magnitude with nothing beyond 2σ — but it is not a pure-noise fit and should not be presented as one. The systematic is unidentified; documented candidates are the ~2 % first-order/MC signal disagreement and the ≤0.050 % device-visibility gap (§24.5).
- **Statistical power holds.** Every row exceeds the 3,000-sifted target (minimum 3,384) and **nothing is clipped** — 122 km used 663e6 pulses against the 1e9 ceiling, where the pre-fix run hit the ceiling exactly. The committed artifact is no longer a 20,000-pulse smoke run with zero sifted bits at 122 km.
- **Nothing is fitted anywhere in the chain.** `P_E` measured, `MU_EFF` derived from the 1.6:1 split, `e_mod` from Fig. 3, drift rate and transfer duration from the text, `afterpulse_prob = 0` from §19.5's three grounds.
- 315 tests pass.

---

### Session: GOBBY-7c — a SPAD detection-probability defect, and the afterpulse double-count

| Change | Files | Rationale |
|---|---|---|
| **`spad.detect` Poisson form corrected** | `src/detectors/spad.py` | Computed `eta*(1 - exp(-mu))` — "at least one photon arrives, *then* one coin flip at eta". Photons are detected independently, so the detected count is Poisson(`eta*mu`) and the exact form is `1 - exp(-eta*mu)`. **Affects every protocol in the repository.** `apd.py` checked and correct — it applies `qe` before the Poisson draw. |
| **`AFTERPULSE_PROB = 0.0`** | `analysis/val_gobby/validate_gobby.py` | §19.5 concluded this from the paper long ago but the code never carried it. With `e_mod` implemented since GOBBY-6, keeping 0.05 double-counted the floor. |
| **Stale documentation reasserted** | `opto-sim-issues-and-fixes.md` | Each claim re-derived against current behaviour rather than the heading rewritten. |
| **Poisson-form and saturation tests** | `tests/test_spad.py`, `tests/test_analytic_gobby.py` | Both defects that drove `S_mc/S_analytic` off unity are now guarded. |

**Key results:**

- **The unattributed residual from GOBBY-7b is found, and it was a real physics defect — not a Gobby artifact.** Chasing it *before* committing an hour to the nine-point sweep is what surfaced it. The old form predicted 0.9625 of the analytic value against the measured **0.9533 ± 0.0123** — 0.75σ, which identified it.
- **The error grows with intensity**, so it was benign where the simulator has mostly been used and badly wrong outside it: −0.05% at `mu` = 0.001, −3.58% at Gobby's operating point, **−20.4% at `mu` = 0.5, −54.8% at `mu` = 2.0**. Classical visibility measurements, alignment runs and any high-`mu` study were substantially wrong.
- **Agreement with the Monte Carlo after both fixes**: 0 km **1.019 ± 0.012**, 40 km 1.011 ± 0.025, 65 km 0.985 ± 0.031 — ~2%, consistent with unity across the range, from 0.954 / 1.012 / 0.994 before. One RNG draw either way, so stream alignment is unchanged and the difference is purely physical.
- **The afterpulse default was double-counting the floor.** §19.5 had already established `afterpulse_prob = 0` on three sourced grounds — Fig. 3's dashed curve starting at ~0, the stated `P_e` having no afterpulse term, and the closing summary naming three mechanisms without it — but the code stayed at the ID230 0.05 because until GOBBY-6 afterpulsing was the only thing supplying a floor. On defaults the 0 km QBER read **5.806%** against Gobby's 3.3%; the 2.5 pp excess is exactly the `p_ap/2` §19.5 predicts. Corrected: **3.088 ± 0.247%**, 0.86σ from the stated value.
- **Ownership was checked, and nothing moved.** `afterpulse_prob` *is* a SPAD parameter (`spad.py`, default 0.05 = ID230); the constant in `validate_gobby` is the replication-level override, which is where a claim about Gobby's apparatus belongs rather than one about SPADs generally.
- **Numbers previously reported are restated with the correction named**, not quietly replaced — floor 3.253 → 3.088%, jitter control 3.579 → 3.369%, `S_mc/S_analytic` 0.954 → 1.019. The `T_INT` pass corrected the analytic helper and the Monte Carlo was bit-identical across it; this pass corrects the physical chain, so movement was expected.
- 308 tests pass.

**Documentation reasserted rather than patched.** §19's heading still read "⚠️ OPEN (top priority)" though GOBBY-2 closed several passes ago; §20.6's "Modulation error is not implemented" was false since GOBBY-6 — and the concern it raised, that afterpulsing and `e_mod` were two mechanisms producing a similar number, turned out to be live and is what this pass acted on. §20.6's arm-equalisation bullet was re-checked and **left standing, still true**.

---

### Session: GOBBY-7b — `signal_click_prob()` was missing `T_INT`; jitter demoted

| Change | Files | Rationale |
|---|---|---|
| **`T_INT` applied in `signal_click_prob()`** | `analysis/val_gobby/validate_gobby.py` | It computed `MU * T_link * ETA_BOB` and never applied the polarisation-multiplexed interferometer transmission `2/(1+r)` = 0.769231, putting it **41% above** what the chain delivers. Imported from `src/analytic/gobby_model.py` rather than restated, so one definition and the existing traceability test keeps covering it. Propagates to `model_qber()` and `predicted_visibility()`, which is the point. |
| **Jitter demoted to a negative result** | `src/channel/phase_modulator.py`, `src/protocols/bb84_time_bin.py`, `validate_gobby.py` | `phase_noise_rad` *is* the drive-voltage-noise knob. Never a default anywhere; `PHASE_NOISE_RAD` marked diagnostic-only. Kept, not deleted — modulators with genuinely random shot-to-shot error exist. |
| **Saturation bound + regression guards** | `tests/test_analytic_gobby.py` | The reason a signal error could never explode the QBER, plus a guard so the `T_INT` omission cannot silently return. |

**Key results:**

- **The 41% was a real omission, not a measurement artifact**, and decomposes cleanly at 0 km: **`T_INT` 23.1 pp**, detector dead time 3.1 pp, unattributed residual ~4.7 pp. `0.769231 × 0.9533 × 0.9689 = 0.7105` closes against the measurement. **Restoring `T_INT` is not fitting** — no value is chosen; one the repository already derives from the stated 1.6:1 split is used where it belongs.
- **Three independent checks improved together**, which is the evidence it was a defect rather than a comparison artifact — no single number was targeted:
  - mean |residual| over Gobby's four measured points **0.378 → 0.221 pp** (122 km −0.89 → +0.36, 100 km −0.80 → −0.28, 80 km −0.38 → −0.17);
  - visibility @ 122 km **0.9058 → 0.8809** against the stated **0.884**, i.e. +0.022 off → −0.003 off — the sharpest test, the paper giving a value there rather than a bound;
  - visibility @ 65 km 0.9903, still above the stated >0.99.
- **The residual is bounded, not attributed.** After `T_INT` the form is a few percent optimistic — 0.954 ± 0.012 at 0 km, 1.012 ± 0.025 at 40 km, 0.994 ± 0.031 at 65 km — with **no significant distance dependence** (0 and 65 km differ by 1.2σ). Ruled out: gate width, pulse width, the `T_INT` machinery itself. Recorded as a quantified docstring bound rather than guessed at. Dead time is deliberately not carried in a first-order optical budget.
- **The 65 km cross-check still fails and is still not tuned**: 0.485% against the stated <0.4%, over by 21%. The bound would allow `P_E <= 7.00e-07` against the 8.50e-07 carried. A regression test asserts it keeps missing, so a later change cannot launder it into a pass.
- **The Monte Carlo did not move** — the check that matters, since this corrected the analytic helper and not the physical chain. Floor still 3.253 ± 0.258% at 0 km; jitter control still 3.579 ± 0.267%.
- **Why the defect stayed survivable:** `P_e/(S + 2*P_e)` saturates at 1/2 as `S → 0`, so a signal error shifts the share within a bounded range (+0.008 pp at 0 km to +1.675 pp at 122 km) and can never blow it up. Now asserted in `TestErroneousCountSaturation` rather than argued in prose.
- 303 tests pass (276 + 27).

**Retracted:** an interim estimate in this session scaled S by a constant `K = 0.7105` at every distance and reported 0.254 pp with a 0.525% miss. Wrong in method — dead time is distance-dependent and vanishes at range, so a constant is not the right correction. The `T_INT` fix is both more defensible and better: 0.221 pp, 0.485%.

---

### Session: GOBBY-7 — arm-length drift implemented, bias in volts, two prose cross-checks

| Change | Files | Rationale |
|---|---|---|
| **`phase_drift_rad_s` on `AsymmetricMZI`** | `src/channel/interferometer.py` | Gobby *measure* interferometer phase drift at <0.05°/s and the code had no such parameter. It is **arm-length** drift — "variations in the relative lengths of the two arms" — so it is a property of the interferometer. New `arm_phase_offset(t)` owns the law; `modulate` gains `t=0.0`. |
| **`bias_offset_v` on `PhaseModulator`** | `src/channel/phase_modulator.py` | "Slight inaccuracies of the phase modulator biases" is literally a bias voltage. Converted through the existing crystal-derived `V_pi` (3.8826 V), so 20.93° reads as 451.5 mV. Supplying it *and* `phase_error_rad` raises — one mechanism, two unit systems. No signature change. |
| **Protocol queries, never restates** | `src/protocols/bb84_time_bin.py` | Supplies only its clock, reusing the `pulse_idx / repetition_rate` the SPAD calls already computed, and subtracts the component's own `t=0` value since the PERF-2 coefficients were extracted at rest. Neither law is written down in the protocol. |
| **Gobby default: jitter → static bias** | `analysis/val_gobby/validate_gobby.py` | The old default was justified by "drive noise is random per pulse" — the one mechanism the arithmetic rules out. The paper says "biases" and "drift"; neither is per-pulse noise. |
| **Two cross-checks from the prose** | same | `print_paper_cross_checks` reports both stated bounds on every run, pass or fail. |
| **`bias_offset_v` wired into Duplinskiy** | `src/protocols/bb84_duplinskiy.py` | Overfitting guard: the test suite requires a **non-Gobby** protocol to exercise the parameter. |

**Key results:**

- **Drive-voltage noise is ruled out as the mechanism, and the negative result is recorded.** Reproducing the 3.3% floor would need 451.5 mV = **11.6% of V_pi**; realistic drive electronics sit at 0.1–1%, worth 0.0002–0.025% QBER.
- **The framing was corrected, not patched.** There is no more-fundamental number behind the 3.3% — it is a measured aggregate from a hand-biased apparatus. The bias magnitude being derived from it is **the correct treatment of an experimental result**, not a shortcoming. GOBBY-6 §23.3's caution is rewritten accordingly; what survives is that it remains not a measurement of their modulator.
- **Drift verified end-to-end at measurable amplitude.** 0.091% is too small to resolve without ~1e8 pulses, so the mechanism was tested by sweeping the phase 0 → 90° across a 2e6-pulse run: **18.342 ± 0.833%** against the analytic time-average **18.169%** (0.21σ), and exactly 0.000% with drift off. The cited rate checks arithmetically at 6.000° over their 2-minute transfer, 0.075° over a 1.5 s run. **The rate is passed unscaled** — rescaling it to fit a shorter run would turn a cited constant into a fitted one.
- **The floor holds under the corrected default.** Static bias measures **3.253 ± 0.258%** at 0 km against Gobby's stated **3.300%** (0.18σ); the old jitter default gave 3.579 ± 0.267%. Both within ~1σ, so this changes provenance rather than the result.
- **One prose cross-check passes, one fails — and the failing one is reported, not fixed.** Device visibility >99.9% is satisfied (balanced arms give 1.0). Erroneous counts at 65 km **miss**: 0.525% against the stated <0.4%, over by 31%. Tuning `P_E` to pass would re-fit exactly what GOBBY-1 unfitted.
- **A larger discrepancy surfaced while checking it.** `signal_click_prob()` runs **~41% above** what the chain delivers (0 km: 4.50e-03 vs 3.20e-03 measured), which is what flips the 65 km verdict. It also underpins `model_qber()` and `predicted_visibility()`, i.e. every analytic column. **Deliberately untouched** — a separate change with its own verification.
- **LiNbO3 bias drift deferred on a physical criterion, not for convenience.** Its relaxation runs minutes to days, so over runs of this length it is indistinguishable from a constant offset and is already absorbed by `bias_offset_v`. Recorded with four sources (including the RC-circuit model form) and a falsifiable trigger: implement when a protocol simulates tens of minutes or longer.
- **Overfitting checked against the codebase.** Five modules construct `PhaseModulator` and the Gobby chain is not one of them. `bias_offset_v` is therefore exercised by `bb84_duplinskiy` in tests — accepting it, staying bit-identical at zero, raising QBER at half-V_pi.
- 295 tests pass (276 + 19). Transmission still 0.7692307692307693; default `'balanced'` path bit-identical at 0/65/122 km with both new parameters present and zero.

**Known omission, with its number:** the gap between our device visibility of 1.0 and the paper's stated >0.999 bound is real residual imperfection worth at most 0.050% QBER, ≤1.5% of the floor. Not modelled, the same treatment linewidth received.

---

## 2026-08-08 — PERF-2 extended to stochastic phases; linewidth and modulation error

### Session: GOBBY-6 — the floor reproduced from stated parameters, nothing fitted

| Change | Files | Rationale |
|---|---|---|
| **PERF-2 extended, not refactored** | `src/protocols/bb84_time_bin.py` | PERF-2's 8 precomputed outcomes are exact only while the response is deterministic; a per-pulse random phase makes it continuous. The same linearity argument extends: `P(delta) = g0 + 2*Re(S*exp(i*delta))`, with `(g0, S)` fixed by three evaluations at `delta = 0, pi/2, pi`. Three field propagations **per point** instead of per pulse; any phase thereafter is two multiplies. Documented in the module docstring — it is why a 1e8-pulse sweep is feasible and is easy to lose in a later refactor. |
| **Coefficients extracted from the chain, not re-derived** | same | No second expression of the physics to drift out of step, and correct for both topologies without special-casing. |
| **Consistency assertion added** | same | See below — the negative control exposed that the closed form had stopped deriving the phase-sign convention from the physics. |
| **`linewidth` / `path_mismatch`** | `src/protocols/bb84_time_bin.py`, `analysis/val_gobby/validate_gobby.py` | The chain had no laser model at all. Every real laser has a linewidth, and a 2026 replicator must be able to express their apparatus. |
| **`phase_error_rad` / `phase_noise_rad`** | `src/channel/phase_modulator.py` | Modulation error as a static offset or per-pulse jitter. |

**Key results:**

- **The floor is reproduced from stated parameters with nothing fitted.** `e_mod + P_e` measures **3.407 ± 0.226 %** at 0 km against Gobby's stated **3.3 %** (a second independent run gave 3.579 ± 0.267 %; both within ~1σ).
- **Contribution budget — what earns its place.** `e_mod` is **85 % of the floor** and load-bearing: without it the chain has no floor at all (0.08 % with afterpulsing correctly off). `P_e` is load-bearing at range, not at the floor (0.03 % at 0 km → 5.96 % at 122 km). **Linewidth is negligible at realistic trim** — 0.0047 % predicted at 3 MHz / 10 ps, 0.14 % of the floor — so it is **kept and documented with its number**, rather than silently included as though it mattered or omitted so a replicator with poor trim cannot express their apparatus.
- **Linewidth couples through the residual mismatch, not the AMZI delay.** The S-L and L-S routes traverse the same total path, so the frequency-noise term cancels — which is why the delay line and stretcher exist, and how >99 % visibility works with an 80 ps pulsed source. Measured against the closed form at four (dnu, d_tau) points, all within ~1σ. Producing the full 3.3 % floor this way would need `d_tau` of 2.17–21.7 ns, larger than the 5.8 ns delay itself. **That is why the paper omits it** — in a path-matched scheme it is irrelevant to the QBER. The omission is correct; no claim about the source is made.
- **Both modulation-error models give 3.299 %** — static `d = 20.93°`, jitter `s = 21.17°`, a 1.1 % difference in angle. The choice describes the hardware rather than matching a number; jitter is the default since drive noise is random per pulse.
- **A negative control caught a hole in itself.** The gate first read "failure" at relative difference 1.0 on the four matched-basis entries — a comparison artifact, since those are the *extinguished* ports where relative error is undefined (baseline ~8.7e-44 by direct interference, closed form ~1e-28 from cancellation, i.e. float64 epsilon × full scale). On absolute error against `P_max` the form is bit-identical at **1.785e-16**, with a 1 mrad offset showing at 5e-4. But **the Bob-sign control gave the same value as the pass** — vacuous, because the coefficients are extracted at `phi_B = 0` and so never observe the sign, leaving the per-pulse formula to *assume* `delta = phi_A - phi_B`. The GOBBY-4 §21.2 phase-arm coupling would no longer have been caught. Fixed by asserting the closed form against the chain at `phi_B != 0`; verified that flipping the sign now trips it.
- 276 tests pass. Transmission still exactly 0.769231; bit-to-port mapping unchanged in both bases; default `'balanced'` path bit-identical.

**Recorded caution:** 21° is a derived consequence of taking Gobby's stated 3.3 % at face value within this model — not a measurement of their modulator, and not checked against era specifications.

**Not usable, not quoted:** the 122 km column of the first budget run produced `e_mod + P_e` *below* `P_e` alone — non-physical, marking the run as under-powered rather than anything being wrong. The long end is characterised analytically; a properly powered long-range budget belongs with the nine-point sweep.

**Still open:** the nine-point sweep, deferred until `e_mod` existed. With it in place and nothing fitted, that run becomes the test of the raw claim.

---

## 2026-08-08 — Retarder global phases removed; phase-arm change confirmed MC-vs-analytic

### Session: GOBBY-5 — the conjecture pays off

**A conjecture that something in the chain was not real-valued is what
produced this fix.** The X basis was the suspected site and came back clean,
but the hypothesis directed an audit that found the actual defect — nobody
inspects retarders no live code calls without a reason to.

| Change | Files | Rationale |
|---|---|---|
| **Retarder global prefactors deleted** | `src/channel/optics.py` | `halfwave` carried `exp(-i*pi/2) = -1j` and `quarterwave` `exp(-i*pi/4)`, multiplying the entire Jones matrix. Unlike the coupler `i` (a choice between valid conventions), a global prefactor has **no physical content** — invisible in `\|E\|^2` for one path, a real relative phase as soon as two interfering paths pass different retarder counts. The standard Jones form omits it, and the module's matrix body was already that form, so this was a deletion rather than a reformulation. |
| **Retardance preserved** | `src/channel/optics.py` | The relative `i` between `quarterwave`'s fast and slow axes stays — that is the retardance and the whole point of the component. |
| Convention note | `src/channel/optics.py` | Added beside the existing coupler note, distinguishing a global prefactor (removed) from a relative axis phase (kept). |
| Phase-pinning test corrected | `tests/test_optics.py` | `test_quarterwave_zero_angle_leaves_h` asserted `exp(-1j*pi/4)` under a docstring reading "up to a global phase" — it asserted the very thing it called incidental. Now `1.0`. |

**Key results:**

- **`halfwave` \|Im\|/\|Re\|: 1.6e+16 -> 0.** A real input had been coming out *purely* imaginary. `quarterwave` H-aligned: 1.0 -> 0. Retardance at 45° preserved (`E_y = -i E_x`).
- **No power-based test moved**, which is the expected signature: a unit-modulus prefactor cannot change `|E|^2`. Every other retarder test was already phase-blind — itself evidence the prefactor carried nothing.
- **Phase-arm change confirmed end to end.** §21.2 verified it at the gate-table level; this pass tests the whole chain against the closed form. Link budget only (`afterpulse_prob=0`): MC/analytic **1.278 / 0.893 / 1.114** at 65/100/122 km, weighted mean **1.061 ± 0.131 — 0.47σ from unity**, every point within 1σ.
- **The `MU_EFF` prediction held.** On identical MC numbers, the old fitted 0.0793 gave a weighted mean of 1.091 and the derived value gives 1.061 — the ratio moved **2.8 % toward 1.0**, matching the 3.12 % mismatch it removed.
- **The 3.12 % provably cannot propagate.** The elasticity is analytic: `(dQBER/QBER)/(dS/S) = -V·(e_counts/QBER) <= V <= 1`. Measured 0.009 / 0.127 / 0.405 / **0.567** at 4.4/65/101/122 km — predicted ΔQBER matching a finite difference to the third decimal. **Worst case 0.567, strictly sub-unity**: a 3.12 % move in `mu_eff` yields at most 1.77 % in QBER. The mapping is a contraction, not an amplification. `ALPHA_DB`, `ETA_BOB`, `P_E`, `E_MOD` are independent literals and do not move at all; the MC never reads `MU_EFF`, which is what keeps the 0.769/0.793 comparison a cross-check rather than a tautology.
- 276 tests pass (was 264). Transmission still exactly 0.769231; bit-to-port mapping unchanged in both bases; default `'balanced'` path bit-identical.

**Still deferred:** `coupler_combine`'s `1j` (GOBBY-4 §21.3) — a legitimate
coupler convention merely inconsistent with the real form this module
adopted, which is a different category from a prefactor carrying no physical
content. Modulation error remains unimplemented (GOBBY-2 §19.9 step 3), so
only the link-budget comparison is currently meaningful.

---

## 2026-08-08 — Analytic model finalised (no fits); encoder phase arm corrected

### Session: GOBBY-4 — remove the last fitted parameter; put the phase on the arm that carries it

**The analytic model contained one fitted parameter and its docstring denied
it.** `src/analytic/gobby_model.py` set `MU_EFF = 0.0793` by *inverting
Gobby's measured fringe visibilities* — fitting to the very data the model
exists to predict — beneath a header reading "Every parameter is from the
source paper. Nothing is fitted." That single number was the entire
"analytic vs Gobby" discrepancy; the formulations and stated parameters were
correct throughout.

It was standing in for the **interferometer transmission**, which is
derivable from stated values:

```
T_int  = 2/(1 + r) = 0.769231      r = 1.6, stated
mu_eff = mu * T_int = 0.076923     mu = 0.1, stated
```

| Change | Files | Rationale |
|---|---|---|
| **`MU_EFF` derived, not fitted** | `src/analytic/gobby_model.py` | `MU`, `SPLIT_RATIO`, `T_INT` added; `MU_EFF = MU * T_INT`. Every factor stated or derived in closed form from a stated value. Per-constant citations. |
| **No-fitting rule written into the module** | `src/analytic/gobby_model.py` | Stated as a rule, not a description: no parameter — physical models or analytic comparisons — may be fitted to data the model exists to reproduce, unless the source states it as fitted. **Every parameter must carry a citable source.** |
| **`phase_arm` on `AsymmetricMZI`** | `src/channel/interferometer.py` | Which arm holds the modulator fixes the *sign* of the relative phase. Default `'long'` preserves all existing behaviour. |
| **Encoder phase moved to the encoded arm** | `src/protocols/bb84_time_bin.py` | Gobby encodes on Alice's **short** arm; the code was applying `phi_A` to the long arm, which carries the reference. Two coupled changes: `phase_arm='short'` **and** Bob's sign `exp(-i*phi_B)` -> `exp(+i*phi_B)`. |

**Key results:**

- **Cost of removing the fit: 0.026 pp.** Residuals against Gobby's four measured points are +0.03 / +0.49 / −0.25 / +0.36, **mean 0.282 pp** against the fitted value's 0.256 pp. The model now predicts rather than reproduces.
- **Independent cross-check retained** (as a check, never an input): inverting Gobby's stated visibilities gives `mu_eff/mu = 0.793` against the geometry's 0.769 — **3 % agreement between two unrelated routes**.
- **The traceability guard is exact equality.** `test_mu_eff_is_derived_not_fitted` asserts `MU_EFF == MU * T_INT`, so writing the literal 0.0793 fails even though the two agree to 3 %. A rule that is only a comment is not a rule.
- **Phase-arm fix verified by gate-table equivalence, captured before any edit.** Production vs baseline across 3 distances × 8 entries × 2 ports: **max relative difference 2.35e-16, nothing changed**. Intensities are untouched because `cos` is even — which is exactly why the bug survived every existing test.
- **With a negative control, because a check that cannot fail proves nothing.** Applying the arm move *without* the sign flip changes the table by **1.00** and does so on **Y0Y and Y1Y only** — the Y basis inverts, X is untouched. That confirms the equivalence check can actually detect a sign error.
- 264 tests pass (was 259). Default `'balanced'` path bit-identical at 0/65/122 km; transmission still exactly 0.769231.

**Convention audit — recorded, deliberately not fixed.** Measured
`max|Im|/max|Re|` through the polmux chain in the X basis: **0.000e+00 at
every stage — real end to end**. The imaginary fields are not from the
time-bin chain. Per component: `coupler_split`/`pbs`/`pbc`/`voa` clean;
`coupler_combine` **1.0** (still the `1j` convention); `circular_analyser`
**1.0** (by design, and live in four protocols); `quarterwave` **1.0**;
`halfwave` **1.6e+16** — a real input comes out *purely* imaginary. Written
up as GOBBY-4 §21.3 with scope for a later pass. `|E|^2` is unaffected so no
published number is wrong, but the simulator is meant to be physically real.

**Numbers superseded:** §19.3's decomposition table and §19.2's derivation
were computed at the fitted 0.0793 and are correct only for that value.
§19.3's sensitivity row "0.0793 / from V / 0.26 pp" becomes
"0.0769 / derived / 0.28 pp".

---


## 2026-08-08 — Generic unbalanced-MZI components; the μ/2 loss was ours

### Session: GOBBY-3 — polarisation-multiplexed topology from generic parts

The Gobby chain gated only `μ/2`. §19.7 resolved that by feeding
`2·μ_eff = 0.1586` — a scalar with no physical referent, putting 1.59× the
paper's stated `μ` into a replication whose claim is source traceability.
The real cause was modelling the wrong interferometer: a balanced 50:50
pair produces four paths and dumps half the light into satellite bins the
gate discards, whereas Gobby's polarising beam combiner/splitter routes
deterministically, so `S-S` and `L-L` never form.

Fixed by adding the missing **generic** capability rather than a
Gobby-shaped class. Nothing in `src/` knows the number 1.6.

| Change | Files | Rationale |
|---|---|---|
| **`coupler_split` fixed** | `src/channel/optics.py` | It applied a bogus Hadamard to the fields and let `ratio` move only the returned *powers* — the two were mutually inconsistent. **The old field behaviour was wrong, not merely different**; any caller relying on it was relying on a bug (only tests did). Now amplitude-splits as `√ratio` / `√(1−ratio)` per **Zeilinger, Am. J. Phys. 49(9), 882–883, 1981**. |
| **`pbc` added** | `src/channel/optics.py` | Polarising beam combiner, inverse of `pbs`, per **Collett, *Field Guide to Polarization*, SPIE 2005**. A non-polarising combiner interferes its inputs and loses half the light; a PBC keeps both because they occupy orthogonal states. This is what makes deterministic routing expressible. |
| **Coupler phase convention documented** | `src/channel/optics.py` | Zeilinger: `[[t, ir],[ir, t]]` and `[[t, r],[r, −t]]` are both unitary, equivalent *only if used consistently*. Mixing them is not cosmetic — the imaginary form makes interference go as `sin Δφ`, and BB84 encodes in `{0, π}` where `sin` vanishes. The module now commits to the real form. |
| **`AsymmetricMZI` generalised** | `src/channel/interferometer.py` | `split_ratio=0.5` (default preserves behaviour exactly); `recombine=False` taps the arms; the decoder accepts a pre-split `(E_a, E_b)` pair. Arm tapping and injection are ordinary requirements for an unbalanced MZI. |
| Topology composed, not hard-coded | `src/protocols/bb84_time_bin.py` | `interferometer='polarisation_multiplexed'` wires split → delay → modulate → `pbc` → fibre → `pbs` → delay → `voa` → recombine. Balancing uses the **existing** `voa` at `10·log10(r)` dB. `split_ratio` is a caller argument. |
| Bespoke class deleted | `src/channel/interferometer.py` | An earlier `PolarizationMultiplexedAMZI` hard-coded the topology and re-implemented coupler algebra inline. Removed. |

**Key results:**

- **The geometry predicts the signal.** Equalising the arms discards the excess reference power: transmission `2·κ_A` = **0.769231 measured, exact**. Inverting Gobby's *stated* visibilities via `V = S/(S + 2·P_e)` independently gives `μ_eff/μ = 0.793`. **3 % agreement, zero free parameters** — and it explains the ~21 % shortfall §19.7 left unattributed.
- **MC/analytic excess gone.** Link budget only: 1.36 / **0.976** / 1.31 at 65/100/122 km, consistent with 1.0. §19.10(2) recorded **1.9–2.5×**. The best-powered point matches the 0.769/0.793 = 0.970 ratio the two independent derivations predict. §19.10(2)'s premise is withdrawn — the excess was structural, not detector-side.
- **Default path bit-identical** at 0/65/122 km; balanced chain still gates exactly 0.500.
- **Two bugs found that the existing tests could not catch**: the `sin Δφ` convention error, and Bob's phase on the encoded arm giving `φ_A + φ_B`, which inverts *only* the Y basis. Both produced ~48 % QBER end to end. Twenty passing tests missed them because each swept phase continuously or checked energy — none asserted bit→port. `TestBB84Mapping` now does.
- 259/259 tests pass (was 242).

**Still open:** modulation error is not implemented (GOBBY-2 §19.9 step 3), so the chain's floor comes from afterpulsing while the analytic carries `e_mod = 3.3 %` — different mechanisms, similar magnitude, exactly the pattern §19.5 warns about. Only the link-budget comparison is currently meaningful. Arm equalisation is an inference, not a paper statement. The nine-point sweep has not been re-run at the corrected topology.

---

## 2026-08-07 — GOBBY-2 step ②: the afterpulse=0 sweep — the floor collapses, and the MC/analytic excess is a 3 dB encoder split, not the detector

### Session: run §19.9(2) at full power, root-cause the ~2× MC/analytic excess

Step ②'s prediction held exactly where it was sharp and failed where the
step was designed to look: **the 3.3 % floor collapses** — 0 km goes
**2.55 % → 0.08 %** with `afterpulse_prob = 0` — so afterpulsing really
was standing in for the missing e_mod. But the MC/analytic excess does
**not** collapse: it holds at **~2×** across the sweep. It was never a
detector-state-machine effect — it is the **50:50 encoder split**.

**Result (9 points, `--target-sifted 3000`, seed 42, ≈1.64×10⁹ pulses,
≈1 h 20 m logged run):**

| z (km) | MC | analytic | MC/analytic | Gobby |
|---|---|---|---|---|
| 0 | 0.08 ± 0.05 | 0.019 | 4.4× | 3.3 (i) |
| 4 | 0.03 ± 0.03 | 0.023 | 1.2× | 3.3 (m@4.4) |
| 10 | 0.08 ± 0.05 | 0.030 | 2.8× | 3.3 (i) |
| 20 | 0.11 ± 0.06 | 0.047 | 2.3× | 3.3 (i) |
| 40 | 0.30 ± 0.09 | 0.119 | 2.5× | 3.3 (i) |
| 65 | 0.87 ± 0.14 | 0.374 | 2.3× | 3.3 (m) |
| 80 | 1.77 ± 0.22 | 0.741 | 2.4× | 4.4 (i) |
| 100 | 3.51 ± 0.27 | 1.82 | 1.9× | 6.0 (m@101) |
| 122 | **8.90 ± 0.40** | 4.71 | 1.9× | 8.9 (m) |

Rows below 40 km are noise: with the floor gone everything reads
0.03–0.11 %, and the ratio column there is meaningless.

### Root cause — the chain's gate captures μ/2, not μ

A balanced double AMZI (50:50 encoder + 50:50 decoder) sends **half** the
launched energy to the central interference peak and a quarter to each
satellite bin, which the gate excludes. The chain therefore delivers
`μ/2` at the gate. Verified twice:

- direct check (2M pulses, dark = 0, afterpulse = 0, V = 1):
  `P(click) = 0.002241` vs `μ·η_Bob = 0.0045` → **ratio 2.008**;
- sifted fraction at 0 km: `1.10e-3` = `(μ·η_Bob)/2` exactly.

The closed form with the halved signal then reproduces the Monte Carlo:
**3.51 % @100 km** (the measured 3.5054, to the digit), **8.61 % @122 km**
(measured 8.90; the 0.29 pp tail is dead-time/tie-break), and — with
GOBBY-1's `a = 0.05` — **10.59 % @122 km** vs GOBBY-1's measured 10.69.

**This corrects the GOBBY-1 attribution.** §18.6's "1.56× at 122 km —
excess in the detector state machine" was the analytic side comparing a
full-`μ` closed form against a chain that gates `μ/2`. There is no
state-machine excess to instrument.

### The 122 km point is a coincidence — do not report it as validation

8.90 % vs the paper's 8.9 % looks like a hit but is two compensating
defects: the missing e_mod (−3.3 pp, flat) cancels the halved-signal
error share (+3.9 pp, distance-dependent) only near 122 km. At 65 km the
same two effects read **0.87 % against Gobby's 3.3 %**.

### Decision recorded for steps ③–④ (§19.7, §19.9)

Feed the chain `mu = 2·μ_eff = 0.1586` so the gate receives
`μ_eff = 0.0793` — QBER-identical (`S/(S+2·P_e)` unchanged) to modelling
Gobby's apparatus directly — chosen over implementing a 1.6:1 encoder or
documenting the 2× offset. e_mod goes on the `PhaseModulator` as
`phase_error_rad = 0.3653 rad` (static offset, §19.6). Code changes for
steps ③–④ were deliberately **not** made this session; this sweep record
and the root cause are the deliverable.

---

## 2026-08-07 — GOBBY-1: the link budget, corrected; Gobby's slope reproduced with no fitted parameters

### Session: implement §18.5, run the nine-point sweep, fix a derived-column defect

The headline: **the simulated curve now has Gobby's slope** — fitted
**+4.99 pp/100 km** against the paper's own **+4.39**, where the previous
sweep fitted **+0.011** (flat). The flatness was never missing physics; it
was a mis-specified link budget. Correcting it removed the flatness
**without adding a mechanism and without fitting anything**, and
supersedes the two-parameter fit recorded on 2026-08-06.

**Result (9 points, `--target-sifted 3000`, seed 42, 4 h 33 min,
1.56×10⁹ pulses):** mean |residual| **1.22 pp** over Gobby's four
genuinely measured points, signed mean −0.33 pp.

| z (km) | Simulated | Gobby | Δ |
|---|---|---|---|
| 4 | 1.96 ± 0.23 | 3.3 (measured @4.4) | −1.34 |
| 65 | 2.86 ± 0.25 | 3.3 (measured) | −0.44 |
| 100 | 4.67 ± 0.35 | 6.0 (measured @101) | −1.33 |
| 122 | **10.69 ± 0.43** | 8.9 (measured) | **+1.79** |

Compare: the previous session scored 0.36 pp with **two** fitted
parameters; this scores 1.22 pp with **zero**. That is the intended
trade — the old fit landed close for the wrong reason, its 1,788 Hz
"dark count rate" being a lumped term absorbing the missing 5 dB of Bob
loss and the wrong α.

### Re-parameterisation — `analysis/val_gobby/validate_gobby.py`

| Parameter | From | To |
|---|---|---|
| `ALPHA_dB` | 0.182 | **0.2** (paper's value) |
| `ETA` → `ETA_BOB` | 0.10 | **0.045** — η_Bob, including Bob's 5 dB apparatus loss, which the bare detector QE silently discarded |
| `GATE_WIDTH` | 1 ns | **3.5 ns** |
| `REP_RATE` / `PULSE_WIDTH` | 2.5 MHz / 100 ps | **2 MHz / 80 ps** |
| `DCR` → `P_E` | 15 Hz | **8.5e-7/clock** — Gobby's *measured* total error probability (dark 3.2e-7 + 1.3 µm clock-laser stray light 5.3e-7, deliberately lumped), 242.9 Hz equivalent |
| `VISIBILITY` | 0.934 injected | **1.0** — V is an OUTPUT, `V = S/(S+2·P_e)` |
| `--dcr` | — | **`--p-e`**, plus new `--afterpulse` |
| `PILOT_BITS` / `CEILING` | 200k / 500M | **2M / 1e9** |

Net effect of the old parameterisation: signal **3.68× too high**, error
rate **57× too low**, error/signal ratio — the only thing QBER measures —
off by **209×**.

### Bugs found

| Bug | Location | Consequence |
|---|---|---|
| **`alpha_dB` never passed to the Monte Carlo** | `simulate_qber()` | The MC used `simulate_bb84_time_bin`'s own 0.182 default while the analytic curve used the module constant. They agreed *by coincidence*. Changing `ALPHA_dB` alone would have moved the curve and left the data untouched, and the figure would have read as a modelling disagreement. |
| **Nine rows of "measured" over four measurements** | table writer | `gobby_measured_qber()` is `np.interp` over four published points, but the column was headed "Gobby et al. measured" at all nine distances. Below 4.4 km `np.interp` clamps rather than interpolates, so a flat 3.3 % across four rows was **one published number repeated four times**, presented as data. Same species as BLOCK-1, one step milder. Rows now carry `(m)` / `(m@x)` / `(i)`; measured rows print the published value at the published distance (`6.0 (m@101)`, not the interpolated 5.9). |
| Stale reference comment | `src/channel/interferometer.py` | Ref [4] asserted "Gobby's 3.3 % floor is interferometer visibility, so V = 0.934" — the claim §18.4 refutes, and the origin of the injected-visibility error. Corrected in place. |

### Verified before the sweep

- `model_qber(V=1, a=0)` reproduces the closed form `P_e/(S+2P_e)` to 1e-9.
- `predicted_visibility()` reproduces the paper's stated visibilities:
  **0.9925 at 65 km** (paper: >0.99) and **0.9058 at 122 km** (paper: 0.884).
- **The gate widening is energy-neutral for signal** — 0.050000 photons at
  1 ns vs 0.050057 at 3.5 ns — while background scales exactly 3.5× and
  lands on **8.500e-7**, Gobby's measured P_e to the digit.
- 203/203 tests pass.

### Open, and deliberately not closed by fitting

- **122 km overshoots by 1.79 pp (4.1σ)**, having previously undershot by
  3.60. The MC runs **1.56× above the closed form** there — the largest
  divergence in the sweep, in the regime the closed form disclaims — so
  the excess is in the detector state machine, not the link budget. The
  encoder split (our 0.05 photons vs Gobby's 0.04) predicts the *opposite*
  sign, so the underlying excess is larger than 1.79 pp.
- **The 3.3 % floor is undershot by 0.4–1.3 pp.** Afterpulsing at the
  ID230 datasheet value supplies ~2.3 pp of the paper's acknowledged third
  error source; the remainder is unmodelled and reported as unmodelled.
- Next step is a decomposition run (`--afterpulse 0`, `--p-e 0`), **not**
  tuning `afterpulse_prob` — which would work, and would turn the whole
  sweep back into a fit.

---

## 2026-08-06 — Gobby replication measured and diagnosed; statistical-power tooling

### Session: OPEN-2/OPEN-3 full-power runs, V/DCR refit, analytic-model fix, two memory bugs

The headline: the 122 km Gobby point was **5.30 %**, not the 8.33 % a
12-sifted-bit run had reported, and the simulated curve was **flat**.
After diagnosis and refitting, the chain reproduces Gobby across 122 km to
a mean residual of **0.36 pp**.

| Change | Files | Rationale |
|---|---|---|
| `--target-sifted N` budgeting | `analysis/val_gobby/validate_gobby.py`, `analysis/val_system_scenarios.py` | The sifted fraction spans 164× across the Gobby sweep (2.3e-3 at 0 km to 1.4e-5 at 122 km), so a flat pulse budget starves the only end anyone scrutinises. Flat 10M leaves 122 km at σ = 2.34 pp against a 3.78 pp effect — 1.6σ, unable to settle anything. Pilots the sifted rate, scales, retries on undershoot. Measured 288k pulses at 0 km vs 4.45M at 65 km for matched 0.85 pp error bars. |
| Write guard | both scripts | `check_statistical_power()` runs **before any file is written** and raises with an actionable message. Caught the documented pathologies exactly: Gobby at 20k pulses (*"122 km: 0 sifted"*), scenarios at 60k (*"+ CD: 4 sifted; + PMD: 4 sifted"* — four rows identical because they held the same single error). `--allow-underpowered` for smoke runs. |
| **Two O(N) memory bugs** | `analysis/val_system.py`, `src/protocols/bb84_time_bin.py` | Both accumulated five per-pulse Python lists then sifted afterwards — multiple GB at the ~1e8 pulses `--target-sifted` needs. **This killed the first scenario sweep** (observed at 2.1 GB, climbing). Now sifted inline with counters, O(1). No RNG draw touched; **bit-identical**, verified at four configurations each. |
| Targeting knife-edge | both scripts | `--target-sifted 3000` defaulted `--min-sifted` to the same 3000, and a one-shot pilot extrapolation undershoots ~half the time. An 8-row sweep completed *and* passed its code-path check, then was rejected at the final write for landing at 2422–2668 — an hour discarded for a 20 % miss. Retry with `1.15 + 3/√n` headroom; `MIN_SIFTED` now below target. |
| **`model_qber()` rewritten** | `analysis/val_gobby/validate_gobby.py` | Disagreed with the Monte Carlo by ~3× on the dark term (ratios 2.5–3.6), and both were plotted in the same figure as if they agreed. Three defects: counted one detector's dark rate when both dark-count; halved the dark error term that should be whole; added misalignment as an unweighted constant rather than a share of clicks. Afterpulsing absent entirely. Replaced with a single weighted ratio; both limits now exact. Mean \|residual\| vs MC: **2.86 → 0.59 pp**. Docstring carries a "do not fit parameters against this function" warning. |
| CD code-path check | `analysis/val_system_scenarios.py`, `analysis/val_system.py` | `simulate_point` gained a `delay` parameter (`delay=None` bit-identical). At the Gobby 5.8 ns delay CD cannot move QBER at all — CD and the AMZI are both LTI and commute, so CD acts identically on both ports and cannot move the port ratio the bit comes from; bin crosstalk would need z ≈ 5,674 km. Criterion is the **sifted rate**, not the QBER: measured Δ(QBER) ≤ 0.7σ but Δ(sifted) 6–12σ at 191/400/800/1500 km. Reports LIVE at 11.6σ. |
| `--dcr` exposed | `analysis/val_gobby/validate_gobby.py` | Detector parameters are a *fitted input* to this replication, not a given. |

**Key results:**

- **§16 hypothesis refuted.** 122 km measured 5.30 ± 0.49 % against Gobby's 8.9 % (7.3σ). The old 8.33 % was 1 error in 12 bits and landed near the published value by chance. Simulated slope +0.011 pp/100 km — flat. Visibility is distance-independent, so it could never have created a distance trend.
- **Two defects, opposite signs, partly cancelling.** (a) V = 0.934 was fitted with noise *off* (3.357 %) then applied with `afterpulse_prob = 0.05`, double-counting the floor — afterpulsing measured in isolation at **2.262 ± 0.176 %**. (b) DCR = 15 Hz is the ID230 spec, a 2020-era part on a 2004 experiment; `P_dark/p_signal = 2.5e-4` even at 122 km, so nothing could produce Gobby's rise.
- **Refit reproduces Gobby.** V = 0.9792, DCR = 1,788 Hz (bisected against the **MC**, not the analytic model): mean \|residual\| **0.36 pp**, 8 of 9 points within 1.4σ, 122 km at **9.02 ± 0.57 %** vs 8.90 %. The three 122 km values across the session — 5.30, 17.49 (analytic-fitted DCR), 9.02 — bracket the diagnosis. **See the retraction below: this fits, but for the wrong reason.**
- **Caveats recorded:** two free parameters; 100 km remains a genuine 3.1σ miss at a *measured* Gobby distance.
- **Scenario nulls now real.** CD, PMD and birefringence rows came out **bit-identical** to attenuation-only (2514 sifted, 74 errors each) versus the original 86 sifted with one shared error, while visibility (6.11 %) and dark counts (13.01 %) move by 6σ and 14σ through the same machinery.
- 203/203 tests pass.

**Scenario sweep (BLOCK-2 / OPEN-3) completed and emitted.** All 8 rows cleared the target: CD, PMD and birefringence came out **bit-identical to attenuation-only at 3758 sifted bits** (2.93 ± 0.27 %, 110 errors each), while visibility (6.05 %) and dark counts (13.24 %) move by 8σ and 19σ through the same code path. The CD code-path check reports **LIVE at 11.6σ** on the sifted rate. `val_system_scenarios--seed42.{csv,tex}` written.

### Retraction, same session — GOBBY-1

The paper's own parameters were then located, and they **supersede the refit above**. Recorded in full as **GOBBY-1** in `opto-sim-issues-and-fixes.md` (specified, deliberately **not implemented**).

Gobby states **α = 0.2 dB/km**, **η_Bob = 0.045** *(including Bob's 5 dB apparatus loss)*, **gate width 3.5 ns**, **2 MHz**, **80 ps**, and a measured **error probability per clock cycle P_e = 8.5×10⁻⁷** — of which only **3.2×10⁻⁷ is detector dark count**, the remaining **5.3×10⁻⁷ being stray light from the 1.3 µm clock laser**, a source this chain does not model at all.

Against those: our **signal is 3.68× too high** (wrong α, and Bob's 5 dB silently discarded) and our **error rate 57× too low**, making the error/signal ratio at 122 km **209× too small**. That — not missing physics — is why the sweep came out flat.

- **The "119× the ID230 spec" characterisation is withdrawn.** Correcting the signal by its 3.68× error puts the required error rate at 4.85×10⁻⁷ against Gobby's measured 8.5×10⁻⁷ — **within 1.75×**. The fitted 1,788 Hz was never a dark count rate; it was a lumped error term absorbing the missing Bob loss and wrong α. Calling it a DCR was a category error.
- **Visibility is an output, not an input.** The paper gives `V = S/(S + 2·P_e)` with `S = mu·10^(−αL/10)·η_Bob`. Verified against its own stated values: **99.25 % at 65 km** (paper > 99 %) and **90.58 % at 122 km** (paper 88.4 %). Injecting `visibility=0.934` therefore double-counts the error physics — the same mistake as the V/afterpulsing double-count, one level up. OPEN-1 and the entire §16 visibility hypothesis rest on it.
- Residual after (1−V)/2 is **2.9–4.2 pp and roughly distance-independent**, matching the paper's acknowledged **third error source** and consistent with the 2.26 % afterpulse floor measured here.

**Still outstanding:** implement GOBBY-1 §18.5 (re-parameterise; drive error counts from the measured P_e instead of fitting; stop injecting visibility), after which the 122 km point becomes a *prediction against measured inputs* rather than a two-parameter fit. `val_gobby_table.tex` must **not** be regenerated at the superseded (V = 0.9792, DCR = 1,788 Hz) values.

---

## 2026-08-06 — Eighth pass (OPEN-3 groundwork): gate-width knob, runtime measurement, budget decision

### Session: OPEN-3 measurement + `gate_width` parameter

OPEN-3 is the scenario-table regeneration; before committing to a run
budget the per-pulse cost was measured on the actual chain.

| Change | Files | Rationale |
|---|---|---|
| `simulate_point()` gains `gate_width` (default 1 ns) controlling both the SPAD gate and the integration window | `analysis/val_system.py` | Enables the OPEN-3 CD positive control: shrink the gate below the CD-broadened pulse width (30 → 78 ps at 100 km, L_D ≈ 41.5 km) so spill-over measurably moves the QBER. Additive; the default path is bit-identical (regression point re-verified: 28 sifted / 1 err / 3.571 % @ 500k, 100 km att-only). |

**Runtime measurement (the point of the exercise):** `simulate_point` runs
at **35.8 µs/pulse** @ 100 km attenuation-on — the ~10 µs/pulse PERF-2
figure applies to the Gobby chain only, not to this SPAD Monte Carlo.
Budget decision (user): **lean — ≥3×10³ sifted per row, all 8 rows**;
σ(QBER) ≈ 0.2 %, >8σ separation between attenuation rows (1.16 %) and the
2.81 % control. ≈7×10⁷ pulses/row ≈ 42 min, ~5 h total (10⁴ sifted would
be ~2.3 h/row ≈ 16 h).

**CD positive-control prototype (gate 150 ps, 100 km, 2M pulses, seed 42):**
DCR 10 kHz → CD 2.74 % vs no-CD 2.27 % (2 errors each — inconclusive);
DCR 50 kHz pair aborted mid-run, pending. Generator update + the ~5 h
regeneration are the remaining OPEN-3 work; OPEN-2 (10M Gobby with
V = 0.934) follows.

---

## 2026-08-06 — Seventh pass: PERF-2 + OPEN-1 — Gobby chain unblocked

### Session: PERF-2 / OPEN-1

The statistics bottleneck was the binding constraint on every remaining
result (the 10M Gobby re-run, the scenario table). PERF-2 removed it by
memoising the field chain; OPEN-1 wired the visibility the replication
was always missing.

| Change | Files | Rationale |
|---|---|---|
| 8-outcome `(P_c, P_d)` gate-power table built once per point; per-pulse loop reduced to a dict lookup + SPAD Monte Carlo. RNG draw order untouched. | `src/protocols/bb84_time_bin.py` | The chain is scalar attenuation + deterministic AMZIs, so 2M pulses performed 4M field propagations to obtain 8 distinct answers. Measured: 75 s → 1.9 s per 200k point (~40×), ~10 µs/pulse → a 10M point ≈ 1.7 min. |
| `VISIBILITY = 0.934` at module scope (cited to Gobby's 3.3 % floor via `e_opt = (1−V)/2`), passed into `simulate_bb84_time_bin`; `--visibility` CLI flag keeps the V = 1.0 control runnable. | `analysis/val_gobby/validate_gobby.py` | OPEN-1: the validator silently ran a perfect decoder (V = 1.0), so its floor was pure afterpulsing, not the misalignment the paper compares against. |

**Verification (the refactor acceptance test).** Old vs new implementation,
200k pulses, distances {0, 4.4, 65, 101, 122} km, seed 42 — `qber`,
`n_sifted`, `n_errors` bit-identical at every distance; `modulate` call
counter = 12 per point regardless of pulse count. OPEN-1 physics check:
V = 0.934, zero noise, 0 km → **3.357 %** (Gobby's 3.3 % floor); V = 0.934,
full noise, 122 km → **8.33 %** vs Gobby's 8.9 %.

**Tests:** 203/203 pass; `run_all.py --seed 42` 8/8 in 95.2 s (Gobby
24.5 s, was 158 s).

---

## 2026-08-06 — Sixth pass: BLOCK-2 — impairment-table generator

### Session: BLOCK-2

The last code-side item: Table 11 (`paperwork/tables/val_system_table.tex`)
was hand-written — no generating script, no CSV, no seed, no commit.
`paperwork/` has since been deleted, but the reproducibility contract is
now met by a named generator on the ARCH-1 time-bin chain.

| Change | Files | Rationale |
|---|---|---|
| `val_system_scenarios.py` — 8 explicit impairment configs, each run at 100 km with recorded seed; prints the exact config dict per row; writes CSV + LaTeX table with script/seed/pulses/commit hash in the caption | `analysis/val_system_scenarios.py` (NEW), `val_system/val_system_scenarios--seed42.{csv,tex}` | Every table row is self-contained and regenerable at any budget (`--bits`, `--distance`). Replaces the hand-written table. |
| `simulate_point()` gains independent `birefringence`/`attenuation`/`cd`/`pmd` toggles (`dispersion` remains the legacy alias) | `analysis/val_system.py` | The scenarios need per-impairment control; defaults unchanged, so the ARCH-1 panels are bit-identical. |

**Results (2M pulses/row, seed 42, 100 km):** no impairments 2.81 % ± 0.24
(4761 sifted); attenuation-only / +CD / +PMD / +birefringence all 1.16 % ±
1.16 (86 sifted — bit-identical: rate-only impairments for time-bin);
+visibility 0.934 → 3.53 % (the (1−V)/2 floor); +DCR 10 kHz → 13.08 %;
full chain 3.53 %. Attenuation-on rows are sample-limited at 2M pulses
(±1.2 %, explicit in the artifacts; ~28 µs/pulse runtime).

**Tests:** 203/203 pass (no component code touched).

---

## 2026-08-06 — Fifth pass: PHYS-5 — phenomenological birefringence model deleted

### Session: PHYS-5 (+ compensation-flag bug found during the model test)

#### Part 1 — The test that decided it

Before touching anything, the quasi-static multi-section model was tested at
long distances (test + results in `opto-sim-issues-and-fixes.md` §14):

- **Matrix level:** sectional `θ(L)` saturates to uniform SU(2) (mean ≈ 2.2 rad,
  Haar expectation) by ~5 km — the random walk scrambles within a few
  correlation cells (per-section retardance δ ≈ 10 rad wraps mod 2π). The
  phenomenological `θ = min(π, √(L/L_char)·π/2)` grows to the π clamp only by
  ~150 km — fitted constants (L₀ = 75 km, Δn₀ = 0.87e-5), no literature source.
  Unitarity err ~1e-13, power conserved to ~1e-12, 0.16 ms/apply at 122 km
  (N = 2440 sections) — no speed advantage for the old model.
- **Protocol level** (`bb84_duplinskiy`, 200k pulses): compensated QBER is
  model-independent (~1.7 % floor — the calibration loop undoes any unitary J);
  uncompensated QBER is realization-dependent with the sectional model
  (38.7 %/23.0 % for seeds 42/7 — a per-installation unknown, correct physics)
  vs the phenomenological deterministic 25.86 %.

#### Part 2 — The bug the test exposed

`--no-compensation` was a silent no-op: the per-pulse compensation was applied
whenever the fibre had a Jones matrix, ignoring the `compensate` parameter
(`if U_comp is not None` instead of `if compensate and U_comp is not None`).
Fixed; the ARCH-3 control value (25.86 % @ seed 42) is re-verified. Regression
test added in `tests/test_protocols.py` (comp 1.5 % vs uncomp 49.6 % @ 10 km,
μ = 2, sectional, seed 42).

#### Part 3 — The deletion

| Change | Files | Rationale |
|---|---|---|
| `_build_jones_phenomenological`, `_apply_birefringence_phenomenological`, dead `_random_su2_rotation_rng`/`_random_su2_rotation` removed; `SECTIONAL_LIMIT` gone | `src/channel/fiber.py` | Single model at all lengths; `model='auto'` ≡ `'sectional'`; `model='phenomenological'` raises ValueError with a pointer to PHYS-5. Docstrings updated. |
| 5 phenomenological tests → long-distance sectional; dispatch test → `auto ≡ sectional at all lengths` + `phenomenological raises` | `tests/test_fiber.py` | The former phenomenological tests (power conservation at 100 km, T/λ/seed dependence, output variation) still cover exactly the right properties, now for the one model that exists. |
| Panels B–D rewritten on the single model | `analysis/validation/validate_birefringence.py` | B: long-distance ensemble plateau (uniform SU(2), |E_x|² → ½ within ~200 m); C/D: temperature & bend sensitivity at 2 m — the single-correlation-cell regime where δ(T)/δ(R) is visible (beyond ~1 km the ensemble mean is saturated and these effects are invisible by physics, not by modelling). Poincaré convergence unchanged: |mean(S)| = 0.944 → 0.034. |
| `--birefringence-model` choices reduced; `compensate` flag fixed | `src/protocols/bb84_duplinskiy.py` | Choices now auto/sectional only. |
| Docs | `AGENTS.md`, `README.md` | Hybrid-dispatch sections → single model; test count 203. |

**Tests:** 203/203 pass (was 201; −6 +7 in `test_fiber.py`, +1 `test_protocols.py`).
Harness: `run_all.py --seed 42 --skip gobby` — 7/7 PASS (114.6 s).

**PHYS-5 outcome:** both temperature coefficients (`-3.0e-9` vs `-5e-7`) and
both clamp values (`5e-10` vs `1e-10`) collapse to the confirmed sectional
values. Manuscript items M2 (Eq. 15 coefficient) and M3 (clamp unification)
now have a single answer: `-3.0e-9 /°C` and `5e-10`. No published number
moves: time-bin (Gobby) never used birefringence; the Duplinskiy compensated
table is model-independent; only the uncompensated control becomes
seed-dependent (25.86 % @ seed 42 is still reproducible).

---

## 2026-08-06 — Fourth pass: ARCH-1 system validation rebuild + ARCH-3 polarization compensation

### Session: ARCH-1 + ARCH-3

#### Part 1 — ARCH-1: Section 5 rebuilt on the genuine time-bin SPAD chain

| Change | Files | Rationale |
|---|---|---|
| `analysis/val_system.py` fully rewritten | `analysis/val_system.py` | The old demo (CWLaser → polarizer → phase modulator → fiber → PBS → APD, 1000 km sweep) was classically invalid: it reported QBER from a CW (unmodulated) carrier. The rebuild is the time-bin chain used by the Duplinskiy/Gobby work — CWLaser → MZM carve → encoder AMZI → FiberRealization → decoder AMZI → 2×SPAD — with the real impairment set (CD, PMD, birefringence, attenuation) and ID230 SPAD statistics (η = 0.10, dead 13 µs, afterpulse 5 %, DCR 15 Hz, 1 ns gate). |
| Linearity shortcut | `analysis/val_system.py` | Per-bit field propagation costs ~16.7 ms, making per-point sweeps infeasible. The early/late basis pulses are propagated once per parameter point, and the per-bit probability is evaluated from the exact linearity identity `P_c = g0_c + 2·Re[S_c·e^(−jδ)]`, with `δ = φB − φA − θ` (θ phase jitter from 2πΔν·delay, θ_σ = 190.9 mrad). ~5 s per 1M-pulse point. |
| MZM X-cut polarisation fix | `analysis/val_system.py` | The MZM is X-cut (modulates Ey only); the source must therefore be Y-polarized (`polarization_azimuth = π/2`), otherwise the unmodulated Ex component forms a CW floor that destroys the μ calibration. Also: `S_c`/`S_d` must be built with `complex()` (a `float()` cast silently zeroed the imaginary part). |
| Panel set replaced | `analysis/val_system.py` | Old panels (QBER vs distance with bit-rate title, temperature, bend radius — flat for time-bin) → A: distance 0–122 km (all impairments vs attenuation-only); B: pulse σ 5–50 ps @ 75 km; C: decoder visibility 0.90–1.00 @ 75 km; D: μ 0.02–2.0 @ 75 km; E: DCR 0–10 kHz @ 122 km. Time-bin is immune to slow birefringence (both bins traverse the identical quasi-static operator); CD/PMD change rate, not QBER. |

**ARCH-1 results (1M pulses/point, seed 42):** 0 km QBER 5.4 %; 75 km floor ≈ 6.4 % (decomposition: visibility (1−V)/2 = 3.3 %, phase jitter ~0.8 %, afterpulse ~1.5 %, double-click ~0.5 %); C: 8.5 %→4.3 % across V; D: 0 % @ μ = 0.02 (sample-limited, 24 sifted) → 6.2 % @ μ = 2 (2022 sifted); E: 0 %→16.7 % @ 10 kHz DCR; 122 km sifted = 17 (sample-limited). Artifacts: `val_system/val_system--seed42.{png,csv}` regenerated.

#### Part 2 — ARCH-3: Duplinskiy polarization compensation

| Change | Files | Rationale |
|---|---|---|
| `FiberRealization.birefringence_matrix()` | `src/channel/fiber.py` | Public accessor returning a copy of the quasi-static SU(2) fibre matrix (None when birefringence is disabled). |
| Compensation in `bb84_duplinskiy.py` | `src/protocols/bb84_duplinskiy.py` | `compensate=True` (default) with `--no-compensation` control: `U_comp = J_channel.conj().T` applied per pulse between fibre and VOA — the fixed point of the paper's three-controller calibration loop for a quasi-static channel. Docstring rewritten (loop ⇔ inverse J). |
| Accessor tests | `tests/test_fiber.py` | `TestBirefringenceMatrixAccessor` (5): unitarity, matches `apply()` with other impairments off, roundtrip J†J restores the field, None when disabled, quasi-static (same matrix across calls). |

**ARCH-3 results (seed 42):** 0 km 2.83 % (460 sifted/13 err, 200k); 50 km compensated 1.72 % (58/1, 200k) and 0.98 % (307/3, 1M); 50 km uncompensated control 25.86 % (58/15, 200k). Paper (Duplinskiy et al., Opt. Express 25(23) 28886, 2017; arXiv 1709.06655): QBER ≈ 2 % at 50 km — replication matches within statistics, and the uncompensated control quantifies the necessity of the calibration loop.

**Tests:** 201/201 pass (was 196; +5 accessor tests). Harness: `run_all.py --seed 42 --skip gobby` — 7/7 PASS, 66.3 s.

---

## 2026-08-06 — Third pass: birefringence literature validation, BB84 script consolidation, README reproducibility

### Session: REPRO-4 + ARCH-2 + REPRO-2

#### Part 1 — REPRO-4: real literature comparisons for birefringence

| Change | Files | Rationale |
|---|---|---|
| `bend_birefringence(bend_radius, r_fiber=62.5e-6)` — public, documented | `src/channel/fiber.py` | The bend-induced birefringence `Δn_bend = 0.135·(r/R)²` (Ulrich [7] Eq. 1) was duplicated inline in both Jones builders with their own copies of the constants. Now a single cited function both builders call. |
| `TestUlrichBendLaw` (10 tests) | `tests/test_fiber.py` | (1) Exact match against the published `0.135·(r/R)²` law across a bend-radius sweep (rtol 1e-12, stated tolerance); (2) channel-level: the bend term recovered from the *single-section retardance* (eigenphase of the SU(2) matrix, axis-independent) and compared to the law, tolerances 5–18 % absorbing the model's 10 % stochastic residual. Fibre length is scaled per radius so the retardance is exactly π/2 — an SU(2) matrix only carries retardance mod 2π, so unwrapped recovery needs δ < π. |
| 13 self-consistency checks moved | `analysis/validation/validate_birefringence.py` → `tests/test_fiber.py` | The old script's 13 `[PASS]`-printing invariants (power conservation, temperature/wavelength/seed dependence, Poincaré variation, zero-length, auto-dispatch, `enabled=False`) are now `TestBirefringenceSelfConsistency` (13 tests). The validator keeps only what it validates: the 6-panel figure, Poincaré convergence figure + assertion, and the CSVs. Table CSV now records that the checks moved. |

#### Part 2 — ARCH-2: five BB84 implementations, three active

| Change | Files | Rationale |
|---|---|---|
| `bb84_ideal.py`, `bb84_high_bitrate.py` → `src/protocols/examples/` | moved files | Nothing in the repo imports them (grep-verified); they are legacy CW-based demos superseded by `bb84_time_bin.py`. The active set is now unambiguous: `bb84_time_bin.py` (Gobby validation), `bb84_test_dispersion.py` (dispersion study), `bb84_duplinskiy.py` (ARCH-3 replication). Smoke-tested: `python -m src.protocols.examples.bb84_ideal --fiber-length 5` → QBER 0.0. |
| Docs updated | `README.md`, `AGENTS.md` | Architecture trees and script descriptions point at the new locations. |

#### Part 3 — REPRO-2: README paths and invocations

| Change | Files | Rationale |
|---|---|---|
| Gobby path fixed | `README.md` | `python analysis/validation/validate_gobby.py` → `python analysis/val_gobby/validate_gobby.py` (the doc's exact complaint — the file never lived in `analysis/validation/`). |
| Test counts 77 → 173 | `README.md` | Two stale counts. |
| "Exact invocations behind published results" table | `README.md` | Documents that the published Gobby table used `--bits 10000000` (10M pulses/point; default 200k), the `run_all.py` equivalents, and that the birefringence self-consistency checks now live in `tests/test_fiber.py`. |

**Tests:** 196/196 pass (was 173; +13 self-consistency, +10 Ulrich bend law).

**Key results:**
- `bend_birefringence` matches Ulrich's law exactly (rtol 1e-12) and flows through the sectional Jones matrix: measured Δn_bend = 1.47e-5 at R = 6 mm vs law 1.46e-5 (within noise tolerance)
- Validator harness: 7/7 validators pass (Birefringence 12.4 s, unchanged runtime); Gobby validator unchanged

---

## 2026-08-06 — Second pass: Gobby table columns, MZI visibility, validation harness, untested-module coverage

### Session: BLOCK-1 + BLOCK-3 + REPRO-1 + REPRO-3 + REPRO-5

Continues the code-side work from `opto-sim-issues-and-fixes.md` (first pass
listed below). Manuscript items remain out of scope.

#### Part 1 — BLOCK-1: Gobby table presented the analytic column as measured

| Change | Files | Rationale |
|---|---|---|
| `analytical_qber()` → `model_qber()`; new `gobby_measured_qber()` | `analysis/val_gobby/validate_gobby.py` | The table's numeric column inherited the internal function name `analytical_qber` while being presented as Gobby's measured data — an integrity exposure. The measured values are now interpolated from the published Fig. 3 points (`GOBBY_DIST_KM = [4.4, 65, 101, 122]`, `GOBBY_QBER = [3.3, 3.3, 6.0, 8.9]`) via `np.interp`, so the columns are provably distinct. |
| Table regenerated | `analysis/val_gobby/val_gobby_table.tex` | Columns now read "This work analytic (%)" and "Gobby et al. measured (%)". |

#### Part 2 — BLOCK-3: `AsymmetricMZI` had no visibility or phase-error knob

| Change | Files | Rationale |
|---|---|---|
| `visibility` (default 1.0) + `phase_error` (default 0.0) | `src/channel/interferometer.py` | Modelled as a combiner amplitude imbalance (r, s with r² + s² = 1, 2rs = V) rather than scalar field mixing, so power is conserved at the fringe: P_c ∝ (1 + V·cos Δφ), P_d ∝ (1 − V·cos Δφ), giving a minimum-error floor e_opt = (1 − V)/2 at Δφ = 0. `phase_error` shifts the fringe peak. `visibility` validated to (0, 1]. |
| Params + CLI flags wired through | `src/protocols/bb84_time_bin.py` | `--visibility`, `--phase-error`. |
| 8 tests (`TestAsymmetricMZIVisibility`) | `tests/test_interferometer.py` | Visibility fringe contrast, error floor, power conservation, phase-error fringe shift, parameter validation, backward compatibility. |

**Verified — the audit's afterpulse hypothesis is now confirmed:**
- 0 km QBER with `afterpulse_prob=0`: **2.786 % → 0.0000 %** — the short-range
  floor is afterpulsing, not misalignment.
- `visibility=0.934`, no noise, 1.5M pulses: **3.2961 %** ≈ Gobby's 3.3 %
  short-range floor.
- Fringe contrast matches V exactly to 4 decimals; `visibility=1.0` is
  byte-identical to the old 50:50 combiner.

#### Part 3 — REPRO-1: eager tkinter import broke headless runs

| Change | Files | Rationale |
|---|---|---|
| `polarimeter()` is now a lazy wrapper; tkinter imported only on call | `src/visualization/__init__.py` | `import opto_sim` no longer needs a display server (CI/headless reproducibility). |

#### Part 4 — REPRO-3: `run_all.py` was UTF-16 and checked only exit codes

| Change | Files | Rationale |
|---|---|---|
| Re-encoded UTF-16LE (BOM) → UTF-8; rewritten as a harness | `run_all.py` | The file died with `SyntaxError: Non-UTF-8 code` before executing a line. The harness now runs all 8 validators (CD, PMD, attenuation, birefringence, APD, CW laser, MZM, Gobby), checks exit codes **and** the expected output file, supports `--skip` and `--gobby-bits`, and exits non-zero on any failure. The output-file check matters because `validate_apd`/`validate_cwlaser`/`validate_mzm` print no `[PASS]` marker and are assert-free. |

**Verified:** full harness green — 96 s without Gobby, 158 s with 20k Gobby bits.

#### Part 5 — REPRO-5: untested modules, and two defects found

| Added | File | Covers |
|---|---|---|
| 20 tests | `tests/test_spad.py` | Dead time, DCR convergence, Poisson click rate 1 − exp(−qe·µ), afterpulse excess-rate ~5 %, afterpulse requires a prior click, PHYS-7 no-leak regression (re-armed detector must not fire a stale pending afterpulse). |
| 35 tests | `tests/test_optics.py` | Unitarity (HWP/QWP/rotator/hadamard/circular_analyser), PBS vs analyser phase-blindness (PHYS-6 regression), projector idempotence, cascaded-polariser double extinction, VOA scaling, coupler split/combine bookkeeping. |
| 20 tests | `tests/test_phase_modulator.py` | V_π formula (Alferness [2], Weis & Gaylord [1]), X-cut/Y-cut modulation axis, DC phase application, RF per-sample phase, V_π caching, parameter validation. |

**Defects found and fixed by the new tests:**

1. **`PhaseModulator` partial-`params` crash** — passing a partial dict
   (e.g. `{'wavelength': ...}`) never merged with the defaults, so `get_vpi()`
   died with `AttributeError: 'n_o'`. Fixed by merging `{**defaults, **params}`.
   All production call sites pass `params=None`, so behaviour there is
   unchanged. (`src/channel/phase_modulator.py`)
2. **`coupler_combine` power doubling (legacy, unused)** — the 2-port combine
   applied `[[1, j], [j, 1]]` without the 1/√2 normalization, doubling output
   power. Both branches now apply the ideal 3 dB coupler scattering matrix
   (`E_out1 = (E1 + j·E2)/√2`, `E_out2 = (j·E1 + E2)/√2`), which is unitary, and
   returned powers are derived from the output fields per the project's
   field-derived-power convention. (`src/channel/optics.py`)

**Tests:** 173/173 pass (was 89 at session start; 97 after BLOCK-3; 171 after
REPRO-5; 173 after the coupler fix — the one `xfail` became a real pass).

**Key results:**
- BLOCK-3 premise confirmed: the 0 km QBER floor is afterpulsing, not misalignment
- Gobby's 3.3 % short-range floor reproduced as a visibility V = 0.934 → 3.2961 %
- V_π(X-cut) = 3.8826 V, matching the crystal-derived value all BB84 scripts use

**Still open (unchanged from the first pass):** BLOCK-2 (Table 11 script),
ARCH-1/2/3, PHYS-5, REPRO-2, REPRO-4 remainder (Ulrich bend law), §10, M1–M14.
Gobby validator still runs at default `visibility=1.0`; its short-range floor
now has a physically attributable knob (V ≈ 0.934 → 3.3 %) if the replication
target is the Gobby floor.

---

## 2026-08-06 — Audit fixes: quasi-static fibre, impairment separation, physics corrections

### Session: ROOT-1 + impairment separation + PHYS-1/2/3/4/6/7 + PERF-1

Addresses the code-side findings in `opto-sim-issues-and-fixes.md`. Manuscript
items (M1–M14) were explicitly out of scope and are untouched.

#### Part 1 — ROOT-1: birefringence was a per-bit depolarising channel

| Change | Files | Rationale |
|---|---|---|
| `FiberRealization` class | `src/channel/fiber.py`, `src/channel/__init__.py` | A real fibre's Jones matrix is quasi-static — it drifts on seconds-to-minutes timescales, not per 4 ns bit. `np.random.uniform(...)` was being called *inside* the propagation function, so every bit traversed an independently sampled random unitary, converting the fibre into a fully depolarising channel. The class builds the matrix once at construction from its own `np.random.Generator` (seeded independently of the global RNG used for bit/basis choices and detector noise) and reuses it for every `apply()`. |
| Jones-matrix builders factored out | `src/channel/fiber.py` | `_build_jones_sectional`, `_build_jones_phenomenological`, `_build_jones_matrix`, `_random_su2_rotation_rng` — all take an `rng` argument, accepting either the `numpy.random` module (stateless per-call API) or a `Generator` (FiberRealization). Keeps the stateless path byte-for-byte identical; the validation scripts legitimately want a fresh draw per ensemble sample. |
| Wired into all per-bit loops | `analysis/val_system.py`, `src/protocols/bb84_ideal.py`, `bb84_high_bitrate.py`, `bb84_test_dispersion.py`, `bb84_duplinskiy.py` | One `FiberRealization` built before the loop, reused per bit. |

**Verified:** five `apply()` calls on one realization are bit-identical; two seeds differ; `FiberRealization(L_m=0)` is the identity. Panel A's QBER-vs-distance no longer saturates smoothly at 50% — it now jumps between near-0% and near-100% depending on where each fibre's *fixed* rotation lands relative to Bob's basis, the expected signature of a compensable fixed unitary rather than a depolarising channel.

#### Part 2 — Impairment separation (requested feature, beyond the audit)

| Change | Files | Rationale |
|---|---|---|
| `FiberRealization` owns all four impairments | `src/channel/fiber.py` | `birefringence`, `cd`, `pmd`, `attenuation` are now independent constructor flags defaulting to `propagate()`'s existing defaults (`True, False, False, True`). Previously the class froze only birefringence while `propagate()` re-applied attenuation/CD/PMD separately on every call. `apply(E, dt=None)` runs the enabled set in the original order: birefringence → CD → PMD → attenuation. |
| PMD split into sampling + application | `src/channel/fiber.py` | `_sample_pmd_dgd(rng, ...)` (Maxwellian DGD + axis-swap draw) and `_apply_pmd_fixed(E, dt, dgd, swap)` (deterministic frequency-domain operator), mirroring the birefringence split. PMD's DGD is now frozen at construction alongside birefringence — same physical origin (fibre asymmetry), same quasi-static argument. Verified bit-identical to the old code under a seeded global RNG. |
| `propagate()` delegates wholesale | `src/channel/fiber.py` | Given a `fiber_realization`, returns `fiber_realization.apply(E, dt=dt)` and ignores its own impairment arguments (documented). Stateless path unchanged. |
| Call sites use `fibre.apply()` directly | `analysis/val_system.py`, 4× `src/protocols/bb84_*.py` | Each passes its impairment flags to the constructor; no more indirection through `propagate()` with dead arguments. |

#### Part 3 — Physics and performance corrections

| Fix | Files | Detail |
|---|---|---|
| **PHYS-1** — PMD coefficient 31.6× too large | `src/channel/fiber.py`, `analysis/val_system.py`, `analysis/validation/validate_pmd.py` | `pm_dispersion` (s/√m) → `pmd_coeff_ps_sqrt_km` (ps/√km, the datasheet convention), converted internally. Cited to new **ref [12] Corning SMF-28 Ultra PI1463 (2021)** — PMD link value ≤ 0.06, typical ≤ 0.1 ps/√km. `validate_pmd.py` now fits **0.09982 ps/√km** vs 0.1 nominal (was 3.16), R² = 0.999977. Panel D sweeps raw ps/√km, matching its axis label. |
| **PHYS-2** — `D_TOTAL` double-counted | `src/channel/fiber.py` | 17.0 was the *total* SMF-28 dispersion, not the material component; subtracting a waveguide term on top left `D_TOTAL` 18% low. `D_MATERIAL` 17.0 → **22.0**, `D_WAVEGUIDE` −3.0 → **−5.0** (**Agrawal [6] Fig. 2.10**), giving `D_TOTAL = 17.0` (**Corning [12] / ITU-T G.652**). Decomposition kept because `validate_cd.py` Panel E plots the components separately. `validate_cd.py` now prints `D = 17.0`, `β₂ = −2.168e-26 s²/m`, `L_D = 41.51 km`. |
| **PHYS-3** — sifting counted coin flips | `analysis/val_system.py` | Non-detections were assigned `random.randint(0,1)` and *kept* in the sifted key, so QBER converged to 50% by construction at long range. Now tracks `has_click` and sifts on `bases match AND detected`, matching `bb84_time_bin.py`/`bb84_duplinskiy.py`. `simulate_bb84_full()` returns `{'qber', 'sifted_fraction'}`; CSV gained a `Sifted_fraction` column. |
| **PHYS-4** — `section_length` conflated two lengths | `src/channel/fiber.py`, `analysis/validation/validate_birefringence.py` | Each section draws a *new random axis*, so the parameter is physically the birefringence correlation length L_c, not a discretisation step. Renamed `section_length` → `correlation_length`, default 1.0 → **50.0 m** (**Menyuk & Wai [10]**, physical range 10–100 m). |
| **PHYS-6** — `optics.pbs()` was not a PBS | `src/channel/optics.py`, `src/channel/__init__.py`, 5 call sites | The matrix `(1/√2)·[[1,−1j],[−1j,1]]` is a circular-basis projector (QWP+PBS), not H/V splitting. Renamed `circular_analyser()` and documented; a real `pbs()` doing H/V projection added alongside. All callers switched explicitly, so **behaviour is unchanged** — for pure-H input `circular_analyser` reproduces the old `(0.707, −0.707j)` while the new `pbs` gives `(1, 0)`. |
| **PHYS-7** — latent afterpulse state leak | `src/detectors/spad.py` | On re-arm the code left `_afterpulse_pending`/`_afterpulse_time` set, so a scheduled-but-never-fired afterpulse could fire on a *later* dead period. Both fields now cleared on re-arm, and `_schedule_afterpulse` always assigns so a failed roll clears rather than preserves. |
| **PERF-1** — O(N) matmul loop | `src/channel/fiber.py` | New `_ordered_product()` computes the ordered Jones product by pairwise tree reduction in O(log N) vectorised steps. Exact, not approximate — matrix multiplication is associative. |

#### Tests

| Added | File | Covers |
|---|---|---|
| `TestPMD` (1) | `tests/test_fiber.py` | Maxwellian mean DGD `2·(σ/√3)·√(2/π)` = 0.921 ps at 0.1 ps/√km over 100 km, against 20 000 samples. |
| `TestBirefringenceDepolarization` (1) | `tests/test_fiber.py` | Analytic depolarisation law `⟨S⟩ = pᴺ·S_in`, `p = (1+2cos α)/3` (**Menyuk & Wai [10]**), 400 realizations at L = 10/50/100/200 m. |
| `TestOrderedProductTreeReduction` (10) | `tests/test_fiber.py` | Tree reduction vs naive left-fold for N = 1…1001, covering odd/even and power-of-two boundaries. |

**Tests:** 89/89 pass (was 77).

**Key results:**
- `validate_cd.py`: D = 17.0 ps/(nm·km), L_D = 41.51 km, max error 3.01e-14 %
- `validate_pmd.py`: fitted 0.09982 ps/√km (nominal 0.1), mean DGD 0.933 ps (expected 0.921), KS p = 0.17
- `validate_attenuation.py`: α = 0.182 dB/km, R² = 1.0000000000, max error 3.21e-14 %
- `validate_birefringence.py`: 13/13 self-consistency checks pass; Poincaré convergence |mean(S)| = 0.944 → 0.034
- PHYS-7 leak, measured against pre-fix code: **2000/2000 trials leaked before, 0/2000 after**; effective afterpulse rate still 4.2%/5.1% vs 5% nominal
- All five protocol scripts and `val_system.py` run clean end-to-end

**Knock-on fix.** `validate_birefringence.py`'s convergence demo carried a hardcoded assertion (`|mean(S)| > 0.15` at L = 10 m) encoding the *old* 1 m parameter, which failed under PHYS-4. With L_c = 50 m > L_B ≈ 31 m, a sub-correlation-length fibre is a single random-axis section whose retardation phase wraps past 2π within a few metres, so `|mean(S)|` oscillates with L rather than decaying monotonically. First demo point moved 10 m → 2 m (well under one beat length), restoring a coherent state. Reasoning documented in-script.

**Deliberately not done — see `opto-sim-issues-and-fixes.md` for the full list.** The phenomenological birefringence model was **retained**, so `SECTIONAL_LIMIT` is still 2 km and every run past 2 km still uses its fitted `Δn₀ = 0.87e-5`. ROOT-1's promised simplification (single `Δn₀`, no fit parameter, no 2 km discontinuity) is therefore not yet realised, and **PHYS-5 did not resolve itself as the audit predicted** — both temperature coefficients (`-3.0e-9` sectional, `-5e-7` phenomenological) are still live. Deleting that model will move all long-distance numbers and warrants an explicit decision. Also still open: BLOCK-1/2/3, REPRO-1/2/3/5, REPRO-4 (partial — `p^N` law added, Ulrich bend law still missing), ARCH-1/2/3.

**Separately noted:** `run_all.py` is saved as UTF-16LE with a BOM and dies with `SyntaxError: Non-UTF-8 code` before executing a line. Pre-existing, unrelated to these changes, and a prerequisite for REPRO-3.

---

## 2026-07-22 — Time-bin phase-encoding BB84 (Gobby et al. 2004 replication)

### Session: Pulsed laser, AsymmetricMZI, time-bin protocol, Gobby validation, manuscript update

| Change | Files | Rationale |
|---|---|---|
| Switchable pulsed mode for CWLaser | `src/lasers/cwlaser.py` | `pulsed=False` (default, backward-compatible). Adds `pulse_width`, `repetition_rate`, `timing_jitter_rms`. Gaussian pulse train with `⟨\|g\|²⟩=1` preserves average power. Literature: Agrawal §3.4, Gobby 2004. |
| AsymmetricMZI encoder/decoder | `src/channel/interferometer.py` (new), `src/channel/__init__.py` | 50:50 Hadamard split/combine, delay via `np.roll`, phase on delayed arm, insertion loss. Encoder: single pulse → two time bins. Decoder: recombine → constructive + destructive ports with `cos²(Δφ/2)` fringe. Literature: Townsend 1993, Bennett & Brassard 1984. |
| Time-bin BB84 protocol | `src/protocols/bb84_time_bin.py` (new) | Full chain: pulsed Gaussian pulse → AMZI encoder → fiber (attenuation only; time-bin immune to birefringence) → AMZI decoder → 2× SPAD. X/Y basis encoding, sifting, QBER computation. |
| Gobby et al. 2004 validation | `analysis/val_gobby/validate_gobby.py` (new) | QBER vs distance sweep (0–122 km). Monte Carlo + analytical overlay against Gobby Fig 3 data. 0 km QBER = 3.19% (paper: ~3.3%). |
| AMZI tests (21) | `tests/test_interferometer.py` (new) | Construction, encoder (output shape, two time bins, delay, power conservation, phase), decoder (interference fringes, constructive/destructive dominance, power sum), roundtrip, insertion loss. |
| Pulsed laser tests (8) | `tests/test_cwlaser.py` | Shape, power conservation, FWHM, rep-rate, zero-power inter-pulse, energy/pulse, jitter. |
| Manuscript updated | `paperwork/manuscript.tex`, `paperwork/tables/val_gobby_table.tex` | New §2.4 (AsymmetricMZI), §2.7 (SPAD), §3.8 (Time-bin BB84 validation). Updated abstract, contributions, reproducibility (48→77 tests), summary table, conclusion. 38 references. 28 pages. |

**Key results:**
- 0 km: QBER = 3.19% (paper: ~3.3%) — baseline validated
- 4–40 km: QBER = 1–3% (consistent with Gobby plateau)
- Longer distances need more pulses for statistical significance
- Performance: ~2600 pulses/s (3rd Gen i5)
- 77/77 tests pass

---

## 2026-07-22 — Gobby validation: 10M-pulse sweep, manuscript update

### Session: Final Gobby sweep with 10M pulses per distance, manuscript refined

| Change | Files | Rationale |
|---|---|---|
| Gobby distance sweep at 10M pulses | `analysis/val_gobby/` — updated table, figure | 9 distances, 10M pulses each. Resolved stochastic noise >80 km. QBER: 2.43% (0 km) → 4.55% (122 km). Sifted bits: 23,546 → 154. |
| Manuscript updated | `paperwork/manuscript.tex` | Removed "statistical noise" language from Gobby section; added discussion of systematic offset vs paper (uncharacterised environmental noise). Updated all QBER references (3.19% → 2.43% at 0 km). |
| Table path fixed | `paperwork/tables/val_gobby_table.tex` | Copied from analysis output directory. |
| Manuscript compiled | `paperwork/manuscript.pdf` | 29 pages, all cross-references resolved. |

## 2026-07-21 — SPAD detector, VOA, Duplinskiy et al. BB84 protocol

### Session: single-photon detection, protocol replication, distance sweep

| Change | Files | Rationale |
|---|---|---|
| Geiger-mode SPAD detector | `src/detectors/spad.py` (new), `src/detectors/__init__.py` | Inherits from APD; same photoelectric physics, different operating point (bias above breakdown). Adds dead time (13 μs), dark count rate (15 Hz), afterpulsing (5%), gated detection (20 ns). Binary click output. ID230 specs. |
| Variable optical attenuator | `src/channel/optics.py` | `voa(E, attenuation_dB)` — applies amplitude scaling `sqrt(10^(-dB/10))` to field. Used for Bob's internal loss in Duplinskiy replication. |
| Duplinskiy et al. BB84 protocol | `src/protocols/bb84_duplinskiy.py` (new) | SPAD-based BB84 example. 0 km back-to-back QBER = 2.6% validates detection chain against paper ~2%. Distance sweep removed — uncompensated birefringence is not comparable to paper's active polarization tracking. |
| Manuscript updates | `paperwork/manuscript.tex`, `paperwork/tables/val_system_table.tex` | Removed old model references, added attenuation dominance explanation, updated system-level scenarios. |
| Reference PDF | `paperwork/Low_loss_QKD_optical_scheme_for_fast_polarization_.pdf`, `paperwork/Image.png` | Duplinskiy et al. paper and optical scheme diagram for reference. |

**Key results:**
- 0 km: QBER = 2.63% (dark count floor, matches paper's ~2% lab baseline)
- 50 km: QBER = 26.45% (paper: ~2% lab spool — discrepancy from phenomenological birefringence model)
- 100 km: QBER = 43.90% (approaching dark-count-limited 50%)
- Loss model validated: perfectly linear at 0.2 dB/km + 2 dB Bob

**Note:** SPAD physics same as APD (photoelectric effect, Poisson statistics) — only voltage threshold changes. No separate physics validation needed.

## 2026-07-20 — Birefringence recalibration (SMF-28), fiber.py → fiber_sectional.py

### Session: Recalibrate Δn₀ to SMF-28 literature, rename file

| Change | Files | Rationale |
|---|---|---|
| Recalibrated birefringence to SMF-28 | `src/channel/fiber.py` → `src/channel/fiber_sectional.py` | Δn₀ = 8.7e-6 → 5.0e-8 (L_B = 31 m, Agrawal §4.1 SMF-28 range 10–100 m). T_coeff = -5e-7 → -3e-9/°C (same ~6%/°C ratio). Clamping floor 1e-8 → 5e-10. Stochastic residual 10% of Δn₀. |
| Renamed fiber.py → fiber_sectional.py | `src/channel/fiber_sectional.py` | Multi-section model now has explicit filename. All imports updated. |
| Updated imports | `src/channel/__init__.py`, 5 validation scripts, `test_fiber.py` | s/src.channel.fiber/src.channel.fiber_sectional/g |
| Updated manuscript parameter values | `paperwork/manuscript.tex` | Δn₀ 8.7e-6 → 5e-8, L_B 0.18 m → 31 m. Removed "compressed correlation length" language — now SMF-28 realistic. |
| Updated validation Panel D units | `analysis/validation/validate_birefringence.py` | Beat length in metres (not mm) for new L_B ≈ 31 m. |
| Updated documentation | `AGENTS.md`, `journal_paper_outline.tex`, `CHANGELOG.md` | Fiber references point to fiber_sectional.py. |

**Tests:** 48/48 pass.

## 2026-07-20 — Rename fiber_sectional → fiber, cable → propagate, individual impairment flags, rework birefringence validation

### Session: renames, independent impairment toggles, split-model validation

**File rename:** `src/channel/fiber_sectional.py` → `src/channel/fiber.py`
- Updated all imports: `__init__.py`, `test_fiber.py`, 4 validation scripts
- Updated references in `journal_paper_outline.tex` and `AGENTS.md`

**Function rename:** `cable()` → `propagate()` in `src/channel/fiber.py`
- Updated all call sites: `analysis/val_system.py`, `tests/test_fiber.py`, all 3 BB84 protocol scripts, `main.py`
- Updated `__init__.py` export
- Updated AGENTS.md and CHANGELOG.md

**Independent impairment flags in `propagate()`:**
- `birefringence` (default True), `cd` (default None → uses `dispersion`), `pmd` (default None → uses `dispersion`), `attenuation` (default True)
- Backward compatible: `dispersion=True` still enables cd+pmd; `dispersion=False` (default) disables them
- New patterns: `propagate(..., birefringence=False, attenuation=False)` → no impairments;
  `propagate(..., birefringence=False)` → attenuation only;
  `propagate(..., attenuation=False)` → birefringence only
- `apply_birefringence()` also gets `enabled=True` parameter
- `dt` required only when `cd=True` or `pmd=True` (not just `dispersion=True`)

**Birefringence validation reworked (`validate_birefringence.py`):**
- Now validates both models explicitly: sectional (L < 2 km) and phenomenological (L ≥ 2 km) in separate test functions
- Added auto-dispatch test: verifies model selection at boundary
- Self-consistency checks now labelled per-model

**Tests:** 48/48 pass, validation clean.

## 2026-07-20 — Hybrid birefringence dispatch (sectional + phenomenological)

### Session: dual-model dispatch, bug fixes, performance tuning

**Hybrid birefringence model:**
- `apply_birefringence()` dispatches automatically via `model='auto'`:
  - **Short fibres** (L < `SECTIONAL_LIMIT` = 2 km, `model='sectional'`): multi-section ordered product of random-axis SU(2) matrices, L_B ≈ 31 m (Agrawal §4.1). For DV-QKD and DPS QKD.
  - **Long fibres** (L >= 2 km, `model='phenomenological'`): single SU(2) rotation with `θ = min(π, √(L/L_char)·π/2)` (Menyuk & Wai 1994), L₀ = 75 km. For long-haul BB84 with distance-dependent QBER.
  - Auto-dispatch required because multi-section model converges to uniform SU(2) within ~1 km regardless of parameters, producing flat ~50% QBER.
- Fixed: duplicate return statement, `section_length` default to 1.0, `propagate()` indentation, unused `L0` in validation

**System demo results (hybrid mode):**
- 0–70 km: 0% QBER → 80 km: 23% → 200 km: 68% (peak) → 500+ km: ~45–53% (dark count floor)
- Temperature: V-shaped null at ~35–50°C; pulse width: 21.7% (5 ps) → 1.7% (30 ps); bend radius: 69.2% (2 mm) → 19.6% (5 cm)

**Performance:** `section_length=1.0` → 67 ms/call (2 hrs sweep), `section_length=100.0` → 1 ms/call (2 min sweep). Sectional model defaults to 100 m sections.

## 2026-07-19 — Random-axis birefringence model, system demo 0–1000 km, manuscript 21 pp

### Session: random birefringence + extended demo + manuscript

**Random-axis birefringence model:**
- `src/channel/fiber.py`: new `_random_su2_rotation()` and updated `apply_birefringence()`:
  - SU(2) rotation around uniformly random Poincaré-sphere axis (per-bit axis varies)
  - Rotation angle follows diffusive walk: `θ = min(π, √(L/L_char)·π/2)`
  - Characteristic length: `L_char = L₀·(Δn₀/|Δn|)²`, `L₀ = 75 km` at base Δn = 1.2e-7
  - Temperature (`T_coeff = -5e-7/°C`) and bend radius (`0.135·(r_clad/R)²`) modulate Δn → change L_char → scrambling rate
  - Literature: Menyuk & Wai (1994), Wai & Menyuk (1996) for diffusive polarization model
- `analysis/validation/validate_birefringence.py`: rewritten for random model — 6 self-consistency checks (power conservation, zero-length identity, temp/wavelength dependence, seeded reproducibility, polarization variance at long distance)

**System-level demo extended:**
- `analysis/val_system.py`: distance sweep 0–200 km (10 km steps) + 250–1000 km (50 km steps) = 37 points
- Three regimes: 0% QBER (0–80 km) → peak ~68% at 200 km (birefringence-dominated) → decay to ~48% dark-count floor (500–1000 km)
- Effective bit rate: 250 MHz (4 ns window per bit, 4000 samples × 1 ps); cited to Takesue (2007, 1.6 GHz BB84 over 200 km) and Dixon (2008, GHz-clock QKD)

**Manuscript updated:**
- Section 3.3 (Birefringence): describes SU(2) random rotation, diffusive angle, per-bit axis variation
- Section 4 (Birefringence Validation): rewritten for new validation figure (6 panels)
- Section 5 (System-Level Demo): three-regime QBER description, 0–1000 km, literature citations
- Model Limitations #1: phenomenological birefringence noted
- 3 new bibliography entries: Gobby (2004, APL), Takesue (2007, Nature Photonics), Dixon (2008, APL); total 37 references
- Manuscript now 21 pages (up from 18), compiles cleanly

**Output files:**
- `val_system/val_system--seed42.png` (265 KB, 200 DPI) — full 1000 km sweep with bit-rate title
- `val_system/val_system--seed42.csv` — QBER data at all 37 distances
- `analysis/val_birefringence/val_birefringence--seed42.png` — updated 6-panel validation figure
- `analysis/val_birefringence/val_birefringence--seed42.csv` — validation data

---

## 2026-07-12 — LaTeX compilation fix: mdframed → colorbox + scope clarification

### Session: LaTeX fix + project scope broadened

**Background:** `journal_paper_outline.tex` failed with "Not in outer par mode" on every `\begin{table}[!ht]`.

**Root cause:** [`mdframed`](https://ctan.org/pkg/mdframed) patches `\@xfloat` using `\color@vbox`, which was removed from the LaTeX kernel in 2024–2025.

**Fix:**
- Removed `mdframed` + `caption` packages
- Replaced abstract box: `\mdframed{...}` → `\colorbox{highlightblue!5}{...}`
- `journal_paper_outline.pdf` now compiles (13 pages, xelatex)

**Scope clarification:** The project is a general-purpose physical-layer fiber-optic simulator, not QKD-specific. BB84 is one protocol on top. AGENTS.md updated accordingly.

---

## 2026-07-08 — Literature verification fixes + Ulrich bend model

### Session: literature verification (commit 25cf9b9)

| Fix | File | Detail |
|-----|------|--------|
| 1 — Stokes S3 clip | `src/visualization/stokes.py` | Added `np.clip(S3, -1.0, 1.0)` before `np.arcsin()` — prevented NaN at near-circular polarization |
| 2 — Symmetric Jones + standard Δβ | `src/channel/fiber.py` | Switched to `diag(exp(±j·Δβ·L/2))` with Δβ = 2π·Δn/λ (Agrawal [6] Eq 4.1.2). Relative phase = Δβ·L = 2π·L/L_B as before. |
| 3 — PMD sign randomized | `src/channel/fiber.py` | `if np.random.rand() < 0.5: Hx, Hy = Hy, Hx` — 50:50 fast/slow axis per realization |
| 4 — Vπ doc + refs | `src/channel/phase_modulator.py` | Added numbered refs [1] Weis & Gaylord 1985, [2] Alferness 1988; clarified Vπ is MZM push-pull effective value |
| 5 — Fields title | `src/visualization/fields.py` | Changed to `plt.suptitle()` spanning all subplots |
| — Yuan ref corrected | `src/channel/fiber.py` | Vol 27, 2019 → vol 24, no. 2, pp. 1062-1071, 2016 + DOI |

### Session: Ulrich bend model (this session)

**Bend model rebuilt — `num_bends` → `bend_radius`:**

| File | Change |
|------|--------|
| `src/channel/fiber.py` | Replaced `num_bends=0` (Yuan stress rod factor 2.4e-4) with `bend_radius=None` using Δn_bend = 0.135·(r_fiber/R)². References [7] Ulrich 1980, [8] Smith 1980, [9] Shibata 1986. |
| `analysis/validation/validate_birefringence.py` | Sweeps `bend_radius` (2 mm–2 cm), fits Δn vs (r_f/R)². 0.0000% error on slope coefficient 0.135. |
| `analysis/validation/validate_cd.py`, `validate_pmd.py`, `validate_attenuation.py` | `num_bends=0` → `bend_radius=None` |
| `tests/test_fiber.py` | `num_bends=0` → `bend_radius=None` |
| `src/protocols/bb84_*.py` | `num_bends=10` → `bend_radius=None` |
| `literature_verification_report.md` | Updated table, issue summary, recommendations (all fixed). |
| `validation_report.md` | Updated Fix 2 description. |
| `AGENTS.md` | Updated parameters, birefringence description, removed known issue. |

**Cumulative status:** All 5 literature verification issues resolved. All 4 Tier-1 validations pass with 0.0000% error. Bend model blocker eliminated.

---

## 2026-07-02 — Tier 1 Channel Validation

### Session: ~11:00–11:30 UTC+5

**Restructured validation into `analysis/validation/` with per-task output dirs:**
- `analysis/validation/validate_cd.py` → outputs to `analysis/val_cd/`
- `analysis/validation/validate_pmd.py` → outputs to `analysis/val_pmd/`
- `analysis/validation/validate_attenuation.py` → outputs to `analysis/val_att/`
- `analysis/validation/validate_birefringence.py` → outputs to `analysis/val_biref/`

**Task 1 — CD validation (Agrawal [6] Fig 2.6, Eq 3.2.6):**
- 30 ps Gaussian pulse at z/L_D = 0.0, 0.5, 1.0, 2.0
- RMS width ratio vs analytic: **0.0000% error** at all points
- Output: `analysis/val_cd/val_cd--seed42.{png,csv}`

**Task 2 — PMD validation (Razavi [5] Fig 2.11):**
- Fixed `fiber.py:145`: `np.random.rayleigh` → `scipy.stats.maxwell.rvs`
  Old: Rayleigh (2D) gave RMS = 1.128 × target
  New: Maxwell (3D) gives RMS = target exactly
- 20k realizations: RMS DGD = 22.36 ps (target 22.36 ps), mean = 20.62 ps
- KS test p = 0.82 — data consistent with Maxwellian
- All 12 fiber tests continue to pass

**Task 3 — Attenuation validation (Keiser [1] Eq 3.6, SMF-28):**
- Sweep 0–200 km, 0.182 dB/km @ 1550 nm
- **0.0000% error** at all distances

**Task 4 — Birefringence validation (Yuan [4] Fig 1):**
- Physics-based bend birefringence: Δn_bend = 8.762e-10 / R²
- L_B vs R for Δn₀ = {0.5, 1.0, 2.0} × 10⁻⁵
- Confirmed L_B ∝ R / Δn₀ scaling in bend-dominated regime
- **Finding**: fiber.py uses fixed `bend_effect_factor = 2.4e-4` (not R-dependent);
  correct Yuan model requires 1/R² scaling — flagged for future update

**Task 5 — Reproducibility pass:**
- `run_all.py` at project root runs all 4 validations with `--seed N`
- All outputs tagged with seed in filename: `{name}--seed{N}.{ext}`
- Total runtime: ~17.5s for all 4 scripts

## 2026-06-24 — QBER vs distance dispersion graph

### Session: ~16:00–16:10 UTC+5

**New analysis script — `analysis/qber_vs_distance_dispersion.py`:**
- Sweeps fiber length 10–200 km with 5 ps MZM-carved pulses (300 bits/point).
- Two curves: dispersion ON (CD+PMD active) and dispersion OFF (baseline).
- Dispersion OFF: flat 0% QBER at all distances.
- Dispersion ON: QBER climbs from ~0% at 10 km to ~42% at 200 km (CD + PMD accumulate with distance).
- Output: `analysis/qber_vs_distance_dispersion.png`

## 2026-06-24 — BB84 migration to `sample_field` + dispersion test

### Session: ~15:30–16:00 UTC+5

**BB84 scripts now use `sample_field()` (not `instantaneous_field`):**
- `bb84_ideal.py`, `bb84_high_bitrate.py`: `alice_laser.instantaneous_field(normalize=False, over_period=True)` → `alice_laser.sample_field(dt=1e-12, n_samples=1000)`
- This returns the complex envelope over 1 ns (one bit at 1 Gbaud), unblocking chromatic dispersion and PMD.
- `--dispersion` CLI flag added to both scripts (default False for backward compatibility).
- `dt=1e-12` is now passed to `cable()` so that the FFT frequency grid is valid.

**New file — `src/protocols/bb84_test_dispersion.py`:**
- MZM-carved Gaussian pulses for broadband field generation (5–30 ps σ).
- Laser initialized Ey-only so X-cut MZM modulates the entire field.
- MZM biased at V_pi (V=0 → null, V=V_pi → full transmission).
- Higher laser power (+10 dBm) compensates for low pulsed duty cycle.
- `dispersion=True` by default (this is the test file's purpose).
- CLI flags: `--pulse-sigma`, `--short-pulse`, `--no-dispersion`.

**CD/PMD now produces measurable QBER (100 km, seed=42):**
| Pulse σ | dispersion | QBER |
|---|---|---|
| 30 ps | False | 0.00 % |
| 30 ps | True  | 0.00 % (PMD < pulse width) |
| 5 ps  | True  | 15.00 % (z/LD ≈ 87, PMD >> pulse) |

## 2026-06-24 — `get_electric_field` → `instantaneous_field`

### Session: ~15:00–15:10 UTC+5

**Rename:**
- `CWLaser.get_electric_field()` renamed to `CWLaser.instantaneous_field()` to make its purpose unambiguous — returns the full optical field over one ~5 fs period for fast single-bit polarisation/phase validation.
- `instantaneous_field` docstring explicitly warns: NOT for CD, PMD, or baud-rate physics. Use `sample_field()` instead.
- `sample_field` docstring updated to note that `instantaneous_field` is available for quick validation.

**All callers updated:**
- `src/protocols/bb84_ideal.py`, `src/protocols/bb84_high_bitrate.py`
- `analysis/laser_characterization.py`
- `tests/test_cwlaser.py`
- `src/channel/fiber.py` (docstring reference)
- `AGENTS.md`, `README.md`

**Unchanged (legacy, still uses `get_electric_field`):**
- `src/deprecated/sslaser.py` — SolidStateLaser's method kept as-is.
- `main.py`, `main.ipynb`, `scripts/` — all use SolidStateLaser, not CWLaser.

---

## 2026-06-24 — Repository restructure & Tier 0 (testing/reproducibility)

### Session: ~10:30–11:45 UTC+5

**Structural changes:**
- `src/__init__.py` removed — `src` is now a namespace root, not a package (PEP 420).
- `src/opto_eq/` → `src/channel/` — clearer name for optical channel components.
- `src/viewers/` → `src/visualization/` — descriptive name for plotting utilities.
- `src/protocols/examples/` flattened → `src/protocols/` — BB84 scripts are the main protocols, not examples.
- `src/lasers/sslaser.py`, `src/lasers/ndyag.py` → `src/deprecated/` — broken/unused lasers, out of sight but kept.
- `src/lasers/__init__.py` now only exports `CWLaser`.
- 7 root-level loose scripts moved to `scripts/` (except `main.py`).
- `opto-sim.rar` deleted.

**Tier 0 — Testing & reproducibility:**
- `tests/` directory created with `conftest.py` (auto-seeds `random` + `np.random`, `--seed` CLI arg).
- 48 unit tests across 4 files:
  - `test_cwlaser.py` (11): power convention, phase noise, RIN scaling, seeded reproducibility.
  - `test_mzm.py` (13): Vpi, null/peak, quadrature bias, push-pull vs single-drive, insertion loss.
  - `test_fiber.py` (10): attenuation, birefringence, temperature, CD power conservation.
  - `test_apd.py` (11): responsivity, noise scaling, detect_photons, thermal floor.
- `--seed` CLI argument added to `bb84_ideal.py` and `bb84_high_bitrate.py`.
- `pytest>=8.0` added to `requirements.txt`.

---

## 2026-06-04 — Physics-informed detector overhaul & bug fixes

### Session: 11:30–11:45 UTC+5

---

### 1. `src/channel/fiber.py`

| Item | Detail |
|---|---|
| **File** | `src/channel/fiber.py` |
| **Total lines** | 119 (was 101) |
| **Change type** | Edit |

**Change A — Literature sources (lines 1–16)**
- **Old:** Lines 1–4: three bare comment lines (`# Gerd Keiser Book Chapter 3`, `# Behzad Razavi`, `# Thorlabs`)
- **New:** Lines 1–16: formatted literature block citing Keiser [1], Hui [2], Keck [3], Yuan [4], Razavi [5], Agrawal [6]

**Change B — Attenuation formula (lines 105–117)**
- **Old** (was line ~94): `pout = pin / (10**(-attenuation_factor * fiber_length / 10))`
  - Bug: `pin / (10^(-αL/10))` = `pin * 10^(+αL/10)`, so power *increased* with loss.
- **New** (lines 105–117):
  ```python
  att_lin = 10 ** (-attenuation_factor * fiber_length / 10)
  pout = pin * att_lin
  E = E * np.sqrt(att_lin)
  ```
  - Fix: `pout = pin * 10^(-αL/10)` for correct physical loss (Keiser [1] Eq 3.6).
  - Addition: electric field scaled by `sqrt(att_lin)` so `|E_out|² / |E_in|² = att_lin`, keeping the field consistent with the power budget.

---

### 2. `src/detectors/apd.py`

| Item | Detail |
|---|---|
| **File** | `src/detectors/apd.py` |
| **Total lines** | 129 (was 96) |
| **Change type** | Rewrite |

**Complete rewrite of the APD detector class:**

| Aspect | Old | New | Reason |
|---|---|---|---|
| **`__init__` signature** | `(self, wavelength, excess_noise_factor, load_resistance, temperature, gain=12, frequency=3e8/1550e-09, quantum_efficiency=0.9, dark_current=10e-6)` | Removed `frequency` parameter | Frequency is derived from `c/λ` — it is not a user-tunable knob. The old code's `frequency=40` in BB84 scripts was a unit-compensation hack, not physical. |
| **Duplicate constants** | `self.charge = 1.6e-19` and `self.q = 1.602e-19` (both e) | Single `self.charge = 1.602e-19`. `self.q` removed. | Prevent numerical inconsistency. |
| **`detect_photons()`** | `mpn = field_energy/h * frequency` → units of `1/m³`; then `mpn * area * exposure_time` gives `s/m` (not dimensionless). Missing factor of `c`. | `power / (h·ν) * exposure_time * η` where `power` is the actual optical power in Watts. Units: `(W / J) · s = 1` (dimensionless). | Old formula was dimensionally wrong. New formula follows Agrawal [2] Eq 4.1.2 and Saleh & Teich [3] Eq 17.1-10. |
| **`detect_photons()` args** | `(self, field_energy, area, exposure_time=1e-9)` | `(self, power, exposure_time)` | Power passed directly in Watts; no need to derive from field energy density. |
| **`calculate_output_current()`** | `I_signal = self.gain * self.R * power` (correct math), accepted unused `frequency` arg | Removed `frequency` arg. Same formula, now `(self, power)`. | `frequency` was dead code. |
| **`calculate_noise()`** | `I_noise = sqrt(i_d² + i_q² + i_th²) * F` — excess noise factor F applied to *all* noise terms including thermal. | `I_noise = sqrt(F·(i_d² + i_q²) + i_th²)` — F applies only to shot-noise terms (Kasap [1] Eq 4.45). | Physical correction: thermal noise is not multiplied by the avalanche excess noise factor. |
| **`output()`** | `(self, E, area, bandwidth, details=False)` — took electric field array, returned `detected_photons` (int) | `(self, power, bandwidth, area=1, exposure_time=None, details=False)` — takes optical power in Watts, returns noisy `I_total` (float, Amperes) | Current-based output is physically meaningful for receiver design. Includes Gaussian noise realization `I_total = normal(I_signal, I_noise)`. |
| **Photon noise sampling** | Clipped negative counts but didn't add noise to output current | Always returns `I_total` with additive Gaussian noise; photon count still available via `details=True` | Realistic receiver output includes both signal and noise. |
| **Literature citations** | None | Lines 3–9: Kasap [1], Agrawal [2], Saleh & Teich [3] | Traceability for every formula. |

---

### 3. `src/detectors/apd_v2.py`

| Item | Detail |
|---|---|
| **File** | `src/detectors/apd_v2.py` |
| **Change type** | Deleted |

- `apd_v2` was generated by an external coding agent, not by the repository owner.
- It was never exported in `src/detectors/__init__.py` (which only lists `apd`).
- The `bb84_high_bitrate.py` script that imported it was updated to use the regular `apd` instead (see below).

---

### 4. `src/visualization/stokes.py`

| Item | Detail |
|---|---|
| **File** | `src/visualization/stokes.py` |
| **Total lines** | 77 (was 66) |
| **Change type** | Edit |

| Lines | Change |
|---|---|
| 4–10 | Added literature block: Collett [1], Hecht [2], Born & Wolf [3] |
| 27–37 | Reordered computations: `S0` computed first, guard `if S0 == 0` moved before division. Added inline citations to Collett [1] Eq 2.12–2.15 |
| 39–43 | **Removed** lines `psi = 0.5 * arctan2(S2, S1)` and `chi = 0.5 * arcsin(S3/S0)`. The `chi` line was dividing the already-normalized `S3` by `S0` again, causing `|S3/S0| >> 1` → `arcsin(NaN)`. This code was dead (computed but never returned). Replaced with a commented block noting the correct formula and the clipping needed to prevent floating-point NaNs. |

---

### 5. `src/visualization/fields.py`

| Item | Detail |
|---|---|
| **File** | `src/visualization/fields.py` |
| **Total lines** | 37 (was 29) |
| **Change type** | Edit |

| Lines | Change |
|---|---|
| 1–6 | Added `import numpy as np` and `import matplotlib.pyplot as plt` (were missing, causing `NameError` at runtime). Added literature citation Hecht [1]. |

---

### 6. `src/lasers/cwlaser.py` (NEW FILE)

| Item | Detail |
|---|---|
| **File** | `src/lasers/cwlaser.py` |
| **Total lines** | 199 |
| **Change type** | Created |

New physics-informed continuous-wave laser model for QKD simulation:

| Feature | Implementation | Source |
|---|---|---|
| Phase noise | Wiener process with diffusion coefficient `D_φ = 2π·Δν` | Henry [1] Eq 18 |
| RIN | White noise low-pass filtered at `rin_bandwidth` | Kikuchi [5], Coldren [2] Ch. 5 |
| Polarization | Jones vector from azimuth `ψ` and ellipticity `χ` | Yariv [3] Ch. 6 |
| Power | User-specified in dBm, converted to Watts | Steady-state CW model |
| Field output | `get_electric_field(t, over_period, normalize)` — same interface as `SolidStateLaser` | Backward compatibility |

**Literature sources in file headers:**
- [1] Henry, C. H., "Theory of the Linewidth of Semiconductor Lasers", IEEE JQE 1982
- [2] Coldren, Corzine & Mashanovitch, "Diode Lasers and Photonic Integrated Circuits", 2nd ed., Wiley 2012
- [3] Yariv, A., "Optical Electronics", 4th ed., Saunders 1991, Ch. 6
- [4] Schawlow & Townes, "Infrared and Optical Masers", Phys. Rev. 1958
- [5] Kikuchi, K., "Characterization of semiconductor-laser phase noise", Opt. Express 2012

---

### 7. `src/lasers/__init__.py`

| Item | Detail |
|---|---|
| **File** | `src/lasers/__init__.py` |
| **Total lines** | 5 (was 4) |
| **Change type** | Edit |

| Lines | Change |
|---|---|
| 3, 5 | Added `from .cwlaser import CWLaser` and `'CWLaser'` to `__all__` |

---

### 8. `src/protocols/examples/bb84_ideal.py`

| Item | Detail |
|---|---|
| **File** | `src/protocols/examples/bb84_ideal.py` |
| **Total lines** | 146 (was 135) |
| **Change type** | Edit |

| Lines | Change |
|---|---|
| 30–33 | Removed `frequency=40` from `apd()` constructor call (APD no longer accepts this parameter) |
| 71–72 | Added `pout = alice_laser.power_out * 1e-3` — convert SolidStateLaser's mW output to Watts for the detector chain |
| 81–85 | `E, _ = cable(...)` → `E, pout = cable(...)` — now captures attenuated power from cable (in Watts) |
| 96–120 | **Rewrote detection section:** |
| 97–106 | New: derive PBS arm power from field ratio × `pout` (calibration-independent) |
| 108–110 | `detector.output(power=power_x, bandwidth=1e6)` — new APD signature |
| 112–120 | Old: ratio-based comparison of photon counts (`if abs(px-py)/(px+py) > 0.001`). New: differential detection with 3σ noise floor threshold (`if I_x > threshold or I_y > threshold: bob_bit = 0 if I_x > I_y else 1`) |

---

### 9. `src/protocols/examples/bb84_high_bitrate.py`

| Item | Detail |
|---|---|
| **File** | `src/protocols/examples/bb84_high_bitrate.py` |
| **Total lines** | 132 (was 135 initially, then rewrites) |
| **Change type** | Edit |

| Lines | Change |
|---|---|
| 6 | `from src.detectors.apd_v2 import APD_v2` → `from src.detectors import apd` (apd_v2 was deleted) |
| 20–23 | Removed `frequency=40` from `apd()` constructor call |
| 53–55 | Added `pout = alice_laser.power_out * 1e-3` — mW → W conversion |
| 66–69 | `E, _ = cable(...)` → `E, pout = cable(...)` — captures attenuated power |
| 82–105 | **Rewrote detection section** matching `bb84_ideal.py`: power from field ratio; current-based differential detection with 3σ threshold |

---

### 10. Generated graphs

| File | Description | Command |
|---|---|---|
| `analysis/qber_vs_distance.png` | QBER vs fiber length (10–200 km), dispersion off, corrected APD | `py run_experiments.py` (with dispersion patch) |
| `analysis/qber_vs_bitrate.png` | QBER vs detector bandwidth (10 MHz–10 GHz), 100 km fiber, corrected APD | `py run_bitrate_experiments.py` (with dispersion patch) |

---

### Summary of bugs fixed

| Bug | File | Impact | Root cause |
|---|---|---|---|
| Output power increased with fiber loss | `fiber.py:96` | All QKD results were wildly wrong | `pin / 10^(-αL/10)` instead of `pin * 10^(-αL/10)` |
| Dark count attribute `self.darkcount` undefined | `apd.py:80` | `details=True` would raise `AttributeError` | Typo: should be `self.dcr` |
| Missing imports | `fields.py:1` | `plot_field()` would raise `NameError` | `import numpy` / `pyplot` were absent |
| `arcsin(S3/S0)` returned NaN | `stokes.py:29` | RuntimeWarning whenever Stokes were computed | `S3` already normalized; `S3/S0` exceeded 1 due to FP noise |
| `frequency` used as unit-compensation hack | `apd.py:20` | Physical meaning of parameter was ambiguous | Missing `c` in `detect_photons()` forced a workaround |
| Excess noise factor applied to thermal noise | `apd.py:68` | Noise current overestimated | F should only scale shot-noise terms (Kasap Eq 4.45) |
| Photon detection units broken | `apd.py:31–33` | Photon count had units of `s/m` | Missing factor of `c` in `detect_photons()` |
| Electric field not scaled by attenuation | `fiber.py:94` | Detector received pre-attenuation field amplitude | Field was never multiplied by `sqrt(att_lin)` |

---

## 2026-06-04 — Field-as-power convention (removed separate pin/pout tracking)

### Session: 13:30–13:45 UTC+5

**Design rationale:** The E-field is the single source of truth for both polarization and optical power. Previously, `fiber.py` and `apd.py` tracked power as a separate float (`pin`/`pout`), creating a dichotomy with `channel/` which already derives power from the field. Now, power is always derived from `mean(|E|²)`.

**Field convention:** `mean(|E|²)` = optical power in **Watts**. The field is calibrated once at laser output in each script (e.g., BB84 scripts scale `E` so that `mean(|E|²) = laser.power_out * 1e-3`).

---

### 1. `src/channel/fiber.py` — cable() returns field only

| Lines | Change |
|---|---|
| 17–18 | Removed `pin` from `cable()` signature. Function is now `cable(fiber_length, E, dispersion, attenuation_factor, temperature, num_bends, pm_dispersion)` |
| 105–117 | Removed `pout = pin * att_lin`. Attenuation applied directly to field: `E = E * sqrt(att_lin)`. Convention: `mean(|E_out|²) = att_lin · mean(|E_in|²)`. Return value is just `E` (was `E, pout`). |

### 2. `src/detectors/apd.py` — output() takes E-field

| Lines | Change |
|---|---|
| 83–84 | `output(self, power, ...)` → `output(self, E, ...)`. Power derived internally: `power = mean(|E|²)`. Convention assumes mean(|E|²) = power in Watts (caller must pre-calibrate). |

### 3. `src/protocols/examples/bb84_ideal.py`

| Lines | Change |
|---|---|
| 70–72 | Added field calibration: `power_W = laser.power_out * 1e-3; E *= sqrt(power_W / mean(|E|²))` |
| 82–85 | `E, pout = cable(... pin=pout ...)` → `E = cable(...)` |
| 99–106 | Removed `field_power`, `power_x`, `power_y` computation (was deriving power from field × pout ratio). PBS arms `Ex, Ey` passed directly to detector. |
| 109–110 | `detector.output(power=power_x, ...)` → `detector.output(E=Ex, ...)` |

### 4. `src/protocols/examples/bb84_high_bitrate.py`

| Lines | Change |
|---|---|
| 53–55 | Same field calibration as bb84_ideal |
| 66–69 | `E, pout = cable(... pin=pout ...)` → `E = cable(...)` |
| 83–86 | Removed field_power/power_x/power_y computation |
| 89–90 | `detector.output(power=power_x, ...)` → `detector.output(E=Ex, ...)` |

### 5. `main.py`

| Lines | Change |
|---|---|
| 28 | Removed `pout = source.power_out` (was only used as cable input) |
| 60 | `E, _ = cable(100, E, pout, dispersion=True)` → `E = cable(100, E, dispersion=True)` |

### Files affected (this session)

| File | Lines changed | Change type |
|---|---|---|
| `src/channel/fiber.py` | 17–18, 105–117 | Edit |
| `src/detectors/apd.py` | 83–84 | Edit |
| `src/protocols/examples/bb84_ideal.py` | 70–72, 82–85, 99–110 | Edit |
| `src/protocols/examples/bb84_high_bitrate.py` | 53–55, 66–69, 83–90 | Edit |
| `main.py` | 28, 60 | Edit |

---

## 2026-06-04 — CWLaser RIN: relaxation-oscillation model (Coldren Eq 5.3.38)

### Session: 14:00–14:30 UTC+5

**`src/lasers/cwlaser.py`** — Completely rewrote `_sample_rin()` and `_generate_rin()`:

| Aspect | Old | New | Reason |
|---|---|---|---|
| **RIN model** | White noise → 1st-order LP filter at `rin_bandwidth` | White noise → shape by sqrt(S_RIN(f)) from linearized rate equations (Coldren [2] Eq 5.3.38) | The old model had no physical basis — missed the relaxation oscillation resonance entirely. |
| **Parameters** | `rin_bandwidth` (removed) | `relaxation_frequency`, `damping_rate` | f_RO and γ are physically meaningful; rin_bandwidth was arbitrary. |
| **Implementation** | Time-domain IIR (first-order) | Frequency-domain: rFFT → shape → irFFT | Avoids bilinear-transform instability at optical sampling rates. |
| **RIN PSD** | Flat LP | `S_RIN(f) ∝ (γ²+ω²) / | ω_R²-ω²+jγω|²` | Correct resonance: flat below f_RO, peak at f_RO, 1/f² roll-off above. |
| **Citations** | Kikuchi [5] only | Coldren [2] §5.3.3 (Eq 5.3.38), Petermann [5] Ch. 7 | Traceable physics. |

**Bug fixed:** `H_sq` was the raw `num/den` without normalization by the DC value, causing H_sq values ~10⁻²¹ and making the RIN output have ~0 variance. Fixed by `H_sq = H_sq_raw / H_sq_raw[0]`.

**Other improvements:**
- Added `_rin_dt_min = 1/(10·f_RO)` to decouple RIN time resolution from optical-time sampling. When `dt` is < `_rin_dt_min`, RIN is generated at the coarser rate and interpolated, preventing FFT artifacts at extreme sample rates (~10¹⁷ Hz for optical periods).
- Single-sample `get_electric_field(t, over_period=False)` now uses `_sample_rin` for the RIN value instead of `np.random.normal(0, sqrt(RIN_lin * 1e9))` which had a hardcoded 1 GHz bandwidth.

---

## 2026-06-04 — Laser characterisation script with eye diagrams

### Session: 15:00–15:15 UTC+5

**Created `analysis/laser_characterization.py`** — comprehensive verification suite for the CWLaser:

| Feature | Method | Verification |
|---|---|---|
| **Power convention** | `mean(|E|²)` vs `_power_w` | Bar chart showing error < 1% |
| **Optical linewidth** | Complex-envelope PSD via Welch, Lorentzian fit (curve_fit) | FWHM from fit vs specified Δν |
| **Phase noise** | Structure function `D_φ(τ) = ⟨[φ(t+τ)-φ(t)]²⟩` | Slope `2π·Δν` for Wiener process (Henry [1]) |
| **RIN spectrum** | Welch PSD of `δP(t)` vs Coldren Eq 5.3.38 | Resonance peak at `f_RO`, DC level at `RIN_0` |
| **Polarisation** | Stokes parameters via `compute_stokes_parameters()`; Poincaré sphere via existing `poincare()` | ψ, χ match laser settings |
| **Eye diagram (NRZ-OOK)** | PRBS → intensity modulation → RIN + phase noise → direct detection | Eye closure at increasing bitrates |

**Output files:**
| File | Content |
|---|---|
| `analysis/laser_characterization.png` | 8-panel dashboard |
| `analysis/poincare_sphere.png` | Poincaré sphere (via `stokes.py:poincare()`) |
| `analysis/eye_diagrams.png` | 5/10/25 Gbaud eyes side-by-side |

**Note on OOK modulation:** A phase modulator alone cannot produce OOK (amplitude modulation). The script uses an idealised intensity-modulator model (`E_mod = E_cw · wfm`). Real OOK requires direct laser current modulation or a Mach-Zehnder modulator.

---

## 2026-06-05 — CWLaser API redesign: `sample_field()`, MZM device model, characterisation rewrite

### Session: 08:30–10:00 UTC+5

**Three interconnected changes:**

| Change | Files | Rationale |
|---|---|---|
| `CWLaser.sample_field(dt, n_samples)` | `src/lasers/cwlaser.py` | New public method returning complex-envelope field (n_samples, 2) with all physical effects (power, phase noise, RIN, polarisation). The laser owns the physics — no more ad-hoc noise generation in characterisation scripts. |
| `get_electric_field` API cleanup | `src/lasers/cwlaser.py` | `t=0` → `dt=1e-12` (descriptive). Hardcoded `1000` → parameter `n_samples=1000`. Backward-compatible for all existing callers that use `over_period=True`. |
| MZM physical model | `src/channel/mzm.py` (new), `src/channel/__init__.py` | Push-pull MZM: `E_out = E_in·cos(π·V/V_pi)·exp(j·π·V_bias/V_pi)`. Replaces idealised `E·wfm` intensity modulator with a physically correct interferometric model (Agrawal [1] §4.2). |
| Characterisation rewritten | `analysis/laser_characterization.py` | Removed `_field_complex_envelope` / `_field_series` helpers. All plots now use `laser.sample_field()`. Eye diagrams use `MZM`. Phase noise and RIN measured from full field output (end-to-end verification). Added `matplotlib.use('Agg')` for headless operation. |

**Impact on call sites:**
- `get_electric_field(dt=..., over_period=True, n_samples=...)` — all existing callers use keyword args `over_period=True, normalize=False` which are unchanged.
- `get_electric_field(t=...)` for single-sample → now `get_electric_field(dt=...)`. No existing callers use single-sample mode, so no breakage.
- `sample_field()` is the recommended API going forward for both characterisation and BB84 scripts.
- MZM is importable as `from src.channel.mzm import MZM` or `from src.channel import MZM`.

### BB84 scripts updated to CWLaser — `src/protocols/examples/bb84_ideal.py` and `bb84_high_bitrate.py`

| Change | Rationale |
|---|---|
| `SolidStateLaser` → `CWLaser` | CWLaser provides physics-informed RIN, phase noise, and correct power scaling. |
| Removed `E *= sqrt(power_W / mean(|E|²))` calibration | CWLaser's `get_electric_field()` already has `mean(|E|²) = P_W` — no rescaling needed. |
| Removed unused `eve_laser` (bb84_ideal) | Dead code — variable was defined but never referenced in the simulation loop. |
| `dispersion=True` → `dispersion=False` (bb84_ideal) | The FFT-based dispersion function (`apply_dispersion` in fiber.py:53) computes `f = (1/100)*(1/(D·L·0.2e-9))` which produces ~10¹⁷ Hz — unphysical and corrupts the 193 THz CWLaser field. This is a pre-existing blocked issue documented in AGENTS.md. The high-bitrate script already used `dispersion=False`. |

**Verified:** 0% QBER at 10–200 km (bb84_ideal), 0% QBER at 1 MHz–5 GHz bandwidth / 100 km (bb84_high_bitrate).

---

## 2026-06-08 — Physics-based MZM rewrite; Vpi calibration fix in BB84

### Session

| Change | Files | Rationale |
|---|---|---|
| MZM rewritten as MZI + PhaseModulator | `src/channel/mzm.py` | MZM now internally uses `PhaseModulator` instances per arm. Y-branch splitter/combiner model with configurable insertion loss and extinction ratio. Supports push-pull (zero chirp) and single-drive (residual chirp) modes. `V_pi` derived from crystal parameters — no more empirical `cos(π·V/V_pi)` with a user-specified V_pi. |
| Hardcoded `Vpi = 3.757` removed from BB84 scripts | `bb84_ideal.py`, `bb84_high_bitrate.py` | Both scripts now derive V_pi from `pm_alice.Vpi` at runtime. The stale hardcoded value (3.757 V) differed by 3.2 % from the PhaseModulator's crystal-computed Vpi (3.8826 V), causing mismatched-basis QBER to drift to ~38 % instead of 50 %. With the fix: sifted QBER = 0 %, total QBER = 25 %. |
| Characterisation script updated | `analysis/laser_characterization.py` | `MZM(V_pi=5.0)` → `MZM()` (uses default PhaseModulator). Docstring updated for new switching voltage convention. |

**Verified:** `analysis/laser_characterization.py` runs clean, all three output files generated.

---

## 2026-06-09 — FFT-based chromatic dispersion, birefringence/PMD fixes in fiber.py

### Session

| Change | Files | Rationale |
|---|---|---|
| Chromatic dispersion rewritten with FFT | `src/channel/fiber.py` | Old `apply_dispersion()` used `f = (1/100)*(1/(D·L·0.2e-9))` producing ~10¹⁷ Hz. Replaced with `H(Ω) = exp(-j·β₂·Ω²·L/2)` via `np.fft.fftfreq` (Agrawal [6] §2.4). Applied to both Ex/Ey (CD is isotropic). Verified against Gaussian pulse: broadening ratio error < 0.06 % at 0.5–2.0× LD. |
| Unit fix in GVD calculation | `src/channel/fiber.py` | D stored as `17e-12` (ps/(nm·km)) but β₂ formula needs s/m². Added conversion `D_SI = D × 1e-6`. Previously β₂ was 6 orders too small. |
| Birefringence Jones matrix fixed | `src/channel/fiber.py` | Removed spurious `del_T = pmd_sd²` factor that zeroed out the beat-length phase (~10⁻¹⁶ rad instead of the correct ~10⁶ rad). Now uses `Δβ = 4π·Δn/λ`, Jones = `diag(exp(j·Δβ·L/2), 1)`, preserving the L/2 beat-length convention (SM vs PM fibre discrimination). |
| PMD rewritten (frequency-domain) | `src/channel/fiber.py` | Old `apply_pmd()` added random per-sample phase (phase noise, not PMD). Replaced with frequency-domain DGD: Maxwellian-distributed Δτ, Jones matrix `diag(exp(∓j·ω·Δτ/2), exp(±j·ω·Δτ/2))`. |
| New parameters | `src/channel/fiber.py` | Added `dt` (required for dispersion), `wavelength` (no longer hardcoded). Default `attenuation_factor` changed from 0.25 → 0.182 dB/km (SMF-28 at 1550 nm). |
| Updated `main.py` | `main.py` | Added `dt=1e-12` to `cable()` call. |

**Verified:** `analysis/laser_characterization.py`, `bb84_ideal.py` (0 % sifted QBER), `bb84_high_bitrate.py` (0 % sifted QBER). Gaussian pulse broadening matches Agrawal theory.
