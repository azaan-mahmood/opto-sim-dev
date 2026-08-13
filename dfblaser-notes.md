# DFB laser model and driver

This branch holds the split step DFB laser model and the driver that
supplies its injection current. The code now imports and runs. It is still
not wired into the rest of the project, and its test is still deferred.

## What is here

`src/lasers/dfblaser.py` is a time domain split step model of a DFB laser
diode, following Kim, Chung and Lee, IEEE JQE 36(7), pages 787 to 794,
2000. The paper works in CGS units and this is converted to SI. The
docstring maps each block of code back to the numbered equation it comes
from, so checking it against the paper is easy.

It tracks four field arrays, forward and backward for each of two
polarisations, plus a carrier density for each section. Every step applies
the coupling matrix, then a per section gain and detuning factor, then
updates the carriers. Spontaneous emission enters as seeded complex
Gaussian noise.

`src/lasers/laser_driver.py` holds the current waveform. `DriveParams` is
a validated dataclass covering CW and gain switched pulsed operation, with
gaussian or trapezoidal pulse shapes. `LaserDriver` puts a device and a
waveform together and runs them.

## What was fixed

### The rename was finished

The file went from `dfblaser_v2.py` to `dfblaser.py`, but `laser_driver.py`
was still importing the old name, so neither file imported at all. Both the
import on line 34 and the docstring reference on line 3 now point at the
new name. Nothing else in the repository referred to the old one.

### Recorded output fields no longer add two polarisations together

`simulate()` used to record `Fx[-1] + Fy[-1]` as the output envelope. Those
are the x and y polarised components, which are orthogonal, so adding them
into one complex number does not give a field and its squared magnitude is
not a power. It carries a cross term that is interference between two
components that cannot interfere. The recorded power on the line above used
the sum of squared magnitudes, so the two disagreed. Measured on a short
run, the worst disagreement was 99.5 per cent.

`E_right` and `E_left` are now `(n_rec, 2)` arrays holding the two
components separately, so `sum(abs(E)**2, axis=1)` reproduces the matching
power column exactly. That also matches the shape every other component in
this project uses, so the output can feed CWLaser, the MZM, the fibre and
the AMZIs without a conversion step.

### The gaussian drive pulse is a whole pulse

`tau` was measured as `t % period`, so the gaussian peaked at the start of
each period with nothing before it and only the falling half survived. The
effective width was half the requested FWHM, measured at 50 ps for a
`width` of 100 ps.

`tau` is now the distance to the nearest period boundary, so the leading
half comes from the tail of the preceding period. Measured FWHM is now
100.00 ps for a `width` of 100 ps. The trapezoidal branch is one sided by
construction and deliberately does not get the wrap.

### The convergence criterion is checked

`__init__` now warns when kappa times dz goes above 0.2, which is the limit
from figure 5 of the paper. It warns rather than raises, because going over
degrades accuracy gradually instead of failing visibly, which is exactly
why it needs saying out loud. The message says how many sections would be
needed to get back under.

The default section count went from 15 to 20. At 15 the product was exactly
0.200, sitting on the limit rather than inside it. At 20 it is 0.150.

### The default run time is long enough to show something

`run_time` went from 5 ps to 5 ns. At 5 ps you got about ten steps, and
with `record_every` at 10 that was a single recorded point. Turn on
behaviour and relaxation oscillations need nanoseconds.

### The stepping loop does less repeated work

`_coupling()` was being called on every step even though it only depends on
kappa and dz, both fixed at construction. It is now computed once in
`__init__` and read once per `simulate()` call.

The four arrays holding the next step were allocated fresh every step. They
are now allocated once and swapped with the state arrays at the end of each
step. Swapping matters. Assigning would make the state and the write target
the same array, and the next sweep would read values it had already
overwritten.

The x and y sweeps were separate loops doing identical work on different
arrays. They are now one loop, with the boundary values read into locals.

Be aware this is a smaller win than expected. It measures between 1.07 and
1.60 times faster depending on section count, best at low section counts
where the hoisted coupling call is the larger share. An earlier version of
these notes estimated two to three times, which was too optimistic. Fusing
the loops removes iteration bookkeeping but not arithmetic, so it does
little once the section count is large. A real speedup still needs numba or
moving the recurrence into C.

## Two changes that move existing results

The default section count and the default run time both changed. Section
count sets dz, and dt follows from dz over the group velocity, so anything
that relied on the old defaults will now produce different numbers. Nothing
in the repository calls this model yet, so nothing broke, but it is worth
knowing before comparing against an older run.

## What is still open

`src/lasers/__init__.py` still exports only `CWLaser`. Both new modules are
reachable only by full path, for example
`from src.lasers.dfblaser import Laser`. Worth sorting out at the same time
as deciding whether the class should be called `Laser` or `DFBLaser`, so
the public surface only changes once.

`tests/test_dfblaser.py` is still deferred. It was written against an API
that does not exist. It imports `DFBLaser`, `LaserParams` and
`_iir_lowpass` from `src.lasers.dfblaser`, and a module called
`src.lasers.drive`. What exists is `Laser` and `SimResult` in
`src.lasers.dfblaser`, and `src.lasers.laser_driver`. It also reads
`res.power_right`, where the result object has `P_right`.

The file fails at import, and pytest abandons the whole run when a
collection fails, so a bare `pytest` at the repository root still collects
nothing. Use `pytest --ignore=tests/test_dfblaser.py`, which passes 332
tests. The other way to park it is a module level skip at the top of the
file, which leaves the test bodies untouched and lets a bare `pytest` work
again.

The finite gain bandwidth filter from section III of the paper is still not
implemented. The module docstring explains that at this gain bandwidth the
filter is wider than the simulation band, and that argument still holds.

## How the fixes were checked

The loop changes touch no arithmetic and no random draw, so they must not
move any number. This was checked by pulling the pre-fix module straight
out of commit 69e6c68 and running it beside the current one at matched
settings. Power and field come out bit for bit identical at 15, 20 and 24
sections, which is the check that nothing got reordered.

The rest was checked directly. Power and field agree exactly now and did
not before. The warning fires at 10 sections and is silent at the default.
The gaussian pulse is symmetric and its measured FWHM matches the `width`
argument, where the old one measured half. The trapezoidal pulse still
starts on the period boundary. `LaserDriver.run` returns a `SimResult` with
finite non-negative power.
