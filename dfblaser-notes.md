# DFB laser work in progress

This branch holds the split step DFB laser model and the driver that
supplies its injection current. Neither is finished, and neither is wired
into the rest of the project yet. These notes record what state the code
is in, so that picking it up later does not mean working it out from
scratch.

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

`tests/test_dfblaser.py` is written against a different API than the one
that exists. See the last section.

## What was done

The file was renamed from `dfblaser_v2.py` to `dfblaser.py`. That is the
only change made to the code so far.

## What is broken right now

`laser_driver.py` still imports the old module name. Line 34 reads
`from .dfblaser_v2 import Laser, SimResult`, and line 3 of the docstring
points at `src.lasers.dfblaser_v2.Laser`. Both need the `_v2` dropped.
Until that happens neither file imports at all, so this is the first thing
to fix. Nothing else in the repository refers to the old name.

## What needs fixing

Listed roughly in the order I would do them.

### 1. The recorded output fields add two polarisations together

In `simulate()` the recorded envelope is built like this:

    e_right[rec] = self.Fx[-1] + self.Fy[-1]

`Fx` and `Fy` are the x and y polarised components. They are orthogonal,
so adding them into a single complex number does not give you a field, and
the square of that number is not the power. The line above it gets power
right, using the sum of the two squared magnitudes, so the two recorded
quantities disagree with each other. The `SimResult` docstring says that
the mean of the squared magnitude is in watts, and that only holds when
one polarisation is zero. Both are driven by independent spontaneous
noise, so in practice neither is.

The fix is to keep them apart. Make `E_right` and `E_left` arrays of shape
`(n_rec, 2)` holding the x and y components separately. That also matches
how every other component in this project passes fields around. `CWLaser`,
the MZM, the fibre and the AMZIs all use `(n, 2)` complex arrays, so the
laser output would then feed straight into them with no conversion step.
The `SimResult` docstring needs updating too, since it describes the
scalar version.

### 2. The gaussian drive pulse is only half a pulse

In `DriveParams.i()`:

    tau = t % self.period
    current = self.i_bias + self.i_peak * np.exp(-0.5 * (tau / self._pulse_sigma) ** 2)

`tau` runs from zero up to `period`, and the gaussian peaks at zero. So
every period starts at full amplitude and decays away, and you only ever
get the falling half. The parameter is documented as a full width at half
maximum, which implies a symmetric pulse, so what you actually get is half
the width you asked for.

The one line fix is to measure distance to the nearest period boundary
instead:

    tau = t % self.period
    tau = np.minimum(tau, self.period - tau)

Then the leading tail comes from the end of the previous period and the
pulse is symmetric. If the pulse needs to sit somewhere other than on the
period boundary, the alternative is an explicit offset parameter, but that
adds API for no obvious gain.

The trapezoidal branch does not have this problem. It already produces a
complete pulse. So at the moment the two waveform options do not behave
the same way as each other.

### 3. The default section count sits exactly on the convergence limit

With the defaults, kappa is 5000 per metre and dz is 600 microns divided
by 15 sections, which is 40 microns. That puts kappa times dz at exactly
0.2, which is the convergence criterion rather than safely inside it.

Two things are worth doing. Warn from `__init__` when kappa times dz goes
above 0.2, so that anyone who lowers the section count or lengthens the
grating finds out immediately instead of getting quietly wrong numbers.
Then raise the default section count from 15 to 20, which puts the product
at 0.15.

Note that changing the section count changes dz, and dt follows from dz
divided by the group velocity, so the default step size and every default
configuration result will move. That is a change in behaviour, not just a
guard, so it should be called out rather than slipped in quietly.

### 4. The default run time is about ten steps

dt works out at roughly 0.49 ps and `run_time` defaults to 5 ps, so you
get about ten steps. With `record_every` at 10 that leaves a single
recorded point. It is fine for checking that the code executes, but far
too short to show anything physical. Turn on behaviour and relaxation
oscillations need nanoseconds.

Raising it to 5 ns gives around ten thousand steps. Worth doing after item
5, because at the moment that is slow enough to be irritating as a
default.

### 5. Speed up the stepping loop

The sweep over sections has to stay a loop. The comment in the code
explains why and it is correct. `F` is advanced using the value that was
just computed, so doing the sections in parallel would break energy
conservation in the stopband. But there is repeated work that does not
need to be repeated.

`_coupling()` is called every step even though it only depends on kappa
and dz, both fixed at construction. Compute it once in `__init__`.

The four arrays holding the next step are allocated fresh every step.
Allocate them once and reuse them.

The x and y loops are separate but do identical work on different arrays.
Fusing them into a single loop over sections, working on a stacked pair,
roughly halves the Python level loop overhead.

That is maybe two to three times faster, with no new dependencies. Going
much beyond that means numba, or moving the recurrence into C, which is a
bigger decision and should not be bundled in with this.

## What was deliberately left alone

`src/lasers/__init__.py` still exports only `CWLaser`. Neither new module
is exported, so both are reachable only by full path, for example
`from src.lasers.dfblaser import Laser`. Worth sorting out, but better
done at the same time as deciding whether the class should be called
`Laser` or `DFBLaser`, so that the public surface only changes once.

`tests/test_dfblaser.py` is staying deferred. It was written against an
API that does not exist yet. It imports `DFBLaser`, `LaserParams` and
`_iir_lowpass` from `src.lasers.dfblaser`, and a module called
`src.lasers.drive`. What actually exists is `Laser` and `SimResult` in
`src.lasers.dfblaser`, and `src.lasers.laser_driver`. It also reads
`res.power_right`, where the result object has `P_right`.

Because of that the file fails at import. Pytest stops the whole run when
a collection fails, so a bare `pytest` at the repository root currently
collects nothing and exits with an error. Running
`pytest --ignore=tests/test_dfblaser.py` gets the rest of the suite back.
The other way to park it is a module level skip at the top of the file,
which leaves every test body untouched and lets a bare `pytest` work
again.

The finite gain bandwidth filter from section III of the paper is still
not implemented. The module docstring explains the reasoning, that at this
gain bandwidth the filter is wider than the simulation band, and that
argument still holds.
