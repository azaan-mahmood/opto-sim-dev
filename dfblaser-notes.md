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

The default section count stays at 15, which is what the paper recommends.
It works the criterion through for its own device and says that with a
coupling coefficient of 50 per cm the subsection length has to be under
about 40 microns, so 15 subsections or more are enough for a 600 micron
device. That puts kappa times dz at exactly 0.2, on the boundary the paper
accepts rather than inside it.

An earlier version of this work raised the default to 20 to sit clear of
the limit. That had no support in the paper, and it moved every default
result, because the subsection length sets the time step. It has been put
back to 15.

The guard carries a small tolerance so that the paper's own configuration
cannot trip the paper's own criterion. At 15 sections the product comes out
as 0.19999999999999998, which happens to fall under the limit, but that is
luck in the last bit rather than something to depend on. Checked: 14
sections and below warn, 15 and above stay silent.

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

## One change that moves existing results

The default run time went from 5 ps to 5 ns. The default section count is
unchanged at 15, so the time step is unchanged as well. Nothing in the
repository calls this model yet, so nothing broke, but it is worth knowing
before comparing against an older run.

## How it compares to the paper so far

Two checks have been done. Both used explicit section counts, so the
default revert above does not affect them.

The first is a section count sweep at 1, 10, 100 and 1000 sections, all run
to 3 ns at 100 mA. Once inside the valid range the answer stops moving: 10
sections gives 2.716 mW at the right facet and 100 sections gives 2.693 mW,
which is under one percent apart, and the ripple falls from 2.1 to 0.7 per
cent. One section gives 1.2 nanowatts with a 64 per cent ripple, which is
the right kind of nonsense for a run at fifteen times the coupling limit.
The 1000 section point is still running.

The second is a light against current curve at 40 sections. Threshold comes
out at 72.8 mA and the slope at 0.192 mW per mA counting both facets. The
two facets agree to within 0.07 per cent, which is expected with both anti
reflection coated, so output power is ambiguous by a factor of two unless
you say which one you mean.

The paper gives 8.75 mW at 100 mA with an AM response of 0.263 mW per mA at
100 MHz. That frequency is far below the relaxation oscillation frequency,
so the response is effectively the slope. Those two numbers together imply
a threshold of 100 minus 8.75 over 0.263, which is 66.7 mA. That ratio does
not care whether the paper means one facet or both, which makes it the
safer thing to compare against.

So the threshold agrees to within about 8 per cent, 72.8 against 66.7. That
is the useful result, because it says the gain, loss and recombination
settings are close to right. What is left is efficiency: our slope is 0.192
against the paper's 0.263, and our total power at 100 mA is 5.40 mW against
8.75 mW.

Two honest caveats on that comparison. The paper's figure is for a complex
coupled device, with both index and gain coupling, while this model has
index coupling only, so some difference is expected rather than wrong. And
at 100 mA the device is only about 27 mA above threshold, sitting on the
knee of the curve, where the power is very sensitive to anything that
shifts the threshold. Going from 100 to 150 mA nearly triples the output.
That makes 100 mA close to the worst point to compare absolute power at,
and it makes the slope and the threshold the more trustworthy comparisons.

One modelling artifact worth recording: both polarisations lase, and split
the output evenly once above threshold. A real DFB lases in one
polarisation only, because the other has far lower confinement and never
reaches threshold. This does not explain the power gap, since gain clamping
fixes the total no matter how many modes share it, but it is not physical.

## The grating reflection spectrum, and how it is derived

This is the strongest check done so far, because it lands on a number that
has a closed form. It reproduces the paper's figure 4: reflectivity and
phase against wavelength from 1548 to 1552 nm.

It matters that this test is built out of the model's own section
operator, the same `sech` and `j tanh` pair the time stepping uses, rather
than out of a textbook grating formula. If it had been built the second
way it would only have tested the formula.

### From the section operator to a transfer matrix

The model advances one section by mapping what enters the section to what
leaves it. Writing `p = exp((G - j delta) dz)`, `c00 = sech(gamma dz)` and
`c01 = j tanh(gamma dz)` with `gamma = |kappa|`:

    F(i+1) = p * ( c00 * F(i) + c01 * R(i+1) )        (1)
    R(i)   = p * ( c01 * F(i) + c00 * R(i+1) )        (2)

That is a scattering form. It takes the wave coming in from the left,
`F(i)`, and the wave coming in from the right, `R(i+1)`, and returns the
two waves leaving. It cannot be cascaded directly, because the inputs and
outputs sit at different planes. What cascades is a transfer matrix, which
relates both fields at one plane to both fields at the next.

Rearrange (2) for `R(i+1)`:

    p * c00 * R(i+1) = R(i) - p * c01 * F(i)

    R(i+1) = R(i) / (p * c00) - (c01 / c00) * F(i)    (3)

Put (3) into (1):

    F(i+1) = p * c00 * F(i) + p * c01 * [ R(i)/(p c00) - (c01/c00) F(i) ]
           = p * c00 * F(i) + (c01/c00) * R(i) - (p c01^2 / c00) * F(i)
           = p * (c00^2 - c01^2)/c00 * F(i) + (c01/c00) * R(i)

The bracket collapses, and this is the step that makes the whole thing
work:

    c00^2 - c01^2 = sech^2(x) - (j tanh(x))^2
                  = sech^2(x) + tanh^2(x)
                  = 1

because `sech^2 = 1 - tanh^2`. So

    F(i+1) = (p / c00) * F(i) + (c01 / c00) * R(i)    (4)

and (3) with (4) give the transfer matrix

            [  p/c00      c01/c00   ]
    T   =   [                       ]
            [ -c01/c00    1/(p c00) ]

Its determinant is exactly 1:

    det T = (p/c00)(1/(p c00)) - (c01/c00)(-c01/c00)
          = (1 + c01^2) / c00^2
          = (1 - tanh^2) / sech^2
          = 1

Note that `p` cancels, so the determinant stays 1 whatever the gain or
loss is. That is a useful invariant to assert if this ever gets a test.

### Getting the reflection coefficient

The grating is uniform, so every section has the same `T` and the whole
device is `T` raised to the power `N`. Shine light in from the left with
nothing coming back from the right, which means `R(N) = 0`:

    [ F(N) ]          [ 1 ]
    [      ]  =  T^N  [   ]
    [  0   ]          [ r ]

The second row gives the answer directly:

    0 = T21 + T22 * r        so        r = -T21 / T22

and the transmission, if wanted, is `t = T11 + T12 * r`. Power
reflectivity is `|r|^2` and the phase is the argument of `r`.

The detuning comes from the model's own expression, evaluated at
transparency so the carrier induced index change is zero:

    delta(lambda) = (2 pi / lambda) * n_eff_0 - pi / Lambda

`Lambda` is fixed by the design wavelength, so `delta` is zero at 1550 nm
and the stopband sits there.

### What came out

    peak power reflectivity   N = 15    0.99013
                              N = 200   0.99013
    tanh^2(kappa L)                     0.99013

`tanh^2(kappa L)` is the closed form peak reflectivity of a uniform
grating. The section operator reproduces it to five decimals. That is the
result worth having, because it says the coupling matrix is correct rather
than merely reasonable.

The stopband is centred on 1550 nm and runs roughly 1549.3 to 1550.7 nm,
which is the `|delta| < kappa` condition, since `delta = kappa` is 0.58 nm
off Bragg. Side lobes fall away either side at 0.31, then 0.14, then 0.02.

The zeros between the side lobes were computed independently, from
`sqrt(delta^2 - kappa^2) L = m pi`, and the model sits between 4e-4 and
1.6e-3 at those wavelengths against a peak of 0.99. They are not exactly
zero only because the sweep samples on a 1 pm grid and the nulls are
sharp, so the true minima fall between grid points.

The unwrapped phase runs from -pi/2 to 3 pi/2. It ramps smoothly across
the stopband and jumps by pi at every zero, which is what it should do,
since the reflection changes sign passing through a null.

### Two things the spectrum says beyond the check

Fifteen sections really is enough. The N = 200 curve lies underneath the
N = 15 curve almost everywhere, including side lobe heights and null
positions. That is the paper's figure 5 conclusion seen directly, and it
supports putting the default back to 15.

The waveguide loss is doing a lot of damage. The same grating with alpha
at 40 per cm drops from a peak of 0.99 to 0.45, and the nulls fill in
completely, because alpha times L is 2.4. Worth holding onto, because the
efficiency gap in the section above points at the same parameter. Our
slope is 0.192 mW per mA against the paper's 0.263. That is now the second
independent place alpha looks high. It should still not be tuned to close
the gap, because the paper's device is complex coupled and this model is
index coupled only, so gain coupling is a live alternative explanation.

The script is `analysis/validation/validate_dfb_reflection.py` and it
writes to `analysis/val_dfb/`.

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
