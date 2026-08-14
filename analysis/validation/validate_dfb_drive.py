"""What the DFB device does under CW and gain-switched drive.

The device is the same in both; only the injection-current waveform
changes.  Nothing measured here is an input to the model -- relative
intensity noise, chirp, the pulse-to-pulse phase and the pulse shape all
come out of the spontaneous-emission seeding and the carrier rate
equation.

What the script fails on, and why these
---------------------------------------
Shape is what validates a model like this, because shape comes from the
structure of the physics.  Absolute output power depends on parameters no
article states in full, so nothing here compares against an absolute
power:

1. **Intensity noise is present** and concentrated in a few-GHz band.
2. **Gain switching produces pulses** -- a large peak-to-mean contrast --
   rather than a modulated CW level.
3. **Pulse-to-pulse phase is near-random.**  Each pulse grows from
   spontaneous emission, so it carries no phase memory of the last one.
   This is the property that matters if the source is ever used for QKD,
   and it is the sharpest single discriminator between gain switching and
   CW.
4. **The chirp runs blue then red across the pulse.**  Frequency rises on
   the leading edge, crosses zero just after the peak and falls down the
   trailing edge.  That ordering is the gain-switching chirp signature.
5. Under ``--rin-scaling``, **RIN falls as roughly 1/P**.

Measuring the chirp: not by differentiating the phase
-----------------------------------------------------
An earlier version took ``np.gradient`` of the unwrapped phase at the
0.4933 ps device step and reported 297 GHz.  **That number was noise.**
The instantaneous frequency it produced alternated sign every sample and
ranged over 1123 GHz, which is above the device Nyquist of 1013 GHz -- the
giveaway that it was measuring the differentiator and not the field.

What is used instead is the single-lag autocorrelation already in
``_carrier_removed``, applied over a sliding boxcar: the product
``conj(E[k]) * E[k+1]`` is weighted by ``|E[k]| |E[k+1]|``, so the
near-zero stretches between pulses -- where the phase is pure noise --
contribute almost nothing.  It gives **70 to 90 GHz** and barely moves
with the window (101.8, 90.7, 86.3 GHz at 8, 20 and 40 samples).

A claim withdrawn: the resonance does NOT fail a scaling law
------------------------------------------------------------
This script used to report that the intensity-noise resonance fails the
sqrt(I - I_th) scaling of Coldren Sec. 5.3.1.  **That overstated what was
measured and is withdrawn.**

At most currents there is no single resolvable resonance line to scale.
There is a broad cluster of comparable peaks -- at 100 mA, five peaks
between 3.46 and 5.81 GHz all within a factor of 1.5 of each other -- so
``argmax`` picks between near-equals and hops, which is what produced the
non-monotonic 4.70, 6.68, 1.98, 2.97 GHz.  A spectral centroid, which is
well defined, does not scale either (r = +0.11), so a better peak-picker
does not rescue it: the quantity the law is about does not exist in this
spectrum.  Withdrawn as unmeasurable, not as a model failure.

What replaces it is stronger, because it uses the whole band: the
integrated relative intensity noise scales as **P^-1.20 (r = -0.89)**
against the textbook 1/P, over an order of magnitude in power, with
nothing fitted.  See ``--rin-scaling``.

Two measured facts drive the run parameters
-------------------------------------------
The device takes about 30 ns to settle from its zero-field initial state,
so every window here starts after a discarded settle.  And the recorded
envelope carries a large mean phase ramp -- the lasing mode's offset from
the Bragg reference the device integrates against, in the hundreds of GHz
and moving with current -- which is removed before any phase, chirp or
spectrum is taken.

References
----------
[1] Kim, Chung & Lee, "An Efficient Split-Step Time-Domain Dynamic
    Modeling of DFB/DBR Laser Diodes", IEEE J. Quantum Electron. 36(7),
    787-794 (2000).
[2] Coldren, Corzine & Mashanovitch, "Diode Lasers and Photonic
    Integrated Circuits", 2nd ed., Wiley, 2012.  Sec. 5.3: relaxation
    oscillation and RIN.
[3] Agrawal, "Fiber-Optic Communication Systems", 5th ed., Wiley, 2021,
    Sec. 3.4 (gain switching) and Sec. 3.5 (chirp).
[4] Koch & Bowers, "Nature of wavelength chirping in directly modulated
    semiconductor lasers", Electron. Lett. 20(25), 1038-1040 (1984).
"""
import argparse
import os
import sys

import numpy as np
from scipy import signal

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from src.lasers import CWLaser, DFBLaser, DriveParams, LaserDriver

OUT_DIR = os.path.join(os.path.dirname(__file__), '..', 'val_dfb')

SETTLE = 40e-9          # discarded turn-on, measured ~30 ns
N_SECTIONS = 15         # the paper's recommendation at 600 um
SEED = 11

# Noise band.  Wide enough not to prejudge the answer, narrow enough to
# exclude the DC drift and the mode structure hundreds of GHz out.
F_LO, F_HI = 1e9, 25e9

# Welch segment length, in device steps.  Plain band-averaging of one
# periodogram is not usable here: with a 25 ns window, averaging 200 bins
# gives 8 GHz per band, so every current returned the same "resonance"
# because it was reporting the first band centre.  16384 device steps is
# 8.1 ns, i.e. 0.12 GHz resolution.
NPERSEG = 16384

# Sliding window for the instantaneous-frequency estimator, in device
# steps.  20 steps is 9.9 ps.  The chirp is stable against this (101.8,
# 90.7, 86.3 GHz at 8, 20, 40); the alpha fit below is not, which is why
# one is a pass/fail check and the other is only reported.
CHIRP_WINDOW = 20

# CWLaser builds RIN on an internal grid no finer than 1/(10*f_RO) = 20 ps
# and interpolates below that, so comparing on a finer grid measures the
# interpolation rather than the model.  See vs_cwlaser().
CW_COMPARE_DT = 20e-12


def _carrier_removed(e, t, dt_dev):
    """Envelope with the lasing-mode offset divided out. Returns (e, f_offset)."""
    z = np.sum(np.conj(e[:-1]) * e[1:])
    if np.abs(z) == 0.0:
        return e, 0.0
    f0 = float(np.angle(z) / dt_dev / (2.0 * np.pi))
    return e * np.exp(-2j * np.pi * f0 * t), f0


def instantaneous_frequency(e, dt, window=CHIRP_WINDOW):
    """Instantaneous frequency (Hz), amplitude-weighted over ``window`` steps.

    Do not differentiate the unwrapped phase instead.  At the device step
    that is dominated by the differentiator, giving a signal that
    alternates sign every sample and exceeds the Nyquist frequency.  The
    lag-one product used here is weighted by the field magnitude, so the
    near-zero gaps between pulses carry almost no weight.
    """
    g = np.conj(e[:-1]) * e[1:]
    kern = np.ones(window) / window
    nu = np.angle(np.convolve(g, kern, mode='same')) / dt / (2.0 * np.pi)
    return np.append(nu, nu[-1])          # keep len(e)


def cw_run(current, t_window, record_every=1):
    """Settled CW window. Returns (t, P_total, e_carrier_removed, f_offset, dt).

    Power is the sum of both facets.  Above about 110 mA the device picks
    one facet over the other and which one wins depends on the noise seed,
    so a single facet is not a usable measure of output there; the total
    is.
    """
    las = DFBLaser(n_sections=N_SECTIONS, seed=SEED)
    drv = LaserDriver(las, DriveParams(mode='cw', i_bias=current), seed=SEED)
    res = drv.run(t_end=SETTLE + t_window, record_every=record_every)
    m = res.t >= SETTLE
    t = res.t[m] - SETTLE
    dt = las.dt * record_every
    e, f0 = _carrier_removed(res.E_right[m, 0], t, dt)
    return t, res.P_right[m] + res.P_left[m], e, f0, dt


def rin_spectrum(P, dt):
    """RIN spectral density. Returns (f, S) with S in 1/Hz.

    ``S`` is the PSD of the *relative* intensity fluctuation
    ``(P - Pbar)/Pbar``, which is the standard definition
    ``RIN(f) = S_dP(f) / Pbar**2``, quoted in dB/Hz.
    """
    x = (P - P.mean()) / P.mean()
    nper = min(NPERSEG, len(x))
    return signal.welch(x, fs=1.0 / dt, nperseg=nper, noverlap=nper // 2)


def rin_summary(f, S, f_lo=F_LO, f_hi=F_HI):
    """Mean RIN in dB/Hz over the band, and the integrated variance."""
    b = (f >= f_lo) & (f <= f_hi)
    if not b.any() or not np.any(S[b] > 0):
        return np.nan, 0.0
    df = f[1] - f[0]
    return float(10.0 * np.log10(np.mean(S[b]))), float(np.sum(S[b]) * df)


def band_peaks(f, S, n=5):
    """Strongest peaks in the band, as (frequency, height above band median).

    Reported instead of a single f_RO because at most currents this is a
    band and not a line: five peaks within a factor of 1.5 of each other
    means an argmax picks between near-equals and hops with current.
    """
    b = (f >= F_LO) & (f <= F_HI)
    if not b.any():
        return []
    fb, Sb = f[b], S[b]
    med = float(np.median(Sb))
    if med <= 0:
        return []
    idx = [j for j in range(1, len(Sb) - 1) if Sb[j] > Sb[j - 1] and Sb[j] > Sb[j + 1]]
    idx = sorted(sorted(idx, key=lambda j: -Sb[j])[:n], key=lambda j: fb[j])
    return [(float(fb[j]), float(Sb[j] / med)) for j in idx]


def phase_variance_by_decade(e, dt, lags=(1, 2, 4, 8, 16, 32, 64, 128, 256, 512)):
    """Phase-variance slope against lag, decade by decade.

    Reported as a set rather than one number because the DFB's phase has
    three regimes -- diffusion below ~16 ps, frequency wander to ~1 ns,
    then saturation -- and a single straight-line fit across them measures
    nothing.  A linewidth is only defined where the slope is 1.
    """
    ph = np.unwrap(np.angle(e))
    lags = np.array([L for L in lags if L < len(ph) // 2])
    var = np.array([np.var(ph[L:] - ph[:-L]) for L in lags])
    slopes = [(lags[i] * dt, float(np.log(var[i] / var[i - 1]) / np.log(lags[i] / lags[i - 1])))
              for i in range(1, len(lags))]
    return lags, var, slopes


def straight_line_linewidth(e, dt, lags=(1, 2, 4, 8, 16, 32, 64, 128, 256, 512)):
    """Linewidth from a straight-line fit of variance against lag, in Hz.

    Only meaningful where the variance really is linear in lag.  Used on
    CWLaser, where it is, as the control for the estimator.
    """
    lags_, var, _ = phase_variance_by_decade(e, dt, lags)
    return float(np.polyfit(lags_ * dt, var, 1)[0] / (2.0 * np.pi))


def gain_switched_run(waveform, i_bias, i_peak, period, width, t_window,
                      t_rise=20e-12):
    """Settled gain-switched window. Returns a dict of arrays and the drive."""
    las = DFBLaser(n_sections=N_SECTIONS, seed=SEED)
    d = DriveParams(mode='gain_switched', waveform=waveform, i_bias=i_bias,
                    i_peak=i_peak, period=period, width=width, t_rise=t_rise)
    res = LaserDriver(las, d, seed=SEED).run(t_end=SETTLE + t_window, record_every=1)
    m = res.t >= SETTLE
    t = res.t[m] - SETTLE
    e, f0 = _carrier_removed(res.E_right[m, 0], t, las.dt)
    return dict(t=t, P=res.P_right[m], i=res.i[m], e=e, f0=f0, dt=las.dt, drive=d)


def pulse_statistics(r):
    """Per-pulse peak, phase and chirp. Returns a dict of measurements."""
    t, P, e, dt, period = r['t'], r['P'], r['e'], r['dt'], r['drive'].period
    nu = instantaneous_frequency(e, dt)
    peaks, phases, idx, chirps, lead, trail = [], [], [], [], [], []
    for k in range(1, int(t[-1] // period)):
        w = (t >= k * period) & (t < (k + 1) * period)
        if w.sum() < 10:
            continue
        i = int(np.argmax(np.where(w, P, -1.0)))
        if P[i] <= 0:
            continue
        a, b = max(0, i - 80), min(len(P), i + 80)
        seg, nn, loc = P[a:b], nu[a:b], i - a
        s = seg > 0.10 * seg.max()
        if s.sum() < 6:
            continue
        peaks.append(P[i])
        phases.append(np.angle(e[i]))
        idx.append(i)
        chirps.append(nn[s].max() - nn[s].min())
        if s[:loc].any():
            lead.append(np.mean(nn[:loc][s[:loc]]))
        if s[loc:].any():
            trail.append(np.mean(nn[loc:][s[loc:]]))
    phases = np.array(phases)
    n = min(len(lead), len(trail))
    blue_then_red = float(np.mean(np.array(lead[:n]) > np.array(trail[:n]))) if n else 0.0
    conc = float(np.abs(np.mean(np.exp(1j * phases)))) if len(phases) else 1.0
    # Rayleigh test for circular uniformity.  A bare threshold on the
    # concentration is the wrong check: for N uniformly random phases
    # |<exp(i phi)>| falls as about 0.886/sqrt(N), so a fixed cut is far too
    # tight on a short window and too loose on a long one.  Z = N * R^2 is
    # distributed as Exp(1) under uniformity whatever N is, and goes to N
    # when the phase locks.
    rayleigh_z = float(len(phases) * conc ** 2) if len(phases) else np.inf
    return dict(
        n=len(peaks),
        peak=float(np.max(peaks)) if peaks else 0.0,
        spread=float(100 * np.std(peaks) / np.mean(peaks)) if peaks else np.nan,
        conc=conc,
        rayleigh_z=rayleigh_z,
        chirp=float(np.mean(chirps)) if chirps else 0.0,
        lead=float(np.mean(lead)) if lead else np.nan,
        trail=float(np.mean(trail)) if trail else np.nan,
        blue_then_red=blue_then_red,
        idx=np.array(idx, dtype=int),
    )


def optical_fwhm_and_delay(r):
    """Optical FWHM and the current-to-optical delay, both in seconds."""
    t, P, period = r['t'], r['P'], r['drive'].period
    k = max(1, int(t[-1] // period) // 2)
    w = (t >= k * period) & (t < (k + 1) * period)
    tt, PP = t[w] - k * period, P[w]
    if PP.max() <= 0:
        return np.nan, np.nan
    above = np.where(PP >= PP.max() / 2)[0]
    fwhm = float(tt[above[-1]] - tt[above[0]]) if len(above) > 1 else np.nan
    # Where the drive peaks inside a period, evaluated on its own function
    # rather than on the recorded samples: the gaussian peaks on the period
    # boundary and the trapezoid at the end of its rising edge, so a bare
    # argmax over one window would compare different reference points.
    tg = np.linspace(0, period, 4001, endpoint=False)
    drive_peak = float(tg[int(np.argmax(r['drive'].i(tg)))])
    return fwhm, float(tt[int(np.argmax(PP))] - drive_peak)


def alpha_from_chirp(r, window=CHIRP_WINDOW):
    """Linewidth enhancement factor implied by the chirp. Returns (alpha, R^2).

    From the standard relation (Koch & Bowers [4], Agrawal [3] Sec. 3.5),
    ``dnu = (alpha/4pi) * [dln(P)/dt + kappa*P]``, so regressing the
    instantaneous frequency on ``dln(P)/dt`` returns ``alpha/(4 pi)``.

    Reported only, never a pass/fail check.  The answer drifts with the
    smoothing window -- 3.37, 4.04, 4.90, 5.68 at 10, 20, 40, 80 steps
    against a model value of 5.0 -- so it passes through the right answer
    without settling there, and picking the window that lands on it would
    be fitting.
    """
    P, e, dt = r['P'], r['e'], r['dt']
    kern = np.ones(window) / window
    nu = instantaneous_frequency(e, dt, window)
    Ps = np.convolve(P, kern, mode='same')
    dlnP = np.gradient(np.log(np.maximum(Ps, 1e-30)), dt)
    sel = Ps > 0.05 * Ps.max()
    if sel.sum() < 50:
        return np.nan, np.nan
    A = np.column_stack([dlnP[sel], np.ones(sel.sum())])
    c, *_ = np.linalg.lstsq(A, nu[sel], rcond=None)
    r2 = 1.0 - np.var(nu[sel] - A @ c) / np.var(nu[sel])
    return float(c[0] * 4 * np.pi), float(r2)


# --- sections ---------------------------------------------------------------

def cw_section(current, window, failures):
    t, P, e, f0, dt = cw_run(current, window)
    f, S = rin_spectrum(P, dt)
    rin_db, var = rin_summary(f, S)
    peaks = band_peaks(f, S)

    print(f"\n  CW at {current * 1e3:.0f} mA over {window * 1e9:.0f} ns")
    print(f"    steady power (both) : {P.mean() * 1e3:.4f} mW")
    print(f"    rms intensity ripple: {100 * P.std() / P.mean():.2f} %")
    print(f"    lasing mode offset  : {f0 / 1e9:+.1f} GHz from the Bragg reference")
    print(f"    RIN, {F_LO / 1e9:.0f}-{F_HI / 1e9:.0f} GHz  : {rin_db:.1f} dB/Hz "
          f"(integrated variance {var:.3e})")

    if not peaks:
        print("    noise band          : NOTHING FOUND")
        failures.append("no intensity noise found in the "
                        f"{F_LO / 1e9:.0f}-{F_HI / 1e9:.0f} GHz band")
    else:
        print("    noise band          : "
              + "  ".join(f"{fp / 1e9:.2f} GHz:{h:.1f}x" for fp, h in peaks))
        print("      -> a band of comparable peaks, not one resonance line. No "
              "single f_RO\n         is quoted, and no scaling law is tested on "
              "one (see the docstring).")
        if max(h for _, h in peaks) < 2.0:
            failures.append("intensity noise is flat across the band "
                            f"(strongest peak only {max(h for _, h in peaks):.1f}x "
                            "the median)")

    print("    what RIN does not cover: it is an amplitude measure. This device")
    print("      has two instabilities it cannot see -- the facet bistability")
    print("      above ~110 mA (total power is identical either way) and the")
    print("      frequency wander below. Low RIN is necessary, not sufficient.")

    _, _, slopes = phase_variance_by_decade(e, dt)
    print("    phase variance slope by lag:")
    print("      " + "  ".join(f"{lag * 1e12:.0f}ps:{s:.2f}" for lag, s in slopes[:8]))
    print("      -> 1.0 would be diffusion. Diffusion holds only below ~16 ps; "
          "beyond\n         that the mode frequency wanders, so no linewidth is "
          "quoted.")
    return t, P, f, S, peaks


def gain_switched_section(waveform, window, failures, check_chirp_shape):
    i_bias, i_peak, period, width = 0.060, 0.140, 500e-12, 100e-12
    r = gain_switched_run(waveform, i_bias, i_peak, period, width, window)
    st = pulse_statistics(r)
    fwhm, delay = optical_fwhm_and_delay(r)
    mean_p = float(np.mean(r['P']))
    contrast = st['peak'] / mean_p if mean_p > 0 else 0.0

    print(f"\n  gain switched, {waveform}: {i_bias * 1e3:.0f} -> "
          f"{(i_bias + i_peak) * 1e3:.0f} mA, {width * 1e12:.0f} ps drive "
          f"at {1e-9 / period:.1f} GHz")
    print(f"    peak power          : {st['peak'] * 1e3:.3f} mW")
    print(f"    mean power          : {mean_p * 1e3:.3f} mW")
    print(f"    peak / mean         : {contrast:.1f}")
    print(f"    pulses measured     : {st['n']}   peak spread {st['spread']:.1f} %")
    print(f"    optical FWHM        : {fwhm * 1e12:.1f} ps from a "
          f"{width * 1e12:.0f} ps drive  -> the pulse comes out narrower")
    print(f"    turn-on delay       : {delay * 1e12:.1f} ps after the drive peak")
    print(f"    chirp across a pulse: {st['chirp'] / 1e9:.1f} GHz "
          f"(leading {st['lead'] / 1e9:+.1f}, trailing {st['trail'] / 1e9:+.1f})")
    print(f"    blue-then-red       : {100 * st['blue_then_red']:.0f} % of pulses")
    print(f"    phase concentration : |<exp(i phi)>| = {st['conc']:.3f}  "
          f"(0 = random per pulse, 1 = locked)")
    print(f"    Rayleigh Z = N R^2  : {st['rayleigh_z']:.2f}   "
          f"(Exp(1) under uniformity; {st['n']} would be fully locked)")

    alpha, r2 = alpha_from_chirp(r)
    las_alpha = DFBLaser(n_sections=N_SECTIONS).alpha_m
    if r2 < 0.20:
        print(f"    implied alpha       : not quoted, R^2 = {r2:.2f}  -- incidental")
        print("      -> dnu = (alpha/4pi) dln(P)/dt does not describe this pulse. "
              "The drive\n         ends before the pulse forms, so the carrier "
              "dynamics behind the\n         relation are not the ones shaping it.")
    else:
        print(f"    implied alpha       : {alpha:.2f} (model {las_alpha:.1f}), "
              f"R^2 = {r2:.2f}  -- incidental")
        print("      -> drifts with the smoothing window (3.37/4.04/4.90/5.68 at "
              "10/20/40/80\n         steps), so it is not a pass/fail check.")

    # Chirp measured at a second window, to show it is not window-sensitive.
    nu_alt = instantaneous_frequency(r['e'], r['dt'], window=40)
    alt = []
    for i in st['idx']:
        a, b = max(0, i - 80), min(len(r['P']), i + 80)
        s = r['P'][a:b] > 0.10 * r['P'][a:b].max()
        if s.sum() >= 6:
            alt.append(nu_alt[a:b][s].max() - nu_alt[a:b][s].min())
    if alt:
        print(f"    chirp at window 40  : {np.mean(alt) / 1e9:.1f} GHz "
              f"(vs {st['chirp'] / 1e9:.1f} at 20) -- stable, unlike the alpha fit")

    if contrast < 5.0:
        failures.append(f"{waveform} gain switching did not produce pulses "
                        f"(peak/mean = {contrast:.1f}, expected > 5)")
    if st['n'] < 5:
        failures.append(f"{waveform}: only {st['n']} pulses measured; window too "
                        "short to judge the phase statistics")
    elif st['rayleigh_z'] > 15.0:
        failures.append(f"{waveform}: pulse-to-pulse phase is not random "
                        f"(Rayleigh Z = {st['rayleigh_z']:.1f} on {st['n']} pulses, "
                        f"|<exp(i phi)>| = {st['conc']:.3f}); pulses are not "
                        "starting from spontaneous emission")
    if not (10e9 < st['chirp'] < 500e9):
        failures.append(f"{waveform}: chirp {st['chirp'] / 1e9:.1f} GHz is outside "
                        "10-500 GHz; check the frequency estimator before the model")
    if check_chirp_shape and st['blue_then_red'] < 0.80:
        failures.append(f"{waveform}: chirp runs blue-then-red in only "
                        f"{100 * st['blue_then_red']:.0f} % of pulses (expected >= 80)")
    return r, st


def period_sweep(failures):
    print("\n  repetition-period sweep (drive width scaled with the period, so "
          "this is\n  the device's limit and not the drive overlapping itself)")
    print("    period      rate   trough/peak   extinction   peak mW   verdict")
    periods = (1000e-12, 500e-12, 350e-12, 300e-12, 250e-12, 220e-12,
               200e-12, 150e-12, 100e-12)
    shortest_clean = None
    for per in periods:
        r = gain_switched_run('gaussian', 0.060, 0.140, per,
                              min(100e-12, per * 0.2), 10e-9,
                              t_rise=min(20e-12, per * 0.05))
        t, P = r['t'], r['P']
        tr = []
        for k in range(1, int(t[-1] // per)):
            w = (t >= k * per) & (t < (k + 1) * per)
            if w.sum() >= 8:
                s = P[w]
                tr.append(s.min() / max(s.max(), 1e-30))
        e = float(np.mean(tr)) if tr else 1.0
        db = -10 * np.log10(max(e, 1e-30))
        clean = e < 1e-2
        if clean:
            shortest_clean = per
        print(f"    {per * 1e12:5.0f}ps  {1e-9 / per:5.2f}GHz   {e:10.2e}   "
              f"{db:7.1f} dB  {P.max() * 1e3:8.3f}   "
              f"{'clean' if clean else 'MERGING'}")
    if shortest_clean is not None:
        print(f"\n    shortest period holding 20 dB extinction: "
              f"{shortest_clean * 1e12:.0f} ps ({1e-9 / shortest_clean:.1f} GHz)")
        print("    The limit is the turn-on delay: the optical pulse peaks 80 to")
        print("    150 ps after the drive and runs 40 to 60 ps FWHM plus a tail, so")
        print("    it needs roughly 200 ps of whatever period it is given. This is a")
        print("    bulk InGaAsP active region (0.2 um, Kim Table I), not a quantum")
        print("    well, so a few GHz is the expected ceiling for this device.")
    else:
        failures.append("no period in the sweep held 20 dB extinction")


def rin_scaling(failures):
    print("\n  RIN against output power")
    print("    RIN(f) = S_dP(f)/Pbar^2 in 1/Hz. Theory for a laser well above")
    print("    threshold is RIN ~ 1/P, so the integrated variance should go as P^-1.")
    print("     I(mA)   P_tot(mW)   RIN 1-25GHz   integrated variance")
    P_list, var_list = [], []
    for I in (0.10, 0.12, 0.15, 0.20, 0.25, 0.30):
        _, P, _, _, dt = cw_run(I, 60e-9)
        f, S = rin_spectrum(P, dt)
        db, var = rin_summary(f, S)
        P_list.append(P.mean() * 1e3)
        var_list.append(var)
        print(f"    {I * 1e3:5.0f} {P.mean() * 1e3:11.3f} {db:13.1f} dB/Hz "
              f"{var:20.3e}")
    x, y = np.log(np.array(P_list)), np.log(np.array(var_list))
    slope = float(np.polyfit(x, y, 1)[0])
    r = float(np.corrcoef(x, y)[0, 1])
    print(f"\n    integrated variance ~ P^{slope:.2f}   r = {r:+.4f}   "
          "(theory: -1)")
    if not (-1.7 <= slope <= -0.7):
        failures.append(f"RIN scales as P^{slope:.2f}, outside -0.7 to -1.7; "
                        "theory for a laser above threshold is 1/P")
    if abs(r) < 0.85:
        failures.append(f"RIN against power is not a power law (|r| = {abs(r):.2f})")


def vs_cwlaser(failures):
    dt, n = CW_COMPARE_DT, 4096
    print(f"\n  DFB in cw mode against CWLaser, matched power, {dt * 1e12:.0f} ps grid")
    print("    The grid is CWLaser's, not a choice: it builds RIN on an internal")
    print("    grid no finer than 1/(10*f_RO) = 20 ps and interpolates below that,")
    print("    so a finer comparison measures the interpolation. Total noise power")
    print("    is right at any dt; only its spectral distribution is not.")

    las = DFBLaser(n_sections=N_SECTIONS, seed=SEED)
    drv = LaserDriver(las, DriveParams(mode='cw', i_bias=0.120), seed=SEED)
    Ed = drv.sample_field(dt, n)
    Pd = np.sum(np.abs(Ed) ** 2, axis=1)
    dbm = 10 * np.log10(Pd.mean() * 1e3)
    cw = CWLaser(wavelength=1550e-9, power_dbm=dbm)
    Ec = cw.sample_field(dt, n)
    Pc = np.sum(np.abs(Ec) ** 2, axis=1)

    def _pol(E):
        return (np.mean(np.abs(E[:, 1]) ** 2)
                / max(np.mean(np.abs(E[:, 0]) ** 2), 1e-30))

    fd, Sd = rin_spectrum(Pd, dt)
    fc, Sc = rin_spectrum(Pc, dt)
    db_d, _ = rin_summary(fd, Sd)
    db_c, _ = rin_summary(fc, Sc)
    pk_d, pk_c = band_peaks(fd, Sd, n=1), band_peaks(fc, Sc, n=1)

    print(f"\n    {'':24s} {'DFB cw':>16s} {'CWLaser':>16s}")
    for lab, a_, b_ in (
        ("mean power (mW)", Pd.mean() * 1e3, Pc.mean() * 1e3),
        ("rms ripple (%)", 100 * Pd.std() / Pd.mean(), 100 * Pc.std() / Pc.mean()),
        ("Py / Px", _pol(Ed), _pol(Ec)),
        ("RIN 1-25GHz (dB/Hz)", db_d, db_c),
        ("strongest peak (GHz)", pk_d[0][0] / 1e9 if pk_d else np.nan,
         pk_c[0][0] / 1e9 if pk_c else np.nan),
        ("  its height (x median)", pk_d[0][1] if pk_d else np.nan,
         pk_c[0][1] if pk_c else np.nan),
    ):
        print(f"    {lab:24s} {a_:16.4f} {b_:16.4f}")
    print("      -> CWLaser is its resonance and nothing else, being built from the")
    print("         analytic PSD; the DFB carries a broadband spontaneous-emission")
    print("         floor underneath, so its peak stands lower above the median.")

    print("\n    phase variance slope by lag -- the control is CWLaser:")
    for E, lab in ((Ed, 'DFB cw '), (Ec, 'CWLaser')):
        _, _, slopes = phase_variance_by_decade(E[:, 0], dt)
        print(f"      {lab}: " + "  ".join(f"{lag * 1e12:.0f}ps:{s:.2f}"
                                           for lag, s in slopes[:7]))
    lw = straight_line_linewidth(Ec[:, 0], dt)
    print(f"\n    CWLaser configured linewidth {cw.linewidth / 1e6:.2f} MHz, "
          f"recovered {lw / 1e6:.2f} MHz")
    if not (0.3e6 < lw < 3e6):
        failures.append(f"the phase estimator recovered {lw / 1e6:.2f} MHz from a "
                        f"CWLaser configured at {cw.linewidth / 1e6:.2f} MHz; "
                        "without that control nothing it says about the DFB stands")
    print("      -> CWLaser is pure diffusion at every lag, as constructed, and the")
    print("         estimator returns its configured value. So the DFB's departure")
    print("         from slope 1 is the DFB's phase, not a broken measurement.")


# --- driver -----------------------------------------------------------------

def run(cw_current=0.120, cw_window=60e-9, gs_window=30e-9, quick=False,
        do_period=False, do_rin=False, do_cwlaser=False):
    os.makedirs(OUT_DIR, exist_ok=True)
    if quick:
        cw_window, gs_window = 25e-9, 12e-9

    failures = []
    print("=" * 70)
    print("DFB laser under CW and gain-switched drive")
    print("=" * 70)
    print(f"  device: N = {N_SECTIONS}, seed = {SEED}, "
          f"settle discarded = {SETTLE * 1e9:.0f} ns")

    t, P, f, S, peaks = cw_section(cw_current, cw_window, failures)
    rg, stg = gain_switched_section('gaussian', gs_window, failures,
                                    check_chirp_shape=True)
    rt, stt = gain_switched_section('trapezoidal', gs_window, failures,
                                    check_chirp_shape=False)
    print("\n    The trapezoidal chirp shape is reported, not checked. Its drive")
    print("    ends before the optical pulse forms (turn-on delay exceeds the")
    print("    drive), so the carriers shaping its chirp are not the ones shaping")
    print("    the gaussian's, and the blue-then-red ordering does not carry over.")

    if do_period:
        period_sweep(failures)
    if do_rin:
        rin_scaling(failures)
    if do_cwlaser:
        vs_cwlaser(failures)

    _write_csv(f, S, cw_current, cw_window, P)
    _figure(t, P, f, S, peaks, rg, stg, rt, cw_current)

    print()
    if failures:
        print("[FAIL]")
        for x in failures:
            print(f"  - {x}")
        return 1
    print("[PASS] noise present in band under CW; both drives produce pulses with")
    print("       near-random phase and a blue-then-red chirp where checked")
    return 0


def _write_csv(f, S, current, window, P):
    path = os.path.join(OUT_DIR, 'val_dfb_drive--cw_rin.csv')
    db, var = rin_summary(f, S)
    with open(path, 'w') as fh:
        fh.write("# DFB relative intensity noise under CW drive, "
                 "validate_dfb_drive.py\n")
        fh.write(f"# I={current:.4f} A N={N_SECTIONS} seed={SEED} "
                 f"settle={SETTLE:.3e} s window={window:.3e} s\n")
        fh.write(f"# mean power {P.mean():.6e} W, RIN {db:.2f} dB/Hz over "
                 f"{F_LO:.3e}-{F_HI:.3e} Hz, integrated variance {var:.6e}\n")
        fh.write("frequency_hz,rin_per_hz\n")
        for k in range(len(f)):
            fh.write(f"{f[k]:.6e},{S[k]:.6e}\n")
    print(f"\n  CSV: {path}")


def _figure(t, P, f, S, peaks, rg, stg, rt, cw_current):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("  matplotlib missing, skipping the figure")
        return

    fig, ax = plt.subplots(2, 2, figsize=(12, 8))

    ax[0, 0].plot(t * 1e9, P * 1e3, lw=0.6, color='tab:blue')
    ax[0, 0].set_xlabel('time after settle (ns)')
    ax[0, 0].set_ylabel('total power, both facets (mW)')
    ax[0, 0].set_title(f'CW at {cw_current * 1e3:.0f} mA, settled')
    ax[0, 0].grid(True, alpha=0.3)

    ax[0, 1].loglog(f[1:] / 1e9, S[1:], lw=0.8, color='tab:blue')
    for fp, h in peaks:
        ax[0, 1].axvline(fp / 1e9, color='tab:red', ls=':', lw=0.8, alpha=0.7)
    ax[0, 1].set_xlabel('frequency (GHz)')
    ax[0, 1].set_ylabel('RIN (1/Hz)')
    ax[0, 1].set_title('intensity noise: a band, not one resonance line')
    ax[0, 1].grid(True, alpha=0.3, which='both')

    # one gaussian pulse with its chirp
    if len(stg['idx']):
        i = int(stg['idx'][len(stg['idx']) // 2])
        a, b = max(0, i - 160), min(len(rg['P']), i + 160)
        tt = (rg['t'][a:b] - rg['t'][i]) * 1e12
        nu = instantaneous_frequency(rg['e'], rg['dt'])[a:b]
        ax[1, 0].plot(tt, rg['P'][a:b] * 1e3, color='tab:blue', lw=1.2,
                      label='power')
        ax2 = ax[1, 0].twinx()
        ax2.plot(tt, nu / 1e9, color='tab:red', lw=1.0, label='chirp')
        ax2.axhline(0, color='0.6', lw=0.6)
        ax2.set_ylabel('instantaneous frequency (GHz)', color='tab:red')
        ax[1, 0].set_ylabel('power (mW)', color='tab:blue')
        h1, l1 = ax[1, 0].get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax[1, 0].legend(h1 + h2, l1 + l2, fontsize=8, loc='upper right')
    ax[1, 0].set_xlabel('time from pulse peak (ps)')
    ax[1, 0].set_title('gaussian drive: blue leading edge, red trailing')
    ax[1, 0].grid(True, alpha=0.3)

    n = min(len(rt['t']), int(2e-9 / rt['dt']))
    ax3 = ax[1, 1].twinx()
    ax3.plot(rt['t'][:n] * 1e9, rt['i'][:n] * 1e3, color='0.7', lw=1.0,
             label='injection current')
    ax3.set_ylabel('current (mA)', color='0.5')
    ax[1, 1].plot(rt['t'][:n] * 1e9, rt['P'][:n] * 1e3, color='tab:blue', lw=0.9,
                  label='optical power')
    ax[1, 1].set_xlabel('time after settle (ns)')
    ax[1, 1].set_ylabel('right-facet power (mW)', color='tab:blue')
    ax[1, 1].set_title('trapezoidal drive: pulse forms after the current ends')
    ax[1, 1].grid(True, alpha=0.3)
    h1, l1 = ax[1, 1].get_legend_handles_labels()
    h2, l2 = ax3.get_legend_handles_labels()
    ax[1, 1].legend(h1 + h2, l1 + l2, fontsize=8, loc='upper right')

    fig.tight_layout()
    png = os.path.join(OUT_DIR, 'val_dfb_drive.png')
    fig.savefig(png, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  figure: {png}")


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--current', type=float, default=0.120,
                    help='CW operating current in A (default 0.120)')
    ap.add_argument('--window', type=float, default=60e-9,
                    help='measured window after the settle, in s')
    ap.add_argument('--quick', action='store_true',
                    help='short windows, for a smoke run')
    ap.add_argument('--period-sweep', action='store_true',
                    help='sweep the repetition period and report the minimum '
                         '(nine runs, about a minute)')
    ap.add_argument('--rin-scaling', action='store_true',
                    help='RIN against output power over six currents, checked '
                         'against the 1/P law (about a minute)')
    ap.add_argument('--vs-cwlaser', action='store_true',
                    help='compare the cw mode against CWLaser at matched power')
    a = ap.parse_args()
    sys.exit(run(cw_current=a.current, cw_window=a.window, quick=a.quick,
                 do_period=a.period_sweep, do_rin=a.rin_scaling,
                 do_cwlaser=a.vs_cwlaser))
