"""Replicate Gobby, Yuan & Shields (2004) — QBER vs distance.

References
----------
[1] Gobby, C., Yuan, Z. L., & Shields, A. J. (2004). Quantum key
    distribution over 122 km of standard telecom fiber. Appl. Phys.
    Lett. 84(19), 3762-3764.
"""
import argparse
import math
import numpy as np
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from src.protocols.bb84_time_bin import simulate_bb84_time_bin
# T_INT is topology-specific -- it is the polarisation-multiplexed
# interferometer's transmission.  Importing it is safe here because this
# module models only that chain (`simulate_qber` hardcodes
# interferometer='polarisation_multiplexed'); it must NOT be applied
# unconditionally in code that also builds the balanced topology.
from src.analytic.gobby_model import T_INT

# Gobby paper data points (from Fig 3)
GOBBY_DIST_KM = np.array([4.4, 65.0, 101.0, 122.0])
GOBBY_QBER = np.array([3.3, 3.3, 6.0, 8.9])  # percent

# --- Gobby's link budget --------------------------------------------------
# Every value in this block is taken from [1].  None of them is fitted.
#
# QBER measures the error/signal ratio and nothing else, so an error in
# either term propagates directly into it.  These five parameters set that
# ratio, and getting them from the paper rather than from convenient
# defaults is what makes the sweep a prediction.
LAM = 1550e-9
ALPHA_dB = 0.2         # paper's specified fibre attenuation
MU = 0.1               # photons per clock cycle leaving Alice

# Alice's reference:encoded intensity ratio.  The paper states 1.6:1,
# chosen so the *encoded* pulse carries 0.04 of the 0.1 photons leaving
# Alice -- a security requirement, not an interferometric one.
#
# This also fixes Bob's recombining coupler, which is NOT a free
# parameter: with Alice asymmetric, a balanced 50:50 Bob would cap the
# interferometric visibility at 0.9730, and the paper reports fringe
# visibility > 0.99 up to 65 km.  The mirrored value kappa_B = 1 - kappa_A
# gives V_int = 1 exactly and is what PolarizationMultiplexedAMZI uses by
# default.  It is derived from a published measurement and must never be
# tuned against QBER -- doing so would make this a fit again.
SPLIT_RATIO = 1.6      # -> encoded fraction 1/2.6 = 0.0385 (paper: 0.04)
GATE_WIDTH = 3.5e-9    # gated InGaAs APD
REP_RATE = 2e6         # 2 MHz clock
PULSE_WIDTH = 80e-12   # DFB pulse width

# ETA_BOB is Bob's *end-to-end* detection efficiency, not the bare
# detector QE.  The paper gives eta_Bob = 0.045, which folds the ~12 %
# APD efficiency together with 5 dB of loss in Bob's apparatus.  Passing
# the bare 0.10 -- as this script used to -- silently discards that 5 dB
# and inflates the signal by 2.22x.
ETA_BOB = 0.045

# P_E is the per-detector background click probability per clock cycle:
# the paper's measured *total* error probability, 8.5e-7, of which the
# detector dark count in the 3.5 ns gate accounts for 3.2e-7 and stray
# light from the 1.3 um clock laser leaking through the WDM filter
# accounts for the remaining 5.3e-7.
#
# The two are deliberately NOT split.  Their sum is what Gobby measured;
# splitting them would add a parameter this replication cannot
# independently constrain.  The larger half is stray light, so this is not
# a dark count rate and must not be called one: a detector-spec name on a
# lumped term invites someone to compare it against a datasheet and
# conclude the detector is anomalous.
P_E = 8.5e-7

# The SPAD model expresses background counts as a rate, so convert once
# here.  This is a lumped error-count term wearing the detector model's
# units -- 242.9 Hz -- and is not a physical dark count rate.
BACKGROUND_RATE_HZ = P_E / GATE_WIDTH

# Dead time is NOT from Gobby -- the paper does not state one.  13 us is
# the ID230 figure, kept so the detector model stays self-consistent.
# Flagged because it is the second of only two non-Gobby inputs here.
DEAD_TIME = 13e-6

# Afterpulsing is the second, and it is set to zero for this replication
# on the paper's own evidence, not on a measurement of ours:
#
#   * Fig. 3's dashed curve starts near zero, where afterpulsing at the
#     ID230 rate would put it at p_ap/2 = 2.5 %;
#   * the stated P_e = 8.5e-7 is dark count plus stray light, with no
#     afterpulse term in it;
#   * the closing summary enumerates three error mechanisms and
#     afterpulsing is not among them.
#
# That is physically consistent with the apparatus: at 2 MHz gating with
# a 13 us dead time the detector is off during most of the interval over
# which trapped carriers release.
#
# Carrying an afterpulse floor *and* the modulation error would double
# count -- both would be supplying the same measured 3.3 %, and the 0 km
# QBER comes out about 2.5 pp high, which is exactly the afterpulse term.
#
# Ownership: `afterpulse_prob` is a SPAD parameter (spad.py, default 0.05
# = ID230).  This constant is the replication-level override, a claim
# about Gobby's apparatus rather than about SPADs in general, and the
# datasheet value stays reachable through --afterpulse.
AFTERPULSE_PROB = 0.0

# Visibility is an OUTPUT of this link budget, not an input:
#
#     V = S / (S + 2*P_e),   S = mu * 10^(-alpha*L/10) * eta_Bob
#
# which reproduces both visibilities the paper states (>0.99 at 65 km ->
# 0.9925 predicted; 0.884 at 122 km -> 0.9058 predicted).  The error
# counts *produce* the visibility degradation, so injecting a measured
# visibility as well applies the same physics twice.  The old
# VISIBILITY = 0.934 did exactly that, and was wrong in kind rather than
# merely in value.  The decoder is therefore ideal and
# --visibility survives only as a diagnostic override.
#
# The paper separately states a DEVICE bound: "the classical interference
# visibility is better than 99.9%", measured with bright light, so it is
# the interferometer's own contrast rather than the link-budget quantity
# above.  It is used here as a CHECK, not an input.  In the
# polarisation-multiplexed topology the arm amplitudes set the device
# visibility -- balanced arms give exactly 1.0 -- which satisfies ">0.999",
# and bb84_time_bin refuses an injected visibility on that path for the
# double-counting reason above.
#
# KNOWN OMISSION, with its number: the gap between our 1.0 and their
# stated bound is real residual imperfection (VOA balance, polarisation
# extinction, residual path mismatch) worth at most (1 - 0.999)/2 = 0.050 %
# QBER, i.e. <=1.5 % of the 3.3 % floor.  Not modelled, for the same reason
# linewidth is carried but reported as negligible.
VISIBILITY = 1.0

# --- Modulation error ----------------------------------------------------
# Gobby decompose the QBER additively: a constant modulation-error floor
# plus the distance-dependent erroneous counts.  Fig. 3's arrow gives the
# floor as 3.3 %; the paper's closing summary names phase modulation as
# one of three mechanisms.  This is a STATED value, cited, not fitted.
#
# The magnitude that reproduces it in this model:
#     static offset d = arccos(1 - 2*e_mod)      = 0.3653 rad = 20.93 deg
#     Gaussian jitter s = sqrt(-2*ln(1 - 2*e_mod)) = 0.3695 rad = 21.17 deg
# They agree to ~1 %, so the choice describes the hardware rather than
# matching a number -- and the paper says which hardware:
#
#   "errors in the phase modulation, resulting from slight inaccuracies of
#    the phase modulator biases, as well as phase drift during the
#    experiment"
#
# Both named mechanisms are static or slow.  Neither is per-pulse random
# noise, so the STATIC offset is the default here and jitter is not used;
# PHASE_NOISE_RAD survives for hardware whose error genuinely is random
# shot to shot.  The offset is carried as a bias VOLTAGE, the unit a
# modulator is actually set in -- PhaseModulator converts it through its
# crystal-derived V_pi = 3.8826 V, giving 451.5 mV for 20.93 deg.
#
# On the magnitude: ~21 deg is derived from the stated 3.3 % within this
# model, and that is the correct treatment rather than a shortcoming.  The
# paper names the mechanism and reports its aggregate size; no separately
# measured value exists to recover, because 3.3 % *is* the measurement of a
# hand-biased apparatus.  It remains not a measurement of their modulator
# and must not be presented as one.
E_MOD = 0.033


def _bias_for_aggregate(e_mod, rate_rad_s, duration_s):
    """Static bias whose time-average WITH drift reproduces `e_mod`.

    The paper's 3.3 % is an aggregate of BOTH named mechanisms -- "slight
    inaccuracies of the phase modulator biases, **as well as phase drift
    during the experiment**".  So the bias is not `arccos(1 - 2*e_mod)`;
    that is the bias-only solution, and using it while ALSO applying drift
    counts the drift twice.

    Solve instead for the `d0` satisfying

        (1/T) * integral_0^T (1 - cos(d0 + rate*t))/2 dt  =  e_mod

    which has the closed form

        1/2 - [sin(d0 + rate*T) - sin(d0)] / (2*rate*T)  =  e_mod

    Every input is stated: `e_mod` from Fig. 3's arrow, `rate` from
    "less than 0.05 deg per second", `T` from "averaged over a 2-minute
    key transfer".  Nothing is fitted -- this is a decomposition of a
    stated aggregate into its two stated mechanisms, not a free parameter.

    At the Gobby values it gives 17.86 deg, drifting to 23.86 deg by the
    end of the transfer -- a ramp centred on the 20.93 deg the bias-only
    reading would have assigned.

    Falls back to the bias-only solution when there is no drift.
    """
    if rate_rad_s <= 0.0 or duration_s <= 0.0:
        return math.acos(1.0 - 2.0 * e_mod)
    d = rate_rad_s * duration_s
    lo, hi = 0.0, math.pi
    for _ in range(200):                      # bisection: mean is monotone in d0
        mid = 0.5 * (lo + hi)
        mean = 0.5 - (math.sin(mid + d) - math.sin(mid)) / (2.0 * d)
        lo, hi = (mid, hi) if mean < e_mod else (lo, mid)
    return 0.5 * (lo + hi)


# NOTE: defined after PHASE_DRIFT_RAD_S / KEY_TRANSFER_S below, which is why
# the assignment sits further down rather than here.

# DIAGNOSTIC ONLY, and not the mechanism for this replication.  Per-pulse
# drive jitter large enough to produce the stated floor would need
# 21.17 deg = 0.457 V = 11.8 % of V_pi, one to two orders above what drive
# electronics deliver.  The paper also attributes the floor to "slight
# inaccuracies of the phase modulator biases", which is a setting that is
# wrong and stays wrong, not noise that is fresh every pulse.
#
# Reachable so the alternative can be run deliberately:
#     simulate_qber(0, N, phase_error_rad=0.0, phase_noise_rad=PHASE_NOISE_RAD)
# Never a default here or in any component.
PHASE_NOISE_RAD = math.sqrt(-2.0 * math.log(1.0 - 2.0 * E_MOD))  # jitter 21.17

# --- Interferometer arm-length drift -------------------------------------
# The second mechanism the paper names, and the one they MEASURE:
#
#   "A drift in the phase of the interferometer, due to variations in the
#    relative lengths of the two arms, could contribute directly to the
#    QBER.  By casing both Alice's and Bob's setups in enclosures to
#    prevent air convection, we found the phase drift rate to be less than
#    0.05 deg per second"
#
# A stated, citable number -- not derived and not fitted.  It is arm-length
# drift, a property of the interferometer, so it lives on AsymmetricMZI and
# not on the modulator despite the paper's prose grouping both under
# "errors in the phase modulation".
#
# It only matters for long runs, since accumulated phase is rate * t:
# their 2-minute key transfer reaches 6.0 deg and contributes 0.091 %,
# while a 3e6-pulse run at 2 MHz lasts 1.5 s and contributes ~0.  RUN_TIME
# below lets a short run be evaluated at the paper's own duration.
PHASE_DRIFT_RAD_S = math.radians(0.05)      # 8.727e-4 rad/s
KEY_TRANSFER_S = 120.0                      # "averaged over a 2-minute key transfer"

# The static bias, solved so that bias + drift together reproduce the stated
# 3.3 % over the stated transfer duration.  See `_bias_for_aggregate`.
#
# Bias and drift must be solved TOGETHER.  The paper reports one aggregate
# floor and names two mechanisms contributing to it, so taking the
# bias-only solution `arccos(1 - 2*E_MOD)` and then applying drift on top
# counts the same measurement twice.  The error is invisible on short runs,
# where accumulated drift is a fraction of a degree, and grows with
# duration: over the paper's own 120 s transfer it puts the 0 km floor
# about 1 pp high.
PHASE_ERROR_RAD = _bias_for_aggregate(E_MOD, PHASE_DRIFT_RAD_S, KEY_TRANSFER_S)

# --- Laser linewidth -----------------------------------------------------
# Gobby state no linewidth, and in a path-matched scheme they do not need
# to: the S-L and L-S routes traverse the same total path, so the
# frequency-noise term cancels and linewidth couples only through the
# RESIDUAL mismatch left after the delay line and fibre stretcher.
#
# Carried anyway, because a real laser has a linewidth and a 2026
# replicator must be able to express their apparatus.  1550 nm DFB diodes
# run from several hundred kHz to ~10 MHz (current parts cite 2-3.2 MHz);
# 2004 hardware sits at the wider end.  A documented ASSUMPTION with a
# cited range -- never derived, never fitted.  Default mismatch 0 keeps it
# inert.
LINEWIDTH = 3.0e6          # Hz, mid-range DFB assumption
PATH_MISMATCH = 0.0        # s, residual after trimming (0 = perfectly trimmed)


# --- Statistical power guards --------------------------------------------
# A flat pulse budget is badly matched to this sweep: the sifted fraction
# runs from ~1.1e-3 at 0 km to ~1e-5 at 122 km, so a flat count
# over-samples the short end and starves the long end -- which is the only
# end anyone scrutinises.  The 122 km point once landed within 0.6 pp of
# Gobby's published value on **12 sifted bits** and was briefly read as
# corroboration; the Clopper-Pearson interval on 1/12 runs from 0.21 % to
# 38.5 %.  At full power the same point read 5.30 %.  That is what these
# guards exist to prevent.
#
# MIN_SIFTED sits below TARGET_SIFTED_DEFAULT on purpose.  Equal values
# make the run a knife-edge where any undershoot discards the whole sweep
# at the final write.
MIN_SIFTED = 500
TARGET_SIFTED_DEFAULT = 3000      # sigma ~ 0.5 pp at QBER ~8 %

# Both raised for the corrected link budget.  Correcting the link
# budget cut the signal by 3.68x, so every pulse count in this sweep grew
# by about the same factor:
#   - PILOT_BITS 200k gave only ~2 sifted bits at 122 km, far too noisy a
#     base to extrapolate a budget from (the retry loop recovered, but by
#     jumping straight to the ceiling and wasting the run).  2M gives ~20.
#   - CEILING 500M would have clipped the 122 km point, which measures a
#     sifted fraction near 1e-5 and therefore needs of order 3-7e8 pulses
#     for 3000 sifted bits.  A silently clipped point looks identical to
#     an honest one in the artifact.
PILOT_BITS = 2_000_000
CEILING_DEFAULT = 1_000_000_000


def signal_click_prob(dist_km):
    """Signal click probability per clock cycle,
    S = mu * T_INT * T_link * eta_Bob.

    `T_INT` is the polarisation-multiplexed interferometer transmission,
    `2/(1 + r)` = 0.769231 for the paper's 1.6:1 split.  It is imported
    from the analytic model rather than restated here so there is one
    definition, already covered by
    `test_analytic_gobby.py::test_mu_eff_is_derived_not_fitted`.

    Omitting it puts this function 41% above what the chain delivers and
    makes every quantity built on it (`model_qber`,
    `predicted_visibility`) optimistic.  It is derived from the stated
    split ratio, not fitted.

    ACCURACY AGAINST THE MONTE CARLO, measured with dead time and
    afterpulsing off:

        0 km   S_mc / S_here = 1.019 +/- 0.012   [7,056 sifted]
        40 km                  1.011 +/- 0.025   [1,664 sifted]
        65 km                  0.985 +/- 0.031   [1,025 sifted]

    Agreement to ~2%, consistent with unity across the range.

    Detector dead time is a further ~3% at 0 km, falling to nothing at
    range as click rates drop.  That is deliberately NOT carried here: a
    first-order link budget describes the optical chain, and dead time is a
    property of the detector.  The Monte Carlo has it.
    """
    return (MU * T_INT * 10.0 ** (-ALPHA_dB * np.asarray(dist_km) / 10.0)
            * ETA_BOB)


def predicted_visibility(dist_km, p_e=P_E):
    """Fringe visibility as Gobby's own closed form, V = S/(S + 2*P_e).

    This is an *output* of the link budget, not a tunable input -- see the
    VISIBILITY note above.  Reproduces both visibilities
    the paper states: >0.99 at 65 km (0.9903 here) and 0.884 at 122 km
    (0.8809 here).

    The 122 km figure is the sharper test, the paper giving a value rather
    than a bound, and it is what confirmed the `T_INT` omission above:
    before the fix this read 0.9058, +0.022 off; after, 0.8809, -0.003 off.
    """
    s = signal_click_prob(dist_km)
    return s / (s + 2.0 * p_e)


def model_qber(dist_km, p_e=P_E, visibility=VISIBILITY,
               afterpulse_prob=AFTERPULSE_PROB):
    """First-order analytic QBER for the time-bin chain, in percent.

    Every click channel is weighted into a single ratio rather than added
    as a standalone offset:

        c_sig  = mu * eta_Bob * T_link        signal click probability
        c_bkg  = 2 * P_e                      BOTH detectors see background
        c_prim = c_sig + c_bkg                primary clicks
        c_tot  = c_prim * (1 + a)             plus afterpulses

        errors = c_sig*(1-V)/2                wrong port, finite visibility
               + c_bkg/2                      background in the wrong detector
               + a*c_prim/2                   afterpulses land at random

        QBER   = errors / c_tot

    With the decoder ideal (V = 1) and afterpulsing off this collapses to
    P_e/(c_sig + 2*P_e), which is exactly e_opt = (1 - V_pred)/2 for the
    emergent visibility above -- the two forms are the same statement.

    Correctness of the limits: as c_sig -> 0 this tends to 1/2 (a
    background-dominated link is a coin flip), and as c_bkg -> 0 it tends
    to [(1-V)/2 + a/2]/(1+a), the misalignment-plus-afterpulse floor.

    ------------------------------------------------------------------
    THIS IS AN APPROXIMATION.  THE MONTE CARLO IS AUTHORITATIVE.
    ------------------------------------------------------------------
    Do NOT fit detector parameters against this function.  It is a
    first-order closed form and it does not track the Monte Carlo closely
    enough to carry a fit: a parameter tuned to make this expression match
    the data will not make the chain match it.  Every parameter here comes
    from the paper, which is what makes the sweep a prediction.

    Weighting every click channel into a single ratio is what makes it
    correct in its limits.  The form it replaced,
    `(P_dark/2)/(p_signal + P_dark) + QBER_opt`, had three defects: it
    counted only one detector's background rate, halved the background
    error term that should be whole, and bolted the misalignment on as an
    unweighted additive constant instead of a share of the clicks.

    Residual error is concentrated at the background-dominated end, where
    dead time, afterpulse chaining off background clicks, and double-click
    tie-breaks all interact in ways this closed form does not attempt to
    capture.  Treat it as a guide curve for figures, not as a predictive
    model.
    """
    c_sig = signal_click_prob(dist_km)
    c_bkg = 2.0 * p_e                  # both detectors see the background
    c_prim = c_sig + c_bkg

    errors = (c_sig * (1.0 - visibility) / 2.0
              + c_bkg / 2.0
              + afterpulse_prob * c_prim / 2.0)
    c_tot = c_prim * (1.0 + afterpulse_prob)
    return errors / c_tot * 100.0


def gobby_measured_qber(dist_km):
    """Gobby et al. (2004) measured QBER, interpolated from Fig 3 data.

    Returns the published value exactly at the four measured distances
    (4.4, 65, 101, 122 km); clamped at the endpoints elsewhere.

    Anything this returns away from those four distances is a *derived*
    number.  Callers that put it in an artifact must say which is which
    -- see gobby_is_measured() and the column marker in the table writer.
    """
    return np.interp(dist_km, GOBBY_DIST_KM, GOBBY_QBER)


# Half-width, in km, within which one of our sweep distances is treated as
# comparable to one of Gobby's published ones (4 km <-> 4.4, 100 <-> 101).
GOBBY_MATCH_TOL_KM = 1.5


def gobby_nearest_measured(dist_km):
    """(published_km, published_qber) if this distance is close enough to
    one Gobby actually measured, else None.

    Gobby published four points.  Our sweep has nine, so five entries in
    the "Gobby et al. measured" column are np.interp output -- and below
    4.4 km np.interp does not even interpolate, it clamps to the endpoint,
    so the flat 3.3 % at 0/10/20/40 km is one published number repeated
    four times rather than four data points.

    Reporting those as "measured" is the same species of defect as
    a derived column presented as experimental data.  The values do come
    from the paper, but a reader comparing nine rows would believe there
    are nine measurements to disagree with, and there are four.

    Returning the *published* pair rather than a bare flag matters: at
    100 km np.interp gives 5.93 while what Gobby actually measured is
    6.0 at 101 km.  Marking a row "measured" while printing the
    interpolated value would be worse than not marking it at all.
    """
    i = int(np.argmin(np.abs(GOBBY_DIST_KM - dist_km)))
    if abs(GOBBY_DIST_KM[i] - dist_km) <= GOBBY_MATCH_TOL_KM:
        return float(GOBBY_DIST_KM[i]), float(GOBBY_QBER[i])
    return None


def gobby_is_measured(dist_km):
    """Boolean form of gobby_nearest_measured(), for the CSV flag column."""
    return gobby_nearest_measured(dist_km) is not None


# --- Cross-checks against the paper's prose ------------------------------
# Both bounds below come from the text rather than Fig. 3, so they test the
# budget against something it was not built from.  They are REPORTED, never
# used to adjust a parameter: tuning P_E until the 65 km check passed would
# make it fitted.

# "the contributions due to detector dark counts and stray light are less
# than 0.4%" for fibre lengths up to 65 km.
#
# READING: the sentence names BOTH terms, so the bound
# covers P_E in full (3.2e-7 dark + 5.3e-7 stray).  Reading it as the dark
# term alone would let us pass comfortably -- 0.184% at 65 km, crossing only
# at 82 km -- but the text does not support that, and taking it would be
# choosing the interpretation that flatters the model.
ERRONEOUS_BOUND_PCT = 0.4    # "less than 0.4%" up to 65 km
ERRONEOUS_BOUND_KM = 65.0
DEVICE_VISIBILITY_BOUND = 0.999   # "better than 99.9%", classical/bright-light


def erroneous_count_share(dist_km, p_e=P_E):
    """Dark-count + stray-light contribution to QBER alone, in percent.

    This is the dashed line of Fig. 3: the error counts weighted against
    the signal, with the modulation error left out.
    """
    s = signal_click_prob(dist_km)
    return 100.0 * p_e / (s + 2.0 * p_e)


def erroneous_bound_check(p_e=P_E, s_65=None):
    """Check the budget against the stated <0.4% erroneous-count bound.

    `s_65` overrides the signal click probability at 65 km.  Pass the
    value the Monte Carlo actually delivers -- `signal_click_prob` is the
    first-order form and runs optimistic (see the note in
    print_paper_cross_checks), and the verdict depends on which is used.

    Returns (value_pct, passes, implied_p_e).  The implied P_E asks what
    the bound would allow at fixed signal, so the size of any miss is
    legible -- not so it can be substituted in.
    """
    s = signal_click_prob(ERRONEOUS_BOUND_KM) if s_65 is None else s_65
    v = 100.0 * p_e / (s + 2.0 * p_e)
    f = ERRONEOUS_BOUND_PCT / 100.0
    implied = f * s / (1.0 - 2.0 * f)
    return v, v <= ERRONEOUS_BOUND_PCT, implied


def print_paper_cross_checks(p_e=P_E, visibility=VISIBILITY, s_65=None):
    """Report both prose bounds, pass or fail, without changing anything."""
    print("  Cross-checks against the paper's text (not Fig. 3):")
    v, ok, implied = erroneous_bound_check(p_e, s_65)
    src = "first-order" if s_65 is None else "MC-measured"
    print(f"    [{'OK  ' if ok else 'MISS'}] dark+stray at "
          f"{ERRONEOUS_BOUND_KM:g} km = {v:.3f}% vs stated "
          f"<{ERRONEOUS_BOUND_PCT}%   [{src} S]")
    if not ok:
        # Where the curve actually crosses, so the miss is legible as an
        # endpoint effect rather than a failure across the whole range.
        lo, hi = 0.0, ERRONEOUS_BOUND_KM
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            share = 100.0 * p_e / (signal_click_prob(mid) + 2.0 * p_e)
            lo, hi = (mid, hi) if share < ERRONEOUS_BOUND_PCT else (lo, mid)
        print(f"           over by {v / ERRONEOUS_BOUND_PCT - 1.0:+.0%} at the "
              f"endpoint; holds to {lo:.1f} km.")
        # Propagate the paper's OWN 122 km visibility back to 65 km using
        # only its own alpha and P_e -- no part of this model enters -- and
        # see where it lands against its own bound.
        v122 = 0.884
        s122_req = v122 * 2.0 * p_e / (1.0 - v122)
        s65_req = s122_req * 10.0 ** (ALPHA_dB * (122.0 - 65.0) / 10.0)
        share_paper = 100.0 * p_e / (s65_req + 2.0 * p_e)
        print(f"           CAUSE: the paper's own V(122) = {v122} "
              f"implies S(122) = {s122_req:.3e},")
        print(f"           hence S(65) = {s65_req:.3e} and a share of "
              f"{share_paper:.3f}% -- above its own bound.")
        print(f"           We read {v:.3f}%, agreeing with that to "
              f"{abs(v / share_paper - 1) * 100:.0f}%.  So we reproduce the "
              f"paper's")
        print(f"           VISIBILITY and inherit its tension with the paper's "
              f"BOUND; the excess is")
        print(f"           not a defect in this budget.  Note '<0.4%' is a "
              f"BOUND, not a measurement --")
        print(f"           the two simply cannot both be tight, and we match "
              f"the measured one.")
        print("           NOT tuned -- moving P_E would make it a fitted "
              "parameter.")
    if s_65 is None:
        print(f"           (first-order S agrees with the chain to ~2% "
              f"since the missing T_INT and")
        print(f"           the SPAD Poisson form are both correct, so the "
              f"number above")
        print(f"           is trustworthy.)")
    dev_ok = visibility >= DEVICE_VISIBILITY_BOUND
    print(f"    [{'OK  ' if dev_ok else 'MISS'}] device visibility = "
          f"{visibility:.4f} vs stated >{DEVICE_VISIBILITY_BOUND}"
          f"   (balanced arms give 1.0)")


def simulate_qber(dist_km, num_bits, seed=42, verbose=False,
                  visibility=VISIBILITY, p_e=P_E,
                  afterpulse_prob=AFTERPULSE_PROB,
                  phase_error_rad=PHASE_ERROR_RAD, phase_noise_rad=0.0,
                  linewidth=LINEWIDTH, path_mismatch=PATH_MISMATCH,
                  phase_drift_rad_s=PHASE_DRIFT_RAD_S,
                  run_duration=KEY_TRANSFER_S):
    """Run Monte Carlo at a given distance.

    Every optical and link parameter comes from the module constants,
    which come from [1].  `p_e` and `afterpulse_prob` are exposed only so
    the diagnostic control runs are reachable (--p-e, --afterpulse), not
    because they are meant to be tuned -- see the P_E and AFTERPULSE_PROB
    notes above.

    `alpha_dB` is passed explicitly.  It previously was not, so the Monte
    Carlo silently used simulate_bb84_time_bin's own 0.182 default while
    the analytic curve used this module's ALPHA_dB.  The two agreed only
    because both happened to be 0.182; correcting one without the other
    would have desynced the curve from the data drawn over it.
    """
    results = simulate_bb84_time_bin(
        num_bits=num_bits,
        fiber_length=dist_km, alpha_dB=ALPHA_dB,
        mu=MU, wavelength=LAM,
        repetition_rate=REP_RATE, pulse_width=PULSE_WIDTH,
        spad_eta=ETA_BOB,
        # Lumped error-count term (detector dark + clock-laser stray
        # light) wearing the SPAD model's units.  Not a dark count rate.
        dark_count_rate=p_e / GATE_WIDTH,
        afterpulse_prob=afterpulse_prob, dead_time=DEAD_TIME,
        gate_width=GATE_WIDTH,
        visibility=visibility,      # ideal by default; V is an output
        phase_error=0.0,
        # Gobby's actual apparatus: a polarising beam combiner/splitter
        # routes Alice-short -> Bob-long and Alice-long -> Bob-short, so
        # only the S-L and L-S paths exist.  A balanced 50:50 pair would
        # instead produce S-S and L-L satellite bins carrying half the
        # launched energy, which the gate discards -- a mu/2 loss that is
        # an artefact of the wrong interferometer, not a property of the
        # experiment.  Measured: balanced gates mu/2 (ratio 2.0 with dead
        # time off), polarisation-multiplexed gates the full mu (1.0).
        interferometer='polarisation_multiplexed',
        split_ratio=SPLIT_RATIO,
        phase_error_rad=phase_error_rad,
        phase_noise_rad=phase_noise_rad,
        linewidth=linewidth, path_mismatch=path_mismatch,
        # Both cited: the rate from "less than 0.05 deg per second", the
        # duration from "averaged over a 2-minute key transfer".  Neither
        # is tuned.
        #
        # `run_duration` is what makes the drift independent of the pulse
        # budget.  Without it the sweep's own statistics set the simulated
        # experiment length -- the 122 km point needs 1e9 pulses = 500 s at
        # 2 MHz, which accumulated 25 deg of drift against their 6 and put
        # the QBER at 13.52 % against a stated 8.9 %.  The pulse count at
        # long range therefore EXCEEDS what the apparatus sent in 120 s,
        # deliberately: these are more samples of the same two-minute
        # experiment, not a longer one.
        phase_drift_rad_s=phase_drift_rad_s,
        run_duration=run_duration,
        seed=seed, verbose=verbose,
    )
    return results


def qber_err_pp(qber_frac, n_sifted):
    """Binomial s.d. of a QBER estimate, in percentage points."""
    if n_sifted <= 0:
        return float('inf')
    return np.sqrt(qber_frac * (1.0 - qber_frac) / n_sifted) * 100.0


def run_to_target(dist_km, target_sifted, ceiling, seed, visibility,
                  pilot_bits=PILOT_BITS, p_e=P_E,
                  afterpulse_prob=AFTERPULSE_PROB):
    """Grow the pulse count at one distance until `target_sifted` sifted
    bits are collected, or `ceiling` pulses are spent.

    Runs a cheap pilot to measure this distance's sifted fraction, then
    scales to the required count and re-runs once from the same seed (so
    the result stays reproducible from the reported pulse count alone --
    it is a single deterministic run, not an accumulation).

    Returns (pulses_used, results_dict).
    """
    pilot_n = int(min(pilot_bits, ceiling))
    r = simulate_qber(dist_km, pilot_n, seed=seed, visibility=visibility,
                      p_e=p_e, afterpulse_prob=afterpulse_prob)

    n = pilot_n
    for _ in range(3):
        if r['n_sifted'] >= target_sifted or n >= ceiling:
            break
        if r['n_sifted'] == 0:
            n_next = ceiling    # no signal in the pilot: spend the budget
        else:
            # Re-estimate from the largest run so far. Headroom covers the
            # sampling error on the sifted count (~1/sqrt(n_sifted)); a
            # flat 15% undershoots about half the time, and an undershoot
            # that trips the write guard throws the whole sweep away.
            frac = r['n_sifted'] / n
            headroom = 1.15 + 3.0 / np.sqrt(max(r['n_sifted'], 1))
            n_next = int(np.ceil(target_sifted / frac * headroom))
        n_next = int(min(max(n_next, n + 1), ceiling))
        if n_next <= n:
            break
        n = n_next
        r = simulate_qber(dist_km, n, seed=seed, visibility=visibility,
                          p_e=p_e, afterpulse_prob=afterpulse_prob)

    return n, r


def check_statistical_power(rows, min_sifted, allow_underpowered):
    """Refuse to emit a table whose rows cannot support a claim.

    `rows` is a sequence of (distance_km, n_sifted).  Raises unless every
    row clears `min_sifted`, or the caller explicitly opts out.
    """
    weak = [(d, s) for d, s, in rows if s < min_sifted]
    if not weak:
        return
    detail = ', '.join(f"{d:g} km: {s} sifted" for d, s in weak)
    if allow_underpowered:
        print(f"\n  !! WARNING: {len(weak)} row(s) below {min_sifted} sifted "
              f"bits ({detail}).")
        print("     Emitting anyway because --allow-underpowered was given.")
        print("     This table is a smoke run. Do NOT cite it as a result.")
        return
    raise RuntimeError(
        f"Refusing to write an under-powered table: {len(weak)} row(s) have "
        f"fewer than {min_sifted} sifted bits ({detail}).\n"
        f"A row below {min_sifted} cannot resolve a Gobby-scale QBER "
        f"difference at 3 sigma, so the artifact would look like a result "
        f"without being one: at such counts a point can land within a fraction of a pp of the published value by chance.\n"
        f"Fix: re-run with --target-sifted {TARGET_SIFTED_DEFAULT} (grows the "
        f"pulse count per distance), or raise --bits.\n"
        f"To emit anyway for a smoke test, pass --allow-underpowered."
    )


def _stem(base, reduced):
    """Smoke runs write to their own files.

    A run that opts out of the statistical-power guard is by definition
    not quotable, so it must not be able to replace an artifact that is.
    Before this, `--allow-underpowered` wrote straight over
    `val_gobby_table.tex`: a
    smoke table sitting in the repository looking like a result.  The
    guard refusing to write was the only thing standing in the way, and
    the flag exists precisely to switch that off.

    `reduced` is driven by `allow_underpowered` rather than by a budget
    threshold, because the budget needed to be quotable is not a fixed
    pulse count -- it depends on distance, and `--target-sifted` is what
    expresses it.  Opting out of the check is the honest signal.

    The marker goes last, before the extension: `.gitignore` matches
    `*--quick.csv`, `*--quick.png` and `*--quick.md`, so a name like
    `val_gobby--quick--seed42.csv` would fall outside the ignore rules.
    """
    return f'{base}--quick' if reduced else base


# Gates for check_results.  Set against the published sweep, which agrees
# with Gobby's four points to between 0.05 and 0.69 pp and has a binomial
# sigma of 0.28 to 0.43 pp, so each gate below has roughly three times the
# margin of the disagreement it is meant to allow.
PAPER_TOL_PP = 2.0      # simulated vs published, at the four measured points
FLOOR_TOL_PP = 1.5      # simulated minus analytic background, vs E_MOD
RISE_SIGMA = 5.0        # 122 km must exceed 80 km by this many sigma


def check_results(dist, qber, sifted, p_e, visibility, afterpulse_prob,
                  allow_underpowered):
    """Check the sweep against the paper and against the analytic model.

    Returns 0 on pass, 1 on fail.  Three claims, and none of them is new:
    the script already computed all three quantities and printed them side
    by side without ever comparing them.

    Skipped entirely on an under-powered run.  `--quick` switches the
    statistical-power guard off and collects a few dozen sifted bits per
    point, where the binomial sigma is larger than every gate here, so
    asserting any of it would fail at random and mean nothing when it
    passed.
    """
    if allow_underpowered:
        print("\n  Checks skipped: this run opted out of the "
              "statistical-power guard,")
        print("  so its counts cannot support them. Use --full for a run "
              "that checks.")
        return 0

    failures = []
    print("\n  Checks")
    print("  " + "-" * 68)

    # 1. The four points Gobby actually published.  The other five rows of
    #    the sweep are np.interp output and cannot agree or disagree with
    #    an experiment, so they are not compared.
    print(f"  {'km':>5} {'simulated':>12} {'published':>10} {'diff':>8}")
    n_measured = 0
    for d, q, s in zip(dist, qber, sifted):
        near = gobby_nearest_measured(d)
        if near is None:
            continue
        n_measured += 1
        zg, qg = near
        diff = q - qg
        print(f"  {d:5.0f} {q:9.2f} +/-{qber_err_pp(q / 100.0, s):4.2f} "
              f"{qg:10.1f} {diff:+8.2f} pp")
        if abs(diff) > PAPER_TOL_PP:
            failures.append(
                f"at {d:g} km the simulation gives {q:.2f} % against a "
                f"published {qg:.1f} % at {zg:g} km, a difference of "
                f"{diff:+.2f} pp past the {PAPER_TOL_PP:g} pp allowed")
    if n_measured < len(GOBBY_DIST_KM):
        failures.append(
            f"only {n_measured} of the {len(GOBBY_DIST_KM)} published "
            f"distances were covered by this sweep, so the replication "
            f"claim is not being tested at full width")

    # 2. The two error sources add, and the modulation floor does not
    #    depend on distance.  `model_qber` is the background term only, so
    #    the simulation minus it should sit at E_MOD everywhere.
    floors = [q - model_qber(d, p_e, visibility, afterpulse_prob)
              for d, q in zip(dist, qber)]
    worst = max(range(len(floors)),
                key=lambda i: abs(floors[i] - 100 * E_MOD))
    print(f"  simulated minus analytic background: "
          f"{min(floors):.2f} to {max(floors):.2f} pp "
          f"against E_MOD = {100 * E_MOD:.1f} pp")
    if abs(floors[worst] - 100 * E_MOD) > FLOOR_TOL_PP:
        failures.append(
            f"at {dist[worst]:g} km the simulation sits {floors[worst]:.2f} "
            f"pp above the analytic background, against a modulation floor "
            f"of {100 * E_MOD:.1f} pp. The two sources add and the floor "
            f"does not depend on distance, so this gap should be flat")

    # 3. Background takes over at long range.  This is the shape of the
    #    paper's figure: flat at the floor, then rising as the signal falls
    #    and the background keeps its rate.
    i_lo = min(range(len(dist)), key=lambda i: abs(dist[i] - 80))
    i_hi = len(dist) - 1
    rise = qber[i_hi] - qber[i_lo]
    s_rise = math.hypot(qber_err_pp(qber[i_hi] / 100.0, sifted[i_hi]),
                        qber_err_pp(qber[i_lo] / 100.0, sifted[i_lo]))
    print(f"  rise from {dist[i_lo]:g} to {dist[i_hi]:g} km: "
          f"{rise:+.2f} pp ({rise / s_rise:.1f} sigma)")
    if rise < RISE_SIGMA * s_rise:
        failures.append(
            f"QBER rose only {rise:+.2f} pp between {dist[i_lo]:g} and "
            f"{dist[i_hi]:g} km ({rise / s_rise:.1f} sigma). The signal "
            f"falls with distance while the background does not, so the "
            f"error rate has to climb at the far end")

    print()
    if failures:
        print("[FAIL]")
        for f_ in failures:
            print(f"  - {f_}")
        return 1
    print(f"[PASS] all {n_measured} published points reproduced within "
          f"{PAPER_TOL_PP:g} pp")
    print("[PASS] the modulation floor is flat across the sweep and the "
          "two error")
    print("       sources add, as the analytic model assumes")
    print(f"[PASS] QBER climbs at long range ({rise / s_rise:.1f} sigma), "
          f"which is the shape of the paper's figure")
    return 0


def run_validation(num_bits=200000, seed=42, distances=None,
                   visibility=VISIBILITY, target_sifted=None,
                   ceiling=CEILING_DEFAULT, allow_underpowered=False,
                   min_sifted=MIN_SIFTED, p_e=P_E,
                   afterpulse_prob=AFTERPULSE_PROB):
    if distances is None:
        distances = [0, 4, 10, 20, 40, 65, 80, 100, 122]
    print("=" * 68)
    print("Gobby et al. 2004 — QBER vs Distance Validation")
    print("=" * 68)
    print("  Link budget (all from the paper — nothing fitted):")
    print(f"    alpha    = {ALPHA_dB} dB/km        mu       = {MU} photons/clock")
    print(f"    eta_Bob  = {ETA_BOB}          (incl. Bob's 5 dB apparatus loss)")
    print(f"    gate     = {GATE_WIDTH*1e9:.1f} ns          clock    = "
          f"{REP_RATE/1e6:g} MHz")
    print(f"    P_e      = {p_e:.2e}/clock  (dark 3.2e-7 + stray light 5.3e-7)")
    print(f"             = {p_e/GATE_WIDTH:.1f} Hz equivalent in the SPAD model")
    print("  Non-Gobby inputs (detector model self-consistency):")
    print(f"    afterpulse = {afterpulse_prob}         dead time = "
          f"{DEAD_TIME*1e6:g} us   [ID230]")
    print(f"    decoder V  = {visibility}"
          + ("          (ideal — V is an OUTPUT here)"
             if visibility == 1.0 else
             "        !! DIAGNOSTIC OVERRIDE — V is normally an output"))
    print_paper_cross_checks(p_e, visibility)
    if target_sifted:
        print(f"  Budget: --target-sifted {target_sifted} per point "
              f"(ceiling {ceiling:,} pulses)")
    else:
        print(f"  Budget: flat {num_bits:,} pulses per point")
    print(f"  Seed: {seed}")
    print("=" * 68)

    sim_dist = []
    sim_qber = []
    sim_sifted = []
    sim_pulses = []

    for d in distances:
        print(f"\n--- Distance: {d} km ---")
        if target_sifted:
            used, r = run_to_target(d, target_sifted, ceiling, seed,
                                    visibility, p_e=p_e,
                                    afterpulse_prob=afterpulse_prob)
        else:
            used, r = num_bits, simulate_qber(
                d, num_bits, seed=seed, verbose=True, visibility=visibility,
                p_e=p_e, afterpulse_prob=afterpulse_prob)
        sim_dist.append(d)
        sim_qber.append(r['qber'] * 100.0)
        sim_sifted.append(r['n_sifted'])
        sim_pulses.append(used)
        if r['n_sifted'] > 0:
            sd = qber_err_pp(r['qber'], r['n_sifted'])
            print(f"  QBER: {r['qber']*100:.2f} +/- {sd:.2f}%  "
                  f"({used:,} pulses, {r['n_sifted']} sifted)")
            print(f"        (Gobby measured: {gobby_measured_qber(d):.1f}%, "
                  f"this work analytic: "
                  f"{model_qber(d, p_e, visibility, afterpulse_prob):.1f}%)")
            print(f"        (emergent V = {predicted_visibility(d, p_e):.4f}, "
                  f"S = {signal_click_prob(d):.3e}/clock)")
        else:
            print(f"  No sifted bits in {used:,} pulses.")

    sim_dist = np.array(sim_dist)
    sim_qber = np.array(sim_qber)
    sim_sifted = np.array(sim_sifted)
    sim_pulses = np.array(sim_pulses)

    # Power guard: refuse to emit an artifact that looks like a result
    # but cannot support one.  Checked before any file is written.
    check_statistical_power(list(zip(sim_dist, sim_sifted)), min_sifted,
                            allow_underpowered)

    # Analytical curve
    dense_dist = np.linspace(0, 130, 131)
    analytic = model_qber(dense_dist, p_e, visibility, afterpulse_prob)

    budget_mode = (f"--target-sifted {target_sifted}" if target_sifted
                   else f"--bits {num_bits}")
    invocation = (f"--seed {seed} {budget_mode}"
                  + ("" if visibility == 1.0 else f" --visibility {visibility}")
                  + ("" if p_e == P_E else f" --p-e {p_e:g}")
                  + ("" if afterpulse_prob == AFTERPULSE_PROB
                     else f" --afterpulse {afterpulse_prob:g}"))

    def provenance(comment, marking=True):
        """Parameter provenance, shared by the .tex and .csv headers.

        Records which numbers came from the paper and which did not,
        because that is the one fact a reader cannot recover from the data.
        Parameters, sources, units.  It does not argue: the reasoning for
        each choice lives beside its constant above, where it can be kept
        correct, rather than being copied into every artifact.

        `marking=False` omits the (m)/(i) paragraph, for callers that print
        their own legend next to the table it describes.
        """
        c = comment
        lines = [
            "Gobby et al. (2004) link budget -- from the paper, not fitted:",
            f"  alpha   {ALPHA_dB} dB/km       mu     {MU} photons/clock",
            f"  eta_Bob {ETA_BOB}         (includes Bob's 5 dB apparatus loss)",
            f"  gate    {GATE_WIDTH*1e9:.1f} ns         clock  {REP_RATE/1e6:g} MHz"
            f"       pulse  {PULSE_WIDTH*1e12:.0f} ps",
            f"  P_e     {p_e:.2e}/clock/detector -- measured total error",
            "          probability: dark count 3.2e-7 + clock-laser stray",
            "          light 5.3e-7, lumped because the sum is what was measured",
            "",
            "Not from Gobby, both at datasheet value and NOT to be fitted --",
            "fitting either would make the 122 km point a fit, not a prediction:",
            f"  afterpulse_prob {afterpulse_prob}   dead_time {DEAD_TIME*1e6:g} us   [ID230]",
            "",
            f"Visibility is an output of this budget, V = S/(S + 2*P_e), and the",
            f"decoder is ideal at V = {visibility}.  The V column is that prediction;",
            "the paper states > 0.99 at 65 km and 0.884 at 122 km.",
        ]
        if marking:
            lines += [
                "",
                "The Gobby column is a MEASUREMENT only on rows marked (m); the",
                f"paper published four points, {', '.join(f'{d:g}' for d in GOBBY_DIST_KM)} km.  Rows marked (i)",
                "are np.interp between them, or below 4.4 km the endpoint",
                "repeated, and cannot disagree with anything.",
            ]
        return ''.join(f"{c} {ln}\n" if ln else f"{c}\n" for ln in lines)

    # Save results table.  An under-powered run writes to its own name so
    # it cannot replace the quotable table -- see _stem.
    table_path = os.path.join(os.path.dirname(__file__),
                              _stem('val_gobby_table', allow_underpowered)
                              + '.tex')
    with open(table_path, 'w') as f:
        f.write(f"% Generated by analysis/val_gobby/validate_gobby.py "
                f"{invocation}\n")
        f.write("% QBER error bars are binomial sqrt(q(1-q)/n_sifted).\n%\n")
        # marking=False: this file prints its own (m)/(m@x)/(i) legend below
        # the tabular, where a reader of the table needs it. The CSV has no
        # such place, so it keeps the paragraph in its header.
        f.write(provenance('%', marking=False))
        f.write("%\n")
        f.write(r"\begin{tabular}{rcccccc}" + "\n")
        f.write(r"  Distance & Pulses & Sifted & QBER & Predicted & This work "
                r"& Gobby et al. \\" + "\n")
        f.write(r"  (km)     &        & bits   & (\%) & $V$ & analytic (\%) "
                r"& (\%) \\" + "\n")
        f.write(r"\hline" + "\n")
        for d, q, s, n in zip(sim_dist, sim_qber, sim_sifted, sim_pulses):
            g_model = model_qber(d, p_e, visibility, afterpulse_prob)
            sd = qber_err_pp(q / 100.0, s)
            # Measurement vs interpolation, explicitly.  Without this the
            # last column reads as nine experimental points when the paper
            # published four.  Measured rows print the PUBLISHED value at
            # the published distance, not np.interp evaluated at ours.
            near = gobby_nearest_measured(d)
            if near is not None:
                zg, qg = near
                at = "" if abs(zg - d) < 1e-9 else f"@{zg:g}"
                cell = f"{qg:.1f}\\,(m{at})"
            else:
                cell = f"{gobby_measured_qber(d):.1f}\\,(i)"
            f.write(f"  {d:.0f} & {n:,} & {s} & ${q:.2f} \\pm {sd:.2f}$ "
                    f"& {predicted_visibility(d, p_e):.4f} "
                    f"& {g_model:.1f} & {cell} \\\\\n")
        f.write(r"\end{tabular}" + "\n")
        f.write("%\n% (m)   = Gobby measured here; the value shown is the "
                "PUBLISHED one.\n"
                "% (m@x) = same, but published at x km rather than our "
                "sweep distance\n"
                "%         (our 4 km vs their 4.4; our 100 km vs their "
                "101).\n"
                "% (i)   = np.interp between published points, or below "
                "4.4 km the\n%         endpoint repeated. NOT a "
                "measurement -- an (i) row cannot\n%         agree or "
                "disagree with the experiment.\n"
                "% Published points: "
                + ", ".join(f"{d:g} km -> {q:g}%"
                            for d, q in zip(GOBBY_DIST_KM, GOBBY_QBER))
                + ".\n")
    print(f"\nTable saved to {table_path}")

    # CSV alongside the .tex so the numbers are machine-readable
    csv_path = os.path.join(os.path.dirname(__file__),
                            _stem(f'val_gobby--seed{seed}',
                                  allow_underpowered) + '.csv')
    with open(csv_path, 'w') as f:
        f.write(f"# validate_gobby.py {invocation}\n#\n")
        f.write(provenance('#'))
        f.write("#\n")
        f.write("distance_km,pulses,sifted_bits,qber_pct,qber_err_pp,"
                "signal_per_clock,predicted_visibility,"
                "analytic_pct,gobby_pct,gobby_is_measured\n")
        for d, q, s, n in zip(sim_dist, sim_qber, sim_sifted, sim_pulses):
            f.write(f"{d:g},{n},{s},{q:.4f},{qber_err_pp(q/100.0, s):.4f},"
                    f"{signal_click_prob(d):.6e},"
                    f"{predicted_visibility(d, p_e):.6f},"
                    f"{model_qber(d, p_e, visibility, afterpulse_prob):.4f},"
                    f"{gobby_measured_qber(d):.4f},"
                    f"{int(gobby_is_measured(d))}\n")
    print(f"CSV saved to {csv_path}")

    # Plot
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(dense_dist, analytic, 'b-',
                label=f'This work, analytic ($P_e$={p_e:.1e}/clock, '
                      f'a={afterpulse_prob:g})')
        mc_label = (f'Monte Carlo (>={target_sifted} sifted/pt)'
                    if target_sifted else f'Monte Carlo ({num_bits:,} pulses)')
        yerr = [qber_err_pp(q / 100.0, s)
                for q, s in zip(sim_qber, sim_sifted)]
        ax.errorbar(sim_dist, sim_qber, yerr=yerr, fmt='ro-', capsize=3,
                    label=mc_label)
        ax.plot(GOBBY_DIST_KM, GOBBY_QBER, 'gs', label='Gobby paper (Fig 3)')
        ax.set_xlabel('Distance (km)')
        ax.set_ylabel('QBER (%)')
        # Self-describing, because a figure travels away from the script
        # that made it. Which paper, which seed, what budget.
        ax.set_title('Replication of Gobby, Yuan & Shields (2004): '
                     'time-bin BB84 QBER vs distance\n'
                     f'seed {seed}, {budget_mode}, link budget from the '
                     'paper with nothing fitted')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.axhline(y=11.0, color='r', linestyle='--', alpha=0.5, label='BB84 threshold (11%)')
        ax.set_xlim(0, 130)
        ax.set_ylim(0, 50)

        png_path = os.path.join(os.path.dirname(__file__),
                                _stem(f'val_gobby--seed{seed}',
                                      allow_underpowered) + '.png')
        fig.savefig(png_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"Figure saved to {png_path}")
    except ImportError:
        print("matplotlib not available — skipping figure")

    rc = check_results(sim_dist, sim_qber, sim_sifted, p_e, visibility,
                       afterpulse_prob, allow_underpowered)
    print("\nDone.")
    return rc


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Gobby et al. 2004 QBER validation")
    parser.add_argument('--bits', type=int, default=200000,
                        help='Flat pulses per point. Ignored when '
                             '--target-sifted is given.')
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--visibility', type=float, default=VISIBILITY,
                        help='DIAGNOSTIC ONLY. Decoder visibility (default '
                             '1.0 = ideal). V is an OUTPUT of this link '
                             'budget, V = S/(S+2*P_e); injecting one as well '
                             'double-counts the error physics. Overriding '
                             'it is a debugging '
                             'aid, not a way to match the paper.')
    parser.add_argument('--target-sifted', type=int, default=None,
                        metavar='N',
                        help=f'Grow the pulse count per distance until N '
                             f'sifted bits are collected, instead of a flat '
                             f'--bits budget. Recommended: '
                             f'{TARGET_SIFTED_DEFAULT} (sigma ~0.5 pp). A '
                             f'flat budget starves the long distances, which '
                             f'are the contested ones.')
    parser.add_argument('--p-e', type=float, default=P_E, dest='p_e',
                        metavar='P',
                        help=f'Background click probability per clock per '
                             f'detector (default {P_E:g} = Gobby\'s measured '
                             f'total error probability: dark count 3.2e-7 + '
                             f'clock-laser stray light 5.3e-7). Replaces the '
                             f'old --dcr, which was a category error -- the '
                             f'larger half of this term is stray light, not '
                             f'dark counts. Pass 0 to '
                             f'isolate the afterpulse floor.')
    parser.add_argument('--afterpulse', type=float, default=AFTERPULSE_PROB,
                        metavar='A',
                        help=f'SPAD afterpulse probability (default '
                             f'{AFTERPULSE_PROB:g}, the ID230 datasheet '
                             f'value). This is the candidate for Gobby\'s '
                             f'acknowledged third error source. Pass 0 for '
                             f'the control run that isolates the link '
                             f'budget. Do NOT fit it.')
    parser.add_argument('--ceiling', type=int, default=CEILING_DEFAULT,
                        help='Hard cap on pulses per point under '
                             '--target-sifted (default 500M)')
    parser.add_argument('--min-sifted', type=int, default=MIN_SIFTED,
                        help=f'Refuse to write tables if any row has fewer '
                             f'sifted bits than this (default {MIN_SIFTED})')
    parser.add_argument('--quick', action='store_true',
                        help='Smoke run: a small flat budget with the '
                             'statistical-power guard switched off, writing '
                             'to --quick artifact names so it cannot replace '
                             'the quotable table. Exercises the whole chain '
                             'end to end in about 20 s. The numbers it '
                             'produces are not citable -- that is the point '
                             'of the separate names.')
    parser.add_argument('--allow-underpowered', action='store_true',
                        help='Emit tables even when rows fall below '
                             '--min-sifted. For smoke runs only; the '
                             'artifact must not be cited.')
    args = parser.parse_args()

    # --quick is a smoke run: small flat budget, guard off, separate names.
    # It overrides --target-sifted rather than combining with it, because
    # growing the budget to hit a sifted target is the expensive thing the
    # flag exists to avoid.
    num_bits = args.bits
    target_sifted = args.target_sifted
    allow_underpowered = args.allow_underpowered
    if args.quick:
        num_bits = min(args.bits, 200_000)
        target_sifted = None
        allow_underpowered = True

    try:
        rc = run_validation(num_bits=num_bits, seed=args.seed,
                            visibility=args.visibility,
                            target_sifted=target_sifted,
                            ceiling=args.ceiling,
                            min_sifted=args.min_sifted,
                            allow_underpowered=allow_underpowered,
                            p_e=args.p_e, afterpulse_prob=args.afterpulse)
    except RuntimeError as exc:
        print(f"\nERROR: {exc}", file=sys.stderr)
        sys.exit(1)
    sys.exit(rc)
