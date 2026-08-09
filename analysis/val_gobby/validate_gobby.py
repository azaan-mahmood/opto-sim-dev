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

# --- Gobby's link budget (GOBBY-1) ---------------------------------------
# Every value in this block is taken from [1].  None of them is fitted.
#
# The previous parameterisation was wrong at the link-budget level, which
# is documented as GOBBY-1 in section 18 of opto-sim-issues-and-fixes.md.
# It ran alpha = 0.182, eta = 0.10 and a 1 ns gate against a 15 Hz DCR:
# the signal came out 3.68x too high and the error rate 57x too low, so
# the error/signal ratio -- the only quantity QBER measures -- was off by
# 209x.  That, not missing physics, is why the 9th-pass sweep came out
# flat (slope +0.011 pp/100 km) and could never have reproduced Gobby's
# rise from 3.3 % to 8.9 %.
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
# independently constrain.  Note that the larger half is stray light, so
# calling this a "dark count rate" is a category error -- the 9th pass
# made exactly that mistake, fitted 1,788 Hz, and described it as "119x
# the ID230 spec".  That characterisation is withdrawn (section 18.3):
# corrected for the 3.68x signal error the required rate is 4.85e-7,
# within 1.75x of the measured 8.5e-7, and no fit is needed at all.
P_E = 8.5e-7

# The SPAD model expresses background counts as a rate, so convert once
# here.  This is a lumped error-count term wearing the detector model's
# units -- 242.9 Hz -- and is not a physical dark count rate.
BACKGROUND_RATE_HZ = P_E / GATE_WIDTH

# Dead time is NOT from Gobby -- the paper does not state one.  13 us is
# the ID230 figure, kept so the detector model stays self-consistent.
# Flagged because it is the second of only two non-Gobby inputs here.
DEAD_TIME = 13e-6

# Afterpulsing is the first, and it is the open question of section 18.5.
# Gobby's QBER carries a third error source beyond the dark count and the
# stray light, worth 2.9-4.2 pp and roughly distance-independent (see the
# residual column in section 18.4).  Afterpulsing is the obvious
# candidate: it is already in the SPAD model and was measured here in
# isolation at e_ap = 2.262 +/- 0.176 % for afterpulse_prob = 0.05.
#
# 0.05 is the ID230 datasheet value, used as-is.  It is NOT fitted, and
# it must not become fitted -- tuning it against Gobby's QBER would
# reintroduce precisely the free parameter GOBBY-1 removes, and would
# make the 122 km point a fit again instead of a prediction.  Pass
# --afterpulse 0 for the control run that isolates the link budget.
#
# RESOLVED (GOBBY-7b): the value above is 0, not the ID230 datasheet 0.05.
#
# §19.5 already settled this from the paper, on three independent grounds:
# Fig. 3's dashed curve starts at ~0 where afterpulsing would put it at
# p_ap/2 = 2.5 %; the stated P_e = 8.5e-7 is dark + stray with no
# afterpulse term; and the closing summary enumerates three mechanisms,
# none of them afterpulsing.  Physically consistent too -- at 2 MHz gating
# with a 13 us dead time the detector is off while trapped carriers
# release.
#
# It was nonetheless left at 0.05, because until GOBBY-6 afterpulsing was
# the only thing in the chain supplying a floor at all.  Now that the
# modulation error it was standing in for is implemented, carrying both
# DOUBLE-COUNTS: on the old default the 0 km QBER read 5.806 % against
# Gobby's 3.3 %, the excess being exactly the 2.5 % §19.5 predicts.  That
# is the "pattern warning, third occurrence" §19.5 records against itself.
#
# Ownership is unchanged and correct: `afterpulse_prob` is a SPAD
# parameter (spad.py, default 0.05 = ID230), and this constant is the
# replication-level override -- a claim about Gobby's apparatus, not about
# SPADs in general.  The datasheet value stays reachable via --afterpulse.
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
# merely in value (section 18.4).  The decoder is therefore ideal and
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
# linewidth is carried but reported as negligible -- see GOBBY-7.
VISIBILITY = 1.0

# --- Modulation error (GOBBY-2 section 19.1 / GOBBY-6) -------------------
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
# and must not be presented as one.  See GOBBY-7.
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

# DIAGNOSTIC ONLY -- the jitter path is a recorded NEGATIVE RESULT, not an
# alternative parameterisation.  Drive noise was tested as the source of the
# 3.3 % floor and ruled out twice over (GOBBY-7 §24.1, §24.4):
#
#   * it would need 21.17 deg = 0.457 V = 11.8 % of V_pi, one to two orders
#     above what drive electronics deliver;
#   * measured head to head at 0 km, static bias gives 3.253 +/- 0.258 %
#     against jitter's 3.579 +/- 0.267 %, with the stated floor at 3.300 %.
#
# Kept so the negative control stays reachable and re-runnable:
#     simulate_qber(0, N, phase_error_rad=0.0, phase_noise_rad=PHASE_NOISE_RAD)
# Never a default here or in any component.
PHASE_NOISE_RAD = math.sqrt(-2.0 * math.log(1.0 - 2.0 * E_MOD))  # jitter 21.17

# --- Interferometer arm-length drift (GOBBY-7) ---------------------------
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
# This corrects a double-count introduced in GOBBY-7: PHASE_ERROR_RAD was
# `arccos(1 - 2*E_MOD)` = 20.93 deg, the bias-ONLY solution, while drift was
# separately applied on top.  Harmless while runs were ~1 s (drift ~0.05 deg),
# but once run_duration was set to the paper's own 120 s it put the 0 km floor
# at 4.30 % against their 3.3 %.  Same category of error as the afterpulse
# double-count in §19.5 -- two mechanisms, one measured aggregate, counted
# twice.  See GOBBY-7c.
PHASE_ERROR_RAD = _bias_for_aggregate(E_MOD, PHASE_DRIFT_RAD_S, KEY_TRANSFER_S)

# --- Laser linewidth (GOBBY-6) -------------------------------------------
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
# inert; see the contribution budget in GOBBY-6 for what it is worth.
LINEWIDTH = 3.0e6          # Hz, mid-range DFB assumption
PATH_MISMATCH = 0.0        # s, residual after trimming (0 = perfectly trimmed)


# --- OPEN-2: statistical power guards ------------------------------------
# A flat pulse budget is badly matched to this sweep: the sifted fraction
# runs from ~1.1e-3 at 0 km to ~1e-5 at 122 km, so a flat count
# over-samples the short end and starves the long end -- which is the only
# end anyone scrutinises.  The 122 km point once landed within 0.6 pp of
# Gobby's published value on **12 sifted bits** and was briefly read as
# corroboration; the Clopper-Pearson interval on 1/12 runs from 0.21 % to
# 38.5 %.  At full power the same point read 5.30 %.  That is what these
# guards exist to prevent.
#
# MIN_SIFTED sits below TARGET_SIFTED_DEFAULT on purpose -- see the same
# note in val_system_scenarios.py.  Equal values make the run a knife-edge
# where any undershoot discards the whole sweep at the final write.
MIN_SIFTED = 500
TARGET_SIFTED_DEFAULT = 3000      # sigma ~ 0.5 pp at QBER ~8 %

# Both raised for the GOBBY-1 parameterisation.  Correcting the link
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

    GOBBY-7 §24.5: this factor was **missing**, which put this function
    41% above what the chain delivers and made every quantity built on it
    (`model_qber`, `predicted_visibility`) optimistic.  Restoring it is not
    fitting -- `T_INT` is derived from the stated split ratio, and the same
    value already sits in `src/analytic/gobby_model.py` as `MU_EFF`.

    ACCURACY AGAINST THE MONTE CARLO, measured with dead time and
    afterpulsing off:

        0 km   S_mc / S_here = 1.019 +/- 0.012   [7,056 sifted]
        40 km                  1.011 +/- 0.025   [1,664 sifted]
        65 km                  0.985 +/- 0.031   [1,025 sifted]

    Agreement to ~2%, consistent with unity across the range.

    It read 0.954 / 1.012 / 0.994 until the SPAD detection probability was
    corrected: `spad.detect` computed `eta*(1 - exp(-mu))` where photons are
    detected independently and the right form is `1 - exp(-eta*mu)`.  That
    was the residual disagreement this docstring previously recorded as
    bounded-but-unattributed, and it is now closed (GOBBY-7b §24.5).

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
    VISIBILITY note above and section 18.4.  Reproduces both visibilities
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
    Do NOT fit detector parameters against this function.  The 8th pass
    did exactly that and it produced 7,593 Hz, which made the Monte Carlo
    overshoot to 17.5 % at 122 km against Gobby's 8.9 %.  Fitting the same
    endpoint against the MC directly gave 1,788 Hz and reproduced the
    sweep to a mean residual of 0.36 pp -- and GOBBY-1 then showed that
    *both* numbers were artifacts of a mis-specified link budget, and that
    the paper's own measured P_e removes the need to fit anything.  The
    accuracy figures previously quoted here were measured against that
    superseded parameterisation and have been dropped rather than
    restated; they described a chain that no longer exists.

    The structural fix in this function -- weighting every click channel
    into a single ratio -- stands on its own and is independent of the
    parameter question (section 18.6).  The form it replaced,
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
    BLOCK-1, one step milder: a derived column presented as experimental
    data.  It is not a copy-paste slip here -- the values do come from the
    paper -- but a reader comparing nine rows would believe there are nine
    measurements to disagree with, and there are four.

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


# --- Cross-checks against the paper's prose (GOBBY-7) --------------------
# Both bounds below come from the text rather than Fig. 3, so they test the
# budget against something it was not built from.  They are REPORTED, never
# used to adjust a parameter: tuning P_E until the 65 km check passed would
# make it fitted, which is exactly what GOBBY-1 removed.

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
        print(f"           over by {v / ERRONEOUS_BOUND_PCT - 1.0:+.0%}; the "
              f"bound would allow P_E <= {implied:.2e} (carrying {p_e:.2e}).")
        print(f"           Left as-is deliberately -- a discrepancy in our "
              f"budget to report, not a knob")
        print(f"           to turn. Tuning P_E here would re-fit what "
              f"GOBBY-1 unfitted. See GOBBY-7.")
    if s_65 is None:
        print(f"           (first-order S agrees with the chain to ~2% "
              f"since the missing T_INT and")
        print(f"           the SPAD Poisson form were fixed (GOBBY-7b "
              f"§24.5), so this verdict is")
        print(f"           reliable -- the miss is in the budget, not in "
              f"the comparison.)")
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
        # See GOBBY-2 section 19.7 in opto-sim-issues-and-fixes.md.
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
        # experiment, not a longer one.  See GOBBY-7c.
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
    """Refuse to emit a table whose rows cannot support a claim (OPEN-2).

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
        f"without being one (see OPEN-2 in opto-sim-issues-and-fixes.md; "
        f"the 122 km point once landed within 0.6 pp of the published "
        f"value on 12 sifted bits, purely by chance).\n"
        f"Fix: re-run with --target-sifted {TARGET_SIFTED_DEFAULT} (grows the "
        f"pulse count per distance), or raise --bits.\n"
        f"To emit anyway for a smoke test, pass --allow-underpowered."
    )


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

    # OPEN-2 guard: refuse to emit an artifact that looks like a result
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

    def provenance(comment):
        """Parameter provenance, shared by the .tex and .csv headers.

        Written into every artifact because the single most important
        fact about these numbers is which of them came from the paper
        (all the link parameters) and which did not (two detector
        parameters) -- see GOBBY-1, section 18.
        """
        c = comment
        lines = [
            "Gobby et al. (2004) link budget -- from the paper, not fitted:",
            f"  alpha = {ALPHA_dB} dB/km   mu = {MU} photons/clock   "
            f"eta_Bob = {ETA_BOB} (incl. Bob's 5 dB)",
            f"  gate = {GATE_WIDTH*1e9:.1f} ns   clock = {REP_RATE/1e6:g} MHz   "
            f"pulse = {PULSE_WIDTH*1e12:.0f} ps",
            f"  P_e = {p_e:.2e}/clock/detector -- Gobby's measured total error",
            "        probability (dark count 3.2e-7 + 1.3 um clock-laser stray",
            "        light 5.3e-7, deliberately lumped: their sum is what was",
            "        measured, and splitting them adds an unconstrained knob).",
            "",
            "Detector parameters NOT from Gobby (the only free inputs here):",
            f"  afterpulse_prob = {afterpulse_prob}   dead_time = "
            f"{DEAD_TIME*1e6:g} us   [ID230 datasheet]",
            "  Afterpulsing is the candidate for the paper's acknowledged",
            "  third error source (section 18.4 residual: 2.9-4.2 pp, roughly",
            "  distance-independent; measured here in isolation at 2.26 %).",
            "  It is used at its datasheet value and MUST NOT be fitted --",
            "  doing so would make the 122 km point a fit rather than a",
            "  prediction, which is the whole point of this parameterisation.",
            "",
            "Visibility is an OUTPUT, V = S/(S + 2*P_e), not an input; the",
            f"decoder is ideal (V = {visibility}).  Injecting a measured",
            "visibility as well would apply the same physics twice -- see",
            "section 18.4.  The V column below is that prediction; the paper",
            "states > 0.99 at 65 km and 0.884 at 122 km.",
            "",
            "Structural difference worth recording: Gobby's 1.6:1 reference-",
            "to-encoded intensity ratio puts 0.04 photons in the encoded",
            "pulse, whereas this encoder AMZI splits 50:50.  That is a",
            "difference in the chain, not a parameter, and may account for",
            "part of any residual.",
            "",
            "The Gobby column is a MEASUREMENT only on rows marked (m).",
            f"The paper published four points -- {', '.join(f'{d:g}' for d in GOBBY_DIST_KM)} km",
            "-- and rows marked (i) are np.interp between them, or below",
            "4.4 km the endpoint value repeated.  Compare against the (m)",
            "rows; the (i) rows cannot disagree with anything.",
        ]
        return ''.join(f"{c} {ln}\n" if ln else f"{c}\n" for ln in lines)

    # Save results table
    table_path = os.path.join(os.path.dirname(__file__), 'val_gobby_table.tex')
    with open(table_path, 'w') as f:
        f.write(f"% Generated by analysis/val_gobby/validate_gobby.py "
                f"{invocation}\n")
        f.write("% QBER error bars are binomial sqrt(q(1-q)/n_sifted).\n%\n")
        f.write(provenance('%'))
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
                            f'val_gobby--seed{seed}.csv')
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
        ax.set_title('Time-bin BB84 — QBER vs Distance')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.axhline(y=11.0, color='r', linestyle='--', alpha=0.5, label='BB84 threshold (11%)')
        ax.set_xlim(0, 130)
        ax.set_ylim(0, 50)

        png_path = os.path.join(os.path.dirname(__file__), f'val_gobby--seed{seed}.png')
        fig.savefig(png_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"Figure saved to {png_path}")
    except ImportError:
        print("matplotlib not available — skipping figure")

    print("\nDone.")


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
                             'double-counts the error physics (GOBBY-1 '
                             'section 18.4). Overriding it is a debugging '
                             'aid, not a way to match the paper.')
    parser.add_argument('--target-sifted', type=int, default=None,
                        metavar='N',
                        help=f'Grow the pulse count per distance until N '
                             f'sifted bits are collected, instead of a flat '
                             f'--bits budget. Recommended: '
                             f'{TARGET_SIFTED_DEFAULT} (sigma ~0.5 pp). A '
                             f'flat budget starves the long distances, which '
                             f'are the contested ones -- see OPEN-2.')
    parser.add_argument('--p-e', type=float, default=P_E, dest='p_e',
                        metavar='P',
                        help=f'Background click probability per clock per '
                             f'detector (default {P_E:g} = Gobby\'s measured '
                             f'total error probability: dark count 3.2e-7 + '
                             f'clock-laser stray light 5.3e-7). Replaces the '
                             f'old --dcr, which was a category error -- the '
                             f'larger half of this term is stray light, not '
                             f'dark counts (GOBBY-1 section 18.3). Pass 0 to '
                             f'isolate the afterpulse floor.')
    parser.add_argument('--afterpulse', type=float, default=AFTERPULSE_PROB,
                        metavar='A',
                        help=f'SPAD afterpulse probability (default '
                             f'{AFTERPULSE_PROB:g}, the ID230 datasheet '
                             f'value). This is the candidate for Gobby\'s '
                             f'acknowledged third error source. Pass 0 for '
                             f'the control run that isolates the link '
                             f'budget. Do NOT fit it -- see section 18.5.')
    parser.add_argument('--ceiling', type=int, default=CEILING_DEFAULT,
                        help='Hard cap on pulses per point under '
                             '--target-sifted (default 500M)')
    parser.add_argument('--min-sifted', type=int, default=MIN_SIFTED,
                        help=f'Refuse to write tables if any row has fewer '
                             f'sifted bits than this (default {MIN_SIFTED})')
    parser.add_argument('--allow-underpowered', action='store_true',
                        help='Emit tables even when rows fall below '
                             '--min-sifted. For smoke runs only; the '
                             'artifact must not be cited.')
    args = parser.parse_args()
    try:
        run_validation(num_bits=args.bits, seed=args.seed,
                       visibility=args.visibility,
                       target_sifted=args.target_sifted,
                       ceiling=args.ceiling,
                       min_sifted=args.min_sifted,
                       allow_underpowered=args.allow_underpowered,
                       p_e=args.p_e, afterpulse_prob=args.afterpulse)
    except RuntimeError as exc:
        print(f"\nERROR: {exc}", file=sys.stderr)
        sys.exit(1)
