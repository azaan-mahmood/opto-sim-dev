"""Replicate Gobby, Yuan & Shields (2004) — QBER vs distance.

References
----------
[1] Gobby, C., Yuan, Z. L., & Shields, A. J. (2004). Quantum key
    distribution over 122 km of standard telecom fiber. Appl. Phys.
    Lett. 84(19), 3762-3764.
"""
import argparse
import numpy as np
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from src.protocols.bb84_time_bin import simulate_bb84_time_bin

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
AFTERPULSE_PROB = 0.05

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
VISIBILITY = 1.0

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
    """Signal click probability per clock cycle, S = mu * T_link * eta_Bob."""
    return MU * 10.0 ** (-ALPHA_dB * np.asarray(dist_km) / 10.0) * ETA_BOB


def predicted_visibility(dist_km, p_e=P_E):
    """Fringe visibility as Gobby's own closed form, V = S/(S + 2*P_e).

    This is an *output* of the link budget, not a tunable input -- see the
    VISIBILITY note above and section 18.4.  Reproduces both visibilities
    the paper states: >0.99 at 65 km (0.9925 here) and 0.884 at 122 km
    (0.9058 here).
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


def simulate_qber(dist_km, num_bits, seed=42, verbose=False,
                  visibility=VISIBILITY, p_e=P_E,
                  afterpulse_prob=AFTERPULSE_PROB):
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
